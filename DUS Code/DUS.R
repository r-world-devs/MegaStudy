
# Define output folder ----
outputFolder <- here::here("storage")   
# Create output folder if it doesn't exist
if (!file.exists(outputFolder)){
  dir.create(outputFolder, recursive = TRUE)}


# Start log ----
log_file <- here::here(outputFolder, paste0(dbName, "_log.txt"))
logger <- create.logger()
logfile(logger) <- log_file
level(logger) <- "INFO"

# Create cdm object ----
info(logger, 'CREATE CDM OBJECT')
cdm <- cdmFromCon(
  con = db,
  cdmSchema = c(schema = cdmSchema),
  writeSchema = c(schema = writeSchema, prefix = writePrefix),
  cdmName = dbName
)


# cdm snapshot ----
info(logger, 'CREATE SNAPSHOT')
write.csv(
  x = snapshot(cdm),
  file = here(outputFolder, paste0("snapshot_", cdmName(cdm), ".csv")),
  row.names = FALSE
)


# if SIDIAP filter cdm$drug_exposure
if (cdmName(cdm) == "SIDIAP") {
  info(logger, "FILTER DRUG EXPOSURE TABLE SIDIAP")
  cdm$drug_exposure <- cdm$drug_exposure %>%
    filter(drug_type_concept_id == 32839) %>%
    compute()
}


## Generate indication concepts ----
info(logger, "GENERATE INDICATION CONCEPTS")
source(here("individual_indications.R"))

## Generate indication cohorts ----
info(logger, "GENERATE INDICATION COHORTS")

covid <- CDMConnector::readCohortSet(
  here::here("json_cohort"))

cdm <- CDMConnector::generateCohortSet(cdm,
                                       covid,
                                       name = "covid",
                                       overwrite = TRUE)

attr(cdm$covid, "cohort_attrition") <- attr(cdm$covid, "cohort_attrition") %>% 
  mutate(
    number_subjects = as.integer(number_subjects),
    number_records = as.integer(number_records)
  )

name_indication <- c(
  "neutropenia",
  "bacteraemia",
  "infection",
  "pneumonia",
  "orthopaedic_surgery",
  "cardiovascular_surgery",
  "gynecologic_surgery",
  "chronic_bronchitis",
  "cataract_surgery",
  "meningitis",
  "urogenital_surgery",
  "gastrointestinal_surgery",
  "colorectal_surgery",
  "endocarditis",
  "sepsis",
  "lyme_disease",
  "copd",
  "syphilis",
  "gonorrhea",
  "otitis_media",
  "osteomyelitis",
  "cellulitis_erysipelas_woundinfection",
  "pyelonephritis",
  "cystitis_bacteriuria",
  "sinusitis",
  "helicobacter_pylori",
  "typhoid",
  "tonsillitis_pharyngitis",
  "rheumatic_fever",
  "yaws_pinta",
  "gingivitis",
  "scarlet_fever",
  "urethritis",
  "cervicitis",
  "fabry",
  "gaucher",
  "assisted_reproduction",
  "implant_transplant",
  "apl",
  "hereditary_angioedema",
  "smoking_cessation",
  "rheumatoid_arthritis",
  "giant_cell_arteritis",
  "jia",
  "pcv_cnv",
  "csc",
  "exudative_amd",
  "catheter_flushing",
  "ischemic_stroke",
  "pe",
  "mi"
)

for (ind in name_indication) {
  
  ind_obj <- get(ind)
  ind_list <- list(as.integer(ind_obj$concept_id))
  names(ind_list) <- ind
  
  info(logger, paste0(ind))
  
  cdm <- generateConceptCohortSet(
    cdm,
    conceptSet = ind_list,
    name = ind,  
    limit = "all",
    requiredObservation = c(0, 0),
    end = 1,
    overwrite = TRUE
  )
  
  attr(cdm[[ind]], "cohort_attrition") <- attr(cdm[[ind]], "cohort_attrition") %>% 
    mutate(
      number_subjects = as.integer(number_subjects),
      number_records = as.integer(number_records)
    )
}

## Generate drug concepts ----
info(logger, "GENERATE DRUG CONCEPTS")
source(here("drugs.R"))


## get denominator population ----
info(logger, "GENERATE DENOMINATOR COHORT")
cdm <- generateDenominatorCohortSet(
  cdm = cdm, 
  name = "denominator",
  cohortDateRange = as.Date(c("2010-01-01",NA)),
  ageGroup = list(
    c(0, 150) 
  ),
  sex = c("Both"), 
  daysPriorObservation = c(0,30),   
  requirementInteractions = TRUE
)


## Generate incident drug cohorts---
info(logger, "GENERATE INCIDENT DRUG COHORTS")

inc_attrition <- tibble::as_tibble(NULL)
inc_pat <- tibble::as_tibble(NULL)
for (i in seq_along(concept_drugs)) {
  
  info(logger, paste0(names(concept_drugs[i])))
  
  cdm <- generateDrugUtilisationCohortSet(
    cdm,
    name = "drug_cohorts",
    conceptSet = concept_drugs[i],
    durationRange = c(1, Inf),
    imputeDuration = "none",
    gapEra = 0,
    priorUseWashout = 0,
    priorObservation = 0,
    cohortDateRange = as.Date(c("2010-01-01",NA)),
    limit = "all"
  )
  
  inc <- estimateIncidence(
    cdm = cdm,
    denominatorTable = "denominator",
    outcomeTable = "drug_cohorts",
    interval = "years",
    completeDatabaseIntervals = FALSE,
    outcomeWashout = 30,
    repeatedEvents = TRUE,
    minCellCount = 10,
    returnParticipants = TRUE
  ) 
  
  inc_attrition <- bind_rows(inc_attrition,attrition(inc))
  
  inc_analysis_id <- inc %>% 
    filter(denominator_days_prior_observation == "30") %>%
    distinct(analysis_id,outcome_cohort_name) %>% 
    mutate(analysis_id = as.integer(analysis_id))
  
  addition <- IncidencePrevalence::participants(result = inc, analysisId = inc_analysis_id$analysis_id) %>%
    filter(!is.na(outcome_start_date)) %>%
    mutate(
      calendar_year = lubridate::year(outcome_start_date)
    ) %>% 
    distinct(subject_id, calendar_year, .keep_all = TRUE)  %>%
    collect() %>%
    mutate(analysis_id = inc_analysis_id$analysis_id) %>%
    inner_join(inc_analysis_id, by = "analysis_id") %>%
    mutate(cohort_definition_id = analysis_id * i) %>%
    select(-"analysis_id") %>% 
    rename(study_end = cohort_end_date ,
           study_start = cohort_start_date ) %>%
    inner_join(cdm[["drug_cohorts"]] %>% select(-"cohort_definition_id"), 
               by = c("subject_id","outcome_start_date" = "cohort_start_date"), copy = TRUE) %>%
    select(!c("study_start","study_end")) %>%
    rename(cohort_start_date = outcome_start_date)
  
  inc_pat <- bind_rows(inc_pat, addition)
  
}

cdm <- insertTable(cdm, "inc_pat", inc_pat)
cdm$inc_pat <- omopgenerics::newCohortTable(cdm[["inc_pat"]]) 


info(logger, "WRITE INCIDENCE ATTRITION")
write.csv(
  inc_attrition,
  here("storage", paste0("incidence_attrition_", cdmName(cdm), ".csv"))
)


## characterisation of demographics of incident patients ---------------------
info(logger, "DEMOGRAPHICS INCIDENT PATIENTS")

inc_demographics <- cdm[["inc_pat"]]  %>%
  addDemographics(
    indexDate = "cohort_start_date",
    age = TRUE,
    ageName = "age",
    ageGroup = list(
      c(0,18),
      c(19,64),
      c(65,150)
    ),
    sex = TRUE,
    sexName = "sex",
    priorObservation = FALSE,
    futureObservation = FALSE,
  ) %>% 
  CohortCharacteristics::summariseCharacteristics(
    strata = list( "Calendar Year" = "calendar_year"),
    demographics = TRUE,
    ageGroup = list(
      c(0,18),
      c(19,64),
      c(65,150)
    ),
    # commented arguments are not used anymore, instead we have more detailed ones, all set to `list()` by default
    #tableIntersect = list(), 
    #cohortIntersect = list(),
    #conceptIntersect = list(),
    otherVariables = character()
  ) %>% 
  suppress(minCellCount = 10) %>%
  mutate(cohort_definition_id = as.integer(str_sub(group_level, 8))) %>%
  left_join(
    cdm[["inc_pat"]] %>% select(cohort_definition_id, outcome_cohort_name) %>% distinct(), 
    by = c("cohort_definition_id"), copy = TRUE
  ) %>%
  select(-"cohort_definition_id") %>%
  filter(variable_name %in% c("Number subjects","Number records","Age","Sex","Age group"))

write.csv(
  inc_demographics, 
  here("storage", paste0("inc_demographics_summary_", cdmName(cdm), ".csv"))
)

## characterisation of indication of incident patients ---------------------
info(logger, "INDICATION INCIDENT PATIENTS")

cov_indication <- c(name_indication,"covid")

source("applyFilterPatforInd.R")

inc_indication_summary <- tibble::as.tibble(NULL)
for (cov_ind in cov_indication){
  indication_drug_cohort <- applyFilterPatforInd(cdm[["inc_pat"]], cov_ind) %>%
    addIndication(
      indicationCohortName = cov_ind,
      indicationGap = c(0,30,Inf),
      unknownIndicationTable = NULL,
      indicationDate = "cohort_start_date"
    ) %>% 
    summariseIndication(
      strata = list("Calendar Year" = "calendar_year")
    ) %>% 
    suppress(minCellCount = 10) %>%
    mutate(cohort_definition_id = as.integer(str_sub(group_level, 8))) %>%
    inner_join(
      cdm[["inc_pat"]] %>% select(cohort_definition_id, outcome_cohort_name) %>% distinct(), 
      by = c("cohort_definition_id"), copy = TRUE
    ) %>%
    select(-"cohort_definition_id") %>% 
    filter(!is.na(estimate_value), estimate_value != "<5", variable_level != "None")
  
  inc_indication_summary <- bind_rows(inc_indication_summary,indication_drug_cohort)
}

write.csv(
  inc_indication_summary, 
  here("storage", paste0("inc_indication_summary_", cdmName(cdm), ".csv"))
)

## large scale characterisation of incident patients ---------------------
info(logger, "LARGE SCALE CHARACTERISATION INCIDENT PATIENTS")
inc_lsc <- cdm[["inc_pat"]]  %>%
  CohortCharacteristics::summariseLargeScaleCharacteristics(
    window = list(
      c(-Inf, -30), c(-Inf, -1),  c(-30, -1), c(0, 0), c(1, 7),c(8,Inf),c(1,Inf)
    ),
    strata = list("Calendar Year" = "calendar_year"),
    eventInWindow = c("condition_occurrence","observation","measurement"),       
    episodeInWindow = "drug_exposure",            
    indexDate = "cohort_start_date",
    censorDate = NULL,
    minimumFrequency = 0.05,                      
    excludedCodes = NULL
  ) %>%                     
  suppress(minCellCount = 10) %>%                
  mutate(cohort_definition_id = as.integer(str_sub(group_level, 8))) %>%
  inner_join(cdm[["inc_pat"]] %>% select(cohort_definition_id, outcome_cohort_name) %>% distinct(), by = c("cohort_definition_id"), copy = TRUE ) %>%
  select(-"cohort_definition_id") %>% 
  filter(!is.na(estimate_value), estimate_value != "<5")


write.csv(
  inc_lsc, 
  here("storage", paste0("inc_lsc_", cdmName(cdm), ".csv"))
)


## characterisation of drug use of incident patients ---------------------
info(logger, "DRUG USE INCIDENT PATIENTS")

source("applyFilterIngforUse.R")

ingredients <- c(
  "tenecteplase" = 19098548L,
  "alteplase" = 1347450L,
  "sarilumab" = 1594587L,
  "verteporfin" = 912803L,
  "tocilizumab" = 40171288L,
  "varenicline" = 780442L,
  "ceftolozane" = 45892599L,
  "tazobactam" = 1741122L,
  "C1 esterase inhibitor" = 45892906L,
  "arsenic trioxide" = 1333379L,
  "belatacept" = 40239665L,
  "ganirelix" = 1536743L,
  "tigecycline" = 1742432L,
  "imiglucerase" = 1348407L,
  "agalsidase beta" = 1525746L,
  "azithromycin" = 1734104L,
  "clarithromycin" = 1750500L,
  "penicillin V" = 1729720L,
  "penicillin G" = 1728416L,
  "amoxicillin" = 1713332L,
  "clavulanate" = 1759842L,
  "ceftriaxone" = 1777806L,
  "cefotaxime" = 1774470L,
  "meropenem" = 1709170L,
  "cefuroxime" = 1778162L,
  "piperacillin" = 1746114L,
  "streptokinase" = 19136187L,
  "urokinase" = 1307515L,
  "abatacept" = 1186087L,
  "tofacitinib" = 42904205L,
  "baricitinib" = 1510627L,
  "upadacitinib" = 1361580L,
  "etanercept" = 1151789L,
  "infliximab" = 937368L,
  "certolizumab pegol" = 912263L,
  "CERTOLIZUMAB" = 36857573L,
  "golimumab" = 19041065L,
  "anakinra" = 1114375L,
  "ranibizumab" = 19080982L,
  "bevacizumab" = 1397141L,
  "nicotine" = 718583L,
  "icatibant" = 40242044L,
  "ecallantide" = 40168938L,
  "Conestat alfa" = 36878937L,
  "lanadelumab" = 35200405L,
  "berotralstat" = 37003361L,
  "cytarabine" = 1311078L,
  "Cytarabine liposomal" = 902730L,
  "CYTARABINE 5'-PHOSPHATE" = 36863408L,
  "daunorubicin" = 1311799L,
  "idarubicin" = 19078097L,
  "tretinoin" = 903643L,
  "mycophenolic acid" = 19012565L,
  "sirolimus" = 19034726L,
  "cyclosporine" = 19010482L,
  "tacrolimus" = 950637L,
  "cetrorelix" = 1503983L,
  "velaglucerase alfa" = 40174604L,
  "taliglucerase alfa" = 42800246L,
  "Agalsidase alfa" = 36878851L
)

cdm$drug_strength <- cdm$drug_strength %>%
  dplyr::mutate(numerator_value = as.numeric(numerator_value))

inc_use_summary <- tibble::as_tibble(NULL)
for (j in seq_along(ingredients)) {
  ingredient <- unname(ingredients)[j]
  name <- names(ingredients)[j]
  
  info(logger, paste0(ingredient))
  
  use_drug_cohort <- applyFilterIngforUse(cdm[["inc_pat"]], name) %>%
    addDrugUse(
      ingredientConceptId = ingredient,
      duration = TRUE,
      quantity = TRUE,
      dose = TRUE,
      gapEra = 0,
      eraJoinMode = "subsequent",
      overlapMode = "subsequent",
      sameIndexMode = "sum",
      imputeDuration = "none",
      imputeDailyDose = "none"
    ) %>%
    summariseDrugUse(
      strata = list("Calendar Year" = "calendar_year")
      ) %>%  suppress(minCellCount = 10) %>%
    mutate(variable_level = name) %>%
    mutate(cohort_definition_id = as.integer(str_sub(group_level, 8))) %>%
    inner_join(
      cdm[["inc_pat"]] %>% select(cohort_definition_id, outcome_cohort_name) %>% distinct(), 
      by = c("cohort_definition_id"), copy = TRUE 
    ) %>%
    select(-"cohort_definition_id") %>%
    filter(variable_name %in% c("number subjects","number records","duration","number_exposures","cumulative_quantity",
                                "initial_quantity","initial_daily_dose_milligram","cumulative_dose_milligram"))
  
  inc_use_summary <- bind_rows(inc_use_summary,use_drug_cohort)
}

write.csv(
  inc_use_summary, 
  here("storage", paste0("inc_use_summary_", cdmName(cdm), ".csv"))
)

## Generate prevalent drug cohorts---------------------------------------------------------------------------------------------------------------------------
info(logger, "GENERATE PREVALENT DRUG COHORTS")

prev_attrition <- tibble::as_tibble(NULL)
prev_pat <- tibble::as_tibble(NULL)
for (i in seq_along(concept_drugs)){

  cdm <- generateDrugUtilisationCohortSet(
    cdm,
    name = "drug_cohorts",
    conceptSet = concept_drugs[i],
    durationRange = c(1, Inf),
    imputeDuration = "none",
    gapEra = 0,
    priorUseWashout = 0,
    priorObservation = 0,
    cohortDateRange = as.Date(c("2010-01-01",NA)),
    limit = "all"
  )
  
  prev <- estimatePeriodPrevalence(
    cdm = cdm,
    denominatorTable = "denominator",
    outcomeTable = "drug_cohorts",
    interval = "years",
    completeDatabaseIntervals = FALSE,
    minCellCount = 10,
    returnParticipants = TRUE    
  )
  
  prev_attrition <- bind_rows(prev_attrition,attrition(prev))
  
  prev_analysis_id <- prev %>% 
    filter(denominator_days_prior_observation == "0") %>%
    distinct(analysis_id,outcome_cohort_name) %>% 
    mutate(analysis_id = as.integer(analysis_id))
  
  addition <- IncidencePrevalence::participants(result = prev, analysisId = prev_analysis_id$analysis_id) %>%
    filter(!is.na(outcome_start_date)) %>%
    mutate(
      calendar_year = lubridate::year(outcome_start_date)
    ) %>% 
    distinct(subject_id, calendar_year, .keep_all = TRUE)  %>%
    collect() %>%
    mutate(analysis_id = prev_analysis_id$analysis_id) %>%
    inner_join(prev_analysis_id, by = "analysis_id") %>%
    mutate(cohort_definition_id = analysis_id * i) %>%
    select(-"analysis_id") %>% 
    rename(study_end = cohort_end_date ,
           study_start = cohort_start_date ) %>%
    inner_join(cdm[["drug_cohorts"]] %>% select(-"cohort_definition_id"), 
               by = c("subject_id","outcome_start_date" = "cohort_start_date"), copy = TRUE) %>%
    select(!c("study_start","study_end")) %>%
    rename(cohort_start_date = outcome_start_date)
  
  prev_pat <- bind_rows(prev_pat, addition)
  
}

info(logger, "WRITE PREVALENCE ATTRITION")
write.csv(
  prev_attrition,
  here("storage", paste0("prevalence_attrition_", cdmName(cdm), ".csv"))
)

## characterisation of prevalent patients ---------------------
info(logger, "CHARACTERISATION PREVALENT PATIENTS")

cdm <- insertTable(cdm, "prev_pat", prev_pat)
cdm$prev_pat <- omopgenerics::newCohortTable(cdm[["prev_pat"]])    


prev_demographics <- cdm[["prev_pat"]]  %>%
  addDemographics(
    indexDate = "cohort_start_date",
    age = TRUE,
    ageName = "age",
    ageGroup = list(
      c(0,18),
      c(19,64),
      c(65,150)
    ),
    sex = TRUE,
    sexName = "sex",
    priorObservation = FALSE,
    futureObservation = FALSE,
  ) %>% 
  CohortCharacteristics::summariseCharacteristics(
    strata = list( "Calendar Year" = "calendar_year"),
    demographics = TRUE,
    ageGroup = list(
      c(0,18),
      c(19,64),
      c(65,150)
    ),
    # commented arguments are not used anymore, instead we have more detailed ones, all set to `list()` by default
    # tableIntersect = list(),
    # cohortIntersect = list(),
    # conceptIntersect = list(),
    otherVariables = character()
  ) %>% suppress(minCellCount = 10) %>%
  mutate(cohort_definition_id = as.integer(str_sub(group_level, 8))) %>%
  left_join(cdm[["prev_pat"]] %>% select(cohort_definition_id, outcome_cohort_name) %>% distinct(), by = c("cohort_definition_id"), copy = TRUE ) %>%
  select(-"cohort_definition_id") %>%
  filter(variable_name %in% c("Number subjects","Number records","Age","Sex","Age group"))

write.csv(
  prev_demographics, 
  here("storage", paste0("prev_demographics_summary_", cdmName(cdm), ".csv"))
)

## characterisation of indication of prevalent patients ---------------------
info(logger, "INDICATION PREVALENT PATIENTS")

prev_indication_summary <- tibble::as.tibble(NULL)
for (cov_ind in cov_indication){
  indication_drug_cohort <- applyFilterPatforInd(cdm[["prev_pat"]], cov_ind) %>%
    addIndication(
      indicationCohortName = cov_ind,
      indicationGap = c(0,30,Inf),
      unknownIndicationTable = NULL,
      indicationDate = "cohort_start_date"
    ) %>% 
    summariseIndication(
      strata = list("Calendar Year" = "calendar_year")
    ) %>% suppress(minCellCount = 10) %>%
    mutate(cohort_definition_id = as.integer(str_sub(group_level, 8))) %>%
    inner_join(cdm[["prev_pat"]] %>% select(cohort_definition_id, outcome_cohort_name) %>% distinct(), by = c("cohort_definition_id"), copy = TRUE ) %>%
    select(-"cohort_definition_id") %>% 
    filter(!is.na(estimate_value), estimate_value != "<5",
           variable_level != "None")
  
  prev_indication_summary <- bind_rows(prev_indication_summary,indication_drug_cohort)
}

write.csv(
  prev_indication_summary, 
  here("storage", paste0("prev_indication_summary_", cdmName(cdm), ".csv"))
)


## large scale characterisation of prevalent patients ---------------------
info(logger, "LARGE SCALE CHARACTERISATION PREVALENT PATIENTS")
prev_lsc <- cdm[["prev_pat"]]  %>%
  CohortCharacteristics::summariseLargeScaleCharacteristics(
    window = list(
      c(-Inf, -30), c(-Inf, -1),  c(-30, -1), c(0, 0), c(1, 7),c(8,Inf),c(1,Inf)
    ),
    strata = list("Calendar Year" = "calendar_year"),
    eventInWindow = c("condition_occurrence","observation","measurement"),   
    episodeInWindow = "drug_exposure",           
    indexDate = "cohort_start_date",
    censorDate = NULL,
    minimumFrequency = 0.05,                    
    excludedCodes = NULL) %>%                  
  suppress(minCellCount = 10) %>%  
  mutate(cohort_definition_id = as.integer(str_sub(group_level, 8))) %>%
  inner_join(cdm[["prev_pat"]] %>% select(cohort_definition_id, outcome_cohort_name) %>% distinct(), by = c("cohort_definition_id"), copy = TRUE ) %>%
  select(-"cohort_definition_id") %>% 
  filter(!is.na(estimate_value), estimate_value != "<5")


write.csv(
  prev_lsc, 
  here("storage", paste0("prev_lsc_", cdmName(cdm), ".csv"))
)


## characterisation of drug use of prevalent patients ---------------------
info(logger, "DRUG USE PREVALENT PATIENTS")

prev_use_summary <- tibble::as_tibble(NULL)
for (j in seq_along(ingredients)) {
  ingredient <- ingredients[j]
  name <- names[j]  
  
  info(logger, paste0(ingredient))
  
  use_drug_cohort <- applyFilterIngforUse(cdm[["prev_pat"]], name) %>%
    addDrugUse(
      ingredientConceptId = ingredient ,
      duration = TRUE,
      quantity = TRUE,
      dose = TRUE,
      gapEra = 0,
      eraJoinMode = "subsequent",
      overlapMode = "subsequent",
      sameIndexMode = "sum",
      imputeDuration = "none",
      imputeDailyDose = "none"
    ) %>%
    summariseDrugUse(
      strata = list("Calendar Year" = "calendar_year")
    ) %>%  suppress(minCellCount = 10) %>% 
    mutate(variable_level = name) %>%
    mutate(cohort_definition_id = as.integer(str_sub(group_level, 8))) %>%
    inner_join(cdm[["prev_pat"]] %>% select(cohort_definition_id, outcome_cohort_name) %>% distinct(), by = c("cohort_definition_id"), copy = TRUE ) %>%
    select(-"cohort_definition_id") %>%
    filter(variable_name %in% c("number subjects","number records","duration","number_exposures","cumulative_quantity",
                                "initial_quantity","initial_daily_dose_milligram","cumulative_dose_milligram"))
  
  prev_use_summary <- bind_rows(prev_use_summary,use_drug_cohort)
}

write.csv(
  prev_use_summary, 
  here("storage", paste0("prev_use_summary_", cdmName(cdm), ".csv"))
)

## zip everything together ---
info(logger, "ZIP EVERYTHING")
zip(
  zipfile = here::here(paste0("Results_DUS_", cdmName(cdm), ".zip")),
  files = list.files(outputFolder),
  root = outputFolder
)

## remove storage, caused confusion ---
info(logger, "REMOVE STORAGE FOLDER")
unlink(here("storage"), recursive = TRUE)

print("Done!")
print("If all has worked, there should now be a zip file with your DUS results in the same level as the .Rproj")
print("Thank you for running the DUS analysis!")