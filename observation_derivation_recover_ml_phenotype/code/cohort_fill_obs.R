


#' Takes input and fills empty observation table
#' 
#' @param records the data for which to fill observation records
#' 
#' @return the returned observation table with records filled in
#' 

.fill_observation_records <- function(records) {
  
  template <- cdm_tbl('observation') %>% filter(0==1)
  template %>% full_join(records,by=tbl_vars(records),
                         copy=TRUE)
  
}


#' Takes cohort and fills observation records when input is 
#' structured as the `measurement` tbl
#' 
#' @param src_tbl the name of the src_tbl input
#' @param obs_concept_id the integer value to fill in for the `observation_concept_id`
#' in the output
#' 
#' @return 
#' 
#' 

create_observation_records_meas <- function(src_tbl,
                                            obs_concept_id) {
  
  
  src_tbl %>%
    mutate(observation_concept_id=obs_concept_id) %>%
    rename(observation_date=measurement_date,
           observation_datetime=measurement_datetime,
           observation_source_concept_id=measurement_id,
           observation_type_concept_id=measurement_type,
           qualifier_concept_id=measurement_concept_id) %>%
    select(person_id,
           visit_occurrence_id,
           observation_concept_id,
           value_as_concept_id,
           observation_date,
           observation_datetime,
           observation_source_concept_id,
           observation_type_concept_id,
           qualifier_concept_id) %>%
    .fill_observation_records()
  
}

#' Takes cohort and fills observation record when input is
#' structured as the `condition_occurrence` tbl
#' 
#' @param src_tbl the name of the src_tbl input
#' @param obs_concept_id the integer value to fill in for the `observation_concept_id`
#' in the output
#' 
#' @return 
#' 

create_observation_records_cond <- function(src_tbl,
                                            obs_concept_id) {
  
  
  src_tbl %>%
    mutate(observation_concept_id=obs_concept_id) %>%
    rename(observation_date=condition_start_date,
           observation_datetime=condition_start_datetime,
           observation_source_concept_id=condition_occurrence_id,
           observation_type_concept_id=condition_type,
           qualifier_concept_id=condition_qualifier_concept_id) %>%
    select(person_id,
           visit_occurrence_id,
           observation_concept_id,
           value_as_concept_id,
           observation_date,
           observation_datetime,
           observation_source_concept_id,
           observation_type_concept_id,
           qualifier_concept_id) %>%
    .fill_observation_records()
  
}


#' Takes cohort and fills observation record when input is
#' structured as the `immunization` tbl
#' 
#' @param src_tbl the name of the src_tbl input
#' @param obs_concept_id the integer value to fill in for the `observation_concept_id`
#' in the output
#' 
#' @return 
#' 

create_observation_records_imm <- function(src_tbl,
                                            obs_concept_id=2000001492L) {
  
  
  src_tbl %>%
    mutate(observation_concept_id=obs_concept_id) %>%
    rename(observation_date=immunization_date,
           observation_datetime=immunization_datetime,
           observation_source_concept_id=immunization_id,
           observation_type_concept_id=immunization_type,
           qualifier_concept_id=immunization_concept_id) %>%
    select(person_id,
           visit_occurrence_id,
           observation_concept_id,
           value_as_concept_id,
           observation_date,
           observation_datetime,
           observation_source_concept_id,
           observation_type_concept_id,
           qualifier_concept_id) %>%
    .fill_observation_records()
  
}

