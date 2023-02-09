


#' Identifies all COVID-19 PCR tests in the measurement table
#' 
#' @param pcr_codeset codeset that contains pcr tests
#' @param measurement_tbl labs from the measurement table
#' 
#' @return measurement table that contains only 
#' pcr tests; 
#' 
#' Function will perform string search on `measurement_source_value`
#' Function will add a column called `measurement_type` to 
#' annotate whether test was derived from structured code or 
#' source value derived
#' 

find_pcr <- function(pcr_codeset,
                     measurement_tbl=cdm_tbl('measurement_labs')) {
  
  
  all_tests <- 
    measurement_tbl %>%
    inner_join(
      select(pcr_codeset,concept_id),
      by=c('measurement_concept_id'='concept_id')
    ) %>% mutate(measurement_type=2000001517L)
  
  string_search <- 
    measurement_tbl %>%
    filter(measurement_concept_id == 0L &
             str_detect(
               lower(measurement_source_value),'sars-cov-2 osl|sars-cov-2, naa')
           ) %>% compute_new(temporary=TRUE,
                             indexes=list('measurement_concept_id',
                                          'person_id')) %>%
    mutate(measurement_type=2000001519L)
  
  combined <- 
    dplyr::union(all_tests,
                 string_search) %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  
}



#' Identifies all COVID-19 serology tests in the measurement table
#' 
#' @param serology_codeset codeset that contains serology tests
#' @param measurement_tbl labs from the measurement table
#' 
#' @return measurement table that contains only 
#' serology tests; 
#' 
#' Function will perform string search on `measurement_source_value`
#' Function will add a column called `measurement_type` to 
#' annotate whether test was derived from structured code or 
#' source value derived
#' 

find_serology <- function(serology_codeset,
                          measurement_tbl=cdm_tbl('measurement_labs')) {
  
  
  all_tests <- 
    measurement_tbl %>%
    inner_join(
      select(serology_codeset,concept_id),
      by=c('measurement_concept_id'='concept_id')
    ) %>% mutate(measurement_type=2000001517L)
  
  string_search <- 
    measurement_tbl %>%
    filter(measurement_concept_id == 0L &
             str_detect(
               lower(measurement_source_value),
               'sars-cov-2 ace2 blocking activity, igg|
               sars-cov-2 igg, spike antigen s/c ratio|
               diasorin sars-cov-2 ab, igg labcorp|
               sars-cov-2 igg antibody titer cutoff, nucleocapsid|
               nucleocapsid igg|
               sars cov 2 igg index|
               sars-cov-2 spike ab interp labcorp|
               sars-cov-2 semi-quant igg ab labcorp|
               sars-cov-2, igg index|
               sars-cov-2 igm level|
               sars-cov-2 igg level'
    )) %>% compute_new(temporary=TRUE,
                      indexes=list('measurement_concept_id',
                                   'person_id')) %>%
    mutate(measurement_type=2000001519L)
  
  combined <- 
    dplyr::union(all_tests,
                 string_search) %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  
}

#' Identifies all COVID-19 antigen tests in the measurement table
#' 
#' @param antigen codeset that contains antigen tests
#' @param measurement_tbl labs from the measurement table
#' 
#' @return measurement table that contains only 
#' antigen tests; 
#' 
#' Function will perform string search on `measurement_source_value`
#' Function will add a column called `measurement_type` to 
#' annotate whether test was derived from structured code or 
#' source value derived
#' 

find_antigen <- function(antigen_codeset,
                          measurement_tbl=cdm_tbl('measurement_labs')) {
  
  
  all_tests <- 
    measurement_tbl %>%
    inner_join(
      select(antigen_codeset,concept_id),
      by=c('measurement_concept_id'='concept_id')
    ) %>% mutate(measurement_type=2000001517L)
  
  string_search <- 
    measurement_tbl %>%
    filter(measurement_concept_id == 0L &
             str_detect(
               lower(measurement_source_value),
               'binaxnow covid-19 result|bd veritor sars-cov2|
               carestart antigen detection covid-19 results'
             )) %>% compute_new(temporary=TRUE,
                                indexes=list('measurement_concept_id',
                                             'person_id')) %>%
    mutate(measurement_type=2000001519L)
  
  combined <- 
    dplyr::union(all_tests,
                 string_search) %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  
}


#' Annotates test results from a measurement_tbl structure
#' 
#' @param measurement_tests Tbl that assumes the measurement_labs
#' structure, with columns `value_as_concept_id` and `value_source_value`
#' 
#' @return Standardizes `value_as_concept_id` by both 
#' choosing just one `concept_id` as well as doing string search
#' in `value_source_value` and mapping to the selected `concept_id`
#' Specifically, 9189 is chosen for the negative tests and 9191 for positives
#' 2000001415 is unknown
#' 

annotate_test_result <- function(measurement_tests) {
  
  results <- 
    measurement_tests %>%
    mutate(
      value_as_concept_id = case_when(
        value_as_concept_id %in% c(9189L,9190L,45878583L,45884153L) ~ 9189L,
        value_as_concept_id %in% c(9191L,4126681L,45884084L,45878745L,4328749L,45876384L,45881666L) ~ 9191L,
        value_as_concept_id == 0L & str_detect(lower(value_source_value),'^positive|^pos|^detected|^present|^presumptive positive') ~ 9191L,
        value_as_concept_id == 0L & str_detect(lower(value_source_value),'negative|neg|^not|^none|undetected|^no.*detected|^absent') ~ 9189L,
        TRUE ~ 2000001415L
      )
    )
  
}


#' Annotates first test 
#' 
#' @param measurement_values Tbl that assumes the measurement_labs
#' structure, with columns `value_as_concept_id` and `measurement_date`
#' 
#' @return returns same input, replacing `value_as_concept_id` with 
#' the `concept_id` for the first test
#' 

annotate_first_test <- function(measurement_values) {
  
  first <- 
    measurement_values %>%
    filter(value_as_concept_id == 9191L) %>%
    group_by(person_id) %>%
    summarise(first_pos_test = min(measurement_date)) %>%
    mutate(first_pos ='yes')
  
  final <- 
    measurement_values %>%
    left_join(first,
              by=c('person_id',
                   'measurement_date'='first_pos_test')) %>%
    mutate(value_as_concept_id=
             case_when(
               value_as_concept_id == 9191L & first_pos == 'yes' ~ 2000001526L,
               TRUE ~ value_as_concept_id
             ))
    
}





