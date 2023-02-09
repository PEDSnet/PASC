

#' Identifies all COVID-19 PCR tests in the measurement table
#' 
#' @param vx_codeset codeset that contains vaccines
#' @param immunization_tbl labs from the measurement table
#' 
#' @return measurement table that contains only 
#' covid19 immunizations; 
#' 
#' Function will perform string search on `immunization_source_value`
#' Function will add a column called `immunization_type` to 
#' annotate whether test was derived from structured code or 
#' source value derived
#' 

find_immunizations <- function(vx_codeset,
                               immunization_tbl=cdm_tbl('immunization')) {
  
  
  all_imm <- 
    immunization_tbl %>%
    inner_join(
      select(vx_codeset,concept_id),
      by=c('immunization_concept_id'='concept_id')
    ) %>% mutate(immunization_type=2000001517L)
  
  string_search <- 
    immunization_tbl %>%
    filter(immunization_concept_id == 0L &
             str_detect(
               lower(immunization_source_value),'sars|covid')
    ) %>% compute_new(temporary=TRUE,
                      indexes=list('immunization_concept_id',
                                   'person_id')) %>%
    mutate(immunization_type=2000001519L)
  
  combined <- 
    dplyr::union(all_imm,
                 string_search) %>%
    mutate(value_as_concept_id = 2000001492L) %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  
}