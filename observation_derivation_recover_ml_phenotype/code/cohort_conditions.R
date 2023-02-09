


#' Identifies all PASC diagnoses in the condition occurrence table
#'
#' @param pasc_codeset codeset that contains pasc diagnoses
#' @param cond_tbl condition occurrence table
#'
#' @return condition occurrence table that contains only
#' pasc diagnoses
#'
#' Function will perform string search on `condition_source_value`
#' Function will add a column called `condition_type` to
#' annotate whether test was derived from structured code or
#' source value derived
#' Function will add a column called `condition_qualifier_concept_id`
#' column which is set to condition_concept_id if derived from SNOMED
#' code or source value and set to condition_source_concept_id if derived
#' from ICD code. (currently, all PASC diagnosis codes are in the ICD vocabulary)
#'

compute_pasc_dx <- function(pasc_codeset,
                      cond_tbl=cdm_tbl('condition_occurrence') %>% filter(! condition_type_concept_id %in% c(2000001423L, 2000001424L))) {

  pasc_conditions<-
    cond_tbl %>%
    inner_join(select(pasc_codeset, concept_id), by=c("condition_source_concept_id"="concept_id")) %>%
    mutate(condition_type=2000001517L) %>%
    mutate(condition_qualifier_concept_id=condition_source_concept_id) %>%
    compute_new(temporary=TRUE)

  string_search <- cond_tbl %>%
    filter(grepl("U09.9", condition_source_value, ignore.case=TRUE)|
             (condition_source_concept_id %in% c(45605252L,35206022L) & grepl("covid", condition_source_value, ignore.case=TRUE))|
             (
               (
                 grepl("post", condition_source_value, ignore.case=TRUE) &
                   grepl("acute", condition_source_value, ignore.case=TRUE) &
                   (grepl("covid", condition_source_value, ignore.case=TRUE)|grepl("sars", condition_source_value, ignore.case=TRUE))
               )|
                 (
                   grepl("late", condition_source_value, ignore.case=TRUE) &
                     grepl("effect", condition_source_value, ignore.case=TRUE) &
                     (grepl("covid", condition_source_value, ignore.case=TRUE)|grepl("sars", condition_source_value, ignore.case=TRUE))
                 )|
                 (
                   grepl("long", condition_source_value, ignore.case=TRUE) &
                     (grepl("covid", condition_source_value, ignore.case=TRUE)|grepl("sars", condition_source_value, ignore.case=TRUE))
                 )|
                 (
                   grepl("after", condition_source_value, ignore.case=TRUE) &
                     (grepl("covid", condition_source_value, ignore.case=TRUE)|grepl("sars", condition_source_value, ignore.case=TRUE))
                 )&
                 !grepl("U09.9", condition_source_value, ignore.case=TRUE)&
                 !grepl("vaccination", condition_source_value, ignore.case=TRUE)
             )
    ) %>% mutate(condition_type=2000001519L) %>%
    mutate(condition_qualifier_concept_id=condition_concept_id) %>%
    compute_new(temporary=TRUE)

  combined <-
    dplyr::union(pasc_conditions %>% mutate(temp_order=0),
                 string_search %>% mutate(temp_order=1)) %>%
    group_by(condition_occurrence_id) %>%
    slice_min(temp_order) %>%
    ungroup %>%
    select(-temp_order) %>%
    mutate(value_as_concept_id=2000001520L) %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
}

#' Identifies all MIS-C diagnoses in the condition occurrence table
#'
#' @param misc_codeset codeset that contains misc diagnoses
#' @param cond_tbl condition occurrence table
#'
#' @return condition occurrence table that contains only
#' misc diagnoses
#'
#' Function will perform string search on `condition_source_value`
#' Function will add a column called `condition_type` to
#' annotate whether test was derived from structured code or
#' source value derived
#' Function will add a column called `condition_qualifier_concept_id`
#' column which is set to condition_concept_id if derived from SNOMED
#' code or source value and set to condition_source_concept_id if derived
#' from ICD code. (currently, all MIS-C diagnosis codes are in the ICD vocabulary)
#'
compute_misc_dx <- function(misc_codeset,
                         cond_tbl=cdm_tbl('condition_occurrence')%>% filter(!condition_type_concept_id %in% c(2000001423L, 2000001424L))) {

  misc_codeset_collect <- misc_codeset %>% collect

  misc_codeset_vector <- as.vector(misc_codeset_collect$concept_id)

  misc_conditions_coded <-
    cond_tbl %>%
    mutate(condition_qualifier_concept_id=
             case_when(condition_concept_id %in% c(misc_codeset_vector) ~ condition_concept_id,
                       condition_source_concept_id %in% c(misc_codeset_vector) ~ condition_source_concept_id,
                       TRUE ~ 0L)) %>%
    filter(! condition_qualifier_concept_id == 0L) %>%
    mutate(condition_type=2000001517L) %>%
    compute_new()

  string_search <-
    cond_tbl %>%
    filter((str_detect(lower(condition_source_value),'mis-c')) |
             (str_detect(lower(condition_source_value),'multisystem') &
                str_detect(lower(condition_source_value), 'covid'))) %>%
    #filter(#condition_concept_id == 0L &
            # str_detect(lower(condition_source_value), 'misc-c|multisystem|inflam')) %>%
    #filter((grepl("mis-c", condition_source_value, ignore.case=TRUE))|
           #  (grepl("multisystem", condition_source_value, ignore.case=TRUE) & grepl("inflam", condition_source_value, ignore.case=TRUE)))%>%
    mutate(condition_type=2000001519L) %>%
    mutate(condition_qualifier_concept_id=condition_concept_id) %>%
    compute_new()

  combined <-# misc_conditions_snomed %>% mutate(temp_order=0) %>%
    dplyr::union(misc_conditions_coded %>% mutate(temp_order=0),
                 string_search %>% mutate(temp_order=1)) %>%
   # dplyr::union(misc_conditions_icd %>% mutate(temp_order=1)) %>%
   # dplyr::union(string_search %>% mutate(temp_order=2)) %>%
    group_by(condition_occurrence_id) %>%
    slice_min(temp_order) %>%
    ungroup %>%
    select(-temp_order) %>%
    mutate(value_as_concept_id=703578L) #%>%
#    compute_new(temporary=TRUE,
#                indexes=list('person_id'))
}




#' Identifies all B94.8 diagnoses in the condition occurrence table
#'
#' @param b94_8_codeset codeset that contains B94.8 diagnoses
#' @param cond_tbl condition occurrence table
#'
#' @return condition occurrence table that contains only
#' B94.8 diagnoses
#'
#'
compute_b94_8_dx <- function(b94_8_codeset, b94_8_exclude_codeset,
                            cond_tbl=cdm_tbl('condition_occurrence') %>% filter(!condition_type_concept_id %in% c(2000001423L, 2000001424L))) {

  b94_8_codeset_collect <- b94_8_codeset %>% collect

  exclude_codeset_collect <- b94_8_exclude_codeset %>% collect


  b94_8_codeset_vector <- as.vector(b94_8_codeset_collect$concept_id)

  exclude_codeset_vector <- as.vector(exclude_codeset_collect$concept_id)

  exclude_cond_tbl<-cond_tbl %>%
    filter((condition_concept_id %in% exclude_codeset_vector)|(condition_source_concept_id %in% exclude_codeset_vector)) %>%
    compute_new(index="condition_occurrence_id")


  b94_8_conditions_coded <-
    cond_tbl %>%
    mutate(condition_qualifier_concept_id=
             case_when(condition_concept_id %in% c(b94_8_codeset_vector) ~ condition_concept_id,
                       condition_source_concept_id %in% c(b94_8_codeset_vector) ~ condition_source_concept_id,
                       TRUE ~ 0L)) %>%
    filter(! condition_qualifier_concept_id == 0L) %>%
    mutate(condition_type=2000001517L) %>%
    anti_join(exclude_cond_tbl, by="condition_occurrence_id") %>%
    compute_new()


  ###### FILL IN ID
  final <-b94_8_conditions_coded %>%
    mutate(value_as_concept_id=2000001533L) %>%
    compute_new(temporary=TRUE, indexes=list('person_id'))
}

#' Identifies all conditions by codeset and source values that are mapped to concept_id's
#'
#' @param covid_dx_codeset codeset with snomed concepts of interest
#' @param covid_src_codeset codeset with non-specific snomed concepts that have source values that indicate covid
#' @param condition_tbl the condition tbl from OMOP cdm
#'
#' @return A condition table with metadata appropriate to fill in the observation table.
#' The `value_as_concept_id` field will classify whether the code is an exposure/dx, complication, or history code
#' The `qualifier_concept_id` is the `concept_id`, whether mapped or not mapped
#' The `observation_type_concept_id` will categorize whether a structured code was used or if the codes were source
#' value derived
#'

compute_covid_dx <- function(covid_dx_codeset,
                             covid_src_codeset,
                             condition_tbl=cdm_tbl('condition_occurrence') %>% filter(!condition_type_concept_id %in% c(2000001423L, 2000001424L)),
                             covid_dx_broad=c(37311061L)) {

  conds <-
    condition_tbl %>%
    inner_join(
      select(covid_dx_codeset,
             concept_id,
             dx_class),
             by=c('condition_concept_id'='concept_id')
      ) %>%
    filter(! condition_concept_id %in% covid_dx_broad) %>%
    filter(! str_detect(lower(condition_source_value), 'multisystem|post|vaccin|mis-c|hauler|long|kawasaki|preprocedure|preoperative')) %>%
    mutate(dx_class=case_when(str_detect(lower(condition_source_value),'history of covid|history of 2019') ~ 'hx',
                              str_detect(lower(condition_source_value),'covid toes') ~ 'comp',
                              TRUE ~ dx_class)) %>%
    mutate(condition_qualifier_concept_id = condition_concept_id) %>%
    mutate(condition_type = 2000001517L)

  conds_broad <- condition_tbl %>%
    filter(condition_concept_id %in% covid_dx_broad) %>%
    mutate(dx_class = 'dx') %>%
    filter(! str_detect(lower(condition_source_value), 'multisystem|post|vaccin|hauler|mis-c|long|kawasaki')) %>%
    mutate(condition_qualifier_concept_id = condition_concept_id) %>%
    mutate(condition_type = 2000001517L)

  conds_src <-
    condition_tbl %>%
    inner_join(
      select(covid_src_codeset,
             concept_id,
             dx_class),
             by=c('condition_concept_id'='concept_id')
    ) %>% filter(str_detect(lower(condition_source_value),'covid|sars-cov-2')) %>%
    filter(! str_detect(lower(condition_source_value), 'multisystem|post|vaccin|mis-c|hauler|long|kawasaki|ab test negative')) %>%
    mutate(dx_class=case_when(str_detect(lower(condition_source_value), 'exposure') ~ 'exp',
                              TRUE ~ dx_class)) %>%
    mutate(condition_qualifier_concept_id = condition_concept_id) %>%
    mutate(condition_type = 2000001519L)

  dplyr::union(conds,
               conds_broad) %>%
    dplyr::union(conds_src) %>%
    mutate(value_as_concept_id =
             case_when(dx_class == 'dx' ~ 2000001525L,
                       dx_class == 'hx' ~ 2000001522L,
                       dx_class == 'comp' ~ 2000001523L,
                       dx_class == 'exp' ~ 2000001524L)) %>%
    compute_new(temporary=TRUE)

}

