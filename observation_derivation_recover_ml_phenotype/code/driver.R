# Top-level code for execution of data request

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))

# Need to do this for assignInNamespace to work
suppressPackageStartupMessages(library(dbplyr))

# Required for execution using Rscript
suppressPackageStartupMessages(library(methods))

#' Set up the execution environment
#'
#' The .load() function sources the R files needed to execute the query
#' and sets up the execution environment.  In particular, all of the base
#' framework files, as well as files inthe code_dir with names matching
#' `cohort_*.R` or `analyze_*.R` will be sourced.
#'
#' This function is usually run automatically when the `run.R` file is sourced
#' to execute the request.  It may also be executed manually during an
#' interactive session to re-source changed code or to re-establish a connection
#' to the database.
#'
#' **N.B.** You will almost never have to edit this function.
#'
#' @param here The name of the top-level directory for the request.  The default
#'   is `config('base_dir')` if the config function has been set up, or the
#'   global variable `base_dir` if not.
#'
#' @return The value of `here`.
#' @md
.load <- function(here = ifelse(typeof(get('config')) == 'closure',
                                config('base_dir'), base_dir)) {
    source(file.path(here, 'code', 'config.R'))
    source(file.path(here, 'code', 'req_info.R'))
    source(config('site_info'))
    source(file.path(here, config('subdirs')$code_dir, 'setup.R'))
    source(file.path(here, config('subdirs')$code_dir, 'codesets.R'))
    for (fn in list.files(file.path(here, config('subdirs')$code_dir),
                        'util_.+\\.R', full.names = TRUE))
      source(fn)
    for (fn in list.files(file.path(here, config('subdirs')$code_dir),
                          'cohort_.+\\.R', full.names = TRUE))
      source(fn)
    for (fn in list.files(file.path(here, config('subdirs')$code_dir),
                        'analyze_.+\\.R', full.names = TRUE))
      source(fn)
    source(file.path(here, config('subdirs')$code_dir, 'cohorts.R'))

    .env_setup()

    for (def in c('retain_intermediates', 'results_schema')) {
      if (is.na(config(def)))
        config(def, config(paste0('default_', def)))
    }

    here
}

#' Execute the request
#'
#' This function presumes the environment has been set up, and executes the
#' steps of the request.
#'
#' In addition to performing queries and analyses, the execution path in this
#' function should include periodic progress messages to the user, and logging
#' of intermediate totals and timing data through [append_sum()].
#'
#' This function is also typically executed automatically, but is separated from
#' the setup done in [.load()] to facilitate direct invocation during
#' development and debugging.
#'
#' @param base_dir The name of the top-level directory for the request.  The default
#'   is `config('base_dir')`, which should always be valid after execution of
#'   [.load()].
#'
#' @return The return value is dependent on the content of the request, but is
#'   typically a structure pointing to some or all of the retrieved data or
#'   analysis results.  The value is not used by the framework itself.
#' @md
.run  <- function(base_dir = config('base_dir')) {

  message('Starting execution with framework version ',
          config('framework_version'))

  init_sum(cohort = 'Start', persons = 0)
  rslt <- list()

  message('Computing PCR tests')
    rslt$pcr_raw <- find_pcr(pcr_codeset = load_codeset('covid_pcr')) 
    rslt$pcr_values <- annotate_test_result(rslt$pcr_raw)
    rslt$pcr_final <- annotate_first_test(rslt$pcr_values)

    append_sum(cohort = 'PCR tests',
               persons = distinct_ct(rslt$pcr_final))
  
  message('Computing antigen tests')
    rslt$antigen_raw <- find_antigen(antigen_codeset = load_codeset('covid_antigen'))
    rslt$antigen_values <- annotate_test_result(rslt$antigen_raw)
    rslt$antigen_final <- annotate_first_test(rslt$antigen_values)
    
    append_sum(cohort = 'Antigen tests',
               persons = distinct_ct(rslt$antigen_final))
    
  message('Computing serology tests')
    rslt$serology_raw <- find_serology(serology_codeset = load_codeset('covid_serology'))
    rslt$serology_values <- annotate_test_result(rslt$serology_raw)
    rslt$serology_final <- annotate_first_test(rslt$serology_values)
    
    append_sum(cohort = 'Serology tests',
              persons = distinct_ct(rslt$serology_final))
    
  message('Computing immunizations')
    rslt$immunization <- find_immunizations(vx_codeset = load_codeset('covid_vx'))

    append_sum(cohort = 'Immunizations',
               persons = distinct_ct(rslt$immunization))
    
  message('Computing long covid dx')
    rslt$pasc_dx <- compute_pasc_dx(pasc_codeset = load_codeset('pasc_dx'))
    
    append_sum(cohort = 'Long Covid dx',
               persons = distinct_ct(rslt$pasc_dx))
    
  message('Computing C19 conditions')
    rslt$c19_dx <- compute_covid_dx(covid_dx_codeset=load_codeset('covid_dx','icccc'),
                                    covid_src_codeset=load_codeset('covid_dx_src_val')) %>%
      anti_join(rslt$pasc_dx, by="condition_occurrence_id")
    
    append_sum(cohort = 'C19 dx',
               persons = distinct_ct(rslt$c19_dx))
    

    

    
  message('Computing MIS-C')
    rslt$misc_dx <- compute_misc_dx(misc_codeset = load_codeset('misc_dx'))
    
    append_sum(cohort = 'MIS-C dx',
               persons = distinct_ct(rslt$misc_dx))
    
  message('Computing B94.8')
    rslt$b94_8_dx <- (compute_b94_8_dx(b94_8_codeset = load_codeset('b94_8_dx'), b94_8_exclude_codeset = load_codeset('b94_8_exclude')) %>%
      anti_join(rslt$pasc_dx, by="condition_occurrence_id"))%>%
      anti_join(rslt$c19_dx, by="condition_occurrence_id") %>% 
      compute_new(temporary=TRUE,indexes=list('person_id'))
    
    append_sum(cohort = 'B94.8 dx',
               persons = distinct_ct(rslt$b94_8_dx))
    
  message('Filling Observation Table: PCR')
  rslt$pcr_obs <- create_observation_records_meas(src_tbl=rslt$pcr_final,
                                                  obs_concept_id=2000001530L)
  output_tbl(rslt$pcr_obs, 'observation_derivation_recover_pcr')
  
  message('Filling Observation Table: antigen')
  rslt$antigen_obs <- create_observation_records_meas(src_tbl=rslt$antigen_final,
                                                      obs_concept_id=2000001529L)
  output_tbl(rslt$antigen_obs, 'observation_derivation_recover_antigen')
  
  message('Filling Observation Table: Serology')
  rslt$serology_obs <- create_observation_records_meas(src_tbl=rslt$serology_final,
                                                       obs_concept_id=2000001528L)
  output_tbl(rslt$serology_obs, 'observation_derivation_recover_serology')
  
  message('Filling Observation Table: Immunization')
  rslt$immunizations_obs <- create_observation_records_imm(src_tbl=rslt$immunization,
                                                       obs_concept_id=2000001492L)
  output_tbl(rslt$immunizations_obs, 'observation_derivation_recover_imm')
  
  message('Filling Observation Table: Covid19 DX')
  rslt$c19_dx_obs <- create_observation_records_cond(src_tbl=rslt$c19_dx,
                                                     obs_concept_id=2000001527L)
  output_tbl(rslt$c19_dx_obs, 'observation_derivation_recover_c19')
  
  message('Filling Observation Table: PASC')
  rslt$pasc_dx_obs <- create_observation_records_cond(src_tbl=rslt$pasc_dx,
                                                  obs_concept_id=2000001527L)
  output_tbl(rslt$pasc_dx_obs, 'observation_derivation_recover_pasc')
  
  message('Filling Observation Table: MISC')
  rslt$misc_dx_obs <- create_observation_records_cond(src_tbl=rslt$misc_dx,
                                                      obs_concept_id=2000001527L)
  
  message('Filling Observation Table: B94.8')
  rslt$b94_8_dx_obs <- create_observation_records_cond(src_tbl=rslt$b94_8_dx,
                                                      obs_concept_id=2000001527L)
  
  output_tbl(rslt$b94_8_dx_obs, 'observation_derivation_recover_b94_8')
  
  
  rslt_final <- rslt[c(15:22)]
  
  rslt$observation_derivation_recover <- 
    reduce(.x=rslt_final,
           .f=dplyr::union)
  
  output_tbl(rslt$observation_derivation_recover,
             'observation_derivation_recover')
  

  message('Done.')

  #invisible(rslt)

}

#' Set up and execute a data request
#'
#' This function encapsulates a "production" run of the data request.  It sets
#' up the environment, executes the request, and cleans up the environment.
#'
#' Typically, the `run.R` file calls run_request() when in a production mode.
#'
#' @param base_dir Path to the top of the data request files.  This is
#'   typically specified in `run.R`.
#'
#' @return The result of [.run()].
#' @md
run_request <- function(base_dir) {
    base_dir <- .load(base_dir)
    on.exit(.env_cleanup())
    .run(base_dir)
}
