#' Simulated dataset simdataHADS
#' 
#' The data mimic the PREDIALA study described and analyzed in Proust-Lima et al (2021 - 
#' https://arxiv.org/abs/2109.13064). 
#' The study aims to describe the trajectories of depressive symptomatology of patients 
#' suffering end-stage renal disease and registered on the renal transplant waiting list.
#' Repeated measures of anxiety and depression (HADS) were simulated at different times
#' of measurement for 561 subjects. Four time-independent covariates were also generated: 
#' group (dialyzed or pre-emptive), sex and age at entry in the cohort and time on the 
#' waiting list at entry in the cohort. 
#' 
#' 
#' @name simdataHADS
#' @docType data
#' @format A data frame with 1140 observations on the following 13 variables.
#' \describe{
#' \item{grp}{group with 0=dialyzed and 1=preemptive}
#' \item{sex}{sex with 0=woman and 1=man}
#' \item{age}{age at entry in the cohort}
#' \item{hads_2}{item 2 of HADS measuring depression}
#' \item{hads_4}{item 4 of HADS measuring depression}
#' \item{hads_6}{item 6 of HADS measuring depression}
#' \item{hads_8}{item 8 of HADS measuring depression}
#' \item{hads_10}{item 10 of HADS measuring depression}
#' \item{hads_12}{item 12 of HADS measuring depression}
#' \item{hads_14}{item 14 of HADS measuring depression}
#' \item{ID}{subject identification number}
#' \item{time}{time of measurement}
#' \item{time_entry}{time on the waiting list at entry in the cohort}
#' }
#' @keywords datasets
NULL
