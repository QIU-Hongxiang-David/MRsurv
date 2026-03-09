#' Example data
#' 
#' An Rda image containing an example of data that works with this package.
#' 
#' @format The image contains a data frame `follow.up.time` and a list of data frames `covariates`.
#' `follow.up time` contains four variables:
#' \describe{
#' \item{id}{id for each individual}
#' \item{X}{length of follow-up, i.e., minimum of time-to-event and time-to-censoring}
#' \item{Delta}{indicator of observing event}
#' \item{wt}{sample weight (This variable need not be present)}
#' }
#' The i-th element of `covariates` corresponds to covariates measured at the i-th visit time for individuals still at risk at the time, including the `id` variable for data linkage.
