% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_surv.R
\name{fit_survival_forest}
\alias{fit_survival_forest}
\title{Wrapper of \code{grf::survival_forest}}
\usage{
fit_survival_forest(
  formula,
  data,
  id.var,
  time.var,
  event.var,
  nfold = 2,
  time.grid.size = 250,
  oob = TRUE,
  option = list(),
  obs.weight.var,
  ...
)
}
\arguments{
\item{formula}{formula used by \code{\link[grf:survival_forest]{grf::survival_forest}}}

\item{data}{data containing all covariates, follow-up time, event indicator and id}

\item{id.var}{see \code{\link{MRsurv}}}

\item{time.var}{see \code{\link{MRsurv}}}

\item{event.var}{see \code{\link{MRsurv}}}

\item{nfold}{number of folds used when fitting survival curves with sample splitting. Default is 2. If \code{nfold=1}, sample is not split.}

\item{time.grid.size}{size of time grid if more than 1 event times, default to 250}

\item{oob}{whether to use out-of-bag (OOB) fitted values from \code{\link[grf:survival_forest]{grf::survival_forest}} when sample splitting is not used (\code{nfold=1}). Default is \code{TRUE}}

\item{option}{a list containing optional arguments passed to \code{\link[grf:survival_forest]{grf::survival_forest}}. We encourage using a named list. Will be passed to \code{\link[grf:survival_forest]{grf::survival_forest}} by running a command like \code{do.call(survival_forest, option)}. The user should not specify \code{X}, \code{Y}, \code{D}, \code{failure.times} and \code{compute.oob.predictions}.}

\item{obs.weight.var}{see \code{\link{MRsurv}}}

\item{...}{ignored}
}
\value{
a \code{\link{pred_surv}} class containing fitted survival curves for individuals in \code{data}
}
\description{
Wrapper of \code{grf::survival_forest}
}
