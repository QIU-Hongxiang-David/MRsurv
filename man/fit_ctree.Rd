% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_surv.R
\name{fit_ctree}
\alias{fit_ctree}
\title{Wrapper of \code{party::ctree}}
\usage{
fit_ctree(
  formula,
  data,
  id.var,
  time.var,
  event.var,
  nfold = 2,
  time.grid.size = 250,
  option = list(),
  obs.weight.var,
  ...
)
}
\arguments{
\item{formula}{formula used by \code{\link[party:ctree]{party::ctree}}}

\item{data}{data containing all covariates, follow-up time, event indicator and id}

\item{id.var}{see \code{\link{MRsurv}}}

\item{time.var}{see \code{\link{MRsurv}}}

\item{event.var}{see \code{\link{MRsurv}}}

\item{nfold}{number of folds used when fitting survival curves with sample splitting. Default is 2. If \code{nfold=1}, sample is not split.}

\item{time.grid.size}{size of time grid if more than 1 event times, default to 250}

\item{option}{a list containing optional arguments passed to \code{\link[party:ctree]{party::ctree}}. We encourage using a named list. Will be passed to \code{\link[party:ctree]{party::ctree}} by running a command like \code{do.call(ctree, option)}. The user should not specify \code{formula} and \code{data}.}

\item{obs.weight.var}{see \code{\link{MRsurv}}}

\item{...}{ignored}
}
\value{
a \code{\link{pred_surv}} class containing fitted survival curves for individuals in \code{data}
}
\description{
Wrapper of \code{party::ctree}
}
