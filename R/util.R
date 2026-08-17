#copy of .SL.require
.requireNamespace<-function (package, message = paste("loading required package (", package, ") failed", sep = "")) 
{
    if (!requireNamespace(package, quietly = FALSE)) {
        stop(message, call. = FALSE)
    }
    invisible(TRUE)
}

.require<-function(package, message = paste("loading and attaching required package (",package, ") failed", sep = "")){
    if(!(package %in% .packages())){
        message(paste("Loading and attaching required package (",package,") to the search list",sep = ""))
        if (!suppressWarnings(require(package, quietly = TRUE,character.only=TRUE))) {
            stop(message, call. = FALSE)
        }
        invisible(TRUE)
    }
}

#find the first index of x that is TRUE
#noTRUE: return value when none of x is TRUE
find.first.TRUE.index<-function(x,noTRUE=length(x)+1){
    match(TRUE,x,nomatch=noTRUE)
}

#find the last index of x that is TRUE (loop from the last entry back to the first to and find the first encountered TRUE)
#noTRUE: return value when none of x is TRUE
find.last.TRUE.index<-function(x,noTRUE=0){
    length(x)+1-match(TRUE,rev(x),nomatch=length(x)+1-noTRUE)
}

#' @title Clip all elements in an array to fall in an interval
#' @name clip_interval
#' @param x numeric array or matrix
#' @param lower lower end of the interval
#' @param upper upper end of the interval
#' @return clipped array
#' @details For each element `X` in `x`, `clip_interval(X)` equals `lower` if `X`<`lower`, `upper` if `X`>`upper`, and `X` otherwise. This function may be useful to clip predictions to fall in, e.g., \eqn{[0,1]}. \code{lower} and \code{upper} may be vectors for entry-wise clipping.
#' @export
clip_interval<-function(x,lower=-Inf,upper=Inf){
    assert_that(is.numeric(x))
    if(!is.vector(x) && !is.matrix(x)){
        stop("x should be a vector or matrix")
    }
    
    out<-pmax(pmin(x,upper),lower)
    if(is.matrix(x)){
        out<-matrix(out,nrow=nrow(x),ncol=ncol(x))
    }
    out
}


#' @title Apply administrative censoring to follow up times at a fixed censor.time
#' @name admin.censor
#' @description Right-censor all observations in `follow.up.time` at `censor.time`. This is useful when fitting survival curves within each time window
#' @param follow.up.time See \code{\link{MRsurv}}
#' @param time.var See \code{\link{MRsurv}}
#' @param event.var See \code{\link{MRsurv}}
#' @param censor.time Administrative right-censoring time. Default is `Inf`, i.e., no censoring
#' @return A data frame with the same shape as `follow.up.time` with times and event indicators modified to reflect right-censoring at `censor.time`.
#' @export
admin.censor<-function(follow.up.time,time.var,event.var,censor.time=Inf){
    if(censor.time==Inf){
        follow.up.time
    }else{
        viewed.censored<-pull(follow.up.time,all_of(time.var))>censor.time
        follow.up.time%>%
            mutate("{time.var}":=ifelse(.env$viewed.censored,.env$censor.time,.data[[time.var]]),
                   "{event.var}":=ifelse(.env$viewed.censored,0,.data[[event.var]]))
    }
}


# #create k folds of a vector id with roughly balanced events
# create.folds<-function(id,Delta,k){
#     id0<-id[Delta==0]
#     id1<-id[Delta==1]
#     
#     order0<-sample.int(length(id0))
#     order1<-sample.int(length(id1))
#     
#     d0<-suppressWarnings(data.frame(cbind(id0[order0],1:k)))[1:length(id0),,drop=FALSE]
#     d1<-suppressWarnings(data.frame(cbind(id1[order1],1:k)))[1:length(id1),,drop=FALSE]
#     d<-rbind(d0,d1)
#     names(d)<-c("id","fold.id")
#     lapply(tapply(d$id,d$fold.id,identity,simplify=FALSE),sort)
# }

#convert a vector to a row matrix and return the input if it is already a matrix
#essentially a copy of as.matrix.default
# as_matrix_rowvec<-function(x){
#     if(is.matrix(x)){
#         x
#     }else{
#         array(x,c(1L,length(x)),
#               if(!is.null(names(x)))
#                   list(names(x), NULL)
#               else
#                   NULL
#         )
#     }
# }


#shift censoring time to the left a tiny bit
#used for survival analysis with event observed iff T<C (rather than traditionally T<=C)
#status: 1 if observed event; 0 if censoring
left.shift.censoring<-function(time,status){
    epsilon<-c(0,time)%>%unique%>%sort%>%diff%>%min
    epsilon<-min(epsilon*.5,1e-5)
    time[status==0]<-time[status==0]-epsilon
    time
}



#' @title Options of SuperLearner for fitting conditional survival probability functions Q and U
#' @description
#' The output is essentially a list to be passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} by running a command like `do.call(SuperLearner, QU.SuperLearner.control)`. The user should not specify `Y` and `X`.
#' @name QU.SuperLearner.control
#' @param family family passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}. Defaults to `gaussian()`. For compatibility with machine learning algorithms in \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, `family$family` must be "gaussian", with a possibly non-identity link (e.g., `gaussian(link=log)`) for more flexible usage of \code{\link[SuperLearner:SL.glm]{SuperLearner::SL.glm}}.
#' @param SL.library library of learners passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}. Defaults to linear regression "SL.lm"
#' @param ... further named arguments passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}
#' @return a `QU.SuperLearner.control` object
#' @section Custom learners:
#' Custom learners may be specified by providing an character vector to `SL.library`. The user may refer to resources such as \url{https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html} for a guide to create custom learners.
#' @export
QU.SuperLearner.control<-function(family=gaussian(),SL.library="SL.lm",...){
    out<-list(family=family,SL.library=SL.library,...)
    if(any(c("Y","X","obsWeights") %in% names(out))){
        stop("Y, X and obsWeights should not be specified")
    }
    class(out)<-"QU.SuperLearner.control"
    out
}
