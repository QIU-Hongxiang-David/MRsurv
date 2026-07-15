#' @title S3 class for predictive models with an intercept only
#' @name intercept_model
#' @param est the predicted mean
#' @return an "`intercept_model`" object, essentially a list with element `est`.
#' @export
intercept_model<-function(est){
    assert_that(is.number(est))
    out<-list(est=est)
    
    class(out)<-"intercept_model"
    out
}

#' @title S3 class for predictive models with an intercept only and influence function information
#' @name intercept_IF_model
#' @param est the predicted mean
#' @param pseudo.outcome named vector of the pseudo-outcome for observation. The name corresponds to each observation's id. Used in GEE for clustered MR estimator of a scalar estimand (rather than a function).
#' @param IF named vector of influence function evaluated at each observation. The name corresponds to each observation's id. Used to calculate standard error and confidence interval for MR estimator (with iid data) of a scalar estimand (rather than a function).
#' @param SE standard error
#' @return an "`intercept_IF_model`" object, essentially a list with element `est` and `IF`.
#' @export
intercept_IF_model<-function(est,pseudo.outcome,IF,SE){
    assert_that(is.number(est))
    assert_that(is.vector(pseudo.outcome,mode="numeric"))
    assert_that(!is.null(names(pseudo.outcome)))
    if(!is.null(IF)){
        assert_that(is.vector(IF,mode="numeric"))
        assert_that(!is.null(names(IF)))
    }
    out<-list(est=est,pseudo.outcome=pseudo.outcome,IF=IF,SE=SE)
    
    class(out)<-c("intercept_IF_model")
    out
}

#' @export
predict.intercept_model<-function(object,...){
    object$est
}

#' @export
fitted.intercept_IF_model<-predict.intercept_model
#' @export
predict.intercept_IF_model<-predict.intercept_model
#' @export
fitted.intercept_model<-predict.intercept_model


#' @export
print.intercept_IF_model<-function(x,...){
    cat("Point estimate: ",x$est,"\n")
    cat("Standard error: ",x$SE,"\n")
}

#' @export
summary.intercept_IF_model<-function(object,conf.level=.95,...){
    out<-list(est=object$est,SE=object$SE,IF=object$IF,CI=confint(object,conf.level=conf.level),conf.level=conf.level)
    class(out)<-"summary.intercept_IF_model"
    out
}

#' @export
print.summary.intercept_IF_model<-function(x,...){
    cat("Point estimate: ",x$est,"\n")
    cat("Standard error: ",x$SE,"\n")
    cat(round(x$conf.level*100,2),"% confidence interval:\n")
    print(x$CI)
}


#' @export
confint.intercept_IF_model<-function(object,parm,level=.95,...){
    qs<-c((1-level)/2,(1+level)/2)
    CI<-object$est+qnorm(qs)*object$SE
    names(CI)<-paste(round(qs*100,1),"%")
    CI
}


#' @export
print.intercept_model<-function(x,...){
    cat("Point estimate: ",x$est,"\n")
}
