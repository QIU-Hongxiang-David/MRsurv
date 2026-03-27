#' @title G-computation transformation
#' @name Gtransform
#' @description
#' Given a `pred_surv` object for time to event in a time window, calculates the G-computation transformation in the time window. The transformation is used as the outcome when estimating the conditional survival probability at the next visit time.
#' @param follow.up.time see \code{\link{MRsurv}}
#' @param pred_event_obj a `pred_surv` object for time to event in the time window of interest
#' @param tvals see \code{\link{MRsurv}}.
#' @param next.visit.time the next visit time. Default is `Inf`, corresponding to the last time window
#' @param id.var see \code{\link{MRsurv}}
#' @param time.var see \code{\link{MRsurv}}
#' @param event.var see \code{\link{MRsurv}}
#' @return a named one-column matrix of transformations used for regression. Each row corresponds to an individual. Row names are elements in `follow.up.time$id.var`
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{Gsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @export
Gtransform<-function(follow.up.time,pred_event_obj,tvals,next.visit.time=Inf,id.var,time.var,event.var){
    tvals.bar<-pmin(tvals,next.visit.time) #tvals truncated at next visit time
    
    output<-matrix(nrow=nrow(pred_event_obj$surv),ncol=length(tvals.bar))
    rownames(output)<-rownames(pred_event_obj$surv)
    colnames(output)<-as.character(tvals)
    
    for(j in 1:ncol(output)){
        if(j>1 && tvals.bar[j]==tvals.bar[j-1]){
            output[,j]<-output[,j-1]
        }else{
            #find Shat at t
            #k.t is the index of the last event time in event times that is <= current tvals (t)
            #will use k.t to get Shat at t
            #if no event time < t, then set Shat.t to be 1
            k.t<-find.last.TRUE.index(pred_event_obj$time<=tvals.bar[j],noTRUE=0)
            if(k.t==0){
                Shat<-1
            }else{
                Shat<-pred_event_obj$surv[,k.t]
            }
            output[,j]<-Shat
        }
    }
    output
}



#' @title Regression based on G-computation transformation of fitted survival and censoring probabilities using SuperLearner
#' @name Greg.SuperLearner
#' @description Apply G-computation transformation on fitted survival and censoring probabilities in each time window and estimate P(T > t | T > truncation time, covariates available at truncation time) with \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}.
#' @param covariates see \code{\link{MRsurv}}
#' @param follow.up.time see \code{\link{MRsurv}}
#' @param pred_event.list list of `pred_surv` objects for time to event. Each `pred_event_censor` object in the list corresponds to a time window in `visit.times` after `truncation.index` in increasing order.
#' @param visit.times see \code{\link{MRsurv}}
#' @param tvals see \code{\link{MRsurv}}. Must be sorted in ascending order.
#' @param truncation.index see \code{\link{MRsurv}}
#' @param id.var see \code{\link{MRsurv}}
#' @param time.var see \code{\link{MRsurv}}
#' @param event.var see \code{\link{MRsurv}}
#' @param U.formula see \code{\link{MRsurv}}
#' @param Q.formula see \code{\link{MRsurv}}
#' @param U.SuperLearner.control see \code{\link{MRsurv}}
#' @param Q.SuperLearner.control see \code{\link{MRsurv}}
#' @param U.folds a list of vectors of id (identified by variable `id.var`) corresponding to each fold for cross-fitting. Set to a list containing one vector for no cross-fitting.
#' @param obs.weight.var see \code{\link{MRsurv}}
#' @return a list of `SuperLearner` models (conditional probability) or \code{\link{intercept_IF_model}} objects (marginal probability) corresponding to `tvals`.
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{MRsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @section Custom learners:
#' Custom learners may be specified by providing an element named `SL.library` in `Q.SuperLearner.control`.The user may refer to resources such as \url{https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html} for a guide to create custom learners.
#' @export
Greg.SuperLearner<-function(
    covariates,
    follow.up.time,
    pred_event.list,
    visit.times,
    tvals,
    truncation.index,
    id.var,
    time.var,
    event.var,
    U.formula=NULL,
    Q.formula=~.,
    U.SuperLearner.control=list(family=gaussian(),SL.library="SL.lm"),
    Q.SuperLearner.control=U.SuperLearner.control,
    U.folds,
    obs.weight.var=NULL
){
    index.shift<-truncation.index-1 #shift for the index of pred_event.list
    
    models<-lapply(seq_along(tvals),function(i){
        #K is the last visit.time that needs to be considered
        K<-find.last.TRUE.index(visit.times<tvals[i])
        
        for(k in K:truncation.index){
            if(k<length(visit.times)){
                pred_event_obj<-truncate_pred_surv(pred_event.list[[k-index.shift]],visit.times[k+1])
            }else{
                pred_event_obj<-pred_event.list[[k-index.shift]]
            }
            if(k==K){
                next.visit.time<-tvals[i]
            }else{
                next.visit.time<-visit.times[k+1]
            }
            
            history<-reduce(covariates[1:ifelse(k==truncation.index,k,k-1)],.f=function(d1,d2){
                right_join(d1,d2,by=id.var)
            })%>%arrange(.data[[id.var]])
            
            Y.G<-Gtransform(follow.up.time,pred_event_obj,tvals[i],
                            next.visit.time=next.visit.time,
                            id.var,time.var,event.var)
            # Y.G<-Y.G[order(rownames(Y.G)),1]
            Y.G<-sort_by(Y.G[,1],rownames(Y.G))
            if(k==K){
                Y<-Y.G
            }else{
                if(length(U.folds)==1){
                    U<-as.numeric(predict(model,newdata=history%>%select(!.data[[id.var]]),onlySL=TRUE)$pred)
                    names(U)<-history%>%pull(.data[[id.var]])
                }else{
                    U.list<-lapply(1:length(U.folds),function(v){
                        newdata<-history%>%filter(.data[[id.var]] %in% U.folds[[v]])
                        U<-as.numeric(predict(models[[v]],newdata=newdata%>%select(!.data[[id.var]]),onlySL=TRUE)$pred)
                        names(U)<-newdata%>%pull(.data[[id.var]])
                        U
                    })
                    U<-do.call(c,U.list)
                    # U<-U[order(names(U))]
                    U<-sort_by(U,names(U))
                }
                
                Y<-U*Y.G
                names(Y)<-names(U)
            }
            
            
            if(k>truncation.index){
                form<-U.formula[[k-truncation.index]]
            }else{
                form<-Q.formula
            }
            train.data<-history%>%filter(.data[[id.var]] %in% names(Y))%>%arrange(.data[[id.var]])
            X<-model.frame(form,train.data%>%select(!.data[[id.var]]))
            
            if(is.null(obs.weight.var)){
                obsWeights<-NULL
            }else{
                obsWeights<-follow.up.time%>%filter(.data[[id.var]] %in% names(Y))%>%arrange(.data[[id.var]])%>%pull(obs.weight.var)
                names(obsWeights)<-names(Y)
            }
            
            if(k==truncation.index && ncol(X)==0){
                if(is.null(obsWeights)){
                    est<-mean(Y)
                }else{
                    # message(paste(obs.weight.var,"might not be correctly accounted for in the standard error due to failure fully account for the sampling scheme."))
                    est<-mean(Y*obsWeights)/mean(obsWeights)
                    # est<-mean(Y*obsWeights)
                }
                model<-intercept_model(est)
                return(model)
            }else{
                if(k==truncation.index){
                    SuperLearner.arg<-c(
                        list(Y=Y,X=X,obsWeights=obsWeights),
                        Q.SuperLearner.control
                    )
                    model<-do.call(SuperLearner,SuperLearner.arg)
                    return(model)
                }else{
                    if(length(U.folds)==1){
                        SuperLearner.arg<-c(
                            list(Y=Y,X=X,obsWeights=obsWeights),
                            U.SuperLearner.control
                        )
                        model<-do.call(SuperLearner,SuperLearner.arg)
                    }else{
                        models<-lapply(U.folds,function(fold){
                            X<-model.frame(form,train.data%>%filter(!(.data[[id.var]] %in% fold))%>%select(!.data[[id.var]]))
                            SuperLearner.arg<-c(
                                list(Y=Y[!(names(Y) %in% fold)],X=X,obsWeights=obsWeights[!(names(obsWeights) %in% fold)]),
                                U.SuperLearner.control
                            )
                            do.call(SuperLearner,SuperLearner.arg)
                        })
                    }
                }
            }
        }
    })
    names(models)<-as.character(tvals)
    models
}
