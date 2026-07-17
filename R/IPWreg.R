#' @title Inverse-probability weighting (IPCW) transformation
#' @name IPCWtransform
#' @description
#' Given a `pred_surv` object for time to censoring in a time window, calculates the IPCW transformation in the time window. The transformation is used as the outcome when estimating the conditional survival probability at the next visit time.
#' @param follow.up.time see \code{\link{MRsurv}}
#' @param pred_censor_obj a `pred_surv` object for time to censoring in the time window of interest
#' @param tvals see \code{\link{MRsurv}}.
#' @param next.visit.time the next visit time. Default is `Inf`, corresponding to the last time window
#' @param id.var see \code{\link{MRsurv}}
#' @param time.var see \code{\link{MRsurv}}
#' @param event.var see \code{\link{MRsurv}}
#' @param denom.survival.trunc see \code{\link{MRsurv}}
#' @return a named one-column matrix of transformations used for regression. Each row corresponds to an individual. Row names are elements in `follow.up.time[[id.var]]`
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{IPCWsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @export
IPCWtransform<-function(follow.up.time,pred_censor_obj,tvals,next.visit.time=Inf,id.var,time.var,event.var,denom.survival.trunc){
    tvals.bar<-pmin(tvals,next.visit.time) #tvals truncated at next visit time
    
    output<-matrix(nrow=nrow(pred_censor_obj$surv),ncol=length(tvals.bar))
    rownames(output)<-rownames(pred_censor_obj$surv)
    colnames(output)<-as.character(tvals)
    
    for(i in 1:nrow(output)){
        id.matching.data<-follow.up.time%>%filter(.data[[id.var]]==rownames(output)[i])
        X<-pull(id.matching.data,all_of(time.var))
        Delta<-pull(id.matching.data,all_of(event.var))
        for(j in 1:ncol(output)){
            if(j>1 && tvals.bar[j]==tvals.bar[j-1]){
                output[i,j]<-output[i,j-1]
            }else{
                if(Delta==1 && X<=tvals.bar[j]){
                    #find Ghat(X-)
                    #k.GX is the index of the first censoring time in censoring times that is >= X
                    #will use k.GX-1 to get Ghat(X-)
                    #if no censoring time >= X, then set k.GX to be the index after the last censoring time
                    #if k.GX=1 (all censoring time >= X), set Ghat(X-) to 1
                    k.GX<-find.first.TRUE.index(pred_censor_obj$time>=X,
                                                noTRUE=length(pred_censor_obj$time)+1)
                    if(k.GX==1){
                        Ghat.Xminus<-1
                    }else{
                        Ghat.Xminus<-pred_censor_obj$surv[i,k.GX-1]
                    }
                    
                    output[i,j]<-1-1/pmax(Ghat.Xminus,denom.survival.trunc)
                }else{
                    output[i,j]<-1
                }
            }
        }
    }
    
    output
}


#' @title Regression based on IPCW transformation of fitted survival and censoring probabilities using SuperLearner
#' @name IPCWreg.SuperLearner
#' @description Apply IPCW transformation on fitted survival and censoring probabilities in each time window and estimate P(T > t | T > truncation time, covariates available at truncation time) with \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}.
#' @param covariates see \code{\link{MRsurv}}
#' @param follow.up.time see \code{\link{MRsurv}}
#' @param pred_censor.list list of `pred_surv` objects for time to censoring. Each `pred_censor.list` object in the list corresponds to a time window in `visit.times` after `truncation.index` in increasing order.
#' @param visit.times see \code{\link{MRsurv}}
#' @param tvals see \code{\link{MRsurv}}.
#' @param truncation.index see \code{\link{MRsurv}}
#' @param id.var see \code{\link{MRsurv}}
#' @param time.var see \code{\link{MRsurv}}
#' @param event.var see \code{\link{MRsurv}}
#' @param Q.formula formula to specify covariates being used for estimating P(T > t | T > `visit.times[truncation.index]`, covariates available at `visit.times[truncation.index]`). Set to include intercept only (`~ 0` or `~ -1`) for marginal survival probability, which is simply the mean of pseudo-outcomes. Default is `~ .`, which includes main effects of all available covariates up to (inclusive) the `truncation.time`.
#' @param Q.SuperLearner.control see \code{\link{MRsurv}}
#' @param cluster.var see \code{\link{MRsurv}}. If provided, this variable is used to split samples when using cross-fitting and/or cross-valiation, as well as inference of marginal survival probability.
#' @param obs.weight.var see \code{\link{MRsurv}}
#' @param corstr see \code{\link{MRsurv}}
#' @param denom.survival.trunc see \code{\link{MRsurv}}
#' @return a list of `SuperLearner` models (conditional probability) or \code{\link{intercept_IF_model}} objects (marginal probability) corresponding to `tvals`.
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{MRsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @section Custom learners:
#' Custom learners may be specified by providing an element named `SL.library` in `Q.SuperLearner.control`.The user may refer to resources such as \url{https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html} for a guide to create custom learners.
#' @export
IPCWreg.SuperLearner<-function(
    covariates,
    follow.up.time,
    pred_censor.list,
    visit.times,
    tvals,
    truncation.index,
    id.var,
    time.var,
    event.var,
    Q.formula=~.,
    Q.SuperLearner.control=list(family=gaussian(),SL.library="SL.lm"),
    cluster.var=NULL,
    obs.weight.var=NULL,
    corstr="independence",
    denom.survival.trunc=1e-3
){
    assert_that(denom.survival.trunc>=0,denom.survival.trunc<=1)
    
    index.shift<-truncation.index-1 #shift for the index of pred_censor.list
    
    models<-lapply(seq_along(tvals),function(i){
        #K is the last visit.time that needs to be considered
        K<-find.last.TRUE.index(visit.times<tvals[i])
        
        history<-reduce(covariates[1:truncation.index],.f=function(d1,d2){
            right_join(d1,d2,by=id.var)
        })%>%arrange(.data[[id.var]])
        
        pred_censor_obj<-pred_censor.list[[K-index.shift]]
        next.visit.time<-tvals[i]
        Y.IPCW<-IPCWtransform(follow.up.time,pred_censor_obj,tvals[i],
                              next.visit.time=next.visit.time,
                              id.var,time.var,event.var,denom.survival.trunc)
        Y.IPCW<-sort_by(Y.IPCW[,1],rownames(Y.IPCW))
        
        IPCW<-rep(1,length(Y.IPCW))
        if(K>truncation.index){
            for(k in (K-1):truncation.index){
                pred_censor_obj<-pred_censor.list[[k-index.shift]]
                k.t<-find.last.TRUE.index(pred_censor_obj$time<=visit.times[k+1],noTRUE=0)
                if(k.t==0){
                    Ghat.t<-1
                }else{
                    Ghat.t<-pred_censor_obj$surv[,k.t]
                    # Ghat.t<-Ghat.t[order(names(Ghat.t))]
                    Ghat.t<-Ghat.t[names(Ghat.t) %in% names(Y.IPCW)]
                    Ghat.t<-sort_by(Ghat.t,names(Ghat.t))
                    Ghat.t<-pmax(Ghat.t,denom.survival.trunc)
                }
                IPCW<-IPCW/Ghat.t
            }
        }
        Y.IPCW<-Y.IPCW*IPCW
        Y<-numeric(nrow(history))
        names(Y)<-pull(history,all_of(id.var))
        Y[names(Y.IPCW)]<-Y.IPCW
        
        X<-model.frame(Q.formula,history%>%filter(.data[[id.var]] %in% names(.env$Y))%>%arrange(.data[[id.var]])%>%select(!all_of(id.var)))
        
        if(is.null(cluster.var)){
            cluster.id<-NULL
        }else{
            cluster.id<-follow.up.time%>%filter(.data[[id.var]] %in% names(.env$Y))%>%arrange(.data[[id.var]])%>%pull(all_of(cluster.var))
            names(cluster.id)<-names(Y)
        }
        
        if(is.null(obs.weight.var)){
            obsWeights<-NULL
        }else{
            obsWeights<-follow.up.time%>%filter(.data[[id.var]] %in% names(.env$Y))%>%arrange(.data[[id.var]])%>%pull(all_of(obs.weight.var))
            names(obsWeights)<-names(Y)
        }
        
        if(ncol(X)==0){
            if(is.null(cluster.var)){
                if(is.null(obsWeights)){
                    est<-mean(Y)
                }else{
                    # message(paste(obs.weight.var,"might not be correctly accounted for in the standard error due to failure fully account for the sampling scheme."))
                    est<-mean(Y*obsWeights)/mean(obsWeights)
                }
            }else{
                .requireNamespace("geepack")
                gee.df<-data.frame(Y=Y,cluster.id=cluster.id)
                if(is.null(obsWeights)){
                    gee.df$weights<-rep(1,length(Y))
                }else{
                    gee.df$weights<-obsWeights
                }
                gee.df<-gee.df%>%arrange(.data$cluster.id) #sort by cluster so that geeglm identifies clusters correctly
                gee<-geepack::geeglm(Y~1,weights=weights,id=cluster.id,family=gaussian(),corstr=corstr,data=gee.df)
                est<-as.numeric(coef(gee))
            }
            model<-intercept_model(est)
            return(model)
        }else{
            SuperLearner.arg<-c(
                list(Y=Y,X=X,obsWeights=obsWeights),
                Q.SuperLearner.control
            )
            if(!is.null(cluster.var)){
                if("cvControl" %in% names(SuperLearner.arg) && 
                   "stratifyCV" %in% names(SuperLearner.arg$cvControl) &&
                   SuperLearner.arg$cvControl$stratifyCV){
                    message("Setting stratifyCV=FALSE in SuperLearner::SuperLearner.CV.control due to cluster.var being specified")
                }
                SuperLearner.arg$cvControl$stratifyCV<-FALSE
            }
            model<-do.call(SuperLearner,SuperLearner.arg)
            return(model)
        }
    })
    names(models)<-as.character(tvals)
    models
}
