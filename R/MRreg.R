#' @title Doubly robust transformation
#' @name DRtransform
#' @description
#' Given a `pred_event_censor` object in a time window, calculates the doubly robust transformation in the time window. The transformation is used as the outcome when estimating the conditional survival probability at the next visit time.
#' @param follow.up.time see \code{\link{MRsurv}}
#' @param pred_event_censor_obj a `pred_event_censor` object in the time window of interest
#' @param tvals see \code{\link{MRsurv}}. Must be greater than the smallest time in `pred_event_censor_obj`
#' @param next.visit.time the next visit time. Default is `Inf`, corresponding to the last time window
#' @param id.var see \code{\link{MRsurv}}
#' @param time.var see \code{\link{MRsurv}}
#' @param event.var see \code{\link{MRsurv}}
#' @param denom.survival.trunc see \code{\link{MRsurv}}
#' @return a named one-column matrix of transformations used for regression. Each row corresponds to an individual. Row names are elements in `follow.up.time[[id.var]]`
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{MRsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @export
DRtransform<-function(follow.up.time,pred_event_censor_obj,tvals,next.visit.time=Inf,id.var,time.var,event.var,denom.survival.trunc){
    tvals.bar<-pmin(tvals,next.visit.time) #tvals truncated at next visit time
    
    #compute Ghat(s-) at s=event times
    Ghat.minus<-pred_event_censor_obj$event$surv
    for(i in 1:length(pred_event_censor_obj$event$time)){
        #j is the index of the first censoring time in censoring times that is >= current event time
        #will use j-1 to get Ghat(s-)
        #if no censoring time >= current event time, then set j to be the index after the last censoring time
        #if j=1 (all censoring time >= current event time), set Ghat(s-) to 1
        j<-find.first.TRUE.index(pred_event_censor_obj$censor$time>=pred_event_censor_obj$event$time[i],
                                 noTRUE=length(pred_event_censor_obj$censor$time)+1)
        if(j==1){
            Ghat.minus[,i]<-1
        }else{
            Ghat.minus[,i]<-pred_event_censor_obj$censor$surv[,j-1]
        }
    }
    
    #compute integrand in the DR transform at each event time
    integrand<-pred_event_censor_obj$event$surv
    for(i in 1:length(pred_event_censor_obj$event$time)){
        if(i==1){
            integrand[,i]<-(1-pred_event_censor_obj$event$surv[,i])/pred_event_censor_obj$event$surv[,i]/pmax(Ghat.minus[,i],denom.survival.trunc)
        }else{
            integrand[,i]<-(pred_event_censor_obj$event$surv[,i-1]-pred_event_censor_obj$event$surv[,i])/
                pred_event_censor_obj$event$surv[,i]/
                pred_event_censor_obj$event$surv[,i-1]/
                pmax(Ghat.minus[,i],denom.survival.trunc)
        }
    }
    #integral in the DR transform at each event time
    integral<-rowCumsums(integrand)
    
    output<-matrix(nrow=nrow(pred_event_censor_obj$event$surv),ncol=length(tvals.bar))
    rownames(output)<-rownames(pred_event_censor_obj$event$surv)
    colnames(output)<-as.character(tvals)
    for(i in 1:nrow(output)){
        id.matching.data<-follow.up.time%>%filter(.data[[id.var]]==rownames(output)[i])
        X<-pull(id.matching.data,all_of(time.var))
        Delta<-pull(id.matching.data,all_of(event.var))
        for(j in 1:ncol(output)){
            if(j>1 && tvals.bar[j]==tvals.bar[j-1]){
                output[i,j]<-output[i,j-1]
            }else{
                #find Shat at t
                #k.t is the index of the last event time in event times that is <= current tvals (t)
                #will use k.t to get Shat at t
                #if no event time < t, then set Shat.t to be 1
                k.t<-find.last.TRUE.index(pred_event_censor_obj$event$time<=tvals.bar[j],noTRUE=0)
                if(k.t==0){
                    Shat.t<-1
                }else{
                    Shat.t<-pred_event_censor_obj$event$surv[i,k.t]
                }
                
                #find Shat at X
                #same logic as above
                k.SX<-find.last.TRUE.index(pred_event_censor_obj$event$time<=X,noTRUE=0)
                if(k.SX==0){
                    Shat.X<-1
                }else{
                    Shat.X<-pred_event_censor_obj$event$surv[i,k.SX]
                }
                
                #compute first IPW term in the bracket
                if(Delta==1 && X<=tvals.bar[j]){
                    #find Ghat(X-)
                    #same logic as above
                    k.GX<-find.first.TRUE.index(pred_event_censor_obj$censor$time>=X,
                                                noTRUE=length(pred_event_censor_obj$censor$time)+1)
                    if(k.GX==1){
                        Ghat.Xminus<-1
                    }else{
                        Ghat.Xminus<-pred_event_censor_obj$censor$surv[i,k.GX-1]
                    }
                    
                    IPW.term<-1/Shat.X/pmax(Ghat.Xminus,denom.survival.trunc)
                }else{
                    IPW.term<-0
                }
                
                #find integral at min(X,t)
                k.int<-min(k.t,k.SX)
                if(k.int==0){
                    integral.Xt<-0
                }else{
                    integral.Xt<-integral[i,k.int]
                }
                
                if(Shat.t!=0){
                    output[i,j]<-Shat.t-Shat.t*(IPW.term-integral.Xt)
                }else{
                    output[i,j]<-Shat.t
                }
            }
        }
    }
    
    output
}


#' @title Regression based on multiply (sequentially) robust transformation of fitted survival and censoring probabilities using SuperLearner
#' @name MRreg.SuperLearner
#' @description Apply sequentially doubly robust transformation on fitted survival and censoring probabilities in each time window and estimate P(T > t | T > tk, covariates available at tk) with \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}.
#' @param covariates see \code{\link{MRsurv}}
#' @param follow.up.time see \code{\link{MRsurv}}
#' @param pred_event_censor.list list of `pred_event_censor` objects (see \code{\link{pred_event_censor}}). Each `pred_event_censor` object in the list corresponds to a time window in `visit.times` after `truncation.index` in increasing order.
#' @param visit.times see \code{\link{MRsurv}}
#' @param tvals see \code{\link{MRsurv}}
#' @param truncation.index see \code{\link{MRsurv}}
#' @param id.var see \code{\link{MRsurv}}
#' @param time.var see \code{\link{MRsurv}}
#' @param event.var see \code{\link{MRsurv}}
#' @param U.formula see \code{\link{MRsurv}}
#' @param Q.formula see \code{\link{MRsurv}}
#' @param U.SuperLearner.control see \code{\link{MRsurv}}
#' @param Q.SuperLearner.control see \code{\link{MRsurv}}
#' @param U.folds a list of vectors of id (identified by variable `id.var`) corresponding to each fold for cross-fitting. Set to a list containing one vector for no cross-fitting.
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
MRreg.SuperLearner<-function(
    covariates,
    follow.up.time,
    pred_event_censor.list,
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
    cluster.var=NULL,
    obs.weight.var=NULL,
    corstr="independence",
    denom.survival.trunc=1e-3
){
    assert_that(denom.survival.trunc>=0,denom.survival.trunc<=1)
    
    index.shift<-truncation.index-1 #shift for the index of pred_event_censor.list
    
    models<-lapply(seq_along(tvals),function(i){
        #K is the last visit.time that needs to be considered
        K<-find.last.TRUE.index(visit.times<tvals[i])
        
        for(k in K:truncation.index){
            if(k<length(visit.times)){
                pred_event_censor_obj<-truncate_pred_event_censor(pred_event_censor.list[[k-index.shift]],visit.times[k+1])
            }else{
                pred_event_censor_obj<-pred_event_censor.list[[k-index.shift]]
            }
            if(k==K){
                next.visit.time<-tvals[i]
            }else{
                next.visit.time<-visit.times[k+1]
            }
            
            history<-reduce(covariates[1:ifelse(k==truncation.index,k,k-1)],.f=function(d1,d2){
                right_join(d1,d2,by=id.var)
            })%>%arrange(.data[[id.var]])
            
            Y.DR<-DRtransform(follow.up.time,pred_event_censor_obj,tvals[i],
                              next.visit.time=next.visit.time,
                              id.var,time.var,event.var,denom.survival.trunc)
            # Y.DR<-Y.DR[order(rownames(Y.DR)),1]
            Y.DR<-sort_by(Y.DR[,1],rownames(Y.DR))
            
            if(k==K){
                Y<-Y.DR
            }else{
                if(length(U.folds)==1){
                    U<-as.numeric(predict(model,newdata=history%>%select(!all_of(id.var)),onlySL=TRUE)$pred)
                    names(U)<-history%>%pull(all_of(id.var))
                }else{
                    U.list<-lapply(1:length(U.folds),function(v){
                        newdata<-history%>%filter(.data[[id.var]] %in% U.folds[[v]])
                        U<-as.numeric(predict(models[[v]],newdata=newdata%>%select(!all_of(id.var)),onlySL=TRUE)$pred)
                        names(U)<-newdata%>%pull(all_of(id.var))
                        U
                    })
                    U<-do.call(c,U.list)
                    # U<-U[order(names(U))]
                    U<-sort_by(U,names(U))
                }
                
                Y1<-numeric(nrow(history))
                names(Y1)<-history%>%pull(all_of(id.var))
                Y1[names(Y)]<-Y-U[names(Y)]
                
                k.t<-find.last.TRUE.index(pred_event_censor_obj$censor$time<=visit.times[k+1],noTRUE=0)
                if(k.t==0){
                    Ghat.t<-1
                }else{
                    Ghat.t<-pred_event_censor_obj$censor$surv[,k.t]
                    # Ghat.t<-Ghat.t[order(names(Ghat.t))]
                    Ghat.t<-sort_by(Ghat.t,names(Ghat.t))
                    Ghat.t<-pmax(Ghat.t,denom.survival.trunc)
                }
                Y1<-Y1/Ghat.t
                Y<-Y1+U*Y.DR
                names(Y)<-names(U)
            }
            
            
            if(k>truncation.index){
                form<-U.formula[[k-truncation.index]]
            }else{
                form<-Q.formula
            }
            train.data<-history%>%filter(.data[[id.var]] %in% names(.env$Y))%>%arrange(.data[[id.var]])
            X<-model.frame(form,train.data%>%select(!all_of(id.var)))
            
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
            
            if(k==truncation.index && ncol(X)==0){
                pseudo.outcome<-Y
                if(is.null(cluster.var)){
                    if(is.null(obsWeights)){
                        est<-mean(Y)
                        IF<-Y-est
                    }else{
                        # message(paste(obs.weight.var,"is normalized to have sample mean 1"))
                        est<-mean(Y*obsWeights)/mean(obsWeights)
                        IF<-(mean(obsWeights)*obsWeights*Y-mean(Y*obsWeights)*obsWeights)/mean(obsWeights)^2
                    }
                    SE<-sqrt(mean(IF^2))/sqrt(length(IF))
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
                    IF<-as.numeric(gee$geese$infls[1,])
                    IF<-IF*length(IF)
                    names(IF)<-unique(gee.df$cluster.id)
                    SE<-sqrt(vcov(gee)[1,1])
                }
                
                model<-intercept_IF_model(est,pseudo.outcome,IF,SE)
                return(model)
            }else{
                if(k==truncation.index){
                    SuperLearner.arg<-c(
                        list(Y=Y,X=X,id=cluster.id,obsWeights=obsWeights),
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
                }else{
                    if(length(U.folds)==1){
                        SuperLearner.arg<-c(
                            list(Y=Y,X=X,id=cluster.id,obsWeights=obsWeights),
                            U.SuperLearner.control
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
                    }else{
                        models<-lapply(U.folds,function(fold){
                            X<-model.frame(form,train.data%>%filter(!(.data[[id.var]] %in% fold))%>%select(!all_of(id.var)))
                            SuperLearner.arg<-c(
                                list(Y=Y[!(names(Y) %in% fold)],X=X,id=cluster.id[!(names(cluster.id) %in% fold)],obsWeights=obsWeights[!(names(obsWeights) %in% fold)]),
                                U.SuperLearner.control
                            )
                            if(!is.null(cluster.var)){
                                if("cvControl" %in% names(SuperLearner.arg) && 
                                   "stratifyCV" %in% names(SuperLearner.arg$cvControl) &&
                                   SuperLearner.arg$cvControl$stratifyCV){
                                    message("Setting stratifyCV=FALSE in SuperLearner::SuperLearner.CV.control due to cluster.var being specified")
                                }
                                SuperLearner.arg$cvControl$stratifyCV<-FALSE
                            }
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
