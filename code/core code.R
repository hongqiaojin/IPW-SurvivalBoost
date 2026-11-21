############################################################
## Case-Cohort IPW-SurvivalBoost
## C-index + IAUC + MISE
############################################################

rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(survival); library(xgboost)
  library(dplyr); library(timeROC)
})

SB_CLIP <- 10

## =========================================================
## 0. Utilities
## =========================================================

safe_exp <- function(x,clip=20) exp(pmax(pmin(x,clip),-clip))
safe_lp <- function(x,clip=20){
  x[!is.finite(x)] <- NA_real_
  pmax(pmin(x,clip),-clip)
}
has_events <- function(x,n=5) sum(x==1,na.rm=TRUE)>=n
safe_mean <- function(x) if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)
safe_sd <- function(x) if(sum(!is.na(x))<=1) NA_real_ else sd(x,na.rm=TRUE)

safe_q <- function(x,p){
  x<-x[is.finite(x)]
  if(!length(x)) return(NA_real_)
  as.numeric(quantile(x,p,na.rm=TRUE,names=FALSE))
}

trapz <- function(x,y){
  ok<-is.finite(x)&is.finite(y); x<-x[ok]; y<-y[ok]
  if(length(x)<2) return(NA_real_)
  o<-order(x); x<-x[o]; y<-y[o]
  sum(diff(x)*(head(y,-1)+tail(y,-1))/2)
}

## =========================================================
## 1. Evaluation
## =========================================================

get_event_surv <- function(time,status,eps=1e-6){
  f<-tryCatch(survfit(Surv(time,status)~1),
              error=function(e) NULL)
  
  if(is.null(f)||!length(f$time))
    return(function(t) rep(1,length(t)))
  
  tt<-c(0,f$time); ss<-c(1,f$surv)
  
  function(t)
    pmax(ss[pmax(findInterval(t,tt),1L)],eps)
}

get_cindex <- function(time,status,risk){
  risk<-safe_lp(risk,20)
  ok<-is.finite(time)&!is.na(status)&is.finite(risk)
  
  if(sum(ok)<10||sum(status[ok]==1)<2) return(NA_real_)
  
  z<-tryCatch(
    concordance(
      Surv(time[ok],status[ok])~risk[ok],
      reverse=TRUE
    ),error=function(e) NULL)
  
  if(is.null(z)) NA_real_ else as.numeric(z$concordance)
}

get_iauc <- function(time,status,risk,times){
  risk<-safe_lp(risk,20)
  ok<-is.finite(time)&!is.na(status)&is.finite(risk)
  
  if(sum(ok)<10||sum(status[ok]==1)<2) return(NA_real_)
  
  z<-tryCatch(timeROC(
    T=time[ok],delta=status[ok],marker=risk[ok],
    cause=1,times=times,iid=FALSE
  ),error=function(e) NULL)
  
  if(is.null(z)) return(NA_real_)
  
  auc<-as.numeric(z$AUC)
  S<-get_event_surv(time,status)(times)
  ok<-is.finite(auc)&is.finite(S)&is.finite(times)
  
  if(sum(ok)<2) return(NA_real_)
  
  auc<-auc[ok]; S<-S[ok]
  dF<-pmax(head(S,-1)-tail(S,-1),0)
  amid<-(head(auc,-1)+tail(auc,-1))/2
  
  if(sum(dF)<=1e-12) mean(auc)
  else sum(amid*dF)/sum(dF)
}

## True survival + MISE
get_true_surv <- function(dat,times,lambda0=.08){
  eta<-dat$x1-.8*cos(dat$x2+1)+
    .7*abs(dat$x3)+.5*dat$x4^2-.6*sin(dat$x5)
  
  exp(-outer(lambda0*exp(eta),times,"*"))
}

get_mise <- function(pred,true,times){
  if(is.null(pred)||all(is.na(pred))) return(NA_real_)
  
  pred<-as.matrix(pred); true<-as.matrix(true)
  
  mse<-sapply(seq_along(times),function(k){
    ok<-is.finite(pred[,k])&is.finite(true[,k])
    if(!sum(ok)) return(NA_real_)
    mean((pred[ok,k]-true[ok,k])^2)
  })
  
  ok<-is.finite(times)&is.finite(mse)
  if(sum(ok)<2) return(NA_real_)
  
  trapz(times[ok],mse[ok])/diff(range(times[ok]))
}

## =========================================================
## 2. Data generation
## =========================================================

sim_rightcens_data <- function(n=1000,p=20,
                               event_rate_target=.15,
                               lambda0=.08,seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  if(p<5) stop("p must be >= 5.")
  
  X<-matrix(rnorm(n*p),n,p)
  colnames(X)<-paste0("x",1:p)
  
  eta<-X[,1]-.8*cos(X[,2]+1)+
    .7*abs(X[,3])+.5*X[,4]^2-.6*sin(X[,5])
  
  Te<--log(runif(n))/(lambda0*exp(eta))
  
  f<-function(lr)
    mean(exp(-exp(lr)*Te))-event_rate_target
  
  lr<-tryCatch(
    uniroot(f,c(-15,15))$root,
    error=function(e) NA_real_
  )
  
  if(!is.finite(lr)) stop("Censoring calibration failed.")
  
  C<-rexp(n,exp(lr))
  
  data.frame(
    time=pmin(Te,C),
    status=as.integer(Te<=C),
    X
  )
}

## =========================================================
## 3. Case-cohort
## =========================================================

make_casecohort_data <- function(dat,q=.20,seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  
  sub<-rbinom(nrow(dat),1,q)
  xi<-as.integer(dat$status==1|sub==1)
  pi<-dat$status+(1-dat$status)*q
  
  omega<-xi/pi
  omega[xi==0]<-0
  
  out<-transform(
    dat,subcohort=sub,xi=xi,
    pi_obs=pi,omega=omega
  )
  
  list(
    full=out,
    observed=out[out$xi==1,,drop=FALSE],
    sub_only=out[out$subcohort==1,,drop=FALSE],
    q=q
  )
}

## =========================================================
## 4. Weighted Cox objective -- unchanged
## =========================================================

weighted_cox_core <- function(preds,time,status,weight,
                              clip=10,hess_floor=1e-6){
  
  n<-length(preds)
  grad<-numeric(n)
  hess<-rep(hess_floor,n)
  
  keep<-is.finite(time)&!is.na(status)&
    is.finite(weight)&weight>0
  
  if(!sum(keep)||!sum(status[keep]==1))
    return(list(grad=grad,hess=hess,loss=NA_real_))
  
  id<-which(keep)
  tt<-time[keep]; dd<-status[keep]
  ww<-weight[keep]; lp<-preds[keep]
  
  lp[!is.finite(lp)]<-0
  lp<-pmax(pmin(lp,clip),-clip)
  
  o<-order(tt); inv<-order(o)
  tt<-tt[o]; dd<-dd[o]; ww<-ww[o]; lp<-lp[o]
  
  r<-ww*exp(lp)
  risk<-rev(cumsum(rev(r)))
  
  grp<-match(tt,unique(tt))
  first<-which(!duplicated(tt))
  denom<-pmax(risk[first],1e-12)
  
  de<-as.numeric(rowsum(ww*dd,grp,reorder=FALSE))
  A<-cumsum(de/denom)
  B<-cumsum(de/denom^2)
  
  g<--ww*dd+r*A[grp]
  h<-r*A[grp]-r^2*B[grp]
  
  g[!is.finite(g)]<-0
  h[!is.finite(h)]<-hess_floor
  h<-pmax(h,hess_floor)
  
  grad[id]<-g[inv]
  hess[id]<-h[inv]
  
  list(
    grad=grad,hess=hess,
    loss=-sum(ww*dd*lp)+sum(de*log(denom))
  )
}

make_obj <- function(time,status,weight){
  force(time); force(status); force(weight)
  
  function(preds,dtrain){
    z<-weighted_cox_core(preds,time,status,weight)
    list(grad=z$grad,hess=z$hess)
  }
}

make_eval <- function(time,status,weight){
  force(time); force(status); force(weight)
  
  function(preds,dtrain){
    z<-weighted_cox_core(preds,time,status,weight)
    list(metric="wcoxloss",value=z$loss)
  }
}

## =========================================================
## 5. Breslow + survival
## =========================================================

breslow_baseline <- function(time,status,lp,weight=NULL,
                             clip=SB_CLIP){
  
  if(is.null(weight)) weight<-rep(1,length(time))
  lp<-safe_lp(lp,clip)
  
  keep<-is.finite(time)&!is.na(status)&
    is.finite(lp)&is.finite(weight)&weight>0
  
  time<-time[keep]; status<-status[keep]
  lp<-lp[keep]; weight<-weight[keep]
  
  if(!length(time)||!sum(status==1))
    return(list(time=0,cumhaz=0))
  
  o<-order(time)
  time<-time[o]; status<-status[o]
  lp<-lp[o]; weight<-weight[o]
  
  r<-weight*safe_exp(lp,clip)
  risk<-rev(cumsum(rev(r)))
  
  grp<-match(time,unique(time))
  first<-which(!duplicated(time))
  denom<-pmax(risk[first],1e-12)
  
  de<-as.numeric(rowsum(weight*status,grp,reorder=FALSE))
  ut<-unique(time); ev<-de>0
  
  list(
    time=ut[ev],
    cumhaz=cumsum(de[ev]/denom[ev])
  )
}

predict_surv_lp <- function(lp,bh,times,clip=SB_CLIP){
  lp<-safe_lp(lp,clip)
  
  H0<-approx(
    c(0,bh$time),c(0,bh$cumhaz),
    xout=times,method="constant",
    rule=2,f=0
  )$y
  
  S<-exp(-outer(safe_exp(lp,clip),H0,"*"))
  S[!is.finite(S)]<-NA_real_
  
  pmin(pmax(S,0),1)
}

## =========================================================
## 6. Cox
## =========================================================

fit_cox <- function(dat,features,weights=NULL,robust=FALSE){
  if(is.null(dat)||!nrow(dat)||!has_events(dat$status))
    return(NULL)
  
  d<-dat[,c("time","status",features),drop=FALSE]
  keep<-complete.cases(d)
  
  if(!is.null(weights)){
    keep<-keep&is.finite(weights)&weights>0
    weights<-weights[keep]
  }
  
  d<-d[keep,,drop=FALSE]
  if(!nrow(d)||!has_events(d$status)) return(NULL)
  
  form<-as.formula(
    paste0("Surv(time,status)~",paste(features,collapse="+"))
  )
  
  tryCatch(suppressWarnings(
    coxph(
      form,data=d,weights=weights,robust=robust,
      x=TRUE,ties="breslow",singular.ok=TRUE,
      control=coxph.control(
        iter.max=50,eps=1e-7,toler.chol=1e-9
      )
    )
  ),error=function(e) NULL)
}

predict_cox_lp <- function(fit,newdat){
  if(is.null(fit)) return(rep(NA_real_,nrow(newdat)))
  
  safe_lp(tryCatch(
    as.numeric(predict(
      fit,newdata=newdat,type="lp",reference="zero"
    )),
    error=function(e) rep(NA_real_,nrow(newdat))
  ),20)
}

predict_cox_surv <- function(fit,risk,times,n){
  if(is.null(fit)||all(is.na(risk)))
    return(matrix(NA_real_,n,length(times)))
  
  bh<-tryCatch(basehaz(fit,centered=FALSE),
               error=function(e) NULL)
  
  if(is.null(bh)||!nrow(bh))
    return(matrix(NA_real_,n,length(times)))
  
  H0<-approx(
    c(0,bh$time),c(0,bh$hazard),
    xout=times,method="constant",
    rule=2,f=0
  )$y
  
  S<-exp(-outer(safe_exp(risk,20),H0,"*"))
  S[!is.finite(S)]<-NA_real_
  
  pmin(pmax(S,0),1)
}

## =========================================================
## 7. SurvivalBoost base model
## =========================================================

train_sb_base <- function(
    x,time,status,weight,
    nrounds=300,eta=.03,
    max_depth=3,min_child_weight=10,
    subsample=.8,colsample_bytree=.8,
    lambda=1,alpha=0,seed=1){
  
  x<-as.matrix(x)
  
  keep<-complete.cases(x)&is.finite(time)&
    !is.na(status)&is.finite(weight)&weight>0
  
  x<-x[keep,,drop=FALSE]
  time<-time[keep]; status<-status[keep]; weight<-weight[keep]
  
  if(!nrow(x)||!has_events(status)) return(NULL)
  
  dtrain<-xgb.DMatrix(x,label=time)
  
  args<-list(
    params=list(
      eta=eta,max_depth=max_depth,
      min_child_weight=min_child_weight,
      subsample=subsample,
      colsample_bytree=colsample_bytree,
      lambda=lambda,alpha=alpha,
      base_score=0,seed=seed,
      disable_default_eval_metric=1
    ),
    data=dtrain,nrounds=nrounds,
    obj=make_obj(time,status,weight),
    maximize=FALSE,verbose=0
  )
  
  ev<-make_eval(time,status,weight)
  
  if("custom_metric"%in%names(formals(xgb.train)))
    args$custom_metric<-ev
  else
    args$feval<-ev
  
  fit<-tryCatch(
    suppressWarnings(do.call(xgb.train,args)),
    error=function(e) NULL
  )
  
  if(is.null(fit)) return(NULL)
  
  lp<-tryCatch(
    as.numeric(predict(fit,dtrain)),
    error=function(e) rep(NA_real_,nrow(x))
  )
  
  list(
    model=fit,
    lp_train=safe_lp(lp,20)
  )
}

predict_sb_risk <- function(fit,x){
  if(is.null(fit))
    return(rep(NA_real_,nrow(x)))
  
  safe_lp(tryCatch(
    as.numeric(predict(
      fit$model,xgb.DMatrix(as.matrix(x))
    )),
    error=function(e) rep(NA_real_,nrow(x))
  ),20)
}

## =========================================================
## 8. SurvivalBoost final model
## =========================================================

train_sb <- function(dat,features,weight,seed=1){

  fit<-train_sb_base(
    dat[,features,drop=FALSE],
    dat$time,dat$status,weight,
    seed=seed+100
  )

  if(is.null(fit)) return(NULL)

  bh<-breslow_baseline(
    dat$time,dat$status,
    fit$lp_train,
    weight
  )

  list(
    model=fit$model,
    basehaz=bh
  )
}

predict_sb <- function(fit,x,times){

  if(is.null(fit))
    return(list(
      risk=rep(NA_real_,nrow(x)),
      surv=matrix(NA_real_,nrow(x),length(times))
    ))

  risk<-safe_lp(tryCatch(
    as.numeric(predict(
      fit$model,xgb.DMatrix(as.matrix(x))
    )),
    error=function(e) rep(NA_real_,nrow(x))
  ),20)

  list(
    risk=risk,
    surv=predict_surv_lp(
      risk,
      fit$basehaz,
      times
    )
  )
}

## =========================================================
## 9. Evaluation grid
## =========================================================

make_eval_grid <- function(dat,lo=.30,hi=.80,n_times=30){
  
  et<-dat$time[dat$status==1]
  t1<-safe_q(et,lo)
  t2<-safe_q(et,hi)
  
  if(!is.finite(t1)||!is.finite(t2)||t2<=t1)
    stop("Unable to construct evaluation grid.")
  
  seq(t1,t2,length.out=n_times)
}

## =========================================================
## 10. Fit/evaluate all models
## =========================================================

fit_all_models <- function(train,test,q=.20,seed=1){
  
  features<-grep("^x",names(train),value=TRUE)
  cc<-make_casecohort_data(train,q,seed)
  
  full<-cc$full
  obs<-cc$observed
  sub<-cc$sub_only
  
  set.seed(seed+100)
  srs<-train[
    sample(seq_len(nrow(train)),nrow(obs)),
    ,drop=FALSE
  ]
  
  times<-make_eval_grid(train)
  trueS<-get_true_surv(test,times)
  
  ## Cox
  cf<-fit_cox(full,features)
  cs<-fit_cox(sub,features)
  cr<-fit_cox(srs,features)
  ci<-fit_cox(obs,features,obs$omega,TRUE)
  
  rcf<-predict_cox_lp(cf,test)
  rcs<-predict_cox_lp(cs,test)
  rcr<-predict_cox_lp(cr,test)
  rci<-predict_cox_lp(ci,test)
  
  scf<-predict_cox_surv(cf,rcf,times,nrow(test))
  scs<-predict_cox_surv(cs,rcs,times,nrow(test))
  scr<-predict_cox_surv(cr,rcr,times,nrow(test))
  sci<-predict_cox_surv(ci,rci,times,nrow(test))
  
  ## SurvivalBoost
  sf<-train_sb(
    full,features,rep(1,nrow(full)),
    seed+1
  )

  ss<-train_sb(
    sub,features,rep(1,nrow(sub)),
    seed+2
  )

  sr<-train_sb(
    srs,features,rep(1,nrow(srs)),
    seed+3
  )

  si<-train_sb(
    obs,features,obs$omega,
    seed+4
  )

  xt<-test[,features,drop=FALSE]
  
  psf<-predict_sb(sf,xt,times)
  pss<-predict_sb(ss,xt,times)
  psr<-predict_sb(sr,xt,times)
  psi<-predict_sb(si,xt,times)
  
  metric<-function(method,risk,surv){
    data.frame(
      method=method,
      cindex=get_cindex(test$time,test$status,risk),
      auc=get_iauc(test$time,test$status,risk,times),
      mise=get_mise(surv,trueS,times)
    )
  }

  res<-bind_rows(
    metric("Cox_Full",rcf,scf),
    metric("Cox_Sub",rcs,scs),
    metric("Cox_SRS",rcr,scr),
    metric("Cox_IPW",rci,sci),

    metric("SB_Full",psf$risk,psf$surv),
    metric("SB_Sub",pss$risk,pss$surv),
    metric("SB_SRS",psr$risk,psr$surv),
    metric("SB_IPW",psi$risk,psi$surv)
  )

  list(
    metric=res,
    extra=list(
      obs_n=nrow(obs),
      sub_n=nrow(sub),
      eval_t1=min(times),
      eval_t2=max(times)
    )
  )
}

## =========================================================
## 11. Simulation
## =========================================================

run_one <- function(n,p,event_rate,
                    train_ratio=.7,q=.20,seed=2026){
  
  dat<-sim_rightcens_data(
    n=n,p=p,
    event_rate_target=event_rate,
    seed=seed
  )
  
  set.seed(seed+1000)
  id<-sample(seq_len(n),floor(train_ratio*n))
  
  fit_all_models(
    dat[id,,drop=FALSE],
    dat[-id,,drop=FALSE],
    q=q,
    seed=seed+2000
  )
}

run_study <- function(
    nsim=100,
    n_levels=c(1000,1500,2000),
    p_levels=c(10,20),
    event_levels=c(.10,.15,.20),
    train_ratio=.7,q=.20,
    base_seed=2026){
  
  out<-list(); k<-1L
  
  for(ni in seq_along(n_levels))
    for(pi in seq_along(p_levels))
      for(ei in seq_along(event_levels)){
        
        n<-n_levels[ni]
        p<-p_levels[pi]
        pe<-event_levels[ei]
        
        cat("n =",n," p =",p," event =",pe,"\n")
        
        for(b in seq_len(nsim)){
          
          one<-tryCatch(
            run_one(
              n,p,pe,train_ratio,q,
              base_seed+
                100000L*ni+
                10000L*pi+
                1000L*ei+b
            ),
            error=function(e) NULL
          )
          
          if(is.null(one)) next
          
          z<-one$metric
          z$rep<-b
          z$n<-n
          z$p<-p
          z$event_target<-pe
          
          out[[k]]<-z
          k<-k+1L
        }
      }
  
  if(!length(out)) stop("All simulations failed.")
  
  detail<-bind_rows(out)
  
  summary<-detail %>%
    group_by(n,p,event_target,method) %>%
    summarise(
      cindex_mean=safe_mean(cindex),
      cindex_sd=safe_sd(cindex),
      auc_mean=safe_mean(auc),
      auc_sd=safe_sd(auc),
      mise_mean=safe_mean(mise),
      mise_sd=safe_sd(mise),
      valid_n=sum(!is.na(cindex)),
      .groups="drop"
    )
  
  list(detail=detail,summary=summary)
}

## =========================================================
## 12. Final run
## =========================================================

set.seed(2026)

sim_out<-run_study(
  nsim=500,                 # confirm trend first
  n_levels=c(1000,1500,2000),
  p_levels=c(10,20),
  event_levels=c(.10,.15,.20),
  train_ratio=.7,
  q=.20,
  base_seed=2026
)

print(sim_out$summary,n=Inf)

write.csv(
  sim_out$summary,
  "C:/Users/lenovo/Desktop/cohort_new/q_0.2_summary.csv",
  row.names=FALSE
)

write.csv(
  sim_out$detail,
  "C:/Users/lenovo/Desktop/cohort_new/q_0.2_detail.csv",
  row.names=FALSE
)
