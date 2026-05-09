library(dplyr)
library(tidyr)
library(survival)


#storage for simulation results

Nsim <- 1
#descriptives
descriptives_scenario <- matrix(nrow=Nsim, ncol=20)
colnames(descriptives_scenario) <- rep("NA",20)
descr_KM_curves_long <- matrix(nrow = 0, ncol = 4)
colnames(descr_KM_curves_long) = c("sim","time","type","risk")

#calibration
#[ , ,1]: predicted
#[ , ,2]: IPCW observed
#[ , ,3]: counterfactual true
#[ , ,4]: subset predicted
#[ , ,5]: subset observed
calib_risk0<-array(dim=c(Nsim,5,5))
calib_risk1<-array(dim=c(Nsim,5,5))
calib_risk0_group<-array(dim=c(Nsim,10,5))
calib_risk1_group<-array(dim=c(Nsim,10,5))

#discrimination
disc_cindex<-matrix(nrow=Nsim,ncol=10)
disc_auct<-matrix(nrow=Nsim,ncol=10)

#overall prediction error
brier_raw <- matrix(nrow=Nsim, ncol=8)
brier_ipa <- matrix(nrow=Nsim, ncol=8)

#---------------
#function for c-index with inverse prob weight option, allowing tied event times

c_index_ties <- function(time, status, risk, weightmatrix=NULL, tau)
  #Input:
  # time: event/censoring times
  # status: event indicators
  # risk: vector of risk probabilities by time horizon (typically tau) for each individual.
  # weightmatrix: matrix with ipc weights. Rows are subjects, columns are unique event time points (ordered)
  # tau: truncation time point, calculate c-index from zero up to (and including) tau
  #Output:
  # c-index
{
  tt <- sort(unique(time[status == 1 & time <= tau]))
  if(is.null(weightmatrix))
  {
    weightmatrix=matrix(1,nrow=length(time),ncol = length(tt))
  }
  nt <- length(tt)            #number of unique event time points
  x <- risk
  numsum <- denomsum <- 0
  for (i in 1:nt)                   #loop over unique event time points
  {
    ti <- tt[i]
    n1 <- intersect(which(time==ti) ,which(status==1)) #indices of cases at this time point
    n0 <- union(which(time>ti), intersect(which(time==ti),which(status==0)))  	#indices controls at this time point (patients with event time >ti plus patients censored at ti)
    nn1 <- length(n1)                 #number of cases
    nn0 <- length(n0)                 #number of controls
    x <- risk
    for (k in 1:nn1)                  #loop over the cases
    {
      xi <- x[n1][k] # risk of the k'th case at time point ti
      numsum  <- numsum + weightmatrix[n1,i][k] * (x[n0]<xi)%*%weightmatrix[n0,i] + 0.5 * weightmatrix[n1,i][k] * (x[n0] == xi) %*% weightmatrix[n0,i]
      denomsum <- denomsum + weightmatrix[n1,i][k]*sum(weightmatrix[n0,i])

    }

  }
  cindex <- numsum/denomsum
  return(cindex)
}


#---------------
#function for cumulative dynamic AUCt with ipc weights

wCD_AUCt <- function (time, status, risk, seq.time, plot = T, weightmatrix=NULL, xlim=5)
  #Input:
  # time: event/censoring times
  # status: event indicators
  # risk: risk is a vector of risks for each subject (most implementations of C/D AUC just allow a vector, eg timeROC). But could be a matrix with predicitons by each unique event time point (our risk_obs_allt does not contain these for all observed event time points, only up to last follow up of each patient)
  # weightsmatrix: rows are subjects, weights are ordered unique event time points plus the evaluation time point appended (t=4.99 in this case, so weight at t=4)
  # seq.time: time vector where you want to calculate C/D AUC (must be ordered and smallest time should be equal or larger than first event time)
  # xlim: plotting parameter
  #Output:
  # C/D AUCt values at seq.times
  # plot
{
  tt <- sort(unique(time[status==1])) #unique event time points
  if(is.null(weightmatrix))
  {
    weightmatrix=matrix(1,nrow=length(time),ncol = length(tt))
  }
  nseq <- length(seq.time)
  AUCt <- rep(NA,nseq)                  #vector to save AUCt in
  for (i in 1:nseq)                   #loop over time points where you want to calculate C/D AUCt
  {
    numsum <- denomsum <- 0
    ti <- seq.time[i]                 #ith unique time point where you want to calcuate C/D AUCt
    #tti <- which(tt<=ti))         #not used could be index among tt of this seq.time (then tt should be expanded by seq.time) and used in the weights of the controls below
    n1 <- intersect(which(time<=ti) ,which(status==1)) #indices of cumulative(!) cases at this time point
    n0 <- union(which(time>ti), intersect(which(time==ti),which(status==0)))  	#indices of controls at this time point (patients with event time >ti plus patients censored at ti)
    nn1 <- length(n1)                 #number of cases
    nn0 <- length(n0)                 #number of controls
    x <- risk                         #x is vector of risks to be evaluated (could be specific to i'th unique time point)
    for (k in 1:nn1)                  #loop over the cases
    {
      n1k <- n1[k]                    #index in original data of the k'th case
      xi <- x[n1k]                    # risk of the k'th case
      ttk <- which(tt==time[n1k])     #index among unique event times for the k'th case
      numsum   <- numsum   + weightmatrix[n1k,ttk] * (x[n0]<xi)%*%weightmatrix[n0,ncol(weightmatrix)] + 0.5 * weightmatrix[n1k,ttk] * (x[n0] == xi) %*% weightmatrix[n0,ncol(weightmatrix)]
      denomsum <- denomsum + weightmatrix[n1k,ttk] * sum(weightmatrix[n0,ncol(weightmatrix)])
      #note that the weight for the case is evaluated at it's own unique event time point
      #the weight for the controls are evaluated at the time point where you want to calculate C/D AUCt
    }
    AUCt[i] <- numsum/denomsum
  }

  if(plot)
  {
    plot(seq.time,AUCt,xlab = "Time t", ylab = "C/D AUC(t)", type='l', xlim =c(0,xlim), ylim=c(0.5,1), main="cumulative dynamic AUCt")
    abline(h=0.5,lty=3)
  }
  return(list(AUCt=data.frame(time=seq.time,AUC=AUCt)))
}

#----------------------
#function for Brier score

Brier <- function(time, status, risk, seq.time, weights=rep(1, length(time)))
{
  uncensored <- time>seq.time | (status>0 & !is.na(risk))
  status[time>seq.time] <- 0  # only evaluate event status up to seq.time

  # original code uses 1/length(time), we use sum(weights) as its slightly more
  # accurate
  print(sum(weights[uncensored]))
  return(1/(sum(weights[uncensored]))*sum((risk[uncensored]-status[uncensored])^2*weights[uncensored]))
}

#----------------------
#function for scaled Brier score (IPA)

ipa <- function(time, status, risk, seq.time, weights=rep(1, length(time)))
{
  #sf <- survfit(Surv(time, status)~1,weights=weights)
  #nullrisk <- 1-min(sf$surv[sf$time<=seq.time])
  uncensored <- time>seq.time | (status>0 & !is.na(risk))
  nullrisk <- 1/(sum(weights[uncensored]))*sum(status[uncensored]*weights[uncensored])
  brier1 <- Brier(time, status, risk, seq.time, weights=weights)
  brier0 <- Brier(time, status, rep(nullrisk,length(time)), seq.time, weights=weights)
  return((brier0-brier1)/brier0*100)
}
n <- 3000
scenario <- 1; cindex_ylim_low <- .48; cindex_ylim_high <- .675; auc_ylim_low <- 0.54; auc_ylim_high <- 0.75; calib_lim_low <- 0.1; calib_lim_high <- 0.9; calib_risk_lim_high <- 0.8 ;calib_events_lim_high <- (n/1000)*1300; model = " Cox model"


# dat_sim_cox_scenario1 ---------------------------------------------------
U=rnorm(n,0,0.1)


dat_sim_cox_scenario1 <- function() {
  n.visit=5

  #----
  #expit function

  expit=function(x){exp(x)/(1+exp(x))}

  #------------------
  #parameter values
  #------------------

  #model for A|L
  gamma.0=-1
  gamma.L=0.5

  #model for hazard
  alpha.0=-2
  alpha.A=-0.5
  alpha.L=0.5
  alpha.U=0.5

  #test scenario with no predictive value of the model (AUC and c-index around 0.5)
  #alpha.A=0
  #alpha.L=0
  #alpha.U=0

  #test scenario with high predictive value of the model (c-index around 0.7)
  #alpha.A=-0.1
  #alpha.L=1.5
  #alpha.U=0.5

  #------------------
  #simulate data
  #------------------

  #----
  #generate U, A, L

  A=matrix(nrow=n,ncol=n.visit)
  L=matrix(nrow=n,ncol=n.visit)

  L[,1]=rnorm(n,U,1)
  A[,1]=rbinom(n,1,expit(gamma.0+gamma.L*L[,1]))
  for(k in 2:n.visit){
    L[,k]=rnorm(n,0.8*L[,k-1]-A[,k-1]+0.1*(k-1)+U,1)
    A[,k]=ifelse(A[,k-1]==1,1,rbinom(n,1,expit(gamma.0+gamma.L*L[,k])))
  }

  #----
  #generate event times T.obs, and event indicators D.obs

  T.obs=rep(NA,n)

  for(k in 1:n.visit){
    u.t=runif(n,0,1)
    haz=exp(alpha.0+alpha.A*A[,k]+alpha.L*L[,k]+alpha.U*U)
    new.t=-log(u.t)/haz
    T.obs=ifelse(is.na(T.obs) & new.t<1,k-1+new.t,T.obs)
  }
  D.obs=ifelse(is.na(T.obs),0,1)
  T.obs=ifelse(is.na(T.obs),5,T.obs)

  #-----
  #Create data frame

  colnames(A)=paste0("A.",0:4)
  colnames(L)=paste0("L.",0:4)
  dat=data.frame(id=1:n,T.obs,D.obs,A,L)

  #-----
  #set A to 0 in time periods after event/censoring

  dat$A.1=ifelse(dat$T.obs<1,0,dat$A.1)
  dat$A.2=ifelse(dat$T.obs<2,0,dat$A.2)
  dat$A.3=ifelse(dat$T.obs<3,0,dat$A.3)
  dat$A.4=ifelse(dat$T.obs<4,0,dat$A.4)

  #------------------
  #some summaries: may be useful if you wish to change the parameter values used above, to consider other scenarios.
  #------------------

  #proportion always treated
  always.treat=A[,1]+A[,2]+A[,3]+A[,4]+A[,5]

  #proportion never treated
  never.treat=(1-A[,1])+(1-A[,2])+(1-A[,3])+(1-A[,4])+(1-A[,5])

  #------------------
  #Reshape data into 'long' format (multiple rows per individual: 1 row for each visit)
  #------------------

  dat.long=reshape(data = dat,varying=c(paste0("A.",0:4),paste0("L.",0:4)),direction="long",idvar="id")
  dat.long=dat.long[order(dat.long$id,dat.long$time),]

  #generate start and stop times for each row
  dat.long$time.stop=dat.long$time+1

  dat.long=dat.long[dat.long$time<dat.long$T.obs,]

  dat.long$time.stop=ifelse(dat.long$time.stop>dat.long$T.obs,dat.long$T.obs,dat.long$time.stop)

  dat.long$event=ifelse(dat.long$time.stop==dat.long$T.obs & dat.long$D.obs==1,1,0)

  #visit number
  dat.long$visit=ave(rep(1,dim(dat.long)[1]),dat.long$id,FUN=cumsum)

  #generate lagged A values
  dat.long=dat.long %>%
    group_by(id) %>%
    mutate(Alag1 = lag(A,n=1),Alag2 = lag(A,n=2),Alag3 = lag(A,n=3),Alag4 = lag(A,n=4)) %>%
    mutate(Alag1 = replace_na(Alag1,0),Alag2=replace_na(Alag2,0),Alag3=replace_na(Alag3,0),Alag4=replace_na(Alag4,0))

  #generate lagged L values
  dat.long=dat.long %>%
    group_by(id) %>%
    mutate(Llag1 = lag(L,n=1),Llag2 = lag(L,n=2),Llag3 = lag(L,n=3),Llag4 = lag(L,n=4)) %>%
    mutate(Llag1=replace_na(Llag1,0),Llag2=replace_na(Llag2,0),Llag3=replace_na(Llag3,0),Llag4=replace_na(Llag4,0))

  #baseline L
  dat.long=dat.long %>%
    group_by(id) %>%
    mutate(L.baseline = first(L))

  return(list(dat, dat.long))
}

listdatdatlong <- dat_sim_cox_scenario1()

dat.dev<-listdatdatlong[[1]]
dat.long.dev<-listdatdatlong[[2]]

listvalvallong <- dat_sim_cox_scenario1()
dat.val <- listvalvallong[[1]]
dat.long.val <-listvalvallong[[2]]

dat_sim_cox_counterfactual_scenario1 <- function() {
  #------------------------------
  #------------------------------
  # Simulates counterfactual longitudinal data on covariates L and event times under 'always treated' and 'never treated, using a Cox model.
  #------------------------------
  #------------------------------

  #----
  #number of visits (K+1)

  n.visit=5

  #----
  #expit function

  expit=function(x){exp(x)/(1+exp(x))}

  #------------------
  #parameter values
  #------------------

  #model for A|L
  gamma.0=-1
  gamma.L=0.5

  #model for hazard
  alpha.0=-2
  alpha.A=-0.5
  alpha.L=0.5
  alpha.U=0.5

  #---------------------
  #---------------------
  #generate event times and event indicators had each person been NEVER TREATED
  #use U and L[,1] from validation data
  #---------------------
  #---------------------

  A=matrix(nrow=n,ncol=n.visit)
  L=matrix(nrow=n,ncol=n.visit)

  L[,1]=dat.val$L.0
  A[,1]=0
  for(k in 2:n.visit){
    L[,k]=rnorm(n,0.8*L[,k-1]-A[,k-1]+0.1*(k-1)+U,1)
    A[,k]=0
  }

  T.A0=rep(NA,n)
  for(k in 1:n.visit){
    u.t=runif(n,0,1)
    haz=exp(alpha.0+alpha.A*A[,k]+alpha.L*L[,k]+alpha.U*U)
    new.t=-log(u.t)/haz
    T.A0=ifelse(is.na(T.A0) & new.t<1,k-1+new.t,T.A0)
  }
  D.A0=ifelse(is.na(T.A0),0,1)
  T.A0=ifelse(is.na(T.A0),5,T.A0)

  #---------------------
  #---------------------
  #generate event times and event indicators had each person been ALWAYS TREATED
  #use U and L[,1] from above
  #---------------------
  #---------------------

  A=matrix(nrow=n,ncol=n.visit)
  L=matrix(nrow=n,ncol=n.visit)

  L[,1]=dat.val$L.0
  A[,1]=1
  for(k in 2:n.visit){
    L[,k]=rnorm(n,0.8*L[,k-1]-A[,k-1]+0.1*(k-1)+U,1)
    A[,k]=1
  }

  T.A1=rep(NA,n)
  for(k in 1:n.visit){
    u.t=runif(n,0,1)
    haz=exp(alpha.0+alpha.A*A[,k]+alpha.L*L[,k]+alpha.U*U)
    new.t=-log(u.t)/haz
    T.A1=ifelse(is.na(T.A1) & new.t<1,k-1+new.t,T.A1)
  }
  D.A1=ifelse(is.na(T.A1),0,1)
  T.A1=ifelse(is.na(T.A1),5,T.A1)

  #-----
  #Create data frame

  dat.cf=data.frame(id=1:n,T.A0,D.A0,T.A1,D.A1,L.0=dat.val$L.0)

  return(dat.cf)
}

dat.cf <- dat_sim_cox_counterfactual_scenario1()




# model -------------------------------------------------------------------

wt.mod=glm(A~L,family="binomial",data=dat.long.dev[dat.long.dev$Alag1==0,])
pred.wt=predict(wt.mod,type = "response",newdata = dat.long.dev)
dat.long.dev$wt=ifelse(dat.long.dev$A==1,pred.wt,1-pred.wt)
dat.long.dev$wt=ifelse(dat.long.dev$Alag1==1,1,dat.long.dev$wt)
dat.long.dev$wt.cum=ave(dat.long.dev$wt,dat.long.dev$id,FUN=cumprod)

#Numerator of stabilised weights
wt.mod.num.L=glm(A~L.baseline*as.factor(visit),family="binomial",data=dat.long.dev[dat.long.dev$Alag1==0,])
pred.wt.num.L=predict(wt.mod.num.L,type = "response",newdata = dat.long.dev)
dat.long.dev$wt.num.L=ifelse(dat.long.dev$A==1,pred.wt.num.L,1-pred.wt.num.L)
dat.long.dev$wt.num.L=ifelse(dat.long.dev$Alag1==1,1,dat.long.dev$wt.num.L)
dat.long.dev$wt.cum.num.L=ave(dat.long.dev$wt.num.L,dat.long.dev$id,FUN=cumprod)

#Stabilized weights
dat.long.dev$ipw.s.L=dat.long.dev$wt.cum.num.L/dat.long.dev$wt.cum

#-----------------
#MSM-IPTW analysis using stabilized weights
#-----------------

cox.msm=coxph(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+Alag4+L.baseline,data=dat.long.dev,weights = dat.long.dev$ipw.s.L)

#baseline cumulative hazard
cumhaz=basehaz(cox.msm,centered=F)$hazard
event.times=basehaz(cox.msm,centered=F)$time

#step function giving baseline cumulative hazard at any time
cumhaz.fun=stepfun(event.times,c(0,cumhaz))


# validation --------------------------------------------------------------

L.baseline.dat=dat.long.val$L[dat.long.val$visit==1]


t.hor=1:5

#risks under NEVER TREATED
risk0_exp=sapply(t.hor,FUN=function(x){1-exp(-(cumhaz.fun(x)*exp(cox.msm$coefficients["L.baseline"]*L.baseline.dat)))})

#risks under ALWAYS TREATED
risk1_exp=sapply(t.hor,FUN=function(x){1-exp(-(
  cumhaz.fun(min(x,1))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)+
    (x>=1)*(cumhaz.fun(min(x,2))-cumhaz.fun(1))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["Alag1"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)+
    (x>=2)*(cumhaz.fun(min(x,3))-cumhaz.fun(2))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["Alag1"]+cox.msm$coefficients["Alag2"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)+
    (x>=3)*(cumhaz.fun(min(x,4))-cumhaz.fun(3))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["Alag1"]+cox.msm$coefficients["Alag2"]+cox.msm$coefficients["Alag3"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)+
    (x>=4)*(cumhaz.fun(min(x,5))-cumhaz.fun(4))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["Alag1"]+cox.msm$coefficients["Alag2"]+cox.msm$coefficients["Alag3"]+cox.msm$coefficients["Alag4"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)
))})

#------------------------------
#------------------------------
#OBSERVERD risks under the two treatment strategies 'ALWAYS TREATED' and 'NEVER TREATED'
#------------------------------
#------------------------------

#-----------------
#'Observed' risks up to times 1:5, under the 'NEVER TREATED' strategy
#Obtained using censoring and weighting
#-----------------

#fit weights model for people untreated at time 0
wt.mod=glm(A~L,family = "binomial",data=dat.long.val[dat.long.val$Alag1==0,])



#predicted probability that A[t]=0 conditional on A[t-1]=0 and conditional on L[t]
pred.wt0=predict(wt.mod,type = "response",newdata = dat.long.val)
dat.long.val$wt0=1-pred.wt0

#Obtain the IPW at each time using cumulative product of 1/wt up to that time
dat.long.val = dat.long.val %>% group_by(id) %>% mutate(ipw0=1/cumprod(wt0))

#Now impose the artificial censoring when people deviate from the 'never treated' strategy

dat.long.val = dat.long.val %>%group_by(id) %>%mutate(A.baseline=first(A))
dat.long.val$in.dat.0 = (dat.long.val$A==0)

#---
#weighted Kaplan-Meier - using unstabilized weights
km.0=survfit(Surv(time,time.stop,event)~1, data=dat.long.val %>% filter(in.dat.0==1), weights = ipw0)
step.risk0.obs=stepfun(km.0$time,c(1,km.0$surv))#step function giving survival probability at any time

#estimated 'observed' risk at times 1:5
risk0_obs=1-step.risk0.obs(1:5)


#OBSERVED RISKS
#'Observed' risks up to times 1:5, under the 'ALWAYS TREATED' strategy
#Note that under our data generating mechanism people always continue treatment after they start
#i.e. there are no transitions from A=1 to A=0
#so we don't need time-dependent weights under the 'always treated' strategy

wt.mod.baseline=glm(A~L,family = "binomial",data=dat.long.val[dat.long.val$visit==1,])


#predicted probability that A[0]=1 conditional on A[-1]=0 (which is true for everyone) and conditional on L[0]
pred.wt1.baseline=predict(wt.mod.baseline,type = "response",newdata = dat.long.val[dat.long.val$visit==1,])

#Obtain the IPW at each time (which is the same at each time here)
dat.long.val$wt1 = 0
dat.long.val$wt1[dat.long.val$visit==1] = 1/pred.wt1.baseline
dat.long.val = dat.long.val %>% group_by(id) %>% mutate(ipw1=sum(wt1))

#Now impose the artificial censoring when people deviate from the 'always treated' strategy
#Note once a person starts treatment they always continue, in this simulation study
dat.long.val$in.dat.1 = (dat.long.val$A.baseline==1)

#---
#weighted Kaplan-Meier
km.1=survfit(Surv(time,time.stop,event)~1,data=dat.long.val %>% filter(in.dat.1==1),weights = ipw1)
step.risk1.obs=stepfun(km.1$time,c(1,km.1$surv))#step function giving survival probability at any time

#estimated 'observed' risk at times 1:5
risk1_obs=1-step.risk1.obs(1:5)

risk1_obs



# ipeval ------------------------------------------------------------------

dat.val_outcome <- data.frame(id = dat.val$id, time = dat.val$T.obs,
                              status = dat.val$D.obs)

dat.val_long <- wide_to_long(dat.val,
                             baseline_variables = c("id"),
                             wide_variables = list(
                               "A" = paste0("A.", 0:4),
                               "L" = paste0("L.", 0:4)
                               ),
                             visit_times = 0:4,
                             outcome_times = dat.val_outcome$time)
dat.val_long <- add_lag_terms(dat.val_long, "A")



ipscore_results_0 <- ip_score_long(
  probabilities = risk0_exp[,5],
  data_outcome = dat.val_outcome,
  data_long = dat.val_long,
  time_horizon = 5,
  treatment_formula = A ~ A_lag_1 * L,
  treatment_of_interest = rep(0, 5),
  null_model = TRUE,
  metrics = c("auc", "brier", "scaled_brier", "oeratio")
)

ipscore_results_1 <- ip_score_long(
  probabilities = risk0_exp[,5],
  data_outcome = dat.val_outcome,
  data_long = dat.val_long,
  time_horizon = 5,
  treatment_formula = A ~ A_lag_1 * L,
  treatment_of_interest = rep(1, 5),
  null_model = FALSE,
  metrics = c("auc", "brier", "scaled_brier", "oeratio")
)


# validation under interventions metrics ----------------------------------


km.cens.stand <- survfit(Surv(T.obs,D.obs==0)~1,data=dat.val)
km.cens.step <- stepfun(km.cens.stand$time,c(1,km.cens.stand$surv))

#-----------------
#-----------------
# 'never treated' strategy
#-----------------
#-----------------

#-----------------------
# artificially censor the event time and event status in dat.val

dat.val$T.cens.0 <- ifelse(dat.val$A.0==1,0,
                           ifelse(dat.val$A.1==1,1,
                                  ifelse(dat.val$A.2==1, 2,
                                         ifelse(dat.val$A.3==1, 3,
                                                ifelse(dat.val$A.4==1,4,dat.val$T.obs)))))
dat.val$D.cens.0 <- ifelse(dat.val$A.0==1,0,
                           ifelse(dat.val$A.1==1,0,
                                  ifelse(dat.val$A.2==1, 0,
                                         ifelse(dat.val$A.3==1, 0,
                                                ifelse(dat.val$A.4==1,0,dat.val$D.obs)))))

#-----------------------
# time-dependent weights for artificial censoring (ipw0 and ipw1) have been made in 'analysis_x_model_validation.R' at visit time points (and at T.obs)
# for cindex and auct we need for each id weights at each event time that occurred in subjects who are not artificially censored

# unique event time points for 'never treated'
event.times.0 <- sort(unique(dat.long.val$time.stop[dat.long.val$in.dat.0==1 & dat.long.val$event==1]))
# we add the evaluation time point t=5 (we actually use 4.999) as we need it for calculation of CD AUCt and Brier score (at t=5)
event.times.0<-c(event.times.0,4.9999)

# select the rows from dat.long.val relevant for never treated scenario (using indicator in.dat.0)
# split on event times relevant for the 'never treated' scenario (note that data are already split on visit times)
# note this step takes few seconds computation time (dependent on n)
dat.0.split <- survSplit( Surv(time,time.stop,event) ~., data = dat.long.val %>% filter(in.dat.0==1), cut = event.times.0)

# calculate weights for standard censoring at the end time points
dat.0.split$ipw.othercens <- 1/km.cens.step(dat.0.split$time.stop)
# combine these weights with the weights for artificial censoring
dat.0.split$ipw.comb <- dat.0.split$ipw0 * dat.0.split$ipw.othercens

#--------------
# construct weights matrix needed for cindex / auct 'never treated'
# select the weights at event time points + make sure subjects who are censored before the first event time point are kept in the dataset
# put weights in wide format
dat.0.wide<- dat.0.split %>%
  filter(time.stop %in% event.times.0 | (time.stop < min(event.times.0))) %>%
  select(c("id","time.stop","ipw.comb")) %>%
  spread(time.stop, ipw.comb)

# subjects with censoring time before first event time should not add a row but not a column to the weightsmatrix (only columns for event time points are needed)
n.cens.before.first.event.0 <- sum(dat.0.split$time.stop < min(event.times.0))
dat.0.wide <- as.matrix(dat.0.wide[,-1:-(1+n.cens.before.first.event.0)]) #one additional column is deleted ("id")

# the above weightsmatrix only contains with at least A0=0. Expand it so that it has rows of NA for people in who are directly censored at t=0
wt_matrix0 <- matrix(nrow=n,ncol=ncol(dat.0.wide))
ids.0<-unique(dat.long.val$id[dat.long.val$A==0])
wt_matrix0[ids.0,] <- dat.0.wide

# remove last column from weightmatrix (weights at 4.9999) as this is not needed in the cindex calculation
wt_matrix0_eventsonly <- wt_matrix0[,-ncol(wt_matrix0)]

#------------
# AUCt for the 'never treated' strategy
#------------
i <- 1 # 'simulation iteration'
#'real' C/D AUCt (using counterfactual data for the never treated strategy)
rCDauct0 <- wCD_AUCt(time=dat.cf$T.A0, status=dat.cf$D.A0, risk=risk0_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,1] <- rCDauct0$AUCt$AUC[rCDauct0$AUCt$time==4.9999]

#weighted C/D AUCt from artificially censored validation data for the never treated strategy (our proposal)
wCDauct0 <- wCD_AUCt(time=dat.val$T.cens.0, status=dat.val$D.cens.0, risk=risk0_exp[,5], plot = FALSE, weightmatrix = wt_matrix0, seq.time=4.9999)
disc_auct[i,4] <- wCDauct0$AUCt$AUC[wCDauct0$AUCt$time==4.9999]


expect_equal(disc_auct[i,4], ipscore_results_0$score$auc[[2]])



# Brier -------------------------------------------------------------------

#'real' Brier and scaled Brier (IPA) using counterfactual data for the never treated strategy
brier_raw[i,1] <- Brier(dat.cf$T.A0, dat.cf$D.A0, risk0_exp[,5], seq.time=4.9999)
brier_ipa[i,1] <- ipa(dat.cf$T.A0, dat.cf$D.A0, risk0_exp[,5], seq.time=4.9999)

#td-weighted Brier and scaled Brier (IPA) from artificial censored data for the never treated strategy (our proposal)
# first derive a vector of weights
# for individuals with event after art cens (T.cens.0==1) - weight at their own event time
# for individuals known without event by tau after art cens (T.cens.0==5) - weight at tau
# weights for others are not needed (also does not mind if NA or other number)
# taking the last ipw0 weight recorded in dat.long.val for each id suffices here as max time is tau in the data
# calculate last ipw0 value in dat.long.val
dat.long.val <- dat.long.val %>% group_by(id) %>% mutate(ipw0.last=last(ipw0))
#  add this value to dat.val
dat.val <- left_join(dat.val, dat.long.val %>%
                       select(c("id","ipw0.last")) %>%
                       filter(row_number()==1),
                     by="id")

# combine with standard censoring
dat.val$ipw.othercens <- 1/km.cens.step(dat.val$T.obs)
dat.val$ipw.othercens[dat.val$T.obs==5] <- 1/km.cens.step(4.9999)
dat.val$ipw0.comb <- dat.val$ipw0.last * dat.val$ipw.othercens

brier_raw[i,3] <- Brier(dat.val$T.cens.0, dat.val$D.cens.0, risk=risk0_exp[,5], seq.time=4.9999, weights=dat.val$ipw0.comb)
brier_ipa[i,3] <- ipa(dat.val$T.cens.0, dat.val$D.cens.0, risk=risk0_exp[,5], seq.time=4.9999, weights=dat.val$ipw0.comb)

brier_raw
expect_equal(brier_raw[i,3], ipscore_results_0$score$brier[[2]])
expect_equal(brier_ipa[i,3], ipscore_results_0$score$scaled_brier[[2]])

brier_ipa[i, 1]


mean(dat.cf$D.A0)
ipscore_results_0$predictions$km0[[1]]
ipscore_results_0$predictions$`null model`[[1]]

# calibration -------------------------------------------------------------

#------------------------------------
#------------------------------------
#ESTIMATED mean risks at times 1:5 under the two treatment strategies: "always treated" (risk1), "never treated" (risk0)
#------------------------------------
#------------------------------------

calib_risk0[i,,1]=sapply(1:5,FUN=function(x){mean(risk0_exp[,x])})
calib_risk1[i,,1]=sapply(1:5,FUN=function(x){mean(risk1_exp[,x])})

#------------------------------------
#------------------------------------
#OBSERVED risks at times 1:5 under the two treatment strategies: "always treated" (risk1), "never treated" (risk0)
#------------------------------------
#------------------------------------

calib_risk0[i,,2]=risk0_obs
calib_risk1[i,,2]=risk1_obs

ratio <- calib_risk0[,5,2]/calib_risk0[,5,1]

expect_equal(ratio, ipscore_results_0$score$oeratio[[1]])


