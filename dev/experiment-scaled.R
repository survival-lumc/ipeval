library(dplyr)
library(tidyr)
library(survival)


n <- 3000

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

n.visit=5

simulate <- function(seed) {
  set.seed(seed)

  U <- rnorm(n,0,0.1)

  dat_sim_cox_scenario1 <- function() {


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
      repeat {
        u.t=runif(n,0,1)
        haz=exp(alpha.0+alpha.A*A[,k]+alpha.L*L[,k]+alpha.U*U)
        new.t=-log(u.t)/haz
        if (all(new.t >= 0.0000001)) {
          break
        } else {
          print("too small t")
        }
      }

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


  # eval --------------------------------------------------------------------


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


  true_results_0 <- observed_score(list(rep(mean(dat.cf$D.A0), n), risk0_exp[,5]),
                                   dat.cf, outcome = D.A0,
                                   metrics = c("auc", "brier", "scaled_brier",
                                               "oeratio"))

  results <- c(
    "auc_cf" = ipscore_results_0$score$auc[[3]],
    "auc_tr" = true_results_0$score$auc[[2]],

    "brier_cf" = ipscore_results_0$score$brier[[3]],
    "brier_tr" = true_results_0$score$brier[[2]],

    "scaled_cf" = ipscore_results_0$score$scaled_brier[[3]],
    "scaled_tr" = true_results_0$score$scaled_brier[[2]],

    "oe_cf" = ipscore_results_0$score$oeratio[[3]],
    "oe_tr" = true_results_0$score$oeratio[[2]],

    "null_cf" = ipscore_results_0$predictions$`null model`[[1]],
    "null_cf_km" = ipscore_results_0$predictions$km0[[1]],
    "null_tr" = mean(dat.cf$D.A0),

    "briernull_cf" = ipscore_results_0$score$brier[[2]],
    "briernull_tr" = true_results_0$score$brier[[1]]
  )

  return(results)
}

results <- lapply_progress(as.list(1:1000), simulate, "")
results <- as.data.frame(do.call(rbind, results))

results$auc_bias <- results$auc_cf - results$auc_tr
results$brier_bias <- results$brier_cf - results$brier_tr
results$scaled_bias <- results$scaled_cf - results$scaled_tr
results$oe_bias <- results$oe_cf - results$oe_tr
results$briernull_bias <- results$briernull_cf - results$briernull_tr

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))


mean(results$auc_bias) # -8.807769e-05
se(results$auc_bias) # 0.0009173214
mean(results$brier_bias) # 7.329222e-05
se(results$brier_bias) # 0.00026801
mean(results$scaled_bias) # -0.1207792
se(results$scaled_bias) # 0.08484403
mean(results$oe_bias) # -0.0009896135
se(results$oe_bias) # 0.001189472
mean(results$briernull_bias) # -0.0002612435
se(results$briernull_bias) # 0.0001291938

t.test(results$auc_bias) #  p-value = 0.9235
t.test(results$brier_bias) # p-value = 0.7845
t.test(results$scaled_bias) # p-value = 0.1549
t.test(results$oe_bias) # p-value = 0.4056
t.test(results$briernull_bias) # p-value = 0.04343
