#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#Using TMLE as implemented in the lmtp package
#######################################################

#------------------------------------------------------
#DATA SET-UP
#------------------------------------------------------

#identify people meeting eligibility criteria: people who have started a second line treatment
dta.analysis<-dta
dta.analysis$treat_status_any<-ifelse((dta.analysis$treat_status_A+dta.analysis$treat_status_B+
                                         dta.analysis$treat_status_C+dta.analysis$treat_status_D)>=1,1,0)
dta.analysis<-dta.analysis[dta.analysis$treat_status_any==1,]

#generate rownumber
dta.analysis<-dta.analysis%>%group_by(id)%>%mutate(rownum=row_number())

#Restrict to the first row
dta.analysis<-dta.analysis[dta.analysis$rownum==1,]

#restrict to people with HbA1c>=7.5 on starting treatment
dta.analysis<-dta.analysis[dta.analysis$hb>=7.5,]

#generate composite event
dta.analysis$event_composite<-ifelse(dta.analysis$event_type==1|dta.analysis$event_type==2,1,0)

#generate new event time that is relative to the time of starting treatment
dta.analysis$event_time_new<-dta.analysis$event_time-dta.analysis$tstart

#generate categorical treatment type variable
dta.analysis$treat_type<-""
dta.analysis$treat_type<-ifelse(dta.analysis$treat_status_A==1,"A",dta.analysis$treat_type)
dta.analysis$treat_type<-ifelse(dta.analysis$treat_status_B==1,"B",dta.analysis$treat_type)
dta.analysis$treat_type<-ifelse(dta.analysis$treat_status_C==1,"C",dta.analysis$treat_type)
dta.analysis$treat_type<-ifelse(dta.analysis$treat_status_D==1,"D",dta.analysis$treat_type)

#------------------------------------------------------
#set data up in the discrete-time wide format needed for lmtp
#------------------------------------------------------

dta.tmle<-dta.analysis[,c("id","event_time_new","event_type","treat_type",
                          "sex","age","smok_former","smok_current",
                          "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc")]
dta.tmle$event_type<-ifelse(dta.tmle$event_type==3,0,dta.tmle$event_type)

dta.tmle.split<-survSplit(Surv(event_time_new/365.25,as.factor(event_type))~.,data=dta.tmle,cut=seq(0,5,0.1))
dta.tmle.split$event<-ifelse(dta.tmle.split$event=="censor",0,dta.tmle.split$event)
dta.tmle.split$event<-ifelse(dta.tmle.split$event==2,1,dta.tmle.split$event)
dta.tmle.split$event<-ifelse(dta.tmle.split$event==3,2,dta.tmle.split$event)

dta.tmle.split$Y<-ifelse(dta.tmle.split$event==1,1,0)
dta.tmle.split$D<-ifelse(dta.tmle.split$event==2,1,0)
dta.tmle.split$C<-ifelse(dta.tmle.split$event==0,1,0)

dta.tmle.split$tstart<-NULL
dta.tmle.split$tstop<-NULL
dta.tmle.split$event<-NULL

dta.tmle.split<-dta.tmle.split%>%group_by(id)%>%mutate(time=row_number())

dta.tmle.split<-as.data.frame(dta.tmle.split)

dta.wide <- pivot_wider(dta.tmle.split, names_from = "time", 
                         values_from = c("C","Y","D"),names_sep=".",names_sort=F)

maxtime<-max(dta.tmle.split$time)
dta.wide<-dta.wide[,c("id","treat_type","sex","age","smok_former","smok_current",
                      "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc",
                      paste(rep(c("C","D","Y"),maxtime),rep(c(1:maxtime),each=3),sep="."))]

paste(rep(c("C","D","Y"),maxtime),rep(c(1:maxtime),each=3),sep=".")

#if Y_k=1 or D_k=1 then C_[k+1]=NA, for all k
#This is already the case

#if Y_k=1 or D_k=1 then C_[k]=1, for all k
for(k in 1:(maxtime-1)){
  eval(parse(text=paste0("dta.wide$C.",k,"=ifelse(dta.wide$Y.",k,"==1 & dta.wide$D.",k,"==0 & 
  !is.na(dta.wide$Y.",k,")  & !is.na(dta.wide$D.",k,"),
                         1,dta.wide$C.",k,")")))
  
  eval(parse(text=paste0("dta.wide$C.",k,"=ifelse(dta.wide$Y.",k,"==0 & dta.wide$D.",k,"==1 & 
  !is.na(dta.wide$Y.",k,")  & !is.na(dta.wide$D.",k,"),
                         1,dta.wide$C.",k,")")))
}

#if Y_k=1 and D_k=0 then Y_[k+1]=1 and D_[k+1]=0, for all k
#if Y_k=0 and D_k=1 then Y_[k+1]=0 and D_[k+1]=1, for all k
for(k in 1:(maxtime-1)){
  eval(parse(text=paste0("dta.wide$Y.",k+1,"=ifelse(dta.wide$Y.",k,"==1 & dta.wide$D.",k,"==0 & 
  !is.na(dta.wide$Y.",k,")  & !is.na(dta.wide$D.",k,"),
                         dta.wide$Y.",k,",dta.wide$Y.",k+1,")")))
  
  
  eval(parse(text=paste0("dta.wide$D.",k+1,"=ifelse(dta.wide$Y.",k,"==1 & dta.wide$D.",k,"==0 & 
  !is.na(dta.wide$Y.",k,")  & !is.na(dta.wide$D.",k,"),
                         dta.wide$D.",k,",dta.wide$D.",k+1,")")))
  
  eval(parse(text=paste0("dta.wide$Y.",k+1,"=ifelse(dta.wide$Y.",k,"==0 & dta.wide$D.",k,"==1 &
  !is.na(dta.wide$Y.",k,")  & !is.na(dta.wide$D.",k,"),
                         dta.wide$Y.",k,",dta.wide$Y.",k+1,")")))
  
  eval(parse(text=paste0("dta.wide$D.",k+1,"=ifelse(dta.wide$Y.",k,"==0 & dta.wide$D.",k,"==1 &
  !is.na(dta.wide$Y.",k,")  & !is.na(dta.wide$D.",k,"),
                         dta.wide$D.",k,",dta.wide$D.",k+1,")")))
}

#if (C_k=1, Y_k=0, D_k=0) and C_[k+1]=NA then set C_k=0, for all k
for(k in 1:(maxtime-1)){
  eval(parse(text=paste0("dta.wide$C.",k,"=ifelse(dta.wide$Y.",k,"==0 & dta.wide$D.",k,"==0 & dta.wide$C.",k,"==1 &
  !is.na(dta.wide$Y.",k,")  & !is.na(dta.wide$D.",k,")  & is.na(dta.wide$C.",k+1,"),
                         0,dta.wide$C.",k,")")))
}

#if C_k=0 (i.e censored) then set Y_k=NA and D_k=NA, for all k
for(k in 1:maxtime){
  eval(parse(text=paste0("dta.wide$Y.",k,"=ifelse(dta.wide$C.",k,"==0 & !is.na(dta.wide$C.",k,"),
                         NA,dta.wide$Y.",k,")")))
  
  eval(parse(text=paste0("dta.wide$D.",k,"=ifelse(dta.wide$C.",k,"==0 & !is.na(dta.wide$C.",k,"),
                         NA,dta.wide$D.",k,")")))
}

#if C_k=0 (i.e censored) then set C_[k+1]=0, for all k
for(k in 1:(maxtime-1)){
  eval(parse(text=paste0("dta.wide$C.",k+1,"=ifelse(dta.wide$C.",k,"==0 & !is.na(dta.wide$C.",k,"),
                         0,dta.wide$C.",k+1,")")))
}

#------------------------------------------------------
#Using LMTP
#------------------------------------------------------

dta.wide<-dta.wide[,c("id","treat_type","sex","age","smok_former","smok_current",
                           "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc",
                           paste(rep(c("C","D","Y"),maxtime),rep(c(1:maxtime),each=3),sep="."))]

policy.A <- function(data, trt) {
     a <- data[[trt]]
     "A"}

policy.B <- function(data, trt) {
  a <- data[[trt]]
  "B"}

policy.C <- function(data, trt) {
  a <- data[[trt]]
  "C"}

policy.D <- function(data, trt) {
  a <- data[[trt]]
  "D"}

lmtp.A.mace<-lmtp_survival(
  data = dta.wide,
  trt = "treat_type",
  cens = paste0("C.", 1:maxtime),
  compete = paste0("D.", 1:maxtime),
  baseline = c("sex","age","smok_former","smok_current",
               "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
  outcome = paste0("Y.", 1:maxtime),
  shift = policy.A,
  folds = 1,
  estimator = "lmtp_tmle"
)

lmtp.B.mace<-lmtp_survival(
  data = dta.wide,
  trt = "treat_type",
  cens = paste0("C.", 1:maxtime),
  compete = paste0("D.", 1:maxtime),
  baseline = c("sex","age","smok_former","smok_current",
               "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
  outcome = paste0("Y.", 1:maxtime),
  shift = policy.B,
  folds = 1,
  estimator = "lmtp_tmle"
)

lmtp.C.mace<-lmtp_survival(
  data = dta.wide,
  trt = "treat_type",
  cens = paste0("C.", 1:maxtime),
  compete = paste0("D.", 1:maxtime),
  baseline = c("sex","age","smok_former","smok_current",
               "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
  outcome = paste0("Y.", 1:maxtime),
  shift = policy.C,
  folds = 1,
  estimator = "lmtp_tmle"
)

lmtp.D.mace<-lmtp_survival(
  data = dta.wide,
  trt = "treat_type",
  cens = paste0("C.", 1:maxtime),
  compete = paste0("D.", 1:maxtime),
  baseline = c("sex","age","smok_former","smok_current",
               "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
  outcome = paste0("Y.", 1:maxtime),
  shift = policy.D,
  folds = 1,
  estimator = "lmtp_tmle"
)

lmtp.A.death<-lmtp_survival(
  data = dta.wide,
  trt = "treat_type",
  cens = paste0("C.", 1:maxtime),
  compete = paste0("Y.", 1:maxtime),
  baseline = c("sex","age","smok_former","smok_current",
               "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
  outcome = paste0("D.", 1:maxtime),
  shift = policy.A,
  folds = 1,
  estimator = "lmtp_tmle"
)

lmtp.B.death<-lmtp_survival(
  data = dta.wide,
  trt = "treat_type",
  cens = paste0("C.", 1:maxtime),
  compete = paste0("Y.", 1:maxtime),
  baseline = c("sex","age","smok_former","smok_current",
               "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
  outcome = paste0("D.", 1:maxtime),
  shift = policy.B,
  folds = 1,
  estimator = "lmtp_tmle"
)

lmtp.C.death<-lmtp_survival(
  data = dta.wide,
  trt = "treat_type",
  cens = paste0("C.", 1:maxtime),
  compete = paste0("Y.", 1:maxtime),
  baseline = c("sex","age","smok_former","smok_current",
               "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
  outcome = paste0("D.", 1:maxtime),
  shift = policy.C,
  folds = 1,
  estimator = "lmtp_tmle"
)

lmtp.D.death<-lmtp_survival(
  data = dta.wide,
  trt = "treat_type",
  cens = paste0("C.", 1:maxtime),
  compete = paste0("Y.", 1:maxtime),
  baseline = c("sex","age","smok_former","smok_current",
               "diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
  outcome = paste0("D.", 1:maxtime),
  shift = policy.D,
  folds = 1,
  estimator = "lmtp_tmle"
)

#save estimates

tidylmtp.A.mace<-tidy(lmtp.A.mace)
tidylmtp.B.mace<-tidy(lmtp.B.mace)
tidylmtp.C.mace<-tidy(lmtp.C.mace)
tidylmtp.D.mace<-tidy(lmtp.D.mace)

tidylmtp.A.death<-tidy(lmtp.A.death)
tidylmtp.B.death<-tidy(lmtp.B.death)
tidylmtp.C.death<-tidy(lmtp.C.death)
tidylmtp.D.death<-tidy(lmtp.D.death)

cinc.A.mace.lmtp<-1-c(1,tidylmtp.A.mace$estimate)
cinc.B.mace.lmtp<-1-c(1,tidylmtp.B.mace$estimate)
cinc.C.mace.lmtp<-1-c(1,tidylmtp.C.mace$estimate)
cinc.D.mace.lmtp<-1-c(1,tidylmtp.D.mace$estimate)

cinc.A.death.lmtp<-1-c(1,tidylmtp.A.death$estimate)
cinc.B.death.lmtp<-1-c(1,tidylmtp.B.death$estimate)
cinc.C.death.lmtp<-1-c(1,tidylmtp.C.death$estimate)
cinc.D.death.lmtp<-1-c(1,tidylmtp.D.death$estimate)

#save 95% CIs

cinc.A.mace.lmtp.CIL<-1-c(1,tidylmtp.A.mace$conf.high)
cinc.B.mace.lmtp.CIL<-1-c(1,tidylmtp.B.mace$conf.high)
cinc.C.mace.lmtp.CIL<-1-c(1,tidylmtp.C.mace$conf.high)
cinc.D.mace.lmtp.CIL<-1-c(1,tidylmtp.D.mace$conf.high)

cinc.A.mace.lmtp.CIU<-1-c(1,tidylmtp.A.mace$conf.low)
cinc.B.mace.lmtp.CIU<-1-c(1,tidylmtp.B.mace$conf.low)
cinc.C.mace.lmtp.CIU<-1-c(1,tidylmtp.C.mace$conf.low)
cinc.D.mace.lmtp.CIU<-1-c(1,tidylmtp.D.mace$conf.low)

cinc.A.death.lmtp.CIL<-1-c(1,tidylmtp.A.death$conf.high)
cinc.B.death.lmtp.CIL<-1-c(1,tidylmtp.B.death$conf.high)
cinc.C.death.lmtp.CIL<-1-c(1,tidylmtp.C.death$conf.high)
cinc.D.death.lmtp.CIL<-1-c(1,tidylmtp.D.death$conf.high)

cinc.A.death.lmtp.CIU<-1-c(1,tidylmtp.A.death$conf.low)
cinc.B.death.lmtp.CIU<-1-c(1,tidylmtp.B.death$conf.low)
cinc.C.death.lmtp.CIU<-1-c(1,tidylmtp.C.death$conf.low)
cinc.D.death.lmtp.CIU<-1-c(1,tidylmtp.D.death$conf.low)




