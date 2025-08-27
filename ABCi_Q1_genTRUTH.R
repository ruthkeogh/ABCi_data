###############################################################
#Generate counterfactual data for Q1. 
#This file generates data as if everyone starting treatment were to start treatment A, B, C or D.
#The treatment is specified by trt.fix in analysis_Q1_TRUTH.R
#
#Data up to treatment initiation are generated in the same way as in the 'observed' data (ABCi_data_generate.R).
#From the point of treatment initiation, the impacts of treatment on other variables are set to be those of the  specified treatment trt.fix.
#This is done by setting impacts of treatment on hazards for mace/death and on the time-dependent covariates to be as for trt.fix.
###############################################################

#----------
#set seed
#----------

set.seed(1)

#----------
#sample size
#----------
n<-100000

#----------
#set up data frame
#----------

#administrative censoring time: 1826.25 (5 years)
admin.cens.time<-1826.25

#list of times from 0 to 1826.25 days (5 years) in increments of 30 days
time_increment<-30 #could be modified
time_increment_obs<-30

times<-seq(0,admin.cens.time,time_increment)
ntimes<-length(times)

obs.times<-seq(0,admin.cens.time,time_increment_obs)

tstart<-rep(times,n)#start times for each row
tstop<-rep(times+time_increment,n)#stop times for each row
tstop<-ifelse(tstop>admin.cens.time,admin.cens.time,tstop)#change the last stop time to admin.cens.time

dta<-data.frame(id=rep(1:n,each=ntimes),tstart=tstart,tstop=tstop,sex=NA,age=NA,
                smok=NA,diabdur=NA,bmi=NA,hb=NA,hyp=NA,dys=NA,cvd=NA,kidney=NA,panc=NA,
                treat_time_A=NA,treat_time_B=NA,treat_time_C=NA,treat_time_D=NA,
                treat_status_A=NA,treat_status_B=NA,treat_status_C=NA,treat_status_D=NA,
                mace_time=NA,death_time=NA,cens_time=NA,
                mace_status=NA,death_status=NA,cens_status=NA)

#--------------------------------------------
#impacts of treatment on different covariates
#Original values used in the 'observed' data generating mechanism
#--------------------------------------------

#hazard of MACE: main effects of treatment (log HRs)
ORIGcoef.mace.A<-log(1.1)
ORIGcoef.mace.B<-log(1.1*0.8)
ORIGcoef.mace.C<-log(1.1*0.3)
ORIGcoef.mace.D<-log(1.1*0.8)

#on hazard of MACE: interactions with cvd and kidney disease
ORIGcoef.mace.A.cvd<-0
ORIGcoef.mace.B.cvd<-log(0.7)
ORIGcoef.mace.C.cvd<-0
ORIGcoef.mace.D.cvd<-0

ORIGcoef.mace.A.kidney<-0
ORIGcoef.mace.B.kidney<-0
ORIGcoef.mace.C.kidney<-log(0.8)
ORIGcoef.mace.D.kidney<-0

#hazard of death: main effects of treatment (log HRs)
ORIGcoef.death.A<-log(1.1)
ORIGcoef.death.B<-log(1.1*0.65)
ORIGcoef.death.C<-log(1.1*0.5)
ORIGcoef.death.D<-log(1.1*0.75)

#bmi: the effect of treatment is a function of time over the first year after treatment initiation 
#and then remains constant after 1 year
ORIGbmi.coef.A.yr1<-0
ORIGbmi.coef.A.longterm<-0
ORIGbmi.coef.B.yr1<--1.6
ORIGbmi.coef.B.longterm<--1.6
ORIGbmi.coef.C.yr1<--2.6
ORIGbmi.coef.C.longterm<--2.6
ORIGbmi.coef.D.yr1<--0.7
ORIGbmi.coef.D.longterm<--0.7

#ORIGhba1c: the effect of treatment is a function of time over the first 2 years after treatment initiation 
#and then remains constant after year 2
ORIGhb.slope.A.yr1<--2.68
ORIGhb.const.A.yr1to2<--2.68
ORIGhb.slope.A.yr1to2<-0
ORIGhb.const.A.yr2plus<--2.68

ORIGhb.slope.B.yr1<--(2.68+0.23)
ORIGhb.const.B.yr1to2<--(2.68+0.23)
ORIGhb.slope.B.yr1to2<--(2.68+0.23+0.24)
ORIGhb.const.B.yr2plus<--(2.68+0.23+0.24)

ORIGhb.slope.C.yr1<--1
ORIGhb.const.C.yr1to2<--1
ORIGhb.slope.C.yr1to2<-0
ORIGhb.const.C.yr2plus<--1

ORIGhb.slope.D.yr1<-0
ORIGhb.const.D.yr1to2<-0
ORIGhb.slope.D.yr1to2<--(2.68+0.13)
ORIGhb.const.D.yr2plus<--(2.68+0.13)

#ORIGkidney disease: log (ORs) for ORIGkidney disease
ORIGkid.coef.A<-log(1.2)
ORIGkid.coef.B<-log(0.62)
ORIGkid.coef.C<-log(0.86)
ORIGkid.coef.D<-(log(1.2)+log(0.72))

#--------------------------------------------
#parameters: impacts of treatments on different covariates
#set these to be those of the treatment of interest (trt.fix)
#--------------------------------------------

#hazard of MACE: main effects of treatment (log HRs)
eval(parse(text=paste0("coef.mace.A<-coef.mace.B<-coef.mace.C<-coef.mace.D<-ORIGcoef.mace.",trt.fix)))

#on hazard of MACE: interactions with cvd and kidney disease
eval(parse(text=paste0("coef.mace.A.cvd<-coef.mace.B.cvd<-coef.mace.C.cvd<-coef.mace.D.cvd<-ORIGcoef.mace.",trt.fix,".cvd")))
eval(parse(text=paste0("coef.mace.A.kidney<-coef.mace.B.kidney<-coef.mace.C.kidney<-coef.mace.D.kidney<-ORIGcoef.mace.",trt.fix,".kidney")))

#hazard of death: main effects of treatment (log HRs)
eval(parse(text=paste0("coef.death.A<-coef.death.B<-coef.death.C<-coef.death.D<-ORIGcoef.death.",trt.fix)))

#bmi
eval(parse(text=paste0("bmi.coef.A.yr1<-bmi.coef.B.yr1<-bmi.coef.C.yr1<-bmi.coef.D.yr1<-ORIGbmi.coef.",trt.fix,".yr1")))
eval(parse(text=paste0("bmi.coef.A.longterm<-bmi.coef.B.longterm<-bmi.coef.C.longterm<-bmi.coef.D.longterm<-ORIGbmi.coef.",trt.fix,".longterm")))

#hba1c
eval(parse(text=paste0("hb.slope.A.yr1<-hb.slope.B.yr1<-hb.slope.C.yr1<-hb.slope.D.yr1<-ORIGhb.slope.",trt.fix,".yr1")))
eval(parse(text=paste0("hb.slope.A.yr1to2<-hb.slope.B.yr1to2<-hb.slope.C.yr1to2<-hb.slope.D.yr1to2<-ORIGhb.slope.",trt.fix,".yr1to2")))

eval(parse(text=paste0("hb.const.A.yr1to2<-hb.const.B.yr1to2<-hb.const.C.yr1to2<-hb.const.D.yr1to2<-ORIGhb.const.",trt.fix,".yr1to2")))
eval(parse(text=paste0("hb.const.A.yr2plus<-hb.const.B.yr2plus<-hb.const.C.yr2plus<-hb.const.D.yr2plus<-ORIGhb.const.",trt.fix,".yr2plus")))

#--------------------------------------------
#--------------------------------------------
#FROM THIS POINT ONWARDS THE DATA GENERATION CODE IS AS IN ABCi_data_generate.R
#--------------------------------------------
#--------------------------------------------

#--------------------------------------------
#random intercepts and slopes for BMI and HbA1C 
#--------------------------------------------

bmi_random_intercept<-rnorm(n,0,6)
bmi_random_slope<-rnorm(n,0,0.3)

hb_random_intercept<-rnorm(n,0,2)
hb_random_slope<-rnorm(n,0,0.1)

#--------------------------------------------
#functions: hazards for time to initiation of treatments A-D
#--------------------------------------------

haz.A<-function(sex,age,bmi,diabdur,hb,hb_change,smok_former,smok_current,cvd,kidney){
  exp(-9-log(0.83)*sex+log(0.98)*age+log(1.01)*bmi+
        log(0.94)*diabdur+log(5)*(hb>7.5)+log(1.7)*(hb-7.5)*(hb>7.5)+log(2)*hb_change+
        log(1.06)*smok_former+log(1.03)*smok_current)
}

haz.B<-function(sex,age,bmi,diabdur,hb,hb_change,smok_former,smok_current,cvd,kidney){
  exp(log(haz.A(sex,age,bmi,diabdur,hb,hb_change,smok_former,smok_current))+
        log(1.42)*(age<75)*(age>=65)+
        log(1.26)*(age>=75)-
        0.03*hb*(hb>8)*(hb<=9)+
        log(0.68)*(hb>9)+
        log(0.65)*cvd+log(0.28)*kidney)
}

haz.C<-function(sex,age,bmi,diabdur,hb,hb_change,smok_former,smok_current,cvd,kidney){
  exp(log(haz.A(sex,age,bmi,diabdur,hb,hb_change,smok_former,smok_current))-1+
        log(2.14)*(bmi-30)/5+
        log(0.44)*cvd+log(0.36)*kidney)
}

haz.D<-function(sex,age,bmi,diabdur,hb,hb_change,smok_former,smok_current,cvd,kidney){
  exp(log(haz.A(sex,age,bmi,diabdur,hb,hb_change,smok_former,smok_current,cvd,kidney))-1.5+
        log(1.31)*(age<75)*(age>=65)+
        log(1.35)*(age>=75)-
        log(1.03)*hb*(hb>8))
}

#--------------------------------------------
#functions: hazards for mace, death, censoring
#--------------------------------------------

haz.mace<-function(age,sex,smok_former,smok_current,bmi,hb,diabdur,dys,hyp,cvd,kidney,
                   treat_status_A,treat_status_B,treat_status_C,treat_status_D){
  exp(-11.5+log(1.14)*(age-60)+
        log(1.3)*(sex==1)+
        log(1.05)*smok_former+log(1.52)*smok_current+
        log(1.15)*(bmi>=30)+log(1.3)*(bmi>=40)+log(1.1)*(bmi<20)+
        log(1.3)*(hb-7)+
        log(1.05)*diabdur+
        log(2)*dys*(sex==0)+log(1.2)*dys*(sex==1)+
        log(1.7)*hyp+
        log(2)*kidney+
        log(3)*cvd+
        coef.mace.A*treat_status_A+
        coef.mace.B*treat_status_B+
        coef.mace.C*treat_status_C+
        coef.mace.D*treat_status_D+
        #treatment*cvd interactions
        coef.mace.A.cvd*treat_status_A*cvd+
        coef.mace.B.cvd*treat_status_B*cvd+
        coef.mace.C.cvd*treat_status_C*cvd+
        coef.mace.D.cvd*treat_status_D*cvd+
        #treatment*kidney interactions
        coef.mace.A.kidney*treat_status_A*kidney+
        coef.mace.B.kidney*treat_status_B*kidney+
        coef.mace.C.kidney*treat_status_C*kidney+
        coef.mace.D.kidney*treat_status_D*kidney)
}

haz.death<-function(age,sex,smok_former,smok_current,bmi,hb,diabdur,dys,hyp,cvd,kidney,
                    treat_status_A,treat_status_B,treat_status_C,treat_status_D){
  exp(-12.5+log(1.07)*(age-60)+
        log(0.83)*(sex==1)+
        log(1.07)*smok_former+log(1.62)*smok_current+
        (-0.334)*(bmi-30)+0.0054*((bmi-30)^2)+
        log(1)*hb+
        0.05*diabdur+
        log(0.9)*dys+
        log(1)*hyp+
        log(1.49)*kidney+
        log(1.25)*cvd+
        coef.death.A*treat_status_A+
        coef.death.B*treat_status_B+
        coef.death.C*treat_status_C+
        coef.death.D*treat_status_D)
}

haz.cens<-0.0001

#--------------------------------------------
#start of data generation loop
#--------------------------------------------

for(j in 1:ntimes)
{
  print(j)
  
  #--------------------------------------------
  #covariates: time-fixed and deterministically changing variables
  #sex, age, smok, diabdur
  #--------------------------------------------
  
  if(j==1){
    #sex
    dta$sex[dta$tstart==times[j]]<-rbinom(n,1,0.45)
    
    #age
    dta$age[dta$tstart==times[j]]<-rtruncnorm(n,a=0, b=95,mean=62,sd=12)
    dta$age_group[dta$tstart==times[j]]<-cut(dta$age[dta$tstart==times[j]],breaks=c(0,60,70,80,200),label=F)
    
    #smoking
    a<-exp(-1+0.5*(dta$age_group[dta$tstart==times[j]]==2)+1*(dta$age_group[dta$tstart==times[j]]==3|dta$age_group[dta$tstart==times[j]]==4))
    b<-exp(-0.5+0.8*(dta$age_group[dta$tstart==times[j]]==2)+1.2*(dta$age_group[dta$tstart==times[j]]==3|dta$age_group[dta$tstart==times[j]]==4))
    
    p_smok_current<-a/(1+a+b)
    p_smok_former<-b/(1+a+b)
    p_smok_never=1-p_smok_current-p_smok_former
    
    dta$smok[dta$tstart==times[j]]<-sapply(1:n,FUN=function(x){
      which(rmultinom(1,1,c(p_smok_never[x],p_smok_former[x],p_smok_current[x]))==1)})
    
    dta$smok_former[dta$tstart==times[j]]<-ifelse(dta$smok[dta$tstart==times[j]]==2,1,0)
    dta$smok_current[dta$tstart==times[j]]<-ifelse(dta$smok[dta$tstart==times[j]]==3,1,0)
    
    #time since diabetes diagnosis (years)
    dta$diabdur[dta$tstart==times[j]]<-rtruncnorm(n, a=0, b=Inf, mean = 0.3-0.05*(dta$age[dta$tstart==times[j]]-60), sd = 0.3)
  }else if(j>1){
    #sex
    dta$sex[dta$tstart==times[j]]<-dta$sex[dta$tstart==times[1]]
    
    #age
    dta$age[dta$tstart==times[j]]<-dta$age[dta$tstart==times[j-1]]+time_increment/365
    dta$age_group[dta$tstart==times[j]]<-cut(dta$age[dta$tstart==times[j]],breaks=c(0,60,70,80,200),label=F)
    
    #smoking
    dta$smok[dta$tstart==times[j]]<-dta$smok[dta$tstart==times[1]]
    dta$smok_former[dta$tstart==times[j]]<-dta$smok_former[dta$tstart==times[1]]    
    dta$smok_current[dta$tstart==times[j]]<-dta$smok_current[dta$tstart==times[1]]
    
    #time since diabetes diagnosis (years)
    dta$diabdur[dta$tstart==times[j]]<-dta$diabdur[dta$tstart==times[j-1]]+time_increment/365
  }
  
  #--------------------------------------------
  #Time-dependent covariates: bmi and HbA1c (hb)
  #These depend on treatment after the first time period
  #--------------------------------------------
  
  if(j==1){
    dta$bmi[dta$tstart==times[j]]<-rnorm(n,(30+bmi_random_intercept)+
                                           0.01*(dta$age[dta$tstart==times[j]]-60)-
                                           0.03*(dta$sex[dta$tstart==times[j]]==1),
                                         1)
    
    
    #hbA1c
    dta$hb[dta$tstart==times[j]]<-rtruncnorm(n, a=7, b=Inf, 
                                             mean = (7.2+hb_random_intercept)+(0.1+hb_random_slope)*0+
                                               0.09*(dta$age[dta$tstart==times[j]]-60)+
                                               0.1*(dta$bmi[dta$tstart==times[j]]-30)+0.13*(dta$sex[dta$tstart==times[j]]==1)+
                                               0.1*(dta$diabdur[dta$tstart==times[j]]-0.5), 
                                             sd = 0.2)
    
    #HbA1c at 30 days (time_increment) before baseline - used to calculate hb_change
    dta$hb_prebaseline[dta$tstart==times[j]]<-rtruncnorm(n, a=7, b=Inf, 
                                                         mean = (7.2+hb_random_intercept)+(0.1+hb_random_slope)*(-time_increment)/365+
                                                           0.09*(dta$age[dta$tstart==times[j]]-time_increment/365-60)+
                                                           0.1*(dta$bmi[dta$tstart==times[j]]-30)+0.13*(dta$sex[dta$tstart==times[j]]==1)+
                                                           0.1*(dta$diabdur[dta$tstart==times[j]]-time_increment/365-0.5), 
                                                         sd = 0.2)
    
    dta$hb_change[dta$tstart==times[j]]<-dta$hb[dta$tstart==times[j]]-
      dta$hb_prebaseline[dta$tstart==times[j]]
  }else if(j>1){
    #bmi: depends on treatment and time since treatment initiation
    dta$bmi[dta$tstart==times[j]]<-rnorm(n,(30+bmi_random_intercept)+(1+bmi_random_slope)*(times[j]/365)-
                                           0.1*((times[j]/365)^2)+
                                           0.01*(dta$age[dta$tstart==times[j]]-60)-
                                           0.03*(dta$sex[dta$tstart==times[j]]==1)-
                                           0.03*(dta$sex[dta$tstart==times[j]]==1)*(times[j]/365),
                                         1)+
      bmi.coef.A.yr1*dta$ts_A[dta$tstart==times[j]]*(dta$ts_A[dta$tstart==times[j]]<1)+
      bmi.coef.A.longterm*(dta$ts_A[dta$tstart==times[j]]>=1)+
      bmi.coef.B.yr1*dta$ts_B[dta$tstart==times[j]]*(dta$ts_B[dta$tstart==times[j]]<1)+
      bmi.coef.B.longterm*(dta$ts_B[dta$tstart==times[j]]>=1)+
      bmi.coef.C.yr1*dta$ts_C[dta$tstart==times[j]]*(dta$ts_C[dta$tstart==times[j]]<1)+
      bmi.coef.C.longterm*(dta$ts_C[dta$tstart==times[j]]>=1)+
      bmi.coef.D.yr1*dta$ts_D[dta$tstart==times[j]]*(dta$ts_D[dta$tstart==times[j]]<1)+
      bmi.coef.D.longterm*(dta$ts_D[dta$tstart==times[j]]>=1)
    
    #hbA1c: depends on treatment and time since treatment initiation
    dta$hb[dta$tstart==times[j]]<-rtruncnorm(n, a=7, b=Inf, 
                                             mean = (7.2+hb_random_intercept)+(0.1+hb_random_slope)*times[j]/365+
                                               0.09*(dta$age[dta$tstart==times[1]]-60)+
                                               0.1*(dta$bmi[dta$tstart==times[j]]-30)+0.13*(dta$sex[dta$tstart==times[j]]==1)+
                                               0.1*(dta$diabdur[dta$tstart==times[j]]-0.5)+
                                               #treatment A
                                               hb.slope.A.yr1*dta$ts_A[dta$tstart==times[j]]*(dta$ts_A[dta$tstart==times[j]]<1)+
                                               hb.const.A.yr1to2*(dta$ts_A[dta$tstart==times[j]]>=1)*(dta$ts_A[dta$tstart==times[j]]<2)+
                                               hb.slope.A.yr1to2*(dta$ts_A[dta$tstart==times[j]]-1)*(dta$ts_A[dta$tstart==times[j]]>=1)*(dta$ts_A[dta$tstart==times[j]]<2)+
                                               hb.const.A.yr2plus*(dta$ts_A[dta$tstart==times[j]]>=2)+
                                               #treatment B
                                               hb.slope.B.yr1*dta$ts_B[dta$tstart==times[j]]*(dta$ts_B[dta$tstart==times[j]]<1)+
                                               hb.const.B.yr1to2*(dta$ts_B[dta$tstart==times[j]]>=1)*(dta$ts_B[dta$tstart==times[j]]<2)+
                                               hb.slope.B.yr1to2*(dta$ts_B[dta$tstart==times[j]]-1)*(dta$ts_B[dta$tstart==times[j]]>=1)*(dta$ts_B[dta$tstart==times[j]]<2)+
                                               hb.const.B.yr2plus*(dta$ts_B[dta$tstart==times[j]]>=2)+
                                               #treatment C
                                               hb.slope.C.yr1*dta$ts_C[dta$tstart==times[j]]*(dta$ts_C[dta$tstart==times[j]]<1)+
                                               hb.const.C.yr1to2*(dta$ts_C[dta$tstart==times[j]]>=1)*(dta$ts_C[dta$tstart==times[j]]<2)+
                                               hb.slope.C.yr1to2*(dta$ts_C[dta$tstart==times[j]]-1)*(dta$ts_C[dta$tstart==times[j]]>=1)*(dta$ts_C[dta$tstart==times[j]]<2)+
                                               hb.const.C.yr2plus*(dta$ts_C[dta$tstart==times[j]]>=2)+
                                               #treatment D
                                               hb.slope.D.yr1*dta$ts_D[dta$tstart==times[j]]*(dta$ts_D[dta$tstart==times[j]]<1)+
                                               hb.const.D.yr1to2*(dta$ts_D[dta$tstart==times[j]]>=1)*(dta$ts_D[dta$tstart==times[j]]<2)+
                                               hb.slope.D.yr1to2*(dta$ts_D[dta$tstart==times[j]]-1)*(dta$ts_D[dta$tstart==times[j]]>=1)*(dta$ts_D[dta$tstart==times[j]]<2)+
                                               hb.const.D.yr2plus*(dta$ts_D[dta$tstart==times[j]]>=2),
                                             sd = 0.2)
    
    #Change in HbA1c since last visit 
    dta$hb_change[dta$tstart==times[j]]<-dta$hb[dta$tstart==times[j]]-dta$hb[dta$tstart==times[j]-time_increment]
    
  }
  
  #--------------------------------------------
  #Time-dependent covariates: comorbidities (hyp, dys, cvd, kidney, panc)
  #These do not depend directly on treatment
  #--------------------------------------------
  if(j==1){
    #hypertension
    dta$hyp[dta$tstart==times[j]]<-rbinom(n,1,
                                          expit(-0.7+0.04*(dta$age[dta$tstart==times[j]]-60)+
                                                  0.0002*((dta$age[dta$tstart==times[j]]-60)^2)+
                                                  0.2*(dta$sex[dta$tstart==times[j]]==1)+
                                                  0.2*(dta$smok[dta$tstart==times[j]]==3)+
                                                  0.05*(dta$smok[dta$tstart==times[j]]==2)+
                                                  0.02*dta$bmi[dta$tstart==times[j]]+
                                                  0.05*dta$hb[dta$tstart==times[j]]))
    
    #dyslipidemia
    dta$dys[dta$tstart==times[j]]<-rbinom(n,1,
                                          expit(-0.5+0.05*(dta$age[dta$tstart==times[j]]-60)+
                                                  0.1*(dta$sex[dta$tstart==times[j]]==1)+
                                                  0.02*(dta$bmi[dta$tstart==times[j]]-30)))
    
    #history of cvd
    dta$cvd[dta$tstart==times[j]]<-rbinom(n,1,
                                          expit(-3.5+0.5*(dta$age_group[dta$tstart==times[j]]==2)+2*(dta$age_group[dta$tstart==times[j]]==3)+3*(dta$age_group[dta$tstart==times[j]]==4)+
                                                  0.05*(dta$bmi[dta$tstart==times[j]]-30)+
                                                  0.5*dta$dys[dta$tstart==times[j]]+1*dta$hyp[dta$tstart==times[j]]))
    
    #history of kidney disease
    dta$kidney[dta$tstart==times[j]]<-rbinom(n,1,
                                             expit(-5.5+0.5*(dta$age_group[dta$tstart==times[j]]==2)+1*(dta$age_group[dta$tstart==times[j]]==3)+1.5*(dta$age_group[dta$tstart==times[j]]==4)+
                                                     0.2*dta$hb[dta$tstart==times[j]]+
                                                     0.5*dta$hyp[dta$tstart==times[j]]))
    
    #history of pancreatitis
    dta$panc[dta$tstart==times[j]]<-rbinom(n,1,
                                           expit(-7.5+0.1*dta$bmi[dta$tstart==times[j]]+
                                                   2*(dta$smok[dta$tstart==times[j]]==3)+
                                                   0.5*(dta$smok[dta$tstart==times[j]]==2)))
    
  }else if(j>1){
    
    #hypertension
    dta$hyp[dta$tstart==times[j]]<-ifelse(dta$hyp[dta$tstart==times[j-1]]==1,1,
                                          rbinom(n,1,
                                                 expit(-7+0.04*(dta$age[dta$tstart==times[j]]-60)+
                                                         0.0002*((dta$age[dta$tstart==times[j]]-60)^2)+
                                                         0.2*(dta$sex[dta$tstart==times[j]]==1)+
                                                         0.2*(dta$smok[dta$tstart==times[j]]==3)+
                                                         0.05*(dta$smok[dta$tstart==times[j]]==2)+
                                                         0.02*dta$bmi[dta$tstart==times[j]]+
                                                         0.05*dta$hb[dta$tstart==times[j]])))
    
    #dyslipidemia
    dta$dys[dta$tstart==times[j]]<-ifelse(dta$dys[dta$tstart==times[j-1]]==1,1,
                                          rbinom(n,1,
                                                 expit(-5.5+0.05*(dta$age[dta$tstart==times[j]]-60)+
                                                         0.1*(dta$sex[dta$tstart==times[j]]==1)+
                                                         0.02*(dta$bmi[dta$tstart==times[j]]-30))))
    
    #history of cvd
    dta$cvd[dta$tstart==times[j]]<-ifelse(dta$cvd[dta$tstart==times[j-1]]==1,1,
                                          rbinom(n,1,expit(-8+0.5*(dta$age_group[dta$tstart==times[j]]==2)+2*(dta$age_group[dta$tstart==times[j]]==3)+3*(dta$age_group[dta$tstart==times[j]]==4)+
                                                             0.05*(dta$bmi[dta$tstart==times[j]]-30)+
                                                             0.5*dta$dys[dta$tstart==times[j]]+1*dta$hyp[dta$tstart==times[j]])))
    
    #history of kidney disease
    dta$kidney[dta$tstart==times[j]]<-ifelse(dta$kidney[dta$tstart==times[j-1]]==1,1,
                                             rbinom(n,1,
                                                    expit(-8+0.5*(dta$age_group[dta$tstart==times[j]]==2)+1*(dta$age_group[dta$tstart==times[j]]==3)+1.5*(dta$age_group[dta$tstart==times[j]]==4)+
                                                            0.2*dta$hb[dta$tstart==times[j]]+
                                                            0.5*dta$hyp[dta$tstart==times[j]]+
                                                            kid.coef.A*(dta$ts_A[dta$tstart==times[j]]>0)+
                                                            kid.coef.B*(dta$ts_B[dta$tstart==times[j]]>0)+
                                                            kid.coef.C*(dta$ts_C[dta$tstart==times[j]]>0)+
                                                            kid.coef.D*(dta$ts_D[dta$tstart==times[j]]>0))))
    
    #history of pancreatitis
    dta$panc[dta$tstart==times[j]]<-ifelse(dta$panc[dta$tstart==times[j-1]]==1,1,
                                           rbinom(n,1,expit(-12+0.1*dta$bmi[dta$tstart==times[j]]+
                                                              2*(dta$smok[dta$tstart==times[j]]==3)+
                                                              0.5*(dta$smok[dta$tstart==times[j]]==2))))
  }
  
  #--------------------------------------------
  #Treatments
  #--------------------------------------------
  
  #Times to starting treatments A-D, and resulting treatment status in a given time interval
  time_A<-times[j]+(-1/haz.A(dta$sex[dta$tstart==times[j]],dta$age[dta$tstart==times[j]],
                             dta$bmi[dta$tstart==times[j]],dta$diabdur[dta$tstart==times[j]],
                             dta$hb[dta$tstart==times[j]],
                             dta$hb_change[dta$tstart==times[j]],
                             (dta$smok[dta$tstart==times[j]]==2),
                             (dta$smok[dta$tstart==times[j]]==3),
                             dta$cvd[dta$tstart==times[j]],
                             dta$kidney[dta$tstart==times[j]]))*log(runif(n,0,1))
  
  time_B<-times[j]+(-1/haz.B(dta$sex[dta$tstart==times[j]],dta$age[dta$tstart==times[j]],
                             dta$bmi[dta$tstart==times[j]],dta$diabdur[dta$tstart==times[j]],
                             dta$hb[dta$tstart==times[j]],
                             dta$hb_change[dta$tstart==times[j]],
                             (dta$smok[dta$tstart==times[j]]==2),
                             (dta$smok[dta$tstart==times[j]]==3),
                             dta$cvd[dta$tstart==times[j]],
                             dta$kidney[dta$tstart==times[j]]))*log(runif(n,0,1))
  
  time_C<-times[j]+(-1/haz.C(dta$sex[dta$tstart==times[j]],dta$age[dta$tstart==times[j]],
                             dta$bmi[dta$tstart==times[j]],dta$diabdur[dta$tstart==times[j]],
                             dta$hb[dta$tstart==times[j]],
                             dta$hb_change[dta$tstart==times[j]],
                             (dta$smok[dta$tstart==times[j]]==2),
                             (dta$smok[dta$tstart==times[j]]==3),
                             dta$cvd[dta$tstart==times[j]],
                             dta$kidney[dta$tstart==times[j]]))*log(runif(n,0,1))
  
  time_D<-times[j]+(-1/haz.D(dta$sex[dta$tstart==times[j]],dta$age[dta$tstart==times[j]],
                             dta$bmi[dta$tstart==times[j]],dta$diabdur[dta$tstart==times[j]],
                             dta$hb[dta$tstart==times[j]],
                             dta$hb_change[dta$tstart==times[j]],
                             (dta$smok[dta$tstart==times[j]]==2),
                             (dta$smok[dta$tstart==times[j]]==3),
                             dta$cvd[dta$tstart==times[j]],
                             dta$kidney[dta$tstart==times[j]]))*log(runif(n,0,1))
  
  time_treat<-cbind(time_A,time_B,time_C,time_D)
  time_treat_whichmin<-sapply(1:n,FUN=function(x){which.min(time_treat[x,])})
  
  if(j==1){
    dta$treat_time_A[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==1 & time_A<=dta$tstop[j],time_A,NA)
    dta$treat_time_B[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==2 & time_B<=dta$tstop[j],time_B,NA)
    dta$treat_time_C[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==3 & time_C<=dta$tstop[j],time_C,NA)
    dta$treat_time_D[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==4 & time_D<=dta$tstop[j],time_D,NA)
    
    dta$treat_status_A[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==1 & time_A<=dta$tstop[j],1,0)
    dta$treat_status_B[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==2 & time_B<=dta$tstop[j],1,0)
    dta$treat_status_C[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==3 & time_C<=dta$tstop[j],1,0)
    dta$treat_status_D[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==4 & time_D<=dta$tstop[j],1,0)
  }else if(j>1){
    not_treated_yet<-I(is.na(dta$treat_time_A[dta$tstart==times[j-1]]) &
                         is.na(dta$treat_time_B[dta$tstart==times[j-1]]) &
                         is.na(dta$treat_time_C[dta$tstart==times[j-1]]) &
                         is.na(dta$treat_time_D[dta$tstart==times[j-1]]))
    
    dta$treat_time_A[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==1 & time_A<=dta$tstop[j] & 
                                                     not_treated_yet,time_A,
                                                   dta$treat_time_A[dta$tstart==times[j-1]])
    
    dta$treat_time_B[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==2 & time_B<=dta$tstop[j] & 
                                                     not_treated_yet,time_B,
                                                   dta$treat_time_B[dta$tstart==times[j-1]])
    
    dta$treat_time_C[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==3 & time_C<=dta$tstop[j] & 
                                                     not_treated_yet,time_C,
                                                   dta$treat_time_C[dta$tstart==times[j-1]])
    
    dta$treat_time_D[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==4 & time_D<=dta$tstop[j] & 
                                                     not_treated_yet,time_D,
                                                   dta$treat_time_D[dta$tstart==times[j-1]])
    
    dta$treat_status_A[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==1 & time_A<=dta$tstop[j] & 
                                                       not_treated_yet,1,
                                                     dta$treat_status_A[dta$tstart==times[j-1]])
    
    dta$treat_status_B[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==2 & time_B<=dta$tstop[j] & 
                                                       not_treated_yet,1,
                                                     dta$treat_status_B[dta$tstart==times[j-1]])
    
    dta$treat_status_C[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==3 & time_C<=dta$tstop[j] & 
                                                       not_treated_yet,1,
                                                     dta$treat_status_C[dta$tstart==times[j-1]])
    
    dta$treat_status_D[dta$tstart==times[j]]<-ifelse(time_treat_whichmin==4 & time_D<=dta$tstop[j] & 
                                                       not_treated_yet,1,
                                                     dta$treat_status_D[dta$tstart==times[j-1]])
  }
  
  #time since treatment start
  if(j==1){
    dta$ts_A[dta$tstart==times[j]]<-0
    dta$ts_B[dta$tstart==times[j]]<-0
    dta$ts_C[dta$tstart==times[j]]<-0
    dta$ts_D[dta$tstart==times[j]]<-0
  }
  if(j<ntimes){
    dta$ts_A[dta$tstart==times[j+1]]<-ifelse(!is.na(dta$treat_time_A[dta$tstart==times[j]]),times[j+1]-dta$treat_time_A[dta$tstart==times[j]],0)/365.25
    dta$ts_B[dta$tstart==times[j+1]]<-ifelse(!is.na(dta$treat_time_B[dta$tstart==times[j]]),times[j+1]-dta$treat_time_B[dta$tstart==times[j]],0)/365.25
    dta$ts_C[dta$tstart==times[j+1]]<-ifelse(!is.na(dta$treat_time_C[dta$tstart==times[j]]),times[j+1]-dta$treat_time_C[dta$tstart==times[j]],0)/365.25
    dta$ts_D[dta$tstart==times[j+1]]<-ifelse(!is.na(dta$treat_time_D[dta$tstart==times[j]]),times[j+1]-dta$treat_time_D[dta$tstart==times[j]],0)/365.25
  }
  
  #--------------------------------------------
  #Events: mace, death, censoring
  #--------------------------------------------
  
  time_mace<-times[j]+(-1/haz.mace(dta$age[dta$tstart==times[j]],dta$sex[dta$tstart==times[j]],
                                   dta$smok_former[dta$tstart==times[j]],dta$smok_current[dta$tstart==times[j]],
                                   dta$bmi[dta$tstart==times[j]],dta$hb[dta$tstart==times[j]],
                                   dta$diabdur[dta$tstart==times[j]],
                                   dta$dys[dta$tstart==times[j]],dta$hyp[dta$tstart==times[j]],
                                   dta$cvd[dta$tstart==times[j]],dta$kidney[dta$tstart==times[j]],
                                   dta$treat_status_A[dta$tstart==times[j]],
                                   dta$treat_status_B[dta$tstart==times[j]],
                                   dta$treat_status_C[dta$tstart==times[j]],
                                   dta$treat_status_D[dta$tstart==times[j]]))*log(runif(n,0,1))
  
  time_death<-times[j]+(-1/haz.death(dta$age[dta$tstart==times[j]],dta$sex[dta$tstart==times[j]],
                                     dta$smok_former[dta$tstart==times[j]],dta$smok_current[dta$tstart==times[j]],
                                     dta$bmi[dta$tstart==times[j]],dta$hb[dta$tstart==times[j]],
                                     dta$diabdur[dta$tstart==times[j]],
                                     dta$dys[dta$tstart==times[j]],dta$hyp[dta$tstart==times[j]],
                                     dta$cvd[dta$tstart==times[j]],dta$kidney[dta$tstart==times[j]],
                                     dta$treat_status_A[dta$tstart==times[j]],
                                     dta$treat_status_B[dta$tstart==times[j]],
                                     dta$treat_status_C[dta$tstart==times[j]],
                                     dta$treat_status_D[dta$tstart==times[j]]))*log(runif(n,0,1))
  
  time_cens<-times[j]+(-1/haz.cens)*log(runif(n,0,1))
  
  time_event<-cbind(time_mace,time_death,time_cens)
  time_event_whichmin<-sapply(1:n,FUN=function(x){which.min(time_event[x,])})
  
  if(j==1){
    dta$mace_time[dta$tstart==times[j]]<-ifelse(time_event_whichmin==1 & time_mace<=dta$tstop[j],time_mace,NA)
    dta$death_time[dta$tstart==times[j]]<-ifelse(time_event_whichmin==2 & time_death<=dta$tstop[j],time_death,NA)
    dta$cens_time[dta$tstart==times[j]]<-ifelse(time_event_whichmin==3 & time_cens<=dta$tstop[j],time_cens,NA)
    
    dta$mace_status[dta$tstart==times[j]]<-ifelse(time_event_whichmin==1 & time_mace<=dta$tstop[j],1,0)
    dta$death_status[dta$tstart==times[j]]<-ifelse(time_event_whichmin==2 & time_death<=dta$tstop[j],1,0)
    dta$cens_status[dta$tstart==times[j]]<-ifelse(time_event_whichmin==3 & time_cens<=dta$tstop[j],1,0)
  }else if(j>1){
    
    no_event_yet<-I(is.na(dta$mace_time[dta$tstart==times[j-1]]) & 
                      is.na(dta$death_time[dta$tstart==times[j-1]]) &
                      is.na(dta$cens_time[dta$tstart==times[j-1]]))
    
    dta$mace_time[dta$tstart==times[j]]<-ifelse(time_event_whichmin==1 & time_mace<=dta$tstop[j] & 
                                                  no_event_yet,time_mace,
                                                dta$mace_time[dta$tstart==times[j-1]])
    dta$death_time[dta$tstart==times[j]]<-ifelse(time_event_whichmin==2 & time_death<=dta$tstop[j] & 
                                                   no_event_yet,time_death,
                                                 dta$death_time[dta$tstart==times[j-1]])
    dta$cens_time[dta$tstart==times[j]]<-ifelse(time_event_whichmin==3 & time_cens<=dta$tstop[j] & 
                                                  no_event_yet,time_cens,
                                                dta$cens_time[dta$tstart==times[j-1]])
    
    dta$mace_status[dta$tstart==times[j]]<-ifelse(time_event_whichmin==1 & time_mace<=dta$tstop[j] & 
                                                    no_event_yet,1,
                                                  dta$mace_status[dta$tstart==times[j-1]])
    dta$death_status[dta$tstart==times[j]]<-ifelse(time_event_whichmin==2 & time_death<=dta$tstop[j] & 
                                                     no_event_yet,1,
                                                   dta$death_status[dta$tstart==times[j-1]])
    dta$cens_status[dta$tstart==times[j]]<-ifelse(time_event_whichmin==3 & time_cens<=dta$tstop[j] & 
                                                    no_event_yet,1,
                                                  dta$cens_status[dta$tstart==times[j-1]])
  }
  
}

#--------------------------------------------
#end of data generation loop
#--------------------------------------------

#--------------------------------------------
#some additional data management
#--------------------------------------------

#---
#exclude rows after event or censoring

#first set mace_time, death_time, cens_time to be the same in every row for an individual
dta<-dta%>%group_by(id)%>%fill(mace_time,.direction="up")
dta<-dta%>%group_by(id)%>%fill(death_time,.direction="up")
dta<-dta%>%group_by(id)%>%fill(cens_time,.direction="up")

#same for treatment times
dta<-dta%>%group_by(id)%>%fill(treat_time_A,.direction="up")
dta<-dta%>%group_by(id)%>%fill(treat_time_D,.direction="up")
dta<-dta%>%group_by(id)%>%fill(treat_time_B,.direction="up")
dta<-dta%>%group_by(id)%>%fill(treat_time_C,.direction="up")

#generate combined event time and status
dta$event_time<-pmin(dta$mace_time,dta$death_time,dta$cens_time,na.rm=T)

dta$event_type<-0
dta$event_type<-ifelse(dta$event_time==dta$mace_time & !is.na(dta$mace_time),1,dta$event_type)
dta$event_type<-ifelse(dta$event_time==dta$death_time & !is.na(dta$death_time),2,dta$event_type)

#set event_time to admin.cens.time for people who do not have an event time for any event type
dta$event_type<-ifelse(is.na(dta$event_time),3,dta$event_type)
dta$event_time<-ifelse(is.na(dta$event_time),admin.cens.time,dta$event_time)

#exclude rows after event or censoring time
dta<-dta[dta$tstart<dta$event_time,]

#set tstop to event_time in the last row (this is for people who have the event in the last time period)
dta$tstop<-ifelse(dta$tstop>dta$event_time,dta$event_time,dta$tstop)

#---
#amending treatment times that are after the event/cens time
dta$treat_time_A<-ifelse(dta$treat_time_A>dta$event_time,NA,dta$treat_time_A)
dta$treat_time_B<-ifelse(dta$treat_time_B>dta$event_time,NA,dta$treat_time_B)
dta$treat_time_C<-ifelse(dta$treat_time_C>dta$event_time,NA,dta$treat_time_C)
dta$treat_time_D<-ifelse(dta$treat_time_D>dta$event_time,NA,dta$treat_time_D)

#---
#Dropping variables not to be included in final data set 
#(they can be derived from the available variables)
dta$age_group<-NULL
dta$hb_change<-NULL
dta$ts_A<-NULL
dta$ts_B<-NULL
dta$ts_C<-NULL
dta$ts_D<-NULL

