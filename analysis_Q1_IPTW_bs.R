#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#Using IPTW
#
#This file bootstraps the analysis and obtains 95% CIs
#######################################################

#------------------------------------------------------
#Store the original data 
#(this means we don't need to alter the analysis code too much because the name of the data is the same as in the main analysis)
#------------------------------------------------------

dta.analysis.orig<-dta.analysis

#------------------------------------------------------
#Storage of bootstrap estimates
#------------------------------------------------------

cinc.A.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.B.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.C.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.D.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))

cinc.A.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.B.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.C.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.D.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))

diff.B.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.C.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.D.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))

diff.B.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.C.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.D.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))

#------------------------------------------------------
#Bootstrap the analysis dataset: dta.analysis.orig
#Perform analysis in each bootstrap sample
#------------------------------------------------------

for(b in 1:nboot){
  
  print(b)

  samp<-sample(1:dim(dta.analysis.orig)[1],replace=T)
  
  dta.analysis<-dta.analysis.orig[samp,]
  
  #------------------------------------------------------
  #Analysis: Estimate the IPTW
  #------------------------------------------------------
  
  #estimate weights
  mod.trtA.num<-glm(treat_status_A~1,data=dta.analysis,family="binomial")
  mod.trtA.denom<-glm(treat_status_A~sex+age+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                        hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")
  
  mod.trtB.num<-glm(treat_status_B~1,data=dta.analysis,family="binomial")
  mod.trtB.denom<-glm(treat_status_B~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                        hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")
  
  mod.trtC.num<-glm(treat_status_C~1,data=dta.analysis,family="binomial")
  mod.trtC.denom<-glm(treat_status_C~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                        hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")
  
  mod.trtD.num<-glm(treat_status_D~1,data=dta.analysis,family="binomial")
  mod.trtD.denom<-glm(treat_status_D~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                        hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")
  
  dta.analysis$prob.trtA.num<-predict(mod.trtA.num,newdata=dta.analysis,type="response")
  dta.analysis$prob.trtB.num<-predict(mod.trtB.num,newdata=dta.analysis,type="response")
  dta.analysis$prob.trtC.num<-predict(mod.trtC.num,newdata=dta.analysis,type="response")
  dta.analysis$prob.trtD.num<-predict(mod.trtD.num,newdata=dta.analysis,type="response")
  
  dta.analysis$prob.trtA.denom<-predict(mod.trtA.denom,newdata=dta.analysis,type="response")
  dta.analysis$prob.trtB.denom<-predict(mod.trtB.denom,newdata=dta.analysis,type="response")
  dta.analysis$prob.trtC.denom<-predict(mod.trtC.denom,newdata=dta.analysis,type="response")
  dta.analysis$prob.trtD.denom<-predict(mod.trtD.denom,newdata=dta.analysis,type="response")
  
  dta.analysis$iptw<-dta.analysis$treat_status_A*(dta.analysis$prob.trtA.num/dta.analysis$prob.trtA.denom)+
    dta.analysis$treat_status_B*(dta.analysis$prob.trtB.num/dta.analysis$prob.trtB.denom)+
    dta.analysis$treat_status_C*(dta.analysis$prob.trtC.num/dta.analysis$prob.trtC.denom)+
    dta.analysis$treat_status_D*(dta.analysis$prob.trtD.num/dta.analysis$prob.trtD.denom)
  
  #------------------------------------------------------
  #Analysis: estimate cumulative incidences using IPTW
  #------------------------------------------------------
  
  cr.trtA.iptw<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis[dta.analysis$treat_type=="A",],weights = dta.analysis$iptw[dta.analysis$treat_type=="A"])
  cr.trtB.iptw<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis[dta.analysis$treat_type=="B",],weights = dta.analysis$iptw[dta.analysis$treat_type=="B"])
  cr.trtC.iptw<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis[dta.analysis$treat_type=="C",],weights = dta.analysis$iptw[dta.analysis$treat_type=="C"])
  cr.trtD.iptw<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis[dta.analysis$treat_type=="D",],weights = dta.analysis$iptw[dta.analysis$treat_type=="D"])
  
  #obtain cumulative incidences for each time in time.grid
  
  cinc.A.times.iptw<-cr.trtA.iptw$time/365.25
  cinc.B.times.iptw<-cr.trtB.iptw$time/365.25
  cinc.C.times.iptw<-cr.trtC.iptw$time/365.25
  cinc.D.times.iptw<-cr.trtD.iptw$time/365.25
  
  cincstep.A.mace.iptw<-stepfun(cr.trtA.iptw$time/365.25,c(0,cr.trtA.iptw$pstate[,2]))
  cincstep.B.mace.iptw<-stepfun(cr.trtB.iptw$time/365.25,c(0,cr.trtB.iptw$pstate[,2]))
  cincstep.C.mace.iptw<-stepfun(cr.trtC.iptw$time/365.25,c(0,cr.trtC.iptw$pstate[,2]))
  cincstep.D.mace.iptw<-stepfun(cr.trtD.iptw$time/365.25,c(0,cr.trtD.iptw$pstate[,2]))
  
  cincstep.A.death.iptw<-stepfun(cr.trtA.iptw$time/365.25,c(0,cr.trtA.iptw$pstate[,3]))
  cincstep.B.death.iptw<-stepfun(cr.trtB.iptw$time/365.25,c(0,cr.trtB.iptw$pstate[,3]))
  cincstep.C.death.iptw<-stepfun(cr.trtC.iptw$time/365.25,c(0,cr.trtC.iptw$pstate[,3]))
  cincstep.D.death.iptw<-stepfun(cr.trtD.iptw$time/365.25,c(0,cr.trtD.iptw$pstate[,3]))
  
  cinc.A.mace.iptw.bs[b,]<-cincstep.A.mace.iptw(time.grid)
  cinc.B.mace.iptw.bs[b,]<-cincstep.B.mace.iptw(time.grid)
  cinc.C.mace.iptw.bs[b,]<-cincstep.C.mace.iptw(time.grid)
  cinc.D.mace.iptw.bs[b,]<-cincstep.D.mace.iptw(time.grid)
  
  cinc.A.death.iptw.bs[b,]<-cincstep.A.death.iptw(time.grid)
  cinc.B.death.iptw.bs[b,]<-cincstep.B.death.iptw(time.grid)
  cinc.C.death.iptw.bs[b,]<-cincstep.C.death.iptw(time.grid)
  cinc.D.death.iptw.bs[b,]<-cincstep.D.death.iptw(time.grid)
  
  #obtain effect estimates
  
  diff.B.mace.iptw.bs[b,]<-cinc.B.mace.iptw-cinc.A.mace.iptw
  diff.C.mace.iptw.bs[b,]<-cinc.C.mace.iptw-cinc.A.mace.iptw
  diff.D.mace.iptw.bs[b,]<-cinc.D.mace.iptw-cinc.A.mace.iptw
  
  diff.B.death.iptw.bs[b,]<-cinc.B.death.iptw-cinc.A.death.iptw
  diff.C.death.iptw.bs[b,]<-cinc.C.death.iptw-cinc.A.death.iptw
  diff.D.death.iptw.bs[b,]<-cinc.D.death.iptw-cinc.A.death.iptw
}

#------------------------------------------------------
#bootstrap 95% CIs for cumulative incidences and differences (percentile method)
#------------------------------------------------------

cinc.A.mace.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.A.mace.iptw.bs[,x],0.05)})
cinc.B.mace.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.B.mace.iptw.bs[,x],0.05)})
cinc.C.mace.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.C.mace.iptw.bs[,x],0.05)})
cinc.D.mace.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.D.mace.iptw.bs[,x],0.05)})

cinc.A.mace.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.A.mace.iptw.bs[,x],0.95)})
cinc.B.mace.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.B.mace.iptw.bs[,x],0.95)})
cinc.C.mace.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.C.mace.iptw.bs[,x],0.95)})
cinc.D.mace.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.D.mace.iptw.bs[,x],0.95)})

cinc.A.death.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.A.death.iptw.bs[,x],0.05)})
cinc.B.death.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.B.death.iptw.bs[,x],0.05)})
cinc.C.death.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.C.death.iptw.bs[,x],0.05)})
cinc.D.death.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.D.death.iptw.bs[,x],0.05)})

cinc.A.death.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.A.death.iptw.bs[,x],0.95)})
cinc.B.death.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.B.death.iptw.bs[,x],0.95)})
cinc.C.death.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.C.death.iptw.bs[,x],0.95)})
cinc.D.death.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.D.death.iptw.bs[,x],0.95)})

diff.B.mace.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.B.mace.iptw.bs[,x],0.05)})
diff.C.mace.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.C.mace.iptw.bs[,x],0.05)})
diff.D.mace.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.D.mace.iptw.bs[,x],0.05)})

diff.B.mace.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.B.mace.iptw.bs[,x],0.95)})
diff.C.mace.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.C.mace.iptw.bs[,x],0.95)})
diff.D.mace.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.D.mace.iptw.bs[,x],0.95)})

diff.B.death.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.B.death.iptw.bs[,x],0.05)})
diff.C.death.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.C.death.iptw.bs[,x],0.05)})
diff.D.death.iptw.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.D.death.iptw.bs[,x],0.05)})

diff.B.death.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.B.death.iptw.bs[,x],0.95)})
diff.C.death.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.C.death.iptw.bs[,x],0.95)})
diff.D.death.iptw.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.D.death.iptw.bs[,x],0.95)})


