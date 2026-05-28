#######################################################
#Q2: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D 
#or never using treatment, in people with hb>=7.5. 
#
#Using IPTW to handle baseline confounding, 
# and using IPCW to handle censoring due to deviation from non-treatment (for the 'never treat' strategy).
#
#This file bootstraps the analysis and obtains 95% CIs
#######################################################

#------------------------------------------------------
#Store the original data 
#(this means we don't need to alter the analysis code in the bootstrap loop because the name of the data is the same as in the main analysis)
#------------------------------------------------------

dta.analysis.orig<-dta.analysis

#------------------------------------------------------
#Storage of bootstrap estimates
#------------------------------------------------------

cinc.A.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.B.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.C.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.D.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.nevertrt.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))

cinc.A.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.B.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.C.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.D.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.nevertrt.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))

diff.A.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.B.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.C.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.D.mace.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))

diff.A.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.B.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.C.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.D.death.iptw.bs<-matrix(nrow=nboot,ncol=length(time.grid))

#------------------------------------------------------
#Bootstrap the analysis dataset: dta.analysis.orig
#Perform analysis in each bootstrap sample
#------------------------------------------------------

for(b in 1:nboot){
  
  print(b)
  
  #bootstrap sample
  dta.analysis<-bootstrap_cluster(dta.analysis.orig,"id_new")
  
  #------------------------------------------------------
  #Fit weights models to handle artificial censoring in the never treated part of the data
  #The censoring occurs due to starting different treatments, and these depend on covariates in different ways
  #Hence we fit separate models for censoring due to starting different treatments.
  #The models are for treatment initiation, among people who have not yet initiated treatment A,B,C, or D.
  #------------------------------------------------------
  
  #Model for probability of starting treatment A in a given time period
  mod.trtA.num<-glm(treat_status_A~1,
                    data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
  mod.trtA.denom<-glm(treat_status_A~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                        hyp+dys+cvd+kidney+panc,
                      data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
  
  #Model for probability of starting treatment B in a given time period
  mod.trtB.num<-glm(treat_status_B~1,
                    data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
  mod.trtB.denom<-glm(treat_status_B~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                        hyp+dys+cvd+kidney+panc,
                      data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
  #Model for probability of starting treatment C in a given time period
  mod.trtC.num<-glm(treat_status_C~1,
                    data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
  mod.trtC.denom<-glm(treat_status_C~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                        hyp+dys+cvd+kidney+panc,
                      data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
  
  #Model for probability of starting treatment D in a given time period
  mod.trtD.num<-glm(treat_status_D~1,
                    data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
  mod.trtD.denom<-glm(treat_status_D~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                        hyp+dys+cvd+kidney+panc,
                      data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
  
  #probability of remaining uncensored in a given time period
  num.A<-1-predict(mod.trtA.num,newdata=dta.analysis[dta.analysis$treat_type=="nevertrt",],type="response")
  num.B<-1-predict(mod.trtB.num,newdata=dta.analysis[dta.analysis$treat_type=="nevertrt",],type="response")
  num.C<-1-predict(mod.trtC.num,newdata=dta.analysis[dta.analysis$treat_type=="nevertrt",],type="response")
  num.D<-1-predict(mod.trtD.num,newdata=dta.analysis[dta.analysis$treat_type=="nevertrt",],type="response")
  denom.A<-1-predict(mod.trtA.denom,newdata=dta.analysis[dta.analysis$treat_type=="nevertrt",],type="response")
  denom.B<-1-predict(mod.trtB.denom,newdata=dta.analysis[dta.analysis$treat_type=="nevertrt",],type="response")
  denom.C<-1-predict(mod.trtC.denom,newdata=dta.analysis[dta.analysis$treat_type=="nevertrt",],type="response")
  denom.D<-1-predict(mod.trtD.denom,newdata=dta.analysis[dta.analysis$treat_type=="nevertrt",],type="response")
  num.anytrt<-num.A*num.B*num.C*num.D
  denom.anytrt<-denom.A*denom.B*denom.C*denom.D
  
  #ipcw
  dta.analysis$ipcw<-1
  dta.analysis[dta.analysis$treat_type=="nevertrt",]$ipcw<-num.anytrt/denom.anytrt
  dta.analysis<-dta.analysis%>%group_by(id_new)%>%mutate(ipcw=cumprod(ipcw))#important to use id_new here
  
  #drop last row of data in dta.nevertrt (this was just needed for estimation of the weights)
  dta.analysis2<-dta.analysis[dta.analysis$treat_type!="nevertrt"|
                                dta.analysis$treat_type=="nevertrt" & dta.analysis$treat_status_any==0,]
  
  #first row of dta.analysis, which is the population for which we will obtain average treatment effects
  dta.analysis2.row1<-dta.analysis2[dta.analysis2$tstart==0,]
  
  #------------------------------------------------------
  #IPTW analysis: estimation of baseline weights
  #------------------------------------------------------
  
  #fit models for each treatment: for model denominators
  mod.trtA<-glm(I(treat_type=="A")~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                  hyp+dys+cvd+kidney+panc,data=dta.analysis2.row1,family="binomial")
  
  mod.trtB<-glm(I(treat_type=="B")~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                  hyp+dys+cvd+kidney+panc,data=dta.analysis2.row1,family="binomial")
  
  mod.trtC<-glm(I(treat_type=="C")~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                  hyp+dys+cvd+kidney+panc,data=dta.analysis2.row1,family="binomial")
  
  mod.trtD<-glm(I(treat_type=="D")~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                  hyp+dys+cvd+kidney+panc,data=dta.analysis2.row1,family="binomial")
  
  mod.nevertrt<-glm(I(treat_type=="nevertrt")~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,data=dta.analysis2.row1,family="binomial")
  
  #fit models for each treatment: for model numerators
  smod.trtA<-glm(I(treat_type=="A")~1,data=dta.analysis2.row1,family="binomial")
  smod.trtB<-glm(I(treat_type=="B")~1,data=dta.analysis2.row1,family="binomial")
  smod.trtC<-glm(I(treat_type=="C")~1,data=dta.analysis2.row1,family="binomial")
  smod.trtD<-glm(I(treat_type=="D")~1,data=dta.analysis2.row1,family="binomial")
  smod.nevertrt<-glm(I(treat_type=="nevertrt")~1,data=dta.analysis2.row1,family="binomial")
  
  #obtain baseline iptw
  dta.analysis2.row1$iptw.trtA<-predict(smod.trtA,newdata=dta.analysis2.row1,type="response")/predict(mod.trtA,newdata=dta.analysis2.row1,type="response")
  dta.analysis2.row1$iptw.trtB<-predict(smod.trtB,newdata=dta.analysis2.row1,type="response")/predict(mod.trtB,newdata=dta.analysis2.row1,type="response")
  dta.analysis2.row1$iptw.trtC<-predict(smod.trtC,newdata=dta.analysis2.row1,type="response")/predict(mod.trtC,newdata=dta.analysis2.row1,type="response")
  dta.analysis2.row1$iptw.trtD<-predict(smod.trtD,newdata=dta.analysis2.row1,type="response")/predict(mod.trtD,newdata=dta.analysis2.row1,type="response")
  dta.analysis2.row1$iptw.nevertrt<-predict(smod.nevertrt,newdata=dta.analysis2.row1,type="response")/predict(mod.nevertrt,newdata=dta.analysis2.row1,type="response")
  
  dta.analysis2.row1$iptw<-NA
  dta.analysis2.row1$iptw<-ifelse(dta.analysis2.row1$treat_type=="A",dta.analysis2.row1$iptw.trtA,dta.analysis2.row1$iptw)
  dta.analysis2.row1$iptw<-ifelse(dta.analysis2.row1$treat_type=="B",dta.analysis2.row1$iptw.trtB,dta.analysis2.row1$iptw)
  dta.analysis2.row1$iptw<-ifelse(dta.analysis2.row1$treat_type=="C",dta.analysis2.row1$iptw.trtC,dta.analysis2.row1$iptw)
  dta.analysis2.row1$iptw<-ifelse(dta.analysis2.row1$treat_type=="D",dta.analysis2.row1$iptw.trtD,dta.analysis2.row1$iptw)
  dta.analysis2.row1$iptw<-ifelse(dta.analysis2.row1$treat_type=="nevertrt",dta.analysis2.row1$iptw.nevertrt,dta.analysis2.row1$iptw)
  
  #add the iptw into the long data set (long for nevertrt), dta.analysis
  dta.analysis2$iptw<-NA
  dta.analysis2$iptw[dta.analysis2$tstart==0]<-dta.analysis2.row1$iptw
  
  dta.analysis2<-dta.analysis2%>%group_by(id_new)%>%fill(iptw,.direction = "down")
  
  #generate combined iptw*ipcw weight
  dta.analysis2$ipw.comb<-dta.analysis2$iptw*dta.analysis2$ipcw
  
  #------------------------------------------------------
  #IPTW analysis: weighted Aalen-Johansen estimates
  #------------------------------------------------------
  
  cr.trtA.iptw<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                        data=dta.analysis2[dta.analysis2$treat_type=="A",],
                        weights = dta.analysis2$ipw.comb[dta.analysis2$treat_type=="A"],
                        id=dta.analysis2$id_new[dta.analysis2$treat_type=="A"])
  
  cr.trtB.iptw<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                        data=dta.analysis2[dta.analysis2$treat_type=="B",],
                        weights = dta.analysis2$ipw.comb[dta.analysis2$treat_type=="B"],
                        id=dta.analysis2$id_new[dta.analysis2$treat_type=="B"])
  
  cr.trtC.iptw<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                        data=dta.analysis2[dta.analysis2$treat_type=="C",],
                        weights = dta.analysis2$ipw.comb[dta.analysis2$treat_type=="C"],
                        id=dta.analysis2$id_new[dta.analysis2$treat_type=="C"])
  
  cr.trtD.iptw<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                        data=dta.analysis2[dta.analysis2$treat_type=="D",],
                        weights = dta.analysis2$ipw.comb[dta.analysis2$treat_type=="D"],
                        id=dta.analysis2$id_new[dta.analysis2$treat_type=="D"])
  
  cr.nevertrt.iptw<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                            data=dta.analysis2[dta.analysis2$treat_type=="nevertrt",],
                            weights = dta.analysis2$ipw.comb[dta.analysis2$treat_type=="nevertrt"],
                            id=dta.analysis2$id_new[dta.analysis2$treat_type=="nevertrt"])
  
  #obtain cumulative incidences for each time in time.grid
  
  cinc.A.times.iptw<-cr.trtA.iptw$time/365.25
  cinc.B.times.iptw<-cr.trtB.iptw$time/365.25
  cinc.C.times.iptw<-cr.trtC.iptw$time/365.25
  cinc.D.times.iptw<-cr.trtD.iptw$time/365.25
  cinc.nevertrt.times.iptw<-cr.nevertrt.iptw$time/365.25
  
  cincstep.A.mace.iptw<-stepfun(cr.trtA.iptw$time/365.25,c(0,cr.trtA.iptw$pstate[,2]))
  cincstep.B.mace.iptw<-stepfun(cr.trtB.iptw$time/365.25,c(0,cr.trtB.iptw$pstate[,2]))
  cincstep.C.mace.iptw<-stepfun(cr.trtC.iptw$time/365.25,c(0,cr.trtC.iptw$pstate[,2]))
  cincstep.D.mace.iptw<-stepfun(cr.trtD.iptw$time/365.25,c(0,cr.trtD.iptw$pstate[,2]))
  cincstep.nevertrt.mace.iptw<-stepfun(cr.nevertrt.iptw$time/365.25,c(0,cr.nevertrt.iptw$pstate[,2]))
  
  cincstep.A.death.iptw<-stepfun(cr.trtA.iptw$time/365.25,c(0,cr.trtA.iptw$pstate[,3]))
  cincstep.B.death.iptw<-stepfun(cr.trtB.iptw$time/365.25,c(0,cr.trtB.iptw$pstate[,3]))
  cincstep.C.death.iptw<-stepfun(cr.trtC.iptw$time/365.25,c(0,cr.trtC.iptw$pstate[,3]))
  cincstep.D.death.iptw<-stepfun(cr.trtD.iptw$time/365.25,c(0,cr.trtD.iptw$pstate[,3]))
  cincstep.nevertrt.death.iptw<-stepfun(cr.nevertrt.iptw$time/365.25,c(0,cr.nevertrt.iptw$pstate[,3]))
  
  cinc.A.mace.iptw.bs[b,]<-cincstep.A.mace.iptw(time.grid)
  cinc.B.mace.iptw.bs[b,]<-cincstep.B.mace.iptw(time.grid)
  cinc.C.mace.iptw.bs[b,]<-cincstep.C.mace.iptw(time.grid)
  cinc.D.mace.iptw.bs[b,]<-cincstep.D.mace.iptw(time.grid)
  cinc.nevertrt.mace.iptw.bs[b,]<-cincstep.nevertrt.mace.iptw(time.grid)
  
  cinc.A.death.iptw.bs[b,]<-cincstep.A.death.iptw(time.grid)
  cinc.B.death.iptw.bs[b,]<-cincstep.B.death.iptw(time.grid)
  cinc.C.death.iptw.bs[b,]<-cincstep.C.death.iptw(time.grid)
  cinc.D.death.iptw.bs[b,]<-cincstep.D.death.iptw(time.grid)
  cinc.nevertrt.death.iptw.bs[b,]<-cincstep.nevertrt.death.iptw(time.grid)
  
  #obtain effect estimates
  
  diff.A.mace.iptw.bs[b,]<-cinc.A.mace.iptw-cinc.nevertrt.mace.iptw
  diff.B.mace.iptw.bs[b,]<-cinc.B.mace.iptw-cinc.nevertrt.mace.iptw
  diff.C.mace.iptw.bs[b,]<-cinc.C.mace.iptw-cinc.nevertrt.mace.iptw
  diff.D.mace.iptw.bs[b,]<-cinc.D.mace.iptw-cinc.nevertrt.mace.iptw
  
  diff.A.death.iptw.bs[b,]<-cinc.A.death.iptw-cinc.nevertrt.death.iptw
  diff.B.death.iptw.bs[b,]<-cinc.B.death.iptw-cinc.nevertrt.death.iptw
  diff.C.death.iptw.bs[b,]<-cinc.C.death.iptw-cinc.nevertrt.death.iptw
  diff.D.death.iptw.bs[b,]<-cinc.D.death.iptw-cinc.nevertrt.death.iptw
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


