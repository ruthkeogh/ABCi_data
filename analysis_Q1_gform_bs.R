#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#Using g-formula (i.e. regression adjustment followed by standardisation)
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

cinc.A.mace.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.B.mace.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.C.mace.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.D.mace.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))

cinc.A.death.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.B.death.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.C.death.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
cinc.D.death.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))

diff.B.mace.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.C.mace.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.D.mace.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))

diff.B.death.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.C.death.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))
diff.D.death.gform.bs<-matrix(nrow=nboot,ncol=length(time.grid))

#------------------------------------------------------
#Bootstrap the analysis dataset: dta.analysis.orig
#Perform analysis in each bootstrap sample
#------------------------------------------------------

for(b in 1:nboot){
  
  print(b)
  
  samp<-sample(1:dim(dta.analysis.orig)[1],replace=T)
  
  dta.analysis<-dta.analysis.orig[samp,]
  
  #------------------------------------------------------
  #Analysis: standardisation, competing events
  #------------------------------------------------------
  
  #fit two cause-specific hazard models
  cox.mod.cr.mace<-coxph(Surv(event_time_new,event_type_cr==1)~as.factor(treat_type)+
                           sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                           hyp+dys+cvd+kidney+panc,data=dta.analysis,id=dta.analysis$id)
  cox.mod.cr.death<-coxph(Surv(event_time_new,event_type_cr==2)~as.factor(treat_type)+
                            sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                            hyp+dys+cvd+kidney+panc,data=dta.analysis,id=dta.analysis$id)
  
  #get baseline cumulative hazards, and baseline hazards (at each time)
  cbhaz.mace<-basehaz(cox.mod.cr.mace,centered = F)#baseline cumulative haz for mace, at each time
  cbhaz.death<-basehaz(cox.mod.cr.death,centered = F)#baseline cumulative haz for death, at each time
  
  bhaz.mace<-c(cbhaz.mace$haz[1],diff(cbhaz.mace$haz))#baseline hazard for mace
  bhaz.death<-c(cbhaz.death$haz[1],diff(cbhaz.death$haz))#baseline hazard for death
  
  #get linear predictors under each treatment
  newdata.A<-dta.analysis%>%mutate(treat_type="A")
  lp.mace.A<-predict(cox.mod.cr.mace,newdata=newdata.A,type="lp",reference="zero")
  lp.death.A<-predict(cox.mod.cr.death,newdata=newdata.A,type="lp",reference="zero")
  
  newdata.B<-dta.analysis%>%mutate(treat_type="B")
  lp.mace.B<-predict(cox.mod.cr.mace,newdata=newdata.B,type="lp",reference="zero")
  lp.death.B<-predict(cox.mod.cr.death,newdata=newdata.B,type="lp",reference="zero")
  
  newdata.C<-dta.analysis%>%mutate(treat_type="C")
  lp.mace.C<-predict(cox.mod.cr.mace,newdata=newdata.C,type="lp",reference="zero")
  lp.death.C<-predict(cox.mod.cr.death,newdata=newdata.C,type="lp",reference="zero")
  
  newdata.D<-dta.analysis%>%mutate(treat_type="D")
  lp.mace.D<-predict(cox.mod.cr.mace,newdata=newdata.D,type="lp",reference="zero")
  lp.death.D<-predict(cox.mod.cr.death,newdata=newdata.D,type="lp",reference="zero")
  
  #get cumulative hazard for each person (columns) at each time (rows)
  chaz.mace.A<-outer(cbhaz.mace$haz,exp(lp.mace.A))
  chaz.mace.B<-outer(cbhaz.mace$haz,exp(lp.mace.B))
  chaz.mace.C<-outer(cbhaz.mace$haz,exp(lp.mace.C))
  chaz.mace.D<-outer(cbhaz.mace$haz,exp(lp.mace.D))
  
  chaz.death.A<-outer(cbhaz.death$haz,exp(lp.death.A))
  chaz.death.B<-outer(cbhaz.death$haz,exp(lp.death.B))
  chaz.death.C<-outer(cbhaz.death$haz,exp(lp.death.C))
  chaz.death.D<-outer(cbhaz.death$haz,exp(lp.death.D))
  
  #get hazard for each person (columns) at each time (rows)
  haz.mace.A<-outer(bhaz.mace,exp(lp.mace.A))
  haz.mace.B<-outer(bhaz.mace,exp(lp.mace.B))
  haz.mace.C<-outer(bhaz.mace,exp(lp.mace.C))
  haz.mace.D<-outer(bhaz.mace,exp(lp.mace.D))
  
  haz.death.A<-outer(bhaz.death,exp(lp.death.A))
  haz.death.B<-outer(bhaz.death,exp(lp.death.B))
  haz.death.C<-outer(bhaz.death,exp(lp.death.C))
  haz.death.D<-outer(bhaz.death,exp(lp.death.D))
  
  #Obtain cumulative incidences for each person setting x=0
  
  surv.A<-exp(-chaz.mace.A-chaz.death.A)
  surv.B<-exp(-chaz.mace.B-chaz.death.B)
  surv.C<-exp(-chaz.mace.C-chaz.death.C)
  surv.D<-exp(-chaz.mace.D-chaz.death.D)
  
  cinc.mace.A.id<-sapply(1:dim(dta.analysis)[1],FUN=function(x){cumsum(haz.mace.A[,x]*surv.A[,x])})
  cinc.mace.B.id<-sapply(1:dim(dta.analysis)[1],FUN=function(x){cumsum(haz.mace.B[,x]*surv.B[,x])})
  cinc.mace.C.id<-sapply(1:dim(dta.analysis)[1],FUN=function(x){cumsum(haz.mace.C[,x]*surv.C[,x])})
  cinc.mace.D.id<-sapply(1:dim(dta.analysis)[1],FUN=function(x){cumsum(haz.mace.D[,x]*surv.D[,x])})
  
  cinc.death.A.id<-sapply(1:dim(dta.analysis)[1],FUN=function(x){cumsum(haz.death.A[,x]*surv.A[,x])})
  cinc.death.B.id<-sapply(1:dim(dta.analysis)[1],FUN=function(x){cumsum(haz.death.B[,x]*surv.B[,x])})
  cinc.death.C.id<-sapply(1:dim(dta.analysis)[1],FUN=function(x){cumsum(haz.death.C[,x]*surv.C[,x])})
  cinc.death.D.id<-sapply(1:dim(dta.analysis)[1],FUN=function(x){cumsum(haz.death.D[,x]*surv.D[,x])})
  
  #average over individuals
  cinc.mace.A.gform<-rowMeans(cinc.mace.A.id)
  cinc.mace.B.gform<-rowMeans(cinc.mace.B.id)
  cinc.mace.C.gform<-rowMeans(cinc.mace.C.id)
  cinc.mace.D.gform<-rowMeans(cinc.mace.D.id)
  
  cinc.death.A.gform<-rowMeans(cinc.death.A.id)
  cinc.death.B.gform<-rowMeans(cinc.death.B.id)
  cinc.death.C.gform<-rowMeans(cinc.death.C.id)
  cinc.death.D.gform<-rowMeans(cinc.death.D.id)
  
  cinc.times.gform<-cbhaz.mace$time/365.25
  
  #obtain cumulative incidences for each time in time.grid
  
  cincstep.A.mace.gform<-stepfun(cinc.times.gform,c(0,cinc.mace.A.gform))
  cincstep.B.mace.gform<-stepfun(cinc.times.gform,c(0,cinc.mace.B.gform))
  cincstep.C.mace.gform<-stepfun(cinc.times.gform,c(0,cinc.mace.C.gform))
  cincstep.D.mace.gform<-stepfun(cinc.times.gform,c(0,cinc.mace.D.gform))
  
  cincstep.A.death.gform<-stepfun(cinc.times.gform,c(0,cinc.death.A.gform))
  cincstep.B.death.gform<-stepfun(cinc.times.gform,c(0,cinc.death.B.gform))
  cincstep.C.death.gform<-stepfun(cinc.times.gform,c(0,cinc.death.C.gform))
  cincstep.D.death.gform<-stepfun(cinc.times.gform,c(0,cinc.death.D.gform))
  
  cinc.A.mace.gform.bs[b,]<-cincstep.A.mace.gform(time.grid)
  cinc.B.mace.gform.bs[b,]<-cincstep.B.mace.gform(time.grid)
  cinc.C.mace.gform.bs[b,]<-cincstep.C.mace.gform(time.grid)
  cinc.D.mace.gform.bs[b,]<-cincstep.D.mace.gform(time.grid)
  
  cinc.A.death.gform.bs[b,]<-cincstep.A.death.gform(time.grid)
  cinc.B.death.gform.bs[b,]<-cincstep.B.death.gform(time.grid)
  cinc.C.death.gform.bs[b,]<-cincstep.C.death.gform(time.grid)
  cinc.D.death.gform.bs[b,]<-cincstep.D.death.gform(time.grid)
  
  #obtain effect estimates
  
  diff.B.mace.gform.bs[b,]<-cinc.B.mace.gform-cinc.A.mace.gform
  diff.C.mace.gform.bs[b,]<-cinc.C.mace.gform-cinc.A.mace.gform
  diff.D.mace.gform.bs[b,]<-cinc.D.mace.gform-cinc.A.mace.gform
  
  diff.B.death.gform.bs[b,]<-cinc.B.death.gform-cinc.A.death.gform
  diff.C.death.gform.bs[b,]<-cinc.C.death.gform-cinc.A.death.gform
  diff.D.death.gform.bs[b,]<-cinc.D.death.gform-cinc.A.death.gform

}

#------------------------------------------------------
#bootstrap 95% CIs for cumulative incidences and differences (percentile method)
#------------------------------------------------------

cinc.A.mace.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.A.mace.gform.bs[,x],0.05)})
cinc.B.mace.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.B.mace.gform.bs[,x],0.05)})
cinc.C.mace.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.C.mace.gform.bs[,x],0.05)})
cinc.D.mace.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.D.mace.gform.bs[,x],0.05)})

cinc.A.mace.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.A.mace.gform.bs[,x],0.95)})
cinc.B.mace.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.B.mace.gform.bs[,x],0.95)})
cinc.C.mace.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.C.mace.gform.bs[,x],0.95)})
cinc.D.mace.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.D.mace.gform.bs[,x],0.95)})

cinc.A.death.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.A.death.gform.bs[,x],0.05)})
cinc.B.death.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.B.death.gform.bs[,x],0.05)})
cinc.C.death.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.C.death.gform.bs[,x],0.05)})
cinc.D.death.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.D.death.gform.bs[,x],0.05)})

cinc.A.death.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.A.death.gform.bs[,x],0.95)})
cinc.B.death.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.B.death.gform.bs[,x],0.95)})
cinc.C.death.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.C.death.gform.bs[,x],0.95)})
cinc.D.death.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(cinc.D.death.gform.bs[,x],0.95)})

diff.B.mace.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.B.mace.gform.bs[,x],0.05)})
diff.C.mace.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.C.mace.gform.bs[,x],0.05)})
diff.D.mace.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.D.mace.gform.bs[,x],0.05)})

diff.B.mace.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.B.mace.gform.bs[,x],0.95)})
diff.C.mace.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.C.mace.gform.bs[,x],0.95)})
diff.D.mace.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.D.mace.gform.bs[,x],0.95)})

diff.B.death.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.B.death.gform.bs[,x],0.05)})
diff.C.death.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.C.death.gform.bs[,x],0.05)})
diff.D.death.gform.CIL<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.D.death.gform.bs[,x],0.05)})

diff.B.death.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.B.death.gform.bs[,x],0.95)})
diff.C.death.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.C.death.gform.bs[,x],0.95)})
diff.D.death.gform.CIU<-sapply(1:length(time.grid),FUN=function(x){quantile(diff.D.death.gform.bs[,x],0.95)})


