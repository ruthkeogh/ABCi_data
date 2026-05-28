#######################################################
#Q3: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D with a grace period of 90 days
#in people with hb>=7.5. 
#
#Using the clone-censor-weight approach
#######################################################

#------------------------------------------------------
#------------------------------------------------------
#Fit weights models to handle artificial censoring in each data set
#The censoring occurs due to starting different treatments, and these depend on covariates in different ways
#Hence we fit separate models for censoring due to starting different treatments
#------------------------------------------------------
#------------------------------------------------------

#generate lagged treatment status
dta.analysis<-dta.analysis%>%group_by(id)%>%mutate(lag_any_treat=lag(treat_status_A,default=0)+
                                                     lag(treat_status_B,default=0)+lag(treat_status_C,default=0)+
                                                     +lag(treat_status_D,default=0))

#restrict to rows up to and including treatment initiation and within the grace period
dta.ipcw<-dta.analysis[dta.analysis$lag_any_treat==0 & dta.analysis$tstart<=gp,]

#Model for probability of starting treatment A within the grace period or in the next period (using tstart<=gp is important, instead of tstart<gp)
mod.trtA.num<-glm(treat_status_A~1,
                  data=dta.ipcw,family="binomial")
mod.trtA.denom<-glm(treat_status_A~as.factor(tstart)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.ipcw,family="binomial")

#Model for probability of starting treatment B before the end of the grace period
mod.trtB.num<-glm(treat_status_B~1,
                  data=dta.ipcw,family="binomial")
mod.trtB.denom<-glm(treat_status_B~as.factor(tstart)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.ipcw,family="binomial")

#Model for probability of starting treatment C before the end of the grace period
mod.trtC.denom<-glm(treat_status_C~1,
                    data=dta.ipcw,family="binomial")
mod.trtC.denom<-glm(treat_status_C~as.factor(tstart)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.ipcw,family="binomial")

#Model for probability of starting treatment D before the end of the grace period
mod.trtD.num<-glm(treat_status_D~1,
                  data=dta.ipcw,family="binomial")
mod.trtD.denom<-glm(treat_status_D~as.factor(tstart)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.ipcw,family="binomial")

#------------------------------------------------------
#Obtain IPCW for strategy A
#------------------------------------------------------

dta.trtA<-dta.trtA%>%group_by(id)%>%mutate(lag_treat_status_A=lag(treat_status_A,default=0))

dta.trtA<-dta.trtA%>%group_by(id)%>%mutate(lag_any_treat=lag(treat_status_A,default=0)+
                                             lag(treat_status_B,default=0)+lag(treat_status_C,default=0)+
                                             +lag(treat_status_D,default=0))

#probability of remaining uncensored at each time (conditional on being observed at the start of the interval)
#Here, this means not starting another treatment up to and including tstart=gp, and then starting A in the next period

dta.trtA$nocens.B.num<-1
dta.trtA$nocens.C.num<-1
dta.trtA$nocens.D.num<-1
dta.trtA$nocens.B.num[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0]<-1-predict(mod.trtB.num,newdata=dta.trtA[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0,],type="response")
dta.trtA$nocens.C.num[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0]<-1-predict(mod.trtC,newdata=dta.trtA[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0,],type="response")
dta.trtA$nocens.D.num[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0]<-1-predict(mod.trtD.num,newdata=dta.trtA[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0,],type="response")
dta.trtA$nocens.othertrt.num<-dta.trtA$nocens.B.num*dta.trtA$nocens.C.num*dta.trtA$nocens.D.num

dta.trtA$nocens.B.denom<-1
dta.trtA$nocens.C.denom<-1
dta.trtA$nocens.D.denom<-1
dta.trtA$nocens.B.denom[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0]<-1-predict(mod.trtB.denom,newdata=dta.trtA[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0,],type="response")
dta.trtA$nocens.C.denom[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0]<-1-predict(mod.trtC,newdata=dta.trtA[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0,],type="response")
dta.trtA$nocens.D.denom[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0]<-1-predict(mod.trtD.denom,newdata=dta.trtA[dta.trtA$tstart<gp & dta.trtA$lag_any_treat==0,],type="response")
dta.trtA$nocens.othertrt.denom<-dta.trtA$nocens.B.denom*dta.trtA$nocens.C.denom*dta.trtA$nocens.D.denom

dta.trtA$startA.num<-1
dta.trtA$startA.num[dta.trtA$tstart==gp & dta.trtA$lag_treat_status_A==0]<-predict(mod.trtA.num,newdata=dta.trtA[dta.trtA$tstart==gp & dta.trtA$lag_treat_status_A==0,],type="response")

dta.trtA$startA.denom<-1
dta.trtA$startA.denom[dta.trtA$tstart==gp & dta.trtA$lag_treat_status_A==0]<-predict(mod.trtA.denom,newdata=dta.trtA[dta.trtA$tstart==gp & dta.trtA$lag_treat_status_A==0,],type="response")

dta.trtA$nocens.prob.num<-dta.trtA$startA.num*dta.trtA$nocens.othertrt.num
dta.trtA<-dta.trtA%>%group_by(id)%>%mutate(nocens.prob.num.prod=cumprod(nocens.prob.num))

dta.trtA$nocens.prob.denom<-dta.trtA$startA.denom*dta.trtA$nocens.othertrt.denom
dta.trtA<-dta.trtA%>%group_by(id)%>%mutate(nocens.prob.denom.prod=cumprod(nocens.prob.denom))

#ipcw
dta.trtA<-dta.trtA%>%group_by(id)%>%mutate(ipcw=nocens.prob.num.prod/nocens.prob.denom.prod)

#------------------------------------------------------
#Obtain IPCW for strategy B
#------------------------------------------------------

dta.trtB<-dta.trtB%>%group_by(id)%>%mutate(lag_treat_status_B=lag(treat_status_B,default=0))

dta.trtB<-dta.trtB%>%group_by(id)%>%mutate(lag_any_treat=lag(treat_status_A,default=0)+
                                             lag(treat_status_B,default=0)+lag(treat_status_C,default=0)+
                                             +lag(treat_status_D,default=0))

#probability of remaining uncensored at each time (conditional on being observed at the start of the interval)
#Here, this means not starting another treatment up to and including tstart=gp, and then starting A in the next period

dta.trtB$nocens.A.num<-1
dta.trtB$nocens.C.num<-1
dta.trtB$nocens.D.num<-1
dta.trtB$nocens.A.num[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0]<-1-predict(mod.trtA.num,newdata=dta.trtB[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0,],type="response")
dta.trtB$nocens.C.num[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0]<-1-predict(mod.trtC.num,newdata=dta.trtB[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0,],type="response")
dta.trtB$nocens.D.num[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0]<-1-predict(mod.trtD.num,newdata=dta.trtB[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0,],type="response")
dta.trtB$nocens.othertrt.num<-dta.trtB$nocens.A.num*dta.trtB$nocens.C.num*dta.trtB$nocens.D.num

dta.trtB$startB.num<-1
dta.trtB$startB.num[dta.trtB$tstart==gp & dta.trtB$lag_treat_status_B==0]<-predict(mod.trtB.num,newdata=dta.trtB[dta.trtB$tstart==gp & dta.trtB$lag_treat_status_B==0,],type="response")

dta.trtB$nocens.prob.num<-dta.trtB$startB.num*dta.trtB$nocens.othertrt.num
dta.trtB<-dta.trtB%>%group_by(id)%>%mutate(nocens.prob.num.prod=cumprod(nocens.prob.num))

dta.trtB$nocens.A.denom<-1
dta.trtB$nocens.C.denom<-1
dta.trtB$nocens.D.denom<-1
dta.trtB$nocens.A.denom[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0]<-1-predict(mod.trtA.denom,newdata=dta.trtB[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0,],type="response")
dta.trtB$nocens.C.denom[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0]<-1-predict(mod.trtC.denom,newdata=dta.trtB[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0,],type="response")
dta.trtB$nocens.D.denom[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0]<-1-predict(mod.trtD.denom,newdata=dta.trtB[dta.trtB$tstart<gp & dta.trtB$lag_any_treat==0,],type="response")
dta.trtB$nocens.othertrt.denom<-dta.trtB$nocens.A.denom*dta.trtB$nocens.C.denom*dta.trtB$nocens.D.denom

dta.trtB$startB.denom<-1
dta.trtB$startB.denom[dta.trtB$tstart==gp & dta.trtB$lag_treat_status_B==0]<-predict(mod.trtB.denom,newdata=dta.trtB[dta.trtB$tstart==gp & dta.trtB$lag_treat_status_B==0,],type="response")

dta.trtB$nocens.prob.denom<-dta.trtB$startB.denom*dta.trtB$nocens.othertrt.denom
dta.trtB<-dta.trtB%>%group_by(id)%>%mutate(nocens.prob.denom.prod=cumprod(nocens.prob.denom))

#ipcw
dta.trtB<-dta.trtB%>%group_by(id)%>%mutate(ipcw=nocens.prob.num.prod/nocens.prob.denom.prod)

#------------------------------------------------------
#Obtain IPCW for strategy C
#------------------------------------------------------

dta.trtC<-dta.trtC%>%group_by(id)%>%mutate(lag_treat_status_C=lag(treat_status_C,default=0))

dta.trtC<-dta.trtC%>%group_by(id)%>%mutate(lag_any_treat=lag(treat_status_A,default=0)+
                                             lag(treat_status_B,default=0)+lag(treat_status_C,default=0)+
                                             +lag(treat_status_D,default=0))

#probability of remaining uncensored at each time (conditional on being observed at the start of the interval)
#Here, this means not starting another treatment up to and including tstart=gp, and then starting A in the next period

dta.trtC$nocens.A.num<-1
dta.trtC$nocens.B.num<-1
dta.trtC$nocens.D.num<-1
dta.trtC$nocens.A.num[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0]<-1-predict(mod.trtA.num,newdata=dta.trtC[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0,],type="response")
dta.trtC$nocens.B.num[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0]<-1-predict(mod.trtB.num,newdata=dta.trtC[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0,],type="response")
dta.trtC$nocens.D.num[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0]<-1-predict(mod.trtD.num,newdata=dta.trtC[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0,],type="response")
dta.trtC$nocens.othertrt.num<-dta.trtC$nocens.A.num*dta.trtC$nocens.B.num*dta.trtC$nocens.D.num

dta.trtC$startC.num<-1
dta.trtC$startC.num[dta.trtC$tstart==gp & dta.trtC$lag_treat_status_C==0]<-predict(mod.trtC.num,newdata=dta.trtC[dta.trtC$tstart==gp & dta.trtC$lag_treat_status_C==0,],type="response")

dta.trtC$nocens.prob.num<-dta.trtC$startC.num*dta.trtC$nocens.othertrt.num
dta.trtC<-dta.trtC%>%group_by(id)%>%mutate(nocens.prob.num.prod=cumprod(nocens.prob.num))

dta.trtC$nocens.A.denom<-1
dta.trtC$nocens.B.denom<-1
dta.trtC$nocens.D.denom<-1
dta.trtC$nocens.A.denom[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0]<-1-predict(mod.trtA.denom,newdata=dta.trtC[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0,],type="response")
dta.trtC$nocens.B.denom[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0]<-1-predict(mod.trtB.denom,newdata=dta.trtC[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0,],type="response")
dta.trtC$nocens.D.denom[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0]<-1-predict(mod.trtD.denom,newdata=dta.trtC[dta.trtC$tstart<gp & dta.trtC$lag_any_treat==0,],type="response")
dta.trtC$nocens.othertrt.denom<-dta.trtC$nocens.A.denom*dta.trtC$nocens.B.denom*dta.trtC$nocens.D.denom

dta.trtC$startC.denom<-1
dta.trtC$startC.denom[dta.trtC$tstart==gp & dta.trtC$lag_treat_status_C==0]<-predict(mod.trtC.denom,newdata=dta.trtC[dta.trtC$tstart==gp & dta.trtC$lag_treat_status_C==0,],type="response")

dta.trtC$nocens.prob.denom<-dta.trtC$startC.denom*dta.trtC$nocens.othertrt.denom
dta.trtC<-dta.trtC%>%group_by(id)%>%mutate(nocens.prob.denom.prod=cumprod(nocens.prob.denom))

#ipcw
dta.trtC<-dta.trtC%>%group_by(id)%>%mutate(ipcw=nocens.prob.num.prod/nocens.prob.denom.prod)

#------------------------------------------------------
#Obtain IPCW for strategy D
#------------------------------------------------------

dta.trtD<-dta.trtD%>%group_by(id)%>%mutate(lag_treat_status_D=lag(treat_status_D,default=0))

dta.trtD<-dta.trtD%>%group_by(id)%>%mutate(lag_any_treat=lag(treat_status_A,default=0)+
                                             lag(treat_status_B,default=0)+lag(treat_status_C,default=0)+
                                             +lag(treat_status_D,default=0))

#probability of remaining uncensored at each time (conditional on being observed at the start of the interval)
#Here, this means not starting another treatment up to and including tstart=gp, and then starting A in the next period

dta.trtD$nocens.A.num<-1
dta.trtD$nocens.B.num<-1
dta.trtD$nocens.C.num<-1
dta.trtD$nocens.A.num[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0]<-1-predict(mod.trtA.num,newdata=dta.trtD[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0,],type="response")
dta.trtD$nocens.B.num[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0]<-1-predict(mod.trtB.num,newdata=dta.trtD[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0,],type="response")
dta.trtD$nocens.C.num[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0]<-1-predict(mod.trtC.num,newdata=dta.trtD[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0,],type="response")
dta.trtD$nocens.othertrt.num<-dta.trtD$nocens.A.num*dta.trtD$nocens.B.num*dta.trtD$nocens.C.num

dta.trtD$startD.num<-1
dta.trtD$startD.num[dta.trtD$tstart==gp & dta.trtD$lag_treat_status_D==0]<-predict(mod.trtC.num,newdata=dta.trtD[dta.trtD$tstart==gp & dta.trtD$lag_treat_status_D==0,],type="response")

dta.trtD$nocens.prob.num<-dta.trtD$startD.num*dta.trtD$nocens.othertrt.num
dta.trtD<-dta.trtD%>%group_by(id)%>%mutate(nocens.prob.num.prod=cumprod(nocens.prob.num))

dta.trtD$nocens.A.denom<-1
dta.trtD$nocens.B.denom<-1
dta.trtD$nocens.C.denom<-1
dta.trtD$nocens.A.denom[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0]<-1-predict(mod.trtA.denom,newdata=dta.trtD[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0,],type="response")
dta.trtD$nocens.B.denom[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0]<-1-predict(mod.trtB.denom,newdata=dta.trtD[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0,],type="response")
dta.trtD$nocens.C.denom[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0]<-1-predict(mod.trtC.denom,newdata=dta.trtD[dta.trtD$tstart<gp & dta.trtD$lag_any_treat==0,],type="response")
dta.trtD$nocens.othertrt.denom<-dta.trtD$nocens.A.denom*dta.trtD$nocens.B.denom*dta.trtD$nocens.C.denom

dta.trtD$startD.denom<-1
dta.trtD$startD.denom[dta.trtD$tstart==gp & dta.trtD$lag_treat_status_D==0]<-predict(mod.trtC.denom,newdata=dta.trtD[dta.trtD$tstart==gp & dta.trtD$lag_treat_status_D==0,],type="response")

dta.trtD$nocens.prob.denom<-dta.trtD$startD.denom*dta.trtD$nocens.othertrt.denom
dta.trtD<-dta.trtD%>%group_by(id)%>%mutate(nocens.prob.denom.prod=cumprod(nocens.prob.denom))

#ipcw
dta.trtD<-dta.trtD%>%group_by(id)%>%mutate(ipcw=nocens.prob.num.prod/nocens.prob.denom.prod)

#------------------------------------------------------
#Obtain IPCW for 'never treat' strategy
#------------------------------------------------------

#Model for probability of starting treatment A in a given time period
mod.trtA.num<-glm(treat_status_A~1,
                  data=dta.analysis[dta.analysis$lag_any_treat==0,],family="binomial")
mod.trtA.denom<-glm(treat_status_A~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.analysis[dta.analysis$lag_any_treat==0,],family="binomial")

#Model for probability of starting treatment B in a given time period
mod.trtB.num<-glm(treat_status_B~1,
                  data=dta.analysis[dta.analysis$lag_any_treat==0,],family="binomial")
mod.trtB.denom<-glm(treat_status_B~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.analysis[dta.analysis$lag_any_treat==0,],family="binomial")

#Model for probability of starting treatment C in a given time period
mod.trtC.num<-glm(treat_status_C~1,
                  data=dta.analysis[dta.analysis$lag_any_treat==0,],family="binomial")
mod.trtC.denom<-glm(treat_status_C~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.analysis[dta.analysis$lag_any_treat==0,],family="binomial")

#Model for probability of starting treatment D in a given time period
mod.trtD.num<-glm(treat_status_D~1,
                  data=dta.analysis[dta.analysis$lag_any_treat==0,],family="binomial")
mod.trtD.denom<-glm(treat_status_D~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.analysis[dta.analysis$lag_any_treat==0,],family="binomial")

#probability of remaining uncensored in a given time period
dta.nevertrt$num.A<-1-predict(mod.trtA.num,newdata=dta.nevertrt,type="response")
dta.nevertrt$num.B<-1-predict(mod.trtB.num,newdata=dta.nevertrt,type="response")
dta.nevertrt$num.C<-1-predict(mod.trtC.num,newdata=dta.nevertrt,type="response")
dta.nevertrt$num.D<-1-predict(mod.trtD.num,newdata=dta.nevertrt,type="response")
dta.nevertrt$denom.A<-1-predict(mod.trtA.denom,newdata=dta.nevertrt,type="response")
dta.nevertrt$denom.B<-1-predict(mod.trtB.denom,newdata=dta.nevertrt,type="response")
dta.nevertrt$denom.C<-1-predict(mod.trtC.denom,newdata=dta.nevertrt,type="response")
dta.nevertrt$denom.D<-1-predict(mod.trtD.denom,newdata=dta.nevertrt,type="response")
dta.nevertrt$num.anytrt<-dta.nevertrt$num.A*dta.nevertrt$num.B*dta.nevertrt$num.C*dta.nevertrt$num.D
dta.nevertrt$denom.anytrt<-dta.nevertrt$denom.A*dta.nevertrt$denom.B*dta.nevertrt$denom.C*dta.nevertrt$denom.D

#ipcw
dta.nevertrt<-dta.nevertrt%>%group_by(id)%>%mutate(ipcw=cumprod(num.anytrt/denom.anytrt))

dta.nevertrt$tstop.new<-dta.nevertrt$tstop

#------------------------------------------------------
#analysis: IPCW Aalen-Johansen estimates
#------------------------------------------------------

cr.trtA.ccw<-survfit(Surv(tstart,tstop.new,event_type_cr)~1,data=dta.trtA,weights=dta.trtA$ipcw,id=dta.trtA$id)
cr.trtB.ccw<-survfit(Surv(tstart,tstop.new,event_type_cr)~1,data=dta.trtB,weights=dta.trtB$ipcw,id=dta.trtB$id)
cr.trtC.ccw<-survfit(Surv(tstart,tstop.new,event_type_cr)~1,data=dta.trtC,weights=dta.trtC$ipcw,id=dta.trtC$id)
cr.trtD.ccw<-survfit(Surv(tstart,tstop.new,event_type_cr)~1,data=dta.trtD,weights=dta.trtD$ipcw,id=dta.trtD$id)
cr.nevertrt.ccw<-survfit(Surv(tstart,tstop.new,event_type_cr)~1,data=dta.nevertrt,weights=dta.nevertrt$ipcw,id=dta.nevertrt$id)

#obtain cumulative incidences for each time in time.grid

cinc.A.times.ccw<-cr.trtA.ccw$time/365.25
cinc.B.times.ccw<-cr.trtB.ccw$time/365.25
cinc.C.times.ccw<-cr.trtC.ccw$time/365.25
cinc.D.times.ccw<-cr.trtD.ccw$time/365.25
cinc.nevertrt.times.ccw<-cr.nevertrt.ccw$time/365.25

cincstep.A.mace.ccw<-stepfun(cr.trtA.ccw$time/365.25,c(0,cr.trtA.ccw$pstate[,2]))
cincstep.B.mace.ccw<-stepfun(cr.trtB.ccw$time/365.25,c(0,cr.trtB.ccw$pstate[,2]))
cincstep.C.mace.ccw<-stepfun(cr.trtC.ccw$time/365.25,c(0,cr.trtC.ccw$pstate[,2]))
cincstep.D.mace.ccw<-stepfun(cr.trtD.ccw$time/365.25,c(0,cr.trtD.ccw$pstate[,2]))
cincstep.nevertrt.mace.ccw<-stepfun(cr.nevertrt.ccw$time/365.25,c(0,cr.nevertrt.ccw$pstate[,2]))

cincstep.A.death.ccw<-stepfun(cr.trtA.ccw$time/365.25,c(0,cr.trtA.ccw$pstate[,3]))
cincstep.B.death.ccw<-stepfun(cr.trtB.ccw$time/365.25,c(0,cr.trtB.ccw$pstate[,3]))
cincstep.C.death.ccw<-stepfun(cr.trtC.ccw$time/365.25,c(0,cr.trtC.ccw$pstate[,3]))
cincstep.D.death.ccw<-stepfun(cr.trtD.ccw$time/365.25,c(0,cr.trtD.ccw$pstate[,3]))
cincstep.nevertrt.death.ccw<-stepfun(cr.nevertrt.ccw$time/365.25,c(0,cr.nevertrt.ccw$pstate[,3]))

cinc.A.mace.ccw<-cincstep.A.mace.ccw(time.grid)
cinc.B.mace.ccw<-cincstep.B.mace.ccw(time.grid)
cinc.C.mace.ccw<-cincstep.C.mace.ccw(time.grid)
cinc.D.mace.ccw<-cincstep.D.mace.ccw(time.grid)
cinc.nevertrt.mace.ccw<-cincstep.nevertrt.mace.ccw(time.grid)

cinc.A.death.ccw<-cincstep.A.death.ccw(time.grid)
cinc.B.death.ccw<-cincstep.B.death.ccw(time.grid)
cinc.C.death.ccw<-cincstep.C.death.ccw(time.grid)
cinc.D.death.ccw<-cincstep.D.death.ccw(time.grid)
cinc.nevertrt.death.ccw<-cincstep.nevertrt.mace.ccw(time.grid)

