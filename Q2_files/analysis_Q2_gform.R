#######################################################
#Q2: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D 
#or never using treatment, in people with hb>=7.5. 
#
#Using g-formula to handle baseline confounding, 
# and using IPCW to handle censoring due to deviation from non-treatment (for the 'never treat' strategy).
#######################################################


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
summary(mod.trtA.denom)

#Model for probability of starting treatment B in a given time period
mod.trtB.num<-glm(treat_status_B~1,
                  data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
mod.trtB.denom<-glm(treat_status_B~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
summary(mod.trtB.denom)

#Model for probability of starting treatment C in a given time period
mod.trtC.num<-glm(treat_status_C~1,
                  data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
mod.trtC.denom<-glm(treat_status_C~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
summary(mod.trtC.denom)

#Model for probability of starting treatment D in a given time period
mod.trtD.num<-glm(treat_status_D~1,
                  data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
mod.trtD.denom<-glm(treat_status_D~rcs(tstart,4)+sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                      hyp+dys+cvd+kidney+panc,
                    data=dta.analysis[dta.analysis$treat_type=="nevertrt",],family="binomial")
summary(mod.trtD.denom)

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
#G-formula/standardisation analysis
#------------------------------------------------------

#fit two cause-specific hazard models, stratified by treatment
#weighted by the IPCW
cox.mod.cr.mace<-coxph(Surv(tstart/365.25,tstop/365.25,event_type_cr==1)~strata(treat_type)+
                         sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                         hyp+dys+cvd+kidney+panc,data=dta.analysis2,weights = dta.analysis2$ipcw)
cox.mod.cr.death<-coxph(Surv(tstart/365.25,tstop/365.25,event_type_cr==2)~strata(treat_type)+
                          sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
                          hyp+dys+cvd+kidney+panc,data=dta.analysis2,weights = dta.analysis2$ipcw)

#get baseline cumulative hazards
cbhaz.mace<-basehaz(cox.mod.cr.mace,centered = F)
cbhaz.mace.A<-cbhaz.mace[cbhaz.mace$strata=="A",]
cbhaz.mace.B<-cbhaz.mace[cbhaz.mace$strata=="B",]
cbhaz.mace.C<-cbhaz.mace[cbhaz.mace$strata=="C",]
cbhaz.mace.D<-cbhaz.mace[cbhaz.mace$strata=="D",]
cbhaz.mace.nevertrt<-cbhaz.mace[cbhaz.mace$strata=="nevertrt",]

cbhaz.death<-basehaz(cox.mod.cr.death,centered = F)
cbhaz.death.A<-cbhaz.death[cbhaz.death$strata=="A",]
cbhaz.death.B<-cbhaz.death[cbhaz.death$strata=="B",]
cbhaz.death.C<-cbhaz.death[cbhaz.death$strata=="C",]
cbhaz.death.D<-cbhaz.death[cbhaz.death$strata=="D",]
cbhaz.death.nevertrt<-cbhaz.death[cbhaz.death$strata=="nevertrt",]

#get baseline hazards
bhaz.mace.A<-c(cbhaz.mace.A$haz[1],diff(cbhaz.mace.A$haz))
bhaz.mace.B<-c(cbhaz.mace.B$haz[1],diff(cbhaz.mace.B$haz))
bhaz.mace.C<-c(cbhaz.mace.C$haz[1],diff(cbhaz.mace.C$haz))
bhaz.mace.D<-c(cbhaz.mace.D$haz[1],diff(cbhaz.mace.D$haz))
bhaz.mace.nevertrt<-c(cbhaz.mace.nevertrt$haz[1],diff(cbhaz.mace.nevertrt$haz))

bhaz.death.A<-c(cbhaz.death.A$haz[1],diff(cbhaz.death.A$haz))
bhaz.death.B<-c(cbhaz.death.B$haz[1],diff(cbhaz.death.B$haz))
bhaz.death.C<-c(cbhaz.death.C$haz[1],diff(cbhaz.death.C$haz))
bhaz.death.D<-c(cbhaz.death.D$haz[1],diff(cbhaz.death.D$haz))
bhaz.death.nevertrt<-c(cbhaz.death.nevertrt$haz[1],diff(cbhaz.death.nevertrt$haz))

#get linear predictors 
#these are the same under each treatment 
#because we allow a stratified baseline haard by treatment group
lp.mace<-predict(cox.mod.cr.mace,newdata=dta.analysis2.row1,type="lp",reference="zero")
lp.death<-predict(cox.mod.cr.death,newdata=dta.analysis2.row1,type="lp",reference="zero")

#get cumulative hazard for each person (columns) at each time (rows)
#under each treatment strategy
chaz.mace.A<-outer(cbhaz.mace.A$haz,exp(lp.mace))
chaz.mace.B<-outer(cbhaz.mace.B$haz,exp(lp.mace))
chaz.mace.C<-outer(cbhaz.mace.C$haz,exp(lp.mace))
chaz.mace.D<-outer(cbhaz.mace.D$haz,exp(lp.mace))
chaz.mace.nevertrt<-outer(cbhaz.mace.nevertrt$haz,exp(lp.mace))

chaz.death.A<-outer(cbhaz.death.A$haz,exp(lp.death))
chaz.death.B<-outer(cbhaz.death.B$haz,exp(lp.death))
chaz.death.C<-outer(cbhaz.death.C$haz,exp(lp.death))
chaz.death.D<-outer(cbhaz.death.D$haz,exp(lp.death))
chaz.death.nevertrt<-outer(cbhaz.death.nevertrt$haz,exp(lp.death))

#get hazard for each person (columns) at each time (rows)
haz.mace.A<-outer(bhaz.mace.A,exp(lp.mace))
haz.mace.B<-outer(bhaz.mace.B,exp(lp.mace))
haz.mace.C<-outer(bhaz.mace.C,exp(lp.mace))
haz.mace.D<-outer(bhaz.mace.D,exp(lp.mace))
haz.mace.nevertrt<-outer(bhaz.mace.nevertrt,exp(lp.mace))

haz.death.A<-outer(bhaz.death.A,exp(lp.death))
haz.death.B<-outer(bhaz.death.B,exp(lp.death))
haz.death.C<-outer(bhaz.death.C,exp(lp.death))
haz.death.D<-outer(bhaz.death.D,exp(lp.death))
haz.death.nevertrt<-outer(bhaz.death.nevertrt,exp(lp.death))

#Obtain cumulative incidences for each person
surv.A<-exp(-chaz.mace.A-chaz.death.A)
surv.B<-exp(-chaz.mace.B-chaz.death.B)
surv.C<-exp(-chaz.mace.C-chaz.death.C)
surv.D<-exp(-chaz.mace.D-chaz.death.D)
surv.nevertrt<-exp(-chaz.mace.nevertrt-chaz.death.nevertrt)

cinc.mace.A.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.mace.A[,x]*surv.A[,x])})
cinc.mace.B.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.mace.B[,x]*surv.B[,x])})
cinc.mace.C.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.mace.C[,x]*surv.C[,x])})
cinc.mace.D.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.mace.D[,x]*surv.D[,x])})
cinc.mace.nevertrt.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.mace.nevertrt[,x]*surv.nevertrt[,x])})

cinc.death.A.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.death.A[,x]*surv.A[,x])})
cinc.death.B.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.death.B[,x]*surv.B[,x])})
cinc.death.C.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.death.C[,x]*surv.C[,x])})
cinc.death.D.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.death.D[,x]*surv.D[,x])})
cinc.death.nevertrt.id<-sapply(1:dim(dta.analysis2.row1)[1],FUN=function(x){cumsum(haz.death.nevertrt[,x]*surv.nevertrt[,x])})

#average over individuals
cinc.mace.A.gform<-rowMeans(cinc.mace.A.id)
cinc.mace.B.gform<-rowMeans(cinc.mace.B.id)
cinc.mace.C.gform<-rowMeans(cinc.mace.C.id)
cinc.mace.D.gform<-rowMeans(cinc.mace.D.id)
cinc.mace.nevertrt.gform<-rowMeans(cinc.mace.nevertrt.id)

cinc.death.A.gform<-rowMeans(cinc.death.A.id)
cinc.death.B.gform<-rowMeans(cinc.death.B.id)
cinc.death.C.gform<-rowMeans(cinc.death.C.id)
cinc.death.D.gform<-rowMeans(cinc.death.D.id)
cinc.death.nevertrt.gform<-rowMeans(cinc.death.nevertrt.id)

#obtain cumulative incidences for each time in time.grid

cinc.times.A.gform<-cbhaz.mace.A$time
cinc.times.B.gform<-cbhaz.mace.B$time
cinc.times.C.gform<-cbhaz.mace.C$time
cinc.times.D.gform<-cbhaz.mace.D$time
cinc.times.nevertrt.gform<-cbhaz.mace.nevertrt$time

cincstep.A.mace.gform<-stepfun(cinc.times.A.gform,c(0,cinc.mace.A.gform))
cincstep.B.mace.gform<-stepfun(cinc.times.B.gform,c(0,cinc.mace.B.gform))
cincstep.C.mace.gform<-stepfun(cinc.times.C.gform,c(0,cinc.mace.C.gform))
cincstep.D.mace.gform<-stepfun(cinc.times.D.gform,c(0,cinc.mace.D.gform))
cincstep.nevertrt.mace.gform<-stepfun(cinc.times.nevertrt.gform,c(0,cinc.mace.nevertrt.gform))

cincstep.A.death.gform<-stepfun(cinc.times.A.gform,c(0,cinc.death.A.gform))
cincstep.B.death.gform<-stepfun(cinc.times.B.gform,c(0,cinc.death.B.gform))
cincstep.C.death.gform<-stepfun(cinc.times.C.gform,c(0,cinc.death.C.gform))
cincstep.D.death.gform<-stepfun(cinc.times.D.gform,c(0,cinc.death.D.gform))
cincstep.nevertrt.death.gform<-stepfun(cinc.times.nevertrt.gform,c(0,cinc.death.nevertrt.gform))

cinc.A.mace.gform<-cincstep.A.mace.gform(time.grid)
cinc.B.mace.gform<-cincstep.B.mace.gform(time.grid)
cinc.C.mace.gform<-cincstep.C.mace.gform(time.grid)
cinc.D.mace.gform<-cincstep.D.mace.gform(time.grid)
cinc.nevertrt.mace.gform<-cincstep.nevertrt.mace.gform(time.grid)

cinc.A.death.gform<-cincstep.A.death.gform(time.grid)
cinc.B.death.gform<-cincstep.B.death.gform(time.grid)
cinc.C.death.gform<-cincstep.C.death.gform(time.grid)
cinc.D.death.gform<-cincstep.D.death.gform(time.grid)
cinc.nevertrt.death.gform<-cincstep.nevertrt.death.gform(time.grid)

#obtain effect estimates

diff.A.mace.gform<-cinc.A.mace.gform-cinc.nevertrt.mace.gform
diff.B.mace.gform<-cinc.B.mace.gform-cinc.nevertrt.mace.gform
diff.C.mace.gform<-cinc.C.mace.gform-cinc.nevertrt.mace.gform
diff.D.mace.gform<-cinc.D.mace.gform-cinc.nevertrt.mace.gform

diff.A.death.gform<-cinc.A.death.gform-cinc.nevertrt.death.gform
diff.B.death.gform<-cinc.B.death.gform-cinc.nevertrt.death.gform
diff.C.death.gform<-cinc.C.death.gform-cinc.nevertrt.death.gform
diff.D.death.gform<-cinc.D.death.gform-cinc.nevertrt.death.gform

