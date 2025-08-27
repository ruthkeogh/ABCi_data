#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#Using g-formula (i.e. regression adjustment followed by standardisation)
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

#generate factor variable for competing events
dta.analysis$event_type_cr<-ifelse(dta.analysis$event_type==3,0,dta.analysis$event_type)
dta.analysis$event_type_cr<-as.factor(dta.analysis$event_type_cr)

#------------------------------------------------------
#Analysis: standardisation, competing events
#This is quite slow, and suggest using the faster (more 'by hand') version below 
#------------------------------------------------------

# #fit adjusted cause-specific Cox model
# cox.mod.cr<-coxph(Surv(event_time_new,event_type_cr)~as.factor(treat_type)+
#                     sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
#                     hyp+dys+cvd+kidney+panc,data=dta.analysis,id=dta.analysis$id)
# 
# #obtain predictions under each treatment strategy
# newdta.A<-dta.analysis
# newdta.A$treat_type<-"A"
# cinc.A<-survfit(cox.mod.cr,newdata=newdta.A)
# 
# newdta.B<-dta.analysis
# newdta.B$treat_type<-"B"
# cinc.B<-survfit(cox.mod.cr,newdata=newdta.B)
# 
# newdta.C<-dta.analysis
# newdta.C$treat_type<-"C"
# cinc.C<-survfit(cox.mod.cr,newdata=newdta.C)
# 
# newdta.D<-dta.analysis
# newdta.D$treat_type<-"D"
# cinc.D<-survfit(cox.mod.cr,newdata=newdta.D)
# 
# #marginal cumulative incidences
# 
# cinc.A.times.gform<-cinc.A$time/365.25
# cinc.B.times.gform<-cinc.B$time/365.25
# cinc.C.times.gform<-cinc.C$time/365.25
# cinc.D.times.gform<-cinc.D$time/365.25
# 
# cinc.A.mace.gform<-rowMeans(cinc.A$pstate[,,2])
# cinc.B.mace.gform<-rowMeans(cinc.B$pstate[,,2])
# cinc.C.mace.gform<-rowMeans(cinc.C$pstate[,,2])
# cinc.D.mace.gform<-rowMeans(cinc.D$pstate[,,2])
# 
# cinc.A.death.gform<-rowMeans(cinc.A$pstate[,,3])
# cinc.B.death.gform<-rowMeans(cinc.B$pstate[,,3])
# cinc.C.death.gform<-rowMeans(cinc.C$pstate[,,3])
# cinc.D.death.gform<-rowMeans(cinc.D$pstate[,,3])

#------------------------------------------------------
#Analysis: standardisation, competing events
#more manual analysis, but faster
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

cinc.A.mace.gform<-cincstep.A.mace.gform(time.grid)
cinc.B.mace.gform<-cincstep.B.mace.gform(time.grid)
cinc.C.mace.gform<-cincstep.C.mace.gform(time.grid)
cinc.D.mace.gform<-cincstep.D.mace.gform(time.grid)

cinc.A.death.gform<-cincstep.A.death.gform(time.grid)
cinc.B.death.gform<-cincstep.B.death.gform(time.grid)
cinc.C.death.gform<-cincstep.C.death.gform(time.grid)
cinc.D.death.gform<-cincstep.D.death.gform(time.grid)

#obtain effect estimates

diff.B.mace.gform<-cinc.B.mace.gform-cinc.A.mace.gform
diff.C.mace.gform<-cinc.C.mace.gform-cinc.A.mace.gform
diff.D.mace.gform<-cinc.D.mace.gform-cinc.A.mace.gform

diff.B.death.gform<-cinc.B.death.gform-cinc.A.death.gform
diff.C.death.gform<-cinc.C.death.gform-cinc.A.death.gform
diff.D.death.gform<-cinc.D.death.gform-cinc.A.death.gform



