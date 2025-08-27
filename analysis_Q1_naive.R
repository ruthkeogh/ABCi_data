#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#Using a naive analysis, with no control for baseline confounders
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
#Analysis: unweighted
#------------------------------------------------------

cr.trtA.naive<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis[dta.analysis$treat_type=="A",])
cr.trtB.naive<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis[dta.analysis$treat_type=="B",])
cr.trtC.naive<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis[dta.analysis$treat_type=="C",])
cr.trtD.naive<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis[dta.analysis$treat_type=="D",])

#------------------------------------------------------
#Store results
#Cumulative incidences and 95% CIs for each time in time.grid
#------------------------------------------------------

cincstep.A.mace.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$pstate[,2]))
cincstep.B.mace.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$pstate[,2]))
cincstep.C.mace.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$pstate[,2]))
cincstep.D.mace.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$pstate[,2]))

cincstep.A.death.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$pstate[,3]))
cincstep.B.death.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$pstate[,3]))
cincstep.C.death.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$pstate[,3]))
cincstep.D.death.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$pstate[,3]))

cinc.A.mace.naive<-cincstep.A.mace.naive(time.grid)
cinc.B.mace.naive<-cincstep.B.mace.naive(time.grid)
cinc.C.mace.naive<-cincstep.C.mace.naive(time.grid)
cinc.D.mace.naive<-cincstep.D.mace.naive(time.grid)

cinc.A.death.naive<-cincstep.A.death.naive(time.grid)
cinc.B.death.naive<-cincstep.B.death.naive(time.grid)
cinc.C.death.naive<-cincstep.C.death.naive(time.grid)
cinc.D.death.naive<-cincstep.D.death.naive(time.grid)

#95% CIs

CILstep.A.mace.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$lower[,2]))
CILstep.B.mace.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$lower[,2]))
CILstep.C.mace.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$lower[,2]))
CILstep.D.mace.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$lower[,2]))

CIUstep.A.mace.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$upper[,2]))
CIUstep.B.mace.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$upper[,2]))
CIUstep.C.mace.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$upper[,2]))
CIUstep.D.mace.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$upper[,2]))

CILstep.A.death.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$lower[,3]))
CILstep.B.death.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$lower[,3]))
CILstep.C.death.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$lower[,3]))
CILstep.D.death.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$lower[,3]))

CIUstep.A.death.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$upper[,3]))
CIUstep.B.death.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$upper[,3]))
CIUstep.C.death.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$upper[,3]))
CIUstep.D.death.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$upper[,3]))

cinc.A.mace.naive.CIL<-CILstep.A.mace.naive(time.grid)
cinc.B.mace.naive.CIL<-CILstep.B.mace.naive(time.grid)
cinc.C.mace.naive.CIL<-CILstep.C.mace.naive(time.grid)
cinc.D.mace.naive.CIL<-CILstep.D.mace.naive(time.grid)

cinc.A.mace.naive.CIU<-CIUstep.A.mace.naive(time.grid)
cinc.B.mace.naive.CIU<-CIUstep.B.mace.naive(time.grid)
cinc.C.mace.naive.CIU<-CIUstep.C.mace.naive(time.grid)
cinc.D.mace.naive.CIU<-CIUstep.D.mace.naive(time.grid)

cinc.A.death.naive.CIL<-CILstep.A.death.naive(time.grid)
cinc.B.death.naive.CIL<-CILstep.B.death.naive(time.grid)
cinc.C.death.naive.CIL<-CILstep.C.death.naive(time.grid)
cinc.D.death.naive.CIL<-CILstep.D.death.naive(time.grid)

cinc.A.death.naive.CIU<-CIUstep.A.death.naive(time.grid)
cinc.B.death.naive.CIU<-CIUstep.B.death.naive(time.grid)
cinc.C.death.naive.CIU<-CIUstep.C.death.naive(time.grid)
cinc.D.death.naive.CIU<-CIUstep.D.death.naive(time.grid)
