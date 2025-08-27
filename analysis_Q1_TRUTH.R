#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#
#Using data simulated as if everyone had followed treatment strategy A, B, C or D.
#The counterfactual data are generated in ABCi_Q1_genTRUTH.R
#######################################################

#------------------------------------------------------
#------------------------------------------------------
#DATA SET-UP
#------------------------------------------------------
#------------------------------------------------------

for(trt.fix in c("A","B","C","D")){
  print(trt.fix)
  
  source("ABCi_Q1_genTRUTH.R")

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
  
  #------------------------------------------------------
  #Analysis: unweighted KM
  #------------------------------------------------------
  
  dta.analysis$event_type_cr<-ifelse(dta.analysis$event_type==3,0,dta.analysis$event_type)
  dta.analysis$event_type_cr<-as.factor(dta.analysis$event_type_cr)
  
  cr.trt.true<-survfit(Surv(event_time_new,event_type_cr)~1,data=dta.analysis)
  
  if(trt.fix=="A"){
    cr.trtA.true<-cr.trt.true
  }else if(trt.fix=="B"){
    cr.trtB.true<-cr.trt.true
  }else if(trt.fix=="C"){
    cr.trtC.true<-cr.trt.true
  }else if(trt.fix=="D"){
    cr.trtD.true<-cr.trt.true
  }
}

#------------------------------------------------------
#Store results
#Cumulative incidences and 95% CIs for each time in time.grid
#------------------------------------------------------

#estimates

cinc.A.times.true<-cr.trtA.true$time/365.25
cinc.B.times.true<-cr.trtB.true$time/365.25
cinc.C.times.true<-cr.trtC.true$time/365.25
cinc.D.times.true<-cr.trtD.true$time/365.25

cincstep.A.mace.true<-stepfun(cr.trtA.true$time/365.25,c(0,cr.trtA.true$pstate[,2]))
cincstep.B.mace.true<-stepfun(cr.trtB.true$time/365.25,c(0,cr.trtB.true$pstate[,2]))
cincstep.C.mace.true<-stepfun(cr.trtC.true$time/365.25,c(0,cr.trtC.true$pstate[,2]))
cincstep.D.mace.true<-stepfun(cr.trtD.true$time/365.25,c(0,cr.trtD.true$pstate[,2]))

cincstep.A.death.true<-stepfun(cr.trtA.true$time/365.25,c(0,cr.trtA.true$pstate[,3]))
cincstep.B.death.true<-stepfun(cr.trtB.true$time/365.25,c(0,cr.trtB.true$pstate[,3]))
cincstep.C.death.true<-stepfun(cr.trtC.true$time/365.25,c(0,cr.trtC.true$pstate[,3]))
cincstep.D.death.true<-stepfun(cr.trtD.true$time/365.25,c(0,cr.trtD.true$pstate[,3]))

cinc.A.mace.true<-cincstep.A.mace.true(time.grid)
cinc.B.mace.true<-cincstep.B.mace.true(time.grid)
cinc.C.mace.true<-cincstep.C.mace.true(time.grid)
cinc.D.mace.true<-cincstep.D.mace.true(time.grid)

cinc.A.death.true<-cincstep.A.death.true(time.grid)
cinc.B.death.true<-cincstep.B.death.true(time.grid)
cinc.C.death.true<-cincstep.C.death.true(time.grid)
cinc.D.death.true<-cincstep.D.death.true(time.grid)

diff.B.mace.true<-cinc.B.mace.true-cinc.A.mace.true
diff.C.mace.true<-cinc.D.mace.true-cinc.A.mace.true
diff.D.mace.true<-cinc.D.mace.true-cinc.A.mace.true

diff.B.death.true<-cinc.B.death.true-cinc.A.death.true
diff.C.death.true<-cinc.D.death.true-cinc.A.death.true
diff.D.death.true<-cinc.D.death.true-cinc.A.death.true

#95% CIs

CILstep.A.mace.true<-stepfun(cr.trtA.true$time/365.25,c(0,cr.trtA.true$lower[,2]))
CILstep.B.mace.true<-stepfun(cr.trtB.true$time/365.25,c(0,cr.trtB.true$lower[,2]))
CILstep.C.mace.true<-stepfun(cr.trtC.true$time/365.25,c(0,cr.trtC.true$lower[,2]))
CILstep.D.mace.true<-stepfun(cr.trtD.true$time/365.25,c(0,cr.trtD.true$lower[,2]))

CIUstep.A.mace.true<-stepfun(cr.trtA.true$time/365.25,c(0,cr.trtA.true$upper[,2]))
CIUstep.B.mace.true<-stepfun(cr.trtB.true$time/365.25,c(0,cr.trtB.true$upper[,2]))
CIUstep.C.mace.true<-stepfun(cr.trtC.true$time/365.25,c(0,cr.trtC.true$upper[,2]))
CIUstep.D.mace.true<-stepfun(cr.trtD.true$time/365.25,c(0,cr.trtD.true$upper[,2]))

CILstep.A.death.true<-stepfun(cr.trtA.true$time/365.25,c(0,cr.trtA.true$lower[,3]))
CILstep.B.death.true<-stepfun(cr.trtB.true$time/365.25,c(0,cr.trtB.true$lower[,3]))
CILstep.C.death.true<-stepfun(cr.trtC.true$time/365.25,c(0,cr.trtC.true$lower[,3]))
CILstep.D.death.true<-stepfun(cr.trtD.true$time/365.25,c(0,cr.trtD.true$lower[,3]))

CIUstep.A.death.true<-stepfun(cr.trtA.true$time/365.25,c(0,cr.trtA.true$upper[,3]))
CIUstep.B.death.true<-stepfun(cr.trtB.true$time/365.25,c(0,cr.trtB.true$upper[,3]))
CIUstep.C.death.true<-stepfun(cr.trtC.true$time/365.25,c(0,cr.trtC.true$upper[,3]))
CIUstep.D.death.true<-stepfun(cr.trtD.true$time/365.25,c(0,cr.trtD.true$upper[,3]))

cinc.A.mace.true.CIL<-CILstep.A.mace.true(time.grid)
cinc.B.mace.true.CIL<-CILstep.B.mace.true(time.grid)
cinc.C.mace.true.CIL<-CILstep.C.mace.true(time.grid)
cinc.D.mace.true.CIL<-CILstep.D.mace.true(time.grid)

cinc.A.mace.true.CIU<-CIUstep.A.mace.true(time.grid)
cinc.B.mace.true.CIU<-CIUstep.B.mace.true(time.grid)
cinc.C.mace.true.CIU<-CIUstep.C.mace.true(time.grid)
cinc.D.mace.true.CIU<-CIUstep.D.mace.true(time.grid)

cinc.A.death.true.CIL<-CILstep.A.death.true(time.grid)
cinc.B.death.true.CIL<-CILstep.B.death.true(time.grid)
cinc.C.death.true.CIL<-CILstep.C.death.true(time.grid)
cinc.D.death.true.CIL<-CILstep.D.death.true(time.grid)

cinc.A.death.true.CIU<-CIUstep.A.death.true(time.grid)
cinc.B.death.true.CIU<-CIUstep.B.death.true(time.grid)
cinc.C.death.true.CIU<-CIUstep.C.death.true(time.grid)
cinc.D.death.true.CIU<-CIUstep.D.death.true(time.grid)
