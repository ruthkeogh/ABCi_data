#######################################################
#Q2: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D 
#or never using treatment, in people with hb>=7.5. 
#
#Using a naive analysis, with no control for baseline confounders, and no IPCW
#######################################################

#first row of dta.analysis, which is the population for which we will obtain average treatment effects
dta.analysis2.row1<-dta.analysis2[dta.analysis2$tstart==0,]

#drop last row of data in dta.nevertrt (this was just needed for estimation of the weights)
dta.analysis2<-dta.analysis[dta.analysis$treat_type!="nevertrt"|
                              dta.analysis$treat_type=="nevertrt" & dta.analysis$treat_status_any==0,]

#------------------------------------------------------
#Unweighted Aalen-Johansen estimates
#------------------------------------------------------

cr.trtA.naive<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                      data=dta.analysis2[dta.analysis2$treat_type=="A",],
                      id=dta.analysis2$id_new[dta.analysis2$treat_type=="A"])

cr.trtB.naive<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                      data=dta.analysis2[dta.analysis2$treat_type=="B",],
                      id=dta.analysis2$id_new[dta.analysis2$treat_type=="B"])

cr.trtC.naive<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                      data=dta.analysis2[dta.analysis2$treat_type=="C",],
                      id=dta.analysis2$id_new[dta.analysis2$treat_type=="C"])

cr.trtD.naive<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                      data=dta.analysis2[dta.analysis2$treat_type=="D",],
                      id=dta.analysis2$id_new[dta.analysis2$treat_type=="D"])

cr.nevertrt.naive<-survfit(Surv(tstart,tstop,event_type_cr)~1,
                          data=dta.analysis2[dta.analysis2$treat_type=="nevertrt",],
                          id=dta.analysis2$id_new[dta.analysis2$treat_type=="nevertrt"])

#------------------------------------------------------
#Store results
#Cumulative incidences and 95% CIs for each time in time.grid
#------------------------------------------------------

cincstep.A.mace.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$pstate[,2]))
cincstep.B.mace.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$pstate[,2]))
cincstep.C.mace.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$pstate[,2]))
cincstep.D.mace.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$pstate[,2]))
cincstep.nevertrt.mace.naive<-stepfun(cr.nevertrt.naive$time/365.25,c(0,cr.nevertrt.naive$pstate[,2]))

cincstep.A.death.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$pstate[,3]))
cincstep.B.death.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$pstate[,3]))
cincstep.C.death.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$pstate[,3]))
cincstep.D.death.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$pstate[,3]))
cincstep.nevertrt.death.naive<-stepfun(cr.nevertrt.naive$time/365.25,c(0,cr.nevertrt.naive$pstate[,3]))

cinc.A.mace.naive<-cincstep.A.mace.naive(time.grid)
cinc.B.mace.naive<-cincstep.B.mace.naive(time.grid)
cinc.C.mace.naive<-cincstep.C.mace.naive(time.grid)
cinc.D.mace.naive<-cincstep.D.mace.naive(time.grid)
cinc.nevertrt.mace.naive<-cincstep.nevertrt.mace.naive(time.grid)

cinc.A.death.naive<-cincstep.A.death.naive(time.grid)
cinc.B.death.naive<-cincstep.B.death.naive(time.grid)
cinc.C.death.naive<-cincstep.C.death.naive(time.grid)
cinc.D.death.naive<-cincstep.D.death.naive(time.grid)
cinc.nevertrt.death.naive<-cincstep.nevertrt.death.naive(time.grid)

#95% CIs

CILstep.A.mace.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$lower[,2]))
CILstep.B.mace.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$lower[,2]))
CILstep.C.mace.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$lower[,2]))
CILstep.D.mace.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$lower[,2]))
CILstep.nevertrt.mace.naive<-stepfun(cr.nevertrt.naive$time/365.25,c(0,cr.nevertrt.naive$lower[,2]))

CIUstep.A.mace.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$upper[,2]))
CIUstep.B.mace.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$upper[,2]))
CIUstep.C.mace.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$upper[,2]))
CIUstep.D.mace.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$upper[,2]))
CIUstep.nevertrt.mace.naive<-stepfun(cr.nevertrt.naive$time/365.25,c(0,cr.nevertrt.naive$upper[,2]))

CILstep.A.death.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$lower[,3]))
CILstep.B.death.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$lower[,3]))
CILstep.C.death.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$lower[,3]))
CILstep.D.death.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$lower[,3]))
CILstep.nevertrt.death.naive<-stepfun(cr.nevertrt.naive$time/365.25,c(0,cr.nevertrt.naive$lower[,3]))

CIUstep.A.death.naive<-stepfun(cr.trtA.naive$time/365.25,c(0,cr.trtA.naive$upper[,3]))
CIUstep.B.death.naive<-stepfun(cr.trtB.naive$time/365.25,c(0,cr.trtB.naive$upper[,3]))
CIUstep.C.death.naive<-stepfun(cr.trtC.naive$time/365.25,c(0,cr.trtC.naive$upper[,3]))
CIUstep.D.death.naive<-stepfun(cr.trtD.naive$time/365.25,c(0,cr.trtD.naive$upper[,3]))
CIUstep.nevertrt.death.naive<-stepfun(cr.nevertrt.naive$time/365.25,c(0,cr.nevertrt.naive$upper[,3]))

cinc.A.mace.naive.CIL<-CILstep.A.mace.naive(time.grid)
cinc.B.mace.naive.CIL<-CILstep.B.mace.naive(time.grid)
cinc.C.mace.naive.CIL<-CILstep.C.mace.naive(time.grid)
cinc.D.mace.naive.CIL<-CILstep.D.mace.naive(time.grid)
cinc.nevertrt.mace.naive.CIL<-CILstep.nevertrt.mace.naive(time.grid)

cinc.A.mace.naive.CIU<-CIUstep.A.mace.naive(time.grid)
cinc.B.mace.naive.CIU<-CIUstep.B.mace.naive(time.grid)
cinc.C.mace.naive.CIU<-CIUstep.C.mace.naive(time.grid)
cinc.D.mace.naive.CIU<-CIUstep.D.mace.naive(time.grid)
cinc.nevertrt.mace.naive.CIU<-CIUstep.nevertrt.mace.naive(time.grid)

cinc.A.death.naive.CIL<-CILstep.A.death.naive(time.grid)
cinc.B.death.naive.CIL<-CILstep.B.death.naive(time.grid)
cinc.C.death.naive.CIL<-CILstep.C.death.naive(time.grid)
cinc.D.death.naive.CIL<-CILstep.D.death.naive(time.grid)
cinc.nevertrt.death.naive.CIL<-CILstep.nevertrt.death.naive(time.grid)

cinc.A.death.naive.CIU<-CIUstep.A.death.naive(time.grid)
cinc.B.death.naive.CIU<-CIUstep.B.death.naive(time.grid)
cinc.C.death.naive.CIU<-CIUstep.C.death.naive(time.grid)
cinc.D.death.naive.CIU<-CIUstep.D.death.naive(time.grid)
cinc.nevertrt.death.naive.CIU<-CIUstep.nevertrt.death.naive(time.grid)
