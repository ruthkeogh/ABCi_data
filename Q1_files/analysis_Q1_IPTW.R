#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#Using IPTW
#######################################################

#------------------------------------------------------
#Estimate the IPTW
#------------------------------------------------------

#estimate weights
mod.trtA.num<-glm(treat_status_A~1,data=dta.analysis,family="binomial")
mod.trtA.denom<-glm(treat_status_A~sex+age+smok_former+smok_current+diabdur+bmi+hb+
                      hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")
# mod.trtA.denom<-glm(treat_status_A~sex+age+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
#                       hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial") #with spline terms

mod.trtB.num<-glm(treat_status_B~1,data=dta.analysis,family="binomial")
mod.trtB.denom<-glm(treat_status_B~sex+age+smok_former+smok_current+diabdur+bmi+hb+
                      hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")
# mod.trtB.denom<-glm(treat_status_B~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
#                       hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")#with spline terms

mod.trtC.num<-glm(treat_status_C~1,data=dta.analysis,family="binomial")
mod.trtC.denom<-glm(treat_status_C~sex+age+smok_former+smok_current+diabdur+bmi+hb+
                      hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")
# mod.trtC.denom<-glm(treat_status_C~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
#                       hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")#with spline terms

mod.trtD.num<-glm(treat_status_D~1,data=dta.analysis,family="binomial")
mod.trtD.denom<-glm(treat_status_D~sex+age+smok_former+smok_current+diabdur+bmi+hb+
                      hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")
# mod.trtD.denom<-glm(treat_status_D~sex+rcs(age,4)+smok_former+smok_current+rcs(diabdur,4)+rcs(bmi,4)+rcs(hb,4)+
#                       hyp+dys+cvd+kidney+panc,data=dta.analysis,family="binomial")#with spline terms

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
#Estimate cumulative incidences using IPTW
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

cinc.A.mace.iptw<-cincstep.A.mace.iptw(time.grid)
cinc.B.mace.iptw<-cincstep.B.mace.iptw(time.grid)
cinc.C.mace.iptw<-cincstep.C.mace.iptw(time.grid)
cinc.D.mace.iptw<-cincstep.D.mace.iptw(time.grid)

cinc.A.death.iptw<-cincstep.A.death.iptw(time.grid)
cinc.B.death.iptw<-cincstep.B.death.iptw(time.grid)
cinc.C.death.iptw<-cincstep.C.death.iptw(time.grid)
cinc.D.death.iptw<-cincstep.D.death.iptw(time.grid)

#obtain effect estimates

diff.B.mace.iptw<-cinc.B.mace.iptw-cinc.A.mace.iptw
diff.C.mace.iptw<-cinc.C.mace.iptw-cinc.A.mace.iptw
diff.D.mace.iptw<-cinc.D.mace.iptw-cinc.A.mace.iptw

diff.B.death.iptw<-cinc.B.death.iptw-cinc.A.death.iptw
diff.C.death.iptw<-cinc.C.death.iptw-cinc.A.death.iptw
diff.D.death.iptw<-cinc.D.death.iptw-cinc.A.death.iptw

