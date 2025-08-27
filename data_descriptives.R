###############################################################
# Descriptive statistics for the ABCi data
###############################################################

#------------------------------------------------------
#packages
#------------------------------------------------------

library(tableone)
library(janitor)
library(survival)

#----------------
#table of baseline characteristics
#----------------

table1<-CreateTableOne(data=dta[dta$tstart==0,],
                       vars=c("sex","age","smok","diabdur","bmi","hb","hyp","dys","cvd","kidney","panc"),
                       factorVars = c("sex","smok","hyp","dys","cvd","kidney","panc"))
print(table1,showAllLevels = T)

#----------------
#summary of time-dependent covariates at different time points
#----------------

summ.tcov<-matrix(nrow=7,ncol=6)

#summarize at 0,1,2,3,4,5 years (closest time prior to this as variables are recorded every 30 days)
summ.times<-c(0,360,720,1080,1440,1800)

for(k in 1:6){
  n.hyp<-table(dta$hyp[dta$tstart==summ.times[k]])
  n.dys<-table(dta$dys[dta$tstart==summ.times[k]])
  n.cvd<-table(dta$cvd[dta$tstart==summ.times[k]])
  n.kidney<-table(dta$kidney[dta$tstart==summ.times[k]])
  n.panc<-table(dta$panc[dta$tstart==summ.times[k]])
  
  p.hyp<-round(100*prop.table(table(dta$hyp[dta$tstart==summ.times[k]])),1)
  p.dys<-round(100*prop.table(table(dta$dys[dta$tstart==summ.times[k]])))
  p.cvd<-round(100*prop.table(table(dta$cvd[dta$tstart==summ.times[k]])))
  p.kidney<-round(100*prop.table(table(dta$kidney[dta$tstart==summ.times[k]])))
  p.panc<-round(100*prop.table(table(dta$panc[dta$tstart==summ.times[k]])))
  
  mean.bmi<-round(mean(dta$bmi[dta$tstart==summ.times[k]]),1)
  mean.hb<-round(mean(dta$hb[dta$tstart==summ.times[k]]),1)
  
  sd.bmi<-round(sd(dta$bmi[dta$tstart==summ.times[k]]),1)
  sd.hb<-round(sd(dta$hb[dta$tstart==summ.times[k]]),1)
  
  summ.tcov[,k]<-c(paste0(mean.bmi," (",sd.bmi,")"),
                   paste0(mean.hb," (",sd.hb,")"),
                   paste0(n.hyp," (",p.hyp,"%)")[2],
                   paste0(n.dys," (",p.dys,"%)")[2],
                   paste0(n.cvd," (",p.cvd,"%)")[2],
                   paste0(n.kidney," (",p.kidney,"%)")[2],
                   paste0(n.panc," (",p.panc,"%)")[2])
}

rownames(summ.tcov)<-c("bmi","hb","hyp","dys","cvd","kidney","panc")
summ.tcov

#----------------
#plots of HbA1c over time: shown for 100 individuals
#----------------

plot(dta$tstart[dta$id==1],dta$hb[dta$id==1],type="l",ylim=c(6,12),xlim=c(0,2000))
for(i in 2:100){
points(dta$tstart[dta$id==i],dta$hb[dta$id==i],type="l")
}

#----------------
#plots of BMI over time: shown for 100 individuals
#----------------

plot(dta$tstart[dta$id==1],dta$bmi[dta$id==1],type="l",ylim=c(15,60),xlim=c(0,2000))
for(i in 2:100){
  points(dta$tstart[dta$id==i],dta$bmi[dta$id==i],type="l")
}

#----------------
#number of people who start each treatment
#----------------

dta<-dta%>%group_by(id)%>%mutate(rownum=row_number())%>%mutate(lastrow=max(rownum))
tabyl(dta$treat_status_A[dta$rownum==dta$lastrow])
tabyl(dta$treat_status_B[dta$rownum==dta$lastrow])
tabyl(dta$treat_status_C[dta$rownum==dta$lastrow])
tabyl(dta$treat_status_D[dta$rownum==dta$lastrow])

#----------------
#summary of treatment status at different time points
#----------------

dta$treat_status_any<-0
dta$treat_status_any<-ifelse(dta$treat_status_A==1|
                               dta$treat_status_B==1|
                               dta$treat_status_C==1|
                               dta$treat_status_D==1,1,dta$treat_status_any)

summ.treat<-matrix(nrow=5,ncol=6)

for(k in 1:6){
  n.any<-table(dta$treat_status_any[dta$tstart==summ.times[k]])
  n.A<-table(dta$treat_status_A[dta$tstart==summ.times[k]])
  n.B<-table(dta$treat_status_B[dta$tstart==summ.times[k]])
  n.C<-table(dta$treat_status_C[dta$tstart==summ.times[k]])
  n.D<-table(dta$treat_status_D[dta$tstart==summ.times[k]])

  p.any<-round(100*prop.table(table(dta$treat_status_any[dta$tstart==summ.times[k]])),1)
  p.A<-round(100*prop.table(table(dta$treat_status_A[dta$tstart==summ.times[k]])),1)
  p.B<-round(100*prop.table(table(dta$treat_status_B[dta$tstart==summ.times[k]])),1)
  p.C<-round(100*prop.table(table(dta$treat_status_C[dta$tstart==summ.times[k]])),1)
  p.D<-round(100*prop.table(table(dta$treat_status_D[dta$tstart==summ.times[k]])),1)

  summ.treat[,k]<-c(paste0(n.any," (",p.any,"%)")[2],
                    paste0(n.A," (",p.A,"%)")[2],
                    paste0(n.B," (",p.B,"%)")[2],
                    paste0(n.C," (",p.C,"%)")[2],
                    paste0(n.D," (",p.D,"%)")[2])
}

rownames(summ.treat)<-c("Any","A","B","C","D")
summ.treat

#----------------
#cumulative incidence plots for MACE and death
#----------------

dta.events<-dta[dta$tstart==0,c("id","tstart","tstop",
                   "mace_time","death_time","cens_time","mace_status",
                   "death_status","cens_status","event_time","event_type")]
tabyl(dta.events$event_type)

dta.events$event_type<-ifelse(dta.events$event_type==3,0,dta.events$event_type)
tabyl(dta.events$event_type)

dta.events$event<-factor(dta.events$event_type,0:2, labels=c("censor", "mace", "death"))

#cumulative incidences for MACE and death

cinc<-survfit(Surv(event_time,event)~1,data=dta.events,id=dta.events$id)

plot(cinc, col=c(1,1), lty=c(1,2),
     mark.time=FALSE, lwd=2,xscale=365.25,
     xlab="Time (years)", ylab="Cumulative incidence",ylim=c(0,0.30))
legend("topleft",c("MACE", "Other deaths"),
       col=c(1,1), lty=c(1,2), lwd=2, bty='n')

#cumulative incidences at time 5: censoring, mace, death
cinc$pstate[dim(cinc$pstate)[1],]

