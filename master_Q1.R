#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#
#This master file for Q1 runs the following analyses:
#(i) Naive analysis, (ii) IPTW, (iii) g-formula, (iv) TMLE.
#It also generates the true cumulative incidences under the 4 treatment strategies. 
#The results are plotted and tabulated.
#######################################################

#------------------------------------------------------
#packages
#------------------------------------------------------

library(survival)
library(rms)
library(tidyverse)
library(lmtp)
library(truncnorm)
library(ggplot2)
library(gridExtra)
library(tidyr)

#------------------------------------------------------
#functions and parameters
#------------------------------------------------------

expit<-function(x){exp(x)/(1+exp(x))}

#time grid on which cumulative incidences are stored
time.grid<-seq(0,5,0.1)

#number of bootstrap samples 
nboot<-1000

#------------------------------------------------------
#generate observed data
#------------------------------------------------------

source("ABCi_data_generate.R")

#------------------------------------------------------
#Analysis: naive (unadjusted for confounding)
#------------------------------------------------------

#estimates and 95% CIs
source("analysis_Q1_naive.R")

#------------------------------------------------------
#Analysis: IPTW
#------------------------------------------------------

#estimates
source("analysis_Q1_IPTW.R")

#bootstrap to get 95% CIs
source("analysis_Q1_IPTW_bs.R")

#------------------------------------------------------
#Analysis: standardisation/gformula
#------------------------------------------------------

#estimates
source("analysis_Q1_gform.R")

#bootstrap to get 95% CIs
source("analysis_Q1_gform_bs.R")

#------------------------------------------------------
#Analysis: doubly robust
#------------------------------------------------------

#estimates and 95% CIs
source("analysis_Q1_lmtp.R")

#------------------------------------------------------
#Obtain true values
#------------------------------------------------------

source("analysis_Q1_TRUTH.R")


#------------------------------------------------------
#Plot results: each method on separate plot (also showing truth)
#------------------------------------------------------

#---
#data sets for plotting

plot.dta.mace<-data.frame(
  time.grid,
  #true
  cinc.A.mace.true,cinc.B.mace.true,cinc.C.mace.true,cinc.D.mace.true,
  cinc.A.mace.true.CIL,cinc.B.mace.true.CIL,cinc.C.mace.true.CIL,cinc.D.mace.true.CIL,
  cinc.A.mace.true.CIU,cinc.B.mace.true.CIU,cinc.C.mace.true.CIU,cinc.D.mace.true.CIU,
  #iptw
  cinc.A.mace.naive,cinc.B.mace.naive,cinc.C.mace.naive,cinc.D.mace.naive,
  cinc.A.mace.naive.CIL,cinc.B.mace.naive.CIL,cinc.C.mace.naive.CIL,cinc.D.mace.naive.CIL,
  cinc.A.mace.naive.CIU,cinc.B.mace.naive.CIU,cinc.C.mace.naive.CIU,cinc.D.mace.naive.CIU,
  #iptw
  cinc.A.mace.iptw,cinc.B.mace.iptw,cinc.C.mace.iptw,cinc.D.mace.iptw,
  cinc.A.mace.iptw.CIL,cinc.B.mace.iptw.CIL,cinc.C.mace.iptw.CIL,cinc.D.mace.iptw.CIL,
  cinc.A.mace.iptw.CIU,cinc.B.mace.iptw.CIU,cinc.C.mace.iptw.CIU,cinc.D.mace.iptw.CIU,
  #gformula
  cinc.A.mace.gform,cinc.B.mace.gform,cinc.C.mace.gform,cinc.D.mace.gform,
  cinc.A.mace.gform.CIL,cinc.B.mace.gform.CIL,cinc.C.mace.gform.CIL,cinc.D.mace.gform.CIL,
  cinc.A.mace.gform.CIU,cinc.B.mace.gform.CIU,cinc.C.mace.gform.CIU,cinc.D.mace.gform.CIU,
  #lmtp
  cinc.A.mace.lmtp,cinc.B.mace.lmtp,cinc.C.mace.lmtp,cinc.D.mace.lmtp,
  cinc.A.mace.lmtp.CIL,cinc.B.mace.lmtp.CIL,cinc.C.mace.lmtp.CIL,cinc.D.mace.lmtp.CIL,
  cinc.A.mace.lmtp.CIU,cinc.B.mace.lmtp.CIU,cinc.C.mace.lmtp.CIU,cinc.D.mace.lmtp.CIU
)

plot.dta.death<-data.frame(
  time.grid,
  #true
  cinc.A.death.true,cinc.B.death.true,cinc.C.death.true,cinc.D.death.true,
  cinc.A.death.true.CIL,cinc.B.death.true.CIL,cinc.C.death.true.CIL,cinc.D.death.true.CIL,
  cinc.A.death.true.CIU,cinc.B.death.true.CIU,cinc.C.death.true.CIU,cinc.D.death.true.CIU,
  #iptw
  cinc.A.death.naive,cinc.B.death.naive,cinc.C.death.naive,cinc.D.death.naive,
  cinc.A.death.naive.CIL,cinc.B.death.naive.CIL,cinc.C.death.naive.CIL,cinc.D.death.naive.CIL,
  cinc.A.death.naive.CIU,cinc.B.death.naive.CIU,cinc.C.death.naive.CIU,cinc.D.death.naive.CIU,
  #iptw
  cinc.A.death.iptw,cinc.B.death.iptw,cinc.C.death.iptw,cinc.D.death.iptw,
  cinc.A.death.iptw.CIL,cinc.B.death.iptw.CIL,cinc.C.death.iptw.CIL,cinc.D.death.iptw.CIL,
  cinc.A.death.iptw.CIU,cinc.B.death.iptw.CIU,cinc.C.death.iptw.CIU,cinc.D.death.iptw.CIU,
  #gformula
  cinc.A.death.gform,cinc.B.death.gform,cinc.C.death.gform,cinc.D.death.gform,
  cinc.A.death.gform.CIL,cinc.B.death.gform.CIL,cinc.C.death.gform.CIL,cinc.D.death.gform.CIL,
  cinc.A.death.gform.CIU,cinc.B.death.gform.CIU,cinc.C.death.gform.CIU,cinc.D.death.gform.CIU,
  #lmtp
  cinc.A.death.lmtp,cinc.B.death.lmtp,cinc.C.death.lmtp,cinc.D.death.lmtp,
  cinc.A.death.lmtp.CIL,cinc.B.death.lmtp.CIL,cinc.C.death.lmtp.CIL,cinc.D.death.lmtp.CIL,
  cinc.A.death.lmtp.CIU,cinc.B.death.lmtp.CIU,cinc.C.death.lmtp.CIU,cinc.D.death.lmtp.CIU
)

#---
#plot set-up

addlinetoplot <- function(dataset, varx, vary,vcol,vlty) { 
  list(
    geom_step(data=dataset, aes_string(x=varx, y=vary,colour=vcol,linetype=vlty),linewidth=1) 
  )
}

cols=c("true"="#000000","A"="#D55E00","B"="#009E73","C"="#CC79A7","D"="#56B4E9")

ltys=c("true"="solid","other"="dashed")

addribbontoplot <- function(dataset, varymin, varymax,vcol,valpha) { 
  list(
    geom_ribbon(data=dataset, aes_string(ymin=varymin, ymax=varymax),fill=vcol,alpha=valpha)
  )
}

#---
#mace plots

plot.mace.iptw=ggplot(plot.dta.mace,aes(x=time.grid,y=cinc.A.mace.true))+theme_bw()+
  #true
  addlinetoplot(plot.dta.mace,"time.grid","cinc.A.mace.true",vcol='"A"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.B.mace.true",vcol='"B"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.C.mace.true",vcol='"C"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.D.mace.true",vcol='"D"',vlty='"true"')+
  #iptw - estimates
  addlinetoplot(plot.dta.mace,"time.grid","cinc.A.mace.iptw",vcol='"A"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.B.mace.iptw",vcol='"B"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.C.mace.iptw",vcol='"C"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.D.mace.iptw",vcol='"D"',vlty='"other"')+
  #iptw - CIs
  addribbontoplot(plot.dta.mace,varymin="cinc.A.mace.iptw.CIL",varymax="cinc.A.mace.iptw.CIU",vcol="#D55E00", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.B.mace.iptw.CIL",varymax="cinc.B.mace.iptw.CIU",vcol="#009E73", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.C.mace.iptw.CIL",varymax="cinc.C.mace.iptw.CIU",vcol="#CC79A7", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.D.mace.iptw.CIL",varymax="cinc.D.mace.iptw.CIU",vcol="#56B4E9", valpha=0.2)+
  #axes
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.5))+
  scale_colour_manual(NULL,values=cols,labels=c(A="Treatment A",B="Treatment B",C="Treatment C",D="Treatment D"), breaks=c("A","B","C","D"))+
  scale_linetype_manual(NULL,values=ltys,labels=c("True","IPTW"), breaks=c("solid","dashed"))+
  ylab("Cumulative incidence")+xlab("Time (years)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))+
  theme(legend.position="inside",legend.position.inside =c(0.2,0.7))+
  ggtitle("Estimation method: IPTW")

plot.mace.gform=ggplot(plot.dta.mace,aes(x=time.grid,y=cinc.A.mace.true))+theme_bw()+
  #true
  addlinetoplot(plot.dta.mace,"time.grid","cinc.A.mace.true",vcol='"A"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.B.mace.true",vcol='"B"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.C.mace.true",vcol='"C"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.D.mace.true",vcol='"D"',vlty='"true"')+
  #gform - estimates
  addlinetoplot(plot.dta.mace,"time.grid","cinc.A.mace.gform",vcol='"A"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.B.mace.gform",vcol='"B"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.C.mace.gform",vcol='"C"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.D.mace.gform",vcol='"D"',vlty='"other"')+
  #gform - CIs
  addribbontoplot(plot.dta.mace,varymin="cinc.A.mace.gform.CIL",varymax="cinc.A.mace.gform.CIU",vcol="#D55E00", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.B.mace.gform.CIL",varymax="cinc.B.mace.gform.CIU",vcol="#009E73", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.C.mace.gform.CIL",varymax="cinc.C.mace.gform.CIU",vcol="#CC79A7", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.D.mace.gform.CIL",varymax="cinc.D.mace.gform.CIU",vcol="#56B4E9", valpha=0.2)+
  #axes
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.5))+
  scale_colour_manual(NULL,values=cols,labels=c(A="Treatment A",B="Treatment B",C="Treatment C",D="Treatment D"), breaks=c("A","B","C","D"))+
  scale_linetype_manual(NULL,values=ltys,labels=c("True","gform"), breaks=c("solid","dashed"))+
  ylab("Cumulative incidence")+xlab("Time (years)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))+
  theme(legend.position="inside",legend.position.inside =c(0.2,0.7))+
  ggtitle("Estimation method: G-formula")

plot.mace.naive=ggplot(plot.dta.mace,aes(x=time.grid,y=cinc.A.mace.true))+theme_bw()+
  #true
  addlinetoplot(plot.dta.mace,"time.grid","cinc.A.mace.true",vcol='"A"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.B.mace.true",vcol='"B"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.C.mace.true",vcol='"C"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.D.mace.true",vcol='"D"',vlty='"true"')+
  #naive - estimates
  addlinetoplot(plot.dta.mace,"time.grid","cinc.A.mace.naive",vcol='"A"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.B.mace.naive",vcol='"B"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.C.mace.naive",vcol='"C"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.D.mace.naive",vcol='"D"',vlty='"other"')+
  #naive - CIs
  addribbontoplot(plot.dta.mace,varymin="cinc.A.mace.naive.CIL",varymax="cinc.A.mace.naive.CIU",vcol="#D55E00", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.B.mace.naive.CIL",varymax="cinc.B.mace.naive.CIU",vcol="#009E73", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.C.mace.naive.CIL",varymax="cinc.C.mace.naive.CIU",vcol="#CC79A7", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.D.mace.naive.CIL",varymax="cinc.D.mace.naive.CIU",vcol="#56B4E9", valpha=0.2)+
  #axes
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.5))+
  scale_colour_manual(NULL,values=cols,labels=c(A="Treatment A",B="Treatment B",C="Treatment C",D="Treatment D"), breaks=c("A","B","C","D"))+
  scale_linetype_manual(NULL,values=ltys,labels=c("True","naive"), breaks=c("solid","dashed"))+
  ylab("Cumulative incidence")+xlab("Time (years)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))+
  theme(legend.position="inside",legend.position.inside =c(0.2,0.7))+
  ggtitle("Estimation method: Naive")

plot.mace.lmtp=ggplot(plot.dta.mace,aes(x=time.grid,y=cinc.A.mace.true))+theme_bw()+
  #true
  addlinetoplot(plot.dta.mace,"time.grid","cinc.A.mace.true",vcol='"A"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.B.mace.true",vcol='"B"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.C.mace.true",vcol='"C"',vlty='"true"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.D.mace.true",vcol='"D"',vlty='"true"')+
  #lmtp - estimates
  addlinetoplot(plot.dta.mace,"time.grid","cinc.A.mace.lmtp",vcol='"A"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.B.mace.lmtp",vcol='"B"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.C.mace.lmtp",vcol='"C"',vlty='"other"')+
  addlinetoplot(plot.dta.mace,"time.grid","cinc.D.mace.lmtp",vcol='"D"',vlty='"other"')+
  #lmtp - CIs
  addribbontoplot(plot.dta.mace,varymin="cinc.A.mace.lmtp.CIL",varymax="cinc.A.mace.lmtp.CIU",vcol="#D55E00", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.B.mace.lmtp.CIL",varymax="cinc.B.mace.lmtp.CIU",vcol="#009E73", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.C.mace.lmtp.CIL",varymax="cinc.C.mace.lmtp.CIU",vcol="#CC79A7", valpha=0.2)+
  addribbontoplot(plot.dta.mace,varymin="cinc.D.mace.lmtp.CIL",varymax="cinc.D.mace.lmtp.CIU",vcol="#56B4E9", valpha=0.2)+
  #axes
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.5))+
  scale_colour_manual(NULL,values=cols,labels=c(A="Treatment A",B="Treatment B",C="Treatment C",D="Treatment D"), breaks=c("A","B","C","D"))+
  scale_linetype_manual(NULL,values=ltys,labels=c("True","lmtp"), breaks=c("solid","dashed"))+
  ylab("Cumulative incidence")+xlab("Time (years)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))+
  theme(legend.position="inside",legend.position.inside =c(0.2,0.7))+
  ggtitle("Estimation method: Doubly robust")


windows(10,7)
grid.arrange(plot.mace.naive,plot.mace.iptw,plot.mace.gform,plot.mace.lmtp,ncol=2)

#---
#death plots

plot.death.iptw=ggplot(plot.dta.death,aes(x=time.grid,y=cinc.A.death.true))+theme_bw()+
  #true
  addlinetoplot(plot.dta.death,"time.grid","cinc.A.death.true",vcol='"A"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.B.death.true",vcol='"B"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.C.death.true",vcol='"C"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.D.death.true",vcol='"D"',vlty='"true"')+
  #iptw - estimates
  addlinetoplot(plot.dta.death,"time.grid","cinc.A.death.iptw",vcol='"A"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.B.death.iptw",vcol='"B"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.C.death.iptw",vcol='"C"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.D.death.iptw",vcol='"D"',vlty='"other"')+
  #iptw - CIs
  addribbontoplot(plot.dta.death,varymin="cinc.A.death.iptw.CIL",varymax="cinc.A.death.iptw.CIU",vcol="#D55E00", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.B.death.iptw.CIL",varymax="cinc.B.death.iptw.CIU",vcol="#009E73", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.C.death.iptw.CIL",varymax="cinc.C.death.iptw.CIU",vcol="#CC79A7", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.D.death.iptw.CIL",varymax="cinc.D.death.iptw.CIU",vcol="#56B4E9", valpha=0.2)+
  #axes
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,0.10,0.05),limits=c(0,0.10))+
  scale_colour_manual(NULL,values=cols,labels=c(A="Treatment A",B="Treatment B",C="Treatment C",D="Treatment D"), breaks=c("A","B","C","D"))+
  scale_linetype_manual(NULL,values=ltys,labels=c("True","IPTW"), breaks=c("solid","dashed"))+
  ylab("Cumulative incidence")+xlab("Time (years)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))+
  theme(legend.position="inside",legend.position.inside =c(0.2,0.7))+
  ggtitle("Estimation method: IPTW")

plot.death.gform=ggplot(plot.dta.death,aes(x=time.grid,y=cinc.A.death.true))+theme_bw()+
  #true
  addlinetoplot(plot.dta.death,"time.grid","cinc.A.death.true",vcol='"A"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.B.death.true",vcol='"B"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.C.death.true",vcol='"C"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.D.death.true",vcol='"D"',vlty='"true"')+
  #gform - estimates
  addlinetoplot(plot.dta.death,"time.grid","cinc.A.death.gform",vcol='"A"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.B.death.gform",vcol='"B"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.C.death.gform",vcol='"C"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.D.death.gform",vcol='"D"',vlty='"other"')+
  #gform - CIs
  addribbontoplot(plot.dta.death,varymin="cinc.A.death.gform.CIL",varymax="cinc.A.death.gform.CIU",vcol="#D55E00", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.B.death.gform.CIL",varymax="cinc.B.death.gform.CIU",vcol="#009E73", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.C.death.gform.CIL",varymax="cinc.C.death.gform.CIU",vcol="#CC79A7", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.D.death.gform.CIL",varymax="cinc.D.death.gform.CIU",vcol="#56B4E9", valpha=0.2)+
  #axes
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,0.10,0.05),limits=c(0,0.10))+
  scale_colour_manual(NULL,values=cols,labels=c(A="Treatment A",B="Treatment B",C="Treatment C",D="Treatment D"), breaks=c("A","B","C","D"))+
  scale_linetype_manual(NULL,values=ltys,labels=c("True","gform"), breaks=c("solid","dashed"))+
  ylab("Cumulative incidence")+xlab("Time (years)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))+
  theme(legend.position="inside",legend.position.inside =c(0.3,0.7))+
  ggtitle("Estimation method: G-formula")

plot.death.naive=ggplot(plot.dta.death,aes(x=time.grid,y=cinc.A.death.true))+theme_bw()+
  #true
  addlinetoplot(plot.dta.death,"time.grid","cinc.A.death.true",vcol='"A"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.B.death.true",vcol='"B"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.C.death.true",vcol='"C"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.D.death.true",vcol='"D"',vlty='"true"')+
  #naive - estimates
  addlinetoplot(plot.dta.death,"time.grid","cinc.A.death.naive",vcol='"A"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.B.death.naive",vcol='"B"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.C.death.naive",vcol='"C"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.D.death.naive",vcol='"D"',vlty='"other"')+
  #naive - CIs
  addribbontoplot(plot.dta.death,varymin="cinc.A.death.naive.CIL",varymax="cinc.A.death.naive.CIU",vcol="#D55E00", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.B.death.naive.CIL",varymax="cinc.B.death.naive.CIU",vcol="#009E73", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.C.death.naive.CIL",varymax="cinc.C.death.naive.CIU",vcol="#CC79A7", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.D.death.naive.CIL",varymax="cinc.D.death.naive.CIU",vcol="#56B4E9", valpha=0.2)+
  #axes
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,0.10,0.05),limits=c(0,0.10))+
  scale_colour_manual(NULL,values=cols,labels=c(A="Treatment A",B="Treatment B",C="Treatment C",D="Treatment D"), breaks=c("A","B","C","D"))+
  scale_linetype_manual(NULL,values=ltys,labels=c("True","naive"), breaks=c("solid","dashed"))+
  ylab("Cumulative incidence")+xlab("Time (years)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))+
  theme(legend.position="inside",legend.position.inside =c(0.2,0.7))+
  ggtitle("Estimation method: Naive")

plot.death.lmtp=ggplot(plot.dta.death,aes(x=time.grid,y=cinc.A.death.true))+theme_bw()+
  #true
  addlinetoplot(plot.dta.death,"time.grid","cinc.A.death.true",vcol='"A"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.B.death.true",vcol='"B"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.C.death.true",vcol='"C"',vlty='"true"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.D.death.true",vcol='"D"',vlty='"true"')+
  #lmtp - estimates
  addlinetoplot(plot.dta.death,"time.grid","cinc.A.death.lmtp",vcol='"A"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.B.death.lmtp",vcol='"B"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.C.death.lmtp",vcol='"C"',vlty='"other"')+
  addlinetoplot(plot.dta.death,"time.grid","cinc.D.death.lmtp",vcol='"D"',vlty='"other"')+
  #lmtp - CIs
  addribbontoplot(plot.dta.death,varymin="cinc.A.death.lmtp.CIL",varymax="cinc.A.death.lmtp.CIU",vcol="#D55E00", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.B.death.lmtp.CIL",varymax="cinc.B.death.lmtp.CIU",vcol="#009E73", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.C.death.lmtp.CIL",varymax="cinc.C.death.lmtp.CIU",vcol="#CC79A7", valpha=0.2)+
  addribbontoplot(plot.dta.death,varymin="cinc.D.death.lmtp.CIL",varymax="cinc.D.death.lmtp.CIU",vcol="#56B4E9", valpha=0.2)+
  #axes
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,0.10,0.05),limits=c(0,0.10))+
  scale_colour_manual(NULL,values=cols,labels=c(A="Treatment A",B="Treatment B",C="Treatment C",D="Treatment D"), breaks=c("A","B","C","D"))+
  scale_linetype_manual(NULL,values=ltys,labels=c("True","lmtp"), breaks=c("solid","dashed"))+
  ylab("Cumulative incidence")+xlab("Time (years)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))+
  theme(legend.position="inside",legend.position.inside =c(0.2,0.7))+
  ggtitle("Estimation method: Doubly robust")


windows(10,7)
grid.arrange(plot.death.naive,plot.death.iptw,plot.death.gform,plot.death.lmtp,ncol=2)

#------------------------------------------------------
#Results tables - for 5 time points
#------------------------------------------------------

table.dta.mace<-plot.dta.mace[plot.dta.mace$time.grid%in%c(1,2,3,4,5),]
table.dta.death<-plot.dta.death[plot.dta.death$time.grid%in%c(1,2,3,4,5),]

#functions for printing results
tabf<-function(x){sprintf("%.3f",x)}
tabf2<-function(est,CIL,CIU){paste0(tabf(est)," (",tabf(CIL),",",tabf(CIU),")")}

#---
#mace

tab.mace.A<-rbind(c(tabf2(table.dta.mace$cinc.A.mace.true,table.dta.mace$cinc.A.mace.true.CIL,table.dta.mace$cinc.A.mace.true.CIU)),
      c(tabf2(table.dta.mace$cinc.A.mace.naive,table.dta.mace$cinc.A.mace.naive.CIL,table.dta.mace$cinc.A.mace.naive.CIU)),
      c(tabf2(table.dta.mace$cinc.A.mace.iptw,table.dta.mace$cinc.A.mace.iptw.CIL,table.dta.mace$cinc.A.mace.iptw.CIU)),
      c(tabf2(table.dta.mace$cinc.A.mace.gform,table.dta.mace$cinc.A.mace.gform.CIL,table.dta.mace$cinc.A.mace.gform.CIU)),
      c(tabf2(table.dta.mace$cinc.A.mace.lmtp,table.dta.mace$cinc.A.mace.lmtp.CIL,table.dta.mace$cinc.A.mace.lmtp.CIU)))
rownames(tab.mace.A)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.mace.A)<-c("Time 1","Time 2","Time 3","Time 4","Time 5")

tab.mace.B<-rbind(c(tabf2(table.dta.mace$cinc.B.mace.true,table.dta.mace$cinc.B.mace.true.CIL,table.dta.mace$cinc.B.mace.true.CIU)),
                  c(tabf2(table.dta.mace$cinc.B.mace.naive,table.dta.mace$cinc.B.mace.naive.CIL,table.dta.mace$cinc.B.mace.naive.CIU)),
                  c(tabf2(table.dta.mace$cinc.B.mace.iptw,table.dta.mace$cinc.B.mace.iptw.CIL,table.dta.mace$cinc.B.mace.iptw.CIU)),
                  c(tabf2(table.dta.mace$cinc.B.mace.gform,table.dta.mace$cinc.B.mace.gform.CIL,table.dta.mace$cinc.B.mace.gform.CIU)),
                  c(tabf2(table.dta.mace$cinc.B.mace.lmtp,table.dta.mace$cinc.B.mace.lmtp.CIL,table.dta.mace$cinc.B.mace.lmtp.CIU)))
rownames(tab.mace.B)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.mace.B)<-c("Time 1","Time 2","Time 3","Time 4","Time 5")

tab.mace.C<-rbind(c(tabf2(table.dta.mace$cinc.C.mace.true,table.dta.mace$cinc.C.mace.true.CIL,table.dta.mace$cinc.C.mace.true.CIU)),
                  c(tabf2(table.dta.mace$cinc.C.mace.naive,table.dta.mace$cinc.C.mace.naive.CIL,table.dta.mace$cinc.C.mace.naive.CIU)),
                  c(tabf2(table.dta.mace$cinc.C.mace.iptw,table.dta.mace$cinc.C.mace.iptw.CIL,table.dta.mace$cinc.C.mace.iptw.CIU)),
                  c(tabf2(table.dta.mace$cinc.C.mace.gform,table.dta.mace$cinc.C.mace.gform.CIL,table.dta.mace$cinc.C.mace.gform.CIU)),
                  c(tabf2(table.dta.mace$cinc.C.mace.lmtp,table.dta.mace$cinc.C.mace.lmtp.CIL,table.dta.mace$cinc.C.mace.lmtp.CIU)))
rownames(tab.mace.C)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.mace.C)<-c("Time 1","Time 2","Time 3","Time 4","Time 5")

tab.mace.D<-rbind(c(tabf2(table.dta.mace$cinc.D.mace.true,table.dta.mace$cinc.D.mace.true.CIL,table.dta.mace$cinc.D.mace.true.CIU)),
                  c(tabf2(table.dta.mace$cinc.D.mace.naive,table.dta.mace$cinc.D.mace.naive.CIL,table.dta.mace$cinc.D.mace.naive.CIU)),
                  c(tabf2(table.dta.mace$cinc.D.mace.iptw,table.dta.mace$cinc.D.mace.iptw.CIL,table.dta.mace$cinc.D.mace.iptw.CIU)),
                  c(tabf2(table.dta.mace$cinc.D.mace.gform,table.dta.mace$cinc.D.mace.gform.CIL,table.dta.mace$cinc.D.mace.gform.CIU)),
                  c(tabf2(table.dta.mace$cinc.D.mace.lmtp,table.dta.mace$cinc.D.mace.lmtp.CIL,table.dta.mace$cinc.D.mace.lmtp.CIU)))
rownames(tab.mace.D)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.mace.D)<-c("Time 1","Time 2","Time 3","Time 4","Time 5")

#---
#death

tab.death.A<-rbind(c(tabf2(table.dta.death$cinc.A.death.true,table.dta.death$cinc.A.death.true.CIL,table.dta.death$cinc.A.death.true.CIU)),
                  c(tabf2(table.dta.death$cinc.A.death.naive,table.dta.death$cinc.A.death.naive.CIL,table.dta.death$cinc.A.death.naive.CIU)),
                  c(tabf2(table.dta.death$cinc.A.death.iptw,table.dta.death$cinc.A.death.iptw.CIL,table.dta.death$cinc.A.death.iptw.CIU)),
                  c(tabf2(table.dta.death$cinc.A.death.gform,table.dta.death$cinc.A.death.gform.CIL,table.dta.death$cinc.A.death.gform.CIU)),
                  c(tabf2(table.dta.death$cinc.A.death.lmtp,table.dta.death$cinc.A.death.lmtp.CIL,table.dta.death$cinc.A.death.lmtp.CIU)))
rownames(tab.death.A)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.death.A)<-c("Time 1","Time 2","Time 3","Time 4","Time 5")

tab.death.B<-rbind(c(tabf2(table.dta.death$cinc.B.death.true,table.dta.death$cinc.B.death.true.CIL,table.dta.death$cinc.B.death.true.CIU)),
                  c(tabf2(table.dta.death$cinc.B.death.naive,table.dta.death$cinc.B.death.naive.CIL,table.dta.death$cinc.B.death.naive.CIU)),
                  c(tabf2(table.dta.death$cinc.B.death.iptw,table.dta.death$cinc.B.death.iptw.CIL,table.dta.death$cinc.B.death.iptw.CIU)),
                  c(tabf2(table.dta.death$cinc.B.death.gform,table.dta.death$cinc.B.death.gform.CIL,table.dta.death$cinc.B.death.gform.CIU)),
                  c(tabf2(table.dta.death$cinc.B.death.lmtp,table.dta.death$cinc.B.death.lmtp.CIL,table.dta.death$cinc.B.death.lmtp.CIU)))
rownames(tab.death.B)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.death.B)<-c("Time 1","Time 2","Time 3","Time 4","Time 5")

tab.death.C<-rbind(c(tabf2(table.dta.death$cinc.C.death.true,table.dta.death$cinc.C.death.true.CIL,table.dta.death$cinc.C.death.true.CIU)),
                  c(tabf2(table.dta.death$cinc.C.death.naive,table.dta.death$cinc.C.death.naive.CIL,table.dta.death$cinc.C.death.naive.CIU)),
                  c(tabf2(table.dta.death$cinc.C.death.iptw,table.dta.death$cinc.C.death.iptw.CIL,table.dta.death$cinc.C.death.iptw.CIU)),
                  c(tabf2(table.dta.death$cinc.C.death.gform,table.dta.death$cinc.C.death.gform.CIL,table.dta.death$cinc.C.death.gform.CIU)),
                  c(tabf2(table.dta.death$cinc.C.death.lmtp,table.dta.death$cinc.C.death.lmtp.CIL,table.dta.death$cinc.C.death.lmtp.CIU)))
rownames(tab.death.C)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.death.C)<-c("Time 1","Time 2","Time 3","Time 4","Time 5")

tab.death.D<-rbind(c(tabf2(table.dta.death$cinc.D.death.true,table.dta.death$cinc.D.death.true.CIL,table.dta.death$cinc.D.death.true.CIU)),
                  c(tabf2(table.dta.death$cinc.D.death.naive,table.dta.death$cinc.D.death.naive.CIL,table.dta.death$cinc.D.death.naive.CIU)),
                  c(tabf2(table.dta.death$cinc.D.death.iptw,table.dta.death$cinc.D.death.iptw.CIL,table.dta.death$cinc.D.death.iptw.CIU)),
                  c(tabf2(table.dta.death$cinc.D.death.gform,table.dta.death$cinc.D.death.gform.CIL,table.dta.death$cinc.D.death.gform.CIU)),
                  c(tabf2(table.dta.death$cinc.D.death.lmtp,table.dta.death$cinc.D.death.lmtp.CIL,table.dta.death$cinc.D.death.lmtp.CIU)))
rownames(tab.death.D)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.death.D)<-c("Time 1","Time 2","Time 3","Time 4","Time 5")

#------------------------------------------------------
#Results tables - for a single specified time point (here time 5)
#------------------------------------------------------

table.dta.mace5<-plot.dta.mace[plot.dta.mace$time.grid==5,]
table.dta.death5<-plot.dta.death[plot.dta.death$time.grid==5,]

#functions for printing results
tabf<-function(x){sprintf("%.3f",x)}
tabf2<-function(est,CIL,CIU){paste0(tabf(est)," (",tabf(CIL),",",tabf(CIU),")")}

#---
#mace

tab.mace5.A<-rbind(c(tabf2(table.dta.mace5$cinc.A.mace.true,table.dta.mace5$cinc.A.mace.true.CIL,table.dta.mace5$cinc.A.mace.true.CIU)),
                  c(tabf2(table.dta.mace5$cinc.A.mace.naive,table.dta.mace5$cinc.A.mace.naive.CIL,table.dta.mace5$cinc.A.mace.naive.CIU)),
                  c(tabf2(table.dta.mace5$cinc.A.mace.iptw,table.dta.mace5$cinc.A.mace.iptw.CIL,table.dta.mace5$cinc.A.mace.iptw.CIU)),
                  c(tabf2(table.dta.mace5$cinc.A.mace.gform,table.dta.mace5$cinc.A.mace.gform.CIL,table.dta.mace5$cinc.A.mace.gform.CIU)),
                  c(tabf2(table.dta.mace5$cinc.A.mace.lmtp,table.dta.mace5$cinc.A.mace.lmtp.CIL,table.dta.mace5$cinc.A.mace.lmtp.CIU)))

tab.mace5.B<-rbind(c(tabf2(table.dta.mace5$cinc.B.mace.true,table.dta.mace5$cinc.B.mace.true.CIL,table.dta.mace5$cinc.B.mace.true.CIU)),
                   c(tabf2(table.dta.mace5$cinc.B.mace.naive,table.dta.mace5$cinc.B.mace.naive.CIL,table.dta.mace5$cinc.B.mace.naive.CIU)),
                   c(tabf2(table.dta.mace5$cinc.B.mace.iptw,table.dta.mace5$cinc.B.mace.iptw.CIL,table.dta.mace5$cinc.B.mace.iptw.CIU)),
                   c(tabf2(table.dta.mace5$cinc.B.mace.gform,table.dta.mace5$cinc.B.mace.gform.CIL,table.dta.mace5$cinc.B.mace.gform.CIU)),
                   c(tabf2(table.dta.mace5$cinc.B.mace.lmtp,table.dta.mace5$cinc.B.mace.lmtp.CIL,table.dta.mace5$cinc.B.mace.lmtp.CIU)))

tab.mace5.C<-rbind(c(tabf2(table.dta.mace5$cinc.C.mace.true,table.dta.mace5$cinc.C.mace.true.CIL,table.dta.mace5$cinc.C.mace.true.CIU)),
                   c(tabf2(table.dta.mace5$cinc.C.mace.naive,table.dta.mace5$cinc.C.mace.naive.CIL,table.dta.mace5$cinc.C.mace.naive.CIU)),
                   c(tabf2(table.dta.mace5$cinc.C.mace.iptw,table.dta.mace5$cinc.C.mace.iptw.CIL,table.dta.mace5$cinc.C.mace.iptw.CIU)),
                   c(tabf2(table.dta.mace5$cinc.C.mace.gform,table.dta.mace5$cinc.C.mace.gform.CIL,table.dta.mace5$cinc.C.mace.gform.CIU)),
                   c(tabf2(table.dta.mace5$cinc.C.mace.lmtp,table.dta.mace5$cinc.C.mace.lmtp.CIL,table.dta.mace5$cinc.C.mace.lmtp.CIU)))

tab.mace5.D<-rbind(c(tabf2(table.dta.mace5$cinc.D.mace.true,table.dta.mace5$cinc.D.mace.true.CIL,table.dta.mace5$cinc.D.mace.true.CIU)),
                   c(tabf2(table.dta.mace5$cinc.D.mace.naive,table.dta.mace5$cinc.D.mace.naive.CIL,table.dta.mace5$cinc.D.mace.naive.CIU)),
                   c(tabf2(table.dta.mace5$cinc.D.mace.iptw,table.dta.mace5$cinc.D.mace.iptw.CIL,table.dta.mace5$cinc.D.mace.iptw.CIU)),
                   c(tabf2(table.dta.mace5$cinc.D.mace.gform,table.dta.mace5$cinc.D.mace.gform.CIL,table.dta.mace5$cinc.D.mace.gform.CIU)),
                   c(tabf2(table.dta.mace5$cinc.D.mace.lmtp,table.dta.mace5$cinc.D.mace.lmtp.CIL,table.dta.mace5$cinc.D.mace.lmtp.CIU)))

tab.mace5<-cbind(tab.mace5.A,tab.mace5.B,tab.mace5.C,tab.mace5.D)

rownames(tab.mace5)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.mace5)<-c("A","B","C","D")

#---
#death

tab.death5.A<-rbind(c(tabf2(table.dta.death5$cinc.A.death.true,table.dta.death5$cinc.A.death.true.CIL,table.dta.death5$cinc.A.death.true.CIU)),
                   c(tabf2(table.dta.death5$cinc.A.death.naive,table.dta.death5$cinc.A.death.naive.CIL,table.dta.death5$cinc.A.death.naive.CIU)),
                   c(tabf2(table.dta.death5$cinc.A.death.iptw,table.dta.death5$cinc.A.death.iptw.CIL,table.dta.death5$cinc.A.death.iptw.CIU)),
                   c(tabf2(table.dta.death5$cinc.A.death.gform,table.dta.death5$cinc.A.death.gform.CIL,table.dta.death5$cinc.A.death.gform.CIU)),
                   c(tabf2(table.dta.death5$cinc.A.death.lmtp,table.dta.death5$cinc.A.death.lmtp.CIL,table.dta.death5$cinc.A.death.lmtp.CIU)))

tab.death5.B<-rbind(c(tabf2(table.dta.death5$cinc.B.death.true,table.dta.death5$cinc.B.death.true.CIL,table.dta.death5$cinc.B.death.true.CIU)),
                   c(tabf2(table.dta.death5$cinc.B.death.naive,table.dta.death5$cinc.B.death.naive.CIL,table.dta.death5$cinc.B.death.naive.CIU)),
                   c(tabf2(table.dta.death5$cinc.B.death.iptw,table.dta.death5$cinc.B.death.iptw.CIL,table.dta.death5$cinc.B.death.iptw.CIU)),
                   c(tabf2(table.dta.death5$cinc.B.death.gform,table.dta.death5$cinc.B.death.gform.CIL,table.dta.death5$cinc.B.death.gform.CIU)),
                   c(tabf2(table.dta.death5$cinc.B.death.lmtp,table.dta.death5$cinc.B.death.lmtp.CIL,table.dta.death5$cinc.B.death.lmtp.CIU)))

tab.death5.C<-rbind(c(tabf2(table.dta.death5$cinc.C.death.true,table.dta.death5$cinc.C.death.true.CIL,table.dta.death5$cinc.C.death.true.CIU)),
                   c(tabf2(table.dta.death5$cinc.C.death.naive,table.dta.death5$cinc.C.death.naive.CIL,table.dta.death5$cinc.C.death.naive.CIU)),
                   c(tabf2(table.dta.death5$cinc.C.death.iptw,table.dta.death5$cinc.C.death.iptw.CIL,table.dta.death5$cinc.C.death.iptw.CIU)),
                   c(tabf2(table.dta.death5$cinc.C.death.gform,table.dta.death5$cinc.C.death.gform.CIL,table.dta.death5$cinc.C.death.gform.CIU)),
                   c(tabf2(table.dta.death5$cinc.C.death.lmtp,table.dta.death5$cinc.C.death.lmtp.CIL,table.dta.death5$cinc.C.death.lmtp.CIU)))

tab.death5.D<-rbind(c(tabf2(table.dta.death5$cinc.D.death.true,table.dta.death5$cinc.D.death.true.CIL,table.dta.death5$cinc.D.death.true.CIU)),
                   c(tabf2(table.dta.death5$cinc.D.death.naive,table.dta.death5$cinc.D.death.naive.CIL,table.dta.death5$cinc.D.death.naive.CIU)),
                   c(tabf2(table.dta.death5$cinc.D.death.iptw,table.dta.death5$cinc.D.death.iptw.CIL,table.dta.death5$cinc.D.death.iptw.CIU)),
                   c(tabf2(table.dta.death5$cinc.D.death.gform,table.dta.death5$cinc.D.death.gform.CIL,table.dta.death5$cinc.D.death.gform.CIU)),
                   c(tabf2(table.dta.death5$cinc.D.death.lmtp,table.dta.death5$cinc.D.death.lmtp.CIL,table.dta.death5$cinc.D.death.lmtp.CIU)))

tab.death5<-cbind(tab.death5.A,tab.death5.B,tab.death5.C,tab.death5.D)

rownames(tab.death5)<-c("True","Naive","IPTW","G-formula","DR")
colnames(tab.death5)<-c("A","B","C","D")
