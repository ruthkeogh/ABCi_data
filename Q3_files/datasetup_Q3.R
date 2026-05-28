#######################################################
#Q3: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D with a grace period of 90 days
#in people with hb>=7.5. 
#
#This file sets up the data in the format needed for analysis using the clone-censor-weight approach
#######################################################

#------------------------------------------------------
#DATA SET-UP
#------------------------------------------------------

dta.analysis<-dta

#identify people with hb>=7.5 at time 0
dta.analysis<-dta.analysis%>%group_by(id)%>%mutate(hb_0=first(hb))
dta.analysis<-dta.analysis[dta.analysis$hb_0>=7.5,]

#Grace period: gp can be changed here but it should be on the grid of times in the data for the code to work without modification
gp<-90
gp.plus<-gp+30 #end of the time period following the end of the grace period

#create factor variable for competing events
dta.analysis$event_type_cr=0
dta.analysis$event_type_cr<-ifelse(dta.analysis$mace_status==1,1,dta.analysis$event_type_cr)
dta.analysis$event_type_cr<-ifelse(dta.analysis$death_status==1,2,dta.analysis$event_type_cr)
dta.analysis$event_type_cr<-as.factor(dta.analysis$event_type_cr)

#------------------------------------------------------
#data set for strategy: grace period of 3 periods (90 days) and then must start treatment A
#------------------------------------------------------

dta.trtA<-dta.analysis
dta.trtA<-filter(dta.trtA,
                 (!is.na(treat_time_A) & treat_time_A<gp.plus)|
                   (is.na(treat_time_A) & is.na(treat_time_B) & is.na(treat_time_C) & is.na(treat_time_D) & tstop<=gp)|
                   (!is.na(treat_time_B) & tstart<=treat_time_B & tstop<=gp)|
                   (!is.na(treat_time_C) & tstart<=treat_time_C & tstop<=gp)|
                   (!is.na(treat_time_D) & tstart<=treat_time_D & tstop<=gp))

#change tstop to the time of starting other treatments, if that occurs within the grace period
dta.trtA<-dta.trtA%>%mutate(tstop.new=tstop)
dta.trtA<-dta.trtA%>%mutate(tstop.new=ifelse(!is.na(treat_time_B) & tstop.new>treat_time_B,treat_time_B,tstop.new))
dta.trtA<-dta.trtA%>%mutate(tstop.new=ifelse(!is.na(treat_time_C) & tstop.new>treat_time_C,treat_time_C,tstop.new))
dta.trtA<-dta.trtA%>%mutate(tstop.new=ifelse(!is.na(treat_time_D) & tstop.new>treat_time_D,treat_time_D,tstop.new))

#generate an indicator of censoring due to deviation from the treatment strategy
dta.trtA<-dta.trtA%>%mutate(cens.gp=0)
dta.trtA<-dta.trtA%>%mutate(cens.gp=ifelse(!is.na(treat_time_B) & tstop.new==treat_time_B,1,cens.gp))
dta.trtA<-dta.trtA%>%mutate(cens.gp=ifelse(!is.na(treat_time_B) & treat_time_B>gp & tstop.new==gp,1,cens.gp))
dta.trtA<-dta.trtA%>%mutate(cens.gp=ifelse(!is.na(treat_time_C) & tstop.new==treat_time_C,1,cens.gp))
dta.trtA<-dta.trtA%>%mutate(cens.gp=ifelse(!is.na(treat_time_C) & treat_time_C>gp & tstop.new==gp,1,cens.gp))
dta.trtA<-dta.trtA%>%mutate(cens.gp=ifelse(!is.na(treat_time_D) & tstop.new==treat_time_D,1,cens.gp))
dta.trtA<-dta.trtA%>%mutate(cens.gp=ifelse(!is.na(treat_time_D) & treat_time_D>gp & tstop.new==gp,1,cens.gp))
dta.trtA<-dta.trtA%>%mutate(cens.gp=ifelse(is.na(treat_time_A) & is.na(treat_time_B) & is.na(treat_time_C) & is.na(treat_time_D) &
                                             tstop.new==gp,1,cens.gp))

dta.trtA$trt.strategy<-"A"

#------------------------------------------------------
#data set for strategy: grace period of 3 periods (90 days) and then must start treatment B
#This is identical to what we did for treatment A above, but now for treatment B.
#------------------------------------------------------

dta.trtB<-dta.analysis
dta.trtB<-filter(dta.trtB,
                 (!is.na(treat_time_B) & treat_time_B<gp.plus)|
                   (is.na(treat_time_A) & is.na(treat_time_B) & is.na(treat_time_C) & is.na(treat_time_D) & tstop<=gp)|
                   (!is.na(treat_time_A) & tstart<=treat_time_A & tstop<=gp)|
                   (!is.na(treat_time_C) & tstart<=treat_time_C & tstop<=gp)|
                   (!is.na(treat_time_D) & tstart<=treat_time_D & tstop<=gp))

#change tstop to the time of starting other treatments, if that occurs within the grace period
dta.trtB<-dta.trtB%>%mutate(tstop.new=tstop)
dta.trtB<-dta.trtB%>%mutate(tstop.new=ifelse(!is.na(treat_time_A) & tstop.new>treat_time_A,treat_time_A,tstop.new))
dta.trtB<-dta.trtB%>%mutate(tstop.new=ifelse(!is.na(treat_time_C) & tstop.new>treat_time_C,treat_time_C,tstop.new))
dta.trtB<-dta.trtB%>%mutate(tstop.new=ifelse(!is.na(treat_time_D) & tstop.new>treat_time_D,treat_time_D,tstop.new))

#and generate an indicator of censoring due to deviation from the treatment strategy
dta.trtB<-dta.trtB%>%mutate(cens.gp=0)
dta.trtB<-dta.trtB%>%mutate(cens.gp=ifelse(!is.na(treat_time_A) & tstop.new==treat_time_A,1,cens.gp))
dta.trtB<-dta.trtB%>%mutate(cens.gp=ifelse(!is.na(treat_time_A) & treat_time_A>gp & tstop.new==gp,1,cens.gp))
dta.trtB<-dta.trtB%>%mutate(cens.gp=ifelse(!is.na(treat_time_C) & tstop.new==treat_time_C,1,cens.gp))
dta.trtB<-dta.trtB%>%mutate(cens.gp=ifelse(!is.na(treat_time_C) & treat_time_C>gp & tstop.new==gp,1,cens.gp))
dta.trtB<-dta.trtB%>%mutate(cens.gp=ifelse(!is.na(treat_time_D) & tstop.new==treat_time_D,1,cens.gp))
dta.trtB<-dta.trtB%>%mutate(cens.gp=ifelse(!is.na(treat_time_D) & treat_time_D>gp & tstop.new==gp,1,cens.gp))
dta.trtB<-dta.trtB%>%mutate(cens.gp=ifelse(is.na(treat_time_A) & is.na(treat_time_B) & is.na(treat_time_C) & is.na(treat_time_D) &
                                             tstop.new==gp,1,cens.gp))

dta.trtB$trt.strategy<-"B"

#------------------------------------------------------
#data set for strategy: grace period of 3 periods (90 days) and then must start treatment C
#This is identical to what we did for treatment A/B above, but now for treatment C.
#------------------------------------------------------

dta.trtC<-dta.analysis
dta.trtC<-filter(dta.trtC,
                 (!is.na(treat_time_C) & treat_time_C<gp.plus)|
                   (is.na(treat_time_A) & is.na(treat_time_B) & is.na(treat_time_C) & is.na(treat_time_D) & tstop<=gp)|
                   (!is.na(treat_time_A) & tstart<=treat_time_A & tstop<=gp)|
                   (!is.na(treat_time_B) & tstart<=treat_time_B & tstop<=gp)|
                   (!is.na(treat_time_D) & tstart<=treat_time_D & tstop<=gp))

#change tstop to the time of starting other treatments, if that occurs within the grace period
dta.trtC<-dta.trtC%>%mutate(tstop.new=tstop)
dta.trtC<-dta.trtC%>%mutate(tstop.new=ifelse(!is.na(treat_time_A) & tstop.new>treat_time_A,treat_time_A,tstop.new))
dta.trtC<-dta.trtC%>%mutate(tstop.new=ifelse(!is.na(treat_time_B) & tstop.new>treat_time_B,treat_time_B,tstop.new))
dta.trtC<-dta.trtC%>%mutate(tstop.new=ifelse(!is.na(treat_time_D) & tstop.new>treat_time_D,treat_time_D,tstop.new))

#and generate an indicator of censoring due to deviation from the treatment strategy
dta.trtC<-dta.trtC%>%mutate(cens.gp=0)
dta.trtC<-dta.trtC%>%mutate(cens.gp=ifelse(!is.na(treat_time_A) & tstop.new==treat_time_A,1,cens.gp))
dta.trtC<-dta.trtC%>%mutate(cens.gp=ifelse(!is.na(treat_time_A) & treat_time_A>gp & tstop.new==gp,1,cens.gp))
dta.trtC<-dta.trtC%>%mutate(cens.gp=ifelse(!is.na(treat_time_B) & tstop.new==treat_time_B,1,cens.gp))
dta.trtC<-dta.trtC%>%mutate(cens.gp=ifelse(!is.na(treat_time_B) & treat_time_B>gp & tstop.new==gp,1,cens.gp))
dta.trtC<-dta.trtC%>%mutate(cens.gp=ifelse(!is.na(treat_time_D) & tstop.new==treat_time_D,1,cens.gp))
dta.trtC<-dta.trtC%>%mutate(cens.gp=ifelse(!is.na(treat_time_D) & treat_time_D>gp & tstop.new==gp,1,cens.gp))
dta.trtC<-dta.trtC%>%mutate(cens.gp=ifelse(is.na(treat_time_A) & is.na(treat_time_B) & is.na(treat_time_C) & is.na(treat_time_D) &
                                             tstop.new==gp,1,cens.gp))

dta.trtC$trt.strategy<-"C"

#------------------------------------------------------
#data set for strategy: grace period of 3 periods (90 days) and then must start treatment D
#This is identical to what we did for treatment A/B/C above, but now for treatment D.
#------------------------------------------------------

dta.trtD<-dta.analysis
dta.trtD<-filter(dta.trtD,
                 (!is.na(treat_time_D) & treat_time_D<gp.plus)|
                   (is.na(treat_time_A) & is.na(treat_time_B) & is.na(treat_time_C) & is.na(treat_time_D) & tstop<=gp)|
                   (!is.na(treat_time_A) & tstart<=treat_time_A & tstop<=gp)|
                   (!is.na(treat_time_B) & tstart<=treat_time_B & tstop<=gp)|
                   (!is.na(treat_time_C) & tstart<=treat_time_C & tstop<=gp))

#change tstop to the time of starting other treatments, if that occurs within the grace period
dta.trtD<-dta.trtD%>%mutate(tstop.new=tstop)
dta.trtD<-dta.trtD%>%mutate(tstop.new=ifelse(!is.na(treat_time_A) & tstop.new>treat_time_A,treat_time_A,tstop.new))
dta.trtD<-dta.trtD%>%mutate(tstop.new=ifelse(!is.na(treat_time_B) & tstop.new>treat_time_B,treat_time_B,tstop.new))
dta.trtD<-dta.trtD%>%mutate(tstop.new=ifelse(!is.na(treat_time_C) & tstop.new>treat_time_C,treat_time_C,tstop.new))

#and generate an indicator of censoring due to deviation from the treatment strategy
dta.trtD<-dta.trtD%>%mutate(cens.gp=0)
dta.trtD<-dta.trtD%>%mutate(cens.gp=ifelse(!is.na(treat_time_A) & tstop.new==treat_time_A,1,cens.gp))
dta.trtD<-dta.trtD%>%mutate(cens.gp=ifelse(!is.na(treat_time_A) & treat_time_A>gp & tstop.new==gp,1,cens.gp))
dta.trtD<-dta.trtD%>%mutate(cens.gp=ifelse(!is.na(treat_time_B) & tstop.new==treat_time_B,1,cens.gp))
dta.trtD<-dta.trtD%>%mutate(cens.gp=ifelse(!is.na(treat_time_B) & treat_time_B>gp & tstop.new==gp,1,cens.gp))
dta.trtD<-dta.trtD%>%mutate(cens.gp=ifelse(!is.na(treat_time_C) & tstop.new==treat_time_C,1,cens.gp))
dta.trtD<-dta.trtD%>%mutate(cens.gp=ifelse(!is.na(treat_time_C) & treat_time_C>gp & tstop.new==gp,1,cens.gp))
dta.trtD<-dta.trtD%>%mutate(cens.gp=ifelse(is.na(treat_time_A) & is.na(treat_time_B) & is.na(treat_time_C) & is.na(treat_time_D) &
                                             tstop.new==gp,1,cens.gp))

dta.trtD$trt.strategy<-"D"

#------------------------------------------------------
#data set for strategy of never starting treatment
#We identify people untreated in the first time interval (and meeting other eligibility criteria)
#People are then censored at the time of treatment initiation, if they start treatment. 
#This is the same as in Q2
#------------------------------------------------------

dta.nevertrt<-dta
dta.nevertrt$trt.strategy<-"nevertrt"

#restrict to people with hb>=7.5 at time 0
dta.nevertrt<-dta.nevertrt%>%group_by(id)%>%mutate(hb_first=first(hb))
dta.nevertrt<-dta.nevertrt[dta.nevertrt$hb_first>=7.5,]
dta.nevertrt$hb_first<-NULL

#generate lagged treatment status
dta.nevertrt<-dta.nevertrt%>%group_by(id)%>%mutate(lag_any_treat=lag(treat_status_A,default=0)+
                                                     lag(treat_status_B,default=0)+lag(treat_status_C,default=0)+
                                                     +lag(treat_status_D,default=0))

#restrict to rows up to and including treatment initiation
#This is needed for estimation of the weights. Later we will drop the row in which the person initiates treatment i.e. the last row)
dta.nevertrt<-dta.nevertrt[dta.nevertrt$lag_any_treat==0,]
dta.nevertrt$lag_any_treat<-NULL