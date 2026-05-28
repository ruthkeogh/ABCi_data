#######################################################
#Q2: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D 
#or never using treatment, in people with hb>=7.5. 
#
#This file sets up the data in the format needed for analysis
#######################################################

#------------------------------------------------------
#DATA SET-UP: for treatment strategies A,B,C,D
#We identify people starting treatment and set the time of treatment start to be their time zero
#------------------------------------------------------

dta.trt<-dta

#identify people meeting eligibility criteria: people who have started a second line treatment
dta.trt<-dta
dta.trt$treat_status_any<-ifelse((dta.trt$treat_status_A+dta.trt$treat_status_B+
                                    dta.trt$treat_status_C+dta.trt$treat_status_D)>=1,1,0)
dta.trt<-dta.trt[dta.trt$treat_status_any==1,]

#generate rownumber
dta.trt<-dta.trt%>%group_by(id)%>%mutate(rownum=row_number())

#Restrict to the first row
dta.trt<-dta.trt[dta.trt$rownum==1,]

#restrict to people with HbA1c>=7.5 on starting treatment
dta.trt<-dta.trt[dta.trt$hb>=7.5,]

#generate new event time that is relative to the time of starting treatment
dta.trt$tstart<-0
dta.trt$tstop<-dta.trt$event_time-dta.trt$tstart

#generate categorical treatment type variable
dta.trt$treat_type<-""
dta.trt$treat_type<-ifelse(dta.trt$treat_status_A==1,"A",dta.trt$treat_type)
dta.trt$treat_type<-ifelse(dta.trt$treat_status_B==1,"B",dta.trt$treat_type)
dta.trt$treat_type<-ifelse(dta.trt$treat_status_C==1,"C",dta.trt$treat_type)
dta.trt$treat_type<-ifelse(dta.trt$treat_status_D==1,"D",dta.trt$treat_type)

#generate factor variable for competing events
dta.trt$event_type_cr<-ifelse(dta.trt$event_type==3,0,dta.trt$event_type)
dta.trt$event_type_cr<-as.factor(dta.trt$event_type_cr)

#------------------------------------------------------
#DATA SET-UP: for treatment strategy of never starting treatment
#We identify people untreated in the first time interval (and meeting other eligibility criteria)
#People are then censored at the time of treatment initiation, if they start treatment. 
#------------------------------------------------------

dta.nevertrt<-dta
dta.nevertrt$treat_type<-"nevertrt"

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

#generate factor variable for competing events
dta.nevertrt$event_type_cr=0
dta.nevertrt$event_type_cr<-ifelse(dta.nevertrt$mace_status==1,1,dta.nevertrt$event_type_cr)
dta.nevertrt$event_type_cr<-ifelse(dta.nevertrt$death_status==1,2,dta.nevertrt$event_type_cr)
dta.nevertrt$event_type_cr<-as.factor(dta.nevertrt$event_type_cr)

#------------------------------------------------------
#Create combined data set
#Noting that dta.trt has 1 row per person and dta.nevertrt has multiple rows per person
#Note it is convenient to do this here because this is the data that gets bootstrapped.
#------------------------------------------------------


dta.analysis<-rbind(dta.trt,dta.nevertrt)

#generate new id number, e.g. for use in the bootstrapping
#because some people feature twice in dta.analysis (e.g. contributing as a treated person and untreated person pre-treatment)
dta.analysis$id_new<-paste0(dta.analysis$id,".",dta.analysis$treat_type)
