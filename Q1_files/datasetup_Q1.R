#######################################################
#Q1: Estimates cumulative incidences of mace and death 
#under the strategies of starting treatment A,B,C or D in people with hb>=7.5
#
#This file sets up the data in the format needed for analysis
#Creates the data set dta.analysis
#######################################################

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
