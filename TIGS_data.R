library(tidyverse)
library(lubridate)
library(ggplot2)
library(corrplot)
library(reshape2)


TIGS<-read.csv("C:/Users/LENOVO/Downloads/TIGS_data.csv")
TIGS<-TIGS %>% select(c("Sample.ID","Date","Final.Result","Viral.Load..copies.ml.","STP"))


# TIGS data cleaning

class(TIGS$Date)
class(TIGS$date_onset)
class(TIGS$Viral.Load..copies.ml.)

#Homogenizing date formats

TIGS$Date<-parse_date_time(x = TIGS$Date,
                orders = c("d-m-y", "d.m.y", "m/d/y"),
                locale = "eng")

TIGS <- TIGS %>% mutate(date_onset=as.Date(TIGS$Date, format = "%y/%m/%d"))

#Removing inconclusive and invalid entries

TIGS1<-TIGS[!(TIGS$Final.Result=="Inconclusive"),]
TIGS1<-TIGS1[!(TIGS1$Final.Result=="Invalid"),]

#marking entries with Negative result and NA  viral loads as zero

TIGS1$Viral.Load..copies.ml.[TIGS1$Final.Result=="Negative" & is.na(TIGS1$Viral.Load..copies.ml.)]<-0
TIGS$Viral.Load..copies.ml.[TIGS$Final.Result=="Negative" & is.na(TIGS$Viral.Load..copies.ml.)]<-0
TIGS1$Viral.Load..copies.ml.[TIGS1$Final.Result=="Negative" & TIGS1$Viral.Load..copies.ml.>5]<-5
TIGS1$Viral.Load..copies.ml.[TIGS1$Final.Result=="Negative" & TIGS1$Viral.Load..copies.ml.==0]<-1

#Calculating EWMA

lamda<-0.70

#calculating the log
TIGS1<-TIGS1 %>% mutate(log_viral_loads=log(Viral.Load..copies.ml.))

TIGS1<- TIGS1 %>% group_by(STP) %>% mutate(ewma = accumulate(log_viral_loads, ~lamda*.y+(1 - lamda)*.x))


#Merging the final and the old data set

TIGS_merged<-merge(x=TIGS,y=TIGS1,by=c("Sample.ID","Date","Final.Result","Viral.Load..copies.ml.","STP","date_onset"),all.x = T)
write.csv(TIGS_merged,"C:/Users/LENOVO/Downloads/exported data/tigs_EWMA_merged_final.csv")

#Renaming variables

TIGS_merged<-TIGS_merged %>% rename(Lab_ID=Lab.ID,
                                    Final_result=Final.Result,
                                    Viral_load_copies_ml=Viral.Load..copies.ml.)
write.csv(TIGS_merged,"C:/Users/LENOVO/Downloads/exported data/tigs_EWMA_merged1.csv")

#trying dcast to convert to wide format

TIGS_merged_wide<- dcast(TIGS1, date_onset+Lab.ID+Rnase.P+E.Gene+RdRp.Gene~STP+Final.Result, value.var = "Viral.Load..copies.ml.")

