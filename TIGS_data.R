library(tidyverse)
library(lubridate)
library(ggplot2)
library(corrplot)
library(reshape2)
install.packages("RColorBrewer")
library(zoo)
library(RColorBrewer)


TIGS<-read.csv("C:/Users/LENOVO/Downloads/TIGS data.csv")



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

#Calculating EWMA

lamda<-0.70



TIGS1<- TIGS1 %>% group_by(STP) %>% mutate(ewma = accumulate(Viral.Load..copies.ml., ~lamda*.y+(1 - lamda)*.x))
rm(TIGS2)

#Merging the final and the old data set

TIGS_merged<-merge(x=TIGS,y=TIGS1,by=c("Date","STP","Lab.ID","Final.Result","Viral.Load..copies.ml.","Rnase.P","E.Gene","RdRp.Gene","N.Gene","date_onset"),all.x = T)
TIGS_merged2<- TIGS_merged %>% group_by(date_onset,STP) %>% mutate(ewma = accumulate(Viral.Load..copies.ml., ~ lamda * .y + (1 - lamda) * .x))
write.csv(TIGS_merged_weekly,"C:/Users/LENOVO/Downloads/exported data/tigs_merged_weekly_cases.csv")

#Renaming variables

TIGS_merged<-TIGS_merged %>% rename(Lab_ID=Lab.ID,
                                    Final_result=Final.Result,
                                    Viral_load_copies_ml=Viral.Load..copies.ml.)
write.csv(corrdata,"C:/Users/LENOVO/Downloads/exported data/tigs_correlations.csv")


#Converting daily to weekly viral loads

TIGS_merged_weekly<-TIGS_merged %>% drop_na(Date) %>% mutate(weekly_cases = floor_date(Date, unit = "10 days"))

# Calculating the weekly emwa for each site

TIGS_merged_weekly<-TIGS_merged_weekly %>% group_by(weekly_cases,STP,Final_result,Lab_ID) %>% summarise(Viral_load_weekly= sum(Viral_load_copies_ml, na.rm = T))

# Decasting the data to wide

TIGS_merged_wide<- dcast(TIGS_merged_weekly, weekly_cases~STP, value.var = "Viral_load_weekly", fun.aggregate = mean)

# replacing NAN with the mean of the coloumn


TIGS_merged_wide_1<- TIGS_merged_wide %>% select("Agaram":"Yele Mallappa Chetti Kere")
TIGS_merged_wide_2<- TIGS_merged_wide %>% select("weekly_cases")
TIGS_merged_wide_1<-na.aggregate(TIGS_merged_wide_1)
#Combining merged 1 and merged 2
TIGS_merged_final<- cbind(TIGS_merged_wide_1, TIGS_merged_wide_2)
m<-cor(TIGS_merged_wide_1)
corrplot(m, type="full", order="hclust",
col=brewer.pal(n=10, name="RdYlBu"))
corrdata<-as.data.frame(m)
