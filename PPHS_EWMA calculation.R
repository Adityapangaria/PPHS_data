library(tidyverse)
library(lubridate)
library(ggplot2)
library(corrplot)
library(reshape2)


PPHS<-read.csv("C:/Users/LENOVO/Downloads/PPHS2.csv")



# PPHS data cleaning

class(PPHS$Date)
class(PPHS$Viral.Load)
PPHS <- PPHS %>% mutate(date_onset = as.Date(Date, format = "%d/%m/%Y"))
PPHS$Viral.Load<-as.numeric(PPHS$Viral.Load)

#Removing inconclusive and invalid entries

PPHS1<-PPHS[!(PPHS$Final.Result=="Invalid"),]

#marking entries with Begative result and NA  viral loads as zero

PPHS1$Viral.Load[PPHS1$Final.Result=="Negative" & is.na(PPHS1$Viral.Load)]<-0
PPHS1$Viral.Load[PPHS1$Final.Result=="Positive" & is.na(PPHS1$Viral.Load)]<-0
PPHS$Viral.Load[PPHS$Final.Result=="Positive" & is.na(PPHS$Viral.Load)]<-0
PPHS$Viral.Load[PPHS$Final.Result=="Negative" & is.na(PPHS$Viral.Load)]<-0

#Calculating EWMA

lamda<-0.70
PPHS1<- PPHS1 %>% group_by(Zone) %>% mutate(ewma = accumulate(Viral.Load, ~ lamda * .y + (1 - lamda) * .x))

#Merging the final and the old data set

PPHS_merged<-merge(x=PPHS,y=PPHS1,by=c("Lab.ID..","Zone", "Date","Final.Result", "Viral.Load", "Viral.Load..EWMA.", "Rnase.P", "E.Gene", "RdRp.Gene", "Orf1ab", "N.GENE", "X","Latitude", "Longitude", "X.1", "Latitude_1","Longitude_1","X.2","X.3","X.4","X.5", "date_onset"),all.x = T)


write.csv(PPHS_merged,"C:/Users/LENOVO/Downloads/exported data/PPHS_ewma_data_merged.csv")


