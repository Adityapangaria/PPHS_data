library(tidyverse)
library(lubridate)
library(ggplot2)
PPHS_data<-read.csv("C:/Users/LENOVO/Downloads/PPHS.csv")
TPR<-read.csv("C:/Users/LENOVO/Downloads/TPR.csv")
TIGS<-read.csv("C:/Users/LENOVO/Downloads/TIGS_EWMA.csv")



# All the 3 datasets loaded

#removing entries with invalid results


class(PPHS_data$Date)
PPHS_data <- PPHS_data %>% mutate(date_onset = as.Date(Date, format = "%d/%m/%Y"))
class(PPHS_data$date_onset)
PPHS_data$Viral.Load..EWMA.<- as.numeric(PPHS_data$Viral.Load..EWMA.)

#Removed invalid entries

PPHS_data<-PPHS_data[!(PPHS_data$Final.Result=="Invalid"),]

#EWMA summed up datewise

PPHS_data1<- PPHS_data %>% group_by(date_onset) %>% summarise(EWMA_datewise= sum(Viral.Load..EWMA., na.rm = T))

#data cleaning for the positivity data till 3rd march

TPR <- TPR %>% mutate(date_onset = as.Date(Date, format = "%d/%m/%Y"))
class(TPR$date_onset)
TPR_new<- TPR %>% group_by(date_onset) %>% summarise(Total_test_daywise= sum(Total.Tests, na.rm = T), Total_positive_daywise= sum(Total.positives, na.rm = T))
TPR_new<- TPR_new %>% mutate(positivity_daywise= Total_positive_daywise/Total_test_daywise)

#merged tpr and ewma data till 

merged_tpr_ewma<- merge(x = PPHS_data1, y = TPR_new, by = "date_onset")

# tested data for normality
shapiro.test(merged_tpr_ewma$EWMA_datewise)$p             # Data is not normally distributed
shapiro.test(merged_tpr_ewma$positivity_daywise)$p        # Data is not normally distributed

# since the data is not normally distributed using the Spearman correlation method

cor.test(merged_tpr_ewma$EWMA_datewise, merged_tpr_ewma$positivity_daywise, method = "spearman")

# Correlation through the perason method

cor.test(merged_tpr_ewma$EWMA_datewise, merged_tpr_ewma$positivity_daywise)

# p-value less than 0.05 we can reject the null hypothesis of zero correlation
# rho of 0.66 shows that there is a moderate correlation at the daily level between positivity rate and EWMA from OD.



# TIGS data cleaning

class(TIGS$Date)
TIGS <- TIGS %>% mutate(date_onset = as.Date(Date, format = "%d/%m/%Y"))
TIGS<-TIGS[!(TIGS$Final.Result=="Inconclusive"),]
TIGS1<- TIGS %>% group_by(date_onset) %>% summarise(EWMA_datewise_TIGS= sum(Viral.Load..copies.ml..EWMA), na.rm = T)

# Merging TIGS Viral loads data and TPR

merged_TIGS_TPR<- merge(x = TIGS1, y = TPR_new, by = "date_onset")
shapiro.test(merged_TIGS_TPR$EWMA_datewise_TIGS)$p             # Data is not normally distributed
shapiro.test(merged_tpr_ewma$positivity_daywise)$p        # Data is not normally distributed

#CAlculating correlations between TIGS viral loads and Positivity through Spearman and Perason respectively

cor.test(merged_TIGS_TPR$EWMA_datewise_TIGS, merged_TIGS_TPR$positivity_daywise, method = "spearman")
cor.test(merged_TIGS_TPR$EWMA_datewise_TIGS, merged_TIGS_TPR$positivity_daywise)

#for the weekly counts of data

PPHS_data_weekly <- PPHS_data %>% drop_na(date_onset) %>% mutate(weekly_cases = floor_date(date_onset, unit = "week"))
PPHS_data_weekly1<- PPHS_data_weekly %>% group_by(weekly_cases) %>% summarise(EWMA_weekly= sum(Viral.Load..EWMA., na.rm = T))

#for weekly count of TPR


TPR_data_weekly <- TPR %>% drop_na(date_onset) %>% mutate(weekly_cases = floor_date(date_onset, unit = "week"))
TPR_new_weekly<- TPR_data_weekly %>% group_by(weekly_cases) %>% summarise(Total_test_weekwise= sum(Total.Tests, na.rm = T), Total_positive_weekwise= sum(Total.positives, na.rm = T))
TPR_new_weekly<- TPR_new_weekly %>% mutate(positivity_weekwise= Total_positive_weekwise/Total_test_weekwise)

#merged data for weekly calculations on OD and TPR

merged_tpr_ewma_weekly<- merge(x = PPHS_data_weekly1, y = TPR_new_weekly, by = "weekly_cases")
shapiro.test(merged_tpr_ewma_weekly$EWMA_weekly)$p    
shapiro.test(merged_tpr_ewma_weekly$positivity_weekwise)$p 

#test for correlation between the positivity rate weekly and EWMA weekly from OD

cor.test(merged_tpr_ewma_weekly$EWMA_weekly, merged_tpr_ewma_weekly$positivity_weekwise, method = "spearman")

# rho of 0.6890 shows moderate correlation between the weekly data for OD which is slightly higher than the daywise data

cor.test(merged_tpr_ewma_weekly$EWMA_weekly, merged_tpr_ewma_weekly$positivity_weekwise)
#Weekly correlation between the number of positives and viral load for OD
cor.test(merged_tpr_ewma_weekly$EWMA_weekly, merged_tpr_ewma_weekly$Total_positive_weekwise, method = "spearman")
#PEarson correlation
cor.test(merged_tpr_ewma_weekly$EWMA_weekly, merged_tpr_ewma_weekly$Total_positive_weekwise)

# For TIGS data

TIGS_weekly<- TIGS %>% drop_na(date_onset) %>% mutate(weekly_cases = floor_date(date_onset, unit = "week"))
TIGS_new_weekly<- TIGS_weekly %>% group_by(weekly_cases) %>% summarise(EWMA_weekly_TIGS= sum(Viral.Load..copies.ml..EWMA, na.rm = T))

#merged TIGS and TPR weekly data

merged_TIGS_TPR_weekly<- merge(x=TIGS_new_weekly, y=TPR_new_weekly, by="weekly_cases")

#correlation between TIGS and TPR weekly data

cor.test(merged_TIGS_TPR_weekly$EWMA_weekly_TIGS, merged_TIGS_TPR_weekly$positivity_weekwise, method = "spearman")

#High degree of correlation found on weekly basis between TIGS and TPR with rho value of 0.73

#Correlation between the total positives and the EMWA fromm the STP

cor.test(merged_TIGS_TPR_weekly$EWMA_weekly_TIGS, merged_TIGS_TPR_weekly$Total_positive_weekwise, method = "spearman")
cor.test(merged_TIGS_TPR_weekly$EWMA_weekly_TIGS, merged_TIGS_TPR_weekly$Total_positive_weekwise)
cor.test(tigs_od_merged$EWMA_weekly_TIGS, tigs_od_merged$EWMA_weekly)


#Checking correlations till Mid Feb for OD vs Positivity 
merged_tpr_ewma_weekly %>% filter(weekly_cases >= as.Date('2021-05-09') & weekly_cases <= as.Date('2022-02-13')) %>%
                          ggplot(aes(y= positivity_weekwise,x= EWMA_weekly))+geom_smooth(method = "lm", col="black")+geom_point(col="red")+
                          labs(title = "Correlation between OD and weekly positivity", y="Weekly Positivity in %", x="Log of EWMA Viral Load(no. of copies/ml)", size=8)+theme_classic()+
                          scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
                          scale_x_continuous(trans = "log10")+
                          annotate("text", x=1000, y=0.15, label="R=0.91 p-value= 2.2e-16  Method= Pearson", size=5, col="black")
                        
#Correlation and graphs till mid feb for OD and number of Positives
merged_tpr_ewma_weekly %>% filter(weekly_cases >= as.Date('2021-05-09') & weekly_cases <= as.Date('2022-02-13')) %>%
                           ggplot(aes(EWMA_weekly,Total_positive_weekwise))+geom_smooth(method="lm", col="black")+geom_point(col="red")+
                           labs(title = "Correlation between OD and number of positives", y="Number of positives", x="Log of EWMA Viral Load(no. of copies/ml)")+
                           theme_classic()+ scale_y_continuous(trans = "log10")+scale_x_continuous(trans = "log10")+
                           annotate("text", x=1000, y=120000, label="R=0.94 p-value= 2.5e-14  Method= Pearson", size=5, col="black")
                            
#Checking correlations till Mid Feb for STP vs Poitive cases
tigs_od_merged<-merge(merged_TIGS_TPR_weekly,merged_tpr_ewma_weekly, by="weekly_cases")
tigs_of_merged_new<-tigs_od_merged %>% filter(weekly_cases >= as.Date('2021-05-09') & weekly_cases <= as.Date('2022-02-13'))
cor.test(tigs_of_merged_new$EWMA_weekly_TIGS, tigs_of_merged_new$Total_positive_weekwise.x)
cor.test(tigs_of_merged_new$EWMA_weekly_TIGS, tigs_of_merged_new$EWMA_weekly)
cor.test(tigs_of_merged_new$EWMA_weekly, tigs_of_merged_new$Total_positive_weekwise.x)
#Checking correlations till Mid Feb for OD vs STP
tigs_od_merged %>% filter(weekly_cases >= as.Date('2021-05-09') & weekly_cases <= as.Date('2022-02-13')) %>%
                  ggplot(aes(EWMA_weekly_TIGS,EWMA_weekly))+geom_smooth(method="lm", col="black")+theme_classic()+
                  labs(title = "Correlation between the viral loads from OD and STP", x="Viral load from STP weekly(number of copies/ml)", y="Viral load from OD(number of copies/ml)")+
                  scale_x_continuous(trans = "log10")+scale_y_continuous(trans="log10")+geom_point(col="red")+
                  annotate("text", x=45, y=100000, label="R=0.94  p-value=2.2e-14   Method= Pearson", size=5, col="black")

tigs_od_merged %>% filter(weekly_cases >= as.Date('2021-05-09') & weekly_cases <= as.Date('2022-02-13')) %>%
  ggplot(aes(EWMA_weekly_TIGS,Total_positive_weekwise.x))+geom_smooth(method="lm", col="black")+theme_classic()+geom_point(col="red")+
  labs(title = "Correlation between the viral loads from STP and Total Positive cases", x="Viral load from STP weekly(number of copies/ml)", y="Number of positive cases")+
  scale_x_continuous(trans = "log10")+scale_y_continuous(trans="log10")+
  annotate("text", x=45, y=100000, label="R=0.98  p-value=2.2e-16   Method= Pearson", size=5)
#merging the weekly dataset TIGS and OD with their TPR

tigs_od_merged<-merge(merged_TIGS_TPR_weekly, merged_tpr_ewma_weekly, by="weekly_cases")

#Calculating the OD and STP viral laods

cor.test(tigs_od_merged$EWMA_weekly_TIGS, tigs_od_merged$EWMA_weekly)

# Calculating correlation through Spearman

cor.test(tigs_od_merged$EWMA_weekly_TIGS, tigs_od_merged$EWMA_weekly, method = "spearman")


od_tpr_daily_feb_17<-merged_tpr_ewma[merged_tpr_ewma$date_onset >= "2021/05/15" & merged_tpr_ewma$date_onset <= "2022/02/17", ]
od_tpr11<-merged_tpr_ewma %>% filter(EWMA_datewise>="2021/05/15" & EWMA_datewise<="2022/02/17")
#correlation between positivity and viral loads till 17th Feb

cor.test(od_tpr_new_daily$EWMA_datewise, od_tpr_new_daily$positivity_daywise, method = "spearman")
cor.test(od_tpr_new_weekly$EWMA_weekly, od_tpr_new_weekly$Total_positive_weekwise, method = "spearman")
cor.test(merged_TIGS_TPR_weekly$EWMA_weekly_TIGS, merged_TIGS_TPR_weekly$Total_positive_weekwise, method = "spearman")
cor.test()
OD_ewma_weekly_corr_jan<- ggplot(od_tpr_new_weekly, aes(EWMA_weekly, positivity_weekwise))+geom_line()+geom_smooth(method="lm")+annotate("text", x = 50000, y = 0.25, 
                                                                                                                                                 label = "R=0.7291 corr OD EWMA positivity rate weekly(Jan)",
                                                                                                                                                 col = "blue", size = 5)
tpr_ewma_daily_jan<-ggplot(od_tpr_new_weekly,aes(x=weekly_cases))+geom_line(aes(y=EWMA_weekly),color="blue")+geom_line(aes(y=positivity_weekwise),color="red")+scale_y_continuous(trans = "log10")
class(od_tpr_new_weekly$weekly_cases)
tpr_ewma_daily_jan
cor.test(od_tpr_new_weekly$EWMA_weekly, od_tpr_new_weekly$Total_positive_weekwise, method="spearman")
cor.test(merged_TIGS_TPR_weekly$EWMA_weekly_TIGS, merged_TIGS_TPR_weekly$Total_positive_weekwise, method="spearman")
cor.test(od_tpr_new_daily$EWMA_datewise, od_tpr_new_daily$Total_positive_daywise, method= "spearman")
sd(od_tpr_new_daily$Total_test_daywise)

tpr_EMWA_daily_corr_graph<- ggplot(od_tpr_new_daily, aes(EWMA_datewise, positivity_daywise))+geom_line()+geom_smooth(method="lm")+ labs(title = "Correlation EWMA from OD and Positivity rate daily")+labs(x= "EWMA daily in number of copies/ml")+labs(y="Positivity daily in %")+theme_bw()
tpr_EMWA_daily_corr_graph+scale_y_continuous(labels = function(x) paste0(x*100, "%"))+annotate("text", x = 2000, y = .20, 
                                                                                               label = "Rho value=0.6582",
                                                                                               col = "blue", size = 5)
tpr_EMWA_weekly_corr_graph<- ggplot(od_tpr_new_weekly, aes(y=positivity_weekwise, x=EWMA_weekly))+geom_line()+geom_smooth(method="lm")+ labs(title = "Correlation EWMA from OD and Positivity rate weekly")+labs(x= "EWMA weekly")+labs(y="Positivity weekly")+theme_bw()+scale_y_continuous(labels = function(x) paste0(x*100, "%"))
tpr_EMWA_weekly_corr_graph+annotate("text", x = 10000, y = .20, 
                                    label = "Rho value=0.7291",
                                    col = "blue", size = 5)
ewma_tpr_daily_graph<- ggplot(od_tpr_new_daily, aes(y=positivity_daywise, x= EWMA_datewise))+geom_line()+geom_smooth(method="lm")+labs(title = "Rho of 0.6582 between EWMA from OD and Positivity rate daily")+labs(x= "positivity daily")+labs(y="log EWMA daily")
ewma_tpr_daily_graph
ewma_plot<-ggplot(od_tpr_new_daily, aes(date_onset, EWMA_datewise))+geom_line()+labs(title = "EWMA viral loads daily", x="Date", y="EWMA viral load(number of copies/ml)")
ewma_plot<-ewma_plot+scale_x_date(date_breaks = "1 month", date_labels="%B")+theme_bw()
ewma_plot+geom_text(x=as.Date("2021-06-14"),y=10000, label="End of second wave", size=5)+geom_vline(xintercept=as.Date("2021-06-14"), col="black")+geom_vline(xintercept = as.Date("2022-01-03",col="black"))+geom_text(x=as.Date("2022-01-03"),y=10000, label="Start of third wave", size=5)
ewma_weekly_plot<-ggplot(od_tpr_new_weekly, aes(weekly_cases, EWMA_weekly))+geom_line()+labs(title = "EWMA viral loads weekly", x="Date", y="EWMA viral load(number of copies/ml)")
ewma_weekly_plot<-ewma_weekly_plot+scale_x_date(date_breaks = "1 month",date_labels="%B")+theme_bw()
ewma_weekly_plot<-ewma_weekly_plot+geom_text(x=as.Date("2021-06-01"),y=35000, label="End of second wave", size=5)+geom_vline(xintercept = as.Date("2021-06-01"), col="black")
ewma_weekly_plot+geom_text(x=as.Date("2021-12-25"),y=35000, label="Start of third wave", size=5)+geom_vline(xintercept =as.Date("2021-12-25") )
positivity_weekly_plot<-ggplot(od_tpr_new_weekly, aes(weekly_cases,positivity_weekwise))+geom_line()+labs(title = "Weekly Positivity rate", x="Date", y="Weekly positivity rates in %")
positivity_weekly_plot<-ggplot(od_tpr_new_weekly, aes(weekly_cases,positivity_weekwise))+geom_line()+labs(title = "Weekly Positivity rate", x="Date", y="Weekly positivity rates in %")ggplot(od_tpr_new_weekly, aes(weekly_cases,positivity_weekwise))+geom_line()+labs(title = "Weekly Positivity rate", x="Date", y="Weekly positivity rates in %")ggplot(od_tpr_new_weekly, aes(weekly_cases,positivity_weekwise))+geom_line()+labs(title = "Weekly Positivity rate", x="Date", y="Weekly positivity rates in %")
positivity_weekly_plot<-positivity_weekly_plot+geom_text(x=as.Date("2021-06-06"),y=0.05, label="End of second wave", size=5)+geom_vline(xintercept =as.Date("2021-06-06"), col="black")
positivity_weekly_plot+geom_text(x=as.Date("2021-12-19"),y=0.05, label="Start of third wave", size=5)+geom_vline(xintercept =as.Date("2022-01-02"), col="black")+theme_bw()
positivity_daily_plot<-ggplot(od_tpr_new_daily, aes(date_onset,positivity_daywise))+geom_line()+labs(title = "Daily Positivity rate", x="Date", y="Daily positivity rates in %")
positivity_daily_plot<-positivity_daily_plot+scale_x_date(date_breaks = "1 month",date_labels="%B")+theme_bw()
positivity_daily_plot<-positivity_daily_plot+scale_y_continuous(labels = function(x) paste0(x*100, "%"))
positivity_daily_plot<-positivity_daily_plot+geom_text(x=as.Date("2021-06-06"),y=0.05, label="End of second wave", size=5)+geom_vline(xintercept =as.Date("2021-06-06"), col="black")
positivity_daily_plot<-positivity_daily_plot+geom_text(x=as.Date("2021-12-19"),y=0.05, label="Start of third wave", size=5)+geom_vline(xintercept =as.Date("2022-01-02"), col="black")
summary(od_tpr_new_daily$Total_test_daywise)
sd(od_tpr_new_daily$Total_test_daywise)
offff<-od_tpr_new_daily %>% group_by(date_onset) %>% count(Total_test_daywise==0)
rm(offff)
cor.test(od_tpr_new_daily$EWMA_datewise, od_tpr_new_daily$Total_positive_daywise, method="spearman", exact = F)
cor.test(od_tpr_new_weekly$EWMA_weekly, od_tpr_new_weekly$Total_positive_weekwise, method="spearman")
cor.test(od_tpr_new_weekly$EWMA_weekly, od_tpr_new_weekly$positivity_weekwise, method="spearman")
cor.test(od_tpr_new_daily$EWMA_datewise, od_tpr_new_daily$positivity_daywise, method="spearman", exact=F)



daily_ewma_tpr_jan_df<-merged_tpr_ewma[merged_tpr_ewma$date_onset >= "2022-01-03" & merged_tpr_ewma$date_onset <= "2022-03-03", ]
weekly_ewma_tpr_jan_df<-merged_tpr_ewma_weekly[merged_tpr_ewma_weekly$weekly_cases >= "2022-01-02" & merged_tpr_ewma_weekly$weekly_cases <= "2022-02-27", ]


positivity_daily_plot_jan<- ggplot(daily_ewma_tpr_jan_df, aes(date_onset,positivity_daywise))+theme_bw()+geom_line()+labs(title = "Daily Positivity Jan-Mar 2022", X="Date", y="Positivity in %")
positivity_daily_plot_jan+scale_x_date(date_breaks = "1 month",date_labels="%B")+scale_y_continuous(labels = function(x) paste0(x*100, "%"))
positivity_weekly_plot_jan<- ggplot(weekly_ewma_tpr_jan_df, aes(weekly_cases, positivity_weekwise))+geom_line()+labs(title = "Weekly positivity Jan-Mar 2022", x="Date", y="Positivity in %")
positivity_weekly_plot_jan+scale_y_continuous(labels = function(x) paste0(x*100, "%"))+scale_x_date(date_breaks = "1 month",date_labels="%B")+theme_bw()
cor.test(weekly_ewma_tpr_jan_df$EWMA_weekly, weekly_ewma_tpr_jan_df$Total_positive_weekwise)
rm(od_tpr11)
rm(od_tpr_daily_feb_17)
rm(od_tpr_Feb20)
