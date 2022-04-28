library(tidyverse)
install.packages("ggstatsplot")
library(ggstatsplot)
library(lubridate)
install.packages("qcc")
install.packages("mdsstat")
library(qcc)
library(mdsstat)
PPHS_data<-read.csv("C:/Users/LENOVO/Downloads/PPHS.csv")
TPR<-read.csv("C:/Users/LENOVO/Downloads/TPR.csv")
TIGS<-read.csv("C:/Users/LENOVO/Downloads/TIGS_EWMA.csv")
class(PPHS_data$Date)
#removing entries with invalid results
as.date

class(PPHS_data$Date)
PPHS_data <- PPHS_data %>% mutate(date_onset = as.Date(Date, format = "%d/%m/%Y"))
class(PPHS_data$date_onset)
PPHS_data$Viral.Load..EWMA.<- as.numeric(PPHS_data$Viral.Load..EWMA.)

#Removed invalid entries

PPHS_data<-PPHS_data[!(PPHS_data$Final.Result=="Invalid"),]






PPHS_data <- PPHS_data %>% mutate(date_onset = as.Date(Date, format = "%d/%m/%Y"))
class(PPHS_data$date_onset)
PPHS_data$Viral.Load..EWMA.<- as.numeric(PPHS_data$Viral.Load..EWMA.)
PPHS_data_new1<- PPHS_data %>% group_by(date_onset) %>% summarise(Viral_load_datewise= sum(Viral.Load, na.rm = T))
lamda<-0.70

#Writing the EWMA function
df <- PPHS_data_new1 %>%
      mutate(ewma = accumulate(Viral_load_datewise, ~ lamda * .y + (1 - lamda) * .x))

#data cleaning for the positivity data till 3rd march

TPR <- TPR %>% mutate(date_onset = as.Date(Date, format = "%d/%m/%Y"))
class(TPR$date_onset)
TPR_new<- TPR %>% group_by(date_onset) %>% summarise(Total_test_daywise= sum(Total.Tests, na.rm = T), Total_positive_daywise= sum(Total.positives, na.rm = T))
TPR_new<- TPR_new %>% mutate(positivity_daywise= Total_positive_daywise/Total_test_daywise)

#merged tpr and ewma data till 

merged_tpr_ewma<- merge(x = df, y = TPR_new, by = "date_onset")
#for the weekly counts of data

PPHS_data_weekly <- PPHS_data %>% drop_na(date_onset) %>% mutate(weekly_cases = floor_date(date_onset, unit = "week"))
PPHS_data_weekly1<- PPHS_data_weekly %>% group_by(weekly_cases) %>% summarise(Viral_load_weekly= sum(Viral.Load, na.rm = T))

#Writing the EWMA function

df_weekly <- PPHS_data_weekly1 %>%
  mutate(ewma_weekl = accumulate(Viral_load_weekly, ~ lamda * .y + (1 - lamda) * .x))

#for weekly count of TPR


TPR_data_weekly <- TPR %>% drop_na(date_onset) %>% mutate(weekly_cases = floor_date(date_onset, unit = "week"))
TPR_new_weekly<- TPR_data_weekly %>% group_by(weekly_cases) %>% summarise(Total_test_weekwise= sum(Total.Tests, na.rm = T), Total_positive_weekwise= sum(Total.positives, na.rm = T))
TPR_new_weekly<- TPR_new_weekly %>% mutate(positivity_weekwise= Total_positive_weekwise/Total_test_weekwise)
#merging weekly datasets

merged_tpr_ewma_weekly<- merge(x = df_weekly, y = TPR_new_weekly, by = "weekly_cases")

write.csv(merged_tpr_ewma_weekly, "C:/Users/LENOVO/Downloads/exported data/merged_tpr_ewma_weekly.csv")

cor.test(merged_tpr_ewma_weekly$ewma_weekl, merged_tpr_ewma_weekly$Total_positive_weekwise, method = "spearman")
cor.test(merged_tpr_ewma$ewma, merged_tpr_ewma$Total_positive_daywise, method = "spearman")
