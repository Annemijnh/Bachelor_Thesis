# Packages ----------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(iglu)
library(data.table)
library(ggcorrplot)
library(caret)
library(runner)
library(Boruta)
library(randomForest)

# Loading data ------------------------------------------------------------

# cgm_data <- read.table('cgm_data.csv',header = TRUE,sep= ',')
# cgm_data$id <- factor(cgm_data$id)
# cgm_data$time <- as.POSIXct(cgm_data$time)
# activity_data <- read.table('activity_data-Rekenmachine.csv',sep = ',',header = TRUE)
# activity_data$id <- factor(activity_data$id)
# activity_data$time <- as.POSIXct(activity_data$time)
# inclusion <- read.csv('inclusion.csv')

# If other scripts were run first, import step can be skipped for inclusion, activity_data and cgm_data
# But do run this part, otherwise skip:
activity_data <- act_data_selection
cgm_data <- cgm_selection

# Import of subject specific data
subject_data <- fread("Datanaarcastorcompleet.csv",select=c(1,20, 21, 28, 34, 35, 36, 39, 40,41,42,44, 45, 172, 275, 277, 280, 282, 283, 285, 286, 287, 288, 290), fill = TRUE,data.table = FALSE)
subject_data <- subject_data[(subject_data$`Participant Id` %in% inclusion),] # Only select the included subjects
colnames(subject_data) <- c('id', 'sex_female', 'age','DM2_years', 'smoking_current','smoking_ever','packyears','alcohol_units_month', 'length', 'weight', 'BMI', 'waist_circumference', 'hip_circumference','HbA1c', 'fast_insulin', 'fast_insulin_dosage', 'mix_insulin', 'mix_insulin_dosage', 'long_insulin', 'long_insulin_dosage','metformin','metformin_dosage', 'SU_derivatives', 'SU_derivatives_dosage')
# BMI was not calculated for each patient
subject_data <- distinct(subject_data) |>
  mutate(BMI = weight/((length/100)^2))


# Feature extraction ------------------------------------------------------

# First, time-based features
cgm_data_time <- cgm_data|>
  mutate(date = date(time),
         year = year(time),
         month = month(time),
         mday = mday(time),
         wday = lubridate::wday(time, label = TRUE, locale = "English_United Kingdom"),
         wday = factor(wday,ordered = FALSE),
         weekend = as.factor(case_when(wday == "Sat" |wday == "Sun" ~ "Weekend",
                                       wday == "Mon" |wday == "Tue" | wday == "Wed" |wday == "Thu" |wday == "Fri" ~ "Week")),
         hour = hour(time),
         minute = minute(time),
         part_of_day = as.factor(case_when( hour >= 6 & hour < 12 ~ "Morning",
                                            hour >= 12 & hour < 18 ~ "Afternoon",
                                            hour >= 18 & hour < 24 ~ "Evening",
                                            hour >= 0 ~ "Night")),
         season = as.factor(case_when(month >= 3 & month < 6 ~ "Spring",
                                      month >= 6 & month < 9 ~ "Summer",
                                      month >= 9 & month < 12 ~ "Autumn",
                                      TRUE ~ "Winter")))

# Next, CGM based features
cgm_data_cgm <- cgm_data_time|>
  arrange(id,time)|>
  group_by(id)|>
  mutate(slope = gl - lag(gl),
         previous_slope = lag(slope),
         gl_minus_1 = lag(gl),
         gl_minus_2 = lag(gl,n=2),
         gl_minus_3 = lag(gl,n=3),
         gl_minus_4 = lag(gl,n=4),
         gl_minus_5 = lag(gl,n=5),
         gl_minus_6 = lag(gl,n=6),
         gl_previous_day = runner(gl,idx = time,k = 'day', f = first, na_pad = TRUE),
         gl_previous_hour = runner(gl,idx = time,k = 'hour', f = first, na_pad = TRUE),
         rolling_av_1hour = runner(x = gl, mean, k= '1 hour', na_pad = TRUE, lag = 1, idx = time),
         rolling_sd_1hour = runner(x = gl, sd, k= '1 hour', na_pad = TRUE, lag = 1, idx = time),
         rolling_av_2hour = runner(x = gl, mean, k= '2 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_sd_2hour = runner(x = gl, sd, k= '2 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_av_3hour = runner(x = gl, mean, k= '3 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_sd_3hour = runner(x = gl, sd, k= '3 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_av_4hour = runner(x = gl, mean, k= '4 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_sd_4hour = runner(x = gl, sd, k= '4 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_av_5hour = runner(x = gl, mean, k= '5 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_sd_5hour = runner(x = gl, sd, k= '5 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_av_6hour = runner(x = gl, mean, k= '6 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_sd_6hour = runner(x = gl, sd, k= '6 hours', na_pad = TRUE, lag = 1, idx = time),
         rolling_av_day = runner(x = gl, mean, k= '1 day', na_pad = TRUE, lag = 1, idx = time),
         rolling_sd_day = runner(x = gl, sd, k= '1 day', na_pad = TRUE, lag = 1, idx = time),
         # The following functions are based on iglu or gluvarpro package functions. I only needed the essence of the function, not all dplyr tidying since I'm doing that myself
         rolling_iqr_day = runner(x = gl, function(x) IQR(x), k = '1 day', na_pad = TRUE, lag = 1, idx = time),
         rolling_min_day = runner(gl, function(x) min(x), idx = time, lag = 1, k = 'day', na_pad = TRUE),
         rolling_max_day = runner(gl, function(x) max(x), idx = time, lag = 1, k = 'day', na_pad = TRUE),
         rolling_TAR_day = runner(x = gl, function(x) {mean(x > 10, na.rm = TRUE) * 100}, k = '1 day', lag = 1, idx = time, na_pad= TRUE),
         rolling_TIR_day = runner(x = gl, function(x) {mean(x >= 3.9 & x <= 10, na.rm = TRUE) * 100}, k = '1 day', lag = 1, idx = time, na_pad = TRUE),
         rolling_TBR_day = runner(x = gl, function(x) {mean(x < 3.9, na.rm = TRUE)*100}, k = '1 day', lag = 1, idx = time, na_pad = TRUE),
         rolling_TAHR_day = runner(x = gl, function(x) {mean(x > 13.9, na.rm = TRUE) * 100}, k = '1 day', lag = 1, idx = time, na_pad= TRUE),
         rolling_TBLR_day = runner(x = gl, function(x) {mean(x < 3, na.rm = TRUE)*100}, k = '1 day', lag = 1, idx = time, na_pad = TRUE),
         rolling_LBGI_day = runner(x = gl, function(x){
           f_vec <- 1.509 * ( log(x*18)^1.084 - 5.381 )
           rl <- ifelse(f_vec < 0, 10 * f_vec^2,0)
           lbgi <- mean(rl,na.rm = TRUE)
           return(lbgi)}
         , k = '1 day', lag = 1, idx = time, na_pad= TRUE),
         rolling_HBGI_day = runner(x = gl, function(x){
           f_vec <- 1.509 * ( log(x*18)^1.084 - 5.381 )
           rh <- ifelse(f_vec > 0, 10 * f_vec^2,0)
           hbgi <- mean(rh,na.rm = TRUE)
           return(hbgi)}, k = '1 day', lag = 1, idx = time, na_pad= TRUE),
         cv_past_day = rolling_sd_day/rolling_av_day*100,
         J_index_day = 0.324*(rolling_av_day + rolling_sd_day)^2,
         rolling_mage = runner(x = gl, function(x){
           times <- time
             
             smoothed_values <- vector()
             smoothed_values[1:4] <- mean(stats::na.omit(x[1:4]))
             for (i in 5:(length(x)-4)){
               smoothed_values[i] <- (1*x[i - 4] + 2*x[i - 3] + 4*x[i - 2] + 8*x[i - 1] + 16*x[i] +
                                        8*x[i + 1] + 4*x[i + 2] + 2*x[i + 3] + 1*x[i + 4]) / 46
             }
             ti <- length(smoothed_values)
             smoothed_values[(ti + 1):(ti + 4)] <- mean(x[(length(x) - 3):length(x)])
             sd <- stats::sd(x)
             
             u <- unique(smoothed_values)
             n <- length(u)
             v1 <- c(0,u[2:n] - u[1:(n - 1)])
             v2 <- c(u[2:n] - u[1:(n - 1)],0)
             
             peaks <- which(v1 > 0 & v2 < 0)
             troughs <- which(v1 < 0 & v2 > 0)
             if (length(peaks) > length(troughs)) troughs <- append(troughs, 1)
             if (length(peaks) < length(troughs)) peaks <- append(peaks, 1)
             
             differences <- u[peaks] - u[troughs]
             mean(stats::na.omit(differences[which(differences > sd)]))
           
         }, k = '1 day', lag = 1, idx = time, na_pad= TRUE)
         
         )|>
select(-c(slope))

activity_data <- activity_data|>
  mutate(hour = hour(time),
         minute = minute(time),
         date = date(time))|> #id, date, hour, and minute so I can merge with the CGM dataset 
  group_by(id)|>
  arrange(id,time)|>
  mutate(steps_past_hour = runner(steps, f = sum, idx = time, k = 'hour',lag = 1,na_pad = TRUE))|>
  mutate(steps_past_day = runner(steps, f = sum, idx = time, k = 'day',lag = 1,na_pad = TRUE))|>
  select('id','date','hour','minute','steps_past_hour','steps_past_day') #steps per minute is useless, since the CGM data usually has /15 minutes

# Combine data ------------------------------------------------------------
# All NA in CGM or activity data are removed and NA in subject data is set to 0. 
total_data <- merge(cgm_data_cgm,activity_data,by=c('id','date','hour','minute'))|>
  arrange(id,time)
no_na_total_data <- total_data|>
  na.omit()|>
  select(-c(time, minute, year, date, mday, month))
no_na_clean_data <- merge(no_na_total_data,subject_data,by='id')
no_na_clean_data <- replace(no_na_clean_data, is.na(no_na_clean_data), 0)


# Personal set ------------------------------------------------------------
# Determine which patients to include:
adrr_cgm <- cgm_data 
adrr_cgm$gl <- adrr_cgm$gl*18 # Convert to mg/dL
adrr_cgm <- adrr_cgm|>
  group_by(id) |>
  group_split()
adrr <- list()
for(i in 1:length(adrr_cgm)){
  adrr[[i]] <- adrr(adrr_cgm[[i]])
}
adrr <- bind_rows(adrr)
low_risk_id <- adrr$id[which(min(adrr$ADRR) == adrr$ADRR)]
moderate_risk_id <- adrr$id[which(median(adrr$ADRR) == adrr$ADRR)]
high_risk_id <- adrr$id[which(max(adrr$ADRR) == adrr$ADRR)]

# The low-risk patient had a lot more data, so the next lowest was chosen: 629                                  
                                  
pers_model_set <- no_na_clean_data
sub_629 <- subset(pers_model_set, id == 629)|> # Low variablility
  select(-id)
dummy <- dummyVars(" ~ .", data=sub_629)
sub_629 <- data.frame(predict(dummy, newdata = sub_629))|>
  select(-c(weekend.Week))
sub_688 <- subset(pers_model_set, id == 688)|> # High variability
  select(-id)
dummy <- dummyVars(" ~ .", data=sub_688)
sub_688 <- data.frame(predict(dummy, newdata = sub_688))|>
  select(-c(weekend.Week))
sub_708 <- subset(pers_model_set, id == 708)|> # Median ADRR
  select(-id)
dummy <- dummyVars(" ~ .", data=sub_708)
sub_708 <- data.frame(predict(dummy, newdata = sub_708))|>
  select(-c(weekend.Week))                                  

# Population set ------------------------------------------------------------                                  
pop_model_set <- no_na_clean_data|>
  select(-id)
dummy <- dummyVars(" ~ .", data=pop_model_set)
pop_model_set <- data.frame(predict(dummy, newdata = pop_model_set))
pop_model_set <- pop_model_set|>
  select(-c(weekend.Week))                                  
