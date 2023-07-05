library(tidyverse)
library(ggcorrplot)

time_based_features <- c('hour',"wday.Sun","wday.Mon","wday.Tue","wday.Wed","wday.Thu","wday.Fri",
                         "wday.Sat","weekend.Weekend","part_of_day.Afternoon", "part_of_day.Evening",
                         "part_of_day.Morning" ,  "part_of_day.Night" ,    "season.Autumn",
                         "season.Spring" ,"season.Summer","season.Winter")
previous_glucose_measurement_features <- c("gl_minus_1", "gl_minus_2","gl_minus_3", "gl_minus_4",
                                           "gl_minus_5", "gl_minus_6" ,"gl_previous_hour","gl_previous_day")
simple_glucose_calculation_features <- c("rolling_av_1hour" ,"rolling_av_2hour","rolling_av_3hour",
                                         "rolling_av_4hour","rolling_av_5hour" ,"rolling_av_6hour","rolling_av_day",
                                         "rolling_sd_1hour","rolling_sd_2hour","rolling_sd_3hour",
                                         "rolling_sd_4hour","rolling_sd_5hour" ,"rolling_sd_6hour",
                                          "rolling_min_day","rolling_max_day", "previous_slope")
glycemic_variability_features <- c("rolling_sd_day","rolling_iqr_day","rolling_TAR_day",
                                   "rolling_TIR_day","rolling_TBR_day","rolling_TAHR_day",
                                   "rolling_TBLR_day","rolling_LBGI_day","rolling_HBGI_day",
                                   "cv_past_day","J_index_day","rolling_mage")
baseline_subject_features <- c("sex_female","age", "DM2_years","smoking_current","smoking_ever",
                               "packyears","alcohol_units_month","length","weight","BMI","waist_circumference",
                               "hip_circumference","HbA1c")
subject_medication_features <- c("fast_insulin","fast_insulin_dosage","mix_insulin",
                                 "mix_insulin_dosage","long_insulin","long_insulin_dosage",
                                 "metformin", "metformin_dosage","SU_derivatives","SU_derivatives_dosage")
activity_features <- c("steps_past_hour","steps_past_day")

corr_time <- pop_model_set|> select(one_of(time_based_features))
cor_matrix_time <- cor(corr_time,use="pairwise.complete.obs")
png('Corr_time.png', height = 450, width = 750)
ggcorrplot(cor_matrix_time,show.diag=F, type="lower", lab_size=2)+
  theme_light(16)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "",title = "Time-based feature correlation")
dev.off()

corr_gluc <- pop_model_set|> select(one_of(previous_glucose_measurement_features))
cor_matrix_gluc <- cor(corr_gluc,use="pairwise.complete.obs")
png('Corr_prevgluc.png', height = 450, width = 750)
ggcorrplot(cor_matrix_gluc,show.diag=F, type="lower", lab_size=2)+
  theme_light(16)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "",title = "Previous glucose measurement features correlation")
dev.off()

corr_basic <- pop_model_set|> select(one_of(simple_glucose_calculation_features))
cor_matrix_basic <- cor(corr_basic,use="pairwise.complete.obs")
png('Corr_basicgluc.png', height = 450, width = 750)
ggcorrplot(cor_matrix_basic,show.diag=F, type="lower", lab_size=2)+
  theme_light(16)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "",title = "Basic glucose calculation features correlation")
dev.off()

corr_glycvar <- pop_model_set|> select(one_of(glycemic_variability_features))
cor_matrix_glycvar <- cor(corr_glycvar,use="pairwise.complete.obs")
png('Corr_glycvar.png', height = 450, width = 750)
ggcorrplot(cor_matrix_glycvar,show.diag=F, type="lower", lab_size=2)+
  theme_light(16)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "",title = "Glycemic variability features correlation")
dev.off()

corr_sub <- pop_model_set|> select(one_of(baseline_subject_features))
cor_matrix_sub <- cor(corr_sub,use="pairwise.complete.obs")
png('Corr_sub.png', height = 450, width = 750)
ggcorrplot(cor_matrix_sub,show.diag=F, type="lower", lab_size=2)+
  theme_light(16)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "",title = "Baseline subject features correlation")
dev.off()

corr_med <- pop_model_set|> select(one_of(subject_medication_features))
cor_matrix_med <- cor(corr_med,use="pairwise.complete.obs")
png('Corr_med.png', height = 450, width = 750)
ggcorrplot(cor_matrix_med,show.diag=F, type="lower", lab_size=2)+
  theme_light(16)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "",title = "Medication features correlation")
dev.off()

corr_act <- pop_model_set|> select(one_of(activity_features))
cor_matrix_act <- cor(corr_act,use="pairwise.complete.obs")
png('Corr_act.png', height = 450, width = 750)
ggcorrplot(cor_matrix_act,show.diag=F, type="lower", lab_size=2)+
  theme_light(16)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "",title = "Activity features correlation")
dev.off()


corr <- pop_model_set
cor_matrix <- cor(corr,use="pairwise.complete.obs")
# Figure below is basically messy, so in the report I'm only reporting the correlation over 90%
png('Corr.png', height = 450, width = 750)
ggcorrplot(cor_matrix,show.diag=F, type="lower", lab_size=2)+
  theme_light(16)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "",title = "All features correlation")
dev.off()

vars <- findCorrelation(cor(pop_model_set),cutoff = 0.90)
colnames(pop_model_set[,vars])
