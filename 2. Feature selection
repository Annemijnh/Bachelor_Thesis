# Packages ----------------------------------------------------------------
library(Boruta)
library(randomForest)

# Make sure all scripts in 1. are run. 

# Run Boruta
seed = 123
set.seed(seed)
boruta <- Boruta(gl ~ ., data = pop_model_set, doTrace = 3, maxRuns = 100)

# Features categorized
time_based_features <- c('hour',"wday.Sun","wday.Mon","wday.Tue","wday.Wed","wday.Thu","wday.Fri",
                         "wday.Sat","weekend.Weekend","part_of_day.Afternoon", "part_of_day.Evening",
                         "part_of_day.Morning" ,  "part_of_day.Night" ,    "season.Autumn",
                         "season.Spring" ,"season.Summer","season.Winter")
previous_glucose_measurement_features <- c("gl_minus_1", "gl_minus_2","gl_minus_3", "gl_minus_4",
                                           "gl_minus_5", "gl_minus_6","gl_previous_day" ,"gl_previous_hour")
simple_glucose_calculation_features <- c("rolling_av_1hour" ,"rolling_sd_1hour","rolling_av_2hour",
                                         "rolling_sd_2hour","rolling_av_3hour","rolling_sd_3hour",
                                         "rolling_av_4hour","rolling_sd_4hour","rolling_av_5hour" ,
                                         "rolling_sd_5hour" ,"rolling_av_6hour","rolling_sd_6hour" ,
                                         "rolling_av_day", "rolling_min_day","rolling_max_day", "previous_slope" )
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

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Categorize Boruta output, get top 30
boruta_output <- attStats(boruta)|>
  rownames_to_column("feature")|>
  mutate(category = case_when((feature %in% time_based_features)~"Time",
                              (feature %in% previous_glucose_measurement_features)~ "Previous glucose measurements",
                              (feature %in% simple_glucose_calculation_features)~"Basic glucose calculations",
                              (feature %in% glycemic_variability_features)~"Glycemic variability indices",
                              (feature %in% baseline_subject_features)~"Baseline subject characteristics",
                              (feature %in% subject_medication_features)~"Medication",
                              (feature %in% activity_features)~"Activity"))|>
  arrange(-meanImp)|>
  head(30)

# Plot
ggplot(data = boruta_output, 
       aes(x = meanImp, y = reorder(feature, meanImp))) +
  geom_bar(stat="identity", aes(fill=category), width = 0.8) + 
  
  scale_fill_manual(#values = c('#FF8585',"#FFC000","#FFFF47","#92D050","#5DE4FF","#0192FF","#C7A1E3"),
    values=cbPalette,
    name="Category",
    breaks=c("Time", "Previous glucose measurements", "Basic glucose calculations", "Glycemic variability indices",
             "Baseline subject characteristics", "Medication","Activity"),
    labels = c("Time", "Previous glucose measurements", "Basic glucose calculations", "Glycemic variability indices",
               "Baseline subject characteristics", "Medication","Activity"))+
  geom_errorbar(aes(xmin=minImp, xmax=maxImp), width=.5,
                position=position_dodge(.9))+
  theme_bw(base_size = 16)+
  scale_x_continuous(expand = c(0,0))+
  labs(x = 'Feature Importance', y = 'Features', title = 'Boruta Output')+
  theme(legend.position="bottom")
