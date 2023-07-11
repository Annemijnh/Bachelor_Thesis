library(ega)
library(tidyverse)
library(caret)

seed <- 123
seed2 <- list(123,456,789,012,345,678,901,234,567,890)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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

grid_RF <- expand.grid(mtry =seq(42,78,12),
                       splitrule = "extratrees",
                       min.node.size = 5)

metrics_pop_RF <- data.frame()



# Rerun from here ---------------------------------------------------------
# This takes a very long time to run, so I recommend adding save states of the workspace inbetween, so quitting R does not make you lose any data
for (i in seed2){
  print(i)
  set.seed(i)
  train_index <- unlist(createDataPartition(pop_model_set$gl, p = 0.8))
  train <- pop_model_set[train_index,]
  test <- pop_model_set[-train_index,]
  
  set.seed(seed)
  RF_default <-  train(gl~., data = train, method = 'ranger', metric = 'RMSE', importance = 'permutation', trControl = trcv_10)
  set.seed(seed)
  RF_random <- train(gl~., data = train, method = 'ranger', metric = 'RMSE', importance = "permutation",tuneLength = 3, trControl = trcv_10_random)
  set.seed(seed)
  RF_grid <- train(gl~., data = train, method = 'ranger', metric = 'RMSE', importance = 'permutation',tuneGrid = grid_RF, trControl = trcv_10_grid)
  RMSE_RF_default <-  RMSE(test$gl,predict(RF_default,test))
  RMSE_RF_random <-  RMSE(test$gl,predict(RF_random,test))
  RMSE_RF_grid <-  RMSE(test$gl,predict(RF_grid,test))
  RMSE_RF <- data.frame(RMSE = rbind(RMSE_RF_default,
                                     RMSE_RF_random,
                                     RMSE_RF_grid))
  RMSE_RF$names <- rownames(RMSE_RF)
  
  if (which.min(RMSE_RF$RMSE) == 1){best_model = RF_default
  }  else if(which.min(RMSE_RF$RMSE) == 2){best_model = RF_random
  }  else {best_model = RF_grid}
  
  RMSE_RF <- RMSE_RF[which.min(RMSE_RF$RMSE),]
  metric <- RMSE_RF|>
    mutate(MAE = MAE(test$gl,predict(best_model,test)),
           type = gsub("RMSE_RF_","",RMSE_RF$names),
           model = "RF",
           dataset = "pop",
           seed = as.character(i))
  zones <- getClarkeZones(test$gl,predict(best_model,test), unit = "mol")
  CEG_zones <- as.data.frame(round(table(zones)/ length(zones) * 100, digits = 2))|>
    pivot_wider(names_from = zones, values_from = Freq)|>
    mutate(A = ifelse(!("A" %in% zones), 0, A),
           B = ifelse(!("B" %in% zones), 0, B),
           C = ifelse(!("C" %in% zones), 0, C),
           D = ifelse(!("D" %in% zones), 0, D),
           E = ifelse(!("E" %in% zones), 0, E))|>
    mutate(CEG_AB = A + B,
           CEG_CDE = C + D + E)|>
    select(-c(A,B,C,D,E))
  metric <- cbind(metric,CEG_zones)
  metrics_pop_RF <- rbind(metrics_pop_RF,metric)
  best_RF <- metrics_pop_RF|>
    arrange(RMSE)
  if(best_RF$seed[1] == i){
    best_RF_model = best_model
    best_test_RF = test
    best_pred_RF = predict(best_RF_model,best_test_RF)
    var_imp_RF_pop = varImp(best_RF_model)$importance
  }
rm(RF_default,RF_grid,RF_random)
}

png('CEG_RF.png', height = 450, width = 750)
plotClarkeGrid_adj(best_test_RF$gl, best_pred_RF, unit = 'mol',title = 'Random Forest')
dev.off()

var_imp_RF_pop$names <- rownames(var_imp_RF_pop)
var_imp_RF <- var_imp_RF_pop|>
  arrange(-Overall)|>
  head(15)|>
  mutate(category = case_when((names %in% time_based_features)~"Time",
                              (names %in% previous_glucose_measurement_features)~ "Previous glucose measurements",
                              (names %in% simple_glucose_calculation_features)~"Basic glucose calculations",
                              (names %in% glycemic_variability_features)~"Glycemic variability indices",
                              (names %in% baseline_subject_features)~"Baseline subject characteristics",
                              (names %in% subject_medication_features)~"Medication",
                              (names %in% activity_features)~"Activity"))
png('RF_pop.png', height = 450, width = 750)
ggplot(data = var_imp_RF,
       aes(x = Overall, y = reorder(names,Overall))) +
  geom_bar(stat="identity", aes(fill=category), width = 0.8) +
  scale_fill_manual(
    values=cbPalette,
    name="Category",
    breaks=c("Time", "Previous glucose measurements", "Basic glucose calculations", "Glycemic variability indices",
             "Baseline subject characteristics", "Medication","Activity"),
    labels = c("Time", "Previous glucose measurements", "Basic glucose calculations", "Glycemic variability indices",
               "Baseline subject characteristics", "Medication","Activity"))+
  theme_bw(base_size = 16)+
  scale_x_continuous(expand = c(0,0))+
  labs(x = 'Feature Importance', y = 'Features', title = 'Random Forest')+
  theme(legend.position="bottom")
dev.off()

write.csv(metrics_pop_RF, "metrics_pop_RF.csv")
