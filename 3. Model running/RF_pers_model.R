library(ega)
library(tidyverse)
library(caret)


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

# Personal model ----------------------------------------------------------

grid_RF629 <- expand.grid(mtry =seq(40,78,5),
                          splitrule = "extratrees",
                          min.node.size =  seq(1,5,2))
grid_RF688 <- expand.grid(mtry =seq(53,78,3),
                          splitrule = "extratrees",
                          min.node.size =  seq(1,5,2))
grid_RF708 <- expand.grid(mtry =seq(53,78,3),
                          splitrule = "extratrees",
                          min.node.size =  seq(1,5,2))
metrics_629_RF <- data.frame()
for (i in seed2){
  print(i)
  set.seed(i)
  train_index_629 <- unlist(createDataPartition(sub_629$gl, p = 0.8))
  train_629 <- sub_629[train_index_629, ]
  test_629 <- sub_629[-train_index_629, ]
  
  set.seed(seed)
  RF_629_default <- train(gl~., data = train_629, method = 'ranger', metric = 'RMSE', importance = 'permutation', trControl = trcv_10)
  set.seed(seed)
  RF_629_random <- train(gl~., data = train_629, method = 'ranger', metric = 'RMSE', importance = 'permutation',tuneLength = 6, trControl = trcv_10_random)
  set.seed(seed)
  RF_629_grid <- train(gl~., data = train_629, method = 'ranger', metric = 'RMSE', importance = 'permutation',tuneGrid = grid_RF629, trControl = trcv_10_grid)
  
  RMSE_RF_629_default <-  RMSE(test_629$gl,predict(RF_629_default,test_629))
  RMSE_RF_629_random <-  RMSE(test_629$gl,predict(RF_629_random,test_629))
  RMSE_RF_629_grid <-  RMSE(test_629$gl,predict(RF_629_grid,test_629))
  RMSE_RF_629 <- data.frame(RMSE = rbind(RMSE_RF_629_default,
                                         RMSE_RF_629_random,
                                         RMSE_RF_629_grid))
  RMSE_RF_629$names <- rownames(RMSE_RF_629)
  
  if (which.min(RMSE_RF_629$RMSE) == 1){best_model = RF_629_default
  }  else if(which.min(RMSE_RF_629$RMSE) == 2){best_model = RF_629_random
  }  else {best_model = RF_629_grid}
  
  RMSE_RF_629 <- RMSE_RF_629[which.min(RMSE_RF_629$RMSE),]
  metric <- RMSE_RF_629|> 
    mutate(MAE = MAE(test_629$gl,predict(best_model,test_629)),
           type = gsub("RMSE_RF_629_","",RMSE_RF_629$names),
           model = "RF",
           dataset = "629",
           seed = as.character(i))
  zones <- getClarkeZones(test_629$gl,predict(best_model,test_629), unit = "mol")
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
  metrics_629_RF <- rbind(metrics_629_RF,metric)
  best_RF_629 <- metrics_629_RF|>
    arrange(RMSE)
  if(best_RF_629$seed[1] == i){ # If the current model is the best, save the test, prediction and varImp!
    best_RF_629_model = best_model
    best_test_RF_629 = test_629
    best_pred_RF_629 = predict(best_RF_629_model,best_test_RF_629)
    var_imp_RF_629 = varImp(best_RF_629_model)$importance
  }
} 
RMSE(best_pred_RF_629,best_test_RF_629$gl) # This is to test whether the above code worked

# Save the CEG of the best model
png('CEG_LR_RF.png', height = 450, width = 750)
plotClarkeGrid_adj(best_test_RF_629$gl, best_pred_RF_629, unit = 'mol',title = 'Low-Risk, Random Forest')
dev.off()

var_imp_RF_629$names <- rownames(var_imp_RF_629)
var_imp_RF_629 <- var_imp_RF_629|>
  mutate(Risk = 'Low')


metrics_708_RF <- data.frame()
for (i in seed2){
  print(i)
  set.seed(i)
  train_index_708 <- unlist(createDataPartition(sub_708$gl, p = 0.8))
  train_708 <- sub_708[train_index_708, ]
  test_708 <- sub_708[-train_index_708, ]
  
  set.seed(seed)
  RF_708_default <- train(gl~., data = train_708, method = 'ranger', metric = 'RMSE', importance = 'permutation', trControl = trcv_10)
  set.seed(seed)
  RF_708_random <- train(gl~., data = train_708, method = 'ranger', metric = 'RMSE', importance = 'permutation',tuneLength = 6, trControl = trcv_10_random)
  set.seed(seed)
  RF_708_grid <- train(gl~., data = train_708, method = 'ranger', metric = 'RMSE', importance = 'permutation',tuneGrid = grid_RF708, trControl = trcv_10_grid)
  
  RMSE_RF_708_default <-  RMSE(test_708$gl,predict(RF_708_default,test_708))
  RMSE_RF_708_random <-  RMSE(test_708$gl,predict(RF_708_random,test_708))
  RMSE_RF_708_grid <-  RMSE(test_708$gl,predict(RF_708_grid,test_708))
  RMSE_RF_708 <- data.frame(RMSE = rbind(RMSE_RF_708_default,
                                         RMSE_RF_708_random,
                                         RMSE_RF_708_grid))
  RMSE_RF_708$names <- rownames(RMSE_RF_708)
  
  if (which.min(RMSE_RF_708$RMSE) == 1){best_model = RF_708_default
  }  else if(which.min(RMSE_RF_708$RMSE) == 2){best_model = RF_708_random
  }  else {best_model = RF_708_grid}
  RMSE_RF_708 <- RMSE_RF_708[which.min(RMSE_RF_708$RMSE),]
  
  metric <- RMSE_RF_708|> 
    mutate(MAE = MAE(test_708$gl,predict(best_model,test_708)),
           type = gsub("RMSE_RF_708_","",RMSE_RF_708$names),
           model = "RF",
           dataset = "708",
           seed = as.character(i))
  zones <- getClarkeZones(test_708$gl,predict(best_model,test_708), unit = "mol")
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
  metrics_708_RF <- rbind(metrics_708_RF,metric)
  best_RF_708 <- metrics_708_RF|>
    arrange(RMSE)
  if(best_RF_708$seed[1] == i){ # If the current model is the best, save the test, prediction and varImp!
    best_RF_708_model = best_model
    best_test_RF_708 = test_708
    best_pred_RF_708 = predict(best_RF_708_model,best_test_RF_708)
    var_imp_RF_708 = varImp(best_RF_708_model)$importance
  }
} 
png('CEG_MR_RF.png', height = 450, width = 750)
plotClarkeGrid_adj(best_test_RF_708$gl, best_pred_RF_708, unit = 'mol',title = 'Moderate-Risk, Random Forest')
dev.off()

var_imp_RF_708$names <- rownames(var_imp_RF_708)
var_imp_RF_708 <- var_imp_RF_708|>
  arrange(-Overall)|>
  mutate(Risk = 'Moderate')

metrics_688_RF <- data.frame()
for (i in seed2){
  print(i)
  set.seed(i)
  train_index_688 <- unlist(createDataPartition(sub_688$gl, p = 0.8))
  train_688 <- sub_688[train_index_688, ]
  test_688 <- sub_688[-train_index_688, ]
  
  set.seed(seed)
  RF_688_default <- train(gl~., data = train_688, method = 'ranger', metric = 'RMSE', importance = 'permutation', trControl = trcv_10)
  set.seed(seed)
  RF_688_random <- train(gl~., data = train_688, method = 'ranger', metric = 'RMSE', importance = 'permutation',tuneLength = 6, trControl = trcv_10_random)
  set.seed(seed)
  RF_688_grid <- train(gl~., data = train_688, method = 'ranger', metric = 'RMSE', importance = 'permutation',tuneGrid = grid_RF688, trControl = trcv_10_grid)
  
  RMSE_RF_688_default <-  RMSE(test_688$gl,predict(RF_688_default,test_688))
  RMSE_RF_688_random <-  RMSE(test_688$gl,predict(RF_688_random,test_688))
  RMSE_RF_688_grid <-  RMSE(test_688$gl,predict(RF_688_grid,test_688))
  RMSE_RF_688 <- data.frame(RMSE = rbind(RMSE_RF_688_default,
                                         RMSE_RF_688_random,
                                         RMSE_RF_688_grid))
  RMSE_RF_688$names <- rownames(RMSE_RF_688)
  
  if (which.min(RMSE_RF_688$RMSE) == 1){best_model = RF_688_default
  }  else if(which.min(RMSE_RF_688$RMSE) == 2){best_model = RF_688_random
  }  else {best_model = RF_688_grid}
  RMSE_RF_688 <- RMSE_RF_688[which.min(RMSE_RF_688$RMSE),]
  
  metric <- RMSE_RF_688|> 
    mutate(MAE = MAE(test_688$gl,predict(best_model,test_688)),
           type = gsub("RMSE_RF_688_","",RMSE_RF_688$names),
           model = "RF",
           dataset = "688",
           seed = as.character(i))
  zones <- getClarkeZones(test_688$gl,predict(best_model,test_688), unit = "mol")
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
  metrics_688_RF <- rbind(metrics_688_RF,metric)
  best_RF_688 <- metrics_688_RF|>
    arrange(RMSE)
  if(best_RF_688$seed[1] == i){ # If the current model is the best, save the test, prediction and varImp!
    best_RF_688_model = best_model
    best_test_RF_688 = test_688
    best_pred_RF_688 = predict(best_RF_688_model,best_test_RF_688)
    var_imp_RF_688 = varImp(best_RF_688_model)$importance
  }
} 
png('CEG_HR_RF.png', height = 450, width = 750)
plotClarkeGrid_adj(best_test_RF_688$gl, best_pred_RF_688, unit = 'mol',title = 'High-Risk, Random Forest')
dev.off()

var_imp_RF_688$names <- rownames(var_imp_RF_688)
var_imp_RF_688 <- var_imp_RF_688|>
  arrange(-Overall)|>
  mutate(Risk = 'High')

var_imp_RF_pers <- rbind(var_imp_RF_629,var_imp_RF_708,var_imp_RF_688)|>
  mutate(Risk = factor(Risk, levels = c('Low', 'Moderate', 'High'),ordered = T))|>
  group_by(names)|>
  mutate(sum = sum(Overall))|>
  arrange(-sum)|>
  head(45)|> 
  mutate(category = case_when((names %in% time_based_features)~"Time",
                              (names %in% previous_glucose_measurement_features)~ "Previous glucose measurements",
                              (names %in% simple_glucose_calculation_features)~"Basic glucose calculations",
                              (names %in% glycemic_variability_features)~"Glycemic variability indices",
                              (names %in% baseline_subject_features)~"Baseline subject characteristics",
                              (names %in% subject_medication_features)~"Medication",
                              (names %in% activity_features)~"Activity"))
png('RF_pers.png', height = 450, width = 750)
ggplot(var_imp_RF_pers, aes(x=Overall, y=reorder(names,sum))) +
  geom_bar(stat="identity", aes(fill=category)) +                                               
  scale_fill_manual(values=cbPalette,
    name="Category",
    breaks=c("Time", "Previous glucose measurements", "Basic glucose calculations", "Glycemic variability indices",
             "Baseline subject characteristics", "Medication","Activity"),
    labels = c("Time", "Previous glucose measurements", "Basic glucose calculations", "Glycemic variability indices",
               "Baseline subject characteristics", "Medication","Activity"))+
  geom_bar(stat="identity", aes(alpha=Risk)) + 
  scale_alpha_manual(values=c(0.1,0.5,0.8),
                     name="Risk",
                     breaks=c("Low", "Moderate", "High"),        
                     labels=c("Low", "Moderate", "High")) +     
  theme_bw(16) +
  scale_x_continuous(expand = c(0,0)) +                                                          
  labs(x = 'Feature Importance', y = 'Features', title = 'Random Forest')+
  theme(legend.position="bottom")
dev.off()

save.image("C:/Users/annem/OneDrive - University of Twente/Bacheloropdracht/Data/RF_pers_models.RData")
write.csv(metrics_629_RF, "metrics_629_RF.csv")
write.csv(metrics_708_RF,"metrics_708_RF.csv")
write.csv(metrics_688_RF, "metrics_688_RF.csv")
