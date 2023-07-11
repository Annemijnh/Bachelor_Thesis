library(tidyverse)
library(caret)
library(ggplot2)

seed = 123
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
grid_lr <- expand.grid(intercept = c(FALSE,TRUE))

metrics_629_LR <- data.frame()
for (i in seed2){
  print(i)
  set.seed(i)
  train_index_629 <- unlist(createDataPartition(sub_629$gl, p = 0.8))
  train_629 <- sub_629[train_index_629, ]
  test_629 <- sub_629[-train_index_629, ]
  
  set.seed(seed)
  LR_629_default <-  train(gl~., data = train_629, method = 'lm', metric = 'RMSE', trControl = trcv_10)
  set.seed(seed)
  LR_629_grid <- train(gl~., data = train_629, method = 'lm', metric = 'RMSE', trControl = trcv_10_grid,tuneGrid  = grid_lr)
  
  RMSE_LR_629_default <-  RMSE(test_629$gl,predict(LR_629_default,test_629))
  RMSE_LR_629_grid <-  RMSE(test_629$gl,predict(LR_629_grid,test_629))
  RMSE_LR_629 <- data.frame(RMSE = rbind(RMSE_LR_629_default,
                                         RMSE_LR_629_grid))
  RMSE_LR_629$names <- rownames(RMSE_LR_629)
  
  if (which.min(RMSE_LR_629$RMSE) == 1){best_model = LR_629_default
  } else {best_model = LR_629_grid}
                        
  RMSE_LR_629 <- RMSE_LR_629[which.min(RMSE_LR_629$RMSE),]
  metric <- RMSE_LR_629|> 
    mutate(MAE = MAE(test_629$gl,predict(best_model,test_629)),
           type = gsub("RMSE_LR_629_","",RMSE_LR_629$names),
           model = "LR",
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
  metrics_629_LR <- rbind(metrics_629_LR,metric)
  best_LR_629 <- metrics_629_LR|>
    arrange(RMSE)
  if(best_LR_629$seed[1] == i){ # If the current model is the best, save the test, prediction and varImp!
    best_LR_629_model = best_model
    best_test_LR_629 = test_629
    best_pred_LR_629 = predict(best_LR_629_model,best_test_LR_629)
    var_imp_LR_629 = varImp(best_LR_629_model)$importance
  }
} 
RMSE(best_pred_LR_629,best_test_LR_629$gl) # This is to test whether the above code worked

metrics_708_LR <- data.frame()
for (i in seed2){
  print(i)
  set.seed(i)
  train_index_708 <- unlist(createDataPartition(sub_708$gl, p = 0.8))
  train_708 <- sub_708[train_index_708, ]
  test_708 <- sub_708[-train_index_708, ]
  
  set.seed(seed)
  LR_708_default <-  train(gl~., data = train_708, method = 'lm', metric = 'RMSE', trControl = trcv_10)
  set.seed(seed)
  LR_708_grid <- train(gl~., data = train_708, method = 'lm', metric = 'RMSE', trControl = trcv_10_grid,tuneGrid  = grid_lr)
  
  RMSE_LR_708_default <-  RMSE(test_708$gl,predict(LR_708_default,test_708))
  RMSE_LR_708_grid <-  RMSE(test_708$gl,predict(LR_708_grid,test_708))
  RMSE_LR_708 <- data.frame(RMSE = rbind(RMSE_LR_708_default,
                                         RMSE_LR_708_grid))
  RMSE_LR_708$names <- rownames(RMSE_LR_708)
  
  if (which.min(RMSE_LR_708$RMSE) == 1){best_model = LR_708_default
  }  else {best_model = LR_708_grid}
  
  RMSE_LR_708 <- RMSE_LR_708[which.min(RMSE_LR_708$RMSE),]
  metric <- RMSE_LR_708|> 
    mutate(MAE = MAE(test_708$gl,predict(best_model,test_708)),
           type = gsub("RMSE_LR_708_","",RMSE_LR_708$names),
           model = "LR",
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
  metrics_708_LR <- rbind(metrics_708_LR,metric)
  best_LR_708 <- metrics_708_LR|>
    arrange(RMSE)
  if(best_LR_708$seed[1] == i){
    best_LR_708_model = best_model
    best_test_LR_708 = test_708
    best_pred_LR_708 = predict(best_LR_708_model,best_test_LR_708)
    var_imp_LR_708 = varImp(best_LR_708_model)$importance
  }
} 

metrics_688_LR <- data.frame()
for (i in seed2){
  print(i)
  set.seed(i)
  train_index_688 <- unlist(createDataPartition(sub_688$gl, p = 0.8))
  train_688 <- sub_688[train_index_688, ]
  test_688 <- sub_688[-train_index_688, ]
  
  set.seed(seed)
  LR_688_default <-  train(gl~., data = train_688, method = 'lm', metric = 'RMSE', trControl = trcv_10)
  set.seed(seed)
  LR_688_grid <- train(gl~., data = train_688, method = 'lm', metric = 'RMSE', trControl = trcv_10_grid,tuneGrid  = grid_lr)
  
  RMSE_LR_688_default <-  RMSE(test_688$gl,predict(LR_688_default,test_688))
  RMSE_LR_688_grid <-  RMSE(test_688$gl,predict(LR_688_grid,test_688))
  RMSE_LR_688 <- data.frame(RMSE = rbind(RMSE_LR_688_default,
                                         RMSE_LR_688_grid))
  RMSE_LR_688$names <- rownames(RMSE_LR_688)
  
  if (which.min(RMSE_LR_688$RMSE) == 1){best_model = LR_688_default
  }  else {best_model = LR_688_grid}
  RMSE_LR_688 <- RMSE_LR_688[which.min(RMSE_LR_688$RMSE),]
  
  metric <- RMSE_LR_688|> 
    mutate(MAE = MAE(test_688$gl,predict(best_model,test_688)),
           type = gsub("RMSE_LR_688_","",RMSE_LR_688$names),
           model = "LR",
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
  metrics_688_LR <- rbind(metrics_688_LR,metric)
  best_LR_688 <- metrics_688_LR|>
    arrange(RMSE)
  if(best_LR_688$seed[1] == i){
    best_LR_688_model = best_model
    best_test_LR_688 = test_688
    best_pred_LR_688 = predict(best_LR_688_model,best_test_LR_688)
    var_imp_LR_688 = varImp(best_LR_688_model)$importance
  }
} 

# Population model --------------------------------------------------------

metrics_pop_LR <- data.frame()
for (i in seed2){
  print(i)
  set.seed(i)
  train_index <- unlist(createDataPartition(pop_model_set$gl, p = 0.8))
  train <- pop_model_set[train_index,]
  test <- pop_model_set[-train_index,]
  
  set.seed(seed)
  LR_default <-  train(gl~., data = train, method = 'lm', metric = 'RMSE', trControl = trcv_10)
  set.seed(seed)
  LR_grid <- train(gl~., data = train, method = 'lm', metric = 'RMSE', trControl = trcv_10_grid,tuneGrid  = grid_lr)
  
  RMSE_LR_default <-  RMSE(test$gl,predict(LR_default,test))
  RMSE_LR_grid <-  RMSE(test$gl,predict(LR_grid,test))
  RMSE_LR <- data.frame(RMSE = rbind(RMSE_LR_default,
                                         RMSE_LR_grid))
  RMSE_LR$names <- rownames(RMSE_LR)
  
  if (which.min(RMSE_LR$RMSE) == 1){best_model = LR_default
  }  else {best_model = LR_grid}
  
  RMSE_LR <- RMSE_LR[which.min(RMSE_LR$RMSE),]
  metric <- RMSE_LR|> 
    mutate(MAE = MAE(test$gl,predict(best_model,test)),
           type = gsub("RMSE_LR_","",RMSE_LR$names),
           model = "LR",
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
  metrics_pop_LR <- rbind(metrics_pop_LR,metric)
  best_LR <- metrics_pop_LR|>
    arrange(RMSE)
  if(best_LR$seed[1] == i){
    best_LR_model = best_model
    best_test_LR = test
    best_pred_LR = predict(best_LR_model,best_test_LR)
    var_imp_LR_pop = varImp(best_LR_model)$importance
  }
} 

# Output --------------------------------------------------------
# Run plotClarkeGrid_adj.R first!
png('CEG_LR_LR.png', height = 450, width = 750)
plotClarkeGrid_adj(best_test_LR_629$gl, best_pred_LR_629, unit = 'mol',title = 'Low-Risk, Linear Regression')
dev.off()

var_imp_LR_629$names <- rownames(var_imp_LR_629)
var_imp_LR_629 <- var_imp_LR_629|>
  add_row(names = 'rolling_TAHR_day', Overall = 0)|>
  arrange(-Overall)|>
  mutate(Risk = 'Low')

png('CEG_MR_LR.png', height = 450, width = 750)
plotClarkeGrid_adj(best_test_LR_708$gl, best_pred_LR_708, unit = 'mol',title = 'Moderate-Risk, Linear Regression')
dev.off()

var_imp_LR_708$names <- rownames(var_imp_LR_708)
var_imp_LR_708 <- var_imp_LR_708|>
  arrange(-Overall)|>
  mutate(Risk = 'Moderate')

png('CEG_HR_LR.png', height = 450, width = 750)
plotClarkeGrid_adj(best_test_LR_688$gl, best_pred_LR_688, unit = 'mol',title = 'High-Risk, Linear Regression')
dev.off()

var_imp_LR_688$names <- rownames(var_imp_LR_688)
var_imp_LR_688 <- var_imp_LR_688|>
  arrange(-Overall)|>
  mutate(Risk = 'High')

var_imp_LR_pers <- rbind(var_imp_LR_629,var_imp_LR_708,var_imp_LR_688)|>
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

png('LR_pers.png', height = 450, width = 750)
ggplot(var_imp_LR_pers, aes(x=Overall, y=reorder(names,sum))) +
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
  labs(x = 'Feature Importance', y = 'Features', title = 'Linear Regression')+
  theme(legend.position="bottom")
dev.off()

png('CEG_LR.png', height = 450, width = 750)
plotClarkeGrid_adj(best_test_LR$gl, best_pred_LR, unit = 'mol',title = 'Linear Regression')
dev.off()

var_imp_LR_pop$names <- rownames(var_imp_LR_pop)
var_imp_LR <- var_imp_LR_pop|>
  arrange(-Overall)|>
  head(15)|>
  mutate(category = case_when((names %in% time_based_features)~"Time",
                              (names %in% previous_glucose_measurement_features)~ "Previous glucose measurements",
                              (names %in% simple_glucose_calculation_features)~"Basic glucose calculations",
                              (names %in% glycemic_variability_features)~"Glycemic variability indices",
                              (names %in% baseline_subject_features)~"Baseline subject characteristics",
                              (names %in% subject_medication_features)~"Medication",
                              (names %in% activity_features)~"Activity"))
png('LR_pop.png', height = 450, width = 750)
ggplot(data = var_imp_LR, 
       aes(x = Overall, y = reorder(names,Overall))) +
  geom_bar(stat="identity", aes(fill=category), width = 0.8) + 
  scale_fill_manual(values=cbPalette,
    name="Category",
    breaks=c("Time", "Previous glucose measurements", "Basic glucose calculations", "Glycemic variability indices",
             "Baseline subject characteristics", "Medication","Activity"),
    labels = c("Time", "Previous glucose measurements", "Basic glucose calculations", "Glycemic variability indices",
               "Baseline subject characteristics", "Medication","Activity"))+
  theme_bw(base_size = 16)+
  scale_x_continuous(expand = c(0,0))+
  labs(x = 'Feature Importance', y = 'Features', title = 'Linear Regression')+
  theme(legend.position="bottom")
dev.off()

# Save metrics for analysis
write.csv(metrics_629_LR, "metrics_629_LR.csv")
write.csv(metrics_708_LR,"metrics_708_LR.csv")
write.csv(metrics_688_LR, "metrics_688_LR.csv")
write.csv(metrics_pop_LR, "metrics_pop_LR.csv")
