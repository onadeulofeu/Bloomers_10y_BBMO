#packages----
library(caret) #random forest
library(forcats) #reorder in ggplot2
library(pdp) #partial dependence analysis
library(tidyverse)
library(stringr)
library(easystats) # check model assumptions
library(ggpmisc) # add pvalue and r2 to ggplots 

##labels ----
## prepare environmental variable labels for this analysis
labs_env_models <- c( "temperature_no_nas" = 'Temperature (ºC)',
                      "day_length_no_nas" = 'Day length (h)',
                      "bacteria_joint_sm" = 'Total bacterial abundance (cells/mL)',
                      "synechococcus_sm" = 'Synechococcus abundance (cells/mL)',
                      "chla_total_sm"  = 'Chl-a (µg/L)' ,
                      "PO4_sm" =   '[PO43-] (µM)',
                      "NH4_sm"      =   '[NH4+] (µM)',
                      "NO2_sm"  =    '[NO2-] (µM)',       
                      "NO3_sm"  =   '[NO3-] (µM)',
                      "Si_sm"   =   '[SiO4-] (µM)',
                      "PNF_sm"   =    'Phototrophic nanoflagellate abundaance (cells mL ^-1)',
                      "cryptomonas_sm" = 'Cryptomonas cells/mL',
                      "micromonas_sm"  = 'Micromonas cells/mL',
                      "HNF_sm"      = 'Heterotrophic nanoflagellate abundaance (cells mL ^-1)',
                      "ev02_sm"   = 'FL community evenness',
                      "ev3_sm"    = 'PA community evenness',        
                      "fitted_rclr" = 'previous rCLR')

## we need to create a wide dataset with environmental variables and each ASV that we consider it is a potential bloomer
## do we need to normalize environmental variables to z-scores? No need
## Is there a way to add the community effect on a blooming event in this test? add bray curtis as an env variable?
## Could we add eukaryotes data here?

## Caret package https://topepo.github.io/caret/index.html

## when performing a RF the two parameters that most likey have the biggest effect on the final accuracy are:
### mtry: numbr of variables randomly sampled as candidates at each split
### ntree: number of trees to grow
set.seed(100)
##metadata to model with interpolated missing variables----
env_data_interpolated_values_all <- read.csv2('data/env_data/env_data_interpolated_values_all.csv') |>
  rename(sample_id_num = X)

env_data_interpolated_values_all |>
  dim()

env_data_interpolated_values_all_red <- env_data_interpolated_values_all |>
  dplyr::select(- decimal_date,  -"PNF2_5um_Micro_no_nas" ,-"PNF_5um_Micro_no_nas",  -"HNF2_5um_Micro_no_nas" ,-"HNF_5um_Micro_no_nas" )

env_data_interpolated_values_all_red  |>
  colnames() 

env_data_interpolated_values_all_red |>
  dim()

## prepare the dataset----
### for the PA bloomers I need the interpolated abundance for those 3 missing samples in this fraction----
bloo_3_inter  <- read.csv('data/wavelet_3_df_deconstand.csv') |>
  as_tibble() |>
  dplyr::select(-X) |>
  dplyr::mutate(fraction = '3') |>
  dplyr::select(abundance_value = rclr, decimal_date, fraction, asv_num_f =asv_num)

## fl fraction 
bloo_02_model <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::select(decimal_date, abundance_value, fraction, asv_num_f)

## all together 
bloo_rclr_model_tb <- bloo_3_inter |>
  bind_rows(bloo_02_model)

### I pick representatives of different clusters and see how it goes

##EXAMPLES FOR EACH TYPE OF BLOOM 
### 'asv17' (PA), 'asv23' (PA), 'asv11' (FL), 'asv62' (FL), 'asv72' (FL), 'asv555' (FL))

## ideas data that I should add here: diatoms, ciliates, dinophlgellates, 

## do we add diversity in the random forest?

## I prepare the data for all my bloomers and then I will filter by the one I'm interested in----
date_tb <- asv_tab_all_bloo_z_tax |>
  dplyr::select(decimal_date, date) |>
  distinct(decimal_date, date)

diff_time_tb_env <- bloo_rclr_model_tb |>
  arrange(decimal_date) |>
  left_join(date_tb) |>
  arrange(decimal_date) |>
  group_by(asv_num_f, fraction) |>
  dplyr::mutate(diff_rclr = c(NA, diff(abundance_value)),
         diff_time = c(NA, diff(date))) |>
  dplyr::mutate(diff_rclr_time = diff_rclr/diff_time) |>
  dplyr::mutate(sample_id_num = row_number()) |>
  left_join(env_data_interpolated_values_all_red) ## add envdata (not smoothed) 

# One by one example
# asv62_rclr <- asv_tab_all_bloo_z_tax |>
#   dplyr::filter(abundance_type == 'rclr') |>
#   dplyr::filter(asv_num == 'asv62') |>
#   dplyr::select(asv_num, abundance_value, date) |>
#   arrange(date)

##difference in abundance from one day to the next one
# diff_abund_tb <- diff(asv62_rclr$abundance_value) |>
#     as_tibble_col(column_name = 'difference_rclr') |>
#     dplyr::mutate(difference_rclr_num = as.numeric(difference_rclr))
#  
# ## calculate from one point to another the increase or decrease
# ## divide it by the difference in time from one point to the next
# diff_day_tb <- asv_tab_all_bloo_z_tax |>
#   dplyr::select(date) |>
#   dplyr::arrange(date) |>
#   distinct(date) %$%
#   diff(date) |>
#   as_tibble_col(column_name = 'difference_day') |>
#   dplyr::mutate(difference_day_num = as.numeric(difference_day))
# 
# ##calculate the difference in time
# diff_time_tb <- diff_abund_tb |>
#   bind_cols(diff_day_tb) |>
#   dplyr::mutate(diff_time = difference_rclr/difference_day_num) |>
#   dplyr::select(diff_time) |>
#   dplyr::mutate(sample_id_num = row_number())

## add the metadata to my target bloomer

## to predict the change in abundance we need the environmental data from my previous datapoint. difference 2-1 gets the env_data in 1
# diff_time_tb_env <- diff_time_tb  |>
#   left_join(env_data_interpolated_values_all_red)

# MODELING WITHOUT SMOOTHING MY VARIABLES ----

diff_time_tb_env |>
  colnames()

## ASV62----
model_tb <- diff_time_tb_env |>
  #dplyr::left_join(community_evenness_all) |> #calculated in the next chunk
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv62' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -decimal_date, -diff_rclr, -diff_time, -date, -fraction, -abundance_value, -BP_FC1.55_no_nas)

##add previous CLR value as explanatory variable
asvs_prev_abund <- diff_time_tb_env |>
  dplyr::filter(asv_num_f == 'asv62' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  ungroup() |>
  dplyr::select(abundance_value)

asvs_prev_abund <- asvs_prev_abund[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_prev_abund)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset----
bloo_final_asv62 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=300, importance=T)
bloo_final_asv62

## once we have trained our model we should not retrain it again (it will lead to different results.)
importance_df_asv62 <- as.data.frame(bloo_final_asv62$finalModel$importance)
importance_df_asv62 <- importance_df_asv62[order(-importance_df_asv62$`%IncMSE`),]
importance_df_asv62

## Linear model
linear_fit <- lm(diff_rclr_time ~.+temperature_no_nas*day_length_no_nas*chla_total_no_nas*PO4_no_nas*NH4_no_nas, data=model_tb)
summary(linear_fit)
anova(linear_fit, type=2)

## asv23----
model_tb <- diff_time_tb_env |>
 # dplyr::left_join(community_evenness_all) |> #calculated in the next chunk
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv23' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -decimal_date, -diff_rclr, -diff_time, -date, -fraction, -abundance_value, -BP_FC1.55_no_nas)

##add previous CLR value as explanatory variable
asvs_prev_abund <- diff_time_tb_env |>
  dplyr::filter(asv_num_f == 'asv23' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  ungroup() |>
  dplyr::select(abundance_value)

asvs_prev_abund <- asvs_prev_abund[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_prev_abund)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset----
bloo_final_asv23 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv23

## once we have trained our model we should not retrain it again (it will lead to different results.)
importance_df_asv23 <- as.data.frame(bloo_final_asv23$finalModel$importance)
importance_df_asv23 <- importance_df_asv23[order(-importance_df_asv23$`%IncMSE`),]
importance_df_asv23

## Linear model
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## asv72----
model_tb <- diff_time_tb_env |>
   dplyr::left_join(community_evenness_all) |> #calculated in the next chunk
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv72' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -decimal_date, -diff_rclr, -diff_time, -date, -fraction, -abundance_value, -BP_FC1.55_no_nas)

##add previous CLR value as explanatory variable
asvs_prev_abund <- diff_time_tb_env |>
  dplyr::filter(asv_num_f == 'asv72' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  ungroup() |>
  dplyr::select(abundance_value)

asvs_prev_abund <- asvs_prev_abund[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_prev_abund)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset----
bloo_final_asv72 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv72

## once we have trained our model we should not retrain it again (it will lead to different results.)
importance_df_asv72 <- as.data.frame(bloo_final_asv72$finalModel$importance)
importance_df_asv72 <- importance_df_asv72[order(-importance_df_asv72$`%IncMSE`),]
importance_df_asv72

## Linear model
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## asv17----
model_tb <- diff_time_tb_env |>
  # dplyr::left_join(community_evenness_all) |> #calculated in the next chunk
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv17' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -decimal_date, -diff_rclr, -diff_time, -date, -fraction, -abundance_value, -BP_FC1.55_no_nas)

##add previous CLR value as explanatory variable
asvs_prev_abund <- diff_time_tb_env |>
  dplyr::filter(asv_num_f == 'asv17' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  ungroup() |>
  dplyr::select(abundance_value)

asvs_prev_abund <- asvs_prev_abund[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_prev_abund)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset----
bloo_final_asv17 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv17

## once we have trained our model we should not retrain it again (it will lead to different results.)
importance_df_asv17 <- as.data.frame(bloo_final_asv17$finalModel$importance)
importance_df_asv17 <- importance_df_asv17[order(-importance_df_asv17$`%IncMSE`),]
importance_df_asv17

## Linear model
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## asv11----
model_tb <- diff_time_tb_env |>
  # dplyr::left_join(community_evenness_all) |> #calculated in the next chunk
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv11' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -decimal_date, -diff_rclr, -diff_time, -date, -fraction, -abundance_value, -BP_FC1.55_no_nas)

##add previous CLR value as explanatory variable
asvs_prev_abund <- diff_time_tb_env |>
  dplyr::filter(asv_num_f == 'asv11' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  ungroup() |>
  dplyr::select(abundance_value)

asvs_prev_abund <- asvs_prev_abund[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_prev_abund)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset----
bloo_final_asv11 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv11

## once we have trained our model we should not retrain it again (it will lead to different results.)
importance_df_asv11 <- as.data.frame(bloo_final_asv11$finalModel$importance)
importance_df_asv11 <- importance_df_asv11[order(-importance_df_asv11$`%IncMSE`),]
importance_df_asv11

## Linear model
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## asv555----
model_tb <- diff_time_tb_env |>
  # dplyr::left_join(community_evenness_all) |> #calculated in the next chunk
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv555' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -decimal_date, -diff_rclr, -diff_time, -date, -fraction, -abundance_value, -BP_FC1.55_no_nas)

##add previous CLR value as explanatory variable
asvs_prev_abund <- diff_time_tb_env |>
  dplyr::filter(asv_num_f == 'asv555' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  ungroup() |>
  dplyr::select(abundance_value)

asvs_prev_abund <- asvs_prev_abund[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_prev_abund)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset----
bloo_final_asv555 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv555

## once we have trained our model we should not retrain it again (it will lead to different results.)
importance_df_asv555 <- as.data.frame(bloo_final_asv555$finalModel$importance)
importance_df_asv555 <- importance_df_asv555[order(-importance_df_asv555$`%IncMSE`),]
importance_df_asv555

## Linear model
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

##plot importance of variables ----
##This importance is a measure of by how much removing a variable decreases accuracy, and vice versa — by how much including a variable increases accuracy. 
## The %IncMSE (Percentage Increase in Mean Squared Error) is a metric used to measure the 
## importance of each predictor (variable) in a random forest model. It indicates how much the mean squared error (MSE) increases when a particular predictor is randomly permuted while keeping all other predictors unchanged.

##dataset from my example ASVs with it's 
im62 <- importance_df_asv62 |>
  dplyr::mutate(asv_num = 'asv62') |>
  rownames_to_column(var = 'variable') |>
  arrange(desc(`%IncMSE`))

im23 <- importance_df_asv23 |>
  dplyr::mutate(asv_num = 'asv23') |>
  rownames_to_column(var = 'variable')

im72 <- importance_df_asv72 |>
  dplyr::mutate(asv_num = 'asv72') |>
  rownames_to_column(var = 'variable')

im17 <- importance_df_asv17 |>
  dplyr::mutate(asv_num = 'asv17') |>
  rownames_to_column(var = 'variable')

im11 <- importance_df_asv11 |>
  dplyr::mutate(asv_num = 'asv11') |>
  rownames_to_column(var = 'variable')

im555 <- importance_df_asv555 |>
  dplyr::mutate(asv_num = 'asv555') |>
  rownames_to_column(var = 'variable')

im_all <- im62 |>
  bind_rows(im23) |>
  bind_rows(im72) |>
  bind_rows(im17) |>
  bind_rows(im11) |>
  bind_rows(im555)|>
  as_tibble() 

im_all |>
  colnames() <- c('variable', 'incMSE', 'IncNodePurity', 'asv_num')

im_all |>
  str()

im_all <- im_all |>
  dplyr::mutate(asv_num = as.factor(asv_num))

im_all$asv_num <- im_all$asv_num |>
  factor(levels = c('asv23', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

importance_example_bloo <- im_all |>
  ggplot(aes(fct_infreq(variable), incMSE))+
  geom_col()+
  facet_wrap(~ asv_num)+
  scale_x_discrete(labels = labs_env_models)+
  labs(y = '%IncMSE')+
  theme_bw()+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(),
        text = element_text(size = 6))

importance_example_bloo
# 
# ggsave(importance_example_bloo, filename = 'importance_example_bloo_no-sm.pdf',
#        path = 'Results/Figures/random_forest/',
#        width = 188, height = 150, units = 'mm')

# To increase the performance of our model we could smooth our environmental variables----
diff_time_tb_env |>
  colnames()

decimal_date_tb <- diff_time_tb_env |>
  ungroup() |>
  dplyr::filter(fraction == '3') |>
  dplyr::select(-abundance_value, -fraction, -asv_num_f, -diff_rclr, -diff_time, -diff_rclr_time, -sample_id_num) |>
  distinct(decimal_date, bacteria_joint, synechococcus, temperature_no_nas, day_length_no_nas, chla_total_no_nas,  PO4_no_nas,
           NH4_no_nas, NO2_no_nas, NO3_no_nas,        
           Si_no_nas, BP_FC1.55_no_nas, PNF_Micro_no_nas, cryptomonas_no_nas, micromonas_no_nas, HNF_Micro_no_nas) |>
  dplyr::select(decimal_date) |>
  arrange(decimal_date)

## dataset with just env variables to smooth them
diff_rclr_time_tb_env <- diff_time_tb_env |>
  ungroup() |>
  dplyr::select(-abundance_value, -fraction, -asv_num_f, -diff_rclr, -diff_time, -diff_rclr_time, -sample_id_num) |>
  distinct(decimal_date, bacteria_joint, synechococcus, temperature_no_nas, day_length_no_nas, chla_total_no_nas,  PO4_no_nas,
           NH4_no_nas, NO2_no_nas, NO3_no_nas,        
           Si_no_nas, BP_FC1.55_no_nas, PNF_Micro_no_nas, cryptomonas_no_nas, micromonas_no_nas, HNF_Micro_no_nas)

## I would like to add evenness as explanatory variable-------
asv_tab_bbmo_10y_w_rar <- read.csv2('data/asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()|> 
  rename('sample_id' = 'X') 

## Rarefied dataset to calculate Community Evenness----
source('~/Documentos/Doctorat/Bloomers/R/community_evenness.R')

community_eveness_02 <- asv_tab_bbmo_10y_w_rar |>
  #as.data.frame() |>
  #tibble::rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar_02 = community_evenness(abundances = reads_rar, index = 'Pielou'))

##for the PA dataset I need the interpolated values and then rarefy
#### #write.csv(asv_tab_bbmo_10y_w_3_inter, 'data/asv_tab_bbmo_10y_w_3_inter_reads.csv')

asv_tab_bbmo_10y_w_3_inter ## this dataframe is from seasonality script where I interpolated reads for the 3 missing samples

#### calculate the minimum read size of all samples to rarefy at that number
min_n_seqs <- asv_tab_bbmo_10y_w_3_inter |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads') |>
  group_by(sample_id) |>
  dplyr::summarize(n_seqs = sum(reads)) |>
  dplyr::summarize(min = min(n_seqs)) |>
  pull(min)

# we use rarefaction (which repeats the subsampling step many times)
# perform this a 1000 times to get an empirical diversity values calculating the mean value for each timepoint.
# asv_tab_bbmo_10y_w_rar_3_inter <- rrarefy.perm(round(asv_tab_bbmo_10y_w_3_inter),
#                                        sample = min_n_seqs,
#                                        n = 1000,
#                                        round.out = T)

#write.csv(asv_tab_bbmo_10y_w_rar_3_inter, file = 'asv_tab_bbmo_10y_w_rar_3_inter.csv')

community_eveness_3 <- asv_tab_bbmo_10y_w_rar_3_inter |>
  as.data.frame() |>
  rownames_to_column(var = 'date') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  #dplyr::filter(str_detect(sample_id, '_3_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(date) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar_3 = community_evenness(abundances = reads_rar, index = 'Pielou'))

## I bring together evenness in 3 and 02 ----
### add decimal date to the evenness datasets
decimal_date_sample_id_tb <- diff_time_tb_env |>
  ungroup() |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::select(-abundance_value, -fraction, -asv_num_f, -diff_rclr, -diff_time, -diff_rclr_time) |>
  distinct(decimal_date, bacteria_joint, synechococcus, temperature_no_nas, day_length_no_nas, chla_total_no_nas,  PO4_no_nas,
           NH4_no_nas, NO2_no_nas, NO3_no_nas,  date,      
           Si_no_nas, BP_FC1.55_no_nas, PNF_Micro_no_nas, cryptomonas_no_nas, micromonas_no_nas, HNF_Micro_no_nas, sample_id_num) |>
  dplyr::select(decimal_date, date, sample_id_num) |>
  arrange(decimal_date)

community_evenness_02_date <- community_eveness_02 |>
  arrange(sample_id) |>
  bind_cols(decimal_date_sample_id_tb)

community_evenness_3_date <- community_eveness_3 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(decimal_date_sample_id_tb, by = 'date')

community_evenness_all <- community_evenness_02_date |>
  left_join(community_evenness_3_date) |>
  dplyr::select(-sample_id, -date, -sample_id_num)

##add community evenness as explanatory variable-----
diff_rclr_time_tb_env <- diff_rclr_time_tb_env |>
  left_join(community_evenness_all)

## day_length_no_nas----
#fit the loess model
fit <- loess(day_length_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")
## synechococcus----
#fit the loess model
fit <- loess( synechococcus ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

synechococcus_sm <- predicted_values |>
  as_tibble_col(column_name = 'synechococcus_sm')

## bacteria_joint----
#fit the loess model
fit <- loess( bacteria_joint ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

bacteria_joint_sm <- predicted_values |>
  as_tibble_col(column_name = 'bacteria_joint_sm')

## temperature_no_nas----
#fit the loess model
# fit <- loess(temperature_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)
# 
# # predict the fited values
# predicted_values <- predict(fit)
# 
# # Obtain the residuals from the loess model
# residuals <- residuals(fit)
# 
# ## check the residuals
# ## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# # Compute the autocorrelation function (ACF) of the residuals
# acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

## chla_total_no_nas----
#fit the loess model
fit <- loess(chla_total_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

chla_total_sm <- predicted_values |>
  as_tibble_col(column_name = 'chla_total_sm')

## PO4_no_nas----
#fit the loess model
fit <- loess(PO4_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residual

PO4_sm <- predicted_values |>
  as_tibble_col(column_name = 'PO4_sm')

## NH4_no_nas----
#fit the loess model
fit <- loess(NH4_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

NH4_sm <- predicted_values |>
  as_tibble_col(column_name = 'NH4_sm')
 
## NO2_no_nas----
#fit the loess model
fit <- loess(NO2_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

NO2_sm <- predicted_values |>
  as_tibble_col(column_name = 'NO2_sm')

## NO3_no_nas----
#fit the loess model
fit <- loess(NO3_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

NO3_sm <- predicted_values |>
  as_tibble_col(column_name = 'NO3_sm')

## Si_no_nas----
#fit the loess model
fit <- loess(Si_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

Si_sm <- predicted_values |>
  as_tibble_col(column_name = 'Si_sm')

## BP_FC1.55_no_nas----
#fit the loess model
fit <- loess(BP_FC1.55_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

BP_FC1.55_sm <- predicted_values |>
  as_tibble_col(column_name = 'BP_FC1.55_sm')

## PNF_no_nas----
#fit the loess model
fit <- loess(PNF_Micro_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

PNF_sm <- predicted_values |>
  as_tibble_col(column_name = 'PNF_sm')

## cryptomonas_no_nas----
#fit the loess model
fit <- loess(cryptomonas_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

cryptomonas_sm <- predicted_values |>
  as_tibble_col(column_name = 'cryptomonas_sm')

## micromonas_no_nas----
#fit the loess model
fit <- loess(micromonas_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

micromonas_sm <- predicted_values |>
  as_tibble_col(column_name = 'micromonas_sm')

## HNF_Micro_no_nas----
#fit the loess model
fit <- loess(HNF_Micro_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

HNF_sm <- predicted_values |>
  as_tibble_col(column_name = 'HNF_sm')

## Community evenness 02----
#fit the loess model
fit <- loess( community_eveness_rar_02 ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

ev02_sm <- predicted_values |>
  as_tibble_col(column_name = 'ev02_sm')

## Community evenness 3----
#fit the loess model
fit <- loess( community_eveness_rar_3 ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals") # Compute the autocorrelation function (ACF) of the residuals

ev3_sm <- predicted_values |>
  as_tibble_col(column_name = 'ev3_sm')

## variables I do not smooth (too much variability in the residuals)----
#### day_length_no_nas, temperature_no_nas
nsm_emv_data <- diff_rclr_time_tb_env |>
  dplyr::select(day_length_no_nas, temperature_no_nas)
## variables I do not smooth (too much variability in the residuals)
#### day_length_no_nas, temperature_no_nas

nsm_emv_data <- diff_rclr_time_tb_env |>
  dplyr::select(day_length_no_nas, temperature_no_nas)

## variables I smooth (they were noisy and the smoothing work well to remove all those noise)----
## bacteria_joint, synechococcus, chla_total, PO4, NH4, NO2, NO3, Si, BP_FC1.55, PNF, cryptomonas, micromonas, HNF
sm_env_data <- bacteria_joint_sm |>
  bind_cols(synechococcus_sm) |>
  bind_cols(chla_total_sm) |>
  bind_cols(PO4_sm) |>
  bind_cols(NH4_sm) |>
  bind_cols(NO2_sm ) |>
  bind_cols(NO3_sm) |>
  bind_cols(Si_sm) |>
  bind_cols(BP_FC1.55_sm) |>
  bind_cols(PNF_sm) |>
  bind_cols(cryptomonas_sm) |>
  bind_cols(micromonas_sm) |>
  bind_cols(HNF_sm) |>
  bind_cols(ev02_sm) |>
  bind_cols(ev3_sm)

##all new environental data
env_data_new <- nsm_emv_data |>
  bind_cols(sm_env_data) |>
  bind_cols(decimal_date_tb)

## observe the new smooth data (did we removed some noise?)
### decimal_date_tibble <- decimal_date_tb[-120,]

env_data_new_l <- env_data_new |>
 # bind_cols(decimal_date_tb) |>
  dplyr::select( -temperature_no_nas, -day_length_no_nas) |>
  pivot_longer(cols = c(!decimal_date), values_to = 'value', names_to = 'env_variable') |>
  separate(env_variable, into = c('env_variable', 'type', 'another'), sep = '_') |>
  dplyr::mutate(type_value = 'smooth')
  
diff_rclr_time_tb_env_l <- diff_rclr_time_tb_env |>
  dplyr::select( -temperature_no_nas, -day_length_no_nas) |>
  pivot_longer(cols = c(!decimal_date), values_to = 'value', names_to = 'env_variable') |>
  separate(env_variable, into = c('env_variable', 'type', 'another'), sep = '_') |>
  dplyr::mutate(type_value = 'original')

community_evenness_all_l <- community_evenness_all |>
  dplyr::select(ev02 = community_eveness_rar_02 , decimal_date, ev3 = community_eveness_rar_3) |>
  pivot_longer(cols = !decimal_date, names_to = 'env_variable', values_to = 'value')|>
  dplyr::mutate(type_value = 'original')

env_data_new_l |>
   bind_rows(diff_rclr_time_tb_env_l) |>
  dplyr::filter(env_variable != 'community') |>
  bind_rows(community_evenness_all_l) |>
  #dplyr::filter(env_variable == 'community') |>
  ggplot(aes(decimal_date, value))+
  geom_point(aes(shape = type_value), alpha = 0.4)+
  geom_line(aes(color = type_value))+
  facet_wrap(vars(env_variable), scales = 'free')+
  theme_bw()

#prepare the env data
env_data_new_red <- env_data_new |>
  dplyr::select(-decimal_date)

env_data_new_red <- env_data_new_red[-120,]

# Model with some smoothed env variables -----
## BLOOMER ASV62----
model_tb <- diff_time_tb_env |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv62' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -diff_rclr, -diff_time, -date, -fraction, -abundance_value) |>
  dplyr::select(-contains('_no_nas'), -bacteria_joint, -synechococcus) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(decimal_date)

model_tb |>
  ggplot(aes(decimal_date, diff_rclr_time))+
  geom_line()
  
### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry
r2_summary_asv62 <- bloo_rf$results |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv62')

r2_summary_asv62 |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()

y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## remove some variables----
model_tb |>
  colnames()

model_tb <- model_tb |>
  dplyr::select( -bacteria_joint_sm, -BP_FC1.55_sm, -chla_total_sm, -PNF_sm,  -cryptomonas_sm,  -micromonas_sm, -HNF_sm, - synechococcus_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(1000)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=200, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final

## Linear model---
# linear_fit <- lm(diff_rclr_time ~., data=model_tb)
# summary(linear_fit)

## BLOOMER ASV11-----
## Model with some smoothed env variables---
model_tb <- diff_time_tb_env |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv11' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -diff_rclr, -diff_time, -date, -fraction, -abundance_value) |>
  dplyr::select(-contains('_no_nas'), -bacteria_joint, -synechococcus) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## remove some variables
model_tb <- model_tb |>
  dplyr::select( -bacteria_joint_sm, -BP_FC1.55_sm, -chla_total_sm, -PNF_sm,  -cryptomonas_sm,  -micromonas_sm, -HNF_sm, - synechococcus_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=400, importance=T)
bloo_final

## Linear model---
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## BLOOMER ASV23-----
model_tb <- diff_time_tb_env |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv23' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -diff_rclr, -diff_time, -date, -fraction, -abundance_value) |>
  dplyr::select(-contains('_no_nas'), -bacteria_joint, -synechococcus) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## remove some variables
model_tb <- model_tb |>
  dplyr::select( -bacteria_joint_sm, -BP_FC1.55_sm, -chla_total_sm, -PNF_sm,  -cryptomonas_sm,  -micromonas_sm, -HNF_sm, - synechococcus_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=400, importance=T)
bloo_final

## Linear model---
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## BLOOMER asv72-----
model_tb <- diff_time_tb_env |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv72' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -diff_rclr, -diff_time, -date, -fraction, -abundance_value) |>
  dplyr::select(-contains('_no_nas'), -bacteria_joint, -synechococcus) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## remove some variables
model_tb <- model_tb |>
  dplyr::select( -bacteria_joint_sm, -BP_FC1.55_sm, -chla_total_sm, -PNF_sm,  -cryptomonas_sm,  -micromonas_sm, -HNF_sm, - synechococcus_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=400, importance=T)
bloo_final

## Linear model---
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## BLOOMER ASV17-----
model_tb <- diff_time_tb_env |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv17' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -diff_rclr, -diff_time, -date, -fraction, -abundance_value) |>
  dplyr::select(-contains('_no_nas'), -bacteria_joint, -synechococcus) |>
  bind_cols(env_data_new_red) |> ## smooth env variables
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## remove some variables
model_tb <- model_tb |>
  dplyr::select( -bacteria_joint_sm, -BP_FC1.55_sm, -chla_total_sm, -PNF_sm,  -cryptomonas_sm,  -micromonas_sm, -HNF_sm, - synechococcus_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=400, importance=T)
bloo_final

## Linear model---
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## BLOOMER ASV555-----
model_tb <- diff_time_tb_env |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv555' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -diff_rclr, -diff_time, -date, -fraction, -abundance_value) |>
  dplyr::select(-contains('_no_nas'), -bacteria_joint, -synechococcus) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## remove some variables
model_tb <- model_tb |>
  #dplyr::select( -bacteria_joint_sm, -Si_sm, -BP_FC1.55_sm) 
  dplyr::select( -bacteria_joint_sm, -BP_FC1.55_sm, -chla_total_sm, -PNF_sm,  -cryptomonas_sm,  -micromonas_sm, -HNF_sm, - synechococcus_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
set.seed(100)
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=400, importance=T)
bloo_final

## Linear model---
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)




#### SMOOTH MY RESPONSE VARIABLE TOO-----
## dataset with just clr variables to smooth them
rclr_time_tb_3 <- diff_time_tb_env |> 
  ungroup() |>
  dplyr::select(abundance_value, fraction, asv_num_f, decimal_date) |>
  dplyr::filter(fraction == '3') |>
  pivot_wider(names_from = 'asv_num_f', id_cols = c('decimal_date', 'fraction'), values_from = abundance_value) |>
  arrange(decimal_date)

rclr_time_tb_02 <- diff_time_tb_env |>
  ungroup() |>
  dplyr::select(abundance_value, fraction, asv_num_f, decimal_date) |>
  dplyr::filter(fraction == '0.2') |>
  pivot_wider(names_from = 'asv_num_f', id_cols = c('decimal_date', 'fraction'), values_from = abundance_value) |>
  arrange(decimal_date)

##mantain the fraction column
fraction_date_02 <- rclr_time_tb_02 |>
  dplyr::select(fraction, decimal_date)

##mantain the fraction column
fraction_date_3 <- rclr_time_tb_3 |>
  dplyr::select(fraction, decimal_date)

## first i try for my examples and they i will decide if I loop over all the variables-----
#### 3 (asv23, as72, asv17)----

#fit the loess model
fit <- loess(asv23 ~ decimal_date, data = rclr_time_tb_3, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

asv23_sm <- predicted_values |>
  as_tibble_col(column_name = 'asv23') |>
  bind_cols(decimal_date_tb) |>
  dplyr::select(decimal_date, 'fitted_rclr' = 'asv23') |>
  dplyr::mutate(asv_num = 'asv23') |>
  dplyr::left_join(fraction_date_3)

rclr_time_tb_3 |>
  dplyr::select(decimal_date, 'rclr' = 'asv23') |>
  left_join(asv23_sm) |>
  pivot_longer(cols = c('rclr', 'fitted_rclr'), names_to = 'approx', values_to = 'abundance_value') |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line(aes(group = approx, color = approx))+
  scale_color_manual(values = c('black', 'grey'))+
  theme_bw()+
  theme(panel.grid = element_blank())

#fit the loess model
fit <- loess(asv72 ~ decimal_date, data = rclr_time_tb_3, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

asv72_sm <- predicted_values |>
  as_tibble_col(column_name = 'asv72') |>
  bind_cols(decimal_date_tb) |>
  dplyr::select(decimal_date, 'fitted_rclr' = 'asv72') |>
  dplyr::mutate(asv_num = 'asv72') |>
  dplyr::left_join(fraction_date_3)

rclr_time_tb_3 |>
  dplyr::select(decimal_date, 'rclr' = 'asv72') |>
  left_join(asv72_sm) |>
  pivot_longer(cols = c('rclr', 'fitted_rclr'), names_to = 'approx', values_to = 'abundance_value') |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line(aes(group = approx, color = approx))+
  scale_color_manual(values = c('black', 'grey'))+
  theme_bw()+
  theme(panel.grid = element_blank())

#fit the loess model
fit <- loess(asv17 ~ decimal_date, data = rclr_time_tb_3, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

asv17_sm <- predicted_values |>
  as_tibble_col(column_name = 'asv17') |>
  bind_cols(decimal_date_tb) |>
  dplyr::select(decimal_date, 'fitted_rclr' = 'asv17') |>
  dplyr::mutate(asv_num = 'asv17') |>
  dplyr::left_join(fraction_date_3)

rclr_time_tb_3 |>
  dplyr::select(decimal_date, 'rclr' = 'asv17') |>
  left_join(asv17_sm) |>
  pivot_longer(cols = c('rclr', 'fitted_rclr'), names_to = 'approx', values_to = 'abundance_value') |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line(aes(group = approx, color = approx))+
  scale_color_manual(values = c('black', 'grey'))+
  theme_bw()+
  theme(panel.grid = element_blank())

### 02 (asv62, asv11, asv555)----

#fit the loess model
fit <- loess(asv62 ~ decimal_date, data = rclr_time_tb_02, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

asv62_sm <- predicted_values |>
  as_tibble_col(column_name = 'asv62') |>
  bind_cols(decimal_date_tb) |>
  dplyr::select(decimal_date, 'fitted_rclr' = 'asv62') |>
  dplyr::mutate(asv_num = 'asv62') |>
  dplyr::left_join(fraction_date_02)

rclr_time_tb_02 |>
  dplyr::select(decimal_date, 'rclr' = 'asv62') |>
  left_join(asv62_sm) |>
  pivot_longer(cols = c('rclr', 'fitted_rclr'), names_to = 'approx', values_to = 'abundance_value') |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line(aes(group = approx, color = approx))+
  scale_color_manual(values = c('black', 'grey'))+
  theme_bw()+
  theme(panel.grid = element_blank())

#fit the loess model
fit <- loess(asv11 ~ decimal_date, data = rclr_time_tb_02, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

asv11_sm <- predicted_values |>
  as_tibble_col(column_name = 'asv11') |>
  bind_cols(decimal_date_tb) |>
  dplyr::select(decimal_date, 'fitted_rclr' = 'asv11') |>
  dplyr::mutate(asv_num = 'asv11') |>
  dplyr::left_join(fraction_date_02)

rclr_time_tb_02 |>
  dplyr::select(decimal_date, 'rclr' = 'asv11') |>
  left_join(asv11_sm) |>
  pivot_longer(cols = c('rclr', 'fitted_rclr'), names_to = 'approx', values_to = 'abundance_value') |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line(aes(group = approx, color = approx))+
  scale_color_manual(values = c('black', 'grey'))+
  theme_bw()+
  theme(panel.grid = element_blank())

#fit the loess model
fit <- loess(asv555 ~ decimal_date, data = rclr_time_tb_02, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")


asv555_sm <- predicted_values |>
  as_tibble_col(column_name = 'asv555') |>
  bind_cols(decimal_date_tb) |>
  dplyr::select(decimal_date, 'fitted_rclr' = 'asv555') |>
  dplyr::mutate(asv_num = 'asv555') |>
  dplyr::left_join(fraction_date_02)

rclr_time_tb_02 |>
  dplyr::select(decimal_date, 'rclr' = 'asv555') |>
  left_join(asv555_sm) |>
  pivot_longer(cols = c('rclr', 'fitted_rclr'), names_to = 'approx', values_to = 'abundance_value') |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line(aes(group = approx, color = approx))+
  scale_color_manual(values = c('black', 'grey'))+
  theme_bw()+
  theme(panel.grid = element_blank())

### bring together all the new smooth asvs-----
asvs_sm <- asv23_sm |>
  bind_rows(asv555_sm) |>
  bind_rows(asv62_sm) |>
  bind_rows(asv72_sm) |>
  bind_rows(asv17_sm) |>
  bind_rows(asv11_sm)

## calculate the increase in the new smoothed variables----
asvs_sm_diff <- asvs_sm |>
  arrange(decimal_date) |>
  left_join(date_tb) |>
  arrange(decimal_date) |>
  group_by(asv_num, fraction) |>
  dplyr::mutate(diff_rclr = c(NA, diff(fitted_rclr)),
                diff_time = c(NA, diff(date))) |>
  dplyr::mutate(diff_rclr_time = diff_rclr/diff_time) |>
  dplyr::mutate(sample_id_num = row_number())

##observe how it is to plot diff than abundances----
asvs_sm_diff |>
  ggplot(aes(decimal_date, diff_rclr_time))+
  geom_line()+
  facet_wrap(vars(asv_num))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = 'transparent'))

asvs_sm_diff |>
  ggplot(aes(decimal_date, fitted_rclr))+
  geom_line()+
  facet_wrap(vars(asv_num))+  
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = 'transparent'))

## blooms and rate of change from one point to the next one (just to observe) ------
x <- asv_tab_all_bloo_z_tax_examples |>
  dplyr::filter(bloom == 'Bloom' &
                  asv_num == 'asv62' &
                  fraction == '0.2') |>
  dplyr::select(decimal_date)

asvs_sm_diff |>
  dplyr::filter(asv_num == 'asv62') |>
  ggplot(aes(decimal_date, diff_rclr_time))+
  geom_line()+
  geom_vline(xintercept = x$decimal_date)

x <- asv_tab_all_bloo_z_tax_examples |>
  dplyr::filter(bloom == 'Bloom' &
                  asv_num == 'asv72' &
                  fraction == '3') |>
  dplyr::select(decimal_date)

asvs_sm_diff |>
  dplyr::filter(asv_num == 'asv72') |>
  ggplot(aes(decimal_date, diff_rclr_time))+
  geom_line()+
  geom_vline( xintercept = x$decimal_date)

x <- asv_tab_all_bloo_z_tax_examples |>
  dplyr::filter(bloom == 'Bloom' &
                  asv_num == 'asv23' &
                  fraction == '3') |>
  dplyr::select(decimal_date)

asvs_sm_diff |>
  dplyr::filter(asv_num == 'asv23') |>
  ggplot(aes(decimal_date, diff_rclr_time))+
  geom_line()+
  geom_vline(data = x, xintercept = x$decimal_date)

x <- asv_tab_all_bloo_z_tax_examples |>
  dplyr::filter(bloom == 'Bloom' &
                  asv_num == 'asv17' &
                  fraction == '3') |>
  dplyr::select(decimal_date)

asvs_sm_diff |>
  dplyr::filter(asv_num == 'asv17') |>
  ggplot(aes(decimal_date, diff_rclr_time))+
  geom_line()+
  geom_vline(data = x, xintercept = x$decimal_date)

x <- asv_tab_all_bloo_z_tax_examples |>
  dplyr::filter(bloom == 'Bloom' &
                  asv_num == 'asv11' &
                  fraction == '0.2') |>
  dplyr::select(decimal_date)

asvs_sm_diff |>
  dplyr::filter(asv_num == 'asv11') |>
  ggplot(aes(decimal_date, diff_rclr_time))+
  geom_line()+
  geom_vline(data = x, xintercept = x$decimal_date)


x <- asv_tab_all_bloo_z_tax_examples |>
  dplyr::filter(bloom == 'Bloom' &
                  asv_num == 'asv555' &
                  fraction == '0.2') |>
  dplyr::select(decimal_date)

asvs_sm_diff |>
  dplyr::filter(asv_num == 'asv555') |>
  ggplot(aes(decimal_date, diff_rclr_time))+
  geom_line()+
  geom_vline(data = x, xintercept = x$decimal_date)




# Model the smooth rCLR and env variables----
### BLOOMER ASV62 (smooth) ----
#env_data_new_red <- env_data_new_red[-120,]
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv62' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv62' &
                               fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv) |>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) ##in this case is ntree 500

# # View summary of the trained model
# summary(custom)
# 
# # Plot the model
# plot(custom)


## retrain model on entire dataset----
bloo_final_asv62 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=300, importance=T)
bloo_final_asv62
## once we have trained our model we should not retrain it again (it will lead to different results.)
importance_df_asv62 <- as.data.frame(bloo_final_asv62$finalModel$importance)
importance_df_asv62 <- importance_df_asv62[order(-importance_df_asv62$`%IncMSE`),]
importance_df_asv62

## remove some variables----
model_tb |>
  colnames()

model_tb <- model_tb |>
  dplyr::select(-BP_FC1.55_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=200, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final

### BLOOMER asv23----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv23' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv23' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) ##100 is best

## retrain model on entire dataset----
bloo_final_asv23 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=100, importance=T)
bloo_final_asv23

importance_df_asv23 <- as.data.frame(bloo_final_asv23$finalModel$importance)
importance_df_asv23 <- importance_df_asv23[order(-importance_df_asv23$`%IncMSE`),]
importance_df_asv23

## remove some variables----
model_tb <- model_tb |>
  dplyr::select( -BP_FC1.55_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=200, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final

### BLOOMER asv72----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv72' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv72' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) ## ntree = 100

## retrain model on entire dataset----
bloo_final_asv72 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=300, importance=T)
bloo_final_asv72

importance_df_asv72 <- as.data.frame(bloo_final_asv72$finalModel$importance)
importance_df_asv72 <- importance_df_asv72[order(-importance_df_asv72$`%IncMSE`),]
importance_df_asv72

## remove some variables----
model_tb |>
  colnames()

model_tb <- model_tb |>
  dplyr::select(  -BP_FC1.55_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=200, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final

### BLOOMER asv17----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv17' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv17' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) ##ntree = 300

## retrain model on entire dataset----
bloo_final_asv17 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=300, importance=T)
bloo_final_asv17

importance_df_asv17 <- as.data.frame(bloo_final_asv17$finalModel$importance)
importance_df_asv17 <- importance_df_asv72[order(-importance_df_asv17$`%IncMSE`),]
importance_df_asv17

## remove some variables----
model_tb |>
  colnames()

model_tb <- model_tb |>
  dplyr::select( -BP_FC1.55_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=200, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final

### BLOOMER asv11----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) ## ntree = 300

## retrain model on entire dataset----
bloo_final_asv11 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=300, importance=T)
bloo_final_asv11

importance_df_asv11 <- as.data.frame(bloo_final_asv11$finalModel$importance)
importance_df_asv11 <- importance_df_asv11[order(-importance_df_asv11$`%IncMSE`),]
importance_df_asv11

## remove some variables----
model_tb |>
  colnames()

model_tb <- model_tb |>
  dplyr::select( -BP_FC1.55_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=200, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final

### BLOOMER asv555----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv555' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv555' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F) #create a stratified random sample of the data into training test sets (75% of my data)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20) ## the trainControl can be used to speficy the type of resampling
bloo_rf <- train( diff_rclr_time ~., method = 'rf', trControl = control, data = train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) # ntree = 300 

## retrain model on entire dataset----
bloo_final_asv555 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=300, importance=T)
bloo_final_asv555

importance_df_asv555 <- as.data.frame(bloo_final_asv555$finalModel$importance)
importance_df_asv555 <- importance_df_asv555[order(-importance_df_asv555$`%IncMSE`),]
importance_df_asv555

## remove some variables----
model_tb |>
  colnames()

model_tb <- model_tb |>
  dplyr::select(-BP_FC1.55_sm)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

## retrain model on entire dataset
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final



### save the results of the models ------
# Function to process results for a given ASV number
process_results <- function(bloo_final_results, asv_num) {
  results <- bloo_final_results$results %>%
    mutate(asv_num = asv_num)
  return(results)
}

# ASV numbers
asv_numbers <- c("asv62", "asv23", "asv72", "asv11", "asv17", "asv555")

# Process results for each ASV number
all_results <- lapply(asv_numbers, function(asv_num) {
  process_results(get(paste0("bloo_final_", asv_num)), asv_num)
})

# Combine all results into a single data frame
results_models <- bind_rows(all_results)

##plot importance of variables ----
##This importance is a measure of by how much removing a variable decreases accuracy, and vice versa — by how much including a variable increases accuracy. 
## The %IncMSE (Percentage Increase in Mean Squared Error) is a metric used to measure the 
## importance of each predictor (variable) in a random forest model. It indicates how much the mean squared error (MSE) increases when a particular predictor is randomly permuted while keeping all other predictors unchanged.

##dataset from my example ASVs with it's 
im62 <- importance_df_asv62 |>
  dplyr::mutate(asv_num = 'asv62') |>
  rownames_to_column(var = 'variable') |>
  arrange(desc(`%IncMSE`))
  
im23 <- importance_df_asv23 |>
  dplyr::mutate(asv_num = 'asv23') |>
  rownames_to_column(var = 'variable')

im72 <- importance_df_asv72 |>
  dplyr::mutate(asv_num = 'asv72') |>
  rownames_to_column(var = 'variable')

im17 <- importance_df_asv17 |>
  dplyr::mutate(asv_num = 'asv17') |>
  rownames_to_column(var = 'variable')

im11 <- importance_df_asv11 |>
  dplyr::mutate(asv_num = 'asv11') |>
  rownames_to_column(var = 'variable')

im555 <- importance_df_asv555 |>
  dplyr::mutate(asv_num = 'asv555') |>
  rownames_to_column(var = 'variable')

im_all <- im62 |>
  bind_rows(im23) |>
  bind_rows(im72) |>
  bind_rows(im17) |>
  bind_rows(im11) |>
  bind_rows(im555)|>
  as_tibble() 

im_all |>
  colnames() <- c('variable', 'incMSE', 'IncNodePurity', 'asv_num')

im_all |>
  str()

im_all <- im_all |>
  dplyr::mutate(asv_num = as.factor(asv_num))

im_all$asv_num <- im_all$asv_num |>
  factor(levels = c('asv23', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

importance_example_bloo <- im_all |>
  ggplot(aes(fct_infreq(variable), incMSE))+
  geom_col()+
  facet_wrap(~ asv_num)+
  scale_x_discrete(labels = labs_env_models)+
  labs(y = '%IncMSE')+
  theme_bw()+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(),
        text = element_text(size = 6))

importance_example_bloo
# 
# ggsave(importance_example_bloo, filename = 'importance_example_bloo_sm084_all.pdf',
#        path = 'Results/Figures/random_forest/',
#        width = 188, height = 150, units = 'mm')

## PARTIAL DEPENDANCE -----
### we explore how do the variables interact with each other
## we need the pdp package
# https://christophm.github.io/interpretable-ml-book/pdp.html

## models for each ASV
# bloo_final_asv62
# bloo_final_asv23
# bloo_final_asv72
# bloo_final_asv17
# bloo_final_asv11
# bloo_final_asv555

# Compute partial dependence for each variable variable
# partial_dep <- partial(bloo_final, pred.var = "temperature_no_nas")
# 
# # Plot the partial dependence
# partial_dep |>
#   ggplot(aes(x = yhat))+
#   geom_line()
# 
# partial_dep <- partial(bloo_final, pred.var = "fitted_rclr")
# 
# # Plot the partial dependence
# plot(partial_dep)

## loop
#### ASV62------
# Get the names of predictor variables
importance_df_asv62 <- importance_df_asv62 |>
  arrange('%IncMSE')

predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv62)))] 

scale_limits <- bloo_final_asv62$trainingData$.outcome |>
  range()

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv62, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    #scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text(size = 6))
}

# Display all partial dependence plots
pdp_asv62 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv62, filename = 'pdp_asv62.pdf',
       path = 'Results/Figures/random_forest/',
       width = 188, height = 180, units = 'mm')

#### asv23------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv23)))] 

scale_limits <- bloo_final_asv23$trainingData$.outcome |>
  range()

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv23, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    #scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text(size = 6))
}

# Display all partial dependence plots
pdp_asv23 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv23, filename = 'pdp_asv23.pdf',
       path = 'Results/Figures/random_forest/',
       width = 188, height = 180, units = 'mm')

#### asv72------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv72)))]

scale_limits <- bloo_final_asv72$trainingData$.outcome |>
  range()

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv72, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    #scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text (size = 6))
}

# Display all partial dependence plots
pdp_asv72 <- gridExtra::grid.arrange(grobs = partial_plots)


ggsave(pdp_asv72, filename = 'pdp_asv72.pdf',
       path = 'Results/Figures/random_forest/',
       width = 188, height = 180, units = 'mm')


#### asv17------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv17)))] 

scale_limits <- bloo_final_asv17$trainingData$.outcome |>
  range()

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv17, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    #scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text (size = 6))
}

# Display all partial dependence plots
pdp_asv17 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv17, filename = 'pdp_asv17.pdf',
       path = 'Results/Figures/random_forest/',
       width = 188, height = 180, units = 'mm')

#### ASV11----
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv11)))] 
scale_limits <- bloo_final_asv11$trainingData$.outcome |>
  range()

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv11, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    ##scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text(size = 6))
}

# Display all partial dependence plots
pdp_asv11 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv11, filename = 'pdp_asv11.pdf',
       path = 'Results/Figures/random_forest/',
       width = 188, height = 180, units = 'mm')

#### asv555------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv555)))]

scale_limits <- bloo_final_asv555$trainingData$.outcome |>
  range()

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv555, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    ##scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text(size = 6))
}

# Display all partial dependence plots
pdp_asv555 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv555, filename = 'pdp_asv555.pdf',
       path = 'Results/Figures/random_forest/',
       width = 188, height = 180, units = 'mm')









## MODEL WITH NON-SMOOTH ENV VARIABLES BUT SMOOTH RESPONSE VARIABLE----
### define a general trainControl 
control <- trainControl(method = 'cv', number = 10, savePredictions = 'final', returnResamp = 'final', returnData = T, 
                        p = 0.8)

## just see how well the model performs:
### BLOOMER ASV62 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv62' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  #dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv62' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb_62 <- model_tb |>
  bind_cols(asvs_sm_asv) |>
  dplyr::select(-sample_id_num, -BP_FC1.55_no_nas)

### Split data to train and test
indices <- createDataPartition(y=model_tb_62$diff_rclr_time, p =0.8, list = F)
train_tb_62 <- model_tb_62[indices,]
test_tb_62 <- model_tb_62[-indices,]

### Train RF
# control <- trainControl(method = 'cv', 
#                         number = 10,   
#                         returnResamp = 'all',
#                         savePredictions = 'all') 
#control <- trainControl(method = 'cv', number = 10)
bloo_rf_62 <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb_62, metric = 'RMSE', ntree=300, importance=T)

# varUsed(bloo_rf$param[importance])
# plot(varImp(bloo_rf))

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry

importance_df_asv11
importance_df_asv62

bloo_rf$results
bloo_rf$resample |>
  dplyr::reframe(R2_mean = mean(Rsquared),
                 R2_sd = sd(Rsquared))

training_results <- tibble(
  Rsquared_mean = mean(bloo_rf$finalModel$rsq),
  Rsquared_sd = sd(bloo_rf$finalModel$rsq),
  mtry  = bloo_rf$finalModel$mtry,
  type = 'training',
  asv_num = 'asv62'
)

validation_results <- print(bloo_rf, showSD = T) |>
  as.tibble() |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv62') |>
  separate(RMSE, into = c("RMSE_mean", "RMSE_sd"), sep = " \\(") |>
  separate(Rsquared, into = c("Rsquared_mean", "Rsquared_sd"), sep = " \\(") |>
  separate(MAE, into = c("MAE_mean", "MAE_sd"), sep = " \\(") |>
  dplyr::mutate(
    RMSE_sd = as.numeric(gsub("\\)$", "", RMSE_sd)),
    Rsquared_sd = as.numeric(gsub("\\)$", "", Rsquared_sd)),
    MAE_sd = as.numeric(gsub("\\)$", "", MAE_sd)),
    RMSE_mean = as.numeric(RMSE_mean),
    Rsquared_mean = as.numeric(Rsquared_mean),
    MAE_mean = as.numeric(MAE_mean),
    type = 'validation'
  ) |>
  dplyr::select(-RMSE_sd, -RMSE_mean, -MAE_mean, -MAE_sd) |>
  dplyr::mutate(mtry = as.numeric(mtry))

training_results |>
bind_rows(validation_results)

r2_summary_asv62 |>
  ggplot(aes(asv_num, Rsquared_mean))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared_mean-Rsquared_sd, ymax = Rsquared_mean +  Rsquared_sd), width = 0.2) +
  theme_bw()

### BLOOMER asv23 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv23' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  #dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv23' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb_23 <- model_tb |>
  bind_cols(asvs_sm_asv)  |>
  dplyr::select(-sample_id_num, -BP_FC1.55_no_nas)

### Split data to train and test
indices <- createDataPartition(y=model_tb_23$diff_rclr_time, p =0.8, list = F)
train_tb_23 <- model_tb_23[indices,]
test_tb_23 <- model_tb_23[-indices,]

### Train RF
#control <- trainControl(method = 'cv', number = 10)
bloo_rf_23 <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb_23, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation----
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry

r2_summary_asv23 <- print(bloo_rf, showSD = T) |>
  as.tibble() |>
  dplyr::filter(mtry == " 9") |>
  dplyr::mutate(asv_num = 'asv23') |>
  separate(RMSE, into = c("RMSE_mean", "RMSE_sd"), sep = " \\(") |>
  separate(Rsquared, into = c("Rsquared_mean", "Rsquared_sd"), sep = " \\(") |>
  separate(MAE, into = c("MAE_mean", "MAE_sd"), sep = " \\(") |>
  dplyr::mutate(
    RMSE_sd = as.numeric(gsub("\\)$", "", RMSE_sd)),
    Rsquared_sd = as.numeric(gsub("\\)$", "", Rsquared_sd)),
    MAE_sd = as.numeric(gsub("\\)$", "", MAE_sd)),
    RMSE_mean = as.numeric(RMSE_mean),
    Rsquared_mean = as.numeric(Rsquared_mean),
    MAE_mean = as.numeric(MAE_mean)
  )

r2_summary_asv23 |>
  ggplot(aes(asv_num, Rsquared_mean))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared_mean-Rsquared_sd, ymax = Rsquared_mean +  Rsquared_sd), width = 0.2) +
  theme_bw()

### BLOOMER asv72 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv72' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  #dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv72' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)  |>
  dplyr::select(-sample_id_num)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
#control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry

r2_summary_asv72 <- print(bloo_rf, showSD = T) |>
  as.tibble() |>
  dplyr::filter(mtry == " 9") |>
  dplyr::mutate(asv_num = 'asv72') |>
  separate(RMSE, into = c("RMSE_mean", "RMSE_sd"), sep = " \\(") |>
  separate(Rsquared, into = c("Rsquared_mean", "Rsquared_sd"), sep = " \\(") |>
  separate(MAE, into = c("MAE_mean", "MAE_sd"), sep = " \\(") |>
  dplyr::mutate(
    RMSE_sd = as.numeric(gsub("\\)$", "", RMSE_sd)),
    Rsquared_sd = as.numeric(gsub("\\)$", "", Rsquared_sd)),
    MAE_sd = as.numeric(gsub("\\)$", "", MAE_sd)),
    RMSE_mean = as.numeric(RMSE_mean),
    Rsquared_mean = as.numeric(Rsquared_mean),
    MAE_mean = as.numeric(MAE_mean)
  )

r2_summary_asv72 |>
  ggplot(aes(asv_num, Rsquared_mean))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared_mean-Rsquared_sd, ymax = Rsquared_mean +  Rsquared_sd), width = 0.2) +
  theme_bw()

### BLOOMER asv17 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv17' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  #dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv17' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)  |>
  dplyr::select(-sample_id_num)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
#control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry

r2_summary_asv17 <- print(bloo_rf, showSD = T) |>
  as.tibble() |>
  dplyr::filter(mtry == " 2") |> ## I need to write it because it does not recognize the number
  dplyr::mutate(asv_num = 'asv17') |>
  separate(RMSE, into = c("RMSE_mean", "RMSE_sd"), sep = " \\(") |>
  separate(Rsquared, into = c("Rsquared_mean", "Rsquared_sd"), sep = " \\(") |>
  separate(MAE, into = c("MAE_mean", "MAE_sd"), sep = " \\(") |>
  dplyr::mutate(
    RMSE_sd = as.numeric(gsub("\\)$", "", RMSE_sd)),
    Rsquared_sd = as.numeric(gsub("\\)$", "", Rsquared_sd)),
    MAE_sd = as.numeric(gsub("\\)$", "", MAE_sd)),
    RMSE_mean = as.numeric(RMSE_mean),
    Rsquared_mean = as.numeric(Rsquared_mean),
    MAE_mean = as.numeric(MAE_mean)
  )

r2_summary_asv17 |>
  ggplot(aes(asv_num, Rsquared_mean))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared_mean-Rsquared_sd, ymax = Rsquared_mean +  Rsquared_sd), width = 0.2) +
  theme_bw()

### BLOOMER asv11 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  #dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)  |>
  dplyr::select(-sample_id_num)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
#control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry

r2_summary_asv11 <- print(bloo_rf, showSD = T) |>
  as.tibble() |>
  dplyr::filter(mtry == " 2") |>
  dplyr::mutate(asv_num = 'asv11') |>
  separate(RMSE, into = c("RMSE_mean", "RMSE_sd"), sep = " \\(") |>
  separate(Rsquared, into = c("Rsquared_mean", "Rsquared_sd"), sep = " \\(") |>
  separate(MAE, into = c("MAE_mean", "MAE_sd"), sep = " \\(") |>
  dplyr::mutate(
    RMSE_sd = as.numeric(gsub("\\)$", "", RMSE_sd)),
    Rsquared_sd = as.numeric(gsub("\\)$", "", Rsquared_sd)),
    MAE_sd = as.numeric(gsub("\\)$", "", MAE_sd)),
    RMSE_mean = as.numeric(RMSE_mean),
    Rsquared_mean = as.numeric(Rsquared_mean),
    MAE_mean = as.numeric(MAE_mean)
  )

r2_summary_asv11 |>
  ggplot(aes(asv_num, Rsquared_mean))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared_mean-Rsquared_sd, ymax = Rsquared_mean +  Rsquared_sd), width = 0.2) +
  theme_bw()

### BLOOMER asv555 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv555' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  #dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv555' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)  |>
  dplyr::select(-sample_id_num)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p = 0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
#control <- trainControl(method = 'cv', number = 10, savePredictions = 'final', returnResamp = 'final', returnData = T)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry

r2_summary_asv555 <- print(bloo_rf, showSD = T) |>
  as.tibble() |>
  dplyr::filter(mtry == " 2") |>
  dplyr::mutate(asv_num = 'asv555') |>
  separate(RMSE, into = c("RMSE_mean", "RMSE_sd"), sep = " \\(") |>
  separate(Rsquared, into = c("Rsquared_mean", "Rsquared_sd"), sep = " \\(") |>
  separate(MAE, into = c("MAE_mean", "MAE_sd"), sep = " \\(") |>
  dplyr::mutate(
    RMSE_sd = as.numeric(gsub("\\)$", "", RMSE_sd)),
    Rsquared_sd = as.numeric(gsub("\\)$", "", Rsquared_sd)),
    MAE_sd = as.numeric(gsub("\\)$", "", MAE_sd)),
    RMSE_mean = as.numeric(RMSE_mean),
    Rsquared_mean = as.numeric(Rsquared_mean),
    MAE_mean = as.numeric(MAE_mean)
  )

r2_summary_asv555 |>
  ggplot(aes(asv_num, Rsquared_mean))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared_mean-Rsquared_sd, ymax = Rsquared_mean +  Rsquared_sd), width = 0.2) +
  theme_bw()

##plot the at the same time----
r2_summary_all <- r2_summary_asv23 |>
  bind_rows(r2_summary_asv62) |>
  bind_rows(r2_summary_asv72) |>
  bind_rows(r2_summary_asv11) |>
  bind_rows(r2_summary_asv17) |>
  bind_rows(r2_summary_asv555)

r2_summary_all$asv_num <- factor(r2_summary_all$asv_num, levels = c('asv23', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

plot_cross_validation_sd <- r2_summary_all |>
  ggplot(aes(asv_num, Rsquared_mean))+
  geom_point(size = 0.9)+
  geom_errorbar(aes(ymin = Rsquared_mean-Rsquared_sd, ymax = Rsquared_mean +  Rsquared_sd), linewidth = 0.2) +
  labs(x = '', y = expression('R'^2))+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 8))

plot_cross_validation_sd

ggsave(plot_cross_validation_sd, filename = 'plot_cross_validation_diff_clr_env_sm.pdf',
       path = 'Results/Figures/random_forest/',
       width = 88, height = 88, units = 'mm')



## CREATE A SCRITP USING LOOPS TO RUN THE RANDOM FOREST ANALYSIS FOR ALL MY BLOOMERS----

### define a general trainControl 
control <- trainControl(method = 'cv', number = 10, savePredictions = 'final', returnResamp = 'final', returnData = T, 
                        p = 0.8)

## general env data 
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,] # last row of env data is not needed since it won't give us information (there is no next timepoint)

### Step 1 create the tb for modeling----
source('src/create_tb_rf.R')

## FL
#bloo_02$value
model_tb_62 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               asv_num = "asv62", fraction = "0.2")

model_tb_11 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv11", "0.2")
model_tb_555 <- create_model_tb(data_diff = asvs_sm_diff,
                                data_previous_ab = asvs_sm, 
                                "asv555", "0.2")
#model_tb_2 <- create_model_tb("asv2", "0.2")

## PA
#bloo_3$value
model_tb_23 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv23", "3")
model_tb_72 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv72", "3")
model_tb_17 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv17", "3")

### Step 2 perform the random forest save the models----

### Split data to train and test
indices <- createDataPartition(y=model_tb_62$diff_rclr_time, p =0.8, list = F)
train_tb_62 <- model_tb_62[indices,]
test_tb_62 <- model_tb_62[-indices,]

bloo_rf_62 <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb_62, metric = 'RMSE', ntree=300, importance=T)

perform_random_forest <- function(model_tb, response_column_name, response_column, train_percentage = 0.8, ntree = 300) {
  indices <- createDataPartition(y = model_tb[['response_column_name']], p = train_percentage, list = FALSE)
  train_data <- model_tb[indices, ]
  test_data <- model_tb[-indices, ]
  
  # Train the random forest model
  rf_model <- train(
    formula = {{response_column}} ~ .,
    method = 'rf',
    trControl = control,
    data = train_data,
    metric = 'RMSE',
    ntree = ntree,
    importance = TRUE
  )
  
  return(list(model = rf_model, train_data = train_data, test_data = test_data))
}

perform_random_forest(model_tb = model_tb_62, 
                      response_column_name = diff_rclr_time,
                      response_column = diff_rclr_time)
str(model_tb_62)

model_subset <- paste0(model_tb, '$', response_column_name)

### Step 3 extract the R2 from the training and testing of the data----

### Step 4 compare predicted and observed values----

### Step 5 extract the importance results for each taxa----




## compare predicted and observed changes in rclr for all the models created----
## I need to do it specifically for each taxa, since the previous abundance value is specific for each taxa. Then the train dataset needs to
## to be adapted 
# predictions_vs_observed_62 <- extractPrediction(
#   list(asv62 = bloo_rf_62),
#   testX = train_tb_62
# )
# 
# predictions_vs_observed_23 <- extractPrediction(
#   list(asv23 = bloo_rf_23),
#   testX = train_tb_23
# )
# 
# prediction_vs_observed_tb <- predictions_vs_observed_23 |>
#   bind_rows(predictions_vs_observed_62) 
# 
# prediction_vs_observed_tb |>
#   as_tibble() |>
#   ggplot(aes(obs, pred))+
#   geom_point(alpha = 0.8)+
#   geom_smooth(method = 'lm', color = 'black')+
#   facet_wrap(vars(object))+
#   stat_poly_eq(use_label(c('p.value', 'rr.label')),
#                p.digits = 2,
#                rr.digits = 2,
#                method = stats::lm)+
#   labs(x = 'Observed', y = 'Predicted')+
#   theme_bw()+
#   theme(text = element_text(size = 8), strip.background = element_rect(fill = 'transparent'),
#         panel.grid = element_blank())

### i do it in a loop for all my ASVs----
# List of ASV numbers
asv_numbers <- c(62, 23)

# Initialize an empty list to store predictions
predictions_list <- list()

# Loop through each ASV number
for (asv_number in asv_numbers) {
  # Extract prediction for current ASV number
  predictions <- extractPrediction(
    list(asv = get(paste0("bloo_rf_", asv_number))),
    testX = get(paste0("train_tb_", asv_number))
  )
  
  # Modify the object column to include ASV number
  predictions$object <- paste0("asv", asv_number)
  
  # Store predictions in the list with ASV number as the list element name
  predictions_list[[paste0("asv", asv_number)]] <- predictions
}

# Combine predictions into a single tibble
prediction_vs_observed_tb <- bind_rows(predictions_list)

prediction_vs_observed_tb |>
  as_tibble() |>
  ggplot(aes(obs, pred))+
  geom_point(alpha = 0.8)+
  geom_smooth(method = 'lm', color = 'black')+
  facet_wrap(vars(object))+
  stat_poly_eq(use_label(c('p.value', 'rr.label')),
               p.digits = 2,
               rr.digits = 2,
               method = stats::lm)+
  labs(x = 'Observed', y = 'Predicted')+
  theme_bw()+
  theme(text = element_text(size = 8), strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank())

### MODEL ABUNDANCES WITH NON-SMOOTH ENV VARIABLES ----
## just see how well the model performs:
### BLOOMER ASV62 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm |>
  dplyr::filter(asv_num == 'asv62' &
                  fraction == '0.2') |>
  arrange(decimal_date) |>
  dplyr::filter(decimal_date != '2004.07') |> #the first abundance does not have previous timepoint
  ungroup() |>
  dplyr::select(-asv_num,  -fraction) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  dplyr::select(-decimal_date) 

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv62' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date, previous_abund = fitted_rclr )

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)  |>
  dplyr::select(-sample_id_num)

### Split data to train and test
indices <- createDataPartition(y=model_tb$fitted_rclr, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( fitted_rclr  ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry
r2_summary_asv62 <- bloo_rf$results |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv62')

r2_summary_asv62 |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()

### BLOOMER asv23 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm |>
  dplyr::filter(asv_num == 'asv23' &
                  fraction == '3') |>
  arrange(decimal_date) |>
  dplyr::filter(decimal_date != '2004.07') |> #the first abundance does not have previous timepoint
  ungroup() |>
  dplyr::select(-asv_num,  -fraction) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv23' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date, previous_abund = fitted_rclr )

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv) 

### Split data to train and test
indices <- createDataPartition(y=model_tb$fitted_rclr, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( fitted_rclr  ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry
r2_summary_asv23 <- bloo_rf$results |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv23')

r2_summary_asv23 |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()

### BLOOMER asv72 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm |>
  dplyr::filter(asv_num == 'asv72' &
                  fraction == '3') |>
  arrange(decimal_date) |>
  dplyr::filter(decimal_date != '2004.07') |> #the first abundance does not have previous timepoint
  ungroup() |>
  dplyr::select(-asv_num,  -fraction) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv72' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date, previous_abund = fitted_rclr )

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv) 

### Split data to train and test
indices <- createDataPartition(y=model_tb$fitted_rclr, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( fitted_rclr  ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry
r2_summary_asv72 <- bloo_rf$results |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv72')

r2_summary_asv72 |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()

### BLOOMER asv17 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm |>
  dplyr::filter(asv_num == 'asv17' &
                  fraction == '3') |>
  arrange(decimal_date) |>
  dplyr::filter(decimal_date != '2004.07') |> #the first abundance does not have previous timepoint
  ungroup() |>
  dplyr::select(-asv_num,  -fraction) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv17' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date, previous_abund = fitted_rclr )

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv) 

### Split data to train and test
indices <- createDataPartition(y=model_tb$fitted_rclr, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( fitted_rclr  ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry
r2_summary_asv17 <- bloo_rf$results |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv17')

r2_summary_asv17 |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()

### BLOOMER asv11 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '0.2') |>
  arrange(decimal_date) |>
  dplyr::filter(decimal_date != '2004.07') |> #the first abundance does not have previous timepoint
  ungroup() |>
  dplyr::select(-asv_num,  -fraction) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date, previous_abund = fitted_rclr )

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv) 

### Split data to train and test
indices <- createDataPartition(y=model_tb$fitted_rclr, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( fitted_rclr  ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry
r2_summary_asv11 <- bloo_rf$results |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv11')

r2_summary_asv11 |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()

### BLOOMER asv555 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm |>
  dplyr::filter(asv_num == 'asv555' &
                  fraction == '0.2') |>
  arrange(decimal_date) |>
  dplyr::filter(decimal_date != '2004.07') |> #the first abundance does not have previous timepoint
  ungroup() |>
  dplyr::select(-asv_num,  -fraction) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv555' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date, previous_abund = fitted_rclr )

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv) 

### Split data to train and test
indices <- createDataPartition(y=model_tb$fitted_rclr, p =0.8, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 10)
bloo_rf <- train( fitted_rclr  ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry
r2_summary_asv555 <- bloo_rf$results |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv555')

r2_summary_asv555 |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()

##plot the at the same time----
r2_summary_all <- r2_summary_asv23 |>
  bind_rows(r2_summary_asv62) |>
  bind_rows(r2_summary_asv72) |>
  bind_rows(r2_summary_asv11) |>
  bind_rows(r2_summary_asv17) |>
  bind_rows(r2_summary_asv555)

r2_summary_all$asv_num <- factor(r2_summary_all$asv_num, levels = c('asv23', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

r2_summary_all |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()







###TRY THE CODE THAT IS IN CHOLLET 2023----

control <- trainControl(method = 'cv', 
                        number = 10,   
                        returnResamp = 'all',
                        savePredictions = 'all') 
bloo_rf <- train( fitted_rclr  ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)

print(bloo_rf, showSD = T)

# # Define the training control
# train_control <- trainControl(
#   method = 'cv',                   # k-fold cross validation
#   number = 3,                      # number of folds
#   index = folds,                   # provide indices computed with groupKFold for the k-fold CV
#   classProbs = T,                  # should class probabilities be returned
#   summaryFunction = stand.dev ,
#   selectionFunction = lowest       # we want to minimize the metric
# )
# 
# model <- train(fitted_rclr, data = train_tb, metric = selec.metric, method = algorithm, trControl = train.control, tuneGrid = tune.grid)
# # model <- train(x = temp.train[,expl.var], y = temp.train[,list.taxa[j]], metric = selec.metric, method = algorithm, trControl = train.control)
# 
# 
# temp.list[["Trained model"]] <- model
# 
# for(n in 1:length(which.set)){
#   # n <- 1
#   # Observation
#   temp.list[[paste(out[1],which.set[n])]] <- temp.sets[[n]]
#   n.obs <- dim(temp.sets[[n]])[1]
#   
#   # Prediction factors
#   temp.list[[paste(out[2],which.set[n])]] <- predict(model, temp.sets[[n]])
#   
#   # Prediction probabilities
#   temp.list[[paste(out[3],which.set[n])]] <- predict(model, temp.sets[[n]], type = 'prob')
#   
#   # Likelihood
#   likeli <- 1:nrow(temp.sets[[n]])
#   for(i in 1:nrow(temp.sets[[n]])){
#     if(temp.sets[[n]][i,list.taxa[j]] == "present"){
#       likeli[i] <- temp.list[[paste(out[3],which.set[n])]][i, "present"]
#     } else if (temp.sets[[n]][i,list.taxa[j]] == "absent" ){
#       likeli[i] <- temp.list[[paste(out[3],which.set[n])]][i, "absent"]
#     }
#   }
#   likeli[which(likeli < 0.0001)] <- 0.0001 # avoid problems when likelihood too small
#   temp.list[[paste(out[4],which.set[n])]] <- likeli
#   # Performance
#   temp.list[[paste(out[5],which.set[n])]] <- -2 * sum(log(likeli)) / nrow(temp.sets[[n]])
#   
# }
# 
# list.outputs[[j]] <- temp.list
# }
# 
# outputs[[k]] <- list.outputs
# }
# 
# return(outputs)
# }

# PENSAR SI TINDRIA SENTIT O NO------
## MODEL WITH DEGREE OF CHANGE BETWEEN ENV VARIABLES AND CLR SMOOTH-----
## calculate the increase in the env variables too----
env_sm_diff_w <- env_data_new |>
  arrange(decimal_date) |>
  pivot_longer(cols = !decimal_date, names_to = 'variable', values_to = 'values') |>
  left_join(date_tb) |>
  arrange(decimal_date) |>
  group_by(variable) |>
  dplyr::mutate(diff_value = c(NA, diff(values)),
                diff_time = c(NA, diff(date))) |>
  dplyr::mutate(diff_env_time = diff_value/diff_time) |>
  dplyr::mutate(sample_id_num = row_number()) |>
  dplyr::select(-date, -values, -diff_value, -diff_time) |>
  pivot_wider(id_cols = decimal_date, values_from = diff_env_time, names_from = variable) |>
  arrange(decimal_date) 

env_sm_diff_w <- env_sm_diff_w[-1,] #the first row is NA because we can't calculate the difference with the previous time point

### BLOOMER ASV62 (smooth) ----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv62' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  left_join(env_sm_diff_w) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv62' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results)

# View summary of the trained model
summary(custom)

# Plot the model
plot(custom)

##in this case is ntree 200

## retrain model on entire dataset----
bloo_final_asv62_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final_asv62_diff

importance_df_asv62_diff <- as.data.frame(bloo_final_asv62_diff$finalModel$importance)
importance_df_asv62_diff <- importance_df_asv62_diff[order(-importance_df_asv62_diff$`%IncMSE`),]

### BLOOMER asv23----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv23' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  left_join(env_sm_diff_w) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv62' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv) |>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results)

##500 is best

## retrain model on entire dataset----
bloo_final_asv23_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv23_diff

importance_df_asv23_diff <- as.data.frame(bloo_final_asv23_diff$finalModel$importance)
importance_df_asv23_diff <- importance_df_asv23_diff[order(-importance_df_asv23_diff$`%IncMSE`),]

### BLOOMER asv72----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv72' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  left_join(env_sm_diff_w) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv72' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results)

## ntree = 100

## retrain model on entire dataset----
bloo_final_asv72_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=100, importance=T)
bloo_final_asv72_diff

importance_df_asv72_diff <- as.data.frame(bloo_final_asv72_diff$finalModel$importance)
importance_df_asv72_diff <- importance_df_asv72_diff[order(-importance_df_asv72_diff$`%IncMSE`),]
importance_df_asv72_diff

### BLOOMER asv17----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv17' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  left_join(env_sm_diff_w) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv17' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) ##ntree = 500

## retrain model on entire dataset----
bloo_final_asv17_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv17_diff

importance_df_asv17_diff <- as.data.frame(bloo_final_asv17_diff$finalModel$importance)
importance_df_asv17_diff <- importance_df_asv72[order(-importance_df_asv17_diff$`%IncMSE`),]
importance_df_asv17_diff

### BLOOMER asv11----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  left_join(env_sm_diff_w) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) ## ntree = 500

## retrain model on entire dataset----
bloo_final_asv11_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv11_diff

importance_df_asv11_diff <- as.data.frame(bloo_final_asv11_diff$finalModel$importance)
importance_df_asv11_diff <- importance_df_asv11_diff[order(-importance_df_asv11_diff$`%IncMSE`),]
importance_df_asv11_diff

### BLOOMER asv555----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv555' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  left_join(env_sm_diff_w) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv555' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv) |>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F) #create a stratified random sample of the data into training test sets (75% of my data)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20) ## the trainControl can be used to speficy the type of resampling
bloo_rf <- train( diff_rclr_time ~., method = 'rf', trControl = control, data = train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) # ntree = 300 

## retrain model on entire dataset----
bloo_final_asv555_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=300, importance=T)
bloo_final_asv555_diff

importance_df_asv555_diff <- as.data.frame(bloo_final_asv555_diff$finalModel$importance)
importance_df_asv555_diff <- importance_df_asv555_diff[order(-importance_df_asv555_diff$`%IncMSE`),]
importance_df_asv555_diff

### save the results of the models ------
# Function to process results for a given ASV number
process_results <- function(bloo_final_results, asv_num) {
  results <- bloo_final_results$results %>%
    mutate(asv_num = asv_num)
  return(results)
}

# ASV numbers
asv_numbers <- c("asv62", "asv23", "asv72", "asv11", "asv17", "asv555")

# Process results for each ASV number
all_results <- lapply(asv_numbers, function(asv_num) {
  process_results(get(paste0("bloo_final_", asv_num, '_diff')), asv_num)
})

# Combine all results into a single data frame
results_models <- bind_rows(all_results)

##plot importance of variables ----
##This importance is a measure of by how much removing a variable decreases accuracy, and vice versa — by how much including a variable increases accuracy. 
## The %IncMSE (Percentage Increase in Mean Squared Error) is a metric used to measure the 
## importance of each predictor (variable) in a random forest model. It indicates how much the mean squared error (MSE) increases when a particular predictor is randomly permuted while keeping all other predictors unchanged.

##dataset from my example ASVs with it's 
im62 <- importance_df_asv62_diff |>
  dplyr::mutate(asv_num = 'asv62') |>
  rownames_to_column(var = 'variable') |>
  arrange(desc(`%IncMSE`))

im23 <- importance_df_asv23_diff |>
  dplyr::mutate(asv_num = 'asv23') |>
  rownames_to_column(var = 'variable')

im72 <- importance_df_asv72_diff |>
  dplyr::mutate(asv_num = 'asv72') |>
  rownames_to_column(var = 'variable')

im17 <- importance_df_asv17_diff |>
  dplyr::mutate(asv_num = 'asv17') |>
  rownames_to_column(var = 'variable')

im11 <- importance_df_asv11_diff |>
  dplyr::mutate(asv_num = 'asv11') |>
  rownames_to_column(var = 'variable')

im555 <- importance_df_asv555_diff |>
  dplyr::mutate(asv_num = 'asv555') |>
  rownames_to_column(var = 'variable')

im_all <- im62 |>
  bind_rows(im23) |>
  bind_rows(im72) |>
  bind_rows(im17) |>
  bind_rows(im11) |>
  bind_rows(im555)|>
  as_tibble() 

im_all |>
  colnames() <- c('variable', 'incMSE', 'IncNodePurity', 'asv_num')

im_all |>
  str()

im_all <- im_all |>
  dplyr::mutate(asv_num = as.factor(asv_num))

im_all$asv_num <- im_all$asv_num |>
  factor(levels = c('asv23', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

importance_example_bloo <- im_all |>
  ggplot(aes(fct_infreq(variable), incMSE))+
  geom_col()+
  facet_wrap(~ asv_num)+
  scale_x_discrete(labels = labs_env_models)+
  labs(y = '%IncMSE')+
  theme_bw()+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank(),
        text = element_text(size = 6))

importance_example_bloo

# ggsave(importance_example_bloo, filename = 'importance_example_bloo_diff.pdf',
#        path = 'Results/Figures/random_forest/',
#        width = 188, height = 150, units = 'mm')

## PARTIAL DEPENDANCE -----
### we explore how do the variables interact with each other
## we need the pdp package
# https://christophm.github.io/interpretable-ml-book/pdp.html

## loop
#### ASV62------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv62_diff)))] 

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv62_diff, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    #scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text(size = 6))
}

# Display all partial dependence plots
pdp_asv62 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv62, filename = 'pdp_asv62_diff.pdf',
       path = 'Results/Figures/random_forest/diff_env_variables/',
       width = 188, height = 180, units = 'mm')

#### asv23------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv23_diff)))] 


# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv23_diff, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    #scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text(size = 6))
}

# Display all partial dependence plots
pdp_asv23 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv23, filename = 'pdp_asv23_diff.pdf',
       path = 'Results/Figures/random_forest/diff_env_variables/',
       width = 188, height = 180, units = 'mm')

#### asv72------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv72_diff)))]

scale_limits <- bloo_final_asv72_diff$trainingData$.outcome |>
  range()

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv72_diff, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    #scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text (size = 6))
}

# Display all partial dependence plots
pdp_asv72 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv72, filename = 'pdp_asv72_diff.pdf',
       path = 'Results/Figures/random_forest/diff_env_variables/',
       width = 188, height = 180, units = 'mm')

#### asv17------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv17_diff)))] 

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv17_diff, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    #scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text (size = 6))
}

# Display all partial dependence plots
pdp_asv17 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv17, filename = 'pdp_asv17_diff.pdf',
       path = 'Results/Figures/random_forest/diff_env_variables/',
       width = 188, height = 180, units = 'mm')

#### ASV11----
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv11_diff)))] 

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv11_diff, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    ##scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text(size = 6))
}

# Display all partial dependence plots
pdp_asv11 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv11, filename = 'pdp_asv11_diff.pdf',
       path = 'Results/Figures/random_forest/diff_env_variables/',
       width = 188, height = 180, units = 'mm')

#### asv555------
# Get the names of predictor variables
predictor_vars <- predictor_vars[order(match(predictor_vars, rownames(importance_df_asv555_diff)))]

# Create an empty list to store the partial dependence plots
partial_plots <- list()

# Loop through each predictor variable
for(var in predictor_vars) {
  # Compute partial dependence for the current variable
  partial_dep <- partial(bloo_final_asv555_diff, pred.var = var)
  
  # Convert the partial dependence plot to a ggplot object
  partial_gg <- autoplot(partial_dep)
  
  # Set x-axis title based on labs_env_models labeller function
  x_axis_title <- labs_env_models[var]
  
  # Store the ggplot object in the list
  partial_plots[[var]] <- partial_gg +
    ##scale_y_continuous(limits = c(scale_limits))+
    theme_bw()+
    labs(y = 'rCLR change (t-1)', x = x_axis_title)+
    theme(panel.grid = element_blank(),
          text = element_text(size = 6))
}

# Display all partial dependence plots
pdp_asv555 <- gridExtra::grid.arrange(grobs = partial_plots)

ggsave(pdp_asv555, filename = 'pdp_asv555_diff.pdf',
       path = 'Results/Figures/random_forest/diff_env_variables/',
       width = 188, height = 180, units = 'mm')






##### WHAT HAPPENS WITH OTHER BLOOMERS (DESIGN A CODE TO DO IT FOR ALL OF THEM) ####-----
bloo_02$value #discard SAR11 clade "asv8"   "asv5"   "asv3"   "asv2"
bloo_3$value

## ASV1 (in PA and FL)
#### PA----
#fit the loess model
fit <- loess(asv1 ~ decimal_date, data = rclr_time_tb_3, span = 0.084)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

asv1_sm <- predicted_values |>
  as_tibble_col(column_name = 'asv1') |>
  bind_cols(decimal_date_tb) |>
  dplyr::select(decimal_date, 'fitted_rclr' = 'asv1') |>
  dplyr::mutate(asv_num = 'asv1') |>
  dplyr::left_join(fraction_date_3)

rclr_time_tb_3 |>
  dplyr::select(decimal_date, 'rclr' = 'asv1') |>
  left_join(asv1_sm) |>
  pivot_longer(cols = c('rclr', 'fitted_rclr'), names_to = 'approx', values_to = 'abundance_value') |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line(aes(group = approx, color = approx))+
  scale_color_manual(values = c('black', 'grey'))+
  theme_bw()+
  theme(panel.grid = element_blank())



### BLOOMER asv1----
model_tb <- asv1_sm |>
  dplyr::filter(!is.na(diff_time)) |>
  #dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv1' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  left_join(env_sm_diff_w) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv1' &
                  fraction == '0.2') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb <- model_tb |>
  bind_cols(asvs_sm_asv)|>
  dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next

### Split data to train and test
indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
train_tb <- model_tb[indices,]
test_tb <- model_tb[-indices,]

### Train RF
control <- trainControl(method = 'cv', number = 20)
bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
y <- predict(bloo_rf, newdata = test_tb[,-1])
test_tb$prediction <- y
fit <- lm(prediction~diff_rclr_time, data=test_tb)
summary(fit)
importance_df <- as.data.frame(bloo_rf$finalModel$importance)
importance_df <- importance_df[order(-importance_df$`%IncMSE`),]

##try to define best parameters----
# Define train control
control <- trainControl(method = 'cv', number = 20)

# Define the tuning parameter grid
sqrt(ncol(train_tb)) ## the mtry should we a number arround this one

tunegrid <- expand.grid(
  mtry = c(4:6) #, # 
  #ntree = c(100, 200)
)

## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
modellist <- list()
for (i in 1:nrow(tunegrid)) {
  set.seed(100)
  mtry_val <- tunegrid$mtry[i]
  
  for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
    fit <- train(diff_rclr_time ~ ., 
                 data = train_tb, 
                 method = 'rf', 
                 metric = 'RMSE', 
                 tuneGrid = tunegrid,
                 trControl = control)
    
    key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
    modellist[[key]] <- fit
  }
}

# compare results
results <- resamples(modellist)
summary(results)
dotplot(results) ## ntree = 500

## retrain model on entire dataset----
bloo_final_asv1_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv1_diff

importance_df_asv1_diff <- as.data.frame(bloo_final_asv1_diff$finalModel$importance)
importance_df_asv1_diff <- importance_df_asv1_diff[order(-importance_df_asv1_diff$`%IncMSE`),]
importance_df_asv1_diff



### WHAT IF I USE THE CLR INSTEAD OF RCLR----
asv_tab_10y_02_clr

clr_df_inter_deconstand 

## filter it for my bloomers (the one's I want to perform the wavelets analysis on)
clr_03 <- clr_df_inter_deconstand  |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_02, by = c('date')) |>
  dplyr::select(decimal_date, asv_num, clr)

## observe them
clr_03  |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(decimal_date, clr, group = asv_num))+
  geom_line()+
  geom_point()+
  facet_wrap(vars(asv_num))+
  theme_bw()

clr_03_w <- clr_03 |>
  pivot_wider(id_cols = decimal_date, names_from = asv_num, values_from = clr)

#fit the loess model
fit <- loess(asv23 ~ decimal_date, data = clr_03_w, span = 0.01)

# predict the fited values
predicted_values <- predict(fit)

# Obtain the residuals from the loess model
residuals <- residuals(fit)

## check the residuals
## this function gets the residuals and we check that we are not oversmoothing our data, so that there is some remaining variability still in the residuals
# Compute the autocorrelation function (ACF) of the residuals
acf_res <- acf(residuals, main = "Autocorrelation Function of Residuals")

asv23_sm <- predicted_values |>
  as_tibble_col(column_name = 'asv23') |>
  bind_cols(decimal_date_tb) |>
  dplyr::select(decimal_date, 'fitted_rclr' = 'asv23') |>
  dplyr::mutate(asv_num = 'asv23') |>
  dplyr::left_join(fraction_date_3)

## Environmental variables overtime----
env_data_new |>
  #bind_cols(decimal_date_tb) |> 
  pivot_longer(cols = !decimal_date, names_to = 'env_variable', values_to = 'value')|>
  ggplot(aes(decimal_date, value))+
  geom_line()+
  facet_wrap(vars(env_variable), scales = 'free')+
  theme_bw()
  



####### NEW APPROACH USING THE MIKROPML PACKAGE----
library(caret)
library(mikropml)
library(tictoc) ##see how long does it take to run a part of the code
library(purrr) #map function

#input
model_tb_11 #wide data

## preprocess data (remove those variables that have near 0 variance)
model_tb_11_pre <-  preprocess_data(model_tb_11, outcome_colname = 'diff_rclr_time',
                collapse_corr_feats = TRUE, ### perfectly correlated env variables do not add inormation to the model #not the case for my data
                group_neg_corr = TRUE, # none of my variables are in this case
                method =  "center") #if we want to normalize the data 


### define a general trainControl 
control <- trainControl(method = 'cv', number = 10, savePredictions = 'final', returnResamp = 'final', returnData = T, 
                        p = 0.8)

##continue without preprocessing
rf_mikro_11_pre <- run_ml(model_tb_11_pre$dat_transformed, method = 'rf',
       outcome_colname = 'diff_rclr_time',
       find_feature_importance = T,
       cross_val = control,
       training_frac = 0.8,
       perf_metric_name = 'RMSE',
       seed = 1030)

test_hp <- list(ntree = c('100', '200', '300', '500'))

rf_mikro_11 <- run_ml(model_tb_11, method = 'rf',
                      outcome_colname = 'diff_rclr_time',
                      find_feature_importance = T,
                      cross_val = control,
                      training_frac = 0.8,
                      hyperparameters = test_hp,
                      perf_metric_name = 'RMSE',
                      seed = 1030)

rf_mikro_11_ml <- run_ml(model_tb_11, method = 'glmnet',
                      outcome_colname = 'diff_rclr_time',
                      find_feature_importance = T,
                      #cross_val = control,
                      training_frac = 0.8,
                      #perf_metric_name = 'RMSE',
                      seed = 1030)


## see if there are difences between one preprocessed and the one without processing 
rf_mikro_11_pre$performance
rf_mikro_11$performance

rf_mikro_11_pre$feature_importance
rf_mikro_11$feature_importance

rf_mikro_62 <- run_ml(model_tb_62, method = 'rf',
       outcome_colname = 'diff_rclr_time',
       find_feature_importance = T,
       cross_val = control,
       #training_frac = 0.8,
       perf_metric_name = 'RMSE',
       seed = 1030)

rf_mikro_23 <- run_ml(model_tb_23, method = 'rf',
                      outcome_colname = 'diff_rclr_time', #response variable
                      find_feature_importance = T,
                      cross_val = control,
                      #training_frac = 0.8,
                      perf_metric_name = 'RMSE',
                      seed = 1030)

##I run the model with different seeds 
get_hyperparams_list(dataset = model_tb_23, method = 'rf') 

hyperparameter <- list(mtry = c(2, 4, 8))

##try to be able to change the model_tb as an argument of the function.
get_results <- function(seed){
  run_ml(model_tb_23, method = 'rf',
         outcome_colname = 'diff_rclr_time', #response variable
         find_feature_importance = T,
         cross_val = control,
         training_frac = 0.8,
         hyperparameters = hyperparameter,
         perf_metric_name = 'RMSE',
         seed = seed)
}

#tic()
iterative_run_ml_results <- map(c(1:3), get_results) #100 iterations
#toc()
## 7 min/ ASV to run it 100 times.
## do i need to paralelize?
## library(furrr)
## plan ('multicore') ## 
## plan('sequential') we need to go back to the sequential after runing the multicore
iterative_run_ml |>
  str()

get_hyperparams_list(dataset = model_tb_23, method = 'rf')
get_hyperparams_list(dataset = model_tb_62, method = 'rf')

## interesed in the trained model for each seed
performance <- iterative_run_ml_results |>
  map(pluck, 'trained_model') |>
  combine_hp_performance()

plot_hp_performance(performance$dat, mtry, RMSE)

performance$dat |>
  group_by(mtry, RMSE) |>
  dplyr::reframe(mean = mean(RMSE),
                 lquartile = quantile(RMSE, prob = 0.25),
                 uquartile = quantile(RMSE, prob = 0.75)) |>
  top_n(n = 3, mean)

## model performance

model <- rf_mikro_23 

prob <- predict(model$trained_model, model$test_data) #type = 'prob' no em funciona
observed <- model$test_data$diff_rclr_time

  bind_cols(predicted =prob, observed = observed) |>
    arrange(desc(predicted))


