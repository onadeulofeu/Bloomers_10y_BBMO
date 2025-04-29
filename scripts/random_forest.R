# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#packages----
library(caret) #random forest
library(forcats) #reorder in ggplot2
library(pdp) #partial dependence analysis
library(tidyverse)
library(stringr)
library(easystats) # check model assumptions
library(ggpmisc) # add pvalue and r2 to ggplots 

theme_set(theme_bw()) ##by default the theme of my plots will be bw

##labels ----
## prepare environmental variable labels for the plots of this analysis
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

labs_env_models_no_nas <- c( "temperature_no_nas" = 'Temperature (ºC)',
                      "day_length_no_nas" = 'Day length (h)',
                      "bacteria_joint" = 'Total bacterial abundance (cells/mL)',
                      "synechococcus" = 'Synechococcus abundance (cells/mL)',
                      "chla_total_no_nas"  = 'Chl-a (µg/L)' ,
                      "PO4_no_nas" =   '[PO43-] (µM)',
                      "NH4_no_nas"      =   '[NH4+] (µM)',
                      "NO2_no_nas"  =    '[NO2-] (µM)',       
                      "NO3_no_nas"  =   '[NO3-] (µM)',
                      "Si_no_nas"   =   '[SiO4-] (µM)',
                      "PNF_no_nas"   =    'Phototrophic nanoflagellate abundaance (cells mL ^-1)',
                      "cryptomonas_no_nas" = 'Cryptomonas cells/mL',
                      "micromonas_no_nas"  = 'Micromonas cells/mL',
                      "HNF_Micro_no_nas"      = 'Heterotrophic nanoflagellate abundaance (cells mL ^-1)',
                      "ev02_no_nas"   = 'FL community evenness',
                      "ev3_no_nas"    = 'PA community evenness',        
                      "fitted_rclr" = 'previous rCLR')

## we need to create a wide dataset with environmental variables and each ASV that we consider it is a potential bloomer
## Do we need to normalize environmental variables to z-scores? No need
## Is there a way to add the community effect on a blooming event in this test? We added community evenness and then we removed it. We do not expect the community from 1 month before to affect blooming events?
## Could we add eukaryotes data here?

## Caret package https://topepo.github.io/caret/index.html

## when performing a RF the two parameters that most likey have the biggest effect on the final accuracy are:
### mtry: numbr of variables randomly sampled as candidates at each split
### ntree: number of trees to grow
set.seed(100)

## Environmental to model with interpolated missing variables----
env_data_interpolated_values_all <- read.csv2('data/env_data/env_data_interpolated_values_all.csv') |>
  rename(sample_id_num = X)

env_data_interpolated_values_all_red <- env_data_interpolated_values_all |>
  dplyr::select(- decimal_date,  -"PNF2_5um_Micro_no_nas" ,-"PNF_5um_Micro_no_nas",  -"HNF2_5um_Micro_no_nas" ,-"HNF_5um_Micro_no_nas" ) #I need to reduce the explanatory variables nº

## prepare the dataset----
### for the PA bloomers I need the interpolated abundance for those 3 missing samples in this fraction----
bloo_3_inter  <- read.csv('data/wavelet_3_df_deconstand.csv') |> # It is here because I prepared this dataframe for the wavelets analysis
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

# EXAMPLES FOR EACH TYPE OF BLOOM 
### 'asv17' (PA), 'asv11' (PA), 'asv11' (FL), 'asv62' (FL), 'asv72' (FL), 'asv555' (FL))

## ideas data that I should add here: diatoms, ciliates, dinophlgellates, 

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

## asv11----
model_tb <- diff_time_tb_env |>
 # dplyr::left_join(community_evenness_all) |> #calculated in the next chunk
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv11' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -decimal_date, -diff_rclr, -diff_time, -date, -fraction, -abundance_value, -BP_FC1.55_no_nas)

##add previous CLR value as explanatory variable
asvs_prev_abund <- diff_time_tb_env |>
  dplyr::filter(asv_num_f == 'asv11' &
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
bloo_final_asv11 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv11

## once we have trained our model we should not retrain it again (it will lead to different results.)
importance_df_asv11 <- as.data.frame(bloo_final_asv11$finalModel$importance)
importance_df_asv11 <- importance_df_asv11[order(-importance_df_asv11$`%IncMSE`),]
importance_df_asv11

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

im11 <- importance_df_asv11 |>
  dplyr::mutate(asv_num = 'asv11') |>
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
  bind_rows(im11) |>
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
  factor(levels = c('asv11', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

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

## BLOOMER ASV11-----
model_tb <- diff_time_tb_env |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv11' &
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
#### 3 (asv11, as72, asv17)----

#fit the loess model
fit <- loess(asv11 ~ decimal_date, data = rclr_time_tb_3, span = 0.084)

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
  dplyr::left_join(fraction_date_3)

rclr_time_tb_3 |>
  dplyr::select(decimal_date, 'rclr' = 'asv11') |>
  left_join(asv11_sm) |>
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
asvs_sm <- asv11_sm |>
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
                  asv_num == 'asv11' &
                  fraction == '3') |>
  dplyr::select(decimal_date)

asvs_sm_diff |>
  dplyr::filter(asv_num == 'asv11') |>
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

### BLOOMER asv11----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_new_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
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
bloo_final_asv11 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=100, importance=T)
bloo_final_asv11

importance_df_asv11 <- as.data.frame(bloo_final_asv11$finalModel$importance)
importance_df_asv11 <- importance_df_asv11[order(-importance_df_asv11$`%IncMSE`),]
importance_df_asv11

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
asv_numbers <- c("asv62", "asv11", "asv72", "asv11", "asv17", "asv555")

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
  
im11 <- importance_df_asv11 |>
  dplyr::mutate(asv_num = 'asv11') |>
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
  bind_rows(im11) |>
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
  factor(levels = c('asv11', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

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
# bloo_final_asv11
# bloo_final_asv72
# bloo_final_asv17
# bloo_final_asv11
# bloo_final_asv555

#Compute partial dependence for each variable variable

partial_dep <- partial(bloo_final, pred.var = "temperature_no_nas")

##partial_dep <- pdp::partial(bloo_final, pred.var = "temperature_no_nas")

# Plot the partial dependence
partial_dep |>
  ggplot(aes(x = yhat))+
  geom_line()

partial_dep <- partial(bloo_final, pred.var = "fitted_rclr")

# Plot the partial dependence
plot(partial_dep)


#Compute partial dependence for each variable variable
partial_dep <- pdp::partial(bloo_final, pred.var = "temperature_no_nas")

# Plot the partial dependence
partial_dep |>
  ggplot(aes(x = yhat))+
  geom_line()

partial_dep <- pdp::partial(bloo_final, pred.var = "fitted_rclr")

# Plot the partial dependence
plot(partial_dep)


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

#### asv11------
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
    #scale_y_continuous(limits = c(scale_limits))+
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
  geom_errorbar(aes(ymin = Rsquared_mean-Rsquared_sd, ymax = Rsquared_mean +  Rsquared_sd), width = 0.2)

### BLOOMER asv11 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '3') |>
  ungroup() |>
  dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
  #dplyr::select(-contains('_no_nas')) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  #left_join(env_data_new) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '3') |>
  dplyr::arrange(decimal_date) |>
  dplyr::select(-asv_num, -fraction, -decimal_date)

asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one

model_tb_11 <- model_tb |>
  bind_cols(asvs_sm_asv)  |>
  dplyr::select(-sample_id_num, -BP_FC1.55_no_nas)

### Split data to train and test
indices <- createDataPartition(y=model_tb_11$diff_rclr_time, p =0.8, list = F)
train_tb_11 <- model_tb_11[indices,]
test_tb_11 <- model_tb_11[-indices,]

### Train RF
#control <- trainControl(method = 'cv', number = 10)
bloo_rf_11 <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb_11, metric = 'RMSE', ntree=300, importance=T)

##extract the sd between the R2 from the training and the validation----
x <- bloo_rf$bestTune$mtry ##this is the chosen mtry

r2_summary_asv11 <- print(bloo_rf, showSD = T) |>
  as.tibble() |>
  dplyr::filter(mtry == " 9") |>
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
r2_summary_all <- r2_summary_asv11 |>
  bind_rows(r2_summary_asv62) |>
  bind_rows(r2_summary_asv72) |>
  bind_rows(r2_summary_asv11) |>
  bind_rows(r2_summary_asv17) |>
  bind_rows(r2_summary_asv555)

r2_summary_all$asv_num <- factor(r2_summary_all$asv_num, levels = c('asv11', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

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
model_tb_11 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv11", "3")
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
# predictions_vs_observed_11 <- extractPrediction(
#   list(asv11 = bloo_rf_11),
#   testX = train_tb_11
# )
# 
# prediction_vs_observed_tb <- predictions_vs_observed_11 |>
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
asv_numbers <- c(62, 11)

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

### BLOOMER asv11 (smooth) ----
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,]
model_tb <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
                  fraction == '3') |>
  arrange(decimal_date) |>
  dplyr::filter(decimal_date != '2004.07') |> #the first abundance does not have previous timepoint
  ungroup() |>
  dplyr::select(-asv_num,  -fraction) |>
  bind_cols(env_data_interpolated_values_all_red) |>
  dplyr::select(-decimal_date)

##add previous CLR value as explanatory variable
asvs_sm_asv <- asvs_sm |>
  dplyr::filter(asv_num == 'asv11' &
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
r2_summary_asv11 <- bloo_rf$results |>
  dplyr::filter(mtry == x) |>
  dplyr::mutate(asv_num = 'asv11')

r2_summary_asv11 |>
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
r2_summary_all <- r2_summary_asv11 |>
  bind_rows(r2_summary_asv62) |>
  bind_rows(r2_summary_asv72) |>
  bind_rows(r2_summary_asv11) |>
  bind_rows(r2_summary_asv17) |>
  bind_rows(r2_summary_asv555)

r2_summary_all$asv_num <- factor(r2_summary_all$asv_num, levels = c('asv11', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

r2_summary_all |>
  ggplot(aes(asv_num, Rsquared))+
  geom_point()+
  geom_errorbar(aes(ymin = Rsquared - RMSE, ymax = Rsquared + RMSE), width = 0.2) +
  theme_bw()







###TRY THE CODE THAT IS IN CHOLLET 2011----

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

### BLOOMER asv11----
model_tb <- asvs_sm_diff |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num == 'asv11' &
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
bloo_final_asv11_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv11_diff

importance_df_asv11_diff <- as.data.frame(bloo_final_asv11_diff$finalModel$importance)
importance_df_asv11_diff <- importance_df_asv11_diff[order(-importance_df_asv11_diff$`%IncMSE`),]

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
asv_numbers <- c("asv62", "asv11", "asv72", "asv11", "asv17", "asv555")

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

im11 <- importance_df_asv11_diff |>
  dplyr::mutate(asv_num = 'asv11') |>
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
  bind_rows(im11) |>
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
  factor(levels = c('asv11', 'asv62', 'asv72', 'asv17', 'asv11', 'asv555'))

importance_example_bloo <- im_all |>
  ggplot(aes(fct_infreq(variable), incMSE))+
  geom_col()+
  facet_wrap(~ asv_num)+
  scale_x_discrete(labels = labs_env_models)+
  labs(y = '%IncMSE')+
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

#### asv11------
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
    #scale_y_continuous(limits = c(scale_limits))+
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
fit <- loess(asv11 ~ decimal_date, data = clr_03_w, span = 0.01)

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
  dplyr::left_join(fraction_date_3)

## Environmental variables overtime----
env_data_new |>
  #bind_cols(decimal_date_tb) |> 
  pivot_longer(cols = !decimal_date, names_to = 'env_variable', values_to = 'value')|>
  ggplot(aes(decimal_date, value))+
  geom_line()+
  facet_wrap(vars(env_variable), scales = 'free')+
  theme_bw()
  



---------- ####### NEW APPROACH USING THE MIKROPML PACKAGE ######### --------------------

## This is the first try with my 5 example bloomers.

## packages----
library(caret)
library(mikropml)
library(tictoc) ##see how long does it take to run a part of the code
library(purrr) #map function
library(pdp) #partial dependance


## upload data----

### environmental data
## Environmental to model with interpolated missing variables----
env_data_interpolated_values_all <- read.csv2('data/env_data/env_data_interpolated_values_all.csv') |>
  rename(sample_id_num = X)

env_data_interpolated_values_all_red <- env_data_interpolated_values_all |>
  dplyr::select(- decimal_date,  -"PNF2_5um_Micro_no_nas" ,-"PNF_5um_Micro_no_nas",  -"HNF2_5um_Micro_no_nas" ,-"HNF_5um_Micro_no_nas" ) #I need to reduce the explanatory variables nº

### ASVs rCLR smooth as response variables
### for the PA bloomers I need the interpolated abundance for those 3 missing samples in this fraction----
bloo_3_inter  <- read.csv('data/wavelet_3_df_deconstand.csv') |> # It is here because I prepared this dataframe for the wavelets analysis
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

### asvs_sm

### rCLR increase t2-t1/time degree of change from one point and the next


### types of blooms
summary_types_of_blooms <- read.csv('results/tables/summary_types_of_blooms.csv') |>
  as_tibble() |>
  dplyr::select(-X) |>
  dplyr::mutate(fraction = as.character(fraction))

## load functions----
source('src/create_tb_rf.R')
source('src/run_ml_and_create_tibble.R')
source('src/run_ml_importance_and_create_tibble.R')
source('src/train_rf_and_get_importance.R')
source('src/train_rf_and_save_pdp.R')
source('src/plot_partia_dependence.R')

### 1. PREPARE THE INPUT DATA FOR MODELING ------
## general env data 
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,] # last row of env data is not needed since it won't give us information (there is no next timepoint)

## FL----
#bloo_02$value
model_tb_62 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               asv_num = "asv62", fraction = "0.2")

model_tb_11 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv11", "0.2")
# model_tb_555 <- create_model_tb(data_diff = asvs_sm_diff,
#                                 data_previous_ab = asvs_sm, 
#                                 "asv555", "0.2") # this taxa is discarded since it has > 80% the rCLR = 0.

#model_tb_2 <- create_model_tb("asv2", "0.2")

## PA----
#bloo_3$value
# model_tb_11 <- create_model_tb(data_diff = asvs_sm_diff,
#                                data_previous_ab = asvs_sm,
#                                "asv11", "3")
model_tb_23 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv23", "3")
model_tb_72 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv72", "3")
model_tb_17 <- create_model_tb(data_diff = asvs_sm_diff,
                               data_previous_ab = asvs_sm,
                               "asv17", "3")

# The inputs are
## model_tb_11, model_tb_62... 


# 2. PREPROCESS DATA (remove those variables that have near 0 variance, variables that perfectly correlate... )-----
## It is not the case for my data so I can skip this step
### examples 
# model_tb_11_pre <- preprocess_data(model_tb_11, outcome_colname = 'diff_rclr_time')
# 
# model_tb_11_pre <-  preprocess_data(model_tb_11, outcome_colname = 'diff_rclr_time',
#                 collapse_corr_feats = TRUE, ### perfectly correlated env variables do not add inormation to the model #not the case for my data
#                 group_neg_corr = TRUE, # none of my variables are in this case
#                 method =  "center") #if we want to normalize the data 

##continue without preprocessing

# 3. MODEL ENGINEERING -----
## In this step I will explore my models, change parammeters, explore the differences to be sure that the model I'm building is correct.
### I do not perform the feature importance so that the model runs faster and then I can evaluate which parameters I should use
### define a general trainControl for the parammeter tunning with not a lot of cv so that it runs faster 
  control_cv <- trainControl(#method = 'cv', number = 10, 
                             savePredictions = 'final', 
                             returnResamp = 'final', 
                             returnData = T, 
                             p = 0.8,
                             method = "cv", 
                             number = 10#, # repeats = 3 tuneGrid = tibble(mtry = 4)
                            )
  
  rf_11 <- run_ml(model_tb_11, method = 'rf',
                   outcome_colname = 'diff_rclr_time',
                   find_feature_importance = F,
                   cross_val = control_cv,
                   training_frac = 0.8,
                  ntree = 1000, # more trees more accuracy 
                  #mtry = 2,
                   perf_metric_name = 'RMSE',
                   seed = 1030)
  
  rf_62 <- run_ml(model_tb_62, method = 'rf',
                  outcome_colname = 'diff_rclr_time',
                  find_feature_importance = F,
                  cross_val = control_cv,
                  training_frac = 0.8,
                  ntree = 1000,
                  #mtry = 2,
                  perf_metric_name = 'RMSE',
                  seed = 1030)
  
  rf_72 <- run_ml(model_tb_72, method = 'rf',
                  outcome_colname = 'diff_rclr_time',
                  find_feature_importance = F,
                  cross_val = control_cv,
                  training_frac = 0.8,
                  ntree = 1000,
                  #mtry = 2,
                  perf_metric_name = 'RMSE',
                  seed = 1030)
  
  rf_17 <- run_ml(model_tb_17, method = 'rf',
                  outcome_colname = 'diff_rclr_time',
                  find_feature_importance = F,
                  cross_val = control_cv,
                  training_frac = 0.8,
                  ntree = 1000,
                  #mtry = 2,
                  perf_metric_name = 'RMSE',
                  seed = 1030)
  
  rf_11 <- run_ml(model_tb_11, method = 'rf',
                  outcome_colname = 'diff_rclr_time',
                  find_feature_importance = F,
                  cross_val = control_cv,
                  training_frac = 0.8,
                  ntree = 1000,
                  #mtry = 2,
                  perf_metric_name = 'RMSE',
                  seed = 1030)
  
  ### explore the model performance and decide which parameters I will use before scaling to 100 cv and 100 seeds----
 #  rf_11$performance
 #  rf_11$trained_model
 #  
 # rf_performance_asv11  <- rf_11$performance |>
 #    as_tibble() |>
 #    dplyr::mutate(asv_num = 'asv11')
 #  rf_performance_asv11 <- rf_62$performance |>
 #    as_tibble() |>
 #    dplyr::mutate(asv_num = 'asv11')
 #  rf_72$performance
 #  rf_17$performance
 #  rf_11$performance
  
  # Create a function to create tibble from rf performance
  # create_rf_tibble <- function(rf, asv_num) {
  #   rf$performance %>%
  #     as_tibble() %>%
  #     mutate(asv_num = asv_num)
  # }
  # # 
  # mean(rf_11$trained_model$resample$Rsquared)
  # rf_11$performance
  # 
  # rf_11$feature_importance
  # rf_11$trained_model$results
  # rf_11$trained_model$bestTune
  
  # bloo_final_asv72$finalModel$importance
  # 
  # varImpPlot(bloo_final_asv72$finalModel)
  # 
  # plot(bloo_final_asv72$finalModel)

# # plot.train(bloo_final_asv11)
#   boot_perf <- bootstrap_performance( rf_11,
#                                      outcome_colname = "diff_rclr_time",
#                                      bootstrap_times = 100, alpha = 0.05
#   )
#   boot_perf
#   

  # # Create tibbles for each random forest
  # rf_performance_asv11 <- create_rf_tibble(rf_11, 'asv11')
  # rf_performance_asv62 <- create_rf_tibble(rf_62, 'asv62')
  # rf_performance_72 <- create_rf_tibble(rf_72, 'asv72')
  # rf_performance_17 <- create_rf_tibble(rf_17, 'asv17')
  # rf_performance_11 <- create_rf_tibble(rf_11, 'asv11')
  # 
  # # Combine all tibbles into one
  # combined_rf_performance <- bind_rows(
  #   rf_performance_asv11,
  #   rf_performance_asv62,
  #   rf_performance_72,
  #   rf_performance_17,
  #   rf_performance_11
  # )
  # 
  # # Plot the data
  # combined_rf_performance |>
  #   pivot_longer(cols = contains('RMSE'), values_to = 'value', names_to = 'type') |>
  #   ggplot( aes(x = asv_num, y = value)) +
  #   geom_point(aes(color = type), position = position_dodge(width = 0.5)) +
  #   scale_color_manual(values = c('darkblue', 'grey'))+
  #   labs(x = "ASV Number", y = "RMSE", color = '')+
  #   theme_bw()+
  #   theme(strip.background = element_rect(fill = 'transparent'),
  #         panel.grid = element_blank(),
  #         text = element_text(size = 6))
  #   
  # combined_rf_performance |>
  #   #pivot_longer(cols = contains('RMSE'), values_to = 'value', names_to = 'type') |>
  #   ggplot( aes(x = asv_num, y = Rsquared)) +
  #   geom_point() +
  #   #scale_color_manual(values = c('darkblue', 'grey'))+
  #   labs(x = "ASV Number", y = "RMSE", color = '')+
  #   theme_bw()+
  #   theme(strip.background = element_rect(fill = 'transparent'),
  #         panel.grid = element_blank(),
  #         text = element_text(size = 6))
  
  
  
# 4. EVALUATION OF THE MODELS -----
  ### I do it 100 different times so that I get different values and I can evaluate the performance. Values from the training and the final model should be similar.

# Run the process 10 times with different seeds and create tibbles
num_runs <- 100
seed_list <- 101:201  # Example list of 100 different seeds
#seed_list <- 101:121

## asv11----
# Initialize a list to store the tibbles
rf_performance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_performance <- run_ml_and_create_tibble(model_tb_11, 'asv11', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_performance_list[[paste0("rf_performance_", i)]] <- rf_performance
}

combined_rf_performance_11 <- bind_rows(rf_performance_list)

## asv62----
# Initialize a list to store the tibbles
rf_performance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_performance <- run_ml_and_create_tibble(model_tb_62, 'asv62', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_performance_list[[paste0("rf_performance_", i)]] <- rf_performance
}

# Bind all tibbles into one
combined_rf_performance_62 <- bind_rows(rf_performance_list)

## asv72----
# Initialize a list to store the tibbles
rf_performance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_performance <- run_ml_and_create_tibble(model_tb_72, 'asv72', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_performance_list[[paste0("rf_performance_", i)]] <- rf_performance
}

# Bind all tibbles into one
combined_rf_performance_72 <- bind_rows(rf_performance_list)

## asv17----
# Initialize a list to store the tibbles
rf_performance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_performance <- run_ml_and_create_tibble(model_tb_17, 'asv17', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_performance_list[[paste0("rf_performance_", i)]] <- rf_performance
}

# Bind all tibbles into one
combined_rf_performance_17 <- bind_rows(rf_performance_list)

## asv11----
# Initialize a list to store the tibbles
rf_performance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_performance <- run_ml_and_create_tibble(model_tb_11, 'asv11', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_performance_list[[paste0("rf_performance_", i)]] <- rf_performance
}

# Bind all tibbles into one
combined_rf_performance_11 <- bind_rows(rf_performance_list)



##bind all tibbles from all the different ASVs----
combined_rf_performance_all <- combined_rf_performance_62 |>
  bind_rows(combined_rf_performance_23) |>
  bind_rows(combined_rf_performance_72) |>
  bind_rows(combined_rf_performance_17) |>
  bind_rows(combined_rf_performance_11)


## add some infromation about my bloomers before ploting
combined_rf_performance_all <- combined_rf_performance_all |>
  dplyr::mutate(fraction = case_when(asv_num %in% c('asv62', 'asv11') ~ '0.2',
                                     asv_num %in% c( 'asv72', 'asv17', 'asv23') ~ '3')) |>
  left_join(summary_types_of_blooms)


##plot RMSE and cv RMSE -----
RMSE_cv <- combined_rf_performance_all %>%
  pivot_longer(cols = contains('RMSE'), values_to = 'value', names_to = 'type') %>%
  ggplot(aes(x = asv_num, y = value, fill = type, shape = recurrency)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 1, color = 'black') +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  scale_fill_manual(values = c('darkblue', 'grey')) +
  labs(x = "ASV Number", y = "RMSE", fill = 'Type', shape = 'Reccurrent taxa') +
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        text = element_text(size = 6))

# ggsave(RMSE_cv, filename = 'RMSE_cv.pdf',
#        path = 'Results/Figures/random_forest/',
#        width = 120, height = 100, units = 'mm')


##plot R2 and resample R2-----
Rsquared_cv <- combined_rf_performance_all |>
  pivot_longer(cols = contains('Rsquared'), values_to = 'value', names_to = 'type') |>
  ggplot( aes(x = asv_num, y = value, fill = type, shape = recurrency)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 1, color = 'black') +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  scale_fill_manual(values = c('darkblue', 'grey')) +
  scale_y_continuous(limits = c(0,1))+
  #scale_color_manual(values = c('darkblue', 'grey'))+
  labs(x = "ASV Number", y = expression('R'^2), color = 'Type', shape = 'Reccurrent taxa')+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        text = element_text(size = 6))
# # 
# ggsave(Rsquared_cv, filename = 'Rsquared_cv.pdf',
#        path = 'Results/Figures/random_forest/',
#        width = 120, height = 100, units = 'mm')





# 5. FEATURE IMPORTANCE ----
### now that I have evaluated my models I compute the feature importance. 
### Understanding the results feature_importance: If feature importances were calculated, a data frame 
### where each row is a feature or correlated group. The columns are the performance metric of the 
### permuted data, the difference between the true performance metric and the performance metric of the 
### permuted data (true - permuted), the feature name, the ML method, the performance metric name, and 
### the seed (if provided). For AUC and RMSE, the higher perf_metric_diff is, the more important that 
### feature is for predicting the outcome. For log loss, the lower perf_metric_diff is, the more 
### important that feature is for predicting the outcome.
# 
# There are several columns:
# perf_metric: The performance value of the permuted feature.
# perf_metric_diff: The difference between the performance for the actual and permuted data (i.e. test performance minus permuted performance). Features with a larger perf_metric_diff are more important.
# pvalue: the probability of obtaining the actual performance value under the null hypothesis.
# lower: the lower bound for the 95% confidence interval of perf_metric.
# upper: the upper bound for the 95% confidence interval of perf_metric.
# feat: The feature (or group of correlated features) that was permuted.
# method: The ML method used.
# perf_metric_name: The name of the performance metric represented by perf_metric & perf_metric_diff.
# seed: The seed (if set).

# Run the process 10 times with different seeds and create tibbles
num_runs <- 10
seed_list <- 101:115  # Example list of 100 different seeds

## asv23----
# Initialize a list to store the tibbles
rf_importance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_importance <- run_ml_importance_and_create_tibble(model_tb_02_asv11, 'asv11', seed = 101)
  
  # Store the tibble in the list with a unique name
  rf_importance_list[[paste0("rf_importance_", i)]] <- rf_importance
}

combined_rf_importance_11 <- bind_rows(rf_importance_list)

## asv62----
# Initialize a list to store the tibbles
rf_importance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_importance <- run_ml_importance_and_create_tibble(model_tb_62, 'asv62', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_importance_list[[paste0("rf_importance_", i)]] <- rf_importance
}

combined_rf_importance_62 <- bind_rows(rf_importance_list)


## asv72----
# Initialize a list to store the tibbles
rf_importance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_importance <- run_ml_importance_and_create_tibble(model_tb_72, 'asv72', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_importance_list[[paste0("rf_importance_", i)]] <- rf_importance
}

combined_rf_importance_72 <- bind_rows(rf_importance_list)

## asv17----
# Initialize a list to store the tibbles
rf_importance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_importance <- run_ml_importance_and_create_tibble(model_tb_17, 'asv17', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_importance_list[[paste0("rf_importance_", i)]] <- rf_importance
}

combined_rf_importance_17 <- bind_rows(rf_importance_list)
# 


## asv11----
# Initialize a list to store the tibbles
rf_importance_list <- list()

# Loop through each run
for (i in 1:num_runs) {
  # Run ml and create tibble
  rf_importance <- run_ml_importance_and_create_tibble(model_tb_11, 'asv11', seed_list[i])
  
  # Store the tibble in the list with a unique name
  rf_importance_list[[paste0("rf_importance_", i)]] <- rf_importance
}

combined_rf_importance_11 <- bind_rows(rf_importance_list)


##bind all tibbles from all the different ASVs----
combined_rf_importance_all <- combined_rf_importance_62 |>
  bind_rows(combined_rf_importance_23) |>
  bind_rows(combined_rf_importance_72) |>
  bind_rows(combined_rf_importance_17) |>
  bind_rows(combined_rf_importance_11)

## add some information about my bloomers before ploting
combined_rf_importance_all <- combined_rf_importance_all |>
  dplyr::mutate(fraction = case_when(asv_num %in% c('asv62', 'asv11') ~ '0.2',
                                     asv_num %in% c('asv72', 'asv17', 'asv23') ~ '3')) |>
  left_join(summary_types_of_blooms)  

##plot
combined_rf_importance_all |>
  dplyr::mutate(feat = as.factor(feat)) |>
  # dplyr::group_by(asv_num, seed) |>
  # slice_max(order_by = abs(perf_metric_diff), n = 10) |>
  ggplot(aes(fct_infreq(feat), perf_metric_diff, shape = recurrency))+
  scale_x_discrete(labels = labs_env_models)+
  geom_point()+
  geom_boxplot(alpha = 0.6)+
  labs(y = 'Perf metric difference', x = 'Explanatory variables')+
  facet_grid(vars(asv_num), scales = 'free')+
  coord_flip()+
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = 'transparent'))

combined_rf_importance_all |>
  dplyr::mutate(feat = as.factor(feat)) |>
  # dplyr::group_by(asv_num, seed) |>
  # slice_max(order_by = abs(perf_metric_diff), n = 5) |>
  ggplot(aes(fct_infreq(feat), perf_metric, shape = recurrency))+
  scale_x_discrete(labels = labs_env_models)+
  geom_point()+
  geom_boxplot(alpha = 0.6)+
  labs(y = 'Perf metric difference', x = 'Explanatory variables')+
  facet_grid(vars(asv_num), scales = 'free')+
  coord_flip()+
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = 'transparent'))


## I run the train with different seeds so that I get different splits of my data and extract the importance of variables ----
control_cv <- trainControl(#method = 'cv', number = 10, 
  savePredictions = 'final', 
  returnResamp = 'final', 
  returnData = T, 
  p = 0.8,
  method = "cv", 
  number = 10#, # repeats = 3 tuneGrid = tibble(mtry = 4)
)

model_tb_list <- list(model_tb_62, model_tb_11, model_tb_17, model_tb_72, model_tb_23)

# Names for the elements
asv_names <- c("asv62", "asv11", "asv17", 'asv72', 'asv23')

# Assign names to the list elements
names(model_tb_list) <- asv_names

# Seed list for cross-validation
seed_list <- 101:201 #100 different splits

# Usage example:
importance_result <- train_rf_and_get_importance(model_tb_list, 'diff_rclr_time', seed_list)

## plot the importance results----
### add some information about my bloomers before ploting
importance_result <- importance_result |>
  dplyr::mutate(fraction = case_when(asv_num %in% c('asv62', 'asv11') ~ '0.2',
                                     asv_num %in% c('asv72', 'asv17', 'asv23') ~ '3')) |>
  left_join(summary_types_of_blooms)  |>
  dplyr::mutate(feat = as.factor(env_variable))

### varIMP: All measures of importance are scaled to have a maximum value of 100, unless the scale argument of varImp.train is set to FALSE.
importance_result |>
  ggplot()

importance_result |>
  colnames()

importance_result |>
  dplyr::mutate(feat = as.factor(env_variable)) |>
  # dplyr::group_by(asv_num, seed) |>
  # slice_max(order_by = abs(perf_metric_diff), n = 10) |>
  ggplot(aes(fct_infreq(env_variable), Overall, shape = recurrency))+
  scale_x_discrete(labels = labs_env_models)+
  geom_point()+
  geom_boxplot(alpha = 0.6)+
  labs(y = 'Scaled Variable Importance', x = 'Explanatory variables', shape = 'Taxa recurrency')+
  facet_grid(vars(asv_num), scales = 'free')+
  coord_flip()+
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = 'transparent'))

### %IncMSE
importance_result_filt <- importance_result |>
  dplyr::group_by(asv_num, env_variable) |>
  slice_max(perc_IncMSE, n = 1) |>
  ungroup()

perc_importance <- importance_result |>
  # dplyr::group_by(asv_num, seed) |>
  # slice_max(order_by = abs(perf_metric_diff), n = 10) |>
  ggplot(aes(fct_infreq(env_variable), perc_IncMSE, shape = recurrency))+
  scale_x_discrete(labels = labs_env_models_no_nas)+
  geom_col(data = importance_result_filt, aes(x = env_variable, y= perc_IncMSE), alpha = 0.4)+
  geom_point(size = 1, alpha = 0.8)+
  geom_boxplot(alpha = 0.6)+
  labs(y = '%IncMSE', x = 'Explanatory variables', shape = 'Taxa recurrency')+
  facet_wrap(vars(asv_num))+
  coord_flip()+
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 6))

perc_importance

# ggsave(perc_importance, filename = 'perc_importance.pdf',
#        path = 'Results/Figures/random_forest/',
#        width = 180, height = 100, units = 'mm')

### ImportanceSD ???



# 6. PARTIAL DEPENDANCE-----
### we explore how do the variables interact with each other
### https://christophm.github.io/interpretable-ml-book/pdp.html
## One important thing of this analysis is to see the direction of the relationship.
## Ecological significance of the results based on the previous knowledge of this taxa.

#### In this case we use the whole dataset to create the model, because we want to gain more confidence on the inferred patterns.

### Save pdp data for plotting
outcome_colname <- "diff_rclr_time"

predictors <- c("bacteria_joint", "synechococcus", "temperature_no_nas", "day_length_no_nas",  "chla_total_no_nas",  "PO4_no_nas" ,       
                 "NH4_no_nas" ,  "NO2_no_nas"  , "NO3_no_nas" , "Si_no_nas" , "PNF_Micro_no_nas",  "cryptomonas_no_nas",
                "micromonas_no_nas" , "HNF_Micro_no_nas"  , "fitted_rclr" )

seed_list <- 1:10

train_rf_and_save_pdp(data = model_tb_list$asv62, outcome_colname = outcome_colname, predictors = predictors, seed_list = seed_list, 
                      grid_resolution = 119, asv_num = 'asv62')

train_rf_and_save_pdp(data = model_tb_list$asv11, outcome_colname = outcome_colname, predictors = predictors, seed_list = seed_list, 
                      grid_resolution = 119, asv_num = 'asv11')

train_rf_and_save_pdp(data = model_tb_list$asv23, outcome_colname = outcome_colname, predictors = predictors, seed_list = seed_list, 
                      grid_resolution = 119, asv_num = 'asv23')

train_rf_and_save_pdp(data = model_tb_list$asv72, outcome_colname = outcome_colname, predictors = predictors, seed_list = seed_list, 
                      grid_resolution = 119, asv_num = 'asv72')

train_rf_and_save_pdp(data = model_tb_list$asv17, outcome_colname = outcome_colname, predictors = predictors, seed_list = seed_list, 
                      grid_resolution = 119, asv_num = 'asv17')

##read RData----
pdp_asv11_data <- readRDS("results/figures/random_forest/pdp/pdp_asv11_data.RDS")
pdp_asv62_data <- readRDS("results/figures/random_forest/pdp/pdp_asv62_data.RDS")
pdp_asv23_data <- readRDS("results/figures/random_forest/pdp/pdp_asv23_data.RDS")
pdp_asv72_data <- readRDS("results/figures/random_forest/pdp/pdp_asv72_data.RDS")
pdp_asv17_data <- readRDS("results/figures/random_forest/pdp/pdp_asv17_data.RDS")

### 2. compute the partial dependence for each seed in the same plot. 

## extraction of the importance variables result to order the partial dependence plots by importance----
importance_order <- combined_rf_importance_all |>
  dplyr::group_by(asv_num, feat) |>
  dplyr::reframe(mean_diff = mean(perf_metric_diff)) |>
  dplyr::arrange(abs(mean_diff)) 

importance_order |>
  dplyr::filter(asv_num == 'asv62') |>
  slice_tail(n = 3)

# Call the function with the pdp data list----
plot_partial_dependence(pdp_asv62_data, importance_df =  importance_order, asv_num = 'asv62', num_plots = 3)
plot_partial_dependence(pdp_asv11_data, importance_df =  importance_order, asv_num = 'asv11', num_plots = 3)
plot_partial_dependence(pdp_asv72_data, importance_df =  importance_order, asv_num = 'asv72', num_plots = 3)
plot_partial_dependence(pdp_asv17_data, importance_df =  importance_order, asv_num = 'asv17', num_plots = 3)
plot_partial_dependence(pdp_asv23_data, importance_df =  importance_order, asv_num = 'asv23', num_plots = 3)

# Composition with all the top 3 effects on each blooming example----
pdp_bloomers_example <- gridExtra::grid.arrange(pdp_all_plots_grid_asv62,
                        pdp_all_plots_grid_asv72,
                        pdp_all_plots_grid_asv23, 
                        pdp_all_plots_grid_asv17, 
                        pdp_all_plots_grid_asv11, 
                        ncol = 1)

ggsave(pdp_bloomers_example, filename = 'pdp_bloomers_example.pdf',
       path = 'Results/Figures/random_forest/pdp/',
       width = 188, height = 220, units = 'mm')


#--------------- ######## RANDOM FOREST FOR THE OTHER POTENTIAL BLOOMERS IN THE BBMO10Y ########### -----------
## Taxa that we can perform the random forest should have less than 80% of 0 in the dataset-------
occurrence_bloo_bbmo |>
  colnames()

bloo_rf <- occurrence_bloo_bbmo  |>
  dplyr::filter(occurrence_perc > 0.2) |>
  dplyr::select(asv_num, fraction)

bloo_02

bloo_rf_02 <- bloo_rf |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |> # SAR11 clade not real bloomers just increasing when the others decrease 
  dplyr::filter(asv_num %in% bloo_02$value)

bloo_rf_3 <- bloo_rf |>
    dplyr::filter(fraction == '3') |>
    dplyr::filter(asv_num %in% bloo_3$value)

##35 bloomers I could build a random forest with------
10+25

#### first I need to smooth my response variable (ASV rCLR)-----

# Define a function to smooth ASV data
smooth_asv <- function(data, asv_num) {
  # Fit the loess model
  fit <- loess(data[[asv_num]] ~ decimal_date, data = data, span = 0.084)
  
  # Predict the fitted values
  predicted_values <- predict(fit)
  
  # Create a tibble with the smoothed values
  smoothed_data <- tibble::tibble(
    decimal_date = data$decimal_date,
    fitted_rclr = predicted_values,
    asv_num = asv_num
  )
  
  return(smoothed_data)
}

# List of ASVs to smooth (FL)----
asv_list <- as.list(bloo_rf_02$asv_num)  # Add more ASVs as needed

# Apply the smoothing function to each ASV
smoothed_data_list <- lapply(asv_list, function(asv) smooth_asv(rclr_time_tb_02, asv))

# Combine the smoothed data into a single tibble
smoothed_combined_02 <- do.call(bind_rows, smoothed_data_list) |>
  dplyr::mutate(fraction = '0.2')

# List of ASVs to smooth (PA)-----
asv_list <- as.list(bloo_rf_3$asv_num)  # Add more ASVs as needed

# Apply the smoothing function to each ASV
smoothed_data_list <- lapply(asv_list, function(asv) smooth_asv(rclr_time_tb_3, asv))

# Combine the smoothed data into a single tibble
smoothed_combined_3 <- do.call(bind_rows, smoothed_data_list) |>
  dplyr::mutate(fraction = '3')

### calculate change in rCLR----
## calculate the increase in the new smoothed variables----
asvs_sm_diff <- smoothed_combined_02 |>
  bind_rows(smoothed_combined_3) |>
  arrange(decimal_date) |>
  left_join(date_tb) |>
  arrange(decimal_date) |>
  group_by(asv_num, fraction) |>
  dplyr::mutate(diff_rclr = c(NA, diff(fitted_rclr)),
                diff_time = c(NA, diff(date))) |>
  dplyr::mutate(diff_rclr_time = diff_rclr/diff_time) |>
  dplyr::mutate(sample_id_num = row_number()) 

##previous rCLR value
 asvs_sm <- asvs_sm_diff |>
   dplyr::select(asv_num, fraction, decimal_date, fitted_rclr) |>
   ungroup()

 
### 1. PREPARE THE INPUT DATA FOR MODELING ------
## general env data 
env_data_interpolated_values_all_red <- env_data_interpolated_values_all_red[-120,] # last row of env data is not needed since it won't give us information (there is no next timepoint)

## FL----
bloo_rf_02

 # Loop through each asv_num
 for (asv_num in bloo_rf_02$asv_num) {
   # Create the asv_num_id with quotes
   asv_num_id <- asv_num
   
   # Create the model with the appropriate name
   assign(paste0('model_tb_02_', asv_num), 
          create_model_tb(data_diff = asvs_sm_diff,
                          data_previous_ab = asvs_sm,
                          asv_num = asv_num_id, fraction = "0.2"))
 }
 
##one by one
# model_tb_15 <- create_model_tb(data_diff = asvs_sm_diff,
#                                data_previous_ab = asvs_sm,
#                                asv_num =  'asv23', fraction = "0.2")
 
 
## PA----
bloo_rf_3

 # Loop through each asv_num
 for (asv_num in bloo_rf_3$asv_num) {
   # Create the asv_num_id with quotes
   asv_num_id <- asv_num
   
   # Create the model with the appropriate name
   assign(paste0('model_tb_3_', asv_num), 
          create_model_tb(data_diff = asvs_sm_diff,
                          data_previous_ab = asvs_sm,
                          asv_num = asv_num_id, fraction = "3"))
 }
 
# The inputs are
## model_tb_02_11, model_tb_3_62... (for example)

 
 # 2. PREPROCESS DATA (remove those variables that have near 0 variance, variables that perfectly correlate... )-----
 ## It is not the case for my data so I can skip this step
 ### examples 
 # model_tb_11_pre <- preprocess_data(model_tb_11, outcome_colname = 'diff_rclr_time')
 # 
 # model_tb_11_pre <-  preprocess_data(model_tb_11, outcome_colname = 'diff_rclr_time',
 #                 collapse_corr_feats = TRUE, ### perfectly correlated env variables do not add inormation to the model #not the case for my data
 #                 group_neg_corr = TRUE, # none of my variables are in this case
 #                 method =  "center") #if we want to normalize the data 
 
 ##continue without preprocessing
 
 
 # 3. MODEL ENGINEERING -----
 ## In this step I will explore my models, change parammeters, explore the differences to be sure that the model I'm building is correct.
 ### I do not perform the feature importance so that the model runs faster and then I can evaluate which parameters I should use
 ### define a general trainControl for the parammeter tunning with not a lot of cv so that it runs faster 
 control_cv <- trainControl(#method = 'cv', number = 10, 
   savePredictions = 'final', 
   returnResamp = 'final', 
   returnData = T, 
   p = 0.8,
   method = "cv", 
   number = 10#, # repeats = 3 tuneGrid = tibble(mtry = 4)
 )
 
 # FL----
 # Define a list of your model_tb_02_ objects
 bloo_rf_02
model_tb_list <-  bloo_rf_02 |>
   dplyr::select(asv_num) |>
   dplyr::mutate(name = paste0('model_tb_02_', asv_num)) |>
   dplyr::select(name) 
   
concatenated_names <- paste0(model_tb_list)

model_names_clean <- gsub("\"|\\\\", "", concatenated_names )


asv_id_clean 

 model_02_list <- list(model_tb_02_asv38, model_tb_02_asv15, model_tb_02_asv27, model_tb_02_asv17,
                    model_tb_02_asv62, model_tb_02_asv58, model_tb_02_asv1, model_tb_02_asv7, 
                       model_tb_02_asv178, model_tb_02_asv11)
 
 model_tb_list_asv <-  bloo_rf_02 |>
   dplyr::select(asv_num) 
 
 concat_asv_id <-  paste0(model_tb_list_asv)
 asv_id_clean <- gsub("\"|\\\\", "'", concat_asv_id )
 
 # Define a list of corresponding ASV IDs
 asv_id_02_list <- c('asv38', 'asv15', 'asv27', 'asv17', 'asv62', 'asv58', 'asv1', 'asv7', 
                     'asv178', 'asv11')
 
 # Loop through each model and ASV ID
 for (i in seq_along(model_02_list)) {
   # Run the model
   rf <- run_ml(model_list[[i]], method = 'rf',
                outcome_colname = 'diff_rclr_time',
                find_feature_importance = F,
                cross_val = control_cv,
                training_frac = 0.8,
                ntree = 1000,
                perf_metric_name = 'RMSE',
                seed = 1030)
   
   # Assign the result to a variable with a suitable name
   assign(paste0("rf_02_", asv_id_02_list[i]), rf)
 }
 
 # PA----
 # Define a list of your model_tb_3_ objects
 bloo_rf_3
 model_tb_list <-  bloo_rf_3 |>
   dplyr::select(asv_num) |>
   dplyr::mutate(name = paste0('model_tb_3_', asv_num)) |>
   dplyr::select(name) 
 
 concatenated_names <- paste0(model_tb_list)
 
 model_names_clean <- gsub("\"|\\\\", "", concatenated_names )
 
 model_3_list <- list(model_tb_3_asv179, model_tb_3_asv15, model_tb_3_asv72, model_tb_3_asv27, model_tb_3_asv17, 
                    model_tb_3_asv192, model_tb_3_asv84, model_tb_3_asv118, model_tb_3_asv23, model_tb_3_asv85, 
                    model_tb_3_asv25, model_tb_3_asv163, model_tb_3_asv80, model_tb_3_asv116, model_tb_3_asv182, 
                    model_tb_3_asv126, model_tb_3_asv105, model_tb_3_asv28, model_tb_3_asv1, model_tb_3_asv7, 
                    model_tb_3_asv4, model_tb_3_asv31, model_tb_3_asv22, model_tb_3_asv11, model_tb_3_asv42)
 
 model_tb_list_asv <-  bloo_rf_3 |>
   dplyr::select(asv_num) 
 
 concat_asv_id <-  paste0(model_tb_list_asv)
 asv_id_clean <- gsub("\"|\\\\", "'", concat_asv_id )
 
 # Define a list of corresponding ASV IDs
 asv_id_3_list <- c('asv179', 'asv15', 'asv72', 'asv27', 'asv17', 'asv192', 'asv84', 'asv118', 'asv23', 'asv85', 
                  'asv25', 'asv163', 'asv80', 'asv116', 'asv182', 'asv126', 'asv105', 'asv28', 'asv1', 'asv7', 
                  'asv4', 'asv31', 'asv22', 'asv11', 'asv42')
 
 # Loop through each model and ASV ID
 for (i in seq_along(model_list)) {
   # Run the model
   rf <- run_ml(model_list[[i]], method = 'rf',
                outcome_colname = 'diff_rclr_time',
                find_feature_importance = F,
                cross_val = control_cv,
                training_frac = 0.8,
                ntree = 1000,
                perf_metric_name = 'RMSE',
                seed = 1030)
   
   # Assign the result to a variable with a suitable name
   assign(paste0("rf_3_", asv_id_list[i]), rf)
 }
 
 
 # 4. EVALUATION OF THE MODELS -----
 ### I do it 100 different times so that I get different values and I can evaluate the performance. Values from the training and the final model should be similar.
 # Run the process 10 times with different seeds and create tibbles
 num_runs <- 100
 seed_list <- 101:201  # Example list of 100 different seeds
 #seed_list <- 101:121
 
 ## FL---- 
 model_02_list
 asv_id_02_list
 
 names(model_02_list) <- asv_id_02_list
 
 summary_types_of_blooms |>
   dplyr::filter(fraction == '0.2') |>
   distinct(asv_num) |>
   as_vector()

 ##all at the same time
 # Initialize an empty list to store the combined results for all ASVs
 combined_rf_02_performance_all <- list()
 
 # Loop through each ASV in the model_tb_list
 for (asv in names(model_02_list)) {
   # Initialize a list to store the tibbles for each ASV
   rf_02_performance_list <- list()
   
   # Loop through each run
   for (i in 1:num_runs) {
     # Run ml and create tibble for the current ASV and run
     rf_02_performance <- run_ml_and_create_tibble(model_02_list[[asv]], asv, seed_list[i])
     
     # Store the tibble in the list with a unique name
     rf_02_performance_list[[paste0("rf_02_performance_", i)]] <- rf_02_performance
   }
   
   # Bind all tibbles for the current ASV into one
   combined_rf_02_performance <- bind_rows(rf_02_performance_list)
   
   # Store the combined tibble for the current ASV in the list
   combined_rf_02_performance_all[[asv]] <- combined_rf_02_performance
 }
 
 # Combine all tibbles from all ASVs into one
 combined_rf_02_performance_all <- do.call(bind_rows, combined_rf_02_performance_all)
 
 combined_rf_02_performance_all |>
   distinct(asv_num)
 
 summary_types_of_blooms |>
   dplyr::filter(fraction == '0.2') |>
   distinct(asv_num)
 
 ## add some information about my bloomers before ploting (FL)-----
 combined_rf_02_performance_all <- combined_rf_02_performance_all |>
   dplyr::mutate(fraction =  '0.2') |>
   left_join(summary_types_of_blooms, by = c('fraction', 'asv_num'))
 
 ## PA-----
 model_3_list
 asv_id_3_list
 
 names(model_3_list) <- asv_id_3_list
 
 ##all at the same time
 # Initialize an empty list to store the combined results for all ASVs
 combined_rf_3_performance_all <- list()
 
 # Loop through each ASV in the model_tb_list
 for (asv in names(model_3_list)) {
   # Initialize a list to store the tibbles for each ASV
   rf_3_performance_list <- list()
   
   # Loop through each run
   for (i in 1:num_runs) {
     # Run ml and create tibble for the current ASV and run
     rf_3_performance <- run_ml_and_create_tibble(model_3_list[[asv]], asv, seed_list[i])
     
     # Store the tibble in the list with a unique name
     rf_3_performance_list[[paste0("rf_3_performance_", i)]] <- rf_3_performance
   }
   
   # Bind all tibbles for the current ASV into one
   combined_rf_3_performance <- bind_rows(rf_3_performance_list)
   
   # Store the combined tibble for the current ASV in the list
   combined_rf_3_performance_all[[asv]] <- combined_rf_3_performance
 }
 
 # Combine all tibbles from all ASVs into one
 combined_rf_3_performance_all <- do.call(bind_rows, combined_rf_3_performance_all)
 
 ## add some information about my bloomers before plotting (FL)-----
 combined_rf_3_performance_all <- combined_rf_3_performance_all |>
   dplyr::mutate(fraction =  '3') |>
   left_join(summary_types_of_blooms)
 
 combined_rf_3_performance_all$recurrency <- factor(combined_rf_3_performance_all$recurrency, levels = c('no', 'yes'))
 
 ##plot RMSE and cv RMSE -----
 labs_recurrent <- as_labeller(c('no' = 'Non-recurrent', 'yes' = 'Recurrent'))
 
 RMSE_cv <- combined_rf_02_performance_all |>
   pivot_longer(cols = contains('RMSE'), values_to = 'value', names_to = 'type') %>%
   ggplot(aes(x = asv_num, y = value, fill = type)) +
   geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 0.2, color = 'black') +
   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
   facet_wrap(vars(recurrency), scales = 'free_x', labeller = labs_recurrent)+
   scale_fill_manual(values = c('darkblue', 'grey')) +
   labs(x = "ASV Number", y = "RMSE", fill = 'Type', shape = 'Reccurrent taxa') +
   theme(strip.background = element_rect(fill = 'transparent'),
         panel.grid = element_blank(),
         legend.position = 'bottom',
         text = element_text(size = 6),
         plot.margin = margin(2,2,2,8))
 
 RMSE_cv 
 
 ggsave(RMSE_cv, filename = 'RMSE_cv_fl.pdf',
        path = 'Results/Figures/random_forest/',
        width = 188, height = 100, units = 'mm')
 
 RMSE_cv <- combined_rf_3_performance_all |>
   pivot_longer(cols = contains('RMSE'), values_to = 'value', names_to = 'type') %>%
   ggplot(aes(x = asv_num, y = value, fill = type)) +
   geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 0.2, color = 'black') +
   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
   facet_wrap(vars(recurrency), scales = 'free_x', labeller = labs_recurrent)+
   scale_fill_manual(values = c('darkblue', 'grey')) +
   labs(x = "ASV Number", y = "RMSE", fill = 'Type', shape = 'Reccurrent taxa') +
   theme(strip.background = element_rect(fill = 'transparent'),
         panel.grid = element_blank(),
         legend.position = 'bottom',
         text = element_text(size = 6),
         plot.margin = margin(2,2,2,8))
 
 RMSE_cv 
 
 # ggsave(RMSE_cv, filename = 'RMSE_cv_pa.pdf',
 #        path = 'Results/Figures/random_forest/',
 #        width = 188, height = 100, units = 'mm')
 
 ##plot R2 and resample R2-----
 Rsquared_cv <- combined_rf_02_performance_all |>
   pivot_longer(cols = contains('Rsquared'), values_to = 'value', names_to = 'type') |>
   ggplot( aes(x = asv_num, y = value, fill = type)) +
   geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 0.2, color = 'black') +
   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
   scale_fill_manual(values = c('darkblue', 'grey')) +
   scale_y_continuous(limits = c(0,1))+
   facet_wrap(vars(recurrency), scales = 'free_x', labeller = labs_recurrent)+
   #scale_color_manual(values = c('darkblue', 'grey'))+
   labs(x = "ASV Number", y = expression('R'^2), fill = 'Type', shape = 'Reccurrent taxa')+
   theme(strip.background = element_rect(fill = 'transparent'),
         panel.grid = element_blank(),
         text = element_text(size = 6),
         legend.position = 'bottom',
         plot.margin = margin(2,2,2,8))
 
 Rsquared_cv 
 
 # ggsave(Rsquared_cv, filename = 'Rsquared_cv_fl.pdf',
 #        path = 'Results/Figures/random_forest/',
 #        width = 188, height = 100, units = 'mm')
 
 Rsquared_cv <- combined_rf_3_performance_all |>
   pivot_longer(cols = contains('Rsquared'), values_to = 'value', names_to = 'type') |>
   ggplot( aes(x = asv_num, y = value, fill = type)) +
   geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 0.2, color = 'black') +
   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
   scale_fill_manual(values = c('darkblue', 'grey')) +
   scale_y_continuous(limits = c(0,1))+
   facet_wrap(vars(recurrency), scales = 'free_x', labeller = labs_recurrent)+
   #scale_color_manual(values = c('darkblue', 'grey'))+
   labs(x = "ASV Number", y = expression('R'^2), fill = 'Type', shape = 'Reccurrent taxa')+
   theme(strip.background = element_rect(fill = 'transparent'),
         panel.grid = element_blank(),
         text = element_text(size = 6),
         legend.position = 'bottom',
         plot.margin = margin(2,2,2,8))
 
Rsquared_cv 

 # ggsave(Rsquared_cv, filename = 'Rsquared_cv_pa.pdf',
 #        path = 'Results/Figures/random_forest/',
 #        width = 188, height = 100, units = 'mm')
 
 
 
 
 
 
 
 # 5. FEATURE IMPORTANCE ----
 ### now that I have evaluated my models I compute the feature importance. 
 ### Understanding the results feature_importance: If feature importances were calculated, a data frame 
 ### where each row is a feature or correlated group. The columns are the performance metric of the 
 ### permuted data, the difference between the true performance metric and the performance metric of the 
 ### permuted data (true - permuted), the feature name, the ML method, the performance metric name, and 
 ### the seed (if provided). For AUC and RMSE, the higher perf_metric_diff is, the more important that 
 ### feature is for predicting the outcome. For log loss, the lower perf_metric_diff is, the more 
 ### important that feature is for predicting the outcome.
 # 
 # There are several columns:
 # perf_metric: The performance value of the permuted feature.
 # perf_metric_diff: The difference between the performance for the actual and permuted data (i.e. test performance minus permuted performance). Features with a larger perf_metric_diff are more important.
 # pvalue: the probability of obtaining the actual performance value under the null hypothesis.
 # lower: the lower bound for the 95% confidence interval of perf_metric.
 # upper: the upper bound for the 95% confidence interval of perf_metric.
 # feat: The feature (or group of correlated features) that was permuted.
 # method: The ML method used.
 # perf_metric_name: The name of the performance metric represented by perf_metric & perf_metric_diff.
 # seed: The seed (if set).
 
 # Run the process 10 times with different seeds and create tibbles
 num_runs <- 10
 seed_list <- 101:111  # Example list of 10 different seeds
 
 ## FL ----
 model_02_list
 asv_id_02_list
 
 model_02_list |>
   names()

 rf_importance_list <- list()
 
 ## run with fixed seed and the we will see if I can compute it for different seeds
 
 # Loop through each run
 for (i in seq_along(asv_id_02_list)) {
   # Run ml and create tibble
   rf_importance <- run_ml_importance_and_create_tibble(model_02_list[[i]], asv_id_02_list[i], seed = 2002)
   
   # Store the tibble in the list with a unique name
   rf_importance_list[[paste0("rf_importance_", i)]] <- rf_importance
 }
 
 # Combine all tibbles into one
 combined_rf_importance <- bind_rows(rf_importance_list)
 
 # Combine all tibbles from all ASVs into one
 combined_rf_02_importance_all <- do.call(bind_rows, combined_rf_importance)
 
 combined_rf_02_importance_all |>
   distinct(asv_num)
 
 ## PA -----
 model_3_list
 asv_id_3_list
 
 model_3_list |>
   names()
 
 rf_importance_list <- list()
 
 ## run with fixed seed and the we will see if I can compute it for different seeds
 
 # Loop through each run
 for (i in seq_along(asv_id_3_list)) {
   # Run ml and create tibble
   rf_importance <- run_ml_importance_and_create_tibble(model_3_list[[i]], asv_id_3_list[i], seed = 203)
   
   # Store the tibble in the list with a unique name
   rf_importance_list[[paste0("rf_importance_", i)]] <- rf_importance
 }
 
 # Combine all tibbles into one
 combined_rf_importance <- bind_rows(rf_importance_list)
 
 # Combine all tibbles from all ASVs into one
 combined_rf_3_importance_all <- do.call(bind_rows, combined_rf_importance) |>
   dplyr::mutate(fraction = '3')
 
 ## add some information about my bloomers before ploting
 combined_rf_importance_all <- combined_rf_02_importance_all |>
   dplyr::mutate(fraction = '0.2') |>
   bind_rows(combined_rf_3_importance_all) |>
   left_join(summary_types_of_blooms)  
 
 combined_rf_importance_all |>
   group_by(asv_num, fraction) |>
   distinct(asv_num, fraction) |>
   group_by(fraction) |>
   reframe(n = n())
 
 ##plot
 combined_rf_importance_all |>
   dplyr::mutate(feat = as.factor(feat)) |>
   # dplyr::group_by(asv_num, seed) |>
   # slice_max(order_by = abs(perf_metric_diff), n = 10) |>
   ggplot(aes(fct_infreq(feat), perf_metric_diff, shape = recurrency))+
   scale_x_discrete(labels = labs_env_models_no_nas)+
   geom_point()+
   geom_boxplot(alpha = 0.6)+
   #facet_wrap(vars(fraction), labeller = labs_fraction, scales = 'free_x')+
   labs(y = 'Perf metric difference', x = 'Explanatory variables')+
   facet_grid(asv_num~fraction, scales = 'free_y')+
   coord_flip()+
   theme(panel.grid = element_blank(), strip.background = element_rect(fill = 'transparent'))
 
 combined_rf_importance_all |>
   dplyr::filter(fraction == '0.2' ) |>
   dplyr::mutate(feat = as.factor(feat)) |>
   # dplyr::group_by(asv_num, seed) |>
   # slice_max(order_by = abs(perf_metric_diff), n = 5) |>
   ggplot(aes(fct_infreq(feat), perf_metric, shape = recurrency))+
   scale_x_discrete(labels = labs_env_models_no_nas)+
   geom_point()+
   geom_boxplot(alpha = 0.6)+
   labs(y = 'Perf metric difference', x = 'Explanatory variables')+
   facet_wrap(vars(asv_num), scales = 'free_y')+
   coord_flip()+
   theme(panel.grid = element_blank(), strip.background = element_rect(fill = 'transparent'))
 
 
 
 
 
 # 6. PARTIAL DEPENDANCE-----
 ### we explore how do the variables interact with each other
 ### https://christophm.github.io/interpretable-ml-book/pdp.html
 ## One important thing of this analysis is to see the direction of the relationship.
 ## Ecological significance of the results based on the previous knowledge of this taxa.
 
 #### In this case we use the whole dataset to create the model, because we want to gain more confidence on the inferred patterns.
 
 ### Save pdp data for plotting
 outcome_colname <- "diff_rclr_time"
 
 predictors <- c("bacteria_joint", "synechococcus", "temperature_no_nas", "day_length_no_nas",  "chla_total_no_nas",  "PO4_no_nas" ,       
                 "NH4_no_nas" ,  "NO2_no_nas"  , "NO3_no_nas" , "Si_no_nas" , "PNF_Micro_no_nas",  "cryptomonas_no_nas",
                 "micromonas_no_nas" , "HNF_Micro_no_nas"  , "fitted_rclr" )
 
 seed_list <- 1:10
 
 asv_id_02_list
 model_02_list
 
 for (asv_num in asv_id_02_list) {
   train_rf_and_save_pdp(data =  model_02_list[[asv_num]], 
                         outcome_colname = outcome_colname, 
                         predictors = predictors, 
                         seed_list = seed_list, 
                         grid_resolution = 119, 
                         asv_num = asv_num)
 }
 
 ##read RData----
 # Initialize a list to store the data
 pdp_data_02_list <- list()
 
 # Loop through each asv_num
 for (asv_num in asv_id_02_list) {
   # Construct the file path
   file_path <- paste0("results/figures/random_forest/pdp/pdp_", asv_num, "_data.RDS")
   
   # Read the RDS file
   pdp_data <- readRDS(file_path)
   
   # Store the data in the list with a unique name
   pdp_data_02_list[[asv_num]] <- pdp_data
 } 
 
 ### 2. compute the partial dependence for each seed in the same plot. 
 
 ## extraction of the importance variables result to order the partial dependence plots by importance----
 combined_rf_02_importance_all |>
   distinct(asv_num) #check that I have the information for all my ASVs
 
 importance_order <- combined_rf_02_importance_all |>
   dplyr::group_by(asv_num, feat) |>
   dplyr::reframe(mean_diff = mean(perf_metric_diff)) |>
   dplyr::arrange(abs(mean_diff)) 
 
 # Call the function with the pdp data list----
 plot_partial_dependence(pdp_asv62_data, importance_df =  importance_order, asv_num = 'asv62', num_plots = 3)
 plot_partial_dependence(pdp_asv11_data, importance_df =  importance_order, asv_num = 'asv11', num_plots = 3)
 plot_partial_dependence(pdp_asv72_data, importance_df =  importance_order, asv_num = 'asv72', num_plots = 3)
 plot_partial_dependence(pdp_asv17_data, importance_df =  importance_order, asv_num = 'asv17', num_plots = 3)
 plot_partial_dependence(pdp_asv23_data, importance_df =  importance_order, asv_num = 'asv23', num_plots = 3)
 
 for (asv_num in asv_id_02_list) {

     data <- pdp_data_02_list[[asv_num]]
     #asv_num <- paste0("'", [[asv_num]], "'")
  
   plot_partial_dependence(data = , importance_df =  importance_order, asv_num = 'asv62', num_plots = 3)
   
 }
 
 plot_partial_dependence(pdp_asv38_data, importance_df =  importance_order, asv_num = 'asv38', num_plots = 3)
 
 # Composition with all the top 3 effects on each blooming example----
 pdp_bloomers_02 <- gridExtra::grid.arrange(pdp_all_plots_grid_asv62,
                                                 pdp_all_plots_grid_asv72,
                                                 pdp_all_plots_grid_asv23, 
                                                 pdp_all_plots_grid_asv17, 
                                                 pdp_all_plots_grid_asv11, 
                                                 ncol = 1)
 
 # ggsave(pdp_bloomers_02, filename = 'pdp_bloomers_02.pdf',
 #        path = 'Results/Figures/random_forest/pdp/',
 #        width = 188, height = 220, units = 'mm')
 
 
 
 
 
##### WHAT HAPPENS WITH OTHER BLOOMERS (DESIGN A CODE TO DO IT FOR ALL OF THEM) ####-----
# bloo_02$value #discard SAR11 clade "asv8"   "asv5"   "asv3"   "asv2"
# bloo_3$value
# 
# ## ASV1 (in PA and FL)
# #### PA----
# #fit the loess model
# fit <- loess(asv1 ~ decimal_date, data = rclr_time_tb_3, span = 0.084)
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
# 
# asv1_sm <- predicted_values |>
#   as_tibble_col(column_name = 'asv1') |>
#   bind_cols(decimal_date_tb) |>
#   dplyr::select(decimal_date, 'fitted_rclr' = 'asv1') |>
#   dplyr::mutate(asv_num = 'asv1') |>
#   dplyr::left_join(fraction_date_3)
# 
# rclr_time_tb_3 |>
#   dplyr::select(decimal_date, 'rclr' = 'asv1') |>
#   left_join(asv1_sm) |>
#   pivot_longer(cols = c('rclr', 'fitted_rclr'), names_to = 'approx', values_to = 'abundance_value') |>
#   ggplot(aes(decimal_date, abundance_value))+
#   geom_line(aes(group = approx, color = approx))+
#   scale_color_manual(values = c('black', 'grey'))+
#   theme_bw()+
#   theme(panel.grid = element_blank())
# 
# 
# 
# ### BLOOMER asv1----
# model_tb <- asv1_sm |>
#   dplyr::filter(!is.na(diff_time)) |>
#   #dplyr::select(-sample_id_num) |>
#   dplyr::filter(asv_num == 'asv1' &
#                   fraction == '3') |>
#   ungroup() |>
#   dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
#   dplyr::select(-contains('_no_nas')) |>
#   left_join(env_sm_diff_w) |>
#   #left_join(env_data_new) |>
#   dplyr::select(-decimal_date)
# 
# ##add previous CLR value as explanatory variable
# asvs_sm_asv <- asvs_sm |>
#   dplyr::filter(asv_num == 'asv1' &
#                   fraction == '0.2') |>
#   dplyr::arrange(decimal_date) |>
#   dplyr::select(-asv_num, -fraction, -decimal_date)
# 
# asvs_sm_asv <- asvs_sm_asv[-120,] #the last abundance won't be necessary since we do not need it to predict the next one
# 
# model_tb <- model_tb |>
#   bind_cols(asvs_sm_asv)|>
#   dplyr::select(-BP_FC1.55_sm) ##BP from the previous month should not have a high impact on the change in abundance of the next
# 
# ### Split data to train and test
# indices <- createDataPartition(y=model_tb$diff_rclr_time, p =0.75, list = F)
# train_tb <- model_tb[indices,]
# test_tb <- model_tb[-indices,]
# 
# ### Train RF
# control <- trainControl(method = 'cv', number = 20)
# bloo_rf <- train( diff_rclr_time ~., method='rf', trControl=control, data=train_tb, metric = 'RMSE', ntree=300, importance=T)
# y <- predict(bloo_rf, newdata = test_tb[,-1])
# test_tb$prediction <- y
# fit <- lm(prediction~diff_rclr_time, data=test_tb)
# summary(fit)
# importance_df <- as.data.frame(bloo_rf$finalModel$importance)
# importance_df <- importance_df[order(-importance_df$`%IncMSE`),]
# 
# ##try to define best parameters----
# # Define train control
# control <- trainControl(method = 'cv', number = 20)
# 
# # Define the tuning parameter grid
# sqrt(ncol(train_tb)) ## the mtry should we a number arround this one
# 
# tunegrid <- expand.grid(
#   mtry = c(4:6) #, # 
#   #ntree = c(100, 200)
# )
# 
# ## the ntree can not be changed using tunegrid that is why i need to loop over different numbers of ntree
# modellist <- list()
# for (i in 1:nrow(tunegrid)) {
#   set.seed(100)
#   mtry_val <- tunegrid$mtry[i]
#   
#   for (ntree_val in c(200, 300, 500, 100)) { # Iterate over ntree values
#     fit <- train(diff_rclr_time ~ ., 
#                  data = train_tb, 
#                  method = 'rf', 
#                  metric = 'RMSE', 
#                  tuneGrid = tunegrid,
#                  trControl = control)
#     
#     key <- paste("mtry_", mtry_val, "_ntree_", ntree_val, sep = "")
#     modellist[[key]] <- fit
#   }
# }
# 
# # compare results
# results <- resamples(modellist)
# summary(results)
# dotplot(results) ## ntree = 500
# 
# ## retrain model on entire dataset----
# bloo_final_asv1_diff <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
# bloo_final_asv1_diff
# 
# importance_df_asv1_diff <- as.data.frame(bloo_final_asv1_diff$finalModel$importance)
# importance_df_asv1_diff <- importance_df_asv1_diff[order(-importance_df_asv1_diff$`%IncMSE`),]
# importance_df_asv1_diff