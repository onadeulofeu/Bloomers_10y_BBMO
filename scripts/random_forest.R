#packages----
library(caret)
library(forcats)

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
  left_join(env_data_interpolated_values_all_red) ## add metadata (not smoothed) 

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
## maybe I need to smooth my variables first. 
diff_time_tb_env |>
  ggplot(aes(sample_id_num, diff_time))+
  geom_line()+
  geom_smooth(method = 'loess', span = 0.02)

model_tb <- diff_time_tb_env |>
  dplyr::filter(!is.na(diff_time)) |>
  dplyr::select(-sample_id_num) |>
  dplyr::filter(asv_num_f == 'asv62' &
                  fraction == '0.2') |>
  ungroup() |>
  dplyr::select(-asv_num_f, -decimal_date, -diff_rclr, -diff_time, -date, -fraction, -abundance_value)

model_tb |>
  colnames()

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

## remove some variables that i believe are not explanatory
model_tb <- model_tb |>
  dplyr::select( -bacteria_joint, -Si_no_nas, -BP_FC1.55_no_nas)

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
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final

## Linear model
linear_fit <- lm(diff_rclr_time ~., data=model_tb)
summary(linear_fit)

## To increase the performance of our model we could smooth our environmental variables----
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
asv_tab_bbmo_10y_w_rar_3_inter <- rrarefy.perm(round(asv_tab_bbmo_10y_w_3_inter),
                                       sample = min_n_seqs,
                                       n = 1000,
                                       round.out = T)

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
fit <- loess(day_length_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.1)

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
fit <- loess( synechococcus ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.15)

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
fit <- loess( bacteria_joint ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.2)

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
# fit <- loess(temperature_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.09)
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
fit <- loess(chla_total_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.2)

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
fit <- loess(PO4_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.2)

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
fit <- loess(NH4_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.2)

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
fit <- loess(NO2_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.2)

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
fit <- loess(NO3_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.15)

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
fit <- loess(Si_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.15)

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
fit <- loess(BP_FC1.55_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.2)

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
fit <- loess(PNF_Micro_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.15)

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
fit <- loess(cryptomonas_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.2)

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
fit <- loess(micromonas_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.15)

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
fit <- loess(HNF_Micro_no_nas ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.18)

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
fit <- loess( community_eveness_rar_02 ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.18)

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
fit <- loess( community_eveness_rar_3 ~ decimal_date, data = diff_rclr_time_tb_env, span = 0.18)

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

env_data_new_l |>
  bind_rows(diff_rclr_time_tb_env_l) |>
  ggplot(aes(decimal_date, value))+
  geom_point(aes(shape = type_value))+
  geom_line(aes(color = type_value))+
  facet_wrap(vars(env_variable), scales = 'free')+
  theme_bw()

## Model with some smoothed env variables-----

#prepare the env data
env_data_new_red <- env_data_new |>
  dplyr::select(-decimal_date)

env_data_new_red <- env_data_new_red[-120,]

### BLOOMER ASV62----
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

### BLOOMER ASV11-----
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

### BLOOMER ASV23-----
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

### BLOOMER asv72-----
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

rclr_time_tb_3 <- diff_time_tb_env |>
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
fit <- loess(asv23 ~ decimal_date, data = rclr_time_tb_3, span = 0.15)

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
  geom_line(aes(group = approx))

#fit the loess model
fit <- loess(asv72 ~ decimal_date, data = rclr_time_tb_3, span = 0.15)

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
  geom_line(aes(group = approx))

#fit the loess model
fit <- loess(asv17 ~ decimal_date, data = rclr_time_tb_3, span = 0.15)

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
  geom_line(aes(group = approx))

### 02 (asv62, asv11, asv555)----

#fit the loess model
fit <- loess(asv62 ~ decimal_date, data = rclr_time_tb_02, span = 0.15)

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
  geom_line(aes(group = approx))

#fit the loess model
fit <- loess(asv11 ~ decimal_date, data = rclr_time_tb_02, span = 0.15)

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
  geom_line(aes(group = approx))

#fit the loess model
fit <- loess(asv555 ~ decimal_date, data = rclr_time_tb_02, span = 0.15)

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
  geom_line(aes(group = approx))

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
  facet_wrap(vars(asv_num))

asvs_sm_diff |>
  ggplot(aes(decimal_date, fitted_rclr))+
  geom_line()+
  facet_wrap(vars(asv_num))

## blooms and rate of change from one point to the next one
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
set.seed(100)

### BLOOMER ASV62 (smooth) ----
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
dotplot(results)

# View summary of the trained model
summary(custom)

# Plot the model
plot(custom)

##in this case is ntree 200

## retrain model on entire dataset----
bloo_final_asv62 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
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
dotplot(results)

##500 is best

## retrain model on entire dataset----
bloo_final_asv23 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv23

importance_df_asv23 <- as.data.frame(bloo_final_asv23$finalModel$importance)
importance_df_asv23 <- importance_df_asv23[order(-importance_df_asv23$`%IncMSE`),]

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
dotplot(results)

## ntree = 300

## retrain model on entire dataset----
bloo_final_asv72 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=300, importance=T)
bloo_final_asv72

importance_df_asv72 <- as.data.frame(bloo_final_asv72$finalModel$importance)
importance_df_asv72 <- importance_df_asv72[order(-importance_df_asv72$`%IncMSE`),]

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
dotplot(results)

##ntree = 500

## retrain model on entire dataset----
bloo_final_asv17 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final_asv17

importance_df_asv17 <- as.data.frame(bloo_final_asv17$finalModel$importance)
importance_df_asv17 <- importance_df_asv72[order(-importance_df_asv17$`%IncMSE`),]

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
dotplot(results)

## ntree = 200

## retrain model on entire dataset----
bloo_final_asv11 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final_asv11

importance_df_asv11 <- as.data.frame(bloo_final_asv11$finalModel$importance)
importance_df_asv11 <- importance_df_asv11[order(-importance_df_asv11$`%IncMSE`),]

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
dotplot(results)

# ntree = 100 pero poso 200 que surt millor

## retrain model on entire dataset----
bloo_final_asv555 <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=200, importance=T)
bloo_final_asv555

importance_df_asv555 <- as.data.frame(bloo_final_asv555$finalModel$importance)
importance_df_asv555 <- importance_df_asv555[order(-importance_df_asv555$`%IncMSE`),]

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

im_all |>
  ggplot(aes(fct_infreq(variable), incMSE))+
  geom_col()+
  facet_wrap(~ asv_num)+
  labs(y = '%IncMSE')+
  theme_bw()+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank())

## check if the predicted values in the 

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

## PARTIAL DEPENDANCE -----
### we explore how do the variables interact with each other
## we need the pdp package
library(pdp)

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
    theme(panel.grid = element_blank())
}

# Display all partial dependence plots
pdp_asv62 <- gridExtra::grid.arrange(grobs = partial_plots)

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
    theme(panel.grid = element_blank())
}

# Display all partial dependence plots
pdp_asv23 <- gridExtra::grid.arrange(grobs = partial_plots)

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
    theme(panel.grid = element_blank())
}

# Display all partial dependence plots
pdp_asv72 <- gridExtra::grid.arrange(grobs = partial_plots)

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
    theme(panel.grid = element_blank())
}

# Display all partial dependence plots
pdp_asv17 <- gridExtra::grid.arrange(grobs = partial_plots)

#### ASV11----
importance_df_asv11 |>
  arrange('%IncMSE')

importance_df_asv555 |>
  arrange('%IncMSE')

importance_df_asv62 |>
  arrange('%IncMSE')

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
    theme(panel.grid = element_blank())
}

# Display all partial dependence plots
pdp_asv11 <- gridExtra::grid.arrange(grobs = partial_plots)


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
    theme(panel.grid = element_blank())
}


# Display all partial dependence plots
pdp_asv555 <- gridExtra::grid.arrange(grobs = partial_plots)







# PENSAR SI TINDRIA SENTIT O NO-----
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

##in this case is ntree 500

## retrain model on entire dataset----
bloo_final <- train( diff_rclr_time ~., method='rf', data=model_tb, metric = 'RMSE', ntree=500, importance=T)
bloo_final

importance_df_asv62_diff <- as.data.frame(bloo_final$finalModel$importance)
importance_df_asv62_diff <- importance_df[order(-importance_df$`%IncMSE`),]


