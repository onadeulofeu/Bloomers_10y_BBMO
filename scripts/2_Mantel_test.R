# Analyse if blooming events respond to changes in environmental conditions
## Mantel tests are correlation tests that determine the correlation between
## two matrices (rather than two variables). When using the test for microbial
## ecology, the matrices are often distance/dissimilarity matrices
## with corresponding positions. In order to calculate the correlation, the matrix
## values of both matrices are 'unfolded' into long column vectors, which are then
## used to determine correlation. Permutations of one matrix are used to determine
## significance.

## Mantel's test is a regression in which the variables are themselves distance or
## dissimilarity matrices summarizing pairwise similarities among sample locations.

## A significant Mantel test will tell you that the distances between samples
## in one matrix are correlated with the distances between samples in the other
## matrix. 

## r falls in the range of -1 to +1, where being close to -1 indicates strong 
## negative correlation and +1 indicates strong positive correlation. 
## iouiAn r value of 0 indicates no correlation.

## We need species abundance dissimilairy matrix and environmental parameter
## distance matrix: created using Euclidean Distance.

## Unless using the ranked Mantel statistic, the Mantel approach is suited to 
## detect linear relationships between (dis)similarity matrices.

library(vegan)
library(tidyverse)
library(ggplot2)

asv_tab_all_bloo_z_tax <- read.csv2('data/asv_tab_all_bloo_z_tax_new.csv')
asv_tab_rar <- read.csv2('data/asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()

#I create two different datasets one for ASVs and the other for the community------
asv_tab_rar |>
  colnames()

asv_tab_rar_l <- asv_tab_rar |>
  pivot_longer(starts_with('asv'), 
               values_to = 'abundance', 
               names_to = 'asv_num')

asv_tab_all_bloo_z_tax |>
  colnames()

sample_id <- asv_tab_all_bloo_z_tax$sample_id

bbmo_env <- asv_tab_all_bloo_z_tax |>
  dplyr::select(sample_id,
                day_length,
                #sampling_time
                temperature,
                secchi,
                salinity,
                chla_total,
                chla_3um,
                PO4,
                NH4, NO2, NO3,
                Si, BP_FC1.55,
                PNF_Micro, PNF2_5um_Micro,
                PNF_5um_Micro, 
                dryptomonas, micromonas,
                HNF_Micro, HNF2_5um_Micro,   
                HNF_5um_Micro ,  LNA,              
                HNA,   prochlorococcus_FC,
                Peuk1,   Peuk2,                 
                bacteria_joint, synechococcus) |>
  dplyr::select(-sample_id) |>
  dplyr::mutate_if( is.character, as.numeric) |>
  bind_cols(sample_id) |>
  rename(sample_id = '...28') |>
  distinct(sample_id,
           day_length,
           #sampling_time
           temperature,
           secchi,
           salinity,
           chla_total,
           chla_3um,
           PO4,
           NH4, NO2, NO3,
           Si, BP_FC1.55,
           PNF_Micro, PNF2_5um_Micro,
           PNF_5um_Micro, 
           dryptomonas, micromonas,
           HNF_Micro, HNF2_5um_Micro,   
           HNF_5um_Micro ,  LNA,              
           HNA,   prochlorococcus_FC,
           Peuk1,   Peuk2,                 
           bacteria_joint, synechococcus)

bbmo_env_02 <- bbmo_env |>
  dplyr::filter(str_detect(sample_id, '0.2'))

bbmo_env_3 <- bbmo_env |>
  dplyr::filter(str_detect(sample_id, '3_'))

##normalization of environmental data using z-scores

class(bbmo_env_sim$day_length)
sum(is.na(bbmo_env_sim$day_length))

bbmo_env_sim <- bbmo_env |>
  dplyr::select(sample_id, day_length) |>
  dplyr::filter(!is.na(day_length)) |>
  dplyr::mutate(day_length = as.numeric(day_length)) |>
  ungroup()

## I created a function for this purpose 
# calculate_z_score <- function(data, col, name = NULL, group = NULL) {
#   stopifnot(is.numeric(data[[col]]))
#   
#   # Check for NAs in the specified column
#   if (anyNA(data[[col]])) {
#     warning("The specified column contains NA values.")
#   }
#   
#   col_name <- ifelse(!is.null(name), paste0("z_score_", name), 'z_score')
#   
#   if (!is.null(group)) {
#     data <- data |>
#       dplyr::group_by({{group}}) |>
#       dplyr::mutate(!!col_name := (!!{{col}} - base::mean(!!{{col}}, na.rm = TRUE)) / base::sd(!!{{col}}, na.rm = TRUE))
#   } else {
#     data <- data |>
#       dplyr::mutate(!!col_name := (!!{{col}} - base::mean(!!{{col}}, na.rm = TRUE)) / base::sd(!!{{col}}, na.rm = TRUE))
#   }
#   
#   return(data)
# }

calculate_z_score <- function(data, col, name = NULL, group = NULL) {
  stopifnot(is.numeric(data[[col]]))
  
  # Check for NAs in the specified column
  if (anyNA(data[[col]])) {
    warning("The specified column contains NA values.")
  }
  
  col_name <- ifelse(!is.null(name), paste0("z_score_", name), 'z_score')
  
  if (!is.null(group)) {
    data <- data |>
      dplyr::group_by(!!sym(group)) |>
      dplyr::mutate(!!col_name := (!!sym(col) - base::mean(!!sym(col), na.rm = TRUE)) / stats::sd(!!sym(col), na.rm = TRUE))
  } else {
    data <- data |>
      dplyr::mutate(!!col_name := (!!sym(col) - base::mean(!!sym(col), na.rm = TRUE)) / stats::sd(!!sym(col), na.rm = TRUE))
  }
  
  return(data)
}

calculate_z_score(bbmo_env_sim, col = 'day_length', name = 'day_length', group = NULL)

bbmo_env_sim |>
  glimpse()

bbmo_env_z <- 
  bbmo_env |>
  as_tibble() |>
  pivot_longer(cols = c( day_length,
                         #sampling_time
                         temperature,
                         secchi,
                         salinity,
                         chla_total,
                         chla_3um,
                         PO4,
                         NH4, NO2, NO3,
                         Si, BP_FC1.55,
                         PNF_Micro, PNF2_5um_Micro,
                         PNF_5um_Micro, 
                         dryptomonas, micromonas,
                         HNF_Micro, HNF2_5um_Micro,   
                         HNF_5um_Micro ,  LNA,              
                         HNA,   prochlorococcus_FC,
                         Peuk1,   Peuk2,                 
                         bacteria_joint, synechococcus), values_to = 'env_values', names_to = 'environmental_variable') |>
    dplyr::filter(!is.na(env_values)) |>
    dplyr::mutate(env_values = as.numeric(env_values)) |>
  calculate_z_score(col = 'env_values', name = 'environmental_variable', group = NULL) 
  #pivot_wider(id_cols = )
  
bbmo_env_z$samples_id |>
  unique()
  
# ## calculate z-scores without using a function---- 
#   bbmo_env_zscore_02  <- bbmo_env_02 |>
#     as_tibble() |>
#     pivot_longer(cols = c(day_length,
#                            #sampling_time
#                            temperature,
#                            secchi,
#                            salinity,
#                            chla_total,
#                            chla_3um,
#                            PO4,
#                            NH4, NO2, NO3,
#                            Si, BP_FC1.55,
#                            PNF_Micro, PNF2_5um_Micro,
#                            PNF_5um_Micro, 
#                            dryptomonas, micromonas,
#                            HNF_Micro, HNF2_5um_Micro,   
#                            HNF_5um_Micro ,  LNA,              
#                            HNA,   prochlorococcus_FC,
#                            Peuk1,   Peuk2,                 
#                            bacteria_joint, synechococcus), values_to = 'env_values', names_to = 'environmental_variable') |>
#     #group_by(environmental_variable) |>
#     dplyr::filter(!is.na(env_values)) |>
#     calculate_z_score(col = 'env_values', name = 'environmental_variable', group = 'environmental_variable') 
#     
#     # dplyr::mutate(sd = sd(env_values),
#     #               mean = mean(env_values)) |>
#     # dplyr::mutate(z_score = ((env_values - mean(env_values))/ sd(env_values))) |>
#     # ungroup()
#   
#   bbmo_env_zscore_3  <- bbmo_env_3 |>
#     as_tibble() |>
#     pivot_longer(cols = c(day_length,
#                           #sampling_time
#                           temperature,
#                           secchi,
#                           salinity,
#                           chla_total,
#                           chla_3um,
#                           PO4,
#                           NH4, NO2, NO3,
#                           Si, BP_FC1.55,
#                           PNF_Micro, PNF2_5um_Micro,
#                           PNF_5um_Micro, 
#                           dryptomonas, micromonas,
#                           HNF_Micro, HNF2_5um_Micro,   
#                           HNF_5um_Micro ,  LNA,              
#                           HNA,   prochlorococcus_FC,
#                           Peuk1,   Peuk2,                 
#                           bacteria_joint, synechococcus), values_to = 'env_values', names_to = 'environmental_variable') |>
#     group_by(environmental_variable) |>
#     dplyr::filter(!is.na(env_values)) |>
#     dplyr::mutate(sd = sd(env_values),
#                   mean = mean(env_values)) |>
#     dplyr::mutate(z_score = ((env_values - mean(env_values))/ sd(env_values))) |>
#     ungroup()
#     
    bbmo_env_zscore_w <-   bbmo_env_zscore_02 |>
      ungroup() |>
      bind_rows(bbmo_env_zscore_3) |>
      dplyr::select(sample_id, environmental_variable, z_score) |>
      pivot_wider(id_cols = sample_id, names_from = environmental_variable, values_from = z_score)
    
# bbmo_env_l |>
#  # group_by(enviornmental_variable) |>
#   calculate_z_score(col = 'env_values', group = 'enviornmental_variable')
#   
# bbmo_env |>
#   class()
# 
# bbmo_env |>
#   colnames()
# bbmo_env |>
#   glimpse()
# 
#   calculate_z_score(bbmo_env, col = day_length)
  
## plot environmental variables normalized by z-scores-----
  bbmo_env_z |>
    left_join(m_bbmo_10y) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    ggplot(aes(date, z_score_environmental_variable))+
    geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                      ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
    geom_point()+
    geom_line()+
    labs(x = 'Time', y = 'z-scores')+
    facet_wrap(vars(environmental_variable), scales = 'free_y')+
    scale_y_continuous()+
    theme_bw()+
    theme(panel.grid.minor.y = element_blank())

##calculate distance matrices
  abund_dist_02 <- asv_tab_rar |>
    dplyr::filter(str_detect(X, '0.2')) |>
    dplyr::select(-X) |>
    vegdist(method = 'bray')
  
  abund_dist_3 <- asv_tab_rar |>
    dplyr::filter(str_detect(X, '3')) |>
    dplyr::select(-X) |>
    vegdist(method = 'bray')

bbmo_env_02_dist <- bbmo_env |>
  dplyr::filter(str_detect(sample_id, '0.2')) |>
  dist(method = 'euclidean')

bbmo_env_3_dist <- bbmo_env |>
  dplyr::filter(str_detect(...28, '3')) |>
  dist(method = 'euclidean')

## here we use environmental data normalized using z_scores
bbmo_envz_02_dist<- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '0.2')) |>
  dist(method = 'euclidean')

bbmo_envz_3_dist <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '3')) |>
  dist(method = 'euclidean')

#Mantel test----
## The test statistic is the correlation coefficient. r falls in the range
## of -1 to +1, where being close to -1 indicates strong negative correlation
## and +1 indicates strong positive correlation. An r value of 0 indicates no
## correlation.
## 02
abund_dist_02 |>
  class()

abund_env_02 <- mantel(abund_dist_02, bbmo_env_02_dist,
                       method = 'spearman', permutations = 9999,
                       na.rm = TRUE)

#with env data normalized with z_scores
abund_envz_02 <- mantel(abund_dist_02, bbmo_envz_02_dist,
                       method = 'spearman', permutations = 9999,
                       na.rm = TRUE)

##3
abund_env_3 <- mantel(abund_dist_3, bbmo_env_3_dist,
                       method = 'spearman', 
                      permutations = 9999,
                       na.rm = TRUE)


## plot the results of the Mantel test----
correlation_coefficient <- abund_env_02$statistic
p_value <- abund_env_02$signif

results_mantel_02 <- tibble(Mantel_correlation = correlation_coefficient,
                            p_value = p_value)

results_mantel_02 |>
  ggplot(aes(Mantel_correlation, p_value))+
  geom_point()


# Option, correlation matrix with bloomers and environmental data
## subset ASV tab by potential bloomers
asv_bloo <- asv_tab_all_bloo_z_tax |>
  dplyr::select(asv_num) |>
  as_vector() |>
  unique() 

asv_tab_all_bloo_z_tax |>
  colnames()

abund_dist_02_bloo <- asv_tab_rar |>
  dplyr::select(X, all_of(asv_bloo)) |>
  dplyr::filter(str_detect(X, '0.2')) |>
  dplyr::select(-X) |>
  vegdist(method = 'bray')

abund_env_02 <- mantel(abund_dist_02_bloo, bbmo_env_02_dist,
                       method = 'spearman', permutations = 9999,
                       na.rm = TRUE)
library(corrplot)
qcorrplot(correlate(bbmo_env))


##visualization of correlation matrix
corr_bloo_02 <- cor(abund_dist_02_bloo, bbmo_envz_02_dist)

ggcorrplot::ggcorrplot(as.matrix(corr_bloo_02))
ggcorrplot::ggcorrplot(as.matrix(abund_dist_02))
ggcorrplot::ggcorrplot(as.matrix(abund_dist_3))
ggcorrplot::ggcorrplot(as.matrix(bbmo_env_02_dist))
ggcorrplot::ggcorrplot(as.matrix(bbmo_env_3_dist))

ggcorrplot::ggcorrplot(as.matrix(bbmo_envz_02_dist))
ggcorrplot::ggcorrplot(as.matrix(bbmo_envz_3_dist))

##
rda1 <- rda(abund_dist_02, bbmo_env3)

## Mantel correlogram
## This function comptues a multivariate Mantel correlogram.
mantel.correlog(abund_dist_02, bbmo_env_02_dist,
                nperm = 9999)

mantel.correlog(abund_dist_02, bbmo_envz_02_dist,
                nperm = 9999)

|>
  print.mantel.correlog()

bbmo_env_02_dist |>
  head()

bbmo_env_02_dist |>
  dim()

bbmo_env |>
  dim()

## other ideas
### perform a mantel test for each environmental variable and the community 


### perform a mantel test for each environmental variable and each potential bloomer


## Visualizing the data-----
asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  distinct(asv_num)

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
                fraction == 3 &
                  abundance_type == 'relative_abundance') |>
  ggplot(aes(temperature, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(temperature, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(secchi, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(salinity, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(chla_total, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(chla_3um, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(PO4, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(NH4, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(NO2, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(NO3, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(dryptomonas, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

## try the same plots but with z_scores
asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(dryptomonas_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(temperature_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(day_length_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(secchi_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(salinity_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(chla_total_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(PO4_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

## Multipatt function:
### This function studies the association between species patterns and combinations of groups of sites.
library(indicspecies)

multipatt(asv_tab_rar)


## multicorrelation
library(psych)
pairs.panels(bbmo_env[2:13], method = 'pearson')

library(ggcorrplot)
cor_pmat(asv_tab_rar)

ggcorr(asv_tab_rar, method = c("everything", "pearson")) 
