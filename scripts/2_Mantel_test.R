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

asv_tab_all_bloo_z_tax <- read.csv2('asv_tab_all_bloo_z_tax_new.csv')
asv_tab_rar <- read.csv2('asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()

#I create two different datasets one for ASVs and the other for the community
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
  distinct(...28,
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
#' #|>
#'   pivot_longer(cols = c('day_length',
#'                         #'sampling_time',
#'                         'temperature',
#'                         'secchi',
#'                         'salinity',
#'                         'chla_total',
#'                         'chla_3um',
#'                         'PO4',
#'                        'NH4', 'NO2', 'NO3',
#'                         'Si', 'BP_FC1.55',
#'                         'PNF_Micro', 'PNF2_5um_Micro',
#'                         'PNF_5um_Micro', 
#'                         'dryptomonas', 'micromonas',
#'                         'HNF_Micro', 'HNF2_5um_Micro',   
#'                         'HNF_5um_Micro' ,  'LNA',              
#'                         'HNA',   'prochlorococcus_FC',
#'                         'Peuk1',  'Peuk2',                 
#'                         'bacteria_joint', 'synechococcus'), 
#'                values_to = 'values', names_to = 'env_variable')

##normalization of environmental data using z-scores
## I created a function for this purpose 
calculate_z_score <- function(data, col, col_name = NULL,  group = NULL) { #,
 browser()  # Insert this line to pause execution and enter the interactive debugging mode
  stopifnot(is.numeric(data[[col]])) 

  col_name <- ifelse(!is.null(col_name), paste0("z_score_", col_name), 'z_score')
  
  if (!is.null(group)) {
    data <- data %>%
      group_by({{group}}) %>%
      mutate({{col_name}} := ({{col}} - mean({{col}}, na.rm = TRUE)) / sd({{col}}, na.rm = TRUE))
  } else {
    data <- data %>%
      mutate({{col_name}} := ({{col}} - mean({{col}}, na.rm = TRUE)) / sd({{col}}, na.rm = TRUE))
  }
  
  return(data)
}

calculate_z_score(bbmo_env, col = "day_length", col_name = 'day_length', group = NULL) ## this function works

#bbmo_env_z <- 
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
    dplyr::mutate(environmental_variable = as.numeric(environmental_variable)) |>
  calculate_z_score(col = 'env_values', group = environmental_variable, col_name = 'env_variable') 
  #pivot_wider(id_cols = )

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
  dplyr::filter(str_detect(...28, '0.2')) |>
  dist(method = 'euclidean')

bbmo_env_3_dist <- bbmo_env |>
  dplyr::filter(str_detect(...28, '3')) |>
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
corr_bloo_02 <- cor(abund_dist_02_bloo, bbmo_env_02_dist)

ggcorrplot::ggcorrplot(as.matrix(corr_bloo_02))
ggcorrplot::ggcorrplot(as.matrix(abund_dist_02))
ggcorrplot::ggcorrplot(as.matrix(abund_dist_3))
ggcorrplot::ggcorrplot(as.matrix(bbmo_env_02_dist))
ggcorrplot::ggcorrplot(as.matrix(bbmo_env_3_dist))

##
rda1 <- rda(abund_dist_02, bbmo_env3)

## Mantel correlogram
## This function comptues a multivariate Mantel correlogram.
mantel.correlog(abund_dist_02, bbmo_env_02_dist,
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

## Multipatt function:
### This function studies the association between species patterns and combinations of groups of sites.
library(indicspecies)

multipatt(asv_tab_rar)


## multicorrelation
library(psych)
pairs.panels(bbmo_env[2:13], method = 'pearson')
