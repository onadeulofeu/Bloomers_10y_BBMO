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

##packages----
library(vegan)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(multipanelfigure) #merge plots with different sizes
# 
# m_bbmo_10y |>
#   colnames()

## functions----
source('src/calculate_z_scores.R')

## environmental variables labs----
labs_env <- as_labeller(c("day_length" = 'Day length' ,
                                "temperature" = 'Temperature',
                                "secchi"  = 'Turbididty (Secchi disck)',    
                                "salinity" = 'Salinity',
                                "chla_total" = 'Total chl-a',
                                "chla_3um" = 'Chl-a 3 um',
                                "PO4"  =   'Phosphate',       
                                "NH4"     = 'Ammonia',
                                "NO2"   =  'Nitrite',          
                                "NO3" = 'Nitrate',
                                "Si"  = 'Silicate',
                                "BP_FC1.55" = 'Bacterial production (µgC l-1 d-1)',     
                                "PNF_Micro" = 'Phototrophic nanoflagellates',
                                "PNF2_5um_Micro"  = 'Phototrophic nanoflagellates (2-5um)',
                                "PNF_5um_Micro"   = 'Phototrophic nanoflagellates (5um)',
                                "dryptomonas"   = 'Dryptomonas',
                                "micromonas"  = 'Micromonas',
                                "HNF_Micro"  =  'Heterotrophic nanoflagellates',   
                                "HNF2_5um_Micro" ='Heterotrophic nanoflagellates (2-5um)',
                                "HNF_5um_Micro"  ='Heterotrophic nanoflagellates (5um)',
                                "LNA" = 'LNA',
                                "HNA"  = 'HNA',
                                "prochlorococcus_FC" = 'Prochlorococcus',
                                "Peuk1"  = 'Picoeukaryotes population 1',       
                                "Peuk2"   = 'Picoeukaryotes population 2',
                                "bacteria_joint" = 'Bacterial abundance',
                                "synechococcus"= 'Synechococcus',
                          'low_vlp' = 'Low viruses',
                          'med_vlp' = 'Mid viruses',
                          'high_vlp' = 'High viruses',
                          'total_vlp' = 'Total viruses'))

##upload data----
asv_tab_all_bloo_z_tax <- read.csv2('data/asv_tab_all_bloo_z_tax_new_assign.csv')
asv_tab_rar <- read.csv2('data/asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()

library(readxl) ##I upload the general metadata from the whole BBMO 20Y
bbmo_20y <- read_xlsx('data/main_databaseMOSTREIGBLANES_March23_od.xlsx', skip = 0 ) |>
  as_tibble()

bbmo_20y |>
  colnames()

bbmo_20y |>
  head()

##I would like to add viruses to the analysis----
bbmo_20y_v <- bbmo_20y |>
  dplyr::select(`NOM MOSTRA`, "Low VLP", "Med VLP" ,                           
                "High VLP" ,"total VLP")

bbmo_20y_v <- bbmo_20y_v |>
  rename(sample_id = `NOM MOSTRA`,
         low_vlp = "Low VLP",
         med_vlp = "Med VLP" ,  
         high_vlp = "High VLP" ,
         total_vlp = "total VLP") 
  
# bbmo_20y_v$sample_id
bbmo_20y_v_red <-  bbmo_20y_v |>
  dplyr::mutate(sample_id_ed = str_remove(sample_id, '_'),
                sample_id_ed2 = substr(sample_id_ed, 1, 8)) |>
  dplyr::filter(!is.na(sample_id_ed2)) |>
  dplyr::mutate(across(contains('_vlp'), as.numeric)) |>
  dplyr::filter(str_detect(sample_id_ed2, 'BL04|BL05|BL06|BL07|BL08|BL09|BL10|BL11|BL12|BL13|BL14')) |>
  dplyr::filter(!is.na(total_vlp)) |>
  dplyr::filter(!sample_id_ed2 == 'BL130709') |> #transecte DEVOTES
  dplyr::distinct(sample_id_ed2, total_vlp, low_vlp, med_vlp, high_vlp)

# Mantel test whole community structure vs environmental data------
### in this case we use the rarefied community to overcome the compositional problem previous to calculate
### the distances
## I create three different datasets one for ASVs, one for environmental data and the other for the community------
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
           bacteria_joint, synechococcus ) |>
  tidyr::separate(sample_id, into = c('sample_id_ed', 'filter', 'sequencing_num'), sep = '_', remove = FALSE) |>
  left_join(bbmo_20y_v_red, by = c('sample_id_ed' = 'sample_id_ed2'), relationship = "many-to-many") |>  ## add virus data to the environmental data 
  #rename('sample_id' = sample_id.x ) |>
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
           bacteria_joint, synechococcus,  low_vlp ,
           med_vlp ,  
           high_vlp ,
           total_vlp)

# bbmo_env$sample_id
# bbmo_20y_v |>
#   dim()
# 
# bbmo_env |>
#   dim()
## in case I need the env data separated by fractions but it's the same data-----
bbmo_env_02 <- bbmo_env |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) 
  # tidyr::separate(sample_id, into = c('sample_id_ed', 'filter', 'sequencing_num'), sep = '_', remove = FALSE) |>
  # left_join(bbmo_20y_v, by = c('sample_id_ed' = 'sample_id_ed2'), relationship = "many-to-many")

bbmo_env_3 <- bbmo_env |>
  dplyr::filter(str_detect(sample_id, '_3_')) 
  # tidyr::separate(sample_id, into = c('sample_id_ed', 'filter', 'sequencing_num'), sep = '_', remove = FALSE) |>
  # left_join(bbmo_20y_v, by = c('sample_id_ed' = 'sample_id_ed2'), relationship = "many-to-many")

##normalization of environmental data using z-scores----

# class(bbmo_env_sim$day_length)
# sum(is.na(bbmo_env_sim$day_length))

# bbmo_env_sim <- bbmo_env |>
#   dplyr::select(sample_id, day_length) |>
#   dplyr::filter(!is.na(day_length)) |>
#   dplyr::mutate(day_length = as.numeric(day_length)) |>
#   ungroup()

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

#calculate_z_score <- function(data, col, name = NULL, group = NULL) {
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
#       dplyr::group_by(!!sym(group)) |>
#       dplyr::mutate(!!col_name := (!!sym(col) - base::mean(!!sym(col), na.rm = TRUE)) / stats::sd(!!sym(col), na.rm = TRUE))
#   } else {
#     data <- data |>
#       dplyr::mutate(!!col_name := (!!sym(col) - base::mean(!!sym(col), na.rm = TRUE)) / stats::sd(!!sym(col), na.rm = TRUE))
#   }
#   
#   return(data)
# }

# calculate_z_score(bbmo_env_sim, col = 'day_length', name = 'day_length', group = NULL)
# 
# bbmo_env_sim |>
#   glimpse()
## Using the function created to normalize environmental data to z-scores----
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
                         bacteria_joint, synechococcus, 
                         low_vlp ,
                         med_vlp ,  
                         high_vlp ,
                         total_vlp), values_to = 'env_values', names_to = 'environmental_variable') |>
  dplyr::filter(!is.na(env_values)) |>
  dplyr::mutate(env_values = as.numeric(env_values)) |>
  calculate_z_score(col = 'env_values', name = 'environmental_variable', group = 'environmental_variable') 
  #pivot_wider(id_cols = )
#   
# bbmo_env_z$sample_id |>
#   unique()
  
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
## for the Mantel test analysis we need the table in a wider format -----
bbmo_env_zscore_w <-   bbmo_env_z |>
  dplyr::select(sample_id, environmental_variable, z_score_environmental_variable) |>
  pivot_wider(id_cols = sample_id, names_from = environmental_variable, values_from = z_score_environmental_variable)

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

### The restoration of the Blanes harbor strated on 24th March 2010 and finished on the 9th of june 2012----
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))
  

## reorder environmental factors ----
bbmo_env_z$environmental_variable <- factor(bbmo_env_z$environmental_variable, levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                                                                                          "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                                                                                          "bacteria_joint", "LNA", "HNA", "Peuk1",             
                                                                                          "Peuk2",
                                                                                          "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                                                                                          "HNF2_5um_Micro", "HNF_5um_Micro",
                                                                                          "low_vlp" ,
                                                                                          "med_vlp" ,  
                                                                                          "high_vlp" ,
                                                                                          "total_vlp"))
## plot environmental variables normalized by z-scores-----
bbmo_env_z |>
    left_join(m_bbmo_10y) |>
    # mutate(type_of_env = case_when(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,      
    #                                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", 
    #                                                              "bacteria_joint", "LNA", "HNA") ~ 'physico_chem',
    #                                environmental_variable %in% c(  "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
    #                                       "HNF2_5um_Micro", "HNF_5um_Micro", "prochlorococcus_FC", "Peuk1",             
    #                                       "Peuk2", "synechococcus") ~ 'biological')) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    ggplot(aes(date, z_score_environmental_variable))+
    geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                      ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
    geom_point(alpha = 0.6, size = 1/2)+
    geom_line()+
    labs(x = 'Time', y = 'z-scores')+
    facet_wrap(vars(environmental_variable), ncol = 2, scales = 'free_y', labeller = labs_env)+
    #facet_grid(vars(environmental_variable), scales = 'free_y', cols = vars(type_of_env), drop = T)+
    scale_y_continuous(labels = function(x) sprintf("%.1f", x),
                       expand = c(0,0))+
    scale_x_datetime(expand = c(0,0))+
    theme_bw()+
    theme(panel.grid.minor.y = element_blank(), 
          strip.background = element_blank(), 
          axis.text.y = element_text(size = 5),
          panel.grid.major.y = element_blank())

##calculate distance matrices----
asv_tab_bbmo_10y_l |> # upload the original table
  head()

####community data distance-----
## If we transform data to z-scores or rCLR we need to use Euclidean distance because we will have negative values, if
## we don't transform our data, then we can use Bray curtis distance

#####Rarefied dataset---
  abund_rar_02 <- asv_tab_rar |>
    dplyr::filter(str_detect(X, '_0.2_'))|>
    as.data.frame() #|>
    # dplyr::select(-X) |>
    # vegdist(method = 'bray')
  
  abund_rar_3 <- asv_tab_rar |>
    dplyr::filter(str_detect(X, '_3_'))|>
    as.data.frame() #|>
    # dplyr::select(-X) |>
    # vegdist(method = 'bray')
  
  # asv_tab_bbmo_10y_l |> # upload the original table
  #   head()
  
  asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_l |> #transform to wider format
    pivot_wider(id_cols = sample_id, names_from = asv_num, values_from = reads) |>
    as.data.frame()
  
  #####Transformed to rCLR----
  asv_tab_bbmo_10y_w_02 <- asv_tab_bbmo_10y_w |>
    dplyr::filter(str_detect(as.character(sample_id), '_0.2'))
  
  asv_tab_bbmo_10y_w_3 <- asv_tab_bbmo_10y_w |>
    dplyr::filter(str_detect(as.character(sample_id), '_3_'))
  
  row.names(abund_rar_02) <- abund_rar_02[,1] 
  row.names(abund_rar_3) <- abund_rar_3[,1] 

  row.names(asv_tab_bbmo_10y_w_02) <- asv_tab_bbmo_10y_w_02[,1]  
  row.names(asv_tab_bbmo_10y_w_3) <- asv_tab_bbmo_10y_w_3[,1] 
  
  abund_rar_02_ed <- abund_rar_02[,-1]
  abund_rar_3_ed <- abund_rar_3[,-1]
  
  asv_tab_bbmo_10y_w_02_ed <- asv_tab_bbmo_10y_w_02[,-1]
  asv_tab_bbmo_10y_w_3_ed <- asv_tab_bbmo_10y_w_3[,-1]
  
  abund_rar_02_rclr <-  abund_rar_02_ed |>
    decostand(method="rclr")
  abund_rar_3_rclr <-  abund_rar_3_ed |>
    decostand(method="rclr")
  
  asv_tab_bbmo_10y_w_02_rclr <-  asv_tab_bbmo_10y_w_02_ed |>
    decostand(method="rclr")
  asv_tab_bbmo_10y_w_3_rclr <-  asv_tab_bbmo_10y_w_3_ed |>
    decostand(method="rclr")
  
  # data.hel <- asv_tab_bbmo_10y_w_02_ed |>
  #   decostand(method="rclr"); str(data.hel)
  # 
  # data.hel |>
  #   dim()
  # 
  # data.dist_02 <- vegdist(data.hel, method="euclidean", na.rm = TRUE)
  # head(data.dist)
  # 
  # data.dist_02 |>
  #   class()
  
  ####env variables-----
# bbmo_env_02_dist <- bbmo_env |>
#   dplyr::filter(str_detect(sample_id, '0.2')) |>
#   dist(method = 'euclidean')
# 
# bbmo_env_3_dist <- bbmo_env |>
#   dplyr::filter(str_detect(...28, '3')) |>
#   dist(method = 'euclidean')

## here we use environmental data normalized using z_scores
# bbmo_envz_02_dist <- bbmo_env_zscore_w |>
#   dplyr::filter(str_detect(sample_id, '_0.2')) |>
#   dist(method = 'euclidean')
# 
# bbmo_envz_3_dist <- bbmo_env_zscore_w |>
#   dplyr::filter(str_detect(sample_id, '_3')) |>
#   dist(method = 'euclidean')

#Mantel test----
## The test statistic is the correlation coefficient. r falls in the range
## of -1 to +1, where being close to -1 indicates strong negative correlation
## and +1 indicates strong positive correlation. An r value of 0 indicates no
## correlation.
## 02
# abund_dist_02 |>
#   class()
# 
# abund_env_02 <- mantel(abund_dist_02, bbmo_env_02_dist,
#                        method = 'spearman', permutations = 9999,
#                        na.rm = TRUE)
# 
# #with env data normalized with z_scores
# abund_envz_02 <- mantel(abund_dist_02, bbmo_envz_02_dist,
#                        method = 'spearman', permutations = 9999,
#                        na.rm = TRUE)
# 
# ##3
# abund_env_3 <- mantel(abund_dist_3, bbmo_env_3_dist,
#                        method = 'spearman', 
#                       permutations = 9999,
#                        na.rm = TRUE)



## MANTEL TEST environmental variable (z-scores) vs community structure (transformed to zCLR)----

# mantel(data.dist_02, bbmo_envz_02_dist,
#        method = 'spearman', 
#        permutations = 9999,
#        na.rm = TRUE)

##one environmental variable at a time
  
##inputs FL (02) ----
  abund_rar_02_rclr
  bbmo_envz_02 <- bbmo_env_zscore_w |>
    dplyr::filter(str_detect(sample_id, '_0.2_')) |>
    as.data.frame()

  row.names(bbmo_envz_02) <-  bbmo_envz_02[,1] 
  bbmo_envz_02_ed <-  bbmo_envz_02[,-1]
# Extract common sample names
common_sample_names <- intersect(rownames( abund_rar_02_rclr), rownames(bbmo_envz_02_ed))

# Subset both matrices based on common sample names
abundance_data <- abund_rar_02_rclr[common_sample_names, ]
environmental_data <- bbmo_envz_02_ed[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_02_ed[, colnames(bbmo_envz_02_ed) != "other_variables"]  # Adjust as needed

# check that the names match
row.names(abund_rar_02_rclr) == row.names(bbmo_envz_02_ed)

# Initialize a list to store Mantel test results
mantel_results_02 <- list()

# Iterate over each environmental variable
for (variable in colnames(environmental_variables)) {
  
  # Extract the current environmental variable
  current_variable <- environmental_variables[, variable]
  
  # Calculate Euclidean distance
  euclidean_distance <- vegdist(current_variable, method = "euclidean", 
                                na.rm = TRUE)
  
  # Assuming you have another matrix or data frame for the second set of variables
  data.dist_02 <- vegdist(abund_rar_02_rclr, method="euclidean", 
                          na.rm = TRUE)
  # Replace 'your_second_set_of_variables' with your actual data
  # second_set_of_variables <-  asv_tab_bbmo_10y_w_02_ed |>
  #   dplyr::select(starts_with('asv'))

  # Calculate Euclidean distance for the second set of variables
  # euclidean_distance_second_set <- dist(second_set_of_variables, method = "euclidean")
  
  # Perform Mantel test
  mantel_test_result_02 <- mantel(euclidean_distance, data.dist_02, method = "pearson", permutations = 999,
                               na.rm = TRUE)
  
  # Store the Mantel test result in the list
  mantel_results_02[[variable]] <- mantel_test_result_02
}

# 'mantel_results_02' now contains Mantel test results for each environmental variable
# You can access the results using mantel_results$variable_name
# For example, mantel_results$"your_variable_name"$

# mantel_results_02$day_length
# mantel_results_02$temperature
# mantel_results_02$secchi
# mantel_results_02$salinity

# Create an empty tibble to store the summary results
results_mantel_02 <- tibble(
  variable_name = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each variable in 'mantel_results' and add the results to the tibble
for (variable_name in names(mantel_results_02)) {
  results <- mantel_results_02[[variable_name]]
  
    mantel_correlation <- results$statistic

  results_mantel_02 <- bind_rows(results_mantel_02, tibble(
    variable_name = variable_name,
    mantel_correlation = mantel_correlation,
    p_value = results$signif
  ))
}

# Print the summary tibble
print(results_mantel_02)

#write.csv2(results_mantel_02, 'results/mantel_community_02.csv')

#plot the results for the Mantel test between environmental variables and the community structure 
results_mantel_02$variable_name <- results_mantel_02$variable_name  |>
  factor(levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                    "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                    "bacteria_joint", "LNA", "HNA", "Peuk1",             
                    "Peuk2", 'total_vlp', 'low_vlp', 'med_vlp', 'high_vlp',
                    "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                    "HNF2_5um_Micro", "HNF_5um_Micro"))


plot_mantel_02_community <- results_mantel_02 |>
  mutate(community = 'Free living (0.2-3 um)') |>
  #left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  ggplot(aes(community, variable_name, fill = mantel_correlation))+
  geom_tile()+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  scale_y_discrete(labels = labs_env)+
  geom_text(aes(label = ifelse(p_value < 0.05, '*', '')))+
  labs(x = 'Community', y = 'Environmental variables', fill = 'Mantel correlation', title = 'Particle Attached')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), text = element_text(size = 5),
        axis.text.x = element_text(angle = 65, hjust = 1))


##inputs second PA (3) ----
abund_rar_3_rclr |>
  row.names()
bbmo_envz_3 <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  as.data.frame()

row.names(bbmo_envz_3) <-  bbmo_envz_3[,1] 
bbmo_envz_3_ed <-  bbmo_envz_3[,-1]
# Extract common sample names
common_sample_names <- intersect(rownames( abund_rar_3_rclr), rownames(bbmo_envz_3_ed))

# Subset both matrices based on common sample names
abundance_data <- abund_rar_3_rclr[common_sample_names, ]
environmental_data <- bbmo_envz_3_ed[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_3_ed[, colnames(bbmo_envz_3_ed) != "other_variables"]  # Adjust as needed

# check that the names match
row.names(abund_rar_3_rclr) == row.names(bbmo_envz_3_ed)

# Initialize a list to store Mantel test results
mantel_results_3 <- list()

# Iterate over each environmental variable
for (variable in colnames(environmental_variables)) {
  
  # Extract the current environmental variable
  current_variable <- environmental_variables[, variable]
  
  # Calculate Euclidean distance
  euclidean_distance <- vegdist(current_variable, method = "euclidean", 
                                na.rm = TRUE)
  
  # Assuming you have another matrix or data frame for the second set of variables
  data.dist_3 <- vegdist(abund_rar_3_rclr, method="euclidean", 
                          na.rm = TRUE)
  # Replace 'your_second_set_of_variables' with your actual data
  # second_set_of_variables <-  asv_tab_bbmo_10y_w_3_ed |>
  #   dplyr::select(starts_with('asv'))
  
  # Calculate Euclidean distance for the second set of variables
  # euclidean_distance_second_set <- dist(second_set_of_variables, method = "euclidean")
  
  # Perform Mantel test
  mantel_test_result_3 <- mantel(euclidean_distance, data.dist_3, method = "pearson", permutations = 999,
                                  na.rm = TRUE)
  
  # Store the Mantel test result in the list
  mantel_results_3[[variable]] <- mantel_test_result_3
}

# 'mantel_results_3' now contains Mantel test results for each environmental variable
# You can access the results using mantel_results$variable_name
# For example, mantel_results$"your_variable_name"$

# mantel_results_3$day_length
# mantel_results_3$temperature
# mantel_results_3$secchi
# mantel_results_3$salinity

# Create an empty tibble to store the summary results
results_mantel_3 <- tibble(
  variable_name = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each variable in 'mantel_results' and add the results to the tibble
for (variable_name in names(mantel_results_3)) {
  results <- mantel_results_3[[variable_name]]
  
  mantel_correlation <- results$statistic
  
  results_mantel_3 <- bind_rows(results_mantel_3, tibble(
    variable_name = variable_name,
    mantel_correlation = mantel_correlation,
    p_value = results$signif
  ))
}

# Print the summary tibble
print(results_mantel_3)

#plot the results for the Mantel test between environmental variables and the community structure 
results_mantel_3$variable_name <- results_mantel_3$variable_name  |>
  factor(levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                    "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                    "bacteria_joint", "LNA", "HNA", "Peuk1",             
                    "Peuk2", 'total_vlp', 'low_vlp', 'med_vlp', 'high_vlp',
                    "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                    "HNF2_5um_Micro", "HNF_5um_Micro"))

results_mantel_02 <-   results_mantel_02 |>
  dplyr::mutate(community = 'Free living (0.2-3 um)')

#write.csv2(results_mantel_3, 'results/mantel_community_3.csv')

##plot results from FL and PA together-----
plot_mantel_3_02_community <- results_mantel_3 |>
  dplyr::mutate(community = 'Particle attached (3-20 um)') |>
  bind_rows(results_mantel_02) |>
  dplyr::filter(p_value < 0.05 &
                  abs(mantel_correlation) > 0.1) |>
  #left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  ggplot(aes(community, variable_name, fill = mantel_correlation))+
  geom_tile()+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  scale_y_discrete(labels = labs_env)+
  #geom_text(aes(label = ifelse(p_value < 0.05, '*', ''), colour = '#757C76'))+
  labs(x = 'Community', y = 'Environmental variables', fill = 'Mantel correlation')+
  guides( colour = 'none')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), text = element_text(size = 5),
        axis.text.x = element_text(angle = 0), axis.ticks.x = element_blank(), panel.border = element_blank())

# ggsave('mantel_test_community_sig.pdf', plot_mantel_3_02_community, path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/results/figures/",
#        width = 88,
#        height = 88,
#        units = 'mm')

## I perform the Mantel test only with my potential bloomers vs each environmental variable----

### Prepare the inputs for the Mantel test
asv_anom_3_tb <- asv_anom_3 |>
  as_tibble()

asv_anom_02_tb <- asv_anom_02 |>
  as_tibble()

abundance_data_bloo02 <- asv_tab_bbmo_10y_w_02_ed |>
  dplyr::select(asv_anom_02_tb$value) # filter by potential bloomers in FL

abundance_data_bloo3 <- asv_tab_bbmo_10y_w_3_ed |>
  dplyr::select(asv_anom_3_tb$value) # filter by potential bloomers in PA

#abundance data transformed to rCLR: 
abundance_data_clr02 <- abundance_data_bloo02 |>
  decostand(method="clr", na.rm = TRUE, pseudocount = 1) 

abundance_data_clr02 |>
  dim()

abundance_data_clr3 <- abundance_data_bloo3 |>
  decostand(method="clr", na.rm = TRUE, pseudocount = 1) 

abundance_data_clr3 |>
  dim()

###FL -----
bbmo_envz_02 <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  as.data.frame()

row.names(bbmo_envz_02) <-  bbmo_envz_02[,1] 
bbmo_envz_02_ed <-  bbmo_envz_02[,-1]

# Extract common sample names
common_sample_names <- intersect(rownames(abundance_data_clr02), rownames(bbmo_envz_02_ed))

# Subset both matrices based on common sample names
abundance_data <- abundance_data_clr02[common_sample_names, ]
environmental_data <- abundance_data_clr02[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_02_ed[, colnames(bbmo_envz_02_ed) != "other_variables"]  # Adjust as needed

# check that the names match
row.names(abundance_data_clr02) == row.names(bbmo_envz_02_ed)

# Create an empty tibble to store the Mantel test results
mantel_results_02 <- tibble(
  taxon_name = character(),
  variable_name_environmental = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each taxon in the abundance matrix
for (taxon_name in colnames(abundance_data_clr02)) {
  
  # Extract the abundance data for the current taxon
  current_variable_abundance <- abundance_data_clr02[, taxon_name, drop = FALSE]
  
  # Calculate distance for the abundance variable
  distance_abundance <- vegdist(current_variable_abundance, method = "euclidean", na.rm = TRUE)
  
  # Iterate over each environmental variable
  for (variable_name_environmental in colnames(bbmo_envz_02_ed)) {
    
    # Extract the environmental data for the current variable
    current_variable_environmental <- bbmo_envz_02_ed[, variable_name_environmental, drop = FALSE]
    
    # Calculate distance for the environmental variable
    distance_environmental <- vegdist(current_variable_environmental, method = "euclidean", na.rm = TRUE)
    
    # Perform Mantel test
    mantel_test_result <- mantel(distance_abundance, distance_environmental, method = "pearson", permutations = 999, na.rm = TRUE)
    
    # Store the Mantel test result in the tibble
    mantel_results_02 <- bind_rows(mantel_results_02, tibble(
      taxon_name = taxon_name,
      variable_name_environmental = variable_name_environmental,
      mantel_correlation = mantel_test_result$statistic,
      p_value = mantel_test_result$signif
    ))
  }
}

# Print the summary tibble
print(mantel_results_02)

#write.csv2(mantel_results_02, 'results/mantel_results_02_bloo.csv')

## reorder environmental factors ----
mantel_results_02$variable_name_environmental <- mantel_results_02$variable_name_environmental  |>
  factor(levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                    "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                    "bacteria_joint", "LNA", "HNA", "Peuk1",             
                    "Peuk2", 'total_vlp', 'low_vlp', 'med_vlp', 'high_vlp',
                    "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                    "HNF2_5um_Micro", "HNF_5um_Micro"))

# plot the results
plot_mantel_02 <- mantel_results_02 |>
  left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  dplyr::mutate(community = 'Free living') |>
  ggplot(aes(interaction(taxon_name, family_f), variable_name_environmental, fill = mantel_correlation))+
  geom_tile()+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  scale_y_discrete(labels = labs_env)+
  facet_wrap(vars(community))+
  labs(x = 'Taxonomy', y = 'Environmental variables', fill = 'Mantel correlation')+ #, title = 'Free living fraction'
  geom_text(aes(label = ifelse(p_value < 0.05, '*', ''), colour = '#757C76'))+
  guides( colour = 'none')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), text = element_text(size = 5),
        axis.text.x = element_text(angle = 55, hjust = 1), strip.background = element_blank(),
        panel.border = element_blank())

####for the PA fraction------
bbmo_envz_3 <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  as.data.frame()

row.names(bbmo_envz_3) <-  bbmo_envz_3[,1] 
bbmo_envz_3_ed <-  bbmo_envz_3[,-1]

# Extract common sample names
common_sample_names <- intersect(rownames(abundance_data_clr3), rownames(bbmo_envz_3_ed))

# Subset both matrices based on common sample names
abundance_data <- abundance_data_clr3[common_sample_names, ]
environmental_data <- abundance_data_clr3[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_3_ed[, colnames(bbmo_envz_3_ed) != "other_variables"]  # Adjust as needed

# check that the names match
row.names(abundance_data_clr3) == row.names(bbmo_envz_3_ed)

# Create an empty tibble to store the Mantel test results
mantel_results_3 <- tibble(
  taxon_name = character(),
  variable_name_environmental = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each taxon in the abundance matrix
for (taxon_name in colnames(abundance_data_clr3)) {
  
  # Extract the abundance data for the current taxon
  current_variable_abundance <- abundance_data_clr3[, taxon_name, drop = FALSE]
  
  # Calculate distance for the abundance variable
  distance_abundance <- vegdist(current_variable_abundance, method = "euclidean", na.rm = TRUE)
  
  # Iterate over each environmental variable
  for (variable_name_environmental in colnames(bbmo_envz_3_ed)) {
    
    # Extract the environmental data for the current variable
    current_variable_environmental <- bbmo_envz_3_ed[, variable_name_environmental, drop = FALSE]
    
    # Calculate distance for the environmental variable
    distance_environmental <- vegdist(current_variable_environmental, method = "euclidean", na.rm = TRUE)
    
    # Perform Mantel test
    mantel_test_result <- mantel(distance_abundance, distance_environmental, method = "pearson", permutations = 999, na.rm = TRUE)
    
    # Store the Mantel test result in the tibble
    mantel_results_3 <- bind_rows(mantel_results_3, tibble(
      taxon_name = taxon_name,
      variable_name_environmental = variable_name_environmental,
      mantel_correlation = mantel_test_result$statistic,
      p_value = mantel_test_result$signif
    ))
  }
}

# Print the summary tibble
print(mantel_results_3)

#write.csv2(mantel_results_3, 'results/mantel_results_3_bloo.csv')

##reorder taxonomy as factors 
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num))

asv_tab_all_bloo_z_tax$class_f <-  factor(asv_tab_all_bloo_z_tax$class_f, 
                                          levels=unique(asv_tab_all_bloo_z_tax$class_f[order(asv_tab_all_bloo_z_tax$phylum_f)]), 
                                          ordered=TRUE)

asv_tab_all_bloo_z_tax$order_f <-  factor(asv_tab_all_bloo_z_tax$order_f, 
                                          levels=unique(asv_tab_all_bloo_z_tax$order_f[order(asv_tab_all_bloo_z_tax$phylum_f,
                                                                                             asv_tab_all_bloo_z_tax$class_f)]), 
                                          ordered=TRUE)

asv_tab_all_bloo_z_tax$family_f <-  factor(asv_tab_all_bloo_z_tax$family_f, 
                                           levels=unique(asv_tab_all_bloo_z_tax$family_f[order(asv_tab_all_bloo_z_tax$phylum_f,
                                                                                               asv_tab_all_bloo_z_tax$class_f,
                                                                                               asv_tab_all_bloo_z_tax$order_f)]), 
                                           ordered=TRUE)


asv_tab_all_bloo_z_tax$asv_num_f <-  factor(asv_tab_all_bloo_z_tax$asv_num_f, 
                                            levels=unique(asv_tab_all_bloo_z_tax$asv_num_f[order(asv_tab_all_bloo_z_tax$phylum_f,
                                                                                                 asv_tab_all_bloo_z_tax$class_f,
                                                                                                 asv_tab_all_bloo_z_tax$order_f,
                                                                                                 asv_tab_all_bloo_z_tax$family_f)]), 
                                            ordered=TRUE)

# Plot the results
tax_factors <- asv_tab_all_bloo_z_tax |>
  dplyr::select(asv_num_f, family_f, order_f, class_f, phylum_f) |>
  distinct(asv_num_f, family_f, order_f, class_f, phylum_f)

asv_tab_all_bloo_z_tax |>
  colnames()

palete_gradient_cb <- c(#"#240023",
  
  "#4db2a2" = 0,
  "#005a47",
  na.value = '#000000') 

## reorder environmental factors ----
mantel_results_3$variable_name_environmental <- mantel_results_3$variable_name_environmental  |>
  factor(levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                    "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                    "bacteria_joint", "LNA", "HNA", "Peuk1",             
                    "Peuk2", 'total_vlp', 'low_vlp', 'med_vlp', 'high_vlp',
                    "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                    "HNF2_5um_Micro", "HNF_5um_Micro"))

## add FL results to plot them together
mantel_results_02_tax <- mantel_results_02 |>
  left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  dplyr::mutate(community = 'Free living (0.2-3 um)')


plot_mantel_all <- mantel_results_3 |>
  left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  dplyr::mutate(community = 'Particle Attached (3-20 um)') |>
  bind_rows(mantel_results_02_tax) |>
  dplyr::filter(p_value < 0.05 &
                  abs(mantel_correlation) > 0.2) |>
  ggplot(aes(interaction(taxon_name, family_f), variable_name_environmental, fill = mantel_correlation))+
  geom_tile()+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  scale_y_discrete(labels = labs_env)+
  #coord_flip()+
  facet_wrap(vars(community), nrow = 2)+
  labs(x = 'Taxonomy', y = 'Environmental variables', fill = 'Mantel\ncorrelation')+ #, title = 'Free living fraction'
  #geom_text(aes(label = ifelse(p_value < 0.05, '*', ''), colour = '#757C76'))+
  guides(colour = 'none')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), text = element_text(size = 5),
        axis.text.x = element_text(angle = 55, hjust = 1), strip.background = element_blank(),
        panel.border = element_blank(), legend.key.size = unit(4, 'mm'), plot.margin = margin(0,0,5,15),
        legend.margin = margin(0,15,0,5))

# mantel_test_all <-  multi_panel_figure(columns = 1, rows = 2, width = 180, height = 200, 
#                                                row_spacing = 0.2, unit = 'mm',
#                                                panel_label_type = 'upper-alpha')
# 
# mantel_test_all  %<>%
#   fill_panel(plot_mantel_02, column = 1, row = 1) %<>%
#   fill_panel(plot_mantel_3, column = 1, row = 2)
# 
# ggsave('mantel_test_all_bloo_sig.pdf', plot_mantel_all ,
#        path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/results/figures/",
#        width = 88,
#        height = 100,
#        units = 'mm')

## TEST: plot the results of the Mantel test-------
# correlation_coefficient <- abund_env_02$statistic
# p_value <- abund_env_02$signif
# 
# results_mantel_02 <- tibble(Mantel_correlation = correlation_coefficient,
#                             p_value = p_value)
# 
# results_mantel_3 |>
#   ggplot(aes(Mantel_correlation, p_value))+
#   geom_point()
# 
# 
# # Option, correlation matrix with bloomers and environmental data
# ## subset ASV tab by potential bloomers
# asv_bloo <- asv_tab_all_bloo_z_tax |>
#   dplyr::select(asv_num) |>
#   as_vector() |>
#   unique() 
# 
# asv_tab_all_bloo_z_tax |>
#   colnames()
# 
# abund_dist_02_bloo <- asv_tab_rar |>
#   dplyr::select(X, all_of(asv_bloo)) |>
#   dplyr::filter(str_detect(X, '0.2')) |>
#   dplyr::select(-X) |>
#   vegdist(method = 'bray')
# 
# abund_env_02 <- mantel(abund_dist_02_bloo, bbmo_env_02_dist,
#                        method = 'spearman', permutations = 9999,
#                        na.rm = TRUE)
# library(corrplot)
# qcorrplot(correlate(bbmo_env))
# 
# 
# ##visualization of correlation matrix
# corr_bloo_02 <- cor(abund_dist_02_bloo, bbmo_envz_02_dist)
# 
# ggcorrplot::ggcorrplot(as.matrix(corr_bloo_02))
# ggcorrplot::ggcorrplot(as.matrix(abund_dist_02))
# ggcorrplot::ggcorrplot(as.matrix(abund_dist_3))
# ggcorrplot::ggcorrplot(as.matrix(bbmo_env_02_dist))
# ggcorrplot::ggcorrplot(as.matrix(bbmo_env_3_dist))
# 
# ggcorrplot::ggcorrplot(as.matrix(bbmo_envz_02_dist))
# ggcorrplot::ggcorrplot(as.matrix(bbmo_envz_3_dist))
# 
# ##
# rda1 <- rda(abund_dist_02, bbmo_env3)
# 
# ## Mantel correlogram
# ## This function comptues a multivariate Mantel correlogram.
# mantel.correlog(abund_dist_02, bbmo_env_02_dist,
#                 nperm = 9999)
# 
# mantel.correlog(abund_dist_02, bbmo_envz_02_dist,
#                 nperm = 9999)
# 
# |>
#   print.mantel.correlog()
# 
# bbmo_env_02_dist |>
#   head()
# 
# bbmo_env_02_dist |>
#   dim()
# 
# bbmo_env |>
#   dim()

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
