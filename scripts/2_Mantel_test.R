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
library(magrittr)
# 
# m_bbmo_10y |>
#   colnames()

## functions
source('src/calculate_z_scores.R')

## environmental variables labs
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
                                "synechococcus"= 'Synechococcus'))

##upload data----
asv_tab_all_bloo_z_tax <- read.csv2('data/asv_tab_all_bloo_z_tax_new_assign.csv')
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

# calculate_z_score(bbmo_env_sim, col = 'day_length', name = 'day_length', group = NULL)
# 
# bbmo_env_sim |>
#   glimpse()

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
## for the Mantel test analysis we need the table in a wider format 

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

### The remodelation of the Blanes harbour strated on 24th March 2010 and finished on the 9th of june 2012
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))
  
## plot environmental variables normalized by z-scores-----
  bbmo_env_z |>
    left_join(m_bbmo_10y) |>
    mutate(type_of_env = case_when(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,      
                                                                 "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", 
                                                                 "bacteria_joint", "LNA", "HNA") ~ 'physico_chem',
                                   environmental_variable %in% c(  "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                                          "HNF2_5um_Micro", "HNF_5um_Micro", "prochlorococcus_FC", "Peuk1",             
                                          "Peuk2", "synechococcus") ~ 'biological')) |>
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

####community-----
#####Rarefied
  abund_dist_02 <- asv_tab_rar |>
    dplyr::filter(str_detect(X, '0.2')) |>
    dplyr::select(-X) |>
    vegdist(method = 'bray')
  
  abund_dist_3 <- asv_tab_rar |>
    dplyr::filter(str_detect(X, '3')) |>
    dplyr::select(-X) |>
    vegdist(method = 'bray')
  
  #####Transformed to rCLR----
  asv_tab_bbmo_10y_l |> # upload the original table
    head()
  
  asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_l |> #transform to wider format
    pivot_wider(id_cols = sample_id, names_from = asv_num, values_from = reads) |>
    as.data.frame()
  
  asv_tab_bbmo_10y_w_02 <- asv_tab_bbmo_10y_w |>
    dplyr::filter(str_detect(as.character(sample_id), '_0.2'))
  
  asv_tab_bbmo_10y_w_3 <- asv_tab_bbmo_10y_w |>
    dplyr::filter(str_detect(as.character(sample_id), '_3_'))
  
  row.names(asv_tab_bbmo_10y_w_02) <- asv_tab_bbmo_10y_w_02[,1]  
  row.names(asv_tab_bbmo_10y_w_3) <- asv_tab_bbmo_10y_w_3[,1] 
  
  asv_tab_bbmo_10y_w_02_ed <- asv_tab_bbmo_10y_w_02[,-1]
  asv_tab_bbmo_10y_w_3_ed <- asv_tab_bbmo_10y_w_3[,-1]
  
  data.hel <- asv_tab_bbmo_10y_w_02_ed |>
    decostand(method="rclr"); str(data.hel)
  
  data.hel |>
    dim()
  
  data.dist_02 <- vegdist(data.hel, method="euclidean", na.rm = TRUE)
  head(data.dist)
  
  data.dist_02 |>
    class()
  
  ####env variables-----
bbmo_env_02_dist <- bbmo_env |>
  dplyr::filter(str_detect(sample_id, '0.2')) |>
  dist(method = 'euclidean')

bbmo_env_3_dist <- bbmo_env |>
  dplyr::filter(str_detect(...28, '3')) |>
  dist(method = 'euclidean')

## here we use environmental data normalized using z_scores
bbmo_envz_02_dist <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_0.2')) |>
  dist(method = 'euclidean')

bbmo_envz_3_dist <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_3')) |>
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



## MANTEL TEST environmental variable (z-scores) vs community structure (transformed to zCLR)----

# mantel(data.dist_02, bbmo_envz_02_dist,
#        method = 'spearman', 
#        permutations = 9999,
#        na.rm = TRUE)

##one environmental variable at a time
bbmo_envz_02 <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_0.2')) |>
  as.data.frame()

bbmo_envz_02 |>
  colnames()

rownames(bbmo_envz_02) <- bbmo_envz_02$sample_id

bbmo_envz_02 <- bbmo_envz_02[,-1] 

asv_tab_bbmo_10y_w_02_ed |>
 row.names()

row.names(asv_tab_bbmo_10y_w_02_ed) <- asv_tab_bbmo_10y_w_02_ed$sample_id

# Extract common sample names
common_sample_names <- intersect(rownames(asv_tab_bbmo_10y_w_02_ed), rownames(bbmo_envz_02))

# Subset both matrices based on common sample names
abundance_data <- asv_tab_bbmo_10y_w_02_ed[common_sample_names, ]
environmental_data <- bbmo_envz_02[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_02[, colnames(bbmo_envz_02) != "other_variables"]  # Adjust as needed

# Initialize a list to store Mantel test results
mantel_results_02 <- list()

# Iterate over each environmental variable
for (variable in colnames(environmental_variables)) {
  
  # Extract the current environmental variable
  current_variable <- environmental_variables[, variable]
  
  # Calculate Euclidean distance
  euclidean_distance <- vegdist(current_variable, method = "euclidean", na.rm = TRUE)
  
  # Assuming you have another matrix or data frame for the second set of variables
  # Transform ASV table to rCLR
  data.hel <- asv_tab_bbmo_10y_w_02_ed |>
    decostand(method="rclr")
  
  data.dist_02 <- vegdist(data.hel, method="euclidean", na.rm = TRUE)
  # Replace 'your_second_set_of_variables' with your actual data
  # second_set_of_variables <-  asv_tab_bbmo_10y_w_02_ed |>
  #   dplyr::select(starts_with('asv'))

  # Calculate Euclidean distance for the second set of variables
  # euclidean_distance_second_set <- dist(second_set_of_variables, method = "euclidean")
  
  # Perform Mantel test
  mantel_test_result <- mantel(euclidean_distance, data.dist_02, method = "pearson", permutations = 999,
                               na.rm = TRUE)
  
  # Store the Mantel test result in the list
  mantel_results_02[[variable]] <- mantel_test_result
}


row.names(asv_tab_bbmo_10y_w_02_ed) == row.names(bbmo_envz_02)
# 'mantel_results' now contains Mantel test results for each environmental variable
# You can access the results using mantel_results$variable_name
# For example, mantel_results$"your_variable_name"$

# mantel_results_02$day_length
# mantel_results_02$temperature
# mantel_results_02$secchi
# mantel_results_02$salinity

# Create an empty tibble to store the summary results
results_mantel <- tibble(
  variable_name = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each variable in 'mantel_results' and add the results to the tibble
for (variable_name in names(mantel_results_02)) {
  results <- mantel_results_02[[variable_name]]
  
  # Check the structure of 'results' to adapt the code
  if (!is.atomic(results$statistic)) {
    mantel_correlation <- results$statistic$observed
  } else {
    mantel_correlation <- results$statistic
  }
  
  results_mantel <- bind_rows(results_mantel, tibble(
    variable_name = variable_name,
    mantel_correlation = mantel_correlation,
    p_value = results$p
  ))
}

# Print the summary tibble
print(results_mantel)

#plot the results


## Testing to perform multiple Mantel tests over different variables vs each taxa -----
# 'Abundance_data' and 'environmental_data' have the same row order---
# Extract common sample names
common_sample_names <- intersect(rownames(asv_tab_bbmo_10y_w_02_ed), rownames(bbmo_envz_02))

# Subset both matrices based on common sample names
abundance_data <- asv_tab_bbmo_10y_w_02_ed[common_sample_names, ]
environmental_data <- bbmo_envz_02[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_02[, colnames(bbmo_envz_02) != "other_variables"]  # Adjust as needed

row.names(asv_tab_bbmo_10y_w_02_ed) == row.names(bbmo_envz_02)

# Create an empty tibble to store the Mantel test results
mantel_results <- tibble(
  variable_name_abundance = character(),
  variable_name_environmental = character(),
  mantel_correlation = double(),
  p_value = double()
)

##input a abundance data.frame and a z_scores environmental data.frame
abundance_data <- asv_tab_bbmo_10y_w_02_ed |>
  dplyr::select(asv1, asv11, asv17)
environmental_data <- bbmo_envz_02

# Iterate over each variable in the abundance matrix
for (variable_name_abundance in colnames(abundance_data)) {
  # Abundance data transform it to rCLR (overcome compositional limitations)
  abundance_data_clr <- abundance_data |>
    decostand(method="rclr")
  
  # Extract the abundance data for the current variable
  current_variable_abundance <- abundance_data_clr[, variable_name_abundance, drop = FALSE]
  
  # Calculate distance for the abundance variable
  distance_abundance <- vegdist(current_variable_abundance, method = "euclidean", na.rm = TRUE)
  
  # Iterate over each variable in the environmental matrix
  for (variable_name_environmental in colnames(environmental_data)) {
    
    # Extract the environmental data for the current variable
    current_variable_environmental <- environmental_data[, variable_name_environmental, drop = FALSE]
    
    # Calculate distance for the environmental variable
    distance_environmental <- vegdist(current_variable_environmental, method = "euclidean", na.rm = TRUE)
    
    # Perform Mantel test
    mantel_test_result <- mantel(distance_abundance, distance_environmental, method = "pearson", permutations = 999,
                                 na.rm = TRUE)
    
    # Store the Mantel test result in the tibble
    mantel_results <- bind_rows(mantel_results, tibble(
      variable_name_abundance = variable_name_abundance,
      variable_name_environmental = variable_name_environmental,
      mantel_correlation = mantel_test_result$statistic,
      p_value = mantel_test_result$p
    ))
  }
}

# Print the summary tibble
print(mantel_results)


###test outside the loop, once I get the analysis that I need then I will create the loop------

##one environmental data vs two ASVs
abundance_data <- asv_tab_bbmo_10y_w_02_ed |>
  dplyr::select(asv1, asv17, asv22)

#abundance data transformed to rCLR: 
abundance_data_clr <- abundance_data |>
  decostand(method="clr", na.rm = TRUE, pseudocount = 1) 

## select one taxa at a time prior to calculating the distances
abundance_data_clr_taxa1 <- abundance_data_clr[,taxa1]

##environmental data
environmental_data_var1 <- bbmo_envz_02[, temp]

## check that rownames are the same in both datasets
row.names(abundance_data) == row.names(abundance_data_clr)

##calculate distance matrix
distance_abund <- vegdist(abundance_data_clr_taxa1, 
                          method = "euclidean", 
                          na.rm = TRUE) #if we work with CLR we need Euclidean distance

distance_environmental <- vegdist(environmental_data_var1, 
                                  method = "euclidean", 
                                  na.rm = TRUE)

##Mantel test
mantel_test_result <- mantel(distance_abund, 
                             distance_environmental, 
                             method = "pearson", 
                             permutations = 999,
                             na.rm = TRUE)


#### another test (it seems to work)-----
abundance_data <- asv_tab_bbmo_10y_w_02_ed |>
  dplyr::select(asv1, asv17, asv22)

#abundance data transformed to rCLR: 
abundance_data_clr <- abundance_data |>
  decostand(method="clr", na.rm = TRUE, pseudocount = 1) 

##check for the same rownames in abundance data and environmental data
bbmo_envz_02 |>
  row.names() == abundance_data_clr |>
  row.names()

bbmo_envz_02 |>
  dim()

abundance_data |>
  dim()

# Create an empty tibble to store the Mantel test results
mantel_results <- tibble(
  taxon_name = character(),
  variable_name_environmental = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each taxon in the abundance matrix
for (taxon_name in colnames(abundance_data_clr)) {
  
  # Extract the abundance data for the current taxon
  current_variable_abundance <- abundance_data_clr[, taxon_name, drop = FALSE]
  
  # Calculate distance for the abundance variable
  distance_abundance <- vegdist(current_variable_abundance, method = "euclidean", na.rm = TRUE)
  
  # Iterate over each environmental variable
  for (variable_name_environmental in colnames(bbmo_envz_02)) {
    
    # Extract the environmental data for the current variable
    current_variable_environmental <- bbmo_envz_02[, variable_name_environmental, drop = FALSE]
    
    # Calculate distance for the environmental variable
    distance_environmental <- vegdist(current_variable_environmental, method = "euclidean", na.rm = TRUE)
    
    # Perform Mantel test
    mantel_test_result <- mantel(distance_abundance, distance_environmental, method = "pearson", permutations = 999, na.rm = TRUE)
    
    # Store the Mantel test result in the tibble
    mantel_results <- bind_rows(mantel_results, tibble(
      taxon_name = taxon_name,
      variable_name_environmental = variable_name_environmental,
      mantel_correlation = mantel_test_result$statistic,
      p_value = mantel_test_result$p
    ))
  }
}

# Print the summary tibble
print(mantel_results)



## I perform the Mantel test only with my potential bloomers vs each environmental variable----
## we need to transform it to 

asv_anom_3_tb <- asv_anom_3 |>
  as_tibble()

asv_anom_02_tb <- asv_anom_02 |>
  as_tibble()

abundance_data_bloo02 <- asv_tab_bbmo_10y_w_02_ed |>
  dplyr::select(asv_anom_02_tb$value)

abundance_data_bloo3 <- asv_tab_bbmo_10y_w_3_ed |>
  dplyr::select(asv_anom_3_tb$value)

#abundance data transformed to rCLR: 
abundance_data_clr02 <- abundance_data_bloo02 |>
  decostand(method="clr", na.rm = TRUE, pseudocount = 1) 

abundance_data_clr3 <- abundance_data_bloo3 |>
  decostand(method="clr", na.rm = TRUE, pseudocount = 1) 

###FL -----
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
  for (variable_name_environmental in colnames(bbmo_envz_02)) {
    
    # Extract the environmental data for the current variable
    current_variable_environmental <- bbmo_envz_02[, variable_name_environmental, drop = FALSE]
    
    # Calculate distance for the environmental variable
    distance_environmental <- vegdist(current_variable_environmental, method = "euclidean", na.rm = TRUE)
    
    # Perform Mantel test
    mantel_test_result <- mantel(distance_abundance, distance_environmental, method = "pearson", permutations = 999, na.rm = TRUE)
    
    # Store the Mantel test result in the tibble
    mantel_results_02 <- bind_rows(mantel_results_02, tibble(
      taxon_name = taxon_name,
      variable_name_environmental = variable_name_environmental,
      mantel_correlation = mantel_test_result$statistic,
      p_value = mantel_test_result$p
    ))
  }
}

# Print the summary tibble
print(mantel_results_02)


####for the PA fraction------
bbmo_envz_3 <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  as.data.frame()

rownames(bbmo_envz_3) <- bbmo_envz_3$sample_id

bbmo_envz_3 <- bbmo_envz_3[,-1] 

asv_tab_bbmo_10y_w_3_ed |>
  row.names()

row.names(asv_tab_bbmo_10y_w_3_ed) <- asv_tab_bbmo_10y_w_3_ed$sample_id


# Extract common sample names
common_sample_names <- intersect(rownames(asv_tab_bbmo_10y_w_3_ed), rownames(bbmo_envz_3))

rownames(asv_tab_bbmo_10y_w_3_ed) == rownames(bbmo_envz_3)


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
  for (variable_name_environmental in colnames(bbmo_envz_3)) {
    
    # Extract the environmental data for the current variable
    current_variable_environmental <- bbmo_envz_3[, variable_name_environmental, drop = FALSE]
    
    # Calculate distance for the environmental variable
    distance_environmental <- vegdist(current_variable_environmental, method = "euclidean", na.rm = TRUE)
    
    # Perform Mantel test
    mantel_test_result <- mantel(distance_abundance, distance_environmental, method = "pearson", permutations = 999, na.rm = TRUE)
    
    # Store the Mantel test result in the tibble
    mantel_results_3 <- bind_rows(mantel_results, tibble(
      taxon_name = taxon_name,
      variable_name_environmental = variable_name_environmental,
      mantel_correlation = mantel_test_result$statistic,
      p_value = mantel_test_result$p
    ))
  }
}

# Print the summary tibble
print(mantel_results)






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
