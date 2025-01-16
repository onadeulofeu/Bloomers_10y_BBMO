# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     EDMs CCM and MDR-S MAP                  ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Results created by Carmen-García Comas          ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# EDM empircal dynamic modeling
## These methods can be used to provide a mechanistic understanding of dynamical systems. 
## These non-linear statistical methods are rooted in state space reconstruction (SSR), 
## i.e. lagged coordinate embedding of time series data
## These methods do not assume any set of equations governing the system but recover the 
## dynamics from time series data, thus called empirical dynamic modeling. 
## Missing data impart an unavoidably negative influence on the performance of EDM.

## upload packages -----
library(scales)
library(rEDM)
library(pheatmap)
library(tidyverse)
library(Bloomers)
library(ggridges)
library(magrittr)
library(ggplot2)
library(ggalign)

# Palettes ----
palette_bloom_frequency <-   c('Non-Seasonal' = 'black', 
                                'No Bloomer' = "white", 
                                'Seasonal' = '#F7CDCD', 
                                'env' = '#4A785C')

palete_gradient_cb2 <- c('#A7FFE6',
                        '#FFBAA6',
                        #'#FFBAA6',
                        '#AA928B',
                        '#8BAA90')

palette_gradient <- c("white", '#DBDBDB', "#BDBEBE","#545454",
                      "#070607")

palette_gradient <- colorRampPalette(c("white",'#DBDBDB', '#BDBEBE', '#545454', '#070607'))

## Preprocessing 

### 1. Time series of variables should always be normalized to zero mean and unit variance 
### to ensure all variables have the same level of magnitude for comparison and to 
### avoid constructing a distorted state space.

### 2. Linear trends should be removed, either by simple regression or taking the first
### difference, to make the time series stationary.

### 3. Unless there is strong mechanistic reason, we recommend that the time series data 
### are not passed through a linear filter (.g., smoothing or moving average), because
### smoothers can destroy the dynamics and make the signal linear.

### 4. We caution that strong cyclic behavior or seasonality may mask the efficacy of EDM; 
### data standardization methods or surrogate data test to account for seasonality have been 
### developed to overcome these problems, although further methodological development 
### is still underway.

#  ------ ####### ------ PREPARE DATA ------ ####### ------
## We need occurrence and relative abundance - rank table-----
asv_tab_all_bloo_z_tax <- read_csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') |>
  as_tibble()

asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  group_by(asv_num) |>
  distinct(asv_num) #61 potential bloomers in my dataset

occurrence_bloo_bbmo <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(sample_id, date, total_reads, fraction, day, day_of_year, abundance_value, asv_num, seq) |>
  group_by(asv_num, fraction) |>
  dplyr::mutate(n_occurrence = case_when(abundance_value > 0 ~ 1,
                                         abundance_value == 0 ~ 0)) |>
  group_by(fraction) |>
  dplyr::mutate(n_samples = n_distinct(sample_id)) |>
  group_by(asv_num, fraction, seq) |>
  dplyr::mutate(occurrence_perc = sum(n_occurrence)/n_samples)

occurrence_bloo_bbmo_summary <-  occurrence_bloo_bbmo |>
  dplyr::select(-sample_id, -total_reads, -day, -day_of_year, -abundance_value, -n_occurrence) |>
  group_by(fraction, asv_num, seq, occurrence_perc) |>
  distinct() |>
  dplyr::filter(occurrence_perc >= 2/3) 

occurrence_bloo_bbmo_summary |>
  group_by(asv_num) |>
  distinct(asv_num) #12 ASVs with an occurrence >66% in the dataset.

#write.csv(occurrence_bloo_bbmo, 'data/occurrence_bloo_bbmo.csv')

## plots -----
occurrence_bloo_bbmo |>
  dplyr::filter(fraction == '3') |>
  dplyr::left_join(tax_bbmo_10y_new, by = c('seq')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(occurrence_perc = as.numeric(occurrence_perc),
                abundance_value = as.numeric(abundance_value)) |>
  ggplot(aes(date, occurrence_perc, color = family))+
  geom_point(aes(size = abundance_value, color = family, alpha = 0.8))+
  scale_y_continuous(labels = percent_format())+
  facet_wrap(vars(phylum))+
  scale_x_datetime()+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

occurrence_bloo_bbmo |>
  dplyr::filter(fraction == '3') |>
  dplyr::left_join(tax_bbmo_10y_new) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(occurrence_perc = as.numeric(occurrence_perc),
                abundance_value = as.numeric(abundance_value)) |>
  group_by(asv_num) |>
  dplyr::slice_max(order_by = abundance_value, n = 1)|>
  ggplot(aes(abundance_value, occurrence_perc, color = order))+
  geom_point(aes(size = abundance_value, color = order, alpha = 0.8))+
  scale_y_continuous(labels = percent_format())+
  scale_x_continuous(labels = percent_format())+
  scale_color_manual(values = palette_order_assigned_bloo)+
  geom_smooth(method = 'loess', aes(abundance_value, occurrence_perc, group = domain), se = F)+
  labs(size = 'Relative\nabundance (%)', color = 'Order', alpha = '', x = 'Relative abundance (%)', y = 'Occurrence')+
  theme_bw()+
  theme(text = element_text(size = 6),
        panel.grid = element_blank(),
        legend.position = 'bottom')

## upload occurrence data -----
occurrence_bloo_bbmo <- read.delim2('data/occurrence_bloo_bbmo.csv', sep = ',')
occurrence_bloo_bbmo |>
  head()

## Preparation of an ASV_tab with all ASVs of the temporal series that are present in 2/3 of the dataset and a tax table for them with bloomers identified ------
setwd("~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/")
bbmo_10y <-readRDS("blphy10years.rds") ##8,052 asv amb totes les mostres, reads / ASV and sample
#str(bbmo_10y)
#bbmo_10y <- prune_samples(sample_sums(bbmo_10y) > 10000, bbmo_10y) in case we want to filter samples by number of reads

bbmo_10y <-
  prune_taxa(taxa_sums(bbmo_10y@otu_table) >0, ## filter ASVs that are 0 in the whole dataset.
             bbmo_10y)

bbmo_10y |>
  nsamples() #237 samples 

bbmo_10y |>
  ntaxa() #8,052 but if we filter for those that are 0 during the whole dataset then we got 7,849 ASVs

## separate datasets by ASV_tab, taxonomy and metadata
asv_tab_bbmo_10y_l <- bbmo_10y@otu_table |>
  as_tibble()

m_bbmo_10y <- bbmo_10y@sam_data |>
  as_tibble()

# new taxonomy created with the database SILVA 138.1
new_tax <-  readRDS('data/03_tax_assignation/devotes_all_assign_tax_assignation_v2.rds') |>
  as_tibble(rownames = 'sequence')

tax_bbmo_10y_new <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax, by = c('seq' = 'sequence'))

## tidy colnames ----
asv_tab_bbmo_10y_l |>
  colnames()

tax_bbmo_10y_old |>
  colnames()

colnames(asv_tab_bbmo_10y_l) <- c('asv_num', "sample_id", 'reads')

colnames(tax_bbmo_10y_old) <- c("asv_num", "kingdom", "phylum", "class", "order", "family", "genus",
                                "species", "curated", "otu_corr","seq")

colnames(m_bbmo_10y) <- c( "project", "location", "code",             
                          "type", "samname", "fraction", "run",               
                          "date", "basics", "julian_day", "day_of_year",       
                          "decimal_date", "position", "sampling_time", "day_length",        
                          "temperature", "secchi", "salinity", "chla_total",     
                          "chla_3um", "PO4" ,"NH4", "NO2" ,              
                          "NO3",  "Si", "BP_FC1.55", "PNF_Micro",         
                          "PNF2_5um_Micro", "PNF_5um_Micro", "dryptomonas", "micromonas",        
                          "HNF_Micro", "HNF2_5um_Micro", "HNF_5um_Micro", "LNA",               
                          "HNA", "prochlorococcus_FC", "Peuk1",  "Peuk2",          
                          "Year", "Month", "Day", "season",            
                          "bacteria_joint", "synechococcus", "depth", "sample_id")

## Divide metadata into FL and PA----
m_02 <- m_bbmo_10y  |>
  dplyr::filter(fraction == 0.2) 

m_02 <- m_02 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02)))

m_3 <- m_bbmo_10y |>
  dplyr::filter(fraction == 3.0)

m_3 <- m_3 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3)))

## we lack of 3 samples in FL fraction which are:----
# m_3 |>
#   group_by(year) |>
#   summarize(n = n()) |>
#   dplyr::filter(n < 12)
# ## in 2004 and 2005, 
#  m_3 |>
#   dplyr::filter(year %in% c('2004', '2005')) |>
#    group_by(month) |>
#    mutate(n = n()) |>
#    dplyr::filter(n < 2) %$%
#    date

##[1] "2004-02-23" "2004-05-25" "2005-03-09" MISSING SAMPLES 

## check even sequencing----
### there is no even sequencing for this dataset, therefore when we calculate diversity we should use a rarefied dataset
# bbmo_m <- m_02 |>
#   bind_rows(m_3)
# 
# bbmo_m |>
#   colnames()
# 
# asv_tab_bbmo_10y_l |>
#   group_by(sample_id) |>
#   dplyr::summarize(n_reads = sum(reads)) |>
#   inner_join(bbmo_m) |>
#   ggplot(aes(date, n_reads))+
#   geom_col()+
#   facet_grid(vars(fraction))

## Calculate relative abundance----
asv_tab_10y_l_rel <- asv_tab_bbmo_10y_l |>
  calculate_rel_abund(group_cols = sample_id)

asv_tab_10y_3_rel <- asv_tab_10y_l_rel %>%
  dplyr::filter(sample_id %in% m_3$sample_id)

asv_tab_10y_02_rel <- asv_tab_10y_l_rel %>%
  dplyr::filter(sample_id %in% m_02$sample_id)

## we could also do it with rCLR add here in case we want to do it----
sample_id <- asv_tab_bbmo_10y_l |>
  pivot_wider(names_from = asv_num, values_from = reads, values_fill = 0) |>
  dplyr::select(sample_id)

asv_tab_bbmo_10y_rclr <- asv_tab_bbmo_10y_l |>
  pivot_wider(names_from = asv_num, values_from = reads, values_fill = 0) |>
  dplyr::select(-sample_id) |>
  decostand( method = 'rclr') |>
  as_tibble() |>
  bind_cols(sample_id)

# Calculate occurrence----
## Occurrence by fraction (0.2 - 3)
asv_tab_10y_l_rel |> 
  colnames()

m_bbmo_10y_sim <- m_bbmo_10y |>
  dplyr::select(fraction, sample_id, Year, Month, Day)

asv_tab_10y_l_rel_occ <- asv_tab_10y_l_rel |>
  dplyr::left_join(m_bbmo_10y_sim) |>
  group_by(fraction, asv_num) |>
  dplyr::mutate(n_occurrence = case_when(relative_abundance > 0 ~ 1,
                                         relative_abundance == 0 ~ 0)) |>
 group_by(fraction) |>
  dplyr::mutate(n_samples = n_distinct(sample_id)) |>
  group_by(asv_num, fraction) |>
  dplyr::mutate(occurrence_perc = sum(n_occurrence)/n_samples) |>
  dplyr::select(Month, Year, Day, asv_num, fraction, relative_abundance, occurrence_perc, n_samples)

asv_tab_10y_l_rel_occ_filt <- 
  asv_tab_10y_l_rel_occ |>
  dplyr::filter(occurrence_perc >= 2/3)

asv_tab_10y_l_rel_occ_filt |>
  dplyr::group_by(fraction, asv_num) |>
  summarize(n = n_distinct(fraction, asv_num)) |> 
  group_by(fraction) |>
  summarize(n = sum(n))## I keep 47 ASVs in FL and 19 in PA data

occurrence_bloo_bbmo_summary <-  asv_tab_10y_l_rel_occ_filt |>
  group_by(fraction, asv_num, occurrence_perc) |>
  distinct() |>
  dplyr::filter(occurrence_perc > 2/3) 

asv_tab_10y_rel_occ_filt_w <- asv_tab_10y_l_rel_occ_filt |>
  dplyr::mutate(fraction_ed = case_when(str_detect(fraction, '0.2') ~ 'bp',
                                     str_detect(fraction, '3') ~ 'bn'),
                f_asv_num = paste0(fraction_ed,'_',asv_num)) |>
  dplyr::select(Year, Month, Day, f_asv_num, relative_abundance) |>
  pivot_wider(id_cols = c(Year, Month, Day), names_from = f_asv_num, values_from = relative_abundance, values_fill = 0 )

asv_tab_10y_rel_occ_filt_w |>
  dim()

## check that I'm NOT filtering by ASV num
asv_tab_10y_l_rel_occ_filt |>
  dplyr::ungroup() |>
  dplyr::select(asv_num) |>
  distinct() |>
  dplyr::filter(asv_num %in% bloo_02$value)

asv_tab_10y_l_rel_occ_filt |>
  dplyr::ungroup() |>
  dplyr::select(asv_num) |>
  distinct() |>
  dplyr::filter(asv_num %in% bloo_3$value)

##add metadata at the last cols of the dataset
m_bbmo_10y |>
  colnames()

m_bbmo_10y_sim |>
  dim()

m_bbmo_10y_sim_ed <- m_bbmo_10y |>
  dplyr::select(Day, Month, Year, day_length, temperature, secchi, salinity, chla_total, chla_3um, PO4, NH4, NO2, NO3, Si, BP_FC1.55, PNF_Micro,         
                PNF2_5um_Micro, PNF_5um_Micro, dryptomonas, micromonas, HNF_Micro, HNF2_5um_Micro, HNF_5um_Micro,
                LNA, HNA, prochlorococcus_FC, Peuk1, Peuk2, bacteria_joint, synechococcus) |>
  distinct(Day, Month, Year , day_length, temperature, secchi, salinity, chla_total, chla_3um, PO4, NH4, NO2, NO3, Si, BP_FC1.55, PNF_Micro,         
           PNF2_5um_Micro, PNF_5um_Micro, dryptomonas, micromonas, HNF_Micro, HNF2_5um_Micro, HNF_5um_Micro,
           LNA, HNA, prochlorococcus_FC, Peuk1, Peuk2, bacteria_joint, synechococcus)

asv_tab_10y_rel_occ_filt_w_env <- asv_tab_10y_rel_occ_filt_w |>
  left_join(m_bbmo_10y_sim_ed, by = c('Year' = 'Year', 'Month' = 'Month', 'Day' = 'Day'))

#write.csv2(asv_tab_10y_rel_occ_filt_w_env, file = '../../EDM_carmen/asv_tab_10y_rel_occ_filt_w_env.csv')

##prepare taxonomy information
asv_tab_10y_l_rel_occ_filt |>
  ungroup() |>
  dplyr::distinct(asv_num) |>
  left_join(tax_bbmo_10y_old, by = 'asv_num') |>
  dplyr::select(asv_num, seq) |>
  left_join(tax_bbmo_10y_new, by = 'seq') |>
  dplyr::mutate(check_tax = (asv_num.x == asv_num.y)) |>
  summary(check_tax) ##being sure that the asv_num is conserved with both taxonomies.

bloo <- occurrence_bloo_bbmo |>
  ungroup() |>
  distinct(asv_num)

tax_occ_filt_bbmo <- asv_tab_10y_l_rel_occ_filt |>
  ungroup() |>
  dplyr::distinct(asv_num) |>
  left_join(tax_bbmo_10y_old, by = 'asv_num') |>
  dplyr::select(asv_num, seq) |>
  left_join(tax_bbmo_10y_new, by = c('seq', 'asv_num')) |>
  dplyr::mutate(bloom = case_when(asv_num %in% bloo$asv_num ~ '1',
                                  !asv_num %in% bloo$asv_num ~ '0'))

#write.csv(tax_occ_filt_bbmo, '../../EDM_carmen/tax_occ_filt_bbmo_ed.csv')
 
occ_asv <- asv_tab_10y_l_rel_occ_filt |>
  ungroup() |>
  distinct(asv_num) |>
  as_vector() ##add indication of bloomer

## Occurrence in general (237 samples) (the results are the same as previously) ---- 
asv_tab_10y_l_rel |> 
  colnames()
m_bbmo_10y_sim <- m_bbmo_10y |>
  dplyr::select(fraction, sample_id, Year, Month, Day)

asv_tab_10y_l_rel_occ <- asv_tab_10y_l_rel |>
  dplyr::left_join(m_bbmo_10y_sim) |>
  group_by(fraction, asv_num) |>
  dplyr::mutate(n_occurrence = case_when(relative_abundance > 0 ~ 1,
                                         relative_abundance == 0 ~ 0)) |>
  dplyr::mutate(n_samples = n_distinct(sample_id)) |>
  group_by(asv_num, fraction) |>
  dplyr::mutate(occurrence_perc = sum(n_occurrence)/n_samples) |>
  dplyr::select(Month, Year, Day, asv_num, fraction, relative_abundance, occurrence_perc, n_samples)

asv_tab_10y_l_rel_occ_filt <- 
  asv_tab_10y_l_rel_occ |>
  group_by(asv_num, occurrence_perc) |>
  dplyr::filter(occurrence_perc >= 2/3)

asv_tab_10y_l_rel_occ_filt |>
  dplyr::group_by(fraction, asv_num) |>
  summarize(n = n_distinct(fraction, asv_num)) |> 
  group_by(fraction) |>
  summarize(n = sum(n))## I keep 47 ASVs in FL and 19 in PA data

occurrence_bloo_bbmo_summary <-  asv_tab_10y_l_rel_occ_filt |>
  group_by(fraction, asv_num, occurrence_perc) |>
  distinct() |>
  dplyr::filter(occurrence_perc > 2/3) 

asv_tab_10y_l_rel_occ_filt |>
  dplyr::mutate(fraction_ed = case_when(str_detect(fraction, '0.2') ~ 'bp',
                                        str_detect(fraction, '3') ~ 'bn'),
                f_asv_num = paste0(fraction_ed,'_',asv_num)) |>
  dplyr::select(Year, Month, Day, f_asv_num, relative_abundance) |>
  pivot_wider(id_cols = c(Year, Month, Day), names_from = f_asv_num, values_from = relative_abundance, values_fill = 0 )

## filter rCLR by those ASVs that were pesent in 2/3 of the dataset -----
asv_occurrence_frac <- asv_tab_10y_l_rel_occ_filt |>
  dplyr::select(asv_num, fraction) |>
  dplyr::distinct() |>
  dplyr::mutate(fraction_asv = paste0(fraction,'_',asv_num))

asv_tab_bbmo_10y_rclr_occ_filt <- asv_tab_bbmo_10y_rclr |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'clr') |>
  dplyr::mutate(fraction = case_when(str_detect(sample_id, '_0.2') ~ '0.2',
                                        str_detect(sample_id, '_3') ~ '3')) |>
  dplyr::mutate(fraction_asv = paste0(fraction,'_',asv_num)) |>
  dplyr::filter(fraction_asv %in% asv_occurrence_frac$fraction_asv) |>
  dplyr::mutate(fraction_ed = case_when(str_detect(fraction, '0.2') ~ 'bp',
                                        str_detect(fraction, '3') ~ 'bn'),
                f_asv_num = paste0(fraction_ed,'_',asv_num)) |>
  dplyr::left_join(m_bbmo_10y_sim) |>
  dplyr::select(Year, Month, Day, f_asv_num, clr) |>
  pivot_wider(id_cols = c(Year, Month, Day), names_from = f_asv_num, values_from = clr, values_fill = 0 ) |>
  left_join(m_bbmo_10y_sim_ed, by = c('Year' = 'Year', 'Month' = 'Month', 'Day' = 'Day')) ## add metadata at the end of the tibble 
  
#write.csv(asv_tab_bbmo_10y_rclr_occ_filt, '../../EDM_carmen/asv_tab_bbmo_10y_rclr_occ_filt.csv')
  
  asv_tab_bbmo_10y_rclr_occ_filt |>
    dim()
  
  asv_tab_10y_rel_occ_filt_w_env |>
    dim()
  
# Creation of a new dataset that has less missing environmental variables and less of them, I send them to Carmen ----- 
## recover the columns of year, month, day to add them to the metadata
  y_m_d <- m_bbmo_10y |>
    dplyr::select(Year, Month, Day, decimal_date) |>
    distinct(Year, Month, Day, decimal_date)
  
env_data_interpolated_values_all_z_score_red <- read.csv2( '../data/env_data/env_data_interpolated_values_all_z_score.csv') |>
  as_tibble() |>
  dplyr::select(-X, -env_values) |>
  dplyr::mutate(environmental_variable = str_replace(environmental_variable, '_no_nas', '')) |>
  pivot_wider(id_cols = 'decimal_date', names_from = 'environmental_variable', values_from = 'z_score_environmental_variable') |>
  right_join(y_m_d ) |>
  dplyr::select(-decimal_date) 
  
asv_tab_10y_rel_occ_filt_w_env_inter <- asv_tab_10y_rel_occ_filt_w |>
  left_join(env_data_interpolated_values_all_z_score_red)
  
#write.csv(asv_tab_10y_rel_occ_filt_w_env_inter, file = '../../EDM_carmen/asv_tab_10y_rel_occ_filt_w_env_inter.csv')
  
asv_tab_bbmo_10y_rclr_occ_filt |>
  colnames()

asv_tab_bbmo_10y_rclr_occ_filt_inter <- asv_tab_bbmo_10y_rclr_occ_filt |>
  left_join(env_data_interpolated_values_all_z_score_red)

#write.csv(asv_tab_bbmo_10y_rclr_occ_filt_inter, '../../EDM_carmen/asv_tab_bbmo_10y_rclr_occ_filt_inter.csv')

### Remove seasonality from the data before the analysis so that we reduce the Moran effect -----

#write.csv(asv_tab_10y_rel_occ_filt_w_env_inter, file = '../../EDM_carmen/asv_tab_10y_rel_occ_filt_w_env_inter.csv')

#### rclr
asv_tab_bbmo_10y_rclr_occ_filt_inter <- read.csv( '../EDM_carmen/marbits/asv_tab_bbmo_10y_rclr_occ_filt_inter.csv') |>
  dplyr::select(-'X')

## remove seasonality 
mean_sd_season <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  pivot_longer(cols = !c('Year', 'Month', 'Day'), values_to = 'rclr', names_to = 'asv_num') |>
  group_by(Month, asv_num) |>
  dplyr::reframe(mean = mean(rclr, na.rm = T), sd = sd(rclr, na.rm = T))
  
asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  pivot_longer(cols = !c('Year', 'Month', 'Day'), values_to = 'rclr', names_to = 'asv_num') |>
  left_join(mean_sd_season ) |>
  ungroup() |>
  dplyr::mutate(rclr_no_seas = case_when(sd != 0.0000 ~ (rclr-mean)/sd,
                                         sd == 0.0000 ~ 0)) |> ## division by 0 is NaN
  ungroup() |>
  dplyr::select(-mean, -sd) |>
  pivot_wider(id_cols = c('Year', 'Month', 'Day'), names_from = 'asv_num', values_from = 'rclr_no_seas',
              values_fill = 0)

#write.csv(asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas, file = '../EDM_carmen/asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas.csv')
#write.csv(asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas, file = '../EDM_carmen/otu_table067_deseasmonthPep.csv', row.names = F) ## change the name of the table for Carmen analysis

#### rel abund ----
asv_tab_bbmo_10y_rel_occ_filt_inter <- read.csv( '../EDM_carmen/asv_tab_10y_rel_occ_filt_w_env_inter.csv') |>
  dplyr::select(-'X')

## remove seasonality 
mean_sd_season <- asv_tab_bbmo_10y_rel_occ_filt_inter |>
  pivot_longer(cols = !c('Year', 'Month', 'Day'), values_to = 'rclr', names_to = 'asv_num') |>
  group_by(Month, asv_num) |>
  dplyr::reframe(mean = mean(rclr), sd = sd(rclr))

asv_tab_bbmo_10y_rel_occ_filt_inter_noseas <- asv_tab_bbmo_10y_rel_occ_filt_inter |>
  pivot_longer(cols = !c('Year', 'Month', 'Day'), values_to = 'rclr', names_to = 'asv_num') |>
  left_join(mean_sd_season ) |>
  dplyr::mutate(rclr_no_seas = (rclr-mean)/sd) |>
  ungroup() |>
  dplyr::select(-mean, -sd) |>
  pivot_wider(id_cols = c('Year', 'Month', 'Day'), names_from = 'asv_num', values_from = 'rclr_no_seas' )

#write.csv(asv_tab_bbmo_10y_rel_occ_filt_inter_noseas, file = '../EDM_carmen/asv_tab_bbmo_10y_rel_occ_filt_inter_noseas.csv')

# Create a table with presence absence of each AS to not infiere relationships when an ASV is not there -----
## In this table 0 will indicate presence of an ASV and 1 absence!

asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  colnames() ==   
  asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas |>
  colnames()

asv_tab_bbmo_10y_rclr_occ_filt_inter$Day == asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas$Day
asv_tab_bbmo_10y_rclr_occ_filt_inter$Month == asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas$Month
asv_tab_bbmo_10y_rclr_occ_filt_inter$Year == asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas$Year

poszero_bloomBlanes_tb <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::mutate(across(-c(Day, Year, Month), ~ ifelse(. == 0.0000, 1, 0)))

#write.csv2(poszero_bloomBlanes_tb, file = '../EDM_carmen/poszero_bloomBlanes.csv', row.names = F)

#### The analysis is computed in MARBTIS


# ------ ####### ------ EXPLORE EDM RESULTS ------ ####### ------
## Boxplots causality -----
#### Si tienes tiempo, puedes hacer boxplots de causalidad intra e inter clase taxonómica por ejemplo y 
#### también ver un poco si las relaciones causales tops (primeros en cada ranking de las scatter plots que te mandé) tienen sentido
  
## CCM convergent cross mapping results computed in MARBITS (it's a bit computationally expensive) ----
### We have three different dataset 1. removed seasonality and rCLR data transformed, 2. removed seasonality and rel abund data, 
### 3. seasonality not removed and rCLR data transformed and 4. seasonality not removed and rel abund data transformed.

### 1. DATASET: We removed seasonality & data is in the format of rCLR -------
#### upload data
# asv_tab_bbmo_10y_rclr_occ_filt_inter <- read.csv( '../EDM_carmen/marbits/asv_tab_bbmo_10y_rclr_occ_filt_inter_noseas_ed.csv',
#                                                   sep = ';') 
# ccm_rho <- read.csv( '../EDM_carmen/marbits/ccm_rclr_no_seas/ccm_rho_Bl_rclr_noseas.csv', 
#                      sep = ',')  # jacobians matrix
# ccm_sig <- read.csv( '../EDM_carmen/marbits/ccm_rclr_no_seas/ccm_sig_Bl_rclr_noseas.csv', 
#                      sep = ',') 

## new results with the updated code 
asv_tab_bbmo_10y_rclr_occ_filt_inter <- read.csv( '../EDM_carmen/otu_table067_deseasmonthPep.csv',
                                                  sep = ',')
ccm_rho <- read.csv( '../EDM_carmen/results_from_carmen_V2/ccm_rho_Bl_nin120_demo.csv',
                     sep = ',')  # jacobians matrix
ccm_sig <- read.csv( '../EDM_carmen/results_from_carmen_V2/ccm_sig_Bl_nin120_demo.csv',
                     sep = ',')

# number of ASVs included in this analysis
num_asv_included <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  colnames() |>
  as_tibble_col() |>
  dplyr::filter(str_detect(value, 'asv')) |>
  dplyr::filter(str_detect(value, 'bn_')) |>
  dplyr::mutate(value = str_replace(value,'bn_', '')) |>
  distinct(value)

num_asv_included <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  colnames() |>
  as_tibble_col() |>
  dplyr::filter(str_detect(value, 'asv')) |>
  dplyr::filter(str_detect(value, 'bp_')) |>
  dplyr::mutate(value = str_replace(value,'bp_', '')) |>
  distinct(value)

ccm_rho |>
  colnames() <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_sig |>
  colnames() <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_rho |>
  dim()

##keep only those Rho values that are significant (in the matrix sig as 1)
ccm_rho |>
  row.names() <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_rho <- ccm_rho |>
  as.matrix()

ccm_rho_filt <- ifelse(ccm_sig == 1, ccm_rho, NA)

ccm_rho_filt |>
  row.names() <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

## i need to remove the cells where the have the same row name and col name because is the result of being predicted by itself and it won't make sense
# replace cells where row name equals column name with 0
ccm_rho_filt[cbind(match(row.names(ccm_rho_filt)[row.names(ccm_rho_filt) %in% colnames(ccm_rho_filt)], row.names(ccm_rho_filt)), 
        match(row.names(ccm_rho_filt)[row.names(ccm_rho_filt) %in% colnames(ccm_rho_filt)], colnames(ccm_rho_filt)))] <- 0

## input for evaluating bloomers causal variables.

#### explore the distribution of rho values ----
ccm_rho_filt |>
  as_tibble() |>
  pivot_longer(cols = everything(), values_to = 'rho', names_to = 'variables') |>
  dplyr::filter(!is.na(rho)) |>
  dplyr::mutate(type = 'ccm_results_rho') |>
  ggplot(aes(x = rho, y = type))+
  geom_density_ridges()+
  theme_ridges()

ccm_rho_filt |>
  as_tibble() |>
  pivot_longer(cols = everything(), values_to = 'rho', names_to = 'variables') |>
  dplyr::filter(!is.na(rho)) |>
  dplyr::summarise(mean = mean(rho),
                   sd = sd(rho))

ccm_rho_filt |>
  as.matrix() |>
  heatmap()

ccm_rho_filt |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  dplyr::filter(rho > 0.5) |> 
  ggplot(aes(variables, variables_2))+
  geom_tile(aes(fill = rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.border = element_blank())

ccm_rho_filt |>
  dim()

# Set a threshold for the maximum number of allowed NA values
na_threshold <- 100 ## to filter by it i need to have a values lower than 85 which is the dimension of my matrix. I don't want to filter in this case.
## 

# Identify columns where the number of NA values is less than or equal to the threshold
cols_below_threshold <- apply(ccm_rho_filt , 2, function(col) sum(is.na(col)) <= na_threshold)

# Remove columns where the number of NA values exceeds the threshold
matrix_clean <- ccm_rho_filt [, cols_below_threshold]
matrix_clean <- ccm_rho_filt [cols_below_threshold, ]

matrix_clean |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  dplyr::filter(rho > 0.4) |> 
  ggplot(aes(variables, variables_2))+
  geom_tile(aes(fill = rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.border = element_blank())

### 2. DATASET: We removed seasonality & data is in the format of relative abundance -----
asv_tab_bbmo_10y_rel_occ_filt_inter <- read.csv( '../EDM_carmen/marbits/asv_tab_bbmo_10y_rel_occ_filt_inter_noseas.csv', sep = ',') 
asv_tab_bbmo_10y_rel_occ_filt_inter  <- asv_tab_bbmo_10y_rel_occ_filt_inter  |>
  dplyr::select(-'X')

ccm_rho <- read.csv( '../EDM_carmen/marbits/ccm_rel_abund_no_seas/ccm_rho_Bl_rel_noseas.csv', sep = ',') 
ccm_sig <- read.csv( '../EDM_carmen/marbits/ccm_rel_abund_no_seas/ccm_sig_Bl_rel_noseas.csv', sep = ',') 

ccm_rho |>
  colnames() <- asv_tab_bbmo_10y_rel_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_sig |>
  colnames() <- asv_tab_bbmo_10y_rel_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

##keep only those Rho values that are significant (in the matrix sig as 1)
ccm_rho |>
  row.names() <- asv_tab_bbmo_10y_rel_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_rho <- ccm_rho |>
  as.matrix()

ccm_rho_filt <- ifelse(ccm_sig == 1, ccm_rho, NA)

ccm_rho_filt |>
  row.names() <- asv_tab_bbmo_10y_rel_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

# replace cells where row name equals column name with 0
ccm_rho_filt[cbind(match(row.names(ccm_rho_filt)[row.names(ccm_rho_filt) %in% colnames(ccm_rho_filt)], row.names(ccm_rho_filt)), 
                   match(row.names(ccm_rho_filt)[row.names(ccm_rho_filt) %in% colnames(ccm_rho_filt)], colnames(ccm_rho_filt)))] <- 0

#### explore the distribution of rho values ----
ccm_rho_filt |>
  as_tibble() |>
  pivot_longer(cols = everything(), values_to = 'rho', names_to = 'variables') |>
  dplyr::filter(!is.na(rho)) |>
  dplyr::mutate(type = 'ccm_results_rho') |>
  ggplot(aes(x = rho, y = type))+
  geom_density_ridges()+
  theme_ridges()

ccm_rho_filt |>
  as_tibble() |>
  pivot_longer(cols = everything(), values_to = 'rho', names_to = 'variables') |>
  dplyr::filter(!is.na(rho)) |>
  dplyr::summarise(mean = mean(rho),
                   sd = sd(rho))

# Set a threshold for the maximum number of allowed NA values
na_threshold <- 100

# Identify columns where the number of NA values is less than or equal to the threshold
cols_below_threshold <- apply(ccm_rho_filt , 2, function(col) sum(is.na(col)) <= na_threshold)

# Remove columns and rows where the number of NA values exceeds the threshold
matrix_clean <- ccm_rho_filt[, cols_below_threshold]
matrix_clean <- matrix_clean[cols_below_threshold, ]

matrix_clean |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  dplyr::filter(rho > 0.4) |> 
  ggplot(aes(variables, variables_2))+
  geom_tile(aes(fill = rho))+
  scale_fill_gradientn(colors = palete_gradient_cb2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.border = element_blank())

### 3. DATASET: We did NOT remove seasonality & data is in the format of rCLR -------
#### upload data 
asv_tab_bbmo_10y_rclr_occ_filt_inter <- read.csv( '../EDM_carmen/marbits/asv_tab_bbmo_10y_rclr_occ_filt_inter.csv', sep = ',') 
asv_tab_bbmo_10y_rclr_occ_filt_inter  <- asv_tab_bbmo_10y_rclr_occ_filt_inter  |>
  dplyr::select(-'X')
ccm_rho <- read.csv( '../EDM_carmen/marbits/ccm_rclr_seas/ccm_rho_Bl_rclr_seas.csv', sep = ',') 
ccm_sig <- read.csv( '../EDM_carmen/marbits/ccm_rclr_seas/ccm_sig_Bl_rclr_seas.csv', sep = ',') 

ccm_rho |>
  colnames() <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_sig |>
  colnames() <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

##filter rho only by the significant ones 
ccm_rho |>
  row.names() <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_rho <- ccm_rho |>
  as.matrix()

ccm_rho_filt <- ifelse(ccm_sig == 1, ccm_rho, NA)

ccm_rho_filt |>
  row.names() <- asv_tab_bbmo_10y_rclr_occ_filt_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

# replace cells where row name equals column name with 0
ccm_rho_filt[cbind(match(row.names(ccm_rho_filt)[row.names(ccm_rho_filt) %in% colnames(ccm_rho_filt)], row.names(ccm_rho_filt)), 
                   match(row.names(ccm_rho_filt)[row.names(ccm_rho_filt) %in% colnames(ccm_rho_filt)], colnames(ccm_rho_filt)))] <- 0

#### explore the distribution of rho values ----
ccm_rho_filt |>
  as_tibble() |>
  pivot_longer(cols = everything(), values_to = 'rho', names_to = 'variables') |>
  dplyr::filter(!is.na(rho)) |>
  dplyr::mutate(type = 'ccm_results_rho') |>
  ggplot(aes(x = rho, y = type))+
  geom_density_ridges()+
  theme_ridges()

ccm_rho_filt |>
  as_tibble() |>
  pivot_longer(cols = everything(), values_to = 'rho', names_to = 'variables') |>
  dplyr::filter(!is.na(rho)) |>
  dplyr::summarise(mean = mean(rho),
                   sd = sd(rho))

# Set a threshold for the maximum number of allowed NA values
na_threshold <- 100

# Identify columns where the number of NA values is less than or equal to the threshold
cols_below_threshold <- apply(ccm_rho_filt , 2, function(col) sum(is.na(col)) <= na_threshold)

# Remove columns where the number of NA values exceeds the threshold
matrix_clean <- ccm_rho_filt [, cols_below_threshold]
matrix_clean <- matrix_clean[cols_below_threshold, ]

dim(matrix_clean)

matrix_clean |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  dplyr::filter(rho > 0.5) |> 
  ggplot(aes(variables, variables_2))+
  geom_tile(aes(fill = rho))+
  scale_fill_gradientn(colors = palete_gradient_cb2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.border = element_blank())

### 4. DATASET: We did not removed seasonality & data is in the format of relative abundance -----
asv_tab_10y_rel_occ_filt_w_env_inter <- read.csv( '../EDM_carmen/marbits/asv_tab_10y_rel_occ_filt_w_env_inter.csv', sep = ',') 
asv_tab_10y_rel_occ_filt_w_env_inter <- asv_tab_10y_rel_occ_filt_w_env_inter  |>
  dplyr::select(-'X')
colnames_vec <- asv_tab_10y_rel_occ_filt_w_env_inter |>
  colnames()

ccm_rho <- read.csv( '../EDM_carmen/marbits/ccm_rel_abund_seas/ccm_rho_Bl_rel_abund_seas.csv', sep = ',') 
ccm_sig <- read.csv( '../EDM_carmen/marbits/ccm_rel_abund_seas/ccm_sig_Bl_rel_abund_seas.csv', sep = ',') 

ccm_rho |>
  colnames() <- asv_tab_10y_rel_occ_filt_w_env_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_sig |>
  colnames() <- asv_tab_10y_rel_occ_filt_w_env_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

##filter rho only by the significant ones 
ccm_rho |>
  row.names() <- asv_tab_10y_rel_occ_filt_w_env_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

ccm_rho <- ccm_rho |>
  as.matrix()

ccm_rho_filt <- ifelse(ccm_sig == 1, ccm_rho, NA)

ccm_rho_filt |>
  row.names() <- asv_tab_10y_rel_occ_filt_w_env_inter |>
  dplyr::select(-'Year', -'Month', -Day) |>
  colnames()

# replace cells where row name equals column name with 0
ccm_rho_filt[cbind(match(row.names(ccm_rho_filt)[row.names(ccm_rho_filt) %in% colnames(ccm_rho_filt)], row.names(ccm_rho_filt)), 
                   match(row.names(ccm_rho_filt)[row.names(ccm_rho_filt) %in% colnames(ccm_rho_filt)], colnames(ccm_rho_filt)))] <- 0

#### explore the distribution of rho values ----
ccm_rho_filt |>
  as_tibble() |>
  pivot_longer(cols = everything(), values_to = 'rho', names_to = 'variables') |>
  dplyr::filter(!is.na(rho)) |>
  dplyr::mutate(type = 'ccm_results_rho') |>
  ggplot(aes(x = rho, y = type))+
  geom_density_ridges()+
  theme_ridges()

ccm_rho_filt |>
  as_tibble() |>
  pivot_longer(cols = everything(), values_to = 'rho', names_to = 'variables') |>
  dplyr::filter(!is.na(rho)) |>
  dplyr::summarise(mean = mean(rho),
                   sd = sd(rho))

# Set a threshold for the maximum number of allowed NA values
na_threshold <- 100

# Identify columns where the number of NA values is less than or equal to the threshold
cols_below_threshold <- apply(ccm_rho_filt , 2, function(col) sum(is.na(col)) <= na_threshold)

# Remove columns where the number of NA values exceeds the threshold
matrix_clean <- ccm_rho_filt [, cols_below_threshold]
matrix_clean <- matrix_clean [cols_below_threshold, ]

matrix_clean |>
  as.matrix() |>
  heatmap()

matrix_clean 

matrix_clean |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  dplyr::filter(rho > 0.5) |> 
  ggplot(aes(variables, variables_2))+
  geom_tile(aes(fill = rho))+
  scale_fill_gradientn(colors = palete_gradient_cb2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.border = element_blank())

#### Exploration of the best way to visualize the results and what do they explain -----
##### play with some values: change the threshold of the NAs allowed per row or column 
##### upload taxonomic information ----
tax_occ_filt_bbmo <- read_csv('../EDM_carmen/tax_occ_filt_bbmo_ed.csv', col_names = T) |>
  dplyr::select(-'...1') |>
  dplyr::select( 'domain' = 'Kingdom',
 'phylum' = 'Phylum', 
                'class' = 'Class', 'order' = 'Order', 'family' = 'Family', 
                'genus' = 'Genus', asv_num, bloom ) |>
  dplyr::mutate(bloom = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ '1',
                                  asv_num %in% c('asv8', 'asv2', 'asv3', 'asv5') ~ '0', ## not real bloomers
                                  TRUE ~ '0'))

tax_occ_filt_bbmo |>
  dplyr::filter(bloom == '1') %$%
  unique(asv_num) ## bloomers that i can study in these dataset ("asv15" "asv17" "asv23" "asv1"  "asv7"  "asv31" "asv22" "asv11")

##### upload frequency information ----
bloo_all_types_summary_tb <- read.csv('results/tables/bloo_all_types_summary_tb_tax_v2.csv')

bloo_all_types_summary_tb_tax <- bloo_all_types_summary_tb |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(fraction = as.character(fraction))

bloo_all_types_summary_tb <- bloo_all_types_summary_tb |>
  dplyr::select(asv_num, recurrency, fraction) |>
  dplyr::mutate(fraction = as.character(fraction))

##### prepare data for plotting ----
### I work with the matrix clean which has only significant rho values and also filtered rows and cols by a threshold based on the number of NAs accepted ---
# Replace NA values (e.g., with 0)
matrix_clean |>
  dim()

ccm_rho_filt |>
  dim()

matrix_clean[is.na(matrix_clean)] <- 0

# Create the heatmap
matrix_clean |>
  as.matrix() |>
  pheatmap() #color = palete_gradient_cb2

matrix_clean |>
  dim()

## pheatmap edited to understand easily the data being visualized ----
matrix_clean |>
  dim()

## for the datasets coming from data transformed to rel abundance I need to replace NA rho values to 0 to be able to plot the phaetmap
data <- matrix_clean |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |> ## variables 2 are causal variables
  separate(col = variables, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                   asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                
                order = case_when(is.na(order) ~ 'env',
                                   TRUE ~ order), 
                
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                
                asv_num = case_when(asv_num == 'synechococcus' ~ 'Synechococcus',
                                    asv_num == 'joint' ~ 'Bacterial Abundance',
                                    variables == 'day_length' ~ 'Day length',
                                    asv_num == 'temperature' ~ 'Temperature',
                                    variables == 'chla_total' ~ 'Chl-a',
                                    variables == 'BP_FC1.55' ~ 'Bacterial Production',
                                    variables == 'PNF_Micro' ~ 'PNF',
                                    variables == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                                    variables == "PNF_5um_Micro" ~ 'PNF 5 um',
                                    asv_num == 'cryptomonas' ~ 'Cryptomonas',
                                    asv_num == 'micromonas' ~ 'Micromonas',
                                    variables == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                                    variables == 'HNF_Micro' ~ 'HNF',
                                    variables == 'HNF_5um_Micro' ~ 'HNF 5 um',
                                    TRUE ~ asv_num)
  ) |>
  left_join(bloo_all_types_summary_tb) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 5 um', 'PNF 2-5 um', 'Cryptomonas', 'Micromonas', 
                   'HFN 2-5 um', 'HNF', 'HFN 5 um', "PO4" , "NH4" ,  
                   "NO2",  "NO3", "Si") ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  ))

# add annotations
anotr1 <- data |>
  distinct(variables, fraction, asv_num, phylum, class, order, family, bloom, frequency)

anotr2 <- data |>
  dplyr::select(-c(variables, 'fraction', 'asv_num', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'bloom', 'frequency', 'recurrency')) |>
  separate(col = variables_2, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  distinct(variables_2, fraction, asv_num, phylum, class, order, family, bloom) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                
                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order), 
                
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                
                asv_num = case_when(
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  variables_2 == 'bacteria_joint' ~ 'Bacterial Abundance',
                  variables_2 == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables_2 == 'chla_total' ~ 'Chl-a',
                  variables_2 == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables_2 == 'PNF_Micro' ~ 'PNF',
                  variables_2 == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables_2 == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables_2 == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables_2 == 'HNF_Micro' ~ 'HNF',
                  variables_2 == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num) 
  )  |>
  dplyr::left_join(bloo_all_types_summary_tb, by = c('asv_num', 'fraction')) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 2-5 um', 'PNF 5 um', 'Cryptomonas', 'Micromonas', 'HFN 2-5 um',
                   'HNF', 'HNF 5 um', "PO4" , "NH4" ,  "NO2",  "NO3", "Si"
                   ) ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  ))

#####
anotr1 |>
  row.names() |>
  length() == anotr2 |>
  row.names() |>
  length() 

anotr_all <- anotr1 |>
  dplyr::select(variables_2 = variables,'fraction', 'asv_num', 'phylum', 'class', 'order', 'family', 'bloom' , frequency) |>
  bind_rows(anotr2)

all.equal(rownames(matrix_clean), anotr1$variables) # true?

annotr <- data.frame(fraction = as.factor(anotr1$fraction),
                     order = as.factor(anotr1$order),
                     bloom = as.factor(anotr2$frequency)) |>
  set_rownames(anotr1$variables)

annotr |>
  dim()

# Define gaps between row groups
all.equal(colnames(matrix_clean), anotr2$variables_2) # true?

annotc <- data.frame(fraction = as.factor(anotr2$fraction),
                     order = as.factor(anotr2$order),
                     bloom = as.factor(anotr2$frequency)) %>%
  set_rownames(anotr2$variables_2)

# Check if colors match anotr1 columns
all(names(colors) %in% colnames(annotr))

colors <- list(
  bloom = c('Non-Seasonal' = 'black', 'No Bloomer' = "white", 'Seasonal' = '#F7CDCD', 'env' = '#4A785C'),
  order = palette_order_assigned_all,
  fraction = c('0.2' = '#00808F', '3' = '#454545', 'env' = '#4A785C'))

## filter colors by only those classes present in the matrix
unique_class_values <- unique(annotr$class)
unique_order_values <- unique(annotr$order)

# Filter the colors list based on the unique annotation values
filtered_colors <- list(
  bloom = c('Non-Seasonal' = 'black', 'No Bloomer' = "white", 'Seasonal' = '#F7CDCD', 'env' = '#4A785C'),
  order = colors$order[names(colors$order) %in% unique_order_values],
  fraction = c('0.2' = '#00808F', '3' = '#3D4C4D', 'env' = '#4A785C')
)

# plot it
set.seed(123)

matrix_clean |>
  dim()

matrix_clean |>
  colnames()

str(matrix_clean)

# Check the structure of the annotation and label data
str(anotr1)

# Check that the labels match the dimensions of matrix_clean
annotr$fraction |>
  unique()

annotc$fraction |>
  unique()

# Verify label lengths match matrix dimensions
if (length(anotr2$asv_num) == nrow(matrix_clean) && length(anotr1$asv_num) == ncol(matrix_clean)) {
  heatmap <- pheatmap(matrix_clean,
                      annotation_col = annotc,
                      annotation_colors = filtered_colors,
                      labels_row = anotr1$asv_num,
                      labels_col = anotr2$asv_num
  )
} else {
  stop("Length of labels does not match matrix dimensions")
}

## I reorder columns by fraction and environmental variables so that i can observe the tendencies better
# Extract the fraction column from the row and column annotations
annotr$fraction <- factor(annotr$fraction, levels = c('3', '0.2', 'env'))
fraction_row_order <- order(annotr$fraction)

# Reorder the matrix_clean based on the fraction order
matrix_clean_reordered <- matrix_clean[fraction_row_order,]

# Define gaps between row groups
gaps_row <- c(62)  # positions after which to insert gaps (separate fraction and then environmental data)

## Remove environmental data from cols (we would like to only cluster the ASVs) ----
annotc_filt <- annotc |>
  dplyr::filter(fraction != 'env') |>
  row.names() |>
  as_tibble_col(column_name = 'variables')

matrix_clean_f <- matrix_clean[, annotc_filt$variables]

annotc_f <- annotc |>
  dplyr::filter(fraction != 'env')

anotr2_f <- anotr2 |>
  dplyr::filter(fraction != 'env')

# Order the rows of matrix_clean based on the fraction factor
fraction_row_order <- order(annotr$fraction)
matrix_clean_f_reor <- matrix_clean_f[fraction_row_order,]

# Reorder the annotr data frame and row_labels 
annotr_reordered <- annotr[fraction_row_order,]
anotr1_reordered <- anotr1[fraction_row_order,] # I reorder label cols in the same way that i have ordered the annotr_reordered

## Check that the rows and columns annotation and labels match ----
# Capture the comparison result for annotr_reordered
result_annotr <- all(row.names(annotr_reordered) == anotr1_reordered$variables)

# Check if the result is FALSE and print a warning
if (!result_annotr) {
  warning("The row names of annotr_reordered do not match anotr1_reordered$variables. Please reorder rows.")
}
print(result_annotr)

# Capture the comparison result for annotc_f
result_annotc <- all(row.names(annotc_f) == anotr2_f$variables_2)
print(result_annotc)

# Check if the result is FALSE and print a warning
if (!result_annotc) {
  warning("The row names of annotc_f do not match anotr2_f$variables_2. Please reorder columns")
}
print(result_annotc)

# Define gaps between row groups
gaps_row <- c(sum(annotr_reordered$fraction == '3'), 
              sum(annotr_reordered$fraction == '3') + sum(annotr_reordered$fraction == '0.2'))

matrix_clean_f_reor <- ifelse(is.na(matrix_clean_f_reor), 0, matrix_clean_f_reor) ## for matrices that come from rel abund data

heatmap_noclustcols <-
  pheatmap(matrix_clean_f_reor,
           annotation_row = annotr_reordered,
           annotation_col = annotc_f,
           annotation_colors = filtered_colors,
           labels_row = anotr1_reordered$asv_num,
           labels_col = anotr2_f$asv_num,
           cluster_cols = T,
           cluster_row = F,
           gaps_row = gaps_row,
           cellwidth = 10, 
           cellheight = 10,
           border_color = 'white',
           fontsize = 10
  )
       # ggsave( plot = heatmap_noclustcols, filename = 'heatmap_ccm_seas_rel_60na_v2.pdf',
       #         path = 'results/figures/',
       #         width = 400, height = 480, units = 'mm')

##### ONLY RELATIONSHIP WITH MY BLOOMERS AND THE COMMUNITY OR ENVIRONMENTAL DATA -----
## pheatmap edited to understand easily the data being visualized ----

ccm_rho_filt |>
  dim()  ## i work with this matrix without filtering by number of NAs in each column/row.
### Notes from Carmen, matrix rho is a not symetric matrix. The effect of i on j is not the same
### as j on i. Columns are causal variables of rows.

## for the datasets coming from data transformed to rel abundance I need to replace NA rho values to 0 to be able to plot the phaetmap

## variables_2 are the causal variables
## variables are the resulting variables. (I need these filtered by my bloomers, the ones that I'm interested in)

bloo_all_types_summary_tb <- bloo_all_types_summary_tb |>
  dplyr::mutate(fraction = as.character(fraction)) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') )

data <- ccm_rho_filt |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  separate(col = variables, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order), 
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                asv_num = case_when(
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  asv_num == 'joint' ~ 'Bacterial Abundance',
                  variables == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables == 'chla_total' ~ 'Chl-a',
                  variables == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables == 'PNF_Micro' ~ 'PNF',
                  variables == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables == 'HNF_Micro' ~ 'HNF',
                  variables == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num)
  ) |>
  left_join(bloo_all_types_summary_tb) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 5 um', 'PNF 2-5 um', 'Cryptomonas', 'Micromonas', 
                   'HFN 2-5 um', 'HNF 5 um','HNF', 'HFN 5 um', "PO4" , "NH4" ,  
                   "NO2",  "NO3", "Si") ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  ))

# add annotations
annotr_labs <- data |>
  distinct(variables, fraction, asv_num, phylum, class, order, family, bloom, frequency)

anotr <- annotr_labs 

annotc_labs <- data |>
  dplyr::select(-c(variables, 'fraction', 'asv_num', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'bloom', 'frequency', 'recurrency')) |>
  separate(col = variables_2, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  distinct(variables_2, fraction, asv_num, phylum, class, order, family, bloom) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                
                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order), 
                
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                
                asv_num = case_when(
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  variables_2 == 'bacteria_joint' ~ 'Bacterial Abundance',
                  variables_2 == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables_2 == 'chla_total' ~ 'Chl-a',
                  variables_2 == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables_2 == 'PNF_Micro' ~ 'PNF',
                  variables_2 == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables_2 == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables_2 == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables_2 == 'HNF_Micro' ~ 'HNF',
                  variables_2 == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num) 
  )  |>
  dplyr::left_join(bloo_all_types_summary_tb, by = c('asv_num', 'fraction')) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 2-5 um', 'PNF 5 um', 'Cryptomonas', 'Micromonas', 'HFN 2-5 um',
                   'HNF', 'HNF 5 um', "PO4" , "NH4" ,  "NO2",  "NO3", "Si"
    ) ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  ))

anotc <- annotc_labs 

#####
annotr_labs |>
  row.names() |>
  length() == anotc |>
  row.names() |>
  length() 

annotr |>
  length() == annotc |>
  length()

anotr_all <- annotr_labs |>
  dplyr::select(variables_2 = variables,'fraction', 'asv_num', 'phylum', 'class', 'order', 'family', 'bloom' , frequency) |>
  bind_rows(annotc_labs)

all.equal(rownames(ccm_rho_filt), annotr_labs$variables) # true?

annotr <- data.frame(fraction = as.factor(anotr$fraction),
                     order = as.factor(anotr$order),
                     bloom = as.factor(anotr$frequency)) |>
  set_rownames(anotr$variables)

all.equal(colnames(ccm_sig), anotc$variables_2) # true?

annotc <- data.frame(fraction = as.factor(anotc$fraction),
                     order = as.factor(anotc$order),
                     bloom = as.factor(anotc$frequency)) %>%
  set_rownames(anotc$variables_2)

# Check if colors match anotr1 columns
all(names(colors) %in% colnames(annotr))

colors <- list(
  bloom = c('Non-Seasonal' = 'black', 'No Bloomer' = "white", 'Seasonal' = '#F7CDCD', 'env' = '#4A785C'),
  order = palette_order_assigned_all,
  fraction = c('0.2' = '#00808F', '3' = '#454545', 'env' = '#4A785C'))

## filter colors by only those classes present in the matrix
unique_class_values <- unique(annotc$class)
unique_order_values <- unique(annotc$order)

# Filter the colors list based on the unique annotation values
filtered_colors <- list(
  bloom = c('Non-Seasonal' = 'black', 'No Bloomer' = "white", 'Seasonal' = '#F7CDCD', 'env' = '#4A785C'),
  order = colors$order[names(colors$order) %in% unique_order_values],
  fraction = c('0.2' = '#00808F', '3' = '#3D4C4D', 'env' = '#4A785C')
)

# plot it
set.seed(123)

## I reorder columns by fraction and environmental variables so that i can observe the tendencies better
# Extract the fraction column from the row and column annotations
annotc$fraction <- factor(annotc$fraction, levels = c('3', '0.2', 'env'))
fraction_column_order <- order(annotc$fraction)

annotc_labs <- annotc_labs[fraction_column_order,] ## reorder the labels the same way than the annotation
annotc_reordered <- annotc[fraction_column_order,]

## Remove environmental data and other ASVs that are not BLOOMERS from rows (we would like to only cluster the bloomers ASVs) ----
annotr_filt <- annotr |>
  dplyr::filter(fraction != 'env') |>
  dplyr::filter(bloom %in% c('Seasonal', 'Non-Seasonal')) |>
  row.names() |>
  as_tibble_col(column_name = 'variables') ## dataset to filter the matrix

ccm_rho_filt_f <- ccm_rho_filt[annotr_filt$variables, ] # we filter the matrix

annotr_labs_f <- annotr_labs |>
  dplyr::filter(fraction != 'env') |>
  dplyr::filter(frequency %in% c('Seasonal', 'Non-Seasonal'))

annotr_f <- annotr |>
  dplyr::filter(fraction != 'env') |>
  dplyr::filter(bloom %in% c('Seasonal', 'Non-Seasonal'))

fraction_row_order <- order(annotr_f$fraction)

ccm_rho_filt_f_reordered <- ccm_rho_filt_f[,fraction_column_order]
ccm_rho_filt_f_reordered  <- ccm_rho_filt_f_reordered[fraction_row_order,]

annotr_labs_f <- annotr_labs_f[fraction_row_order,]
annotr_f_reorder <- annotr_f[fraction_row_order,]

## Check that the rows and columns annotation and labels match ----
# Capture the comparison result for annotr_reordered
result_annotr <- all(row.names(annotr_f_reorder) == annotr_labs_f$variables)

# Check if the result is FALSE and print a warning
if (!result_annotr) {
  warning("The row names of annotr_reordered do not match annotr_labs_filt$variables. Please reorder rows.")
}
print(result_annotr)

# Capture the comparison result for annotc_f
result_annotc <- all(row.names(annotc_reordered) == annotc_labs$variables_2)

# Check if the result is FALSE and print a warning
if (!result_annotc) {
  warning("The row names of annotc do not match anotr2_f$variables_2. Please reorder columns")
}
print(result_annotc)

# Define gaps between row groups
gaps_col <- c(sum(annotc_reordered$fraction == '3'), 
              sum(annotc_reordered$fraction == '3') + sum(annotc_reordered$fraction == '0.2'))

#### COMMENT THIS LINE WHEN IS DATA RCLR TRANSFORMED!!! Or when it plots the heatmap without complaining because of NAs
#ccm_rho_filt_reordered_f <- ifelse(is.na(ccm_rho_filt_reordered_f), 0, ccm_rho_filt_reordered_f) ## for matrices that come from rel abund data

annotc_labs$variables_2 == annotc_reordered |>
  row.names()

ccm_rho_filt_f_reordered |>
  colnames() == annotc_reordered |>
  row.names()

annotc_labs$variables_2 == ccm_rho_filt_f_reordered |>
  colnames() 

heatmap_noclustcols <-
  pheatmap(ccm_rho_filt_f_reordered,
           color = palette_gradient(5),
           annotation_row = annotr_reordered,
           annotation_col = annotc_reordered,
           annotation_colors = filtered_colors,
           labels_row = annotr_labs_f$asv_num,
           labels_col = annotc_labs$asv_num,
           cluster_cols = F,
           cluster_row = T,
           cellwidth = 10, 
           cellheight = 10,
           border_color = 'white',
           fontsize = 10
           #legend_breaks = c("0", "0.2", "0.4", "0.6", '0.8', '1')
           #filename = 'results/figures/heatmap_ccm.pdf',
           #height = 200,
           #width = 180
  )

# ggsave( plot = heatmap_noclustcols, filename = 'heatmap_ccm_bloo_seas_rel2_v2.pdf',
#         path = 'results/figures/ccm/',
#         width = 480, height = 380, units = 'mm')

# ggsave( plot = heatmap_noclustcols, filename = 'heatmap_ccm_bloo_noseas_rclr.svg',
#         path = 'results/figures/poster_svg_format/',
#         width = 480, height = 380, units = 'mm')

# BOXPLOTS of causality intra e inter taxonomic class ----
## Columns are the causal variables on rows. 
ccm_rho_filt ## dataset with the matrix filtered for only the Rho significant values 

data <- ccm_rho_filt  |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  separate(col = variables, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                
                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order), 
                
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                
                asv_num = case_when(
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  asv_num == 'joint' ~ 'Bacterial Abundance',
                  variables == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables == 'chla_total' ~ 'Chl-a',
                  variables == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables == 'PNF_Micro' ~ 'PNF',
                  variables == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables == 'HNF_Micro' ~ 'HNF',
                  variables == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num)
  ) |>
  left_join(bloo_all_types_summary_tb) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 5 um', 'PNF 2-5 um', 'Cryptomonas', 'Micromonas', 
                   'HFN 2-5 um', 'HNF', 'HFN 5 um', "PO4" , "NH4" ,  
                   "NO2",  "NO3", "Si") ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  ))

data2 <- ccm_rho_filt  |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  separate(col = variables_2, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                
                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order), 
                
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                
                asv_num = case_when(
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  asv_num == 'joint' ~ 'Bacterial Abundance',
                  variables == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables == 'chla_total' ~ 'Chl-a',
                  variables == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables == 'PNF_Micro' ~ 'PNF',
                  variables == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables == 'HNF_Micro' ~ 'HNF',
                  variables == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num)
  ) |>
  left_join(bloo_all_types_summary_tb) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 5 um', 'PNF 2-5 um', 'Cryptomonas', 'Micromonas', 
                   'HFN 2-5 um', 'HNF', 'HFN 5 um', "PO4" , "NH4" ,  
                   "NO2",  "NO3", "Si") ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  ))

## intra taxonomic class (within the same taxonomic class) ----
data |>
  distinct(class)

data |>
  colnames()

data_red <- data |>
  dplyr::select(variables, variables_2, fraction, asv_num, frequency, class, order, family, rho)

data2_red <- data2 |>
  dplyr::select(variables1 = variables, variables_2_2 = variables_2, fraction2 = fraction, asv_num2 = asv_num, 
                frequency2 = frequency, class2 = class, order2 = order, family2 = family)

data_intrainter_class <- data_red |>
  bind_cols(data2_red)

## check that the bind cols is correct
length(data_intrainter_class$variables) == length(data_intrainter_class$variables1) 
comparison <- data_intrainter_class$variables == data_intrainter_class$variables1 
summary_result <- table(comparison) # Summarize the counts of TRUE and FALSE
print(summary_result) # Print the summary result

## filter causal relationships inside the same taxonomic class 
data_intrainter_class |>
  colnames() 

intra_causal_relationships_plot <- data_intrainter_class |>
  dplyr::filter(class != 'env') |>
  dplyr::filter(variables != variables_2) |> # filter relationships between the same ASVs
  dplyr::filter(class == class2) |> # relationships inside the same class
  ggplot(aes(x = class, y = rho))+
  geom_point(aes(color = order, shape = fraction), position = position_jitter(width = 0.1))+
  geom_violin(alpha = 0.1)+
  scale_color_manual(values = palette_order_assigned_all)+
  scale_shape_discrete(labels = labs_fraction)+
  labs(x='Class', y = 'Rho Intra Taxonomic Classes', color = 'Order', shape = 'Fraction')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        panel.grid = element_blank(), text = element_text(size = 10),
        axis.text.x = element_text(size = 7),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm')
        )+
  guides(color = guide_legend(ncol = 1),
         alpha = 'none',
         shape = guide_legend(ncol = 1))

# ggsave( plot = intra_causal_relationships_plot, filename = 'intra_causal_relationships_plot_noseas_rclr.pdf',
#         path = 'results/figures/ccm/',
#         width = 180, height = 160, units = 'mm')

## inter taxonomic class (with other taxonomic classes) -----
intrer_causal_relationships_plot <- data_intrainter_class |>
  dplyr::filter(class != 'env') |>
  dplyr::filter(variables != variables_2) |> # filter relationships between the same ASVs
  dplyr::filter(class != class2)  |>  # relationships between different classes
  dplyr::filter(class2 != 'env') |>
  dplyr::mutate(interclass = paste0(class2, '.', class)) |>
  group_by(interclass) |>
  dplyr::filter(n() > 10) |> # filter by only those classes that had more than x significant causal relationships
  ungroup() |>
  dplyr::filter(!is.na(rho)) |>
  ggplot(aes(x = interclass, y = rho))+
  geom_point(aes(color = order, shape = fraction), position = position_jitter(width = 0.1))+
  geom_violin(alpha = 0.1)+
  coord_flip()+
  scale_color_manual(values = palette_order_assigned_all)+
  labs(x='Class', y = 'Rho Class causal relationship with the other Class', color = 'Order', shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend(ncol = 3),
         alpha = 'none',
         shape = guide_legend(ncol = 1))

intrer_causal_relationships_plot

# ggsave( plot = intrer_causal_relationships_plot, filename = 'inter_causal_relationships_plot_noseas_rclr.pdf',
#         path = 'results/figures/ccm/',
#         width = 180, height = 160, units = 'mm')

# TOP causal relationships for my BLOOMERS ------
## variables are the variables being predicted by variables_2 (columns from initial matrix) which are the causal variables.
data <- ccm_rho_filt  |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  separate(col = variables, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                
                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order), 
                
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                
                asv_num = case_when
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  asv_num == 'joint' ~ 'Bacterial Abundance',
                  variables == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables == 'chla_total' ~ 'Chl-a',
                  variables == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables == 'PNF_Micro' ~ 'PNF',
                  variables == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables == 'HNF_Micro' ~ 'HNF',
                  variables == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num)
  ) |>
  left_join(bloo_all_types_summary_tb) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 5 um', 'PNF 2-5 um', 'Cryptomonas', 'Micromonas', 
                   'HFN 2-5 um', 'HNF', 'HFN 5 um', "PO4" , "NH4" ,  
                   "NO2",  "NO3", "Si") ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  )) |>
  dplyr::filter(asv_num %in% bloo_02_filt$value & fraction == '0.2' |
                  asv_num %in% bloo_3$value & fraction == '3' )
  #dplyr::filter(frequency != 'No Bloomer') # I would only like to explore the causal relationships of my bloomers 

## This dataset is the one with the causal variables, (variables_2: columns in the matrix), i would like to keep only those that have an effect ON my bloomers
# Combine the patterns into a single regex pattern using '|' (OR) to detect any of the patterns
###patterns <- paste(bloo_taxonomy$asv_num, collapse = "|")
data2 <- ccm_rho_filt  |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables'), names_to = 'variables_2', values_to = 'rho') |>
  separate(col = variables_2, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),

                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),

                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order),

                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),

                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),

                asv_num = case_when(
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  asv_num == 'joint' ~ 'Bacterial Abundance',
                  variables == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables == 'chla_total' ~ 'Chl-a',
                  variables == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables == 'PNF_Micro' ~ 'PNF',
                  variables == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables == 'HNF_Micro' ~ 'HNF',
                  variables == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num)
  ) |>
  left_join(bloo_all_types_summary_tb) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length',
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 5 um', 'PNF 2-5 um', 'Cryptomonas', 'Micromonas',
                   'HFN 2-5 um', 'HNF', 'HFN 5 um', "PO4" , "NH4" ,
                   "NO2",  "NO3", "Si") ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  )) |>
  #dplyr::filter(str_detect(variables, patterns)) # I would only like to explore the causal relationships of my bloomers
  separate(col = variables, into = c('fraction_rows', 'asv_num_rows'), remove = F) |>
  dplyr::mutate(asv_num_rows = case_when(is.na(asv_num_rows) ~ fraction_rows,
                                    asv_num_rows == '5um' ~ fraction_rows,
                                    TRUE ~ asv_num_rows)) |>
  dplyr::mutate(fraction_rows = case_when(fraction_rows == 'bp' ~ '0.2',
                                     fraction_rows == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::filter(asv_num_rows %in% bloo_02_filt$value & fraction_rows == '0.2' |
                  asv_num_rows %in% bloo_3$value & fraction_rows == '3' )

# ## Causal relationships dataset to see the effect of variables on my bloomers (CAUSAL RELATIONSHIP) ----
# data_red <- data |>
#   dplyr::select(variables, variables_2, fraction, asv_num, frequency, phylum, class, order, family, rho)
# 
# data2_red <- data2 |>
#   dplyr::select(variables1 = variables, variables_2_2 = variables_2, fraction2 = fraction, asv_num2 = asv_num, 
#                 frequency2 = frequency, phylum2 = phylum, class2 = class, order2 = order, family2 = family)
# 
# data_causal_bloomers <- data_red |>
#   bind_cols(data2_red)
# 
# ## check that the bind cols is correct
# length(data_causal_bloomers$variables) == length(data_causal_bloomers$variables1) 
# comparison <- data_causal_bloomers$variables == data_causal_bloomers$variables1 
# summary_result <- table(comparison) # Summarize the counts of TRUE and FALSE
# print(summary_result) # Print the summary result
# 
# data_causal_bloomers |>
#   distinct(asv_num)
# 
# data_causal_bloomers |>
#   distinct(asv_num2)
# 
# data_causal_bloomers <- data_causal_bloomers |>
#   dplyr::mutate(phylum_f = as_factor(phylum2),
#                 family_f = as_factor(family2),
#                 order_f = as_factor(order2),
#                 class_f = as_factor(class2),
#                 asv_num_f = as_factor(asv_num2))
# 
# data_causal_bloomers$class_f <-  factor(data_causal_bloomers$class_f, 
#                                           levels=unique(data_causal_bloomers$class_f[order(data_causal_bloomers$phylum_f)]), 
#                                           ordered=TRUE)
# 
# data_causal_bloomers$order_f <-  factor(data_causal_bloomers$order_f, 
#                                           levels=unique(data_causal_bloomers$order_f[order(data_causal_bloomers$phylum_f,
#                                                                                              data_causal_bloomers$class_f)]), 
#                                           ordered=TRUE)
# 
# data_causal_bloomers$family_f <-  factor(data_causal_bloomers$family_f, 
#                                            levels=unique(data_causal_bloomers$family_f[order(data_causal_bloomers$phylum_f,
#                                                                                                data_causal_bloomers$class_f,
#                                                                                                data_causal_bloomers$order_f)]), 
#                                            ordered=TRUE)
# 
# data_causal_bloomers$asv_num_f <-  factor(data_causal_bloomers$asv_num_f, 
#                                             levels=unique(data_causal_bloomers$asv_num_f[order(data_causal_bloomers$phylum_f,
#                                                                                                  data_causal_bloomers$class_f,
#                                                                                                  data_causal_bloomers$order_f,
#                                                                                                  data_causal_bloomers$family_f)]), 
#                                             ordered=TRUE)
# 
# data_causal_bloomers <- data_causal_bloomers |>
#   dplyr::mutate(phylum_f_bloo = as_factor(phylum),
#                 family_f_bloo = as_factor(family),
#                 order_f_bloo = as_factor(order),
#                 class_f_bloo = as_factor(class),
#                 asv_num_f_bloo = as_factor(asv_num))
# 
# data_causal_bloomers$class_f_bloo <-  factor(data_causal_bloomers$class_f_bloo, 
#                                         levels=unique(data_causal_bloomers$class_f_bloo[order(data_causal_bloomers$phylum_f_bloo)]), 
#                                         ordered=TRUE)
# 
# data_causal_bloomers$order_f_bloo <-  factor(data_causal_bloomers$order_f_bloo, 
#                                         levels=unique(data_causal_bloomers$order_f_bloo[order(data_causal_bloomers$phylum_f_bloo,
#                                                                                          data_causal_bloomers$class_f_bloo)]), 
#                                         ordered=TRUE)
# 
# data_causal_bloomers$family_f_bloo <-  factor(data_causal_bloomers$family_f_bloo, 
#                                          levels=unique(data_causal_bloomers$family_f_bloo[order(data_causal_bloomers$phylum_f_bloo,
#                                                                                            data_causal_bloomers$class_f_bloo,
#                                                                                            data_causal_bloomers$order_f_bloo)]), 
#                                          ordered=TRUE)
# 
# data_causal_bloomers$asv_num_f_bloo <-  factor(data_causal_bloomers$asv_num_f_bloo, 
#                                           levels=unique(data_causal_bloomers$asv_num_f_bloo[order(data_causal_bloomers$phylum_f_bloo,
#                                                                                              data_causal_bloomers$class_f_bloo,
#                                                                                              data_causal_bloomers$order_f_bloo,
#                                                                                              data_causal_bloomers$family_f_bloo)]), 
#                                           ordered=TRUE)
# 
# top10_causal_bloomers_plot <- data_causal_bloomers |>
#   dplyr::filter(asv_num != asv_num2) |> # filter relationships between the same
#   group_by(asv_num, fraction) |>
#   dplyr::filter(rho > 0.3) |>
#   slice_max(order_by = rho, n = 10, na_rm = T) |>
#   ungroup() |>
#   ggplot(aes(asv_num_f_bloo, interaction(family_f, asv_num_f)))+
#   labs(x = 'Bloomers', y = 'Causal variables', fill = 'Rho')+
#   geom_tile(aes(fill = rho))+
#   scale_fill_gradientn(colours = palete_gradient_cb2)+
#   facet_wrap(vars(fraction), labeller = labs_fraction, dir = 'v')+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         text = element_text(size = 10),
#         axis.text.x = element_text(size = 6))
#   
# top10_causal_bloomers_plot

# ggsave( plot = top10_causal_bloomers_plot, filename = 'top10_causal_bloomers_plot_noseas_rclr.pdf',
#         path = 'results/figures/ccm/',
#         width = 180, height = 160, units = 'mm')

## Evaluation of the role of the environment and the community on bloomers (seasonal and chaotic) 
## and other ASVs (non-bloomers) -------
### boxplot 

### Columns are causal variables of rows.

## for the datasets coming from data transformed to rel abundance I need to replace NA rho values to 0 to be able to plot the phaetmap
## variables_2 are the causal variables
## variables are the resulting variables. (I need these filtered by my bloomers, the ones that I'm interested in)
bloo_all_types_summary_tb <- bloo_all_types_summary_tb |>
  dplyr::mutate(fraction = as.character(fraction)) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') )

data <- ccm_rho_filt |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = 'variables_predicted') |>
  as_tibble() |>
  pivot_longer(cols = !c('variables_predicted'), names_to = 'variables_causal', values_to = 'rho') |>
  separate(col = variables_predicted, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order), 
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                asv_num = case_when(
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  asv_num == 'joint' ~ 'Bacterial Abundance',
                  variables_predicted == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables_predicted == 'chla_total' ~ 'Chl-a',
                  variables_predicted == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables_predicted == 'PNF_Micro' ~ 'PNF',
                  variables_predicted == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables_predicted == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables_predicted == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables_predicted == 'HNF_Micro' ~ 'HNF',
                  variables_predicted == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num)
  ) |>
  left_join(bloo_all_types_summary_tb) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 5 um', 'PNF 2-5 um', 'Cryptomonas', 'Micromonas', 
                   'HFN 2-5 um', 'HNF 5 um','HNF', 'HFN 5 um', "PO4" , "NH4" ,  
                   "NO2",  "NO3", "Si") ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  ))

data$variables_predicted |>
 unique()

data$variables_causal == data$variables_predicted

data_ed <- data |>
  dplyr::select(variables_predicted, asv_num_predicted = asv_num, fraction_predicted = fraction, phylum_predicted = phylum,
                domain_predicted = domain, class_predicted = class, order_predicted = order, family_predicted = family, 
                frequency_predicted = frequency, variables_causal, rho) |>
  separate(col = variables_causal, into = c('fraction', 'asv_num'), remove = F) |>
  dplyr::mutate(asv_num = case_when(is.na(asv_num) ~ fraction,
                                    asv_num == '5um' ~ fraction,
                                    TRUE ~ asv_num)) |>
  dplyr::left_join(tax_occ_filt_bbmo, by = 'asv_num') |>
  dplyr::mutate(fraction = case_when(fraction == 'bp' ~ '0.2',
                                     fraction == 'bn' ~ '3',
                                     TRUE ~ 'env')) |>
  dplyr::mutate(phylum = case_when(is.na(phylum) ~ 'env',
                                   TRUE ~ phylum),
                class = case_when(is.na(class) ~ 'env',
                                  TRUE ~ class),
                order = case_when(is.na(order) ~ 'env',
                                  TRUE ~ order), 
                family = case_when(is.na(family) ~ 'env',
                                   TRUE ~ family),  
                bloom = case_when(is.na(bloom) ~ 'env',
                                  TRUE ~ bloom),
                asv_num = case_when(
                  asv_num == 'synechococcus' ~ 'Synechococcus',
                  asv_num == 'joint' ~ 'Bacterial Abundance',
                  variables_causal == 'day_length' ~ 'Day length',
                  asv_num == 'temperature' ~ 'Temperature',
                  variables_causal == 'chla_total' ~ 'Chl-a',
                  variables_causal == 'BP_FC1.55' ~ 'Bacterial Production',
                  variables_causal == 'PNF_Micro' ~ 'PNF',
                  variables_causal == 'PNF2_5um_Micro' ~ 'PNF 2-5 um',
                  variables_causal == "PNF_5um_Micro" ~ 'PNF 5 um',
                  asv_num == 'cryptomonas' ~ 'Cryptomonas',
                  asv_num == 'micromonas' ~ 'Micromonas',
                  variables_causal == 'HNF2_5um_Micro' ~ 'HFN 2-5 um',
                  variables_causal == 'HNF_Micro' ~ 'HNF',
                  variables_causal == 'HNF_5um_Micro' ~ 'HNF 5 um',
                  TRUE ~ asv_num)
  ) |>
  left_join(bloo_all_types_summary_tb) |>
  dplyr::mutate(frequency = case_when(
    recurrency == 'non-recurrent' ~ 'Non-Seasonal',
    recurrency == 'recurrent' ~ 'Seasonal',
    asv_num %in% c('Synechococcus', 'Bacterial Abundance', 'Day length', 
                   'Temperature', 'Chl-a', 'Bacterial Production', 'PNF',
                   'PNF 5 um', 'PNF 2-5 um', 'Cryptomonas', 'Micromonas', 
                   'HFN 2-5 um', 'HNF 5 um','HNF', 'HFN 5 um', "PO4" , "NH4" ,  
                   "NO2",  "NO3", "Si") ~ 'env',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  )) |>
  dplyr::select(variables_predicted, asv_num_predicted, fraction_predicted, phylum_predicted ,
                domain_predicted, class_predicted, order_predicted , family_predicted , 
                frequency_predicted , variables_causal, rho,
                 asv_num_causal = asv_num, fraction_causal = fraction, phylum_causal = phylum,
                domain_causal = domain, class_causal = class, order_causal = order, family_causal = family, 
                frequency_causal = frequency)

data_ed$frequency_predicted |>
  unique()

data_ed$frequency_causal |>
  unique()

labs_fraction_env <- as_labeller(c('0.2' = 'Free living\n(0.2-3 um)',
                               '3' = 'Particle attached\n(3-20 um)',
                               'env' = 'Environmental\nvariables'))

labs_fraction_bloo <- as_labeller(c('0.2' = 'Free living\n(0.2-3 um)',
                                   '3' = 'Particle attached\n(3-20 um)',
                                   'bloomer' = 'Bloomer',
                                   'no_bloomer' = 'No Bloomer'))

## EFFECT OF THE COMMUNITY ON BLOOMERS AND ON THE COMMUNTY -----
data_ed |>
  colnames()

select_bloomers <- data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::distinct(asv_num_predicted, frequency_predicted, fraction_predicted) |>
  dplyr::filter(frequency_predicted != 'No Bloomer')

## add statistics (by group)
### anova to compare mean of the different groups 
data_ed2 <- data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::mutate(fraction_bloom = paste0(fraction_predicted,bloom))

data_ed2 |>
  str()

data_ed2 |>
  dplyr::filter(bloom == 'no_bloomer') %$%
  unique(asv_num_predicted)

library(car)  ## homocedasticity 
library(dunn.test) ## significant groups in Kruskal Test
library(ggpubr) ## ggqqplot

# I need to perform an anova for each group in case that each group has a normal distribution and that they have homocedasticity ---- 
## I group: 0.2 no bloomers 
data_ed2_f <- data_ed2 |>
  dplyr::filter(fraction_bloom == '0.2no_bloomer')

## check normality 
shapiro.test(as.numeric(data_ed2 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(data_ed2 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed2_f) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  177.26 < 2.2e-16 ***
#       4085                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test

# Perform one-way ANOVA to compare the means of Rho across the three groups
## anova_result <- aov(rho ~ fraction_causal, data = data_ed2_f)  # Replace `bloom` with your actual group variable
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed2_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed2_f$rho,                            # The values you are comparing
  g = as.factor(data_ed2_f$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)

# Col Mean-|
#   Row Mean |        0.2          3
# -------------------------------
#   3 |  -0.144353
# |     1.0000
# |
#   env |   15.66140   13.97643
# |    0.0000*    0.0000*
  
## II group: 0.2 no bloomers -----
data_ed2_f <- data_ed2 |>
  dplyr::filter(fraction_bloom == '0.2bloomer')

data_ed2_f |>
  colnames()

## check normality 
shapiro.test(as.numeric(data_ed2 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(data_ed2 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed2_f) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  177.26 < 2.2e-16 ***
#       4085                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed2_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed2_f$rho,                            # The values you are comparing
  g = as.factor(data_ed2_f$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)

# # Perform one-way ANOVA to compare the means of Rho across the three groups
# anova_result <- aov(rho ~ fraction_causal, data = data_ed2_f)  # Replace `bloom` with your actual group variable
# 
# # View the ANOVA summary
# summary(anova_result)
# 
# # If ANOVA shows a significant difference, perform post-hoc tests to determine which groups differ
# TukeyHSD(anova_result)

data_ed2$fraction_bloom |>
  unique()

## III group: 3no_bloomer -----
data_ed2_f <- data_ed2 |>
  dplyr::filter(fraction_bloom == '3no_bloomer')

data_ed2_f |>
  colnames()

## check normality 
shapiro.test(as.numeric(data_ed2 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(data_ed2 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed2_f) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  177.26 < 2.2e-16 ***
#       4085                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed2_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed2_f$rho,                            # The values you are comparing
  g = as.factor(data_ed2_f$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)

# # Perform one-way ANOVA to compare the means of Rho across the three groups
# anova_result <- aov(rho ~ fraction_causal, data = data_ed2_f)  # Replace `bloom` with your actual group variable
# 
# # View the ANOVA summary
# summary(anova_result)
# 
# # If ANOVA shows a significant difference, perform post-hoc tests to determine which groups differ
# TukeyHSD(anova_result)

## IV group: 3bloomer -----
data_ed2_f <- data_ed2 |>
  dplyr::filter(fraction_bloom == '3bloomer')

data_ed2_f |>
  colnames()

## check normality 
shapiro.test(as.numeric(data_ed2 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(data_ed2 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed2_f) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  177.26 < 2.2e-16 ***
#       4085                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed2_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed2_f$rho,                            # The values you are comparing
  g = as.factor(data_ed2_f$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)

# # Perform one-way ANOVA to compare the means of Rho across the three groups
# anova_result <- aov(rho ~ fraction_causal, data = data_ed2_f)  # Replace `bloom` with your actual group variable
# 
# # View the ANOVA summary
# summary(anova_result)
# 
# # If ANOVA shows a significant difference, perform post-hoc tests to determine which groups differ
# TukeyHSD(anova_result)

## without color (non informative in this case) ----
causal_effects_on_community_plot <- data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  ggplot(aes(x = interaction(fraction_causal), y = rho))+
  geom_point(aes(
    alpha = 0.6),position = position_jitter(width = 0.1))+
  geom_boxplot(alpha = 0.5, notch = T)+
  facet_wrap(fraction_predicted~bloom, labeller = labs_fraction, dir = 'v'_bloo)+
  scale_x_discrete(labels = labs_fraction_env)+
  theme_bw()+
  labs(y = 'Prediction (Rho) ', x= 'Causal effects'#, color = 'Order predicted ASVs'
       )+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid.major.y  = element_blank(),
        panel.grid.minor  = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(alpha = 'none',
         color = guide_legend(ncol = 4))

causal_effects_on_community_plot
# 
# ggsave(plot = causal_effects_on_community_plot, filename = 'causal_effects_on_community_plot_v2.pdf',
#                 path = 'results/figures/',
#                 width = 180, height = 180, units = 'mm')

## statistics (separate bloomers into those that are seasonal and non-seasonal) ----
data_ed2_f <- data_ed2 |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::filter(bloom != 'no_bloomer') |>
  dplyr::filter(frequency_predicted != 'No Bloomer')

## group I -----
data_ed2_f1 <- data_ed2_f |>
  dplyr::filter(fraction_predicted == '0.2' &
                bloom == 'bloomer' &
                frequency_predicted == 'Non-Seasonal')

## check normality 
shapiro.test(as.numeric(data_ed2_f1 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value = 0.1119 (NORMALITY)
ggqqplot(as.numeric(data_ed2_f1 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

shapiro.test(as.numeric(data_ed2_f1 |>
                          dplyr::filter(fraction_causal == '3') %$%
                          rho)) # => p-value = 0.836 (NORMALITY)
ggqqplot(as.numeric(data_ed2_f1 |>
                      dplyr::filter(fraction_causal == '3') %$%
                      rho))

shapiro.test(as.numeric(data_ed2_f1 |>
                          dplyr::filter(fraction_causal == 'env') %$%
                          rho)) # => p-value = 0.1229 (NORMALITY)
ggqqplot(as.numeric(data_ed2_f1 |>
                      dplyr::filter(fraction_causal == 'env') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed2_f1) #0.1711 p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Perform one-way ANOVA to compare the means of Rho across the three groups
 anova_result <- aov(rho ~ fraction_causal, data = data_ed2_f1)  # Replace `bloom` with your actual group variable

# View the ANOVA summary
summary(anova_result)
#summary(kruskal_result)

TukeyHSD(anova_result)

## group II ----
data_ed2_f2 <- data_ed2_f |>
  dplyr::filter(fraction_predicted == '0.2' &
                  bloom == 'bloomer' &
                  frequency_predicted == 'Seasonal')

## check normality 
shapiro.test(as.numeric(data_ed2_f2 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value = 8.073e-05 (NO NORMALITY)
ggqqplot(as.numeric(data_ed2_f2 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed2_f2) #0.002456 p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed2_f2) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed2_f2$rho,                            # The values you are comparing
  g = as.factor(data_ed2_f2$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)

## group III ----
data_ed2_f3 <- data_ed2_f |>
  dplyr::filter(fraction_predicted == '3' &
                  bloom == 'bloomer' &
                  frequency_predicted == 'Non-Seasonal')

## check normality 
shapiro.test(as.numeric(data_ed2_f3 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value = 0.001159 (NO NORMALITY)

ggqqplot(as.numeric(data_ed2_f2 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed2_f3) #0.0001744 p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed2_f3) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed2_f3$rho,                            # The values you are comparing
  g = as.factor(data_ed2_f3$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)

## Group IV ----
data_ed2_f4 <- data_ed2_f |>
  dplyr::filter(fraction_predicted == '3' &
                  bloom == 'bloomer' &
                  frequency_predicted == 'Seasonal')

## check normality 
shapiro.test(as.numeric(data_ed2_f4 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value = 0.2179 (NORMALITY)

ggqqplot(as.numeric(data_ed2_f4 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

shapiro.test(as.numeric(data_ed2_f4 |>
                          dplyr::filter(fraction_causal == '3') %$%
                          rho)) # => p-value = 6.122e-06 (NO NORMALITY)


## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed2_f4) #0.5131p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed2_f4) 

summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed2_f4$rho,                            # The values you are comparing
  g = as.factor(data_ed2_f4$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)

## separate bloomers by frequency ----
causal_effects_on_community_plot <- data_ed2 |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::filter(bloom != 'no_bloomer') |>
  dplyr::filter(frequency_predicted != 'No Bloomer') |>
  ggplot(aes(x = interaction(fraction_causal), y = rho))+
  geom_point(aes( color = order_predicted, alpha = 0.8),position = position_jitter(width = 0.1))+
  geom_boxplot(alpha = 0.1, notch = T)+
  facet_wrap(fraction_predicted~bloom~frequency_predicted)+ #, labeller = labs_fraction, dir = 'v'_bloo
  scale_color_manual(values = palette_order_assigned_all)+
  scale_x_discrete(labels = labs_fraction_env)+
  theme_bw()+
  labs(y = 'Rho', x= 'Causal effects', color = 'Order predicted ASVs')+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid.major.y  = element_blank(),
        panel.grid.minor  = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(alpha = 'none',
         color = guide_legend(ncol = 4))

causal_effects_on_community_plot

# ggsave(plot = causal_effects_on_community_plot, filename = 'causal_effects_on_community_plot_recurrency_v2.pdf',
#                 path = 'results/figures/',
#                 width = 180, height = 180, units = 'mm')

### BLOOMERS EFFECT ON THE COMMUNITY -------
### statistics ---
data_stats <- data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_causal %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_causal %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::filter(bloom == 'bloomer') 

## group I
data_stats_f1 <- data_stats |>
  dplyr::filter(fraction_predicted == '0.2')

## check normality
shapiro.test(as.numeric(data_stats_f1 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value = 3.455e-06 (NO NORMALITY)

ggqqplot(as.numeric(data_ed2_f1 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))
## Welch test ----
group1 <- data_stats_f1$rho[data_stats_f1$fraction_causal == "3"]
group2 <- data_stats_f1$rho[data_stats_f1$fraction_causal == "0.2"]

# Perform Welch's t-test
welch_result <- t.test(group1, group2, var.equal = FALSE)

# View the result
welch_result
## t = -0.36271, df = 386.21, p-value = 0.717

## group II ---
## group I
data_stats_f2 <- data_stats |>
  dplyr::filter(fraction_predicted == '3')

## check normality
shapiro.test(as.numeric(data_stats_f2 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # =>p-value = 0.01061 (NO NORMALITY)

ggqqplot(as.numeric(data_ed2_f2 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))
## Welch test ----
group1 <- data_stats_f2$rho[data_stats_f2$fraction_causal == "3"]
group2 <- data_stats_f2$rho[data_stats_f2$fraction_causal == "0.2"]

# Perform Welch's t-test
welch_result <- t.test(group1, group2, var.equal = FALSE)

# View the result
welch_result
## t = 13.128, df = 226.2, p-value < 2.2e-16

### plot ----
bloomers_effect_on_the_community_ccm_based_plot <- data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_causal %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_causal %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::filter(bloom == 'bloomer') |>
  ggplot(aes(x = interaction(fraction_causal), y = rho))+
  geom_point(aes( color = order_causal),position = position_jitter(width = 0.1))+
  geom_boxplot(alpha = 0.1, notch = T)+
  facet_wrap(fraction_predicted~bloom, labeller = labs_fraction, dir = 'v'_bloo, ncol = 1)+
  scale_color_manual(values = palette_order_assigned_all)+
  scale_x_discrete(labels = labs_fraction_env)+
  guides(color = guide_legend(ncol = 2))+
  theme_bw()+
  labs(y = 'Prediction (Rho)', x= 'Causal effects', color = 'Order')+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid.major.y  = element_blank(),
        panel.grid.minor  = element_blank(), text = element_text(size = 7),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 7),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

bloomers_effect_on_the_community_ccm_based_plot

# ggsave(plot = bloomers_effect_on_the_community_ccm_based_plot, filename = 'bloomers_effect_on_the_community_ccm_based_plot_v2.pdf',
#        path = 'results/figures/',
#        width = 88, height = 120, units = 'mm')

### ARE BLOOMERS KEYSTONE TAXA? -----
### They need to have strong causal interactions with other taxa 
bloomers_effect_on_the_community_ccm_based_plot <- data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_causal %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_causal %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::filter(bloom == 'bloomer') |>
  ggplot(aes(x = interaction(fraction_causal), y = rho))+
  geom_point(aes( color = order_causal),position = position_jitter(width = 0.1))+
  geom_boxplot(alpha = 0.1, notch = T)+
  facet_wrap(fraction_predicted~bloom, labeller = labs_fraction, dir = 'v'_bloo, ncol = 1)+
  scale_color_manual(values = palette_order_assigned_all)+
  scale_x_discrete(labels = labs_fraction_env)+
  guides(color = guide_legend(ncol = 2))+
  theme_bw()+
  labs(y = 'Prediction (Rho)', x= 'Causal effects', color = 'Order')+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid.major.y  = element_blank(),
        panel.grid.minor  = element_blank(), text = element_text(size = 7),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 7),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

relative_rho_value <- data_ed |>
  dplyr::filter(frequency_causal != 'env') |>
  dplyr::filter(fraction_causal != 'env') |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |> 
  dplyr::filter(rho >0) |>
  dplyr::group_by(#fraction_causal, 
    asv_num_causal, #fraction_predicted, 
    order_causal, family_causal) |>
  dplyr::reframe(rho = sum(rho),
               n = n()) |>
  dplyr::mutate(relative_rho = rho/n)

bloo_taxonomy <- bloo_taxonomy |>
  dplyr::filter(!asv_num_f %in% c('asv2', 'asv3', 'asv5', 'asv8')) 

data_causal_relationships <- data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::filter(frequency_causal != 'env') |>
  dplyr::filter(fraction_causal != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_causal %in%  unique(select_bloomers$asv_num_predicted) ~ 'no_bloomer',
                                  asv_num_causal %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |> 
  dplyr::filter(rho >0) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_causal' = 'asv_num')) |>
  dplyr::mutate(family_asv_num = case_when(!is.na(genus) ~ paste0(genus, ' ', asv_num_causal),
                                           is.na(genus) ~ paste0(family, ' ', asv_num_causal))) |>
  dplyr::mutate(order_family_asv_num = case_when(!is.na(genus) ~ paste0(order, ' ', family, ' ', 
                                                                         genus, ' ', asv_num_causal),
                                                 is.na(genus) ~ paste0(order, ' ',
                                                                                       family, ' ', asv_num_causal))) |>
  dplyr::mutate(order_family_asv_num = case_when(asv_num_causal %in% bloo_taxonomy$asv_num_f ~ paste0(order_family_asv_num,'*'),
                                                 TRUE ~ order_family_asv_num)) 

relative_rho_value |>
  colnames()

plot_ccm_asvs_causally_affecting_other_asvs_ordered <- relative_rho_value  |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_causal' = 'asv_num')) |>
  dplyr::mutate(family_asv_num = case_when(!is.na(genus) ~ paste0(genus, ' ', asv_num_causal),
                                           is.na(genus) ~ paste0(family, ' ', asv_num_causal))) |>
  dplyr::mutate(order_family_asv_num = case_when(!is.na(genus) ~ paste0(order, ' ', family, ' ', 
                                                                        genus, ' ', asv_num_causal),
                                                 is.na(genus) ~ paste0(order, ' ',
                                                                       family, ' ', asv_num_causal))) |>
  dplyr::mutate(order_family_asv_num = case_when(asv_num_causal %in% bloo_taxonomy$asv_num_f ~ paste0(order_family_asv_num,'*'),
                                                 TRUE ~ order_family_asv_num)) |>
  ggplot(aes(fct_reorder(order_family_asv_num, relative_rho), relative_rho))+
  geom_point(data = data_causal_relationships, aes(order_family_asv_num, rho
             #,color = fraction_predicted
             ,
             alpha = if_else((asv_num_causal %in% bloo_taxonomy$asv_num_f), 0.6, 0.1)), position = position_jitter(width = 0.25),  size = 0.6)+
  geom_boxplot(data = data_causal_relationships, aes(order_family_asv_num, rho,
                                                     alpha = if_else((asv_num_causal %in% bloo_taxonomy$asv_num_f), 0.8, 0.3)),
               notch = T, outliers = F)+ # I already plot the with the geom point
  geom_point(shape = 18)+
  coord_flip()+
  labs(y = 'Rho', x= 'ASVs causally effectting other ASVs', color = 'Fraction ASVs being causally effected')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        #panel.grid.major.y  = element_blank(),
        panel.grid.minor  = element_blank(), text = element_text(size = 5),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 7),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.2, "mm"))+
  guides(alpha = 'none')

plot_ccm_asvs_causally_affecting_other_asvs_ordered

# ggsave(plot = plot_ccm_asvs_causally_affecting_other_asvs_ordered, 
#         filename = 'plot_ccm_asvs_causally_affecting_other_asvs_ordered.pdf',
#        path = 'results/figures/',
#        width = 130, height = 150, units = 'mm')
  
  occurrence_bloo_bbmo |>
  colnames()

occurrence_bloo_bbmo_sim <- occurrence_bloo_bbmo |>
  distinct(asv_num, fraction, occurrence_perc) # I have only the occurrrence for my bloomers

data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::filter(frequency_causal != 'env') |>
  dplyr::filter(fraction_causal != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_causal %in%  unique(select_bloomers$asv_num_predicted) ~ 'no_bloomer',
                                  asv_num_causal %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |> 
  dplyr::filter(rho >0) |>
  left_join(relative_rho_value) |> 
  ggplot(aes(fct_reorder(as.factor(asv_num_causal), desc(relative_rho)), rho))+
  geom_point(aes(color = fraction_causal), position = position_jitter(width = 0.25))+
  geom_boxplot(notch = F, alpha = 0.7)+
  facet_wrap(vars(bloom))+
  coord_flip()+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid.major.y  = element_blank(),
        panel.grid.minor  = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(alpha = 'none',
         color = guide_legend(ncol = 4))

# -------------- #### MDR RESULTS #### ----------------------------------------- 

## Important definitions: MDR: multiview distance regularized S-map
## Jacobian interaction: quantifies the net effects of abundance changes in nodes j (between two consecutive observations) on the abundance of species i, 
## which is more consistent with the findings of empirical addition/removal experiments
## The interaction Jacobian is time-varying as it integrates two sources of temporal variability: (i) species abundance changed with time; 
## and (ii) the interaction coefficient is time-varying (i.e., αij  =  αij(t);
## upload data ----
mdr_tb <- read.csv2('../EDM_carmen/MDR/Bl_nin120_cvunit0.025_aenet_jcof_Nmvx_Rallx_demo_v2.csv', header = T ) |>
  as_tibble() ## new data

mdr_tb |>
  dim()

mdr_tb |>
  distinct(variable)

bloo_all_types_summary_tb <- read.csv('results/tables/bloo_all_types_summary_tb_tax_v2.csv')

bloo_all_types_summary_tb_tax <- bloo_all_types_summary_tb |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(fraction = as.character(fraction))

bloo_all_types_summary_tb <- bloo_all_types_summary_tb |>
  dplyr::select(asv_num, recurrency, fraction) |>
  dplyr::mutate(fraction = as.character(fraction))

# Bloomers are affected by other ASVs  -----
## prepare data for plotting  ----
### The tibble is structured such that each row corresponds to a variable that receives an interaction, and each column corresponds to a variable that performs (or realizes) the interaction.
### columns ASVs that affect other ASVs, whereas rows are the ones that receive the interaction
variable_num <- mdr_tb |>
  colnames() |>
  as_tibble_col(column_name = 'variable_name') |>
  dplyr::filter(!variable_name %in% c('Insample', 'time', 'variable', 'j0'))

variable_num_tb <- mdr_tb |>
  colnames() |>
  as_tibble_col(column_name = 'variable_name') |>
  dplyr::filter(!variable_name %in% c('Insample', 'time', 'variable', 'j0')) |>
  dplyr::mutate(variable_num = 1:nrow(variable_num))

m_02_red2 <- m_02 |>
  dplyr::select(sample_id, date, fraction, sample_id_num) |>
  dplyr::select(-sample_id)

m_3_red2 <- m_02 |>
  dplyr::select(sample_id, date, fraction, sample_id_num) |>
  dplyr::select(-sample_id)

time_num_02 <- mdr_tb %$%
  unique(time) |>
  as_tibble_col(column_name = 'time_num') |>
  dplyr::mutate(time_num = as.character(time_num)) |>
  left_join(m_02_red2, by = c('time_num' = 'sample_id_num'))

time_num_3 <- mdr_tb %$%
  unique(time) |>
  as_tibble_col(column_name = 'time_num') |>
  dplyr::mutate(time_num = as.character(time_num)) |>
  left_join(m_3_red2, by = c('time_num' = 'sample_id_num')) |>
  dplyr::mutate(fraction = '3')

time_num <- time_num_02 |>
  bind_rows(time_num_3)

mdr_tb_m <- mdr_tb |>
  left_join(variable_num_tb , by = c('variable' = 'variable_num')) |>
  dplyr::mutate(fraction = case_when(str_detect(variable_name, 'bp') ~ '0.2',
                                     str_detect(variable_name, 'bn') ~ '3')) |>
  separate(variable_name, sep = '_', into = c('fraction_b', 'asv_num')) |>
  left_join(bloo_all_types_summary_tb, by = c('asv_num', 'fraction')) |>
  dplyr::mutate(recurrency = case_when(
    recurrency == 'non.recurrent' ~ 'Chaotic',
    recurrency == 'recurrent' ~ 'Seasonal',
    is.na(recurrency) ~ 'No Bloomer',
    TRUE ~ recurrency  # Preserve other values if any
  )) |>
  dplyr::mutate(time = as.character(time)) |>
  left_join(time_num, by = c('fraction', 'time' = 'time_num'))

mdr_tb_m |>
  colnames()

time_num |>
  dplyr::filter(is.na(date))

mdr_tb_m |>
  dplyr::filter(is.na(date))

## how do these interactions change in time (small pheatmaps to have an overview)
matrix_mdr_t1 <- mdr_tb |>
  dplyr::filter(time == '1') |>
  dplyr::select(-Insample, -time, -j0) |>
  left_join(variable_num_tb, by = c('variable' = 'variable_num')) |>
  dplyr::select(-variable) |>
  column_to_rownames(var = 'variable_name') |>
  tibble::as_tibble(rownames = "row_name") |> 
  rowwise() |> 
  dplyr::mutate(across(-row_name, ~ if_else(row_name == cur_column(), NA_real_, as.numeric(.)))) |> 
  ungroup() |> 
  column_to_rownames(var = "row_name") |>
  dplyr::mutate(across(everything(), as.numeric)) |>
  as.matrix() 

mdr_t1 <- pheatmap(matrix_mdr_t1, 
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = 47,
         cellwidth = 10, 
         cellheight = 10,
         border_color = 'white',
         color = palette_gradient,
         breaks = breaks,
         fontsize = 10)

# ggsave( plot = mdr_t1, filename = 'mdr_t1.pdf',
#         path = 'results/figures/MDR/v2/',
#         width = 350, height = 350, units = 'mm')

## loop to observe how did the interactions changed over time ---

# Define a color palette and breaks for the scale
palette_gradient <- c('#0A6260',"#BDBEBE","#545454",
                      "#F2A218")

palette_gradient <- colorRampPalette(c('#0A6260',"#BDBEBE", '#ffffff', "#545454",
                                       "#F2A218"))(10)

breaks <- seq(-1, 1, length.out = 10) # Ensure breaks match the number of colors + 1

# List of unique time points
time_points <- unique(mdr_tb$time)

# Loop over time points
for (t in time_points) {
  
  # Filter data for the current time point and preprocess
  matrix_mdr <- mdr_tb |> 
    dplyr::filter(time == t) |> 
    dplyr::select(-Insample, -time, -j0) |> 
    left_join(variable_num_tb, by = c('variable' = 'variable_num')) |> 
    dplyr::select(-variable) |> 
    column_to_rownames(var = 'variable_name') |> 
    tibble::as_tibble(rownames = "row_name") |> 
    rowwise() |> 
    dplyr::mutate(across(-row_name, ~ if_else(row_name == cur_column(), NA_real_, as.numeric(.)))) |> 
    ungroup() |> 
    column_to_rownames(var = "row_name") |> 
    dplyr::mutate(across(everything(), as.numeric)) |> 
    as.matrix()
  
  annotation_row_tb <- mdr_tb |> 
    dplyr::filter(time == t) |> 
    dplyr::select(-Insample, -time, -j0) |> 
    left_join(variable_num_tb, by = c('variable' = 'variable_num')) |> 
    dplyr::select(-variable) |> 
    column_to_rownames(var = 'variable_name') |> 
    tibble::as_tibble(rownames = "row_name") |> 
    rowwise() |> 
    dplyr::mutate(across(-row_name, ~ if_else(row_name == cur_column(), NA_real_, as.numeric(.)))) |> 
    ungroup() |> 
    column_to_rownames(var = "row_name") |> 
    dplyr::mutate(across(everything(), as.numeric)) |> 
    as.matrix() |>
    row.names() |>
    as_tibble_col(column_name = 'variable_name') |>
    dplyr::mutate(annotation_rows = str_replace(variable_name, 'bp_', '')) |>
    dplyr::mutate(annotation_rows = str_replace(annotation_rows, 'bn_', '')) |>
    column_to_rownames(var = 'variable_name')
  
  annotation_col_tb <- mdr_tb |> 
    dplyr::filter(time == t) |> 
    dplyr::select(-Insample, -time, -j0) |> 
    left_join(variable_num_tb, by = c('variable' = 'variable_num')) |> 
    dplyr::select(-variable) |> 
    column_to_rownames(var = 'variable_name') |> 
    tibble::as_tibble(rownames = "row_name") |> 
    rowwise() |> 
    dplyr::mutate(across(-row_name, ~ if_else(row_name == cur_column(), NA_real_, as.numeric(.)))) |> 
    ungroup() |> 
    column_to_rownames(var = "row_name") |> 
    dplyr::mutate(across(everything(), as.numeric)) |> 
    as.matrix() |>
    colnames() |>
    as_tibble_col(column_name = 'annotation_cols') |>
    dplyr::mutate(annotation_cols_ed = str_replace(annotation_cols, 'bp_', '')) |>
    dplyr::mutate(annotation_cols_ed = str_replace(annotation_cols_ed, 'bn_', '')) 
  
  annotation_row_tb <- annotation_row_tb[!is.na(annotation_row_tb$annotation_rows), ]
  annotation_col_tb <- annotation_col_tb[!is.na(annotation_col_tb$annotation_cols_ed), ]
  
   # Check row names
  all(rownames(matrix_mdr) %in% rownames(annotation_row_tb)) # Should return TRUE
  all(rownames(annotation_row_tb) %in% rownames(matrix_mdr)) # Should return TRUE
  
  # Check column names
  all(colnames(matrix_mdr) %in% annotation_col_tb$annotation_cols) # Should return TRUE
  all(annotation_col_tb$annotation_cols %in% colnames(matrix_mdr))
  
  title = time_num |>
    dplyr::filter(time_num == t ) |>
    distinct(date)
  
  # Create heatmap
  heatmap_plot <- pheatmap(
    matrix_mdr, 
    color = palette_gradient,      # Use the defined color palette
    breaks = breaks,            # Use the fixed breaks
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = 47,
    gaps_col = 47,
    cellwidth = 10, 
    cellheight = 10,
    border_color = 'white',
    fontsize = 10,
    main = title$date
  )
  
  # Save heatmap to file

  
  )
}

## prepare bloomers abundance data ----
asv_tab_all_bloo_z_tax_red <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(asv_num, fraction, date, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 
                               'asv23', 'asv31', 'asv15', 'asv17')) 

## plot interactions of the community on my bloomers (the effect that the community has on bloomers)-----
### important, the community is affecting each ASV but PA and FL are afecting the same ASV PA and FL. 4 different types of interactions!!
### ASV11 (0.2) ----
mdr_tb_m |>
  colnames()

#### I need to group positive and negative interactions for each timepoint 
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv11') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv11' & fraction == '0.2')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                     str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv11_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on ASV11 (0.2)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv11_affectedby_community_plot 

asv11_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_line()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction, dir = 'v') + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.75), units = "lines"))

asv11_rclr_plot

# Create an empty plot for the placeholder
#empty_plot <- ggplot() + theme_void()

# Arrange the plots with an empty spot on the right side for the first plot
effect_community_on_ASV11 <- plot_grid(
# First row with empty space
  asv11_rclr_plot,
  asv11_affectedby_community_plot, # Second row
  rel_heights = c(1, 1.75), # Heights of the rows
  labels = c('A', 'B'), # Labels for the plots
  label_size = 10,
  label_fontface = 'plain',
  ncol = 1 # Arrange in one column
)

ggsave(
  plot = effect_community_on_ASV11, 
  filename = 'effect_community_on_ASV11.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 150, 
  units = 'mm'
)

### ASV23 ----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv23') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv23')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv23_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(fraction~fraction_eff_ed, labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv23 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv23_affectedby_community_plot 

asv_tab_all_bloo_z_tax_red_asv |>
  colnames()

asv23_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_line()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

effect_community_on_ASV23 <- plot_grid(asv23_rclr_plot, 
          asv23_affectedby_community_plot,
          rel_heights = c(1, 2),
          labels = c('A', 'B'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

effect_community_on_ASV23

ggsave(
  plot = effect_community_on_ASV23, 
  filename = 'effect_community_on_ASV23.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 180, 
  units = 'mm'
)
 
### ASV7 ----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv7') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv7')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv7_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(fraction~fraction_eff_ed, labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv7 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv7_affectedby_community_plot 

asv_tab_all_bloo_z_tax_red_asv |>
  colnames()

asv7_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_line()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv7_rclr_plot

effect_community_on_ASV7 <- plot_grid(asv7_rclr_plot, 
          asv7_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('A', 'B'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

ggsave(
  plot = effect_community_on_ASV7, 
  filename = 'effect_community_on_ASV7.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 180, 
  units = 'mm'
)

### ASV1 ----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv1') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv1')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv1_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(fraction~fraction_eff_ed, labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv1 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv1_affectedby_community_plot 

asv_tab_all_bloo_z_tax_red_asv |>
  colnames()

asv1_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_line()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv1_rclr_plot

effect_community_on_ASV1 <- plot_grid(asv1_rclr_plot, 
          asv1_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('A', 'B'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

ggsave(
  plot = effect_community_on_ASV1, 
  filename = 'effect_community_on_ASV1.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 180, 
  units = 'mm'
)

### ASV22 -----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv22') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv22')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv22_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(fraction~fraction_eff_ed, labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv22 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv22_affectedby_community_plot 

asv_tab_all_bloo_z_tax_red_asv |>
  colnames()

asv22_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  dplyr::filter(fraction == '3') |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_line()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv22_rclr_plot

effect_community_on_ASV22 <- plot_grid(asv22_rclr_plot, 
          asv22_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('A', 'B'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

ggsave(
  plot = effect_community_on_ASV22, 
  filename = 'effect_community_on_ASV22.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 180, 
  units = 'mm'
)

### ASV15 -----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv15') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv15')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv15_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(fraction~fraction_eff_ed, labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv15 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv15_affectedby_community_plot 

asv_tab_all_bloo_z_tax_red_asv |>
  colnames()

asv15_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_line()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv15_rclr_plot

effect_community_on_ASV15 <- plot_grid(asv15_rclr_plot, 
          asv15_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('A', 'B'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

 ggsave(
  plot = effect_community_on_ASV15, 
  filename = 'effect_community_on_ASV15.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 180, 
  units = 'mm'
)

### ASV31 -----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv31') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv31')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv31_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(fraction~fraction_eff_ed, labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv31 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv31_affectedby_community_plot 

asv_tab_all_bloo_z_tax_red_asv |>
  colnames()

asv31_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  dplyr::filter(fraction == '3') |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_line()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction, dir = 'v') + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv31_rclr_plot

effect_community_on_ASV31 <- plot_grid(asv31_rclr_plot, 
          asv31_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('A', 'B'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

ggsave(
  plot = effect_community_on_ASV31, 
  filename = 'effect_community_on_ASV31.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 180, 
  units = 'mm'
)

### ASV17 -----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv17') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv17')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv17_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(fraction~fraction_eff_ed, labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv17 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv17_affectedby_community_plot 

asv_tab_all_bloo_z_tax_red_asv |>
  colnames()

asv17_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_line()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv17_rclr_plot

effect_community_on_ASV17 <- plot_grid(asv17_rclr_plot, 
          asv17_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('A', 'B'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

ggsave(
  plot = effect_community_on_ASV17, 
  filename = 'effect_community_on_ASV17.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 180, 
  units = 'mm'
)

## ASV1 and ASV7 closely related taxa how are they interacting? -----


## ASV7 and ASV15 interactions for the MAIN text ----
### asv15 is affected by ASV7? -----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv15') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  dplyr::filter(asv_num_eff == 'asv7') |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv15')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, asv_num_eff, fraction_eff, fraction_b, fraction) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv15_affectedby_asv7_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_hline(yintercept = 0, alpha = 0.5)+
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color, shape = fraction_eff_ed), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  scale_color_identity()+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw() + 
  labs(x = 'Date', 
       y = 'Interaction Strength', 
       color = 'Order', 
       title = 'The ASV7 effect on asv15',
       shape = 'ASV7 fraction\naffecting ASV15') + 
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), title = element_text(size = 6),
        legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(5,15,5,15))

asv15_affectedby_asv7_plot 

asv7_asv15_cluster_interaction_plot <- plot_grid( asv7_asv15_cluster,
           asv15_affectedby_asv7_plot,
           labels = c('A', 'B'),
           label_size = 10,
           label_fontface = 'plain',
                   ncol = 1)

ggsave(
  plot = asv7_asv15_cluster_interaction_plot, 
  filename = 'asv7_asv15_cluster_interaction_plot.pdf',
  path = 'results/figures/main_df3/',
  width = 180, 
  height = 150, 
  units = 'mm'
)

## ASV7 is affected by ASV15-----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv7') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  dplyr::filter(asv_num_eff == 'asv15') |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv7')

data <- mdr_tb_m_rclr |> 
  dplyr::mutate(value = as.numeric(value)) |> 
  dplyr::group_by(date, asv_num_eff, fraction_eff, fraction_b, fraction) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                            str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv7_affectedby_asv15_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_hline(yintercept = 0, alpha = 0.5)+
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color, shape = fraction_eff_ed), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  scale_color_identity()+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw() + 
  labs(x = 'Date', 
       y = 'Interaction Strength', 
       color = 'Order', 
       title = 'The ASV15 effect on asv7',
       shape = 'ASV15 fraction\naffecting ASV7') + 
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), title = element_text(size = 6),
        legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(5,15,5,15))
asv7_affectedby_asv15_plot 

interaction_asv7_and_asv15 <- plot_grid(asv15_affectedby_asv7_plot,
          asv7_affectedby_asv15_plot,
          ncol = 1,
          labels = c('A', 'B'))

interaction_asv7_and_asv15

ggsave(
  plot = interaction_asv7_and_asv15 , 
  filename = 'interaction_asv7_and_asv15.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

ggsave(
  plot = asv7_affectedby_asv15_plot , 
  filename = 'asv7_affectedby_asv15_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 100, 
  units = 'mm'
)

#### DOWN HERE I THINK THAT WE SHOULD REMOVE IT ____ ------------

## Role of total bacterial abundance and blooming events 
mdr_tb_m |>
  colnames()

### relationship between interactions and rCLR ----
asv_tab_all_bloo_z_tax_red <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::select(asv_num, fraction, date, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(
    (asv_num == 'asv1'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv15'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv11' & fraction %in% c('0.2'))       |
      (asv_num == 'asv23' & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv31' & fraction %in% c('3'))         |
      (asv_num == 'asv7'  & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv22' & fraction %in% c('3')) |
      (asv_num == 'asv17' & fraction %in% c('3', '0.2')))

bloom_event_red <- bloom_event |>
  dplyr::filter(
    (asv_num == 'asv1'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv11' & fraction %in% c('0.2'))       |
      (asv_num == 'asv23' & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv31' & fraction %in% c('3'))         |
      (asv_num == 'asv7'  & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv22' & fraction %in% c('3')) |
      (asv_num == 'asv17' & fraction %in% c('3', '0.2')))

asv_tab_all_bloo_z_tax_red_ed <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(bloom_event_red, by = c('date', 'asv_num', 'fraction')) |>
  dplyr::mutate(date = as.Date(date))

mdr_abund <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17')) |> 
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::mutate(date = as.Date(date)) |>
  dplyr::mutate(date = as.Date(date)) |>
  left_join(asv_tab_all_bloo_z_tax_red_ed, by = c('date', 'fraction', 'asv_num'))

mdr_abund_ed <- mdr_abund |>
  left_join(tax_bbmo_10y_new, by = c('asv_num' = 'asv_num')) ## add affecting ASVs to the tb

mdr_abund_ed$bloom_event

mdr_abund |>
  colnames()

mdr_abund |>
  dplyr::select(asv_num, asv_num_eff, order.x, order.y)

mdr_abund |>
  colnames()

mdr_abund |>
  dplyr::filter(abs(as.numeric(value)) > 0.05) |>
  ggplot(aes(abundance_value, as.numeric(value)))+
  geom_smooth(aes(group = fraction), method = 'loess',  color = 'darkgrey')+
  geom_point(aes(shape = fraction, alpha = if_else(bloom_event == 'bloom', 0.9, 0.1)))+
  facet_wrap(vars(asv_num))+
  scale_color_manual(values = palette_order_assigned_all)+
  theme_bw()+
  guides(alpha = 'none')+
  labs(x = 'rCLR', y = 'Interaction Strength', color = 'Order')+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank())

## relative abundances ----
## prepare bloomers abundance data
asv_tab_all_bloo_z_tax_red <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(asv_num, fraction, date, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17')) |>
  left_join(bloom_event)

asv_tab_all_bloo_z_tax_red_ed <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(date = as.Date(date))

mdr_abund <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17')) |> 
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::mutate(date = as.Date(date)) |>
  left_join(asv_tab_all_bloo_z_tax_red_ed, by = c('date', 'fraction', 'asv_num'))

mdr_abund %$% date
asv_tab_all_bloo_z_tax_red_ed %$% date

mdr_abund$abundance_value

mdr_abund |>
  ggplot(aes(abundance_value, as.numeric(value)))+
  geom_point(aes(shape = fraction, alpha = if_else(bloom_event == 'bloom', 0.6, 0.2)))+
  facet_wrap(vars(asv_num))+
  scale_color_manual(values = palette_order_assigned_all)+
  theme_bw()+
  guides(alpha = 'none')+
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order')+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank())

### Plotting the summary of interactions on each ASV and their abundances ----
palette_interactions <- c('interaction_negative' = '#304C57',
                           'interaction_positive' = '#A1D6C1')

labs_interactions <- as_labeller(c('interaction_negative' = 'Negative',
                                 'interaction_positive' = 'Positive'))

labs_fraction <- as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                               '3' = 'Particle attached (3-20 um)',
                               'asv11' = 'Alteromonadaceae asv11',
                              'asv1' =  'Cyanobiaceae asv1', 
                              'asv7' =  'Cyanobiaceae asv7', 
                              'asv31' =  'Cyanobiaceae asv31', 
                              'asv23' = 'Rubritaleaceae asv23',
                              'asv22' = 'Alcanivoracaceae1 asv22',
                               'asv17' = 'Sphingomonadaceae asv17'))

mdr_interactions_rclr |>
  distinct(asv_num, family)

#### rCLR -----
bloom_event <- bloom_event |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 
  
asv_tab_all_bloo_z_tax_red <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::select(asv_num, fraction, date, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17')) |>
  left_join(bloom_event, by = c('date', 'asv_num', 'fraction'))

asv_tab_all_bloo_z_tax_red$bloom_event |>
  unique()

mdr_tb_m_positive <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(as.numeric(interaction_strength) > 0) |>
  dplyr::group_by(date, fraction, asv_num) |>
  dplyr::reframe(interaction_positive = sum(as.numeric(interaction_strength)))

mdr_tb_m_negative <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(as.numeric(interaction_strength) < 0) |>
  dplyr::group_by(date, fraction, asv_num) |>
  dplyr::reframe(interaction_negative = sum(as.numeric(interaction_strength)))

mdr_tb_m_interactions <- mdr_tb_m_positive |>
  left_join(mdr_tb_m_negative) |>
  dplyr::mutate(date = as.Date(date))

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value.x))|>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  dplyr::filter(
    (asv_num == 'asv1'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv11' & fraction %in% c('0.2'))       |
      (asv_num == 'asv23' & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv31' & fraction %in% c('3'))         |
      (asv_num == 'asv7'  & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv22' & fraction %in% c('3')) |
      (asv_num == 'asv17' & fraction %in% c('3', '0.2'))) |>
  dplyr::mutate(date = as.Date(date))

asv_tab_all_bloo_z_tax_red_asv$date
mdr_tb_m_interactions$date

mdr_interactions_rclr <- mdr_tb_m_interactions |>
  left_join(asv_tab_all_bloo_z_tax_red_asv) |>
  pivot_longer(cols = c('interaction_positive', 'interaction_negative'), 
               values_to = 'value', names_to = 'variable') |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'))

plot_interactions_summary_rclr <- mdr_interactions_rclr  |>
  ggplot(aes(date, as.numeric(value)))+
  scale_y_continuous(sec.axis = sec_axis(~./1 , name = 'rCLR'))+
  geom_smooth(aes(date, y = as.numeric(abundance_value.x)*1), span = 0.1, 
              color = '#D6AAA1',  fill = '#D6AAA1', se = T, linetype = 'dotted')+
  geom_smooth(method = 'loess', span = 0.1, aes(color = variable, group = variable, fill = variable))+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 0.9, 0),color = variable, group = variable, fill = variable))+
  facet_wrap(asv_num~fraction, labeller = labs_fraction, dir = 'v', ncol = 3, scales = 'free_y')+#, labeller = labs_fraction, dir = 'v'
  scale_color_manual(values = palette_interactions, labels = labs_interactions)+
  scale_fill_manual(values = palette_interactions, labels = labs_interactions)+
  theme_bw()+
  labs(x = 'Date', y = 'Interaction Strength', color = 'Interactions', fill = 'Interactions')+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 8),
        panel.border = element_blank(),
        legend.position = 'bottom')

plot_interactions_summary_rclr

# ggsave(filename = 'plot_interactions_summary_rclr_v3.pdf',
#        plot = plot_interactions_summary_rclr,
#        path = 'results/figures/',
#        width = 180, height = 300, units = 'mm')

### When there is a bloom event are the interactions stronger? ---- 
mdr_abund |>
  colnames()

mdr_tb_m_positive <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(as.numeric(interaction_strength) > 0) |>
  dplyr::group_by(date, fraction, asv_num) |>
  dplyr::reframe(interaction_positive = sum(as.numeric(interaction_strength)))

mdr_tb_m_negative <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(as.numeric(interaction_strength) < 0) |>
  dplyr::group_by(date, fraction, asv_num) |>
  dplyr::reframe(interaction_negative = sum(as.numeric(interaction_strength)))

mdr_tb_m_interactions <- mdr_tb_m_positive |>
  left_join(mdr_tb_m_negative) |>
  dplyr::mutate(date = as.Date(date))

bloom_event_red <-  bloom_event_red |>
  dplyr::mutate(date = as.character(date)) 

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(
    (asv_num == 'asv1'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv11' & fraction %in% c('0.2'))       |
      (asv_num == 'asv23' & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv31' & fraction %in% c('3'))         |
      (asv_num == 'asv7'  & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv22' & fraction %in% c('3')) |
      (asv_num == 'asv17' & fraction %in% c('3', '0.2'))) |>
  dplyr::mutate(date = as.Date(date)) |>
  dplyr::mutate(date = as.character(date)) |>
  left_join(bloom_event_red, by = c('asv_num', 'fraction', 'date'))

unique(asv_tab_all_bloo_z_tax_red_asv$date)
unique(bloom_event_red$date)
unique(mdr_tb_m_interactions$date)

asv_tab_all_bloo_z_tax_red_asv$date
mdr_tb_m_interactions$date

mdr_tb_m_interactions |>
  str()

asv_tab_all_bloo_z_tax_red_asv |>
  str()

asv_tab_all_bloo_z_tax_red_asv$bloom_event |>
  is.na() |>
  summary() 

test <-  asv_tab_all_bloo_z_tax_red_asv |>
dplyr::filter(date %in% c('2004-03-22', '2005-02-15', '2005-05-10')) 

test2 <- mdr_tb_m_interactions |>
  dplyr::mutate(date = as.Date(date)) |>
  dplyr::mutate(date = as.character(date)) |>
  dplyr::filter(date %in% c('2004-03-22', '2005-02-15', '2005-05-10')) 

test |>
  left_join(test2)

test2 |>
  left_join(test)

mdr_interactions_rclr <- mdr_tb_m_interactions |>
  dplyr::mutate(date = as.Date(date)) |>
  dplyr::mutate(date = as.character(date)) |>
  left_join(asv_tab_all_bloo_z_tax_red_asv, by = c('asv_num', 'fraction', 'date' = 'date')) |>
  left_join(bloo_taxonomy, by = 'asv_num') |>
  dplyr::mutate(
    interaction_negative = case_when(is.na(interaction_negative) ~ 0,  # Numeric '0'
                                     TRUE ~ as.numeric(interaction_negative)),
    interaction_positive = case_when(is.na(interaction_positive) ~ 0,  # Numeric '0'
                                     TRUE ~ as.numeric(interaction_positive))
  ) |>
  dplyr::mutate(interaction_balance = interaction_positive+interaction_negative)

mdr_interactions_rclr |>
  dplyr::filter(is.na(bloom_event))

mdr_interactions_rclr |>
  colnames()

mdr_interactions_rclr |>
  dplyr::filter(!is.na(bloom_event)) |>
  ggplot(aes(bloom_event, as.numeric(interaction_balance)))+
  geom_point(aes(shape = fraction), position = position_dodge(width = 0.25))+ #alpha = if_else(bloom_event == 'bloom', 0.9, 0.1)
  geom_boxplot(aes(), notch = T, alpha = 0.2)+
  facet_wrap(vars(asv_num), ncol = 2)+
  scale_x_discrete(labels = c('bloom' = 'Bloom', 'no-bloom' = 'No Bloom'))+
  scale_color_manual(values = palette_order_assigned_all)+
  theme_bw()+
  guides(alpha = 'none')+
  labs(x = '', y = 'Interaction Balance', color = '')+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank())

#### Relative abundance ----
asv_tab_all_bloo_z_tax_red <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(asv_num, fraction, date, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17'))

mdr_tb_m_positive <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(as.numeric(interaction_strength) > 0) |>
  dplyr::group_by(date, fraction, asv_num) |>
  dplyr::reframe(interaction_positive = sum(as.numeric(interaction_strength)))

mdr_tb_m_negative <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(as.numeric(interaction_strength) < 0) |>
  dplyr::group_by(date, fraction, asv_num) |>
  dplyr::reframe(interaction_negative = sum(as.numeric(interaction_strength)))

mdr_tb_m_interactions <- mdr_tb_m_positive |>
  left_join(mdr_tb_m_negative) |>
  dplyr::mutate(date = as.Date(date))

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  dplyr::filter(
    (asv_num == 'asv1'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv11' & fraction %in% c('0.2'))       |
      (asv_num == 'asv23' & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv31' & fraction %in% c('3'))         |
      (asv_num == 'asv7'  & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv22' & fraction %in% c('3')) |
      (asv_num == 'asv17' & fraction %in% c('3', '0.2'))) |>
  dplyr::mutate(date = as.Date(date))

# asv_tab_all_bloo_z_tax_red_asv$date
# mdr_tb_m_interactions$date

mdr_interactions_rclr <- mdr_tb_m_interactions |>
  left_join(asv_tab_all_bloo_z_tax_red_asv) |>
  #dplyr::mutate(abundance_value = abundance_value) |>
  pivot_longer(cols = c('interaction_positive', 'interaction_negative'), 
               values_to = 'value', names_to = 'variable') |>
  left_join(bloo_taxonomy, by = 'asv_num')

plot_interactions_summary_relative_abundance <- mdr_interactions_rclr  |>
  ggplot(aes(date, as.numeric(value)))+
  scale_y_continuous(sec.axis = sec_axis(~./10 , name = 'Relative abundance'))+
  # geom_line(aes(date, y = as.numeric(abundance_value)*3), 
  #             color = '#D6AAA1', linetype = 'dotted')+
  # geom_smooth(aes(date, y = as.numeric(abundance_value)*10), span = 0.1,
  #             color = '#D6AAA1',  fill = '#D6AAA1', se = T, linetype = 'dotted')+
  geom_smooth(method = 'loess', span = 0.1, aes(color = variable, group = variable, fill = variable), alpha = 0.5)+
  geom_smooth(aes(date, y = as.numeric(abundance_value)*10), span = 0.1,
              color = '#D6AAA1',  fill = '#D6AAA1', se = T, linetype = 'dotted', alpha = 0.5)+
  #geom_line( aes(color = variable, group = variable))+
  facet_wrap(asv_num~fraction, labeller = labs_fraction, dir = 'v', ncol = 3, scales = 'free_y')+#, labeller = labs_fraction, dir = 'v'
  scale_color_manual(values = palette_interactions, labels = labs_interactions)+
  scale_fill_manual(values = palette_interactions, labels = labs_interactions)+
  #scale_linetype_manual()+
  theme_bw()+
  labs(x = 'Date', y = 'Interaction Strength', color = 'Interactions', fill = 'Interactions')+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 8),
        panel.border = element_blank(),
        legend.position = 'bottom')

plot_interactions_summary_relative_abundance
# 
# ggsave(filename = 'plot_interactions_summary_relative_abundance_v2.pdf',
#        plot = plot_interactions_summary_relative_abundance,
#        path = 'results/figures/',
#        width = 180, height = 280, units = 'mm')

## Relationship between interactions and abundances ----
### rCLR ----
asv_tab_all_bloo_z_tax_red <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::select(asv_num, fraction, date, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17'))

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num %in% c('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |> 
  dplyr::filter(
    (asv_num == 'asv1'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv11' & fraction %in% c('0.2'))       |
      (asv_num == 'asv23' & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv31' & fraction %in% c('3'))         |
      (asv_num == 'asv7'  & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv22' & fraction %in% c('3')) |
      (asv_num == 'asv17'  & fraction %in% c('0.2', '3')) ) |>
  dplyr::mutate(date = as.Date(date))

mdr_interactions_rclr <- mdr_tb_m_interactions |>
  left_join(asv_tab_all_bloo_z_tax_red_asv) |>
  #dplyr::mutate(abundance_value = abundance_value) |>
  # pivot_longer(cols = c('interaction_positive', 'interaction_negative'), 
  #              values_to = 'value', names_to = 'variable') |>
  left_join(bloo_taxonomy, by = 'asv_num') |>
  dplyr::mutate(
    interaction_negative = case_when(is.na(interaction_negative) ~ 0,  # Numeric '0'
                                     TRUE ~ as.numeric(interaction_negative)),
    interaction_positive = case_when(is.na(interaction_positive) ~ 0,  # Numeric '0'
                                     TRUE ~ as.numeric(interaction_positive))
  )

##check normality before correlation
shapiro.test(as.numeric(mdr_interactions_rclr$interaction_positive)) # => p-value < 2.2e-16 ( NO NORMALITY)
ggqqplot(as.numeric(mdr_interactions_rclr$interaction_positive))

shapiro.test(as.numeric(mdr_interactions_rclr$interaction_negative)) # => p-value < 2.2e-16 ( NO NORMALITY)
ggqqplot(as.numeric(mdr_interactions_rclr$interaction_negative))

## for non normal data
cor_spearman <- rcorr(as.matrix(mdr_interactions_rclr[,c(4,5)]), type = 'spearman')

plot_interactions_summary_rclr_relationship <- mdr_interactions_rclr  |>
  dplyr::mutate(abundance_value = case_when(
    is.na(as.numeric(abundance_value)) ~ 0,   # Use numeric 0 instead of '0'
    TRUE ~ as.numeric(abundance_value)        # Ensure the value is numeric
  )) |>
  ggplot(aes(interaction_positive, abs(interaction_negative)))+
  geom_abline(slope = 1, intercept = 0, color = '#0F3331', linetype = 'dashed')+
  #scale_y_continuous(sec.axis = sec_axis(~./10 , name = 'rCLR'))+
  geom_point(aes(size = as.numeric(abundance_value), shape = fraction, color = fraction, fill = fraction, alpha = 0.7))+
  stat_poly_line(aes(group = fraction, color = fraction, fill = fraction))+
  stat_cor(aes(group = fraction, color = fraction,  label = paste0(..p.label..)),label.x = 0.2, label.y = 7, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  stat_cor(aes(group = fraction, color = fraction,  label = paste0(..r.label..)),label.x = 0.2, label.y = 8, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #facet_wrap(vars(asv_num), labeller = labs_fraction, dir = 'v')+#, labeller = labs_fraction, dir = 'v'
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
  scale_shape_discrete(labels = labs_fraction)+
  scale_size_continuous(range = c(0, 5), limits = c(0, 10), breaks = c(0, 2.5, 5, 10) )+
  #scale_linetype_manual()+
  guides(alpha = 'none')+
  theme_bw()+
  labs(x = 'Interaction positive', y = 'Interaction negative', color = 'Interactions', fill = 'Interactions',
       shape = 'Fraction', size = 'rCLR')+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        aspect.ratio = 10/10)

plot_interactions_summary_rclr_relationship

# ggsave(filename = 'plot_interactions_summary_rclr_relationship_v2.pdf',
#        plot = plot_interactions_summary_rclr_relationship,
#        path = 'results/figures/',
#        width = 180, height = 180, units = 'mm')

### Relative abundance ----
asv_tab_all_bloo_z_tax_red <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(asv_num, fraction, date, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17'))

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  dplyr::filter(
    (asv_num == 'asv1'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv11' & fraction %in% c('0.2'))       |
      (asv_num == 'asv23' & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv31' & fraction %in% c('3'))         |
      (asv_num == 'asv7'  & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv22' & fraction %in% c('3'))|
      (asv_num == 'asv17'  & fraction %in% c('0.2', '3')) ) |>
  dplyr::mutate(date = as.Date(date))

mdr_interactions_rclr <- mdr_tb_m_interactions |>
  left_join(asv_tab_all_bloo_z_tax_red_asv) |>
  #dplyr::mutate(abundance_value = abundance_value) |>
  # pivot_longer(cols = c('interaction_positive', 'interaction_negative'), 
  #              values_to = 'value', names_to = 'variable') |>
  left_join(bloo_taxonomy, by = 'asv_num') |>
  dplyr::mutate(
    interaction_negative = case_when(is.na(interaction_negative) ~ 0,  # Numeric '0'
                                     TRUE ~ as.numeric(interaction_negative)),
    interaction_positive = case_when(is.na(interaction_positive) ~ 0,  # Numeric '0'
                                     TRUE ~ as.numeric(interaction_positive))
  )

##check normality before correlation
shapiro.test(as.numeric(mdr_interactions_rclr$interaction_positive)) # => p-value < 2.2e-16 ( NO NORMALITY)
ggqqplot(as.numeric(mdr_interactions_rclr$interaction_positive))

shapiro.test(as.numeric(mdr_interactions_rclr$interaction_negative)) # => p-value < 2.2e-16 ( NO NORMALITY)
ggqqplot(as.numeric(mdr_interactions_rclr$interaction_negative))

## for non normal data
cor_spearman <- rcorr(as.matrix(mdr_interactions_rclr[,c(4,5)]), type = 'spearman')

plot_interactions_summary_relative_abundance_relationship <- mdr_interactions_rclr  |>
  dplyr::mutate(abundance_value = case_when(
    is.na(as.numeric(abundance_value)) ~ 0,   # Use numeric 0 instead of '0'
    TRUE ~ as.numeric(abundance_value)        # Ensure the value is numeric
  )) |>
  ggplot(aes(interaction_positive, abs(interaction_negative)))+
  geom_abline(slope = 1, intercept = 0, color = '#0F3331', linetype = 'dashed')+
  #scale_y_continuous(sec.axis = sec_axis(~./10 , name = 'Relative abundance'))+
  geom_point(aes(size = abundance_value, shape = fraction, color = fraction, fill = fraction, alpha = 0.7))+
  stat_poly_line(aes(group = fraction, color = fraction, fill = fraction))+
  stat_cor(aes(group = fraction, color = fraction,  label = paste0(..p.label..)),label.x = 0.2, label.y = 7, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  stat_cor(aes(group = fraction, color = fraction,  label = paste0(..r.label..)),label.x = 0.2, label.y = 8, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #facet_wrap(vars(asv_num), labeller = labs_fraction, dir = 'v')+#, labeller = labs_fraction, dir = 'v'
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
  scale_shape_discrete(labels = labs_fraction)+
  scale_size_continuous(range = c(0, 5), limits = c(0, 0.5), breaks = c(0, 0.25, 0.5) )+
  #scale_linetype_manual()+
  guides(alpha = 'none')+
  theme_bw()+
  labs(x = 'Interaction positive', y = 'Interaction negative', color = 'Interactions', fill = 'Interactions',
       shape = 'Fraction', size = 'Relative Abundance')+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        aspect.ratio = 10/10)

plot_interactions_summary_relative_abundance_relationship

# ggsave(filename = 'plot_interactions_summary_relative_abundance_relationship_v2.pdf',
#        plot = plot_interactions_summary_relative_abundance_relationship,
#        path = 'results/figures/',
#        width = 180, height = 180, units = 'mm')

## Bloomers effects ON other ASVs -----
### I would like to observe what is the effect that some bloomers have on other community members -----
asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::select(asv_num, fraction, date, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  #dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num %in% c ('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  dplyr::filter(
    (asv_num == 'asv1'  & fraction %in% c('0.2', '3')) |
      (asv_num == 'asv11' & fraction %in% c('0.2'))       |
      (asv_num == 'asv23' & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv31' & fraction %in% c('3'))         |
      (asv_num == 'asv7'  & fraction %in% c('0.2', '3'))  |
      (asv_num == 'asv22' & fraction %in% c('3')) |
      (asv_num == 'asv17' & fraction %in% c('3', '0.2')),
    ) |>
  dplyr::mutate(date = as.Date(date))

mdr_tb_m_positive <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(time, variable, bp_asv1, bn_asv1, bp_asv11, bp_asv23, bn_asv23,
                bn_asv31, bn_asv7, bp_asv7, bn_asv7) |>
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  #left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  left_join(variable_num_tb, by = c('variable' = 'variable_num')) |>
  separate(variable_name, sep = '_', into = c('fraction_aff', 'asv_num_aff')) |>
  dplyr::filter(asv_num_aff != asv_num_eff) |>
  dplyr::filter(as.numeric(interaction_strength) > 0) |>
  dplyr::group_by(asv_num_aff, time, fraction_aff) |>
  dplyr::reframe(interaction_strength = sum(as.numeric(interaction_strength)))

mdr_tb_m_negative <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(time, variable, bp_asv1, bn_asv1, bp_asv11, bp_asv23, bn_asv23,
                bn_asv31, bn_asv7, bp_asv7, bn_asv7) |>
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  #left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  left_join(variable_num_tb, by = c('variable' = 'variable_num')) |>
  separate(variable_name, sep = '_', into = c('fraction_aff', 'asv_num_aff')) |>
  dplyr::filter(asv_num_aff != asv_num_eff) |>
  dplyr::filter(as.numeric(interaction_strength) < 0) |>
  dplyr::group_by(asv_num_aff, time, fraction_aff) |>
  dplyr::reframe(interaction_strength = sum(as.numeric(interaction_strength)))

mdr_tb_m_interactions <- mdr_tb_m_positive |>
  full_join(mdr_tb_m_negative) |>
  dplyr::mutate(fraction = case_when(str_detect(fraction_aff, 'bp') ~ '0.2',
                                     str_detect(fraction_aff, 'bn') ~ '3')) |>
  left_join(time_num, by = c('time' = 'time_num', 'fraction')) |>
  dplyr::mutate(date = as.Date(date)) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_aff' = 'asv_num'))

asv_tab_all_bloo_z_tax_red_asv$date
mdr_tb_m_interactions$date

# mdr_interactions_rclr <- mdr_tb_m_interactions |>
#   #left_join(asv_tab_all_bloo_z_tax_red_asv) |>
#   pivot_longer(cols = c('interaction_positive', 'interaction_negative'), 
#                values_to = 'value', names_to = 'variable') |>
#   left_join(bloo_taxonomy, by = 'asv_num')
mdr_tb_m_interactions$order |>
  unique()

mdr_tb_m_interactions$order <- factor(mdr_tb_m_interactions$order, levels = c('Synechococcales',
                                                                              'Enterobacterales',
                                                                              'Pseudomonadales',
                                                                              'SAR11 clade',
                                                                              'Puniceispirillales',
                                                                              'Rhodobacterales',
                                                                              'Sphingomonadales',
                                                                              'Rhodospirillales',
                                                                              'Flavobacteriales',
                                                                              'Verrucomicrobiales'))

mdr_tb_m_interactions |>
  colnames()

plot_interactions_summary_rclr <- mdr_tb_m_interactions  |>
  left_join(asv_tab_all_bloo_z_tax_red_asv, by = c('asv_num_aff' = 'asv_num')) |>
  ggplot(aes(date, as.numeric(interaction_strength)))+
  coord_flip()+
  geom_point(aes(shape = fraction, color = order, fill = order))+
  #scale_y_continuous(sec.axis = sec_axis(~./1 , name = 'rCLR'))+
  # geom_smooth(aes(date, y = as.numeric(abundance_value)*1), span = 0.1, 
  #             color = '#D6AAA1',  fill = '#D6AAA1', se = T, linetype = 'dotted')+
  #geom_smooth(method = 'loess', span = 0.1, aes(color = variable, group = variable, fill = variable))+
  geom_line( aes(x = as.numeric(abundance_value), color = order, group = asv_num))+
  facet_wrap(vars(order))+#, labeller = labs_fraction, dir = 'v'
  geom_violin(alpha = 0.1)+
  scale_color_manual(values = palette_order_assigned_all)+
  scale_fill_manual(values = palette_order_assigned_all)+
  theme_bw()+
  labs(x = 'Date', y = 'Interaction Strength', color = 'Interactions', fill = 'Interactions')+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 10),
        panel.border = element_blank(),
        legend.position = 'bottom')

plot_interactions_summary_rclr

# ggsave(filename = 'plot_interactions_summary_rclr.pdf',
#        plot = plot_interactions_summary_rclr,
#        path = 'results/figures/',
#        width = 180, height = 280, units = 'mm')

## -----  MDR visualization of the ASVs that affect my bloomers (in a boxplot) ----- ##
mdr_abund |>
  colnames()

mdr_abund %$%
  asv_num |>
  unique()

max_interactions <- mdr_abund |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::group_by(asv_num) |>
  slice_max(order_by = as.numeric(value), n = 10) |>
  distinct(asv_num_eff)

mdr_abund |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num_eff %in% max_interactions$asv_num_eff) |>
  ggplot(aes(fct_infreq(asv_num_eff), as.numeric(value)))+
  facet_wrap(vars(asv_num), scales = 'free_x', ncol = 1)+
  geom_point(aes(shape = fraction, color = order.y, alpha = ifelse(bloom_event == 'bloom', 1, '0.2')), 
             position = position_jitter(width =  0.25))+
  scale_color_manual(values = palette_order_assigned_all)+
  geom_boxplot( alpha = 0.2)+ #notch = T
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+  
  theme(strip.background = element_rect(fill = 'transparent'),
                     text = element_text(size = 12),
                     axis.text.x = element_text(size = 8),
                     strip.text = element_text(size = 10),
                     panel.border = element_blank(),
                     legend.position = 'bottom')

## only bloomers interaction between them ----
mdr_abund |>
  colnames()

mdr_abund$bloom_event |>
  unique()

mdr_abund |>
  dplyr::filter(bloom_event == 'bloom')
  
mdr_abund |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num_eff %in% c('asv23', 'asv1', 'asv11', 'asv22', 'asv31', 'asv7', 'asv17')) |>
  ggplot(aes(fct_infreq(asv_num_eff), as.numeric(value)))+
  facet_grid(vars(asv_num), scales = 'free')+
  geom_point(aes(shape = fraction, color = order.y, alpha = ifelse(bloom_event == 'bloom', 1, '0.2')), position = position_jitter(width =  0.25))+
  scale_color_manual(values = palette_order_assigned_all)+
  #geom_boxplot( alpha = 0.2)+ #notch = T
  geom_violin(alpha = 0.2)+
  guides(alpha = 'none')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+  
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 10),
        panel.border = element_blank(),
        legend.position = 'bottom')

## MDR jacobians -1, 0, 1 ----
### additionally normalize by richness 
asv_tab_bbmo_10y_w_rar <- read.csv2('data/asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()|> 
  rename('sample_id' = 'X')  ## richness should always be calculated with the rarefied dataset

asv_tab_bbmo_10y_w_rar |>
  colnames()

## calculate shannon index normalize interactions by diversity ----
asv_tab_bbmo_10y_w_rar_t <- asv_tab_bbmo_10y_w_rar |>
  t() |>
  as_tibble(rownames = NULL)

asv_tab_bbmo_10y_w_rar_t |>
  colnames() <- asv_tab_bbmo_10y_w_rar_t[1,]

asv_tab_bbmo_10y_w_rar_t <- asv_tab_bbmo_10y_w_rar_t[-1,]

shannon_rar <- asv_tab_bbmo_10y_w_rar_t |>
  dplyr::mutate(across(everything(), as.numeric)) |>
  microbiome::alpha('shannon') |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble() |>
  left_join(m_bbmo_10y) 

observed_rar <- asv_tab_bbmo_10y_w_rar_t |>
  dplyr::mutate(across(everything(), as.numeric)) |>
  microbiome::alpha(index = 'observed') |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble() |>
  left_join(m_bbmo_10y) 

shannon_rar |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, diversity_shannon))+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_linetype_manual( labels = labs_fraction, values = c('0.2' = 1, '3' = 2))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(data = shannon_rar |>
              dplyr::filter(fraction == '3'), aes(date, diversity_shannon), alpha = 0.3)+
  labs(x = 'Date', y = 'Shannon Index', color = '', linetype = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside',
    axis.ticks.length.y =  unit(0.2, "mm"))

observed_rar  |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, observed))+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_linetype_manual( labels = labs_fraction, values = c('0.2' = 1, '3' = 2))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(data = observed_rar |>
              dplyr::filter(fraction == '3'), aes(date, observed), alpha = 0.3)+
  labs(x = 'Date', y = 'Observed Richness', color = '', linetype = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside',
    axis.ticks.length.y =  unit(0.2, "mm"))

mdr_tb_m |>
  colnames()

mdr_tb_m_bin <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17')) |> 
  pivot_longer(cols = starts_with(c('bp', 'bn')), values_to = 'interaction_strength') |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::mutate(interaction_strength_ed = case_when(as.numeric(interaction_strength) > 0 ~ '1',
                                                    as.numeric(interaction_strength) < 0 ~ '-1',
                                                    as.numeric(interaction_strength) == 0 ~ '0')) |>
  left_join(m_bbmo_10y)

mdr_tb_m_bin |>
  colnames()

# abundance data 
abund_data <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |> 
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17')) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::select(sample_id, abundance_value, asv_num, fraction) 

data <- mdr_tb_m_bin |>
  group_by(fraction, asv_num, date, sample_id, recurrency) |>
  dplyr::reframe(interactions = sum(as.numeric(interaction_strength_ed))) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))|>
  left_join(abund_data, by = c('sample_id', 'fraction', 'asv_num'))

data |>
  dplyr::filter(fraction == '3') |>
  ggplot(aes(date, interactions))+
  geom_line(aes(group = fraction))+
  facet_wrap(vars(asv_num), ncol = 2)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  geom_line(data = data |>
              dplyr::filter(fraction == '3'), aes(date, abundance_value), color = 'darkred')+
  # geom_line(data =  data |>
  #             dplyr::filter(fraction == '0.2'), aes(date, abundance_value), color = 'purple')+
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_linetype_manual( labels = labs_fraction, values = c('0.2' = 1, '3' = 2))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(data = data |>
              dplyr::filter(fraction == '3'), aes(date, interactions), alpha = 0.3)+
  labs(x = 'Date', y = 'Observed Richness', color = '', linetype = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside',
    axis.ticks.length.y =  unit(0.2, "mm"))

data |>
  dplyr::filter(fraction == '0.2') |>
  ggplot(aes(date, interactions))+
  geom_line(aes(group = fraction))+
  facet_wrap(vars(asv_num), ncol = 2)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  # geom_line(data = data |>
  #             dplyr::filter(fraction == '3'), aes(date, abundance_value), color = 'darkred')+
  geom_line(data =  data |>
              dplyr::filter(fraction == '0.2'), aes(date, abundance_value), color = 'grey')+
  geom_line(data =  data |>
              dplyr::filter(fraction == '0.2'), aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_linetype_manual( labels = labs_fraction, values = c('0.2' = 1, '3' = 2))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Date', y = 'Interaction', color = '', linetype = '', shape = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside',
    axis.ticks.length.y =  unit(0.2, "mm"))

observed_rar_red <- observed_rar |>
  dplyr::select(sample_id, observed)

shannon_rar_red <- shannon_rar |>
  dplyr::select(sample_id, diversity_shannon)

interactions_binnary_mdr_plot <- data |>
  left_join(shannon_rar_red) |>
  #dplyr::filter(fraction == '0.2') |>
  ggplot(aes(abundance_value, interactions/diversity_shannon))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_point(aes(shape = fraction, color = fraction))+
  facet_wrap(vars(interaction(asv_num, recurrency)), ncol = 3)+
  scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+

  geom_smooth(method = 'loess', aes(group = fraction, color = fraction))+
  labs(x = 'Abundance robust - Centered Log Ratio', y = 'Interactions summary', color = '', linetype = '', shape = '')+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside',
    axis.ticks.length.y =  unit(0.2, "mm"))

interactions_binnary_mdr_plot

# ggsave(plot = interactions_binnary_mdr_plot, filename = 'interactions_binnary_mdr_plot_norm_shannon.pdf',
#        path = 'results/figures/',
#        width = 180, height = 180, units = 'mm')

## relative abundance ----
# abundance data 
abund_data <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |> 
  dplyr::filter(asv_num %in% c('asv11', 'asv7', 'asv1', 'asv22', 'asv23', 'asv31', 'asv17')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(sample_id, abundance_value, asv_num, fraction) 

data <- mdr_tb_m_bin |>
  group_by(fraction, asv_num, date, sample_id, recurrency) |>
  dplyr::reframe(interactions = sum(as.numeric(interaction_strength_ed))) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))|>
  left_join(abund_data, by = c('sample_id', 'fraction', 'asv_num'))

interactions_binnary_mdr_plot <- data |>
  left_join(shannon_rar_red) |>
  #dplyr::filter(fraction == '0.2') |>
  ggplot(aes(abundance_value, interactions/diversity_shannon))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_point(aes(shape = fraction, color = fraction))+
  facet_wrap(vars(interaction(asv_num, recurrency)), ncol = 3)+
  scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  
  geom_smooth(method = 'loess', aes(group = fraction, color = fraction))+
  labs(x = 'Relative Abundance', y = 'Interactions summary', color = '', linetype = '', shape = '')+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside',
    axis.ticks.length.y =  unit(0.2, "mm"))

interactions_binnary_mdr_plot

# ggsave(plot = interactions_binnary_mdr_plot, filename = 'interactions_binnary_mdr_plot__relabund_norm_shannon.pdf',
#        path = 'results/figures/',
#        width = 180, height = 180, units = 'mm')

interactions_binnary_mdr_plot <- data |>
  left_join(shannon_rar_red) |>
  #dplyr::filter(fraction == '0.2') |>
  ggplot(aes(abundance_value, interactions))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_point(aes(shape = fraction, color = fraction))+
  facet_wrap(vars(interaction(asv_num, recurrency)), ncol = 3)+
  scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  
  geom_smooth(method = 'loess', aes(group = fraction, color = fraction))+
  labs(x = 'Relative Abundance', y = 'Interactions summary', color = '', linetype = '', shape = '')+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside',
    axis.ticks.length.y =  unit(0.2, "mm"))

interactions_binnary_mdr_plot

# ggsave(plot = interactions_binnary_mdr_plot, filename = 'interactions_binnary_mdr_plot__relabund_norm_shannon.pdf',
#        path = 'results/figures/',
#        width = 180, height = 180, units = 'mm')
