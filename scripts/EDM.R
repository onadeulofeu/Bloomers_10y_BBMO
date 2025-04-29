# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     EDMs CCM and MDR-S MAP                  ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Results created by Carmen-García Comas          ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# EDM empirical dynamic modeling
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
bbmo_10y <-readRDS("~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/blphy10years.rds") ##8,052 asv amb totes les mostres, reads / ASV and sample
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

colnames(tax_bbmo_10y_new) <- c("asv_num",'seq', "kingdom", "phylum", "class", "order", "family", "genus")

colnames(m_bbmo_10y) <- c("sample_id", "project", "location", "code",             
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

## we could also do it with rCLR add here in case we want to do it ----
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

asv_tab_10y_l_rel <- asv_tab_10y_l_rel |>
  dplyr::mutate(sample_id = as.character(sample_id))

m_bbmo_10y_sim <- m_bbmo_10y |>
  dplyr::select(fraction, sample_id,  Year = year, Month = month, Day = day) |>
  dplyr::mutate(sample_id = as.character(sample_id))

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

# Create a table with presence absence of each ASV to not infiere relationships when an ASV is not there -----
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
bloo_taxonomy <- asv_tab_all_bloo_z_tax %>%
  dplyr::select(phylum_f, class_f, order_f, family_f, genus_f, asv_num_f) |>
  unique()

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
           #color = palette_gradient(5),
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
        panel.grid = element_blank(), text = element_text(size = 8),
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
        panel.grid = element_blank(), text = element_text(size = 8),
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
                  TRUE ~ asv_num)) |>
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
  dplyr::filter(asv_num %in% bloo_02$value & fraction == '0.2' |
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
  dplyr::filter(asv_num_rows %in% bloo_02$value & fraction_rows == '0.2' |
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
#         text = element_text(size = 8),
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
                               'env' = 'Environmental\nvariables',
                               'bio' = 'Biological Variables',
                               'physicochem' = 'Physicochemical Variables'))

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
  facet_wrap(fraction_predicted~bloom, labeller = labs_fraction, dir = 'v')+
  scale_x_discrete(labels = labs_fraction_env)+
  theme_bw()+
  labs(y = 'Prediction (Rho) ', x= 'Causal effects'#, color = 'Order predicted ASVs'
       )+
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid.major.y  = element_blank(),
        panel.grid.minor  = element_blank(), text = element_text(size = 8),
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
        panel.grid.minor  = element_blank(), text = element_text(size = 8),
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

### stats for environmental variables separated in physicochemical and biological (MAIN) --------
data_ed2 |>
  colnames()

select_bloomers <- data_ed2 |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::distinct(asv_num_predicted, frequency_predicted, fraction_predicted) |>
  dplyr::filter(frequency_predicted != 'No Bloomer')

## add statistics (by group)
### anova to compare mean of the different groups 
data_ed3 <- data_ed2 |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::mutate(fraction_bloom = paste0(fraction_predicted,bloom))

data_ed3 |>
  str()

data_ed3 |>
  dplyr::filter(bloom == 'no_bloomer') %$%
  unique(asv_num_predicted)

library(car)  ## homocedasticity 
library(dunn.test) ## significant groups in Kruskal Test
library(ggpubr) ## ggqqplot

# I need to perform an anova for each group in case that each group has a normal distribution and that they have homocedasticity ---- 
## I group: 0.2 no bloomers 
data_ed3_f <- data_ed3 |>
  dplyr::filter(fraction_bloom == '0.2no_bloomer')

## check normality 
shapiro.test(as.numeric(data_ed3 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(data_ed3 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed3_f) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    3  74.161 < 2.2e-16 ***
#       2464                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test

# Perform one-way ANOVA to compare the means of Rho across the three groups
## anova_result <- aov(rho ~ fraction_causal, data = data_ed3_f)  # Replace `bloom` with your actual group variable
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed3_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed3_f$rho,                            # The values you are comparing
  g = as.factor(data_ed3_f$fraction_causal),    # The grouping variable
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

### group II ------
## II group: 0.2 bloomers 
data_ed3_f <- data_ed3 |>
  dplyr::filter(fraction_bloom == '0.2bloomer')

## check normality 
shapiro.test(as.numeric(data_ed3 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(data_ed3 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed3_f) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   3  4.8133 0.002743 **
#       296                        
# # ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test

# Perform one-way ANOVA to compare the means of Rho across the three groups
## anova_result <- aov(rho ~ fraction_causal, data = data_ed3_f)  # Replace `bloom` with your actual group variable
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed3_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed3_f$rho,                            # The values you are comparing
  g = as.factor(data_ed3_f$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)
# 
# Comparison of x by group                            
# (Bonferroni)                                  
# Col Mean-|
#   Row Mean |        0.2          3        bio
# ---------+--------------------------------- #
#   3 |  -1.039124
# |     0.8962
# |
#   bio |   2.550914   3.012838
# |     0.0322    0.0078*
#   |
#   phch |   0.453806   1.053441  -1.391741
# |     1.0000     0.8764     0.4920

## group III ----
## I group: 0.2 bloomers 
data_ed3_f <- data_ed3 |>
  dplyr::filter(fraction_bloom == '3no_bloomer')

## check normality 
shapiro.test(as.numeric(data_ed3 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(data_ed3 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed3_f) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
# group   3  21.837 1.236e-13 ***
#       859                      
# ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test

# Perform one-way ANOVA to compare the means of Rho across the three groups
## anova_result <- aov(rho ~ fraction_causal, data = data_ed3_f)  # Replace `bloom` with your actual group variable
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed3_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed3_f$rho,                            # The values you are comparing
  g = as.factor(data_ed3_f$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)
# 
# Comparison of x by group                            
# (Bonferroni)                                  
# Col Mean-|
# Col Mean-|
#   Row Mean |        0.2          3        bio
# ---------+--------------------------------- #
#   3 |  -16.43611
# |    0.0000*
#   |
#   bio |   3.514799   13.49203
# |    0.0013*    0.0000*
#   |
#   phch |  -0.043184   9.947252  -2.677929
# |     1.0000    0.0000*    0.0222*

## Group IV -----
## IV group: 3 bloomers 
data_ed3_f <- data_ed3 |>
  dplyr::filter(fraction_bloom == '3bloomer')

## check normality 
shapiro.test(as.numeric(data_ed3 |>
                          dplyr::filter(fraction_causal == '0.2') %$%
                          rho)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(data_ed3 |>
                      dplyr::filter(fraction_causal == '0.2') %$%
                      rho))

## check homocedasticity 
leveneTest(rho ~ fraction_causal, data = data_ed3_f) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
# group   3  9.0201 8.058e-06 ***
#       478     
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ## Since we don't have normality nor homocedasticity we change to a non parametric test kruskal.test

# Perform one-way ANOVA to compare the means of Rho across the three groups
## anova_result <- aov(rho ~ fraction_causal, data = data_ed3_f)  # Replace `bloom` with your actual group variable
kruskal_result <- kruskal.test(rho ~ fraction_causal, data = data_ed3_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_ed3_f$rho,                            # The values you are comparing
  g = as.factor(data_ed3_f$fraction_causal),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)
# 
# Comparison of x by group                            
# (Bonferroni)                                  
# Comparison of x by group                            
# (Bonferroni)                                  
# Col Mean-|
#   Row Mean |        0.2          3        bio
# ---------+---------------------------------#
#   3 |  -12.47094
# |    0.0000*
#   |
#   bio |   3.737682   11.62683
# |    0.0006*    0.0000*
#   |
#   phch |  -0.336589   6.414832  -2.794449
# |     1.0000    0.0000*    0.0156*



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
  facet_wrap(fraction_predicted~bloom, labeller = labs_fraction, dir = 'v', ncol = 1)+
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
labs_bloom <- as_labeller(c('no_bloomer' = 'No Bloomer',
                            'bloomer' = 'Bloomer'))

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
  dplyr::mutate(bloom = case_when(asv_num_causal %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                  TRUE ~ 'no_bloomer')) |>
  ggplot(aes(fct_reorder(order_family_asv_num, relative_rho), relative_rho))+
  # geom_point(data = data_causal_relationships, aes(order_family_asv_num, rho
  #            #,color = fraction_predicted
  #            ,
  #            alpha = if_else((asv_num_causal %in% bloo_taxonomy$asv_num_f), 0.6, 0.1)), position = position_jitter(width = 0.25),  size = 0.6)+
  geom_boxplot(data = data_causal_relationships, aes(order_family_asv_num, rho,
                                                     alpha = if_else((asv_num_causal %in% bloo_taxonomy$asv_num_f), 0.8, 0.3)),
               notch = T, outliers = F)+ # I already plot the with the geom point
  geom_point(shape = 18, alpha = 0.1)+
  facet_wrap(vars(bloom), labeller = labs_bloom)+
  coord_flip()+
  labs(y = 'Rho', x= 'ASVs causally affectting other ASVs', color = 'Fraction ASVs being causally affected')+
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
#         filename = 'plot_ccm_asvs_causally_affecting_other_asvs_ordered_v2.pdf',
#        path = 'results/figures/',
#        width = 170, height = 150, units = 'mm')
  
## Not ordered by relative rho value simply by the distribution ----
data_causal_relationships |>
  colnames()

plot_ccm_asvs_causally_affecting_other_asvs_ordered <- data_causal_relationships |>
  dplyr::group_by(order_family_asv_num) |>
  dplyr::mutate(mean = mean(rho)) |>
  dplyr::ungroup() |>
  dplyr::mutate(bloom = case_when(asv_num_causal %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                  TRUE ~ 'no_bloomer')) |>
  ggplot(aes(fct_reorder(order_family_asv_num, mean), rho)) +
  # geom_point(aes(
  #                 alpha = if_else((asv_num_causal %in% bloo_taxonomy$asv_num_f), 
  #                                 0.6, 0.1)), position = position_jitter(width = 0.25),  
  #            size = 0.6)+
  geom_boxplot( aes(
                    alpha = if_else((asv_num_causal %in% bloo_taxonomy$asv_num_f), 0.8, 0.3)),
                notch = T, outliers = F)+ # I already plot the with the geom point
  coord_flip()+
  labs(y = 'Predictive skill (rho)',
       x = 'ASVs causally affectting other ASVs', 
       color = 'Fraction ASVs being causally affected')+
  facet_wrap(vars(bloom), labeller = labs_bloom)+
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
#         filename = 'plot_ccm_asvs_causally_affecting_other_asvs_ordered_v3.pdf',
#        path = 'results/figures/',
#        width = 170, height = 150, units = 'mm')

##
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
        panel.grid.minor  = element_blank(), text = element_text(size = 8),
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
### Important for the MDR analysis we need to mantain seasonality otherwise it could uncover interactions
## upload data ----
mdr_tb <- read.csv2('../EDM_carmen/MDR/Bl_nin120_cvunit0.025_aenet_jcof_Nmvx_Rallx_demo_seas.csv', header = T ) |> 
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
                                       "#F2A218"))(6)

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
  dplyr::select(asv_num, fraction, date, abundance_value, z_score_ra) |>
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
  dplyr::group_by(date, #fraction_eff, 
                  fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) #|>
  # dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
  #                                    str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv11_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on ASV11 (0.2)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv11_rclr_plot

# Create an empty plot for the placeholder
#empty_plot <- ggplot() + theme_void()

# Arrange the plots with an empty spot on the right side for the first plot
effect_community_on_ASV11 <- plot_grid(
# First row with empty space
  asv11_rclr_plot,
  asv11_affectedby_community_plot, # Second row
  rel_heights = c(1, 1.75), # Heights of the rows
  labels = c('F'), # Labels for the plots
  label_size = 10,
  label_fontface = 'plain',
  ncol = 1 # Arrange in one column
)

# ggsave(
#   plot = effect_community_on_ASV11, 
#   filename = 'effect_community_on_ASV11_v2.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 150, 
#   height = 150, 
#   units = 'mm'
# )

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
  dplyr::group_by(date, #fraction_eff, 
                  fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) #|>
  # dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
  #                                           str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv23_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction)#~fraction_eff_ed
             , labeller = labs_fraction#, dir = 'v'
             ) + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv23 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

effect_community_on_ASV23 <- plot_grid(asv23_rclr_plot, 
          asv23_affectedby_community_plot,
          rel_heights = c(1, 2),
          labels = c('B'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

effect_community_on_ASV23

# ggsave(
#   plot = effect_community_on_ASV23, 
#   filename = 'effect_community_on_ASV23_v2.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 180, 
#   units = 'mm'
# )
 
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
  dplyr::group_by(date, #fraction_eff, 
                  fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) #|>
  # dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
  #                                           str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv7_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction)
             #~fraction_eff_ed
             , labeller = labs_fraction#, dir = 'v'
             ) + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv7 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_line()+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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

# ggsave(
#   plot = effect_community_on_ASV7, 
#   filename = 'effect_community_on_ASV7_v2.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 180, 
#   units = 'mm'
# )

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
  dplyr::group_by(date, #fraction_eff, 
                  fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) #|>
  # dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
  #                                           str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv1_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction)
             #~fraction_eff_ed
             , labeller = labs_fraction
             #, dir = 'v'
             ) + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv1 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv1_rclr_plot

effect_community_on_ASV1 <- plot_grid(asv1_rclr_plot, 
          asv1_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('E'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

# ggsave(
#   plot = effect_community_on_ASV1, 
#   filename = 'effect_community_on_ASV1_v2.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 180, 
#   units = 'mm'
# )

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
  dplyr::group_by(date, #fraction_eff, 
                  fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg)# |>
  # dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
  #                                           str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv22_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction)
             #~fraction_eff_ed
             , labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv22 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_line()+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv22_rclr_plot

effect_community_on_ASV22 <- plot_grid(asv22_rclr_plot, 
          asv22_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('D'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

# ggsave(
#   plot = effect_community_on_ASV22, 
#   filename = 'effect_community_on_ASV22_V2.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 150, 
#   height = 180, 
#   units = 'mm'
# )

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
  dplyr::group_by(date, #fraction_eff, 
                  fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) #|>
  # dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
  #                                           str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv15_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction),
             #~fraction_eff_ed,
             labeller = labs_fraction#, dir = 'v'
             ) + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv15 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv15_rclr_plot

effect_community_on_ASV15 <- plot_grid(asv15_rclr_plot, 
          asv15_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('C', 'D'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

#  ggsave(
#   plot = effect_community_on_ASV15, 
#   filename = 'effect_community_on_ASV15_v2.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 180, 
#   units = 'mm'
# )

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
  dplyr::group_by(date, #fraction_eff, 
                  fraction_b, fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) #|>
  # dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
  #                                           str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv31_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction),
             #~fraction_eff_ed, 
             labeller = labs_fraction#, dir = 'v'
             ) + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv31 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_line()+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction, dir = 'v') + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv31_rclr_plot

effect_community_on_ASV31 <- plot_grid(asv31_rclr_plot, 
          asv31_affectedby_community_plot,
          rel_heights = c(1, 1.75),
          labels = c('C'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column
# 
# ggsave(
#   plot = effect_community_on_ASV31, 
#   filename = 'effect_community_on_ASV31_V2.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 180, 
#   units = 'mm'
# )

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
  dplyr::group_by(date, #fraction_eff, # account for all taxa affecting my bloomer
                  fraction_b, 
                  fraction, asv_num) |> 
  dplyr::filter(value != 0) |> 
  dplyr::reframe(n  = n(),
                 interactions_pos = sum(value[value > 0], na.rm = TRUE),
                 interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  dplyr::mutate(balance = interactions_pos+interaction_neg) #|>
  # dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
  #                                           str_detect(fraction_eff, 'bn') ~ '3'))
data <- data |> 
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))

asv17_affectedby_community_plot <- data |>     
  ggplot(aes(date, interactions_pos)) + 
  geom_col(aes(y = interactions_pos), fill = '#000000') + 
  geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
  facet_wrap(vars(fraction)
             #~fraction_eff_ed
             , labeller = labs_fraction#, dir = 'v'
             ) + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on asv17 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
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
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_line()+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv17_rclr_plot

effect_community_on_ASV17 <- plot_grid(asv17_rclr_plot, 
          asv17_affectedby_community_plot,
          rel_heights = c(1, 1),
          labels = c('A'), # Labels for the plots
          label_size = 10,
          label_fontface = 'plain',
          ncol = 1) # Arrange in one column

# ggsave(
#   plot = effect_community_on_ASV17, 
#   filename = 'effect_community_on_ASV17_v2.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 180, 
#   units = 'mm'
# )

### summary of interactions during the time series -----
composition_bloom_events_interactions_plots <- plot_grid(effect_community_on_ASV17,
                                                         effect_community_on_ASV31,
                                                         effect_community_on_ASV23,
                                                         effect_community_on_ASV22,
                                                         effect_community_on_ASV1,
                                                         effect_community_on_ASV11,
                                                         rel_widths = c(1, 0.6),
                                                         ncol = 2)
  
composition_bloom_events_interactions_plots

 ggsave(
    plot = composition_bloom_events_interactions_plots, 
    filename = 'composition_bloom_events_interactions_plots_v2.pdf',
    path = 'results/figures/MDR/v2/',
    width = 180, 
    height = 260, 
    units = 'mm'
  )
 
 composition_bloom_events_interactions_sup_plots <- plot_grid(effect_community_on_ASV7,
                                                          effect_community_on_ASV15,
                                                          rel_widths = c(1, 1),
                                                          ncol = 2)
 
 composition_bloom_events_interactions_sup_plots
 
 ggsave(
   plot = composition_bloom_events_interactions_sup_plots, 
   filename = 'composition_bloom_events_interactions_sup_plots_v2.pdf',
   path = 'results/figures/MDR/v2/',
   width = 180, 
   height = 150, 
   units = 'mm'
 )

### HOW ARE BLOOMERS AFFECTING THE COMMUNITY? ----
 ### ASV15 -----
 mdr_tb_m_rclr <- mdr_tb_m |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
   #dplyr::filter(asv_num == 'asv15') |>
   pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
   separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
   left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
   dplyr::filter(asv_num_eff == 'asv15') |>
   dplyr::filter(asv_num_eff != asv_num) 
 
 asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
   dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
   dplyr::filter(asv_num == 'asv15')
 
 data <- mdr_tb_m_rclr |> 
   dplyr::mutate(value = as.numeric(value)) |> 
   dplyr::group_by(date, fraction_eff, 
                   #fraction_b, fraction, 
                   asv_num_eff) |> 
   dplyr::filter(value != 0) |> 
   dplyr::reframe(n  = n(),
                  interactions_pos = sum(value[value > 0], na.rm = TRUE),
                  interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
   dplyr::mutate(balance = interactions_pos+interaction_neg) |>
 dplyr::mutate(fraction_eff_ed = case_when(str_detect(fraction_eff, 'bp') ~ '0.2',
                                           str_detect(fraction_eff, 'bn') ~ '3'))
 data <- data |> 
   dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009'))
 
 asv15_affectting_community_plot <- data |>     
   ggplot(aes(date, interactions_pos)) + 
   geom_col(aes(y = interactions_pos), fill = '#000000') + 
   geom_col(aes(y = interaction_neg), fill = '#8C0009') + 
   geom_point(aes(y = balance, color = balance_color), size = 0.75, alpha = 0.6) + 
   facet_wrap(vars(fraction_eff_ed),
              #~fraction_eff_ed,
              labeller = labs_fraction#, dir = 'v'
   ) + 
   scale_color_identity()+
   theme_bw() + 
   labs(x = 'Date', y = 'Interaction Strength', color = 'Order', title = 'The microbial community being affected by ASV15 (0.2 & 3)') + 
   theme(strip.background = element_rect(fill = 'transparent'),
         text = element_text(size = 8),
         axis.text.x = element_text(size = 10),
         panel.border = element_blank(),
         title = element_text(size = 8))
 
 asv15_affectting_community_plot 
 
 asv_tab_all_bloo_z_tax_red_asv |>
   colnames()
 
 asv15_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
   left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
   ggplot(aes(date, abundance_value, color = order_f)) + 
   #geom_point( size = 0.75, alpha = 0.6) + 
   geom_line()+
   geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
   geom_area(aes(fill = order_f), alpha = 0.5)+
   geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
   scale_shape_identity()+
   scale_fill_manual(values = palette_order_assigned_bloo)+
   scale_color_manual(values = palette_order_assigned_bloo)+
   facet_wrap(vars(fraction), labeller = labs_fraction) + 
   theme_bw() + 
   labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
   theme(strip.background = element_rect(fill = 'transparent'),
         text = element_text(size = 8),
         axis.text.x = element_text(size = 10),
         panel.border = element_blank(),
         title = element_text(size = 8),
         legend.position = 'none')
 
 asv15_rclr_plot
 
 effect_community_on_ASV15 <- plot_grid(asv15_rclr_plot, 
                                        asv15_affectting_community_plot,
                                        rel_heights = c(1, 1.75),
                                        labels = c('C', 'D'), # Labels for the plots
                                        label_size = 10,
                                        label_fontface = 'plain',
                                        ncol = 1) # Arrange in one column
 
 
 effect_community_on_ASV15
 
 
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
  scale_y_continuous(limits = c(-0.2, 0.5))+
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
  filename = 'interaction_asv7_and_asv15_v2.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

ggsave(
  plot = asv7_affectedby_asv15_plot , 
  filename = 'asv7_affectedby_asv15_plot_v2.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 100, 
  units = 'mm'
)

## ----- Different interactions during harbor restoration ? ## -----
### ASV11 -----
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
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data$harbor_restoration <- factor(data$harbor_restoration,
                                  levels = c('pre_perturbation',
                                             'perturbation',
                                             'after_perturbation'))

asv11_affectedby_community_harbor_plot <- data |>   
  ggplot(aes(harbor_restoration, interactions_pos)) + 
  geom_boxplot(aes(), width = 0.25)+
  geom_point(aes(y = balance, color = balance_color), size = 0.75, 
             alpha = 0.6, position = position_jitter(0.25)) + 
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order', title = 'The microbial community effect on ASV11 (0.2)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv11_affectedby_community_harbor_plot 

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
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data$harbor_restoration <- factor(data$harbor_restoration,
                                  levels = c('pre_perturbation',
                                             'perturbation',
                                             'after_perturbation'))

asv23_affectedby_community_harbor_plot <- data |>     
  ggplot(aes(harbor_restoration, interactions_pos)) + 
  geom_boxplot(aes(), width = 0.25)+
  geom_point(aes(y = balance, color = balance_color), size = 0.75, 
             alpha = 0.6, position = position_jitter(0.25)) + 
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV23 (0.2)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv23_affectedby_community_harbor_plot 

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
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data$harbor_restoration <- factor(data$harbor_restoration,
                                  levels = c('pre_perturbation',
                                             'perturbation',
                                             'after_perturbation'))


asv7_affectedby_community_harbor_plot <- data |>     
  ggplot(aes(harbor_restoration, interactions_pos)) + 
  geom_boxplot(aes(), width = 0.25)+
  geom_point(aes(y = balance, color = balance_color), size = 0.75, 
             alpha = 0.6, position = position_jitter(0.25)) + 
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV7 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv7_affectedby_community_harbor_plot 

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
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data$harbor_restoration <- factor(data$harbor_restoration,
                                  levels = c('pre_perturbation',
                                             'perturbation',
                                             'after_perturbation'))

asv1_affectedby_community_harbor_plot <- data |>     
  ggplot(aes(harbor_restoration, interactions_pos)) + 
  geom_boxplot(aes(), width = 0.25)+
  geom_point(aes(y = balance, color = balance_color), size = 0.75, 
             alpha = 0.6, position = position_jitter(0.25)) + 
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV1 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv1_affectedby_community_harbor_plot 

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
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data$harbor_restoration <- factor(data$harbor_restoration,
                                  levels = c('pre_perturbation',
                                             'perturbation',
                                             'after_perturbation'))

asv22_affectedby_community_harbor_plot <- data |>     
  ggplot(aes(harbor_restoration, interactions_pos)) + 
  geom_boxplot(aes(), width = 0.25)+
  geom_point(aes(y = balance, color = balance_color), size = 0.75, 
             alpha = 0.6, position = position_jitter(0.25)) + 
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV22 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))
asv22_affectedby_community_harbor_plot 

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
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data$harbor_restoration <- factor(data$harbor_restoration,
                                  levels = c('pre_perturbation',
                                             'perturbation',
                                             'after_perturbation'))

asv15_affectedby_community_harbor_plot <- data |>     
  ggplot(aes(harbor_restoration, interactions_pos)) + 
  geom_boxplot(aes(), width = 0.25)+
  geom_point(aes(y = balance, color = balance_color), size = 0.75, 
             alpha = 0.6, position = position_jitter(0.25)) + 
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV15 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv15_affectedby_community_harbor_plot 

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
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data$harbor_restoration <- factor(data$harbor_restoration,
                                  levels = c('pre_perturbation',
                                             'perturbation',
                                             'after_perturbation'))

asv31_affectedby_community_harbor_plot <- data |>     
  ggplot(aes(harbor_restoration, interactions_pos)) + 
  geom_boxplot(aes(), width = 0.25)+
  geom_point(aes(y = balance, color = balance_color), size = 0.75, 
             alpha = 0.6, position = position_jitter(0.25)) + 
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV31 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv31_affectedby_community_harbor_plot 

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
  dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data$harbor_restoration <- factor(data$harbor_restoration,
                                  levels = c('pre_perturbation',
                                             'perturbation',
                                             'after_perturbation'))

asv17_affectedby_community_harbor_plot <- data |>     
  ggplot(aes(harbor_restoration, interactions_pos)) + 
  geom_point(aes(y = balance, color = balance_color), size = 0.75, 
             alpha = 0.6, position = position_jitter(0.25)) + 
  geom_boxplot(aes(), width = 0.25)+
  facet_wrap(vars(fraction_eff_ed), labeller = labs_fraction, dir = 'v') + 
  scale_color_identity()+
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV17 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8))

asv17_affectedby_community_harbor_plot 

## composition harbor interactions ----
harbor_interactions_plot <- plot_grid(asv11_affectedby_community_harbor_plot,
                                      asv23_affectedby_community_harbor_plot,
                                      asv7_affectedby_community_harbor_plot,
                                      asv31_affectedby_community_harbor_plot,
                                      asv22_affectedby_community_harbor_plot,
                                      asv1_affectedby_community_harbor_plot,
                                      asv15_affectedby_community_harbor_plot,
                                      asv17_affectedby_community_harbor_plot,
                                       rel_heights = c(1, 1),
                                       labels = c('A', 'B', 'C', 'D',
                                                  'E', 'F', 'G', 'H'), # Labels for the plots
                                       label_size = 10,
                                       label_fontface = 'plain',
                                       ncol = 4) # Arrange in one column
harbor_interactions_plot

ggsave(
  plot = harbor_interactions_plot, 
  filename = 'harbor_interactions_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 180, 
  units = 'mm'
)

## Different approach during the harbor restoration ----
### Did the community interactions changed in the PA and FL fractions? ----
#### Plus, is there any difference between bloomers and non-blooming taxa?
### Mean interaction strength ------
#### +0 to sample id_num abundance 
mdr_tb_m_rclr <- mdr_tb_m |>
  # dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
  #                              'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/n,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/n,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
 dplyr::mutate(bloomer = case_when(asv_num %in% c('asv1', 'asv7', 'asv11',
                                                  'asv22', 'asv23', 'asv17',
                                                  'asv31', 'asv15') ~ 'bloomer',
                                   TRUE ~ 'no-bloomer')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

data <- data |>
  pivot_longer(cols = starts_with('mean')) 

summary_interactions_boxplot_community_during_harbor_restoration <- data |>
  dplyr::filter(!is.na(fraction)) |>
  ggplot(aes(interaction( name, harbor_restoration), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name), alpha = 0.7)+
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(fraction~bloomer)+ #, labeller = labs_fraction
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'The microbial community effect on blooming taxa
       during the harbor restoration') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_during_harbor_restoration

# ggsave(
#   plot = summary_interactions_boxplot_community_during_harbor_restoration , 
#   filename = 'summary_interactions_boxplot_community_during_harbor_restoration.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 100, 
#   height = 150, 
#   units = 'mm'
# )

## in 0 and 1 
mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ 1,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/1,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ 1,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/1,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

data |>
  group_by(date, harbor_restoration, fraction) |>
  distinct(date, harbor_restoration, fraction) |>
  group_by(harbor_restoration, fraction) |>
  dplyr::reframe(n = n())

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% c('asv1', 'asv7', 'asv11',
                                                   'asv22', 'asv23', 'asv17',
                                                   'asv31', 'asv15') ~ 'bloomer',
                                    TRUE ~ 'no-bloomer')) |>
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  )) |>
  dplyr::mutate(samples_period = case_when(harbor_restoration == 'pre_perturbation' ~ 75,
                                           harbor_restoration == 'perturbation' ~ 26 ,
                                           harbor_restoration == 'after_perturbation' ~ 18)) 

data <- data |>
  pivot_longer(cols = starts_with('mean')) |>
  dplyr::mutate(value = value/samples_period)

summary_interactions_boxplot_community_during_harbor_restoration <- data |>
  dplyr::filter(!is.na(fraction)) |>
  ggplot(aes(interaction( harbor_restoration, name), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name), alpha = 0.7)+
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(fraction~bloomer)+ #, labeller = labs_fraction
  theme_bw() + 
  labs(x = '', y = 'Interaction (0,1) / nº samples considered', color = 'Interaction',
       fill = 'Interaction',
       title = 'The microbial community effect on blooming taxa
       during the harbor restoration') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_during_harbor_restoration

# ggsave(
#   plot = summary_interactions_boxplot_community_during_harbor_restoration , 
#   filename = 'summary_interactions_boxplot_community_during_harbor_restoration_binary.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 100, 
#   height = 150, 
#   units = 'mm'
# )

## I do not separate bloomers and non bloomers ----

summary_interactions_boxplot_community_during_harbor_restoration <- data |>
  dplyr::filter(!is.na(fraction)) |>
  ggplot(aes(interaction( harbor_restoration, name), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name), alpha = 0.7)+
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(fraction), labeller = labs_fraction)+ #
  theme_bw() + 
  labs(x = '', y = 'Interaction (0,1) / nº samples considered', color = 'Interaction',
       fill = 'Interaction',
       title = 'The microbial community interactions
       during the harbor restoration') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_during_harbor_restoration

ggsave(
  plot = summary_interactions_boxplot_community_during_harbor_restoration , 
  filename = 'summary_interactions_boxplot_community_during_harbor_restoration_binary_grouped.pdf',
  path = 'results/figures/MDR/v2/',
  width = 120, 
  height = 100, 
  units = 'mm'
)

## MDR differences between bloom and no-bloom ----
### number of interactions 
### 1, 0 , -1 
## sum + and - interactions 
### ASV17 -----
mdr_tb_m |>
  dim()

mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv17') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr$date |>
  unique()

mdr_tb_m_rclr  |>
  colnames()

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = case_when(value > 0 ~ 1,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(pos_interactions = sum(value))

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = case_when(value < 0 ~ 1,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value))

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)
  
asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv17')

bloom_event <- asv_tab_all_bloo_z_tax_red_asv |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                          abundance_value > 0.1) ~ 'bloom',
                TRUE ~ 'no-bloom')) |>
  #dplyr::filter(bloom_event == 'bloom') |>
  dplyr::select(bloom_event, date, fraction) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

 data <- mdr_tb_m_rclr_int |> 
#   dplyr::mutate(value = as.numeric(value))# |> 
  #dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  #dplyr::filter(value != 0) |> 
  # dplyr::reframe(n  = n(),
  #                interactions_pos = sum(value[value > 0], na.rm = TRUE),
  #                interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  # #dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                            str_detect(fraction_b, 'bn') ~ '3'))

# Convert bloom_event$date to POSIXct
bloom_event$date <- as.POSIXct(bloom_event$date, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")

# Convert data$date to POSIXct
data$date <- as.POSIXct(data$date, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")

palette_interactions <- c('neg_interactions' = '#8C0009',
  'pos_interactions' = '#000000')

unique(data$name)

data <- data |>
  #dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d")) |>
  #dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  left_join(bloom_event, by = c('date', 'fraction'))

data <- data |>
  dplyr::mutate(neg_interactions = -neg_interactions) |>
  pivot_longer(cols = matches('pos|neg'))

interactions_01_asv17_plot <- data |>
  ggplot(aes(date, value, group = name, color = name)) + 
  # geom_point(#aes(y = balance, color = balance_color), 
  #   size = 0.75, 
  #            alpha = 0.6) + 
  geom_line(aes(group = name))+
  scale_color_manual(values = c('pos_interactions' = '#000000',
                                'neg_interactions' = '#8C0009'))+
  #geom_boxplot(aes(), width = 0.25)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV17 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

interactions_01_asv17_plot

asv17_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_line()+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv17_rclr_plot

plot_grid(interactions_01_asv17_plot, 
          asv17_rclr_plot,
          cols = 1)

### is the j0 the eigenvalue? ---
mdr_tb_m |>
  dplyr::select(j0, time) |>
  ggplot(aes(as.numeric(time), as.numeric(j0)))+
  geom_point()

### the same without transforming to 1, 0 interactions -----
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num == 'asv17') |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/n,
                sd_pos_int = sd(mean_pos_int))

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/n,
                sd_neg_int = sd(mean_neg_int))

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

asv_tab_all_bloo_z_tax_red_asv <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(abundance_value_scaled = scale(abundance_value))|>
  dplyr::filter(asv_num == 'asv17')

bloom_event <- asv_tab_all_bloo_z_tax_red_asv |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

data <- mdr_tb_m_rclr_int |> 
  #   dplyr::mutate(value = as.numeric(value))# |> 
  #dplyr::group_by(date, fraction_eff, fraction_b, fraction, asv_num) |> 
  #dplyr::filter(value != 0) |> 
  # dplyr::reframe(n  = n(),
  #                interactions_pos = sum(value[value > 0], na.rm = TRUE),
  #                interaction_neg = sum(value[value < 0], na.rm = TRUE)) |>
  # #dplyr::mutate(balance = interactions_pos+interaction_neg) |>
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3'))

# Convert bloom_event$date to POSIXct
bloom_event$date <- as.POSIXct(bloom_event$date, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")

# Convert data$date to POSIXct
data$date <- as.POSIXct(data$date, 
                        format = "%Y-%m-%d %H:%M:%S", 
                        tz = "Europe/Berlin")

palette_interactions <- c('neg_interactions' = '#8C0009',
                          'pos_interactions' = '#000000')

data <- data |>
  #dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d")) |>
  #dplyr::mutate(balance_color = if_else(balance > 0, '#000000', '#8C0009')) |>
  left_join(bloom_event, by = c('date', 'fraction'))

data <- data |>
  dplyr::mutate(mean_neg_int = -mean_neg_int) |>
  pivot_longer(cols = starts_with('mean'))

data |>
  colnames()

data$name |>
  unique()

interactions_01_asv17_plot <- data |>
  ggplot(aes(date, value, group = name, color = name)) + 
  geom_point(#aes(y = balance, color = balance_color),
    size = 0.75,
             alpha = 0.6) +
  geom_line(aes(group = name))+
  scale_color_manual(values = c('pos_interactions' = '#000000',
                                'neg_interactions' = '#8C0009'))+
  #geom_boxplot(aes(), width = 0.25)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = '', y = 'Interaction Strength', color = 'Order',
       title = 'The microbial community effect on ASV17 (0.2 & 3)') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

interactions_01_asv17_plot

asv17_rclr_plot <- asv_tab_all_bloo_z_tax_red_asv |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  ggplot(aes(date, abundance_value, color = order_f)) + 
  #geom_point( size = 0.75, alpha = 0.6) + 
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.5)+
  geom_area(aes(fill = order_f), alpha = 0.5)+
  geom_line()+
  geom_point(aes(shape = if_else((abundance_value > 0.1 & z_score_ra > 2), 4, NA )))+
  scale_shape_identity()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  labs(x = 'Date', y = 'Relative Abundance', color = 'Order') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        panel.border = element_blank(),
        title = element_text(size = 8),
        legend.position = 'none')

asv17_rclr_plot

plot_grid(interactions_01_asv17_plot, 
          asv17_rclr_plot,
          cols = 1)

# HOW IS THE COMMUNITY INTERACTING WITH MY APPARENT BLOOMERS OVER THE TIME SERIES? -----
### abundance at t and interaction t -> t+1
### in this case I'm considering as bloom event in the conditions fulfilled in t
## ---- EVALUATE THE ABUNDANCE BEFORE THE INTERACTION ## ----
### Mean interaction strength ------
#### t+0 to sample id_num abundance 
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/n,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/n,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  left_join(bloom_event , by = c('fraction', 'asv_num', 'time' = 'sample_id_num'))

n_bloom_events <- bloom_event |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::reframe(n = paste0('n=', n()))

data <- data |>
  pivot_longer(cols = starts_with('mean')) |>
  left_join(n_bloom_events)

interactions_pre_mean_plot <- data |>
  ggplot(aes(as.numeric(time), interaction( name,fraction,  asv_num), fill = value)) + 
  geom_tile()+
  geom_point(data = data |>
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), interaction( name,fraction,  asv_num)),
             size = 0.75,
             alpha = 0.6, 
             shape = 4) +
  scale_fill_gradient2(low = "#8C0009", mid = "white", 
                       high = "#000000", midpoint = 0, na.value = 'white')+
  theme_bw() + 
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9),
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'))+
  labs(x = 'Time', y  = 'Blooming ASVs being affected by the community', fill = 'Mean Interaction Strength',
       title = 'The microbial community effect on blooming taxa 
       and the pre abundance after interaction') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 5),
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 4),
        axis.ticks.length = unit(0.2, "mm"))

interactions_pre_mean_plot

ggsave(
  plot = interactions_pre_mean_plot , 
  filename = 'community_interactions_pre_mean_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 100, 
  units = 'mm'
)

## group all those data points in which there is a blooming event vs. those that there is not 
genus_data <- tax_bbmo_10y_new |>
  dplyr::filter(asv_num %in% unique(bloo_taxonomy$asv_num_f)) |>
  dplyr::distinct(asv_num, family, genus)

data <- data |>
  dplyr::filter(!is.na(bloom_event)) |> # missing and non interpolated values
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) 

custom_labels <- data |>
  distinct(asv_num, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num),
      is.na(genus) ~ paste0(family_f, " ", asv_num),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num)
    )
  ) |> 
  pull(label, name = asv_num) 

summary_interactions_boxplot_community_pre_int <- data |>
  ggplot(aes(interaction( bloom_event, name, fraction), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name))+
  geom_text(aes(label = n), 
            position = position_nudge(y = max(abs(data$value), na.rm = TRUE) * 0.1), 
            size = 1, vjust = 0.5,hjust = 1,  check_overlap = T) + 
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(asv_num), labeller = as_labeller(custom_labels)#)+
             ,ncol = 2)+
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'The microbial community effect on blooming taxa
       pre interaction') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_pre_int

ggsave(
  plot = summary_interactions_boxplot_community_pre_int , 
  filename = 'community_summary_interactions_boxplot_community_pre_int.pdf',
  path = 'results/figures/MDR/v2/',
  width = 100, 
  height = 150, 
  units = 'mm'
)

### Sum interactions transformed to 0 and 1 -----
#### +0 to sample id_num abundance 
mdr_tb_m_rclr <- mdr_tb_m |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ 1,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/1,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 1 ~ -1,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/1,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  left_join(bloom_event , by = c('fraction', 'asv_num', 'time' = 'sample_id_num'))

n_bloom_events <- bloom_event |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::reframe(n = paste0('n=', n()),
                 n_blooms = n()) 

data <- data |>
  left_join(n_bloom_events) |>
  dplyr::mutate(dates_considered = case_when(bloom_event == 'bloom' ~ n_blooms,
                                             TRUE ~ 120-n_blooms)) |>
  pivot_longer(cols = starts_with('mean')) |>
  dplyr::mutate(value = value/dates_considered)

interactions_pre_mean_plot <- data |>
  dplyr::filter(!is.na(bloom_event)) |>
  ggplot(aes(as.numeric(time), interaction(name,fraction,  asv_num), fill = value)) + 
  geom_tile()+
  geom_point(data = data |>
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), interaction( name,fraction,  asv_num)),
             size = 0.75,
             alpha = 0.6, 
             shape = 4) +
  scale_fill_gradient2(low = "#8C0009", mid = "white", 
                       high = "#000000", midpoint = 0, na.value = 'white')+
  theme_bw() + 
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9),
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'))+
  labs(x = 'Time', y  = 'Blooming ASVs being affected by the community', fill = 'Mean Interaction Strength',
       title = 'The microbial community effect on blooming taxa 
       and the pre abundance after interaction') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 5),
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 4),
        axis.ticks.length = unit(0.2, "mm"))

interactions_pre_mean_plot

ggsave(
  plot = interactions_pre_mean_plot , 
  filename = 'community_interactions_pre_mean_binary_norm_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 100, 
  units = 'mm'
)

## group all those data points in which there is a blooming event vs. those that there is not 
data$name |>
  unique()

genus_data <- tax_bbmo_10y_new |>
  dplyr::filter(asv_num %in% unique(bloo_taxonomy$asv_num_f)) |>
  dplyr::distinct(asv_num, family, genus)

data <- data |>
  dplyr::filter(!is.na(bloom_event)) |> # missing and non interpolated values
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) 

custom_labels <- data |>
  distinct(asv_num, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num),
      is.na(genus) ~ paste0(family_f, " ", asv_num),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num)
    )
  ) |> 
  pull(label, name = asv_num) 

summary_interactions_boxplot_community_pre_int <- data |>
  ggplot(aes(interaction( bloom_event, name, fraction), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name))+
  geom_text(aes(label = n), 
            position = position_nudge(y = max(abs(data$value), na.rm = TRUE) * 0.1), 
            size = 1, vjust = 0.5,hjust = 1,  check_overlap = T) + 
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(asv_num), labeller = as_labeller(custom_labels)#)+
             ,ncol = 2,
             scales = 'free_x')+
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'The microbial community effect on blooming taxa
       pre interaction') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_pre_int

ggsave(
  plot = summary_interactions_boxplot_community_pre_int , 
  filename = 'community_summary_interactions_boxplot_community_norm_pre_int.pdf',
  path = 'results/figures/MDR/v2/',
  width = 100, 
  height = 150, 
  units = 'mm'
)

## ---- EVALUATE RESULTING ABUNDANCE AFTER INTERACTION ##  ----
### Mean interaction strength ------
#### +1 to sample id_num abundance 
mdr_tb_m_rclr <- mdr_tb_m |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/n,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/n,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)+1) |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  left_join(bloom_event , by = c('fraction', 'asv_num', 'time' = 'sample_id_num'))

n_bloom_events <- bloom_event |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::reframe(n = paste0('n=', n()))

data <- data |>
  pivot_longer(cols = starts_with('mean')) |>
  left_join(n_bloom_events)

interactions_resulting_mean_plot <- data |>
  ggplot(aes(as.numeric(time), interaction( name,fraction,  asv_num), fill = value)) + 
  geom_tile()+
  geom_point(data = data |>
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), interaction( name,fraction,  asv_num)),
             size = 0.75,
             alpha = 0.6, 
             shape = 4) +
  scale_fill_gradient2(low = "#8C0009", mid = "white", 
                       high = "#000000", midpoint = 0, na.value = 'white')+
  theme_bw() + 
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9),
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'))+
  labs(x = 'Time', y  = 'Blooming ASVs being affected by the community', fill = 'Mean Interaction Strength',
       title = 'The microbial community effect on blooming taxa 
       and the resulting abundance after interaction') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 5),
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 4),
        axis.ticks.length = unit(0.2, "mm"))

interactions_resulting_mean_plot

ggsave(
  plot = interactions_resulting_mean_plot , 
  filename = 'community_interactions_resulting_mean_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 100, 
  units = 'mm'
)

 ## group all those data points in which there is a blooming event vs. those that there is not 
data$name |>
  unique()

genus_data <- tax_bbmo_10y_new |>
  dplyr::filter(asv_num %in% unique(bloo_taxonomy$asv_num_f)) |>
  dplyr::distinct(asv_num, family, genus)

data <- data |>
  dplyr::filter(!is.na(bloom_event)) |> # missing and non interpolated values
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) 

custom_labels <- data |>
  distinct(asv_num, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num),
      is.na(genus) ~ paste0(family_f, " ", asv_num),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num)
    )
  ) |> 
  pull(label, name = asv_num) 

data |>
  str()

summary_interactions_boxplot_community_resulting_int <- data |>
  ggplot(aes(interaction( bloom_event, name, fraction), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name))+
  geom_text(aes(label = n), 
            position = position_nudge(y = max(abs(data$value), na.rm = TRUE) * 0.1), 
            size = 1, vjust = 0.5,hjust = 1,  check_overlap = T) + 
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(asv_num), labeller = as_labeller(custom_labels)#)+
               ,ncol = 2)+
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'The microbial community effect on blooming taxa
       resulting interaction') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_resulting_int

ggsave(
  plot = summary_interactions_boxplot_community_resulting_int , 
  filename = 'community_summary_interactions_boxplot_community_resulting_int.pdf',
  path = 'results/figures/MDR/v2/',
  width = 100, 
  height = 150, 
  units = 'mm'
)

### Sum interactions transformed to 0 and 1 -----
#### +1 to sample id_num abundance 
mdr_tb_m_rclr <- mdr_tb_m |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ 1,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/1,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 1 ~ -1,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/1,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)+1) |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  left_join(bloom_event , by = c('fraction', 'asv_num', 'time' = 'sample_id_num'))

n_bloom_events <- bloom_event |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::reframe(n = paste0('n=', n()))

data <- data |>
  pivot_longer(cols = starts_with('mean')) |>
  left_join(n_bloom_events)

interactions_resulting_mean_plot <- data |>
  dplyr::filter(!is.na(bloom_event)) |>
  ggplot(aes(as.numeric(time), interaction(name,fraction,  asv_num), fill = value)) + 
  geom_tile()+
  geom_point(data = data |>
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), interaction( name,fraction,  asv_num)),
             size = 0.75,
             alpha = 0.6, 
             shape = 4) +
  scale_fill_gradient2(low = "#8C0009", mid = "white", 
                       high = "#000000", midpoint = 0, na.value = 'white')+
  theme_bw() + 
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9),
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'))+
  labs(x = 'Time', y  = 'Blooming ASVs being affected by the community', fill = 'Mean Interaction Strength',
       title = 'The microbial community effect on blooming taxa 
       and the resulting abundance after interaction') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 5),
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 4),
        axis.ticks.length = unit(0.2, "mm"))

interactions_resulting_mean_plot

ggsave(
  plot = interactions_resulting_mean_plot , 
  filename = 'community_interactions_resulting_mean_binary_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 100, 
  units = 'mm'
)

## group all those data points in which there is a blooming event vs. those that there is not 
data$name |>
  unique()

genus_data <- tax_bbmo_10y_new |>
  dplyr::filter(asv_num %in% unique(bloo_taxonomy$asv_num_f)) |>
  dplyr::distinct(asv_num, family, genus)

data <- data |>
  dplyr::filter(!is.na(bloom_event)) |> # missing and non interpolated values
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) 

custom_labels <- data |>
  distinct(asv_num, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num),
      is.na(genus) ~ paste0(family_f, " ", asv_num),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num)
    )
  ) |> 
  pull(label, name = asv_num) 

data |>
  str()

summary_interactions_boxplot_community_resulting_int <- data |>
  ggplot(aes(interaction( bloom_event, name, fraction), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name))+
  geom_text(aes(label = n), 
            position = position_nudge(y = max(abs(data$value), na.rm = TRUE) * 0.1), 
            size = 1, vjust = 0.5,hjust = 1,  check_overlap = T) + 
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(asv_num), labeller = as_labeller(custom_labels)#)+
             ,ncol = 2)+
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'The microbial community effect on blooming taxa
       resulting interaction') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_resulting_int

ggsave(
  plot = summary_interactions_boxplot_community_resulting_int , 
  filename = 'community_summary_interactions_boxplot_community_resulting_int.pdf',
  path = 'results/figures/MDR/v2/',
  width = 100, 
  height = 150, 
  units = 'mm'
)

# HOW ARE MY APPARENT BLOOMERS INTERACTING WITH THE COMMUNITY OVER THE TIME SERIES? -----
## ---- EVALUATE RESULTING ABUNDANCE PRE INTERACTION ##  ----
### abundance at t interaction t -> t+1
### Mean interaction strength ------
#### +0 to sample id_num abundance 
mdr_tb_m_rclr <- mdr_tb_m |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(asv_num_eff %in% c('asv17', 'asv22' , 'asv23' , 'asv11', 
                                   'asv1', 'asv31', 'asv7', 'asv15')) |>
  dplyr::filter(fraction_b %in% c('bn', 'bp'))

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/n,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/n,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |> ## add sample_number
  dplyr::select(-date)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  full_join(bloom_event , by = c('fraction', 'asv_num_eff' = 'asv_num', 'time' = 'sample_id_num')) 

data <- data |>
  pivot_longer(cols = starts_with('mean')) 

interactions_pre_mean_plot <- data |>
  ggplot(aes(as.numeric(time), interaction( name, fraction,  asv_num_eff), fill = value)) + 
  geom_tile()+
  geom_point(data = data |>
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), interaction( name,fraction,  asv_num_eff)),
             size = 0.75,
             alpha = 0.6, 
             shape = 4) +
  scale_fill_gradient2(low = "#8C0009", mid = "white", 
                       high = "#000000", midpoint = 0, na.value = 'white')+
  theme_bw() + 
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9),
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'))+
  labs(x = 'Time', y  = 'Blooming ASVs affectting the community',
       fill = 'Mean Interaction Strength',
       title = 'Blooming taxa effect on the community
       and the pre abundance before interaction') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 5),
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 4),
        axis.ticks.length = unit(0.2, "mm"))

interactions_pre_mean_plot

ggsave(
  plot = interactions_pre_mean_plot , 
  filename = 'interactions_pre_mean_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 100, 
  units = 'mm'
)

## group all those data points in which there is a blooming event vs. those that there is not 
genus_data <- tax_bbmo_10y_new |>
  dplyr::filter(asv_num %in% unique(bloo_taxonomy$asv_num_f)) |>
  dplyr::distinct(asv_num, family, genus)

data <- data |>
  dplyr::filter(!is.na(bloom_event)) |> # missing and non interpolated values
  left_join(bloo_taxonomy, by = c('asv_num_eff' = 'asv_num_f')) 

custom_labels <- data |>
  distinct(asv_num_eff, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num_eff),
      is.na(genus) ~ paste0(family_f, " ", asv_num_eff),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff) 

summary_interactions_boxplot_community_pre_int <- data |>
  ggplot(aes(interaction( bloom_event, name, fraction), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name))+
  # geom_text(aes(label = n), 
  #           position = position_nudge(y = max(abs(data$value), na.rm = TRUE) * 0.1), 
  #           size = 1, vjust = 0.5,hjust = 1,  check_overlap = T) + 
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(asv_num_eff), labeller = as_labeller(custom_labels)#)+
             ,ncol = 2)+
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'Blooming taxa effect on the community
       and the pre abundance before interaction') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_pre_int

ggsave(
  plot = summary_interactions_boxplot_community_pre_int , 
  filename = 'summary_interactions_boxplot_community_pre_int.pdf',
  path = 'results/figures/MDR/v2/',
  width = 100, 
  height = 150, 
  units = 'mm'
)

### Sum interactions transformed to 0 and 1 ----
mdr_tb_m_rclr <- mdr_tb_m |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(asv_num_eff %in% c('asv17', 'asv22' , 'asv23' , 'asv11', 
                                   'asv1', 'asv31', 'asv7', 'asv15')) |>
  dplyr::filter(fraction_b %in% c('bn', 'bp'))

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ 1,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/1,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ 1,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/1,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |> ## add sample_number
  dplyr::select(-date)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  full_join(bloom_event , by = c('fraction', 'asv_num_eff' = 'asv_num', 'time' = 'sample_id_num')) 

data <- data |>
  dplyr::mutate(mean_neg_int = -mean_neg_int) |>
  pivot_longer(cols = starts_with('mean')) 

interactions_pre_mean_plot <- data |>
  ggplot(aes(as.numeric(time), interaction( name, fraction,  asv_num_eff), fill = value)) + 
  geom_tile()+
  geom_point(data = data |>
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), interaction( name,fraction,  asv_num_eff)),
             size = 0.75,
             alpha = 0.6, 
             shape = 4) +
  scale_fill_gradient2(low = "#8C0009", mid = "white", 
                       high = "#000000", midpoint = 0, na.value = 'white')+
  theme_bw() + 
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9),
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'))+
  labs(x = 'Time', y  = 'Blooming ASVs affectting the community',
       fill = 'Mean Interaction Strength',
       title = 'Blooming taxa effect on the community
       and the pre abundance before interaction') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 5),
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 4),
        axis.ticks.length = unit(0.2, "mm"))

interactions_pre_mean_plot

ggsave(
  plot = interactions_pre_mean_plot , 
  filename = 'interactions_pre_mean_binary_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 100, 
  units = 'mm'
)

## group all those data points in which there is a blooming event vs. those that there is not 
genus_data <- tax_bbmo_10y_new |>
  dplyr::filter(asv_num %in% unique(bloo_taxonomy$asv_num_f)) |>
  dplyr::distinct(asv_num, family, genus)

data <- data |>
  dplyr::filter(!is.na(bloom_event)) |> # missing and non interpolated values
  left_join(bloo_taxonomy, by = c('asv_num_eff' = 'asv_num_f')) 

custom_labels <- data |>
  distinct(asv_num_eff, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num_eff),
      is.na(genus) ~ paste0(family_f, " ", asv_num_eff),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff) 

summary_interactions_boxplot_community_pre_int <- data |>
  ggplot(aes(interaction( bloom_event, name, fraction), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name))+
  # geom_text(aes(label = n), 
  #           position = position_nudge(y = max(abs(data$value), na.rm = TRUE) * 0.1), 
  #           size = 1, vjust = 0.5,hjust = 1,  check_overlap = T) + 
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(asv_num_eff), labeller = as_labeller(custom_labels)#)+
             ,ncol = 2)+
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'Blooming taxa effect on the community
       and the pre abundance before interaction') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_pre_int

ggsave(
  plot = summary_interactions_boxplot_community_pre_int , 
  filename = 'summary_interactions_boxplot_community_binary_pre_int.pdf',
  path = 'results/figures/MDR/v2/',
  width = 100, 
  height = 150, 
  units = 'mm'
)

## ---- EVALUATE RESULTING ABUNDANCE AFTER INTERACTION ##  ----
### abundance at t+1 interaction t -> t+1
### Mean interaction strength ------
#### +1 to sample id_num abundance 
mdr_tb_m_rclr <- mdr_tb_m |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(asv_num_eff %in% c('asv17', 'asv22' , 'asv23' , 'asv11', 
                                   'asv1', 'asv31', 'asv7', 'asv15')) |>
  dplyr::filter(fraction_b %in% c('bn', 'bp'))

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/n,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/n,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)+1) |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |> ## add sample_number
  dplyr::select(-date)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  full_join(bloom_event , by = c('fraction', 'asv_num_eff' = 'asv_num', 'time' = 'sample_id_num')) 

data <- data |>
  pivot_longer(cols = starts_with('mean')) 

interactions_resulting_mean_plot <- data |>
  ggplot(aes(as.numeric(time), interaction( name, fraction,  asv_num_eff), fill = value)) + 
  geom_tile()+
  geom_point(data = data |>
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), interaction( name,fraction,  asv_num_eff)),
             size = 0.75,
             alpha = 0.6, 
             shape = 4) +
  scale_fill_gradient2(low = "#8C0009", mid = "white", 
                       high = "#000000", midpoint = 0, na.value = 'white')+
  theme_bw() + 
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9),
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'))+
  labs(x = 'Time', y  = 'Blooming ASVs affectting the community',
       fill = 'Mean Interaction Strength',
       title = 'Blooming taxa effect on the community
       and the resulting abundance pos interaction') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 5),
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 4),
        axis.ticks.length = unit(0.2, "mm"))

interactions_resulting_mean_plot

ggsave(
  plot = interactions_resulting_mean_plot , 
  filename = 'bloomers_interactions_resulting_mean_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 100, 
  units = 'mm'
)

## group all those data points in which there is a blooming event vs. those that there is not 
genus_data <- tax_bbmo_10y_new |>
  dplyr::filter(asv_num %in% unique(bloo_taxonomy$asv_num_f)) |>
  dplyr::distinct(asv_num, family, genus)

data <- data |>
  dplyr::filter(!is.na(bloom_event)) |> # missing and non interpolated values
  left_join(bloo_taxonomy, by = c('asv_num_eff' = 'asv_num_f')) 

custom_labels <- data |>
  distinct(asv_num_eff, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num_eff),
      is.na(genus) ~ paste0(family_f, " ", asv_num_eff),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff) 

summary_interactions_boxplot_community_resulting_int <- data |>
  ggplot(aes(interaction( bloom_event, name, fraction), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name))+
  # geom_text(aes(label = n), 
  #           position = position_nudge(y = max(abs(data$value), na.rm = TRUE) * 0.1), 
  #           size = 1, vjust = 0.5,hjust = 1,  check_overlap = T) + 
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(asv_num_eff), labeller = as_labeller(custom_labels)#)+
             ,ncol = 2)+
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'Blooming taxa effect on the community
       and the resulting abundance post interaction') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_resulting_int

ggsave(
  plot = summary_interactions_boxplot_community_resulting_int , 
  filename = 'bloomers_summary_interactions_boxplot_community_resulting_int.pdf',
  path = 'results/figures/MDR/v2/',
  width = 100, 
  height = 150, 
  units = 'mm'
)

### Sum interactions transformed to 0 and 1 ----
mdr_tb_m_rclr <- mdr_tb_m |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(asv_num_eff %in% c('asv17', 'asv22' , 'asv23' , 'asv11', 
                                   'asv1', 'asv31', 'asv7', 'asv15')) |>
  dplyr::filter(fraction_b %in% c('bn', 'bp'))

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ 1,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_pos_int = pos_interactions/1,
                sd_pos_int = sd(mean_pos_int)) |>
  dplyr::select(-n, -pos_interactions)

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ 1,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) |>
  dplyr::mutate(mean_neg_int = neg_interactions/1,
                sd_neg_int = sd(mean_neg_int))|>
  dplyr::select(-n, -neg_interactions)

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)+1) |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |> ## add sample_number
  dplyr::select(-date)

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  full_join(bloom_event , by = c('fraction', 'asv_num_eff' = 'asv_num', 'time' = 'sample_id_num')) 

data <- data |>
  dplyr::mutate(mean_neg_int = -mean_neg_int) |>
  pivot_longer(cols = starts_with('mean')) 

interactions_resulting_mean_plot <- data |>
  ggplot(aes(as.numeric(time), interaction( name, fraction,  asv_num_eff), fill = value)) + 
  geom_tile()+
  geom_point(data = data |>
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), interaction( name,fraction,  asv_num_eff)),
             size = 0.75,
             alpha = 0.6, 
             shape = 4) +
  scale_fill_gradient2(low = "#8C0009", mid = "white", 
                       high = "#000000", midpoint = 0, na.value = 'white')+
  theme_bw() + 
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9),
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'))+
  labs(x = 'Time', y  = 'Blooming ASVs affectting the community',
       fill = 'Mean Interaction Strength',
       title = 'Blooming taxa effect on the community
       and the resulting abundance pos interaction') + 
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 5),
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 4),
        axis.ticks.length = unit(0.2, "mm"))

interactions_resulting_mean_plot

ggsave(
  plot = interactions_resulting_mean_plot , 
  filename = 'bloomers_interactions_resulting_mean_binary_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 150, 
  height = 100, 
  units = 'mm'
)

## group all those data points in which there is a blooming event vs. those that there is not 
genus_data <- tax_bbmo_10y_new |>
  dplyr::filter(asv_num %in% unique(bloo_taxonomy$asv_num_f)) |>
  dplyr::distinct(asv_num, family, genus)

data <- data |>
  dplyr::filter(!is.na(bloom_event)) |> # missing and non interpolated values
  left_join(bloo_taxonomy, by = c('asv_num_eff' = 'asv_num_f')) 

custom_labels <- data |>
  distinct(asv_num_eff, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num_eff),
      is.na(genus) ~ paste0(family_f, " ", asv_num_eff),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff) 

summary_interactions_boxplot_community_resulting_int <- data |>
  ggplot(aes(interaction( bloom_event, name, fraction), abs(value))) + 
  geom_boxplot(notch = F, aes(fill = name, color = name))+
  # geom_text(aes(label = n), 
  #           position = position_nudge(y = max(abs(data$value), na.rm = TRUE) * 0.1), 
  #           size = 1, vjust = 0.5,hjust = 1,  check_overlap = T) + 
  scale_color_manual(values = c('mean_pos_int' = '#000000',
                                'mean_neg_int' = '#8C0009'),
                     labels = c('Negative', 'Positive'))+
  scale_fill_manual(values = c('mean_pos_int' = '#000000',
                               'mean_neg_int' = '#8C0009'),
                    labels = c('Negative', 'Positive'))+
  facet_wrap(vars(asv_num_eff), labeller = as_labeller(custom_labels)#)+
             ,ncol = 2)+
  theme_bw() + 
  labs(x = '', y = 'Mean Interaction Strength', color = 'Interaction',
       fill = 'Interaction',
       title = 'Blooming taxa effect on the community
       and the resulting abundance post interaction') + 
  coord_flip()+
  theme(strip.background = element_rect(fill = 'transparent'),
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        legend.position = 'bottom',
        axis.ticks.length = unit(0.2, "mm"))

summary_interactions_boxplot_community_resulting_int

ggsave(
  plot = summary_interactions_boxplot_community_resulting_int , 
  filename = 'bloomers_summary_interactions_boxplot_community_binary_resulting_int.pdf',
  path = 'results/figures/MDR/v2/',
  width = 100, 
  height = 150, 
  units = 'mm'
)

# HOW IS THE COMMUNITY INTERACTING WITH MY APPARENT BLOOMERS OVER THE TIME SERIES? -----
## Individual plots for each ASV and fraction 
## interaction t -> t+1 and abundance at t ----
#### t+0 to sample id_num abundance 
#### sum total (+) and total (-)
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) 

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) 

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date) |>
  dplyr::filter(!(asv_num == 'asv22' & fraction == '0.2') &
                !(asv_num == 'asv31' & fraction == '0.2') &
                !(asv_num == 'asv23' & fraction == '0.2') &
                !(asv_num == 'asv11' & fraction == '3')) 

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  full_join(bloom_event , by = c('fraction', 'asv_num', 'time' = 'sample_id_num')) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  dplyr::filter(!(asv_num == 'asv23' & fraction == '0.2'))

n_bloom_events <- bloom_event |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::reframe(n = paste0('n=', n()))

data <- data |>
  pivot_longer(cols = matches('pos|neg')) |>
  dplyr::select(-n) |>
  left_join(n_bloom_events)

data$asv_num <- factor(data$asv_num,
                       levels = c('asv17','asv22',
                                  'asv1', 'asv31',
                                  'asv7', 'asv11',
                                   'asv15','asv23'))

custom_labels <- data |>
  distinct(asv_num, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num),
      is.na(genus) ~ paste0(family_f, " ", asv_num),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num)
    )
  ) |> 
  pull(label, name = asv_num)

labs_fraction <-  (c('0.2' = 'Free living (0.2-3 um)',
                     '3' = 'Particle attached (3-20 um)'))
   
# Merge both labellers into one named vector
merged_labels <- c(custom_labels, labs_fraction)

# Use as_labeller
final_labeller <- as_labeller(merged_labels)

interactions_timeseries_pre_plot <- data |> 
  ggplot(aes(as.numeric(time), as.numeric(value))) +  
  geom_hline(yintercept = 4, alpha = 0.7) + 
  geom_hline(yintercept = 0, alpha = 0.7) + 
  scale_fill_manual(values = palette_interactions) + 
  scale_color_manual(values = palette_interactions) + 
  facet_wrap(asv_num ~ fraction, ncol = 3, labeller = final_labeller) + 
  geom_area(aes(y = abundance_value * 16 + 4), fill = '#371215', alpha = 0.2) +  
  geom_line(aes(y = abundance_value * 16 + 4), color = '#371215', alpha = 0.8) + 
  geom_rect(aes(
    xmin = 0, xmax = 120, 
    ymin = 0, ymax = 4
  ), fill = "#ffffff", color = NA, alpha = 1) + 
  geom_col(aes(group = name, fill = name, color = name), alpha = 0.6) +  # Bars are drawn first
  scale_y_continuous(expand = c(0,0), 
                     sec.axis = sec_axis(~./16 , name = 'Relative Abundance', 
                                         breaks = c(0.00, 0.25, 0.50, 0.75),  
                                         labels = c('-0.25', '0.00', '0.25', '0.50'))) + 
  geom_point(data = data |> 
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), y = abundance_value * 16 + 4), 
             size = 0.75, 
             alpha = 0.6, 
             shape = 4) + 
  theme_bw() +  
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9), 
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013')) + 
  labs(x = 'Time', y  = 'Blooming ASVs being affected by the community', 
       fill = 'Sum Interaction Strength', 
       color = 'Sum Interaction Strength', 
       title = 'The microbial community effect on blooming taxa  
       and the pre abundance (t) during interaction') +  
  theme(strip.background = element_rect(fill = 'transparent'), 
        text = element_text(size = 8), 
        axis.text.x = element_text(size = 5), 
        axis.text.y = element_text(size = 6), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 2), 
        title = element_text(size = 5), 
        legend.position = 'bottom', 
        legend.key.size = unit(3, 'mm'), 
        legend.text = element_text(size = 4), 
        axis.ticks.length = unit(0.2, "mm"))

interactions_timeseries_pre_plot

ggsave(
  plot = interactions_timeseries_pre_plot , 
  filename = 'interactions_timeseries_pre_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

## interaction t -> t+1 and abundance at t+1 ----
#### t+0 to sample id_num abundance 
#### sum total (+) and total (-)
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n()) 

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(value),
                 n = n()) 

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)+1) |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date) |>
  dplyr::filter(!(asv_num == 'asv22' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '0.2') &
                  !(asv_num == 'asv23' & fraction == '0.2') &
                  !(asv_num == 'asv11' & fraction == '3')) 

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  full_join(bloom_event , by = c('fraction', 'asv_num', 'time' = 'sample_id_num')) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  dplyr::filter(!(asv_num == 'asv23' & fraction == '0.2'))

n_bloom_events <- bloom_event |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::reframe(n = paste0('n=', n()))

data <- data |>
  pivot_longer(cols = matches('pos|neg')) |>
  dplyr::select(-n) |>
  left_join(n_bloom_events)

data |>
  colnames()

data$asv_num <- factor(data$asv_num,
                       levels = c('asv17','asv22',
                                  'asv1', 'asv31',
                                  'asv7', 'asv11',
                                  'asv15','asv23'
                       ))

custom_labels <- data |>
  distinct(asv_num, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num),
      is.na(genus) ~ paste0(family_f, " ", asv_num),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num)
    )
  ) |> 
  pull(label, name = asv_num)

labs_fraction <-  (c('0.2' = 'Free living (0.2-3 um)',
                     '3' = 'Particle attached (3-20 um)'))


# Merge both labellers into one named vector
merged_labels <- c(custom_labels, labs_fraction)

# Use as_labeller
final_labeller <- as_labeller(merged_labels)

interactions_timeseries_post_plot <- data |> 
  # dplyr::filter(!is.na(value)) |>
  # dplyr::filter(!is.na(abundance_value)) |>
  ggplot(aes(as.numeric(time), as.numeric(value))) +  
  geom_hline(yintercept = 4, alpha = 0.7) + 
  geom_hline(yintercept = 0, alpha = 0.7) + 
  scale_fill_manual(values = palette_interactions) + 
  scale_color_manual(values = palette_interactions) + 
  facet_wrap(asv_num ~ fraction, ncol = 3, labeller = final_labeller) + 
  geom_area(aes(y = abundance_value * 16 + 4), fill = '#371215', alpha = 0.2) +  
  geom_line(aes(y = abundance_value * 16 + 4), color = '#371215', alpha = 0.8) + 
  geom_rect(aes(
    xmin = 0, xmax = 125, 
    ymin = 0, ymax = 4
  ), fill = "#ffffff", color = NA, alpha = 1) + 
  geom_col(aes(group = name, fill = name, color = name), alpha = 0.7) +  # Bars are drawn first
  scale_y_continuous(expand = c(0,0), 
                     sec.axis = sec_axis(~./16 , name = 'Relative Abundance', 
                                         breaks = c(0.00, 0.25, 0.50, 0.75),  
                                         labels = c('-0.25', '0.00', '0.25', '0.50'))) + 
  geom_point(data = data |> 
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), y = abundance_value * 16 + 4), 
             size = 0.75, 
             alpha = 0.6, 
             shape = 4) + 
  theme_bw() +  
  scale_x_continuous(expand = c(0,0),
    breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9), 
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013')) + 
  labs(x = 'Time', y  = 'Blooming ASVs being affected by the community', 
       fill = 'Sum Interaction Strength', 
       color = 'Sum Interaction Strength', 
       title = 'The microbial community effect on blooming taxa  
       and the post abundance (t+1) during interaction') +  
  theme(strip.background = element_rect(fill = 'transparent'), 
        text = element_text(size = 8), 
        axis.text.x = element_text(size = 5), 
        axis.text.y = element_text(size = 6), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 2), 
        title = element_text(size = 5), 
        legend.position = 'bottom', 
        legend.key.size = unit(3, 'mm'), 
        legend.text = element_text(size = 4), 
        axis.ticks.length = unit(0.2, "mm"))

interactions_timeseries_post_plot

ggsave(
  plot = interactions_timeseries_post_plot , 
  filename = 'interactions_timeseries_post_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

### BLOOMERS EFFECT ON THE COMMUNITY (abundance at (t)) ----
#### after a bloom event how do they interact with the community
### abundance at t interaction t -> t+1
mdr_tb_m_rclr <- mdr_tb_m |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  dplyr::filter(asv_num != asv_num_eff) |>
  dplyr::filter(asv_num_eff %in% c('asv17', 'asv22' , 'asv23' , 'asv11', 
                                   'asv1', 'asv31', 'asv7', 'asv15')) |>
  dplyr::filter(fraction_b %in% c('bn', 'bp'))

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::filter(value > 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value > 0 ~ value,
                                  value < 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time) |>
  dplyr::reframe(pos_interactions = sum(value),
                 n = n())

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::filter(value < 0) |>
  dplyr::mutate(value = as.numeric(value)) |>
  dplyr::mutate(value = case_when(value < 0 ~ value,
                                  value > 0 ~ 0)) |>
  dplyr::group_by(asv_num_eff, fraction_b, date, time ) |>
  dplyr::reframe(neg_interactions = sum(abs(value)),
                 n = n())

mdr_tb_m_rclr_int <- mdr_tb_m_rclr_pos |>
  left_join(mdr_tb_m_rclr_neg)

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date) |>
  dplyr::filter(!(asv_num == 'asv22' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '0.2') &
                  !(asv_num == 'asv23' & fraction == '0.2') &
                  !(asv_num == 'asv11' & fraction == '3')) 

data <- mdr_tb_m_rclr_int |> 
  dplyr::mutate(fraction = case_when(str_detect(fraction_b, 'bp') ~ '0.2',
                                     str_detect(fraction_b, 'bn') ~ '3')) |>
  full_join(bloom_event , by = c('fraction', 'asv_num_eff' = 'asv_num', 'time' = 'sample_id_num')) |>
  left_join(bloo_taxonomy, by = c('asv_num_eff' = 'asv_num_f')) |>
  dplyr::filter(!(asv_num_eff == 'asv22' & fraction == '0.2') &
                  !(asv_num_eff == 'asv31' & fraction == '0.2') &
                  !(asv_num_eff == 'asv23' & fraction == '0.2') &
                  !(asv_num_eff == 'asv11' & fraction == '3'))

data <- data |>
  dplyr::mutate(neg_interactions = -neg_interactions) |>
  pivot_longer(cols = matches('pos|neg')) |>
  dplyr::select(-n)

data |>
  colnames()

data$asv_num_eff <- factor(data$asv_num_eff,
                       levels = c('asv17','asv22',
                                  'asv1', 'asv31',
                                  'asv7', 'asv11',
                                  'asv15','asv23'
                       ))

custom_labels <- data |>
  distinct(asv_num_eff, family_f, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f," ", genus, " ", asv_num_eff),
      is.na(genus) ~ paste0(family_f, " ", asv_num_eff),  # Case with missing species
      !is.na(genus) ~ paste0(family_f, " ", genus, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff)

labs_fraction <-  (c('0.2' = 'Free living (0.2-3 um)',
                     '3' = 'Particle attached (3-20 um)'))


# Merge both labellers into one named vector
merged_labels <- c(custom_labels, labs_fraction)

# Use as_labeller
final_labeller <- as_labeller(merged_labels)

interactions_timeseries_pre_plot <- data |> 
  ggplot(aes(as.numeric(time), as.numeric(value))) +  
  geom_hline(yintercept = 4, alpha = 0.7) + 
  geom_hline(yintercept = 0, alpha = 0.7) + 
  scale_fill_manual(values = palette_interactions) + 
  scale_color_manual(values = palette_interactions) + 
  facet_wrap(asv_num_eff ~ fraction, ncol = 3, labeller = final_labeller) + 
  geom_area(aes(y = abundance_value * 16 + 4), fill = '#371215', alpha = 0.2) +  
  geom_line(aes(y = abundance_value * 16 + 4), color = '#371215', alpha = 0.8) + 
  geom_rect(aes(
    xmin = 0, xmax = 120, 
    ymin = 0, ymax = 4
  ), fill = "#ffffff", color = NA, alpha = 1) + 
  geom_col(aes(group = name, fill = name, color = name), alpha = 0.6) +  # Bars are drawn first
  scale_y_continuous(expand = c(0,0), 
                     sec.axis = sec_axis(~./16 , name = 'Relative Abundance', 
                                         breaks = c(0.00, 0.25, 0.50, 0.75),  
                                         labels = c('-0.25', '0.00', '0.25', '0.50'))) + 
  geom_point(data = data |> 
               dplyr::filter(bloom_event == 'bloom'), 
             aes(as.numeric(time), y = abundance_value * 16 + 4), 
             size = 0.75, 
             alpha = 0.6, 
             shape = 4) + 
  theme_bw() +  
  scale_x_continuous(breaks = c(12, 12*2, 12*3, 12*4, 12*5, 12*6, 12*7, 12*8, 12*9), 
                     labels = c('2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013')) + 
  labs(x = 'Time', y  = 'Blooming ASVs affecting the community', 
       fill = 'Sum Interaction Strength', 
       color = 'Sum Interaction Strength', 
       title = 'Blooming taxa effect on the community
       abundance (t) during interaction') +  
  theme(strip.background = element_rect(fill = 'transparent'), 
        text = element_text(size = 8), 
        axis.text.x = element_text(size = 5), 
        axis.text.y = element_text(size = 6), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 2), 
        title = element_text(size = 5), 
        legend.position = 'bottom', 
        legend.key.size = unit(3, 'mm'), 
        legend.text = element_text(size = 4), 
        axis.ticks.length = unit(0.2, "mm"))

interactions_timeseries_pre_plot

ggsave(
  plot = interactions_timeseries_pre_plot , 
  filename = 'bloomers_interactions_timeseries_pre_plot.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

### Community interactions at t-1 (bloom) and t ----
data_blooms <- data |>
  dplyr::filter(bloom_event == 'bloom') 

data_pre_blooms <- data_blooms |>
  dplyr::mutate(time = as.numeric(time)) |>
  dplyr::mutate(time = time-1) |>
  dplyr::filter(time > 0) |>
  dplyr::mutate(time = as.character(time)) |>
  dplyr::select(asv_num, time, fraction, abundance_value) |>
  distinct(asv_num, time, fraction) |>
  dplyr::mutate(bloom_event = 'pre-bloom') |>
  full_join(data_blooms)

data_pre_blooms <- data |>
  right_join(data_pre_blooms) |>
  dplyr::mutate(bloom_event = 'pre-bloom')

data_pre_blooms |>
  dim()

data |>
  dim()

data_blooms_preblooms <- data_blooms |>
  bind_rows(data_pre_blooms)

data_blooms_preblooms$bloom_event <- factor(data_blooms_preblooms$bloom_event,
                                            levels = c('pre-bloom', 'bloom'))

blooms_pre_blooms <- data_blooms_preblooms|>
  ggplot(aes(bloom_event, value)) +
  geom_boxplot(notch = T)+
  facet_wrap(asv_num~fraction)

blooms_pre_blooms

# ggsave(
#   plot = blooms_pre_blooms , 
#   filename = 'blooms_pre_bloom_plot.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 100, 
#   height = 120, 
#   units = 'mm'
# )

# ----- ### ----- MDR REPEAT ANALYSIS without seasonality ----- ### -----

palette_gradient <- c('#0A6260',"#BDBEBE","#545454",
                      "#F2A218")

palette_gradient <- colorRampPalette(c('#ffffff', "#BDBEBE", #"#545454",
                                       '#0A6260', "#F2A218"))(10)

## Is the community affecting differently my blooming ASVs during blooming events or not?
### 1, 0 bloom y no bloom
## heatmap + and - 
### abundance at t
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::mutate(value = if_else(as.numeric(value) > 0, 1, NA_real_)) 

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::mutate(value = if_else(as.numeric(value) < 0, 1, NA_real_))

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date) |>
  dplyr::filter(!(asv_num == 'asv22' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '0.2') &
                  !(asv_num == 'asv23' & fraction == '0.2') &
                  !(asv_num == 'asv11' & fraction == '3')) 

n_dates_fraction_asv_tb <- mdr_tb_m_rclr_pos |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  group_by(asv_num, bloom_event, fraction, date) |>
  distinct(asv_num, bloom_event, fraction, date) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::reframe(n_dates = n()) |>
  dplyr::filter(!is.na(bloom_event))

# mdr_tb_m_rclr_pos_bloo <- mdr_tb_m_rclr_pos |>
#   left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
#   dplyr::filter(!is.na(bloom_event)) |>
#   group_by(asv_num, bloom_event, fraction, asv_num_eff) |>
#   dplyr::reframe(pos_interactions = sum(as.numeric(value), na.rm = T)) |>
#   group_by(asv_num, bloom_event, fraction) |>
#   dplyr::mutate(n = n()) |>
#   group_by(asv_num, bloom_event, fraction) |>
#   dplyr::mutate(n_interactions = sum(pos_interactions)) |>
#   dplyr::mutate(rel_interactions = pos_interactions/(n_interactions)) |>
#   dplyr::mutate(type_interaction = 'positive')

mdr_tb_m_rclr_neg_bloo <- mdr_tb_m_rclr_neg |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  group_by(asv_num, bloom_event, fraction, asv_num_eff) |>
  dplyr::reframe(neg_interactions = sum(as.numeric(value), na.rm = T)) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(n = n()) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(n_interactions = sum(neg_interactions)) |>
  dplyr::mutate(rel_interactions = neg_interactions/(n_interactions)) |>
  dplyr::mutate(type_interaction = 'negative')

rel_interactions_bloo_non_bloom <- mdr_tb_m_rclr_neg_bloo |>
  bind_rows(mdr_tb_m_rclr_pos_bloo) |>
  dplyr::mutate(rel_interactions = as.numeric(rel_interactions)) |>
  left_join(tax_occ_filt_bbmo, by = c('asv_num_eff' = 'asv_num'), suffix = c('', '_tax_eff')) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'), suffix = c('', 'bloo')) |>
  rename( phylum_bloo = phylum_f, classbloo = class_f, orderbloo = order_f, familybloo = family_f, genusbloo = genus_f) |>
  dplyr::mutate(fraction_asv_num = paste0(fraction, ' ', asv_num)) |>
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event))|>
  left_join(n_dates_fraction_asv_tb)

rel_interactions_bloo_non_bloom |>
  colnames()

rel_interactions_bloo_non_bloom$asv_num <- factor(rel_interactions_bloo_non_bloom$asv_num,
                                                  levels = c('asv31', 'asv1', 'asv7',
                                                             'asv15', 'asv22', 'asv23', 'asv11', 'asv17'))

rel_interactions_bloo_non_bloom <- rel_interactions_bloo_non_bloom |>
  dplyr::mutate(phylum_f_eff = as_factor(phylum),
                family_f_eff = as_factor(family),
                order_f_eff = as_factor(order),
                class_f_eff = as_factor(class),
                asv_num_f_eff = as_factor(asv_num_eff))

rel_interactions_bloo_non_bloom$class_f_eff <-  factor(rel_interactions_bloo_non_bloom$class_f_eff, 
                                          levels=unique(rel_interactions_bloo_non_bloom$class_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff)]), 
                                          ordered=TRUE)

rel_interactions_bloo_non_bloom$order_f_eff <-  factor(rel_interactions_bloo_non_bloom$order_f_eff, 
                                          levels=unique(rel_interactions_bloo_non_bloom$order_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                             rel_interactions_bloo_non_bloom$class_f_eff)]), 
                                          ordered=TRUE)

rel_interactions_bloo_non_bloom$family_f_eff <-  factor(rel_interactions_bloo_non_bloom$family_f_eff, 
                                           levels=unique(rel_interactions_bloo_non_bloom$family_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                               rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                               rel_interactions_bloo_non_bloom$order_f_eff)]), 
                                           ordered=TRUE)

rel_interactions_bloo_non_bloom$asv_num_f_eff <-  factor(rel_interactions_bloo_non_bloom$asv_num_f_eff, 
                                            levels=unique(rel_interactions_bloo_non_bloom$asv_num_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                 rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                 rel_interactions_bloo_non_bloom$order_f_eff,
                                                                                                 rel_interactions_bloo_non_bloom$family_f_eff)]), 
                                            ordered=TRUE)

custom_labels_bloo <- rel_interactions_bloo_non_bloom |>
  distinct(asv_num, familybloo, genusbloo, fraction, n_dates) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genusbloo) ~ paste0(familybloo," ", genusbloo, " ", asv_num, " ", fraction, '(', n_dates, ')'),
      is.na(genusbloo) ~ paste0(familybloo, " ", asv_num, '(', n_dates, ')'),  # Case with missing species
      !is.na(genusbloo) ~ paste0(familybloo, " ", genusbloo, " ", asv_num," ", fraction, '(', n_dates, ')')
    )
  ) |> 
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event))|>
  pull(label, name = asv_num_fraction)
  

labs_bloom_events <- as_labeller(c('bloom' = 'Bloom', 'no-bloom' = 'No-Bloom', 'negative' = 'Negative Interactions', 
                                   'positive' = 'Positive Interactions'))

rel_interactions_bloo_non_bloom$asv_num_fraction <- factor(rel_interactions_bloo_non_bloom$asv_num_fraction,
                                                           levels = c('asv1 3 bloom', 'asv1 0.2 bloom', 'asv31 0.2 bloom',
                                                                      'asv7 3 bloom',
                                                                      'asv7 0.2 bloom','asv15 3 bloom', 'asv15 0.2 bloom',
                                                                      'asv17 3 bloom', 'asv22 3 bloom', 'asv23 3 bloom', 'asv11 0.2 bloom',
                                                                      #no-bloom
                                                                      'asv1 3 no-bloom', 'asv1 0.2 no-bloom', 'asv31 0.2 no-bloom',
                                                                      'asv7 3 no-bloom',
                                                                      'asv7 0.2 no-bloom','asv15 3 no-bloom', 'asv15 0.2 no-bloom',
                                                                      'asv17 3 no-bloom', 'asv22 3 no-bloom', 'asv23 3 no-bloom', 'asv11 0.2 no-bloom'))

## COMMUNITY INTERACTION WITH BLOOMERS -----
## edit community labels ----
custom_labels_community <- rel_interactions_bloo_non_bloom |>
  distinct(asv_num_eff, family_f_eff, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f_eff, " ", genus, " ", asv_num_eff),
      TRUE ~ paste0(family_f_eff, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff)

setdiff(unique(rel_interactions_bloo_non_bloom$asv_num_fraction), names(custom_labels_bloo))

rel_interactions_bloo_non_bloom_t <- rel_interactions_bloo_non_bloom |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) ## they only have one bloom event!

rel_interactions_mdr_bloo_nonbloo_abund_t_plot <- rel_interactions_bloo_non_bloom_t |>
  dplyr::filter(!is.na(bloom_event)) |>
  ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
             scales = 'free_x')+
  labs(fill = 'Relative Interactions', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs')+
  scale_x_discrete(label = as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.text.y = element_text( size = 2),
        axis.text.x = element_text(angle = 90, size = 2))

rel_interactions_mdr_bloo_nonbloo_abund_t_plot 

# ggsave(
#   plot = rel_interactions_mdr_bloo_nonbloo_abund_t_plot  , 
#   filename = 'rel_interactions_mdr_bloo_nonbloo_abund_t_plot.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 180, 
#   units = 'mm'
# )

### 1, 0 bloom y no bloom ----
## heatmap + and - 
### abundance at t+1
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 
                               'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::mutate(value = if_else(as.numeric(value) > 0, 1, NA_real_)) 

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::mutate(value = if_else(as.numeric(value) < 0, 1, NA_real_)) 

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)+1) |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |> ## add sample_number
  dplyr::select(-date) |>
  dplyr::filter(!(asv_num == 'asv22' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '0.2') &
                  !(asv_num == 'asv23' & fraction == '0.2') &
                  !(asv_num == 'asv11' & fraction == '3')) 

n_dates_fraction_asv_tb <- mdr_tb_m_rclr_pos |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  group_by(asv_num, bloom_event, fraction, date) |>
  distinct(asv_num, bloom_event, fraction, date) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::reframe(n_dates = n()) |>
  dplyr::filter(!is.na(bloom_event))

mdr_tb_m_rclr_pos_bloo <- mdr_tb_m_rclr_pos |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  group_by(asv_num, bloom_event, fraction, asv_num_eff) |>
  dplyr::reframe(pos_interactions = sum(as.numeric(value), na.rm = T)) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(n = n()) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(n_interactions = sum(pos_interactions)) |>
  dplyr::mutate(rel_interactions = pos_interactions/(n_interactions)) |>
  dplyr::mutate(type_interaction = 'positive')

# mdr_tb_m_rclr_pos_bloo <- mdr_tb_m_rclr_pos |>
#   left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
#   group_by(asv_num, bloom_event, fraction, asv_num_eff) |>
#   dplyr::reframe(pos_interactions = sum(as.numeric(value), na.rm = T)) |>
#   group_by(asv_num, bloom_event, fraction) |>
#   dplyr::mutate(n = n()) |>
#   left_join(n_dates_fraction_asv_tb) |>
#   dplyr::mutate(rel_interactions = pos_interactions/n) |>
#   dplyr::mutate(rel_interactions = rel_interactions/n_dates) |>
#   dplyr::mutate(type_interaction = 'positive')

mdr_tb_m_rclr_pos_bloo |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::reframe(sum_rel = sum(rel_interactions))

mdr_tb_m_rclr_neg_bloo <- mdr_tb_m_rclr_neg |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  group_by(asv_num, bloom_event, fraction, asv_num_eff) |>
  dplyr::reframe(pos_interactions = sum(as.numeric(value), na.rm = T)) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(n = n()) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(n_interactions = sum(pos_interactions)) |>
  dplyr::mutate(rel_interactions = pos_interactions/(n_interactions)) |>
  dplyr::mutate(type_interaction = 'negative')

rel_interactions_bloo_non_bloom <- mdr_tb_m_rclr_neg_bloo |>
  bind_rows(mdr_tb_m_rclr_pos_bloo) |>
  dplyr::mutate(rel_interactions = as.numeric(rel_interactions)) |>
  left_join(tax_occ_filt_bbmo, by = c('asv_num_eff' = 'asv_num'), suffix = c('', '_tax_eff')) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'), suffix = c('', 'bloo')) |>
  rename( phylum_bloo = phylum_f, classbloo = class_f, orderbloo = order_f, familybloo = family_f, genusbloo = genus_f) |>
  dplyr::mutate(fraction_asv_num = paste0(fraction, ' ', asv_num)) |>
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event)) |>
  left_join(n_dates_fraction_asv_tb) 

rel_interactions_bloo_non_bloom$asv_num <- factor(rel_interactions_bloo_non_bloom$asv_num,
                                                  levels = c('asv31', 'asv1', 'asv7',
                                                             'asv15', 'asv22', 'asv23', 'asv11', 'asv17'))

rel_interactions_bloo_non_bloom <- rel_interactions_bloo_non_bloom |>
  dplyr::mutate(phylum_f_eff = as_factor(phylum),
                family_f_eff = as_factor(family),
                order_f_eff = as_factor(order),
                genus_f_eff = as_factor(genus),
                class_f_eff = as_factor(class),
                asv_num_f_eff = as_factor(asv_num_eff))

rel_interactions_bloo_non_bloom$class_f_eff <-  factor(rel_interactions_bloo_non_bloom$class_f_eff, 
                                                       levels=unique(rel_interactions_bloo_non_bloom$class_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff)]), 
                                                       ordered=TRUE)

rel_interactions_bloo_non_bloom$order_f_eff <-  factor(rel_interactions_bloo_non_bloom$order_f_eff, 
                                                       levels=unique(rel_interactions_bloo_non_bloom$order_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                       rel_interactions_bloo_non_bloom$class_f_eff)]), 
                                                       ordered=TRUE)

rel_interactions_bloo_non_bloom$family_f_eff <-  factor(rel_interactions_bloo_non_bloom$family_f_eff, 
                                                        levels=unique(rel_interactions_bloo_non_bloom$family_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$order_f_eff)]), 
                                                        ordered=TRUE)

rel_interactions_bloo_non_bloom$genus_f_eff <-  factor(rel_interactions_bloo_non_bloom$family_f_eff, 
                                                        levels=unique(rel_interactions_bloo_non_bloom$genus_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$order_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$family_f_eff)]), 
                                                        ordered=TRUE)

rel_interactions_bloo_non_bloom$asv_num_f_eff <-  factor(rel_interactions_bloo_non_bloom$asv_num_f_eff, 
                                                         levels=unique(rel_interactions_bloo_non_bloom$asv_num_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$order_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$family_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$genus_f_eff)]), 
                                                         ordered=TRUE)

custom_labels_bloo <- rel_interactions_bloo_non_bloom |>
  distinct(asv_num, familybloo, genusbloo, fraction, n_dates, bloom_event ) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genusbloo) ~ paste0(familybloo," ", genusbloo, " ", asv_num, " ", fraction, '(', n_dates, ')'),
      is.na(genusbloo) ~ paste0(familybloo, " ", asv_num, '(', n_dates, ')'),  # Case with missing species
      !is.na(genusbloo) ~ paste0(familybloo, " ", genusbloo, " ", asv_num," ", fraction, '(', n_dates, ')')
    )
  ) |> 
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event))|>
  pull(label, name = asv_num_fraction)

custom_labels_bloo_inter <- rel_interactions_bloo_non_bloom |> 
  distinct(asv_num, familybloo, genusbloo, n_dates, bloom_event) |> # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genusbloo) ~ paste0(familybloo, " ", genusbloo, " ", asv_num),
      is.na(genusbloo) ~ paste0(familybloo, " ", asv_num)  # Case with missing genus
    )
  ) |> 
  dplyr::distinct(asv_num, label )

labs_bloom_events <- as_labeller(c('bloom' = 'Bloom', 
                                   'no-bloom' = 'No-Bloom', 
                                   'negative' = 'Negative Interactions', 
                                   'positive' = 'Positive Interactions'))

bloom_events <- tribble(
  ~asv_num, ~label,
  "bloom", "Bloom",
  "no-bloom", "No-Bloom",
  "negative", "Negative Interactions",
  "positive", "Positive Interactions"
) 

custom_labels_bloo_inter <- bloom_events |>
  bind_rows(custom_labels_bloo_inter) |>
  pull(label, name = asv_num)

## edit community labels ----
custom_labels_community <- rel_interactions_bloo_non_bloom |>
  distinct(asv_num_eff, family_f_eff, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f_eff, " ", genus, " ", asv_num_eff),
      TRUE ~ paste0(family_f_eff, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff)

setdiff(unique(rel_interactions_bloo_non_bloom$asv_num_fraction), names(custom_labels_bloo))

rel_interactions_bloo_non_bloom %$%
  range(rel_interactions)

rel_interactions_bloo_non_bloom$asv_num_fraction <- factor(rel_interactions_bloo_non_bloom$asv_num_fraction,
                                                  levels = c('asv1 3 bloom', 'asv1 0.2 bloom', 'asv31 0.2 bloom',
                                                             'asv7 3 bloom',
                                                             'asv7 0.2 bloom','asv15 3 bloom', 'asv15 0.2 bloom',
                                                             'asv17 3 bloom', 'asv22 3 bloom', 'asv23 3 bloom', 'asv11 0.2 bloom',
                                                             #no-bloom
                                                             'asv1 3 no-bloom', 'asv1 0.2 no-bloom', 'asv31 0.2 no-bloom',
                                                             'asv7 3 no-bloom',
                                                             'asv7 0.2 no-bloom','asv15 3 no-bloom', 'asv15 0.2 no-bloom',
                                                             'asv17 3 no-bloom', 'asv22 3 no-bloom', 'asv23 3 no-bloom', 'asv11 0.2 no-bloom'))

rel_interactions_mdr_bloo_nonbloo_abund_t1_plot <- rel_interactions_bloo_non_bloom |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) |> ## they only have one bloom event!
  ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events, scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs')+
  scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.text.y = element_text( size = 2),
        axis.text.x = element_text(angle = 90, size = 2))

rel_interactions_mdr_bloo_nonbloo_abund_t1_plot 

# ggsave(
#   plot = rel_interactions_mdr_bloo_nonbloo_abund_t1_plot  , 
#   filename = 'rel_interactions_mdr_bloo_nonbloo_abund_t1_plot .pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 180, 
#   units = 'mm'
# )

## combination plot of t and t+1 -----
### I need the t without axis y text 
rel_interactions_mdr_bloo_nonbloo_abund_t_plot <- rel_interactions_bloo_non_bloom_t |>
  dplyr::filter(!is.na(bloom_event)) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  facet_grid(type_interaction~asv_num, labeller =  as_labeller(custom_labels_bloo_inter),
             #labs_bloom_events,
              # as_labeller(custom_labels_bloo),
             scales = 'free_x')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2), na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = '', 
       x = 'Bloomers ASVs',
       title = 'Post Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(size = 2),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        axis.text.y = element_text( size = 0),
        axis.text.x = element_text(angle = 90, size = 2),
        axis.ticks.y = element_blank())

rel_interactions_mdr_bloo_nonbloo_abund_t_plot 

### I need t1 without legend and strip labels y 
rel_interactions_mdr_bloo_nonbloo_abund_t1_plot <- rel_interactions_bloo_non_bloom |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) |> ## they only have one bloom event!
  dplyr::filter(!is.na(bloom_event)) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  facet_grid(type_interaction~asv_num, labeller = as_labeller(custom_labels_bloo_inter),
             scales = 'free_x')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2),
  #                      na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = 'Previous Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

rel_interactions_mdr_bloo_nonbloo_abund_t1_plot 

plot_interactions_post_pre_blooms <- plot_grid(
  rel_interactions_mdr_bloo_nonbloo_abund_t1_plot,
  rel_interactions_mdr_bloo_nonbloo_abund_t_plot,
  rel_widths = c(1.2, 1))

plot_interactions_post_pre_blooms

ggsave(
  plot = plot_interactions_post_pre_blooms, 
  filename = 'plot_interactions_post_pre_blooms_v6.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

## supplementary table nº bloom events / bloomer and size fraction ----
n_bloom_events <- bloom_event |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::reframe(n = paste0('n=', n()))

# write.csv(n_bloom_events,
#           'results/tables/n_bloom_events.csv', 
#           row.names = F)

# Top 10 interactions ----
## combination plot of t and t+1 -----
### I need the t without axis y text 
rel_interactions_mdr_bloo_nonbloo_abund_t_plot <- rel_interactions_bloo_non_bloom_t |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::slice_max(n = 10, order_by = rel_interactions) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  facet_grid(type_interaction~asv_num, labeller =  as_labeller(custom_labels_bloo_inter),
             #labs_bloom_events,
             # as_labeller(custom_labels_bloo),
             scales = 'free_x')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2), na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = '', 
       x = 'Bloomers ASVs',
       title = 'Post Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(size = 2),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        axis.text.y = element_text( size = 0),
        axis.text.x = element_text(angle = 90, size = 2),
        axis.ticks.y = element_blank())

rel_interactions_mdr_bloo_nonbloo_abund_t_plot 

### I need t1 without legend and strip labels y 
rel_interactions_mdr_bloo_nonbloo_abund_t1_plot <- rel_interactions_bloo_non_bloom |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) |> ## they only have one bloom event!
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::slice_max(n = 10, order_by = rel_interactions) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  facet_grid(type_interaction~asv_num, labeller = as_labeller(custom_labels_bloo_inter),
             scales = 'free_x')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2),
  #                      na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = 'Previous Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

rel_interactions_mdr_bloo_nonbloo_abund_t1_plot 

plot_interactions_post_pre_blooms <- plot_grid(
  rel_interactions_mdr_bloo_nonbloo_abund_t1_plot,
  rel_interactions_mdr_bloo_nonbloo_abund_t_plot,
  rel_widths = c(1.2, 1))

plot_interactions_post_pre_blooms

ggsave(
  plot = plot_interactions_post_pre_blooms, 
  filename = 'plot_interactions_post_pre_blooms_top10_v7.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

## only those that are more relevant the others to supplementary -----
palette_gradient <- colorRampPalette(c("#F2A218", '#ffffff',  #"#545454",
                                       '#0A6260'))(10)

## combination plot of t and t+1 -----
rel_interactions_bloo_non_bloom_t |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::filter(!asv_num %in% c('asv11', 'asv22')) |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::slice_max(n = 10, order_by = rel_interactions) %$% 
  range( rel_interactions)

### I need the t without axis y text 
rel_interactions_mdr_bloo_nonbloo_abund_t_plot <- rel_interactions_bloo_non_bloom_t |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::filter(!asv_num %in% c('asv11', 'asv22')) |>
  # dplyr::group_by(asv_num, fraction, bloom_event) |>
  # dplyr::slice_max(n = 10, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-0.25, 0.25))+
  facet_grid(type_interaction~asv_num, labeller =  as_labeller(custom_labels_bloo_inter),
             #labs_bloom_events,
             # as_labeller(custom_labels_bloo),
             scales = 'free_x')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2), na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = '', 
       x = 'Bloomers ASVs',
       title = 'Post Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(size = 2),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        axis.text.y = element_text( size = 0),
        axis.ticks.length = unit(0.2, "mm"),
        axis.text.x = element_text(angle = 90, size = 2),
        axis.ticks.y = element_blank())

rel_interactions_mdr_bloo_nonbloo_abund_t_plot 

### I need t1 without legend and strip labels y 
rel_interactions_mdr_bloo_nonbloo_abund_t1_plot <- rel_interactions_bloo_non_bloom |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) |> ## they only have one bloom event!
  dplyr::filter(!asv_num %in% c('asv11', 'asv22')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  # dplyr::group_by(asv_num, fraction, bloom_event) |>
  # dplyr::slice_max(n = 10, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-0.25, 0.25))+
  facet_grid(type_interaction~asv_num, labeller =  as_labeller(custom_labels_bloo_inter),
             #labs_bloom_events,
             # as_labeller(custom_labels_bloo),
             scales = 'free_x')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2), na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = 'Previous Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid  = element_blank(),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        strip.text.y = element_blank(),
        axis.ticks.length = unit(0.2, "mm"),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

rel_interactions_mdr_bloo_nonbloo_abund_t1_plot 

plot_interactions_post_pre_blooms <- plot_grid(
  rel_interactions_mdr_bloo_nonbloo_abund_t1_plot,
  rel_interactions_mdr_bloo_nonbloo_abund_t_plot,
  rel_widths = c(1.2, 1))

plot_interactions_post_pre_blooms

ggsave(
  plot = plot_interactions_post_pre_blooms, 
  filename = 'plot_interactions_post_pre_blooms_top10_red_v7.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

## supplementary ---
rel_interactions_mdr_bloo_nonbloo_abund_t_supplot <- rel_interactions_bloo_non_bloom_t |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv22')) |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::slice_max(n = 10, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-0.25, 0.25))+
  facet_grid(type_interaction~asv_num, labeller =  as_labeller(custom_labels_bloo_inter),
             #labs_bloom_events,
             # as_labeller(custom_labels_bloo),
             scales = 'free_x')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2), na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = '', 
       x = 'Bloomers ASVs',
       title = 'Post Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(size = 2),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        axis.text.y = element_text( size = 0),
        axis.ticks.length = unit(0.2, "mm"),
        axis.text.x = element_text(angle = 90, size = 2),
        axis.ticks.y = element_blank())

rel_interactions_mdr_bloo_nonbloo_abund_t_supplot 

### I need t1 without legend and strip labels y 
rel_interactions_mdr_bloo_nonbloo_abund_t1_supplot <- rel_interactions_bloo_non_bloom |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) |> ## they only have one bloom event!
  dplyr::filter(asv_num %in% c('asv11', 'asv22')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event) |>
  dplyr::slice_max(n = 10, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-0.25, 0.25))+
  facet_grid(type_interaction~asv_num, labeller =  as_labeller(custom_labels_bloo_inter),
             #labs_bloom_events,
             # as_labeller(custom_labels_bloo),
             scales = 'free_x')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2), na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nFrequency', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = 'Previous Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid  = element_blank(),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        strip.text.y = element_blank(),
        axis.ticks.length = unit(0.2, "mm"),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

rel_interactions_mdr_bloo_nonbloo_abund_t1_supplot 

plot_interactions_post_pre_blooms <- plot_grid(
  rel_interactions_mdr_bloo_nonbloo_abund_t1_supplot,
  rel_interactions_mdr_bloo_nonbloo_abund_t_supplot,
  rel_widths = c(1.2, 1))

plot_interactions_post_pre_blooms

ggsave(
  plot = plot_interactions_post_pre_blooms, 
  filename = 'plot_interactions_post_pre_blooms_top10_red_supl_v1.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

### 1 when it is always interacting with my ASV and 0 never interacting with the ASV ----
### abundance at t
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::mutate(value = if_else(as.numeric(value) > 0, 1, NA_real_)) 

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::mutate(value = if_else(as.numeric(value) < 0, 1, NA_real_))

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num)

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |>## add sample_number
  dplyr::select(-date) |>
  dplyr::filter(!(asv_num == 'asv22' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '0.2') &
                  !(asv_num == 'asv23' & fraction == '0.2') &
                  !(asv_num == 'asv11' & fraction == '3')) 

n_dates_fraction_asv_tb <- mdr_tb_m_rclr_pos |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  group_by(asv_num, bloom_event, fraction, date) |>
  distinct(asv_num, bloom_event, fraction, date) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::reframe(n_dates = n()) |>
  dplyr::filter(!is.na(bloom_event))

mdr_tb_m_rclr_pos_bloo <- mdr_tb_m_rclr_pos |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  group_by(asv_num, bloom_event, fraction, asv_num_eff) |>
  dplyr::reframe(pos_interactions = sum(as.numeric(value), na.rm = T)) |>
  left_join(n_dates_fraction_asv_tb ) |>
  #group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(rel_interactions = pos_interactions/n_dates) |>
  # group_by(asv_num, bloom_event, fraction) |>
  # dplyr::mutate(n_interactions = sum(pos_interactions)) |>
  # dplyr::mutate(rel_interactions = pos_interactions/(n_interactions)) |>
  dplyr::mutate(type_interaction = 'positive')

mdr_tb_m_rclr_neg_bloo <- mdr_tb_m_rclr_neg |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  group_by(asv_num, bloom_event, fraction, asv_num_eff) |>
  dplyr::reframe(neg_interactions = sum(as.numeric(value), na.rm = T)) |>
  left_join(n_dates_fraction_asv_tb ) |>
  #group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(rel_interactions = neg_interactions/n_dates) |>
  # group_by(asv_num, bloom_event, fraction) |>
  # dplyr::mutate(n_interactions = sum(pos_interactions)) |>
  # dplyr::mutate(rel_interactions = pos_interactions/(n_interactions)) |>
  dplyr::mutate(type_interaction = 'negative')

rel_interactions_bloo_non_bloom <- mdr_tb_m_rclr_neg_bloo |>
  bind_rows(mdr_tb_m_rclr_pos_bloo) |>
  dplyr::mutate(rel_interactions = as.numeric(rel_interactions)) |>
  left_join(tax_occ_filt_bbmo, by = c('asv_num_eff' = 'asv_num'), suffix = c('', '_tax_eff')) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'), suffix = c('', 'bloo')) |>
  rename( phylum_bloo = phylum_f, classbloo = class_f, orderbloo = order_f, familybloo = family_f, genusbloo = genus_f) |>
  dplyr::mutate(fraction_asv_num = paste0(fraction, ' ', asv_num)) |>
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event))|>
  left_join(n_dates_fraction_asv_tb)

rel_interactions_bloo_non_bloom |>
  colnames()

rel_interactions_bloo_non_bloom$asv_num <- factor(rel_interactions_bloo_non_bloom$asv_num,
                                                  levels = c('asv31', 'asv1', 'asv7',
                                                             'asv15', 'asv22', 'asv23', 'asv11', 'asv17'))

rel_interactions_bloo_non_bloom <- rel_interactions_bloo_non_bloom |>
  dplyr::mutate(phylum_f_eff = as_factor(phylum),
                family_f_eff = as_factor(family),
                order_f_eff = as_factor(order),
                class_f_eff = as_factor(class),
                asv_num_f_eff = as_factor(asv_num_eff))

rel_interactions_bloo_non_bloom$class_f_eff <-  factor(rel_interactions_bloo_non_bloom$class_f_eff, 
                                                       levels=unique(rel_interactions_bloo_non_bloom$class_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff)]), 
                                                       ordered=TRUE)

rel_interactions_bloo_non_bloom$order_f_eff <-  factor(rel_interactions_bloo_non_bloom$order_f_eff, 
                                                       levels=unique(rel_interactions_bloo_non_bloom$order_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                       rel_interactions_bloo_non_bloom$class_f_eff)]), 
                                                       ordered=TRUE)

rel_interactions_bloo_non_bloom$family_f_eff <-  factor(rel_interactions_bloo_non_bloom$family_f_eff, 
                                                        levels=unique(rel_interactions_bloo_non_bloom$family_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$order_f_eff)]), 
                                                        ordered=TRUE)

rel_interactions_bloo_non_bloom$asv_num_f_eff <-  factor(rel_interactions_bloo_non_bloom$asv_num_f_eff, 
                                                         levels=unique(rel_interactions_bloo_non_bloom$asv_num_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$order_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$family_f_eff)]), 
                                                         ordered=TRUE)

custom_labels_bloo <- rel_interactions_bloo_non_bloom |>
  distinct(asv_num, familybloo, genusbloo, fraction, n_dates) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genusbloo) ~ paste0(familybloo," ", genusbloo, " ", asv_num, " ", fraction, '(', n_dates, ')'),
      is.na(genusbloo) ~ paste0(familybloo, " ", asv_num, '(', n_dates, ')'),  # Case with missing species
      !is.na(genusbloo) ~ paste0(familybloo, " ", genusbloo, " ", asv_num," ", fraction, '(', n_dates, ')')
    )
  ) |> 
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event))|>
  pull(label, name = asv_num_fraction)


labs_bloom_events <- as_labeller(c('bloom' = 'Bloom', 'no-bloom' = 'No-Bloom', 'negative' = 'Negative Interactions', 
                                   'positive' = 'Positive Interactions'))

rel_interactions_bloo_non_bloom$asv_num_fraction <- factor(rel_interactions_bloo_non_bloom$asv_num_fraction,
                                                           levels = c('asv1 3 bloom', 'asv1 0.2 bloom', 'asv31 0.2 bloom',
                                                                      'asv7 3 bloom',
                                                                      'asv7 0.2 bloom','asv15 3 bloom', 'asv15 0.2 bloom',
                                                                      'asv17 3 bloom', 'asv22 3 bloom', 'asv23 3 bloom', 'asv11 0.2 bloom',
                                                                      #no-bloom
                                                                      'asv1 3 no-bloom', 'asv1 0.2 no-bloom', 'asv31 0.2 no-bloom',
                                                                      'asv7 3 no-bloom',
                                                                      'asv7 0.2 no-bloom','asv15 3 no-bloom', 'asv15 0.2 no-bloom',
                                                                      'asv17 3 no-bloom', 'asv22 3 no-bloom', 'asv23 3 no-bloom', 'asv11 0.2 no-bloom'))

## COMMUNITY INTERACTION WITH BLOOMERS -----
## edit community labels ----
custom_labels_community <- rel_interactions_bloo_non_bloom |>
  distinct(asv_num_eff, family_f_eff, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f_eff, " ", genus, " ", asv_num_eff),
      TRUE ~ paste0(family_f_eff, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff)

setdiff(unique(rel_interactions_bloo_non_bloom$asv_num_fraction), names(custom_labels_bloo))

rel_interactions_bloo_non_bloom_t <- rel_interactions_bloo_non_bloom |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) ## they only have one bloom event!

### 1, 0 bloom y no bloom ----
## heatmap + and - 
### abundance at t+1
mdr_tb_m_rclr <- mdr_tb_m |>
  dplyr::filter(asv_num %in% c('asv17', 'asv22', 'asv23',
                               'asv11', 'asv1', 'asv31', 
                               'asv7', 'asv15')) |>
  pivot_longer(cols = starts_with(c('bp', 'bn'))) |>
  separate(name, sep = '_', into = c('fraction_eff', 'asv_num_eff')) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num_eff' = 'asv_num')) |>
  dplyr::filter(asv_num != asv_num_eff)

mdr_tb_m_rclr_pos <- mdr_tb_m_rclr |>
  dplyr::mutate(value = if_else(as.numeric(value) > 0, 1, NA_real_)) 

mdr_tb_m_rclr_neg <- mdr_tb_m_rclr |>
  dplyr::mutate(value = if_else(as.numeric(value) < 0, 1, NA_real_)) 

m_sample_num <- m_02 |>
  dplyr::select(date, sample_id_num) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)+1) |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

bloom_event <- asv_tab_all_bloo_z_tax_red |>
  dplyr::mutate(bloom_event = case_when((z_score_ra > 2 &
                                           abundance_value > 0.1) ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(bloom_event, date, fraction, asv_num, abundance_value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_sample_num) |> ## add sample_number
  dplyr::select(-date) |>
  dplyr::filter(!(asv_num == 'asv22' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '0.2') &
                  !(asv_num == 'asv23' & fraction == '0.2') &
                  !(asv_num == 'asv11' & fraction == '3')) 

n_dates_fraction_asv_tb <- mdr_tb_m_rclr_pos |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  group_by(asv_num, bloom_event, fraction, date) |>
  distinct(asv_num, bloom_event, fraction, date) |>
  group_by(asv_num, bloom_event, fraction) |>
  dplyr::reframe(n_dates = n()) |>
  dplyr::filter(!is.na(bloom_event))

mdr_tb_m_rclr_pos_bloo <- mdr_tb_m_rclr_pos |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  group_by(asv_num, bloom_event, fraction, asv_num_eff, fraction_eff) |>
  dplyr::reframe(pos_interactions = sum(as.numeric(value), na.rm = T)) |>
  left_join(n_dates_fraction_asv_tb ) |>
  #group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(rel_interactions = pos_interactions/n_dates) |>
  # group_by(asv_num, bloom_event, fraction) |>
  # dplyr::mutate(n_interactions = sum(pos_interactions)) |>
  # dplyr::mutate(rel_interactions = pos_interactions/(n_interactions)) |>
  dplyr::mutate(type_interaction = 'positive')

mdr_tb_m_rclr_neg_bloo <- mdr_tb_m_rclr_neg |>
  left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  group_by(asv_num, bloom_event, fraction, asv_num_eff, fraction_eff) |>
  dplyr::reframe(neg_interactions = sum(as.numeric(value), na.rm = T)) |>
  left_join(n_dates_fraction_asv_tb, by = c('asv_num', 'bloom_event', 'fraction') ) |>
  #group_by(asv_num, bloom_event, fraction) |>
  dplyr::mutate(rel_interactions = neg_interactions/n_dates) |>
  # group_by(asv_num, bloom_event, fraction) |>
  # dplyr::mutate(n_interactions = sum(pos_interactions)) |>
  # dplyr::mutate(rel_interactions = pos_interactions/(n_interactions)) |>
  dplyr::mutate(type_interaction = 'negative')

# mdr_tb_m_rclr_neg_bloo <- mdr_tb_m_rclr_neg |>
#   left_join(bloom_event, by = c('asv_num', 'time' = 'sample_id_num', 'fraction')) |>
#   dplyr::filter(!is.na(bloom_event)) |>
#   group_by(asv_num, bloom_event, fraction, asv_num_eff) |>
#   dplyr::reframe(pos_interactions = sum(as.numeric(value), na.rm = T)) |>
#   group_by(asv_num, bloom_event, fraction) |>
#   dplyr::mutate(n = n()) |>
#   group_by(asv_num, bloom_event, fraction) |>
#   dplyr::mutate(n_interactions = sum(pos_interactions)) |>
#   dplyr::mutate(rel_interactions = pos_interactions/(n_interactions)) |>
#   dplyr::mutate(type_interaction = 'negative')

rel_interactions_bloo_non_bloom <- mdr_tb_m_rclr_neg_bloo |>
  bind_rows(mdr_tb_m_rclr_pos_bloo) |>
  dplyr::mutate(rel_interactions = as.numeric(rel_interactions)) |>
  left_join(tax_occ_filt_bbmo, by = c('asv_num_eff' = 'asv_num'), suffix = c('', '_tax_eff')) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'), suffix = c('', 'bloo')) |>
  rename( phylum_bloo = phylum_f, classbloo = class_f, orderbloo = order_f, familybloo = family_f, genusbloo = genus_f) |>
  dplyr::mutate(fraction_asv_num = paste0(fraction, ' ', asv_num)) |>
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event)) |>
  left_join(n_dates_fraction_asv_tb) 

rel_interactions_bloo_non_bloom$asv_num <- factor(rel_interactions_bloo_non_bloom$asv_num,
                                                  levels = c('asv31', 'asv1', 'asv7',
                                                             'asv15', 'asv22', 'asv23', 'asv11', 'asv17'))

rel_interactions_bloo_non_bloom <- rel_interactions_bloo_non_bloom |>
  dplyr::mutate(phylum_f_eff = as_factor(phylum),
                family_f_eff = as_factor(family),
                order_f_eff = as_factor(order),
                genus_f_eff = as_factor(genus),
                class_f_eff = as_factor(class),
                asv_num_f_eff = as_factor(asv_num_eff))

rel_interactions_bloo_non_bloom$class_f_eff <-  factor(rel_interactions_bloo_non_bloom$class_f_eff, 
                                                       levels=unique(rel_interactions_bloo_non_bloom$class_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff)]), 
                                                       ordered=TRUE)

rel_interactions_bloo_non_bloom$order_f_eff <-  factor(rel_interactions_bloo_non_bloom$order_f_eff, 
                                                       levels=unique(rel_interactions_bloo_non_bloom$order_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                       rel_interactions_bloo_non_bloom$class_f_eff)]), 
                                                       ordered=TRUE)

rel_interactions_bloo_non_bloom$family_f_eff <-  factor(rel_interactions_bloo_non_bloom$family_f_eff, 
                                                        levels=unique(rel_interactions_bloo_non_bloom$family_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                                         rel_interactions_bloo_non_bloom$order_f_eff)]), 
                                                        ordered=TRUE)

rel_interactions_bloo_non_bloom$genus_f_eff <-  factor(rel_interactions_bloo_non_bloom$family_f_eff, 
                                                       levels=unique(rel_interactions_bloo_non_bloom$genus_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                       rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                                       rel_interactions_bloo_non_bloom$order_f_eff,
                                                                                                                       rel_interactions_bloo_non_bloom$family_f_eff)]), 
                                                       ordered=TRUE)

rel_interactions_bloo_non_bloom$asv_num_f_eff <-  factor(rel_interactions_bloo_non_bloom$asv_num_f_eff, 
                                                         levels=unique(rel_interactions_bloo_non_bloom$asv_num_f_eff[order(rel_interactions_bloo_non_bloom$phylum_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$class_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$order_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$family_f_eff,
                                                                                                                           rel_interactions_bloo_non_bloom$genus_f_eff)]), 
                                                         ordered=TRUE)

custom_labels_bloo <- rel_interactions_bloo_non_bloom |>
  distinct(asv_num, familybloo, genusbloo, fraction, n_dates, bloom_event ) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genusbloo) ~ paste0(familybloo," ", genusbloo, " ", asv_num, " ", fraction, '(', n_dates, ')'),
      is.na(genusbloo) ~ paste0(familybloo, " ", asv_num, '(', n_dates, ')'),  # Case with missing species
      !is.na(genusbloo) ~ paste0(familybloo, " ", genusbloo, " ", asv_num," ", fraction, '(', n_dates, ')')
    )
  ) |> 
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event))|>
  pull(label, name = asv_num_fraction)

custom_labels_bloo_inter <- rel_interactions_bloo_non_bloom |> 
  distinct(asv_num, familybloo, genusbloo, n_dates, bloom_event) |> # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genusbloo) ~ paste0(familybloo, " ", genusbloo, " ", asv_num),
      is.na(genusbloo) ~ paste0(familybloo, " ", asv_num)  # Case with missing genus
    )
  ) |> 
  dplyr::distinct(asv_num, label )

labs_bloom_events <- as_labeller(c('bloom' = 'Bloom', 
                                   'no-bloom' = 'No-Bloom', 
                                   'negative' = 'Negative Interactions', 
                                   'positive' = 'Positive Interactions'))

bloom_events <- tribble(
  ~asv_num, ~label,
  "bloom", "Bloom",
  "no-bloom", "No-Bloom",
  "negative", "Negative Interactions",
  "positive", "Positive Interactions"
) 

custom_labels_bloo_inter <- bloom_events |>
  bind_rows(custom_labels_bloo_inter) |>
  pull(label, name = asv_num)

## edit community labels ----
custom_labels_community <- rel_interactions_bloo_non_bloom |>
  distinct(asv_num_eff, family_f_eff, genus) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = case_when(
      !is.na(genus) ~ paste0(family_f_eff, " ", genus, " ", asv_num_eff),
      TRUE ~ paste0(family_f_eff, " ", asv_num_eff)
    )
  ) |> 
  pull(label, name = asv_num_eff)

setdiff(unique(rel_interactions_bloo_non_bloom$asv_num_fraction), names(custom_labels_bloo))

rel_interactions_bloo_non_bloom %$%
  range(rel_interactions)

rel_interactions_bloo_non_bloom$asv_num_fraction <- factor(rel_interactions_bloo_non_bloom$asv_num_fraction,
                                                           levels = c('asv1 3 bloom', 'asv1 0.2 bloom', 'asv31 0.2 bloom',
                                                                      'asv7 3 bloom',
                                                                      'asv7 0.2 bloom','asv15 3 bloom', 'asv15 0.2 bloom',
                                                                      'asv17 3 bloom', 'asv22 3 bloom', 'asv23 3 bloom', 'asv11 0.2 bloom',
                                                                      #no-bloom
                                                                      'asv1 3 no-bloom', 'asv1 0.2 no-bloom', 'asv31 0.2 no-bloom',
                                                                      'asv7 3 no-bloom',
                                                                      'asv7 0.2 no-bloom','asv15 3 no-bloom', 'asv15 0.2 no-bloom',
                                                                      'asv17 3 no-bloom', 'asv22 3 no-bloom', 'asv23 3 no-bloom', 'asv11 0.2 no-bloom'))

rel_interactions_bloo_non_bloom$asv_num_fraction <- factor(rel_interactions_bloo_non_bloom$asv_num_fraction,
                                                           levels = c( 'asv1 3 no-bloom', 'asv1 3 bloom',  'asv1 0.2 no-bloom', 'asv1 0.2 bloom', 
                                                                       'asv31 0.2 no-bloom', 'asv31 0.2 bloom',
                                                                       'asv7 3 no-bloom', 'asv7 3 bloom',
                                                                       'asv7 0.2 no-bloom', 'asv7 0.2 bloom',
                                                                       'asv15 3 no-bloom', 'asv15 3 bloom',  'asv15 0.2 no-bloom', 'asv15 0.2 bloom',
                                                                       'asv17 3 no-bloom', 'asv17 3 bloom', 'asv22 3 no-bloom','asv22 3 bloom', 
                                                                       'asv23 3 no-bloom', 'asv23 3 bloom',   'asv11 0.2 no-bloom', 'asv11 0.2 bloom'#no-bloom,
                                                                      ))

rel_interactions_bloo_non_bloom_t$asv_num_fraction <- factor(rel_interactions_bloo_non_bloom_t$asv_num_fraction,
                                                           levels = c( 'asv1 3 no-bloom', 'asv1 3 bloom',  'asv1 0.2 no-bloom', 'asv1 0.2 bloom', 
                                                                       'asv31 0.2 no-bloom', 'asv31 0.2 bloom',
                                                                       'asv7 3 no-bloom', 'asv7 3 bloom',
                                                                       'asv7 0.2 no-bloom', 'asv7 0.2 bloom',
                                                                       'asv15 3 no-bloom', 'asv15 3 bloom',  'asv15 0.2 no-bloom', 'asv15 0.2 bloom',
                                                                       'asv17 3 no-bloom', 'asv17 3 bloom', 'asv22 3 no-bloom','asv22 3 bloom', 
                                                                       'asv23 3 no-bloom', 'asv23 3 bloom',   'asv11 0.2 no-bloom', 'asv11 0.2 bloom'#no-bloom,
                                                           ))

## combination plot of t and t+1 -----
### I need the t without axis y text 
rel_interactions_mdr_bloo_nonbloo_abund_t_plot <- rel_interactions_bloo_non_bloom_t |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  #scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  facet_grid(.~asv_num, labeller =  as_labeller(custom_labels_bloo_inter),
             #labs_bloom_events,
             # as_labeller(custom_labels_bloo),
             scales = 'free')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2), na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = '', 
       x = 'Bloomers ASVs',
       title = 'Post Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(size = 2),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        axis.text.y = element_text( size = 0),
        axis.ticks.length = unit(0.2, "mm"),
        axis.text.x = element_text(angle = 90, size = 2),
        axis.ticks.y = element_blank())

rel_interactions_mdr_bloo_nonbloo_abund_t_plot 

### I need t1 without legend and strip labels y 
rel_interactions_mdr_bloo_nonbloo_abund_t1_plot <- rel_interactions_bloo_non_bloom |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) |> ## they only have one bloom event!
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(interaction(bloom_event, asv_num_fraction), 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  #scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  #geom_tile()+
  #scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'))+
  facet_grid(.~asv_num, labeller = as_labeller(custom_labels_bloo_inter),
             scales = 'free')+
  # ggplot(aes(asv_num_fraction, asv_num_f_eff, fill = rel_interactions))+
  # geom_tile()+
  # scale_fill_gradientn(colors = palette_gradient, limits = c(0,0.2),
  #                      na.value = '#ffffff')+
  # facet_grid(type_interaction~bloom_event, labeller = labs_bloom_events,
  #            scales = 'free_x')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = 'Previous Interactions')+
  #scale_x_discrete(label =   as_labeller(custom_labels_bloo))+
  scale_y_discrete(labels = as_labeller(custom_labels_community))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

rel_interactions_mdr_bloo_nonbloo_abund_t1_plot 

plot_interactions_post_pre_blooms <- plot_grid(
  rel_interactions_mdr_bloo_nonbloo_abund_t1_plot,
  rel_interactions_mdr_bloo_nonbloo_abund_t_plot,
  rel_widths = c(1.2, 1))

plot_interactions_post_pre_blooms

ggsave(
  plot = plot_interactions_post_pre_blooms, 
  filename = 'plot_interactions_post_pre_blooms_relativeinteractions_top5_events_v2.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 200, 
  units = 'mm'
)

## Plot individual por cada ASV ----
rel_interactions_bloo_non_bloom  |>
  colnames()

rel_interactions_bloo_non_bloom_t <- rel_interactions_bloo_non_bloom_t |>
  dplyr::mutate(interaction_moment = 't')

rel_interactions_bloo_non_bloom_all <- rel_interactions_bloo_non_bloom |>
  dplyr::mutate(interaction_moment = 't_1') |>
  bind_rows(rel_interactions_bloo_non_bloom_t)

rel_interactions_bloo_non_bloom_all$interaction_moment <- factor(rel_interactions_bloo_non_bloom_all$interaction_moment,
                                                                 levels = c('t_1', 't'), labels = c('Community Interactions t-1',
                                                                                                    'Community Interactions t'))

rel_interactions_bloo_non_bloom_all <- rel_interactions_bloo_non_bloom_all |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) ## they only have one bloom event!

custom_labels_bloo_v2 <- rel_interactions_bloo_non_bloom_all |>
  distinct(asv_num, familybloo, genusbloo, fraction, n_dates, bloom_event ) |>    # Ensure one label per `asv_num`
  dplyr::mutate(
    label = paste0(asv_num," ", fraction, bloom_event, '(', n_dates, ')')
  ) |> 
  dplyr::mutate(asv_num_fraction = paste0(asv_num, ' ', fraction, ' ', bloom_event))|>
  pull(label, name = asv_num_fraction)

title_text <- custom_labels_bloo_inter |>
  as_tibble()

title_text <- title_text |>
  dplyr::filter(str_detect(value, 'asv7'))

asv7_interactions_rel_plot <- rel_interactions_bloo_non_bloom_all |>
  dplyr::filter(asv_num == 'asv7') |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction, interaction_moment) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(asv_num_fraction, 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  facet_grid(.~interaction_moment, labeller = ,
             scales = 'free')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = title_text$value)+
  scale_y_discrete(labels = custom_labels_community)+
  scale_x_discrete(label = custom_labels_bloo_v2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

asv7_interactions_rel_plot

title_text <- custom_labels_bloo_inter |>
  as_tibble()

title_text <- title_text |>
  dplyr::filter(str_detect(value, 'asv1'))

asv1_interactions_rel_plot <- rel_interactions_bloo_non_bloom_all |>
  dplyr::filter(asv_num == 'asv1') |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2') &
                  !(asv_num == 'asv31' & fraction == '3')) |> ## they only have one bloom event!
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction, interaction_moment) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(asv_num_fraction, 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  facet_grid(.~interaction_moment, labeller = ,
             scales = 'free')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = title_text$value)+
  scale_y_discrete(labels = custom_labels_community)+
  scale_x_discrete(label = custom_labels_bloo_v2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

asv1_interactions_rel_plot

title_text <- custom_labels_bloo_inter |>
  as_tibble()

title_text <- title_text |>
  dplyr::filter(str_detect(value, 'asv15'))

asv15_interactions_rel_plot <- rel_interactions_bloo_non_bloom_all |>
  dplyr::filter(asv_num == 'asv15') |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction, interaction_moment) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(asv_num_fraction, 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  facet_grid(.~interaction_moment, labeller = ,
             scales = 'free')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = title_text$value)+
  scale_y_discrete(labels = custom_labels_community)+
  scale_x_discrete(label = custom_labels_bloo_v2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

asv15_interactions_rel_plot

title_text <- custom_labels_bloo_inter |>
  as_tibble()

title_text <- title_text |>
  dplyr::filter(str_detect(value, 'asv17'))

asv17_interactions_rel_plot <- rel_interactions_bloo_non_bloom_all |>
  dplyr::filter(asv_num == 'asv17') |>
  dplyr::filter(!(asv_num == 'asv17' & fraction == '0.2')) |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction, interaction_moment) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(asv_num_fraction, 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  facet_grid(.~interaction_moment, labeller = ,
             scales = 'free')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = title_text$value)+
  scale_y_discrete(labels = custom_labels_community)+
  scale_x_discrete(label = custom_labels_bloo_v2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

asv17_interactions_rel_plot

title_text <- custom_labels_bloo_inter |>
  as_tibble()

title_text <- title_text |>
  dplyr::filter(str_detect(value, 'asv22'))

asv22_interactions_rel_plot <- rel_interactions_bloo_non_bloom_all |>
  dplyr::filter(asv_num == 'asv22') |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction, interaction_moment) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(asv_num_fraction, 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  facet_grid(.~interaction_moment, labeller = ,
             scales = 'free')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = title_text$value)+
  scale_y_discrete(labels = custom_labels_community)+
  scale_x_discrete(label = custom_labels_bloo_v2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

asv22_interactions_rel_plot

title_text <- custom_labels_bloo_inter |>
  as_tibble()

title_text <- title_text |>
  dplyr::filter(str_detect(value, 'asv23'))

asv23_interactions_rel_plot <- rel_interactions_bloo_non_bloom_all |>
  dplyr::filter(asv_num == 'asv23') |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction, interaction_moment) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(asv_num_fraction, 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  facet_grid(.~interaction_moment, labeller = ,
             scales = 'free')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = title_text$value)+
  scale_y_discrete(labels = custom_labels_community)+
  scale_x_discrete(label = custom_labels_bloo_v2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'none')

asv23_interactions_rel_plot

title_text <- custom_labels_bloo_inter |>
  as_tibble()

title_text <- title_text |>
  dplyr::filter(str_detect(value, 'asv11'))

asv11_interactions_rel_plot <- rel_interactions_bloo_non_bloom_all |>
  dplyr::filter(asv_num == 'asv11') |>
  dplyr::filter(!is.na(bloom_event)) |>
  dplyr::group_by(asv_num, fraction, bloom_event, type_interaction, interaction_moment) |>
  dplyr::slice_max(n = 5, order_by = rel_interactions) |>
  dplyr::mutate(rel_interactions = case_when(type_interaction == 'negative' ~ -rel_interactions,
                                             TRUE ~ rel_interactions)) |>
  dplyr::mutate(rel_interactions = case_when(rel_interactions == '0' ~ NA_real_,
                                             TRUE ~ rel_interactions)) |>
  ggplot(aes(asv_num_fraction, 
             asv_num_f_eff, fill = rel_interactions))+
  geom_tile()+
  scale_fill_gradientn(colors = palette_gradient, na.value = ('#ffffff'),
                       limits = c(-1, 1))+
  facet_grid(.~interaction_moment, labeller = ,
             scales = 'free')+
  labs(fill = 'Relative\nInteraction\nDates', 
       y = 'ASVs interacting with bloomers', 
       x = 'Bloomers ASVs',
       title = title_text$value)+
  scale_y_discrete(labels = custom_labels_community)+
  scale_x_discrete(label = custom_labels_bloo_v2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        title = element_text(size = 6),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2),
        legend.position = 'bottom')

asv11_interactions_rel_plot

plot <- rel_interactions_bloo_non_bloom_all |>
  ggplot(aes(asv_num_fraction, fill = rel_interactions))+
  scale_fill_gradientn(colors = palette_gradient,
                       limits = c(-1, 1))+
  scale_color_manual(guide = 'none')+
  theme(
        panel.grid = element_blank(),
        title = element_text(size = 6),
        strip.text = element_text(size = 2),
        axis.text.y = element_text( size = 2),
        axis.ticks.length = unit(0.2, "mm"),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 2))+
  guides(fill = guide_legend(title = "Relative\nInteraction\nDates"), 
         color = "none", 
         shape = 'none')

legend_interactions <- cowplot::get_legend(plot)

interactions_rel_pots <- plot_grid(asv1_interactions_rel_plot, 
          asv7_interactions_rel_plot,
          asv15_interactions_rel_plot,
          asv17_interactions_rel_plot,
          asv22_interactions_rel_plot,
          asv23_interactions_rel_plot,
          asv11_interactions_rel_plot,
          rel_heights = c(1,1,1.25),
          labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'),
          ncol = 3)

interactions_rel_pots

ggsave(
  plot = interactions_rel_pots, 
  filename = 'interactions_rel_pots.pdf',
  path = 'results/figures/MDR/v2/',
  width = 180, 
  height = 180, 
  units = 'mm'
)

## BLOOMERS INTERACTION WITH THE COMMUNITY (TO BE PREPARED!) -----

