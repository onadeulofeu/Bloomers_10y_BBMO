# EDM empircal dynamic modeling
## These methods can be used to provide a mechanistic understanding of dynamical systems. 
## These non-linear statistical methods are rooted in state space reconstruction (SSR), 
## i.e. lagged coordinate embedding of time series data
## These methods do not assume any set of equations governing the system but recover the 
## dynamics from time series data, thus called empirical dynamic modeling. 
## Missing data impart an unavoidably negative influence on the performance of EDM.

library(rEDM)

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

### We need occurrence and relative abundance - rank table-----

asv_tab_all_bloo_z_tax <- read_csv2('asv_tab_all_bloo_z_tax.csv') |>
  as_tibble()

asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  group_by(asv_num) |>
  distinct(asv_num) #58 potential bloomers in my dataset

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
  dplyr::filter(occurrence_perc > 2/3) 

occurrence_bloo_bbmo_summary |>
  group_by(asv_num) |>
  distinct(asv_num) #13 ASVs with an occurrence >66% in the dataset.

write.csv(occurrence_bloo_bbmo, 'occurrence_bloo_bbmo.csv')

###plots -----
library(scales)
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
  ggplot(aes(abundance_value, occurrence_perc, color = family))+
  geom_point(aes(size = abundance_value, color = family, alpha = 0.8))+
  scale_y_continuous(labels = percent_format())+
  #facet_wrap(vars(phylum))+
  #scale_x_datetime()+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()



##upload occurrence data -----
occurrence_bloo_bbmo <- read.delim2('occurrence_bloo_bbmo.csv', sep = ',')
occurrence_bloo_bbmo |>
  head()

## Preparation of a ASV_tab with all ASVs of the temporal series that are present in 2/3 of the dataset and a tax table for them with bloomers identified------
setwd("~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/")
bbmo_10y <-readRDS("blphy10years.rds") ##8052 asv amb totes les mostres, no está en % aquesta taula
#str(bbmo_10y)
#bbmo_10y <- prune_samples(sample_sums(bbmo_10y) > 10000, bbmo_10y) in case we want to filter samples by number of reads

bbmo_10y <-
  prune_taxa(taxa_sums(bbmo_10y@otu_table) >0, ##filtrar per les asv que son 0 tot el dataset
             bbmo_10y)

bbmo_10y |>
  nsamples() #237 samples 

bbmo_10y |>
  ntaxa() #8052 but if we filter for those that are 0 during the whole dataset then we got 7,849 ASVs

## separate datasets by ASV_tab, taxonomy and metadata
asv_tab_bbmo_10y_l <- bbmo_10y@otu_table |>
  as_tibble()

m_bbmo_10y <- bbmo_10y@sam_data |>
  as_tibble()

#new taxonomy created with the database SILVA 138
new_tax <-  readRDS('03_tax_assignation/devotes_all_tax_assignation.rds') |>
  as_tibble(rownames = 'sequence')

tax_bbmo_10y_new <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax, by = c('seq' = 'sequence'))


## tidy colnames----
asv_tab_bbmo_10y_l |>
  colnames()

tax_bbmo_10y_old |>
  colnames()

colnames(asv_tab_bbmo_10y_l) <- c('asv_num', "sample_id", 'reads')

colnames(tax_bbmo_10y_old) <- c("asv_num", "kingdom", "phylum", "class", "order", "family", "genus",
                                "species", "curated", "otu_corr","seq")

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
                          "bacteria_joint", "synechococcus", "depth", "name_complete")


## Divide metadata into FL and PA----
m_02 <- m_bbmo_10y  |>
  dplyr::filter(fraction == 0.2) 

m_02 <- m_02 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02)))

m_3 <- m_bbmo_10y |>
  dplyr::filter(fraction == 3.0)

m_3 <- m_3 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3)))

##we lack of 3 samples in FL fraction which are:----
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

#check even sequencing----
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

# Calculate relative abundance----
asv_tab_10y_l_rel <- asv_tab_bbmo_10y_l |>
  calculate_rel_abund(group_cols = sample_id)

asv_tab_10y_3_rel <- asv_tab_10y_l_rel %>%
  dplyr::filter(sample_id %in% m_3$sample_id)

asv_tab_10y_02_rel <- asv_tab_10y_l_rel %>%
  dplyr::filter(sample_id %in% m_02$sample_id)


#Calculate occurrence----
## Occurence by fraction (0.2 - 3)
asv_tab_10y_l_rel |> 
  colnames()
m_bbmo_10y_sim <- m_bbmo_10y |>
  dplyr::select(fraction, sample_id, Year, Month, Day)

asv_tab_10y_l_rel_occ <- asv_tab_10y_l_rel |>
  dplyr::left_join(m_bbmo_10y_sim) |>
  # dplyr::mutate(fraction = case_when(str_detect(sample_id, '0.2') ~ 0.2,
  #                                    str_detect(sample_id, '_3_') ~ 3)) |>
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
 # group_by(asv_num, occurrence_perc) |>
  dplyr::filter(occurrence_perc >= 2/3)

asv_tab_10y_l_rel_occ_filt |>
  dplyr::group_by(fraction, asv_num) |>
  summarize(n = n_distinct(fraction, asv_num)) |> 
  group_by(fraction) |>
  summarize(n = sum(n))## I keep 47 ASVs in FL and 19 in PA data

occurrence_bloo_bbmo_summary <-  asv_tab_10y_l_rel_occ_filt |>
  #dplyr::select(-sample_id, -total_reads, -day, -day_of_year, -abundance_value, -n_occurrence) |>
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

##add metadata at the last cols of the dataset
m_bbmo_10y |>
  colnames()
m_bbmo_10y_sim |>
  dim()

m_bbmo_10y_sim <- m_bbmo_10y |>
  dplyr::select(Day, Month, Year, day_length, temperature, secchi, salinity, chla_total, chla_3um, PO4, NH4, NO2, NO3, Si, BP_FC1.55, PNF_Micro,         
                PNF2_5um_Micro, PNF_5um_Micro, dryptomonas, micromonas, HNF_Micro, HNF2_5um_Micro, HNF_5um_Micro,
                LNA, HNA, prochlorococcus_FC, Peuk1, Peuk2, bacteria_joint, synechococcus) |>
  distinct(Day, Month, Year , day_length, temperature, secchi, salinity, chla_total, chla_3um, PO4, NH4, NO2, NO3, Si, BP_FC1.55, PNF_Micro,         
           PNF2_5um_Micro, PNF_5um_Micro, dryptomonas, micromonas, HNF_Micro, HNF2_5um_Micro, HNF_5um_Micro,
           LNA, HNA, prochlorococcus_FC, Peuk1, Peuk2, bacteria_joint, synechococcus)

asv_tab_10y_rel_occ_filt_w_env <- asv_tab_10y_rel_occ_filt_w |>
  left_join(m_bbmo_10y_sim, by = c('Year' = 'Year', 'Month' = 'Month', 'Day' = 'Day'))

write.csv2(asv_tab_10y_rel_occ_filt_w_env, file = '../../EDM_carmen/asv_tab_10y_rel_occ_filt_w_env.csv')

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

write.csv(tax_occ_filt_bbmo, '../../EDM_carmen/tax_occ_filt_bbmo.csv')
  
occ_asv <- asv_tab_10y_l_rel_occ_filt |>
  ungroup() |>
  distinct(asv_num) |>
  as_vector()

  ##add indication of bloomer

## Occurrence in general (237 samples) SURT EL MATEIX QUE ANTERIORMENT ---- 
asv_tab_10y_l_rel |> 
  colnames()
m_bbmo_10y_sim <- m_bbmo_10y |>
  dplyr::select(fraction, sample_id, Year, Month, Day)

asv_tab_10y_l_rel_occ <- asv_tab_10y_l_rel |>
  dplyr::left_join(m_bbmo_10y_sim) |>
  # dplyr::mutate(fraction = case_when(str_detect(sample_id, '0.2') ~ 0.2,
  #                                    str_detect(sample_id, '_3_') ~ 3)) |>
  group_by(fraction, asv_num) |>
  dplyr::mutate(n_occurrence = case_when(relative_abundance > 0 ~ 1,
                                         relative_abundance == 0 ~ 0)) |>
  #group_by(fraction) |>
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
  #dplyr::select(-sample_id, -total_reads, -day, -day_of_year, -abundance_value, -n_occurrence) |>
  group_by(fraction, asv_num, occurrence_perc) |>
  distinct() |>
  dplyr::filter(occurrence_perc > 2/3) 

asv_tab_10y_l_rel_occ_filt |>
  dplyr::mutate(fraction_ed = case_when(str_detect(fraction, '0.2') ~ 'bp',
                                        str_detect(fraction, '3') ~ 'bn'),
                f_asv_num = paste0(fraction_ed,'_',asv_num)) |>
  dplyr::select(Year, Month, Day, f_asv_num, relative_abundance) |>
  pivot_wider(id_cols = c(Year, Month, Day), names_from = f_asv_num, values_from = relative_abundance, values_fill = 0 )
