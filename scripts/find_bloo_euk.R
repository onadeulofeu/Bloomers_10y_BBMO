# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# upload packages ----
library(readxl)
library(tidyverse)
library(janitor)
library(Bloomers)
library(magrittr)
library(scales)
library(phyloseq)
library(speedyseq)
library(vegan)
library(EcolUtils) #rarefaction
library(ggplot2)
library(scales)
library(zCompositions)

##upload data-----
euk_def <- readRDS('data/18S/Blanes_phyloseq.rds')

euk_def |>
  nsamples()

euk_def_sam_data <- euk_def@sam_data |>
  as_tibble()

euk_def_sam_data |>
  colnames()

euk_def_sam_data %$%
  Classified |>
  str_detect('-Nano') |>
  summary() #we are missing 29 samples from the nano fraction (3-20).

euk_def_sam_data |>
  dplyr::filter(str_detect(Classified, '-Nano'))

euk_def@otu_table |>
  dim() #11318 ASVs and 211 samples

euk_def <- prune_taxa(taxa_sums(euk_def@otu_table) >0, ##filtrar per les asv que son 0 tot el dataset
             euk_def) #they are already filtered

## we want to recover some samples from these dataset----
euk_raw <- read.table('data/18S/BLN_nanoASV_counts_v3.txt', header = T) |>
  as_tibble() ##in these table we have different replicates from the different samples.

tax_euk_raw <- read.csv2('data/18S/BLN_nanoASV_tax_v3.csv') |>
  as_tibble()

euk_raw |>
  colnames()

euk_raw |>
  dim()

euk_raw %>%
  dplyr::filter(rowSums(dplyr::select(., starts_with('BL'))) > 0)

euk_raw |>
  dplyr::select(-ASV) |>
  rowSums() |>
  as.tibble() |>
  dplyr::filter(value < 0) ##there are no 0 in this dataset

samples_raw_18 <- euk_raw |>
  colnames() |>
  as.tibble() |>
  dplyr::filter(value != 'ASV') |>
  separate(value,into = c('sample_id', 'fraction', 'replicate'), sep = '_') |> #looking for those samples that were sequenced many times
  group_by(sample_id) |>
  dplyr::mutate(n = n()) |>
  dplyr::filter(n >= 3)

## I try to find blooms in the data that Ramon considered in his previous work and see if we observe blooms or not in eukaryotes -----

## separate datasets by ASV_tab, taxonomy and metadata
euk_def_l <- euk_def@otu_table |>
  as_tibble()

tax_euk_def <- euk_def@tax_table |> 
  #mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(bbmo_10y@otu_table))) |> #ja té uns nº d'OTUs
  as_tibble()

tax_euk_def |>
  colnames()

m_euk_def <- euk_def@sam_data |>
  as_tibble()

m_euk_def |>
  colnames()

## tidy colnames----
colnames(euk_def_l) <- c('asv_num', "sample_id", 'reads')

colnames(tax_euk_def) <- c("asv_num",'seq', 'length', 'group', 'supergroup', 'match_eukV4',
                                'pident_eukV4', 'query_eukV4', "genus", 'species', 'pident_species',
                                "query_species", 'AN', "strain", "complete_name", "MAST_clade",
                                'MAST_ref', 'pident_MAST', 'fungi', 'occurrence', 'reads', 'pr2match',
                                'pr2sim', 'pr2query')

colnames(m_euk_def) <- c("sample_id", "original_name", "project", "station",             
                          "lat", "lon", "date", "depth",               
                          "depth_zone", "fraction", "extract", "temperature",       
                          "chla", "proj", "st", "dep",        
                          "frac", "ext", "classified")

colnames(tax_euk_raw) <- c('asv_num', 'seq', 'supergroup', 'group', 'genus', 'species', 'species_id', 'species_title', 'perc_identity')

## Divide metadata into FL and PA----
m_02 <- m_euk_def  |>
  dplyr::filter(fraction == 'Pico') 

m_02 <- m_02|>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02)))

m_3 <- m_euk_def |>
  dplyr::filter(fraction == 'Nano')

m_3 <- m_3 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3)))

##we lack some samples in nano fraction which are:----
m_3 |>
  separate(date, into = c('day', 'month', 'year'), sep = "/", remove = F) |>
  group_by(year) |>
  summarize(n = n()) |>
  dplyr::filter(n < 12)

## # A tibble: 4 × 2
# year      n
# <chr> <int>
#   1 2004     11
# 2 2005     11
# 3 2010      4
# 4 2012      5

## they have filtered some samples from year 2010 and 2012. 2011 doesn't have any sample.
m_02 |>
  separate(date, into = c('day', 'month', 'year'), sep = "/", remove = F) |>
  group_by(year) |>
  summarize(n = n()) |>
  dplyr::filter(n < 12) #they kept all samples from the small fraction

## I recover this group of samples from the 'raw' dataset
##check which samples were discarded and which where kept from this dataset. Try to recover those that were not considered.

## calculate relative abundances for the raw samples and observe how different they are between replicates
euk_raw_rel_tax <- euk_raw |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'reads') |>
  calculate_rel_abund(group_cols = sample_id) |>
  dplyr::filter(str_detect(sample_id, '_nano_')) |> # we want to compare only the nano fraction
  left_join(tax_euk_raw, by = c('ASV' = 'asv_num' )) |>
  separate(sample_id, into = c('sample_id_code', 'fraction', 'replicate'), remove = F) |>
  ungroup() |>
  dplyr::mutate(
    project = substr(sample_id_code, 1, 2),
    year = substr(sample_id_code, 3, 4),
    month = substr(sample_id_code, 5, 6),
    day = substr(sample_id_code, 7, 8),
    date = paste0(day,'-',month,'-',year, sep = ''))

 euk_raw_rel_tax$sample_id
 
euk_raw_rel_tax |>
  ggplot(aes(replicate, as.numeric(relative_abundance), fill = supergroup))+
  geom_col()+
  facet_wrap(vars(date)) ##we observe if the replicates are different between them or similar. Is seems that they are mainly equal. 
##I'll keep r3 in case there are two of them. 

euk_raw_rel_tax_rep <- euk_raw_rel_tax |>
  group_by(date) |>
  dplyr::filter(case_when(
    'r3' %in% replicate ~ replicate == 'r3',
    'r2' %in% replicate ~ replicate == 'r2',
    'r1' %in% replicate ~ replicate == 'r1',
    TRUE ~ n() == 1  # Keep if there's only one replicate
  ))


euk_raw_rel_tax_rep |>
  dim()
# 
# samples_raw_18_nano <- samples_raw_18 |>
#   dplyr::filter(fraction == 'nano')

# I'll add these samples to the general dataset and observe what do we see in the eukaryotic commmunity. Do we trust these samples?
## I need the same columns to add the samples that were eliminated to the dataset
euk_raw_rel_tax_rep |>
  colnames()

euk_raw_rel_tax_rep %$%
  sample_id

euk_def_l_rel <- euk_def_l |>
  calculate_rel_abund(group_cols = sample_id) |>
  left_join(tax_euk_def, by = 'asv_num') |>
  ungroup() |>
  dplyr::mutate(
    project = substr(sample_id, 1, 3),
    year = substr(sample_id, 5, 6),
    month = substr(sample_id, 7, 8),
    day = substr(sample_id, 9, 10),
    date = paste0(day,'-',month,'-',year, sep = ''))

euk_def_l_rel |>
  colnames()
  
euk_def_l_rel_f <- euk_def_l_rel |>
  dplyr::select(date, year, month, day, relative_abundance, asv_num, group, supergroup, sample_id, reads.x, seq) |>
  rename(reads = reads.x)
  
euk_raw_rel_tax_rep_f <- euk_raw_rel_tax_rep |>
  #rename(group = Group, genus = Genus, species = Species, seq = ASV, supergroup = Supergroup, asv_num = ASVid) |>
  dplyr::select(date, year, month, day, relative_abundance, seq = ASV, sample_id, reads) 

euk_raw_rel_tax_rep_f |>
  colnames()

euk_all <- euk_def_l_rel_f |>
  bind_rows(euk_raw_rel_tax_rep_f)

View(euk_all)

##I create a common taxonomy table with the taxonomy of all ASVs----
### check if the sequences are the same 
tax_euk_raw_f <- tax_euk_raw |>
  dplyr::select(asv_num, seq, group, supergroup, genus)

tax_euk_all <- tax_euk_def |>
  dplyr::select(asv_num, seq, group, supergroup, genus) |>
  bind_rows(tax_euk_raw) |>
  distinct(seq, asv_num, group, supergroup, genus)

# Transform data to CLR we use the zCompositions package which helps us to deal with 0 before applying this transformation.----
## transform asv_tab into wider format to go into vegan package
## I separate nano and pico dataframes because otherwise it seems that the memory is exhausted and divided in two different groups
euk_all <- euk_all |>
  dplyr::mutate(fraction = case_when(str_detect(sample_id, '_P_') ~ 'pico',
                                    str_detect(sample_id, '_N_') ~ 'nano',
                                     str_detect(sample_id, 'nano') ~ 'nano'))

# for those samples that I already had them in the definitive group then I don't need to recover those that were in the raw file
euk_all |>
  dim()

euk_all |>
  group_by(sample_id) |>
  summarize(n = n()) |>
  distinct(n)

euk_all |>
  dplyr::filter(fraction == 'nano') |>
  dplyr::select(date, sample_id) |>
  group_by(date, sample_id) |>
  dplyr::summarise(n = n()) # I confirm that there are duplicates

euk_all_nano_n <- euk_all %>%
  dplyr::filter(fraction == 'nano') %>%
  #dplyr::filter(date == '01-08-06') |> #060801
  dplyr::filter(str_detect(sample_id, '_N_')) |>
  distinct(date, sample_id)

euk_all_nano_no_n <- euk_all %>%
  dplyr::filter(fraction == 'nano') %>%
  #dplyr::filter(date == '01-08-06') |> #060801
  dplyr::filter(!date %in%  euk_all_nano_n$date) |>
  distinct(date, sample_id)

test <- euk_all_nano_n |>
  bind_rows(euk_all_nano_no_n) |>
  group_by(date, sample_id) %>%
  dplyr::filter(n() > 1) %>%
  ungroup()

euk_all_nano_n_data <- euk_all %>%
  dplyr::filter(fraction == 'nano') |>
  dplyr::filter(sample_id %in% euk_all_nano_n$sample_id) 

euk_all_nano_no_n_data <- euk_all %>%
  dplyr::filter(fraction == 'nano') |>
  dplyr::filter(sample_id %in% euk_all_nano_no_n$sample_id) 

euk_all_nano <- euk_all_nano_n_data |>
  bind_rows(euk_all_nano_no_n_data)

euk_all_nano |>
  group_by(date, sample_id) %>%
  dplyr::filter(n() > 1) %>%
  ungroup()

euk_all_nano |>
  distinct(sample_id) |>
  dim() #116
  
euk_all_pico <- euk_all |>
  dplyr::filter(fraction == 'pico') |>
  group_by(date) |>
  dplyr::filter(
    case_when( #when there's a _N_ that was chosen I want that, but when there's not then I'll keep the raw one
      str_detect(sample_id, '_P_') ~ TRUE,
      str_detect(sample_id, '_pico_') ~ TRUE,
      TRUE ~ n() == 1  # Keep if there's only one sample for the date
    )
  )

euk_all_pico |>
  distinct(sample_id) |>
  dim() #120

euk_all_w_nano <- euk_all_nano |>
  ungroup() |>
  # dplyr::filter(fraction == 'nano') |>
  # group_by(asv_num) |>
  # dplyr::filter(sum(reads) > 0) |>
  # ungroup() |>
  # dplyr::filter(year %in% c('04', '05', '06', '07', '08')) |>
  dplyr::select(sample_id, seq, reads) |>
  dplyr::mutate(reads = as.numeric(reads),
                sample_id = as.character(sample_id),
                seq = as.character(seq)) |>
  as_tibble() |>
  pivot_wider(names_from = 'seq', values_from = 'reads') |>
  as.data.frame()

euk_all_w_pico <- euk_all_pico |>
  ungroup() |>
  # dplyr::filter(fraction == 'pico') |>
  # group_by(asv_num) |>
  # dplyr::filter(sum(reads) > 0) |>
  # ungroup() |>
  # dplyr::filter(year %in% c('04', '05', '06', '07', '08')) |>
  dplyr::select(sample_id, seq, reads) |>
  pivot_wider(names_from = 'seq', values_from = 'reads', values_fill = 0) |>
  as.data.frame()

##check if I have more than one sample for the same date
euk_all_nano |>
  dplyr::select(date, sample_id) |>
  distinct(date, sample_id) |>
  dplyr::filter(duplicated(date))

# euk_all_w_nano |>
#   left_join(euk_all_nano, by = 'sample_id')

# euk_all_w_pico_sub2 <- euk_all |>
#   dplyr::filter(fraction == 'pico') |>
#   dplyr::filter(year %in% c('09', '10', '11', '12')) |>
#   group_by(asv_num) |>
#   dplyr::filter(sum(reads) > 0) |>
#   ungroup() |>
#   dplyr::select(date, asv_num, reads) |>
#   pivot_wider(names_from = 'asv_num', values_from = 'reads', values_fill = 0) |>
#   as.data.frame()
# 
# euk_all_w_pico_sub3 <- euk_all |>
#   dplyr::filter(fraction == 'pico') |>
#   dplyr::filter(year %in% c('13')) |>
#   group_by(asv_num) |>
#   dplyr::filter(sum(reads) > 0) |>
#   ungroup() |>
#   dplyr::select(date, asv_num, reads) |>
#   pivot_wider(names_from = 'asv_num', values_from = 'reads', values_fill = 0) |>
#   as.data.frame()

rownames(euk_all_w_nano) <- euk_all_w_nano$sample_id
rownames(euk_all_w_pico) <- euk_all_w_pico$sample_id
# 
# # asv_tab_bbmo_10y_w |>
# #   dim()
# 
euk_all_w_nano <- euk_all_w_nano[,-1]
euk_all_w_pico <- euk_all_w_pico[,-1]
# 
# #geometric mean
# gm <- function(x){
#   exp(mean(log(x[x>0])))
# }

# ## with this transformation I'm losing samples (due to too much 0 in some samples, z.warning set up to 0.99 to keep all samples)
# ### at 0.8 (default) I lose 30 samples which belonged to the years corresponding to harbour remodelation 
# ### I don't lose samples but I lose ASVs.
# zclr_df_nano <- cmultRepl(euk_all_w_nano, method = 'CZM', output = 'p-count', z.warning = 0.99
#                      #adjust = 0.2,    t = 237, s = 7849
# ) |>
#   as_tibble(rownames = "sample_id") %>%
#   pivot_longer(-sample_id) %>%
#   group_by(sample_id) %>%
#   dplyr::mutate(zclr = log(value/gm(value))) %>%
#   ungroup() %>%
#   dplyr::select(-value) %>%
#   pivot_wider(names_from = name, values_from = zclr, values_fill = 0) %>%
#   column_to_rownames("sample_id")
# 
# zclr_df_pico <- cmultRepl(euk_all_w_pico, method = 'CZM', output = 'p-count', z.warning = 0.99
#                           #adjust = 0.2,    t = 237, s = 7849
# ) |>
#   as_tibble(rownames = "sample_id") %>%
#   pivot_longer(-sample_id) %>%
#   group_by(sample_id) %>%
#   dplyr::mutate(zclr = log(value/gm(value))) %>%
#   ungroup() %>%
#   dplyr::select(-value) %>%
#   pivot_wider(names_from = name, values_from = zclr, values_fill = 0) %>%
#   column_to_rownames("sample_id")
# 
# ## with the deconstant function from vegan using pseudocunt we don't lose samples but, in Coenen 2020 they say that
# ## adding a pseudocount disproportionately affects rare taxa, where the magnitude of differences between samples may 
# ## be similar to the magnitude of the added pseudocount and therefore obscured.
# # zclr_df <- decostand(asv_tab_bbmo_10y_w, method = 'rclr' )
# # 
# # zclr_df |>
# #   dim()
# # 
# # asv_tab_bbmo_10y_w |>
# #   rownames()
# 
# ## I try the deconstant function with pseudoconunt and observe what do we get (go to line 881)----
# library(vegan)
# euk_all_w_nano_clr <- decostand(euk_all_w_nano, method = 'clr', pseudocount = 1) |>
#   rownames_to_column(var = 'sample_id') |>
#   pivot_longer(cols = starts_with(c('A','C','G','T')), names_to = 'seq', values_to = 'clr')
# 
# euk_all_w_pico_clr <- decostand(euk_all_w_pico, method = 'clr', pseudocount = 1)|>
#   rownames_to_column(var = 'sample_id') |>
#   pivot_longer(cols = starts_with(c('A','C','G','T')), names_to = 'seq', values_to = 'clr')
# 
# ##we create two datasets one for FL and one for PA--- this i the continuation if we work with zclr transformation
# asv_tab_10y_nano_zclr <- zclr_df_nano |>
#   rownames_to_column(var = 'sample_id') |>
#   pivot_longer(cols = starts_with(c('A','C','G','T')), names_to = 'seq', values_to = 'zclr')
# 
# asv_tab_10y_pico_zclr <- zclr_df_pico |>
#   rownames_to_column(var = 'sample_id') |>
#   pivot_longer(cols = starts_with(c('A','C','G','T')), names_to = 'seq', values_to = 'zclr')
# 
# ##dimensions are different from relative and pseudo dataset and zclr understand why----
# ### i have less ASVs in the zclr dataset specifically 5282 less
# #7849-2567
# 
# asv_tab_10y_nano_zclr |>
#   group_by(seq) |>
#   summarize(n = n()) |>
#   dim() #4499
# 
# asv_tab_10y_pico_zclr  |>
#   group_by(seq) |>
#   summarize(n = n()) |>
#   dim() #6868
# 
# # Discover anomalies----
# ## For each ASVs based on relative abundances and pseudoabundances -----
# ### those ASVs that are present > 50% of the sampless
# asv_tab_10y_nano_zclr %$%
#   sample_id |>
#   unique() |>
#   summary() #195 samples
# 
# asv_tab_10y_nano_zclr |>
#   colnames()
# 
# ### I filter by relative abundance > 10% at some point so that the computer can compute this easily
# 
# #x <- 120*0.75 ##percentage of ASV present in the dataset that we want to subset by (occurrence)
# 
# z_nano <- asv_tab_10y_nano_zclr |>
#   inner_join(euk_all_nano, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
#   group_by(seq) |>
#   #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
#   #dplyr::filter(num_0 <= x) |>
#   #as_tibble() |>
#   dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
#   #group_by(asv_num) |>
#   dplyr::reframe(anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
#                  anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = zclr, plotting = FALSE)[c(1,2,3)])
# 
# asv_tab_10y_nano_zclr %$%
#   sample_id |>
#   unique() |>
#   summary() #195
# 
# z_pico <- asv_tab_10y_pico_zclr |>
#   inner_join(euk_all_pico, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
#   group_by(seq) |>
#   #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
#   #dplyr::filter(num_0 <= x) |>
#   #as_tibble() |>
#   dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
#   #group_by(asv_num) |>
#   dplyr::reframe(anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
#                  anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = zclr, plotting = FALSE)[c(1,2,3)])
# 
# asv_tab_10y_pico_zclr %$%
#   sample_id |>
#   unique() |>
#   summary() #119
# 
# ##check if I'm finding the correct number of potential bloomers ASVs----
# n_bloomers_pico <-  asv_tab_10y_pico_zclr |>
#   inner_join(euk_all_pico, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
#   group_by(seq) |>
#   #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
#   #dplyr::filter(num_0 <= x) |>
#   #as_tibble() |>
#   dplyr::filter(any(relative_abundance >=  0.1)) |>
#   dplyr::filter(zclr >= 1.96) |>
#   dplyr::distinct(seq) |>
#   dplyr::summarize(n = n()) |>
#   dplyr::summarize(sum = sum(n))
# 
# bloo_pico <- asv_tab_10y_pico_zclr |>
#   inner_join(euk_all_pico, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
#   group_by(seq) |>
#   #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
#   #dplyr::filter(num_0 <= x) |>
#   #as_tibble() |>
#   dplyr::filter(any(relative_abundance >=  0.1)) |>
#   dplyr::filter(zclr >= 1.96) |>
#   dplyr::distinct(seq) |>
#   as_vector()
# 
# #write_csv2(as_tibble(bloo_02), 'data/bloo_02.csv')
# 
# n_bloomers_nano <-  asv_tab_10y_nano_zclr |>
#   inner_join(euk_all_nano, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
#   group_by(seq) |>
#   #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
#   #dplyr::filter(num_0 <= x) |>
#   #as_tibble() |>
#   dplyr::filter(any(relative_abundance >=  0.1)) |>
#   dplyr::filter(zclr >= 1.96) |>
#   dplyr::distinct(seq) |>
#   dplyr::summarize(n = n()) |>
#   dplyr::summarize(sum = sum(n))
# 
# bloo_nano <- asv_tab_10y_nano_zclr |>
#   inner_join(euk_all_nano, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
#   group_by(seq) |>
#   #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
#   #dplyr::filter(num_0 <= x) |>
#   #as_tibble() |>
#   dplyr::filter(any(relative_abundance >=  0.1)) |>
#   dplyr::filter(zclr >= 1.96) |>
#   dplyr::distinct(seq) |>
#   as_vector()
# 
# #write_csv2(as_tibble(bloo_3), 'data/bloo_3.csv')
# 
# # Filter the ASV_tab by only those ASVs that have an anomaly at some point of the dataset ----
# ##function
# # find_asv_with_anomalies <- function(anomalies_result, anomaly_in1, 
# #                                     anomaly_in2, anomaly_in3 = NULL, 
# #                                     logic1 = TRUE, 
# #                                     logic2 = NULL,
# #                                     logic3 = NULL, 
# #                                     asv_col = seq) {
# #   # if(is.list(anomalies_result) == FALSE){
# #   #   stop("Function stopped: anomalies_result needs to be a list form the get_anomalies function")
# #   # }
# #   # if(is.logical({{anomaly_in1}}) == FALSE){
# #   #   stop("Function stopped: anomaly_in1 needs to be logical (TRUE/FALSE)")
# #   # }
# #   
# #   asv_potential_bloomers <-
# #     anomalies_result |>
# #     dplyr::filter(if (!is.null(logic1)) {{anomaly_in1}} %in% logic1 else TRUE) |>
# #     dplyr::filter(if (!is.null(logic2)) {{anomaly_in2}} %in% logic2 else TRUE) |>
# #     dplyr::filter(if (!is.null(logic3)) {{anomaly_in3}} %in% logic3 else TRUE) |>
# #     dplyr::select({{asv_col}}) |>
# #     as_vector()
# #   
# #   return(asv_potential_bloomers)
# # }
# 
# ## I only filter for those anomalies in relative abundance because pseudoabundance I can only use it for fl not for PA and zclr has a problem with 
# ## dealing with many 0.
# 
# asv_anom_pico <- find_asv_with_anomalies(anomalies_result = z_pico, anomaly_in1 = anomalies_ra, 
#                                        anomaly_in2 = NULL,
#                                        anomaly_in3 = anomalies_clr, 
#                                        logic1 = 'TRUE', logic2 = NULL, 
#                                        logic3 = NULL,
#                                        asv_col = seq)
# ##me'n surten 34
# 
# asv_anom_nano <- find_asv_with_anomalies(anomalies_result = z_nano, anomaly_in1 = anomalies_ra, 
#                                       anomaly_in2 = anomalies_clr, 
#                                       anomaly_in3 = NULL, #anomalies_ps
#                                       logic1 = 'TRUE', logic2 = NULL, 
#                                       logic3 = NULL, ##per 3 no és representatiu la pseudoabund
#                                       asv_col = seq)
# 
# ##48 ASVS
# ## número d'ASVs comuns i unics a cada fracció
# asv_anom_nano_tb <- asv_anom_nano |>
#   as_tibble()
# 
# asv_anom_pico_tb <- asv_anom_pico |>
#   as_tibble()
# 
# common_bloomers_tax_euk <- asv_anom_nano_tb |>
#   bind_rows(asv_anom_pico_tb) |>
#   unique() |> ##96 totals >50% del dataset
#   left_join(tax_euk_def, by = c('value' = 'seq'))
# 
# asv_anom_nano_tb |>
#   anti_join(asv_anom_pico_tb) #40 ASV only in 3
# 
# asv_anom_pico_tb |>
#   anti_join(asv_anom_nano_tb) #10 ASVs only in 0.2
# 
# asv_tab_10y_pico_zclr_bloo <-  asv_tab_10y_pico_zclr |>
#   left_join(tax_euk_all, by = c('seq')) |> ##add taxonomy
#   left_join(euk_all_pico, by = c('sample_id', 'seq')) |>
#   #group_by(seq) |>
#   #dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
#   pivot_longer(cols = c(zclr, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
#   dplyr::filter(seq %in% asv_anom_pico |
#                   seq %in% asv_anom_nano) ##recover ASVs that presented anomalies in 02 or 3
# 
# asv_tab_10y_pico_zclr_bloo |>
#   group_by(seq) |>
#   dplyr::summarize(n = n()) |>
#   dplyr::summarize(n_num = n())
# 
# asv_tab_10y_nano_zclr_bloo <-  asv_tab_10y_nano_zclr |>
#   left_join(tax_euk_all, by = c('seq')) |> ##add taxonomy
#   left_join(euk_all_nano, by = c('sample_id', 'seq')) |>
#   group_by(seq) |>
#   dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
#   pivot_longer(cols = c(zclr, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
#   dplyr::filter(seq %in% asv_anom_nano |
#                   seq %in% asv_anom_pico) ##recover ASVs that presented anomalies in 02 or 3
# 
# asv_tab_10y_nano_zclr_bloo |>
#   group_by(seq) |>
#   dplyr::summarize(n = n()) |>
#   dplyr::summarize(n_num = n())
# 
# asv_tab_all_bloo <- asv_tab_10y_pico_zclr_bloo |>
#   bind_rows(asv_tab_10y_nano_zclr_bloo)
# 
# asv_tab_all_bloo |>
#   colnames()
# 
# asv_tab_all_bloo |>
#   group_by(seq) |>
#   dplyr::summarize(n = n()) |>
#   dplyr::summarize(n_num = n()) #17 ASVs
# 
# # asv_tab_all_perc_filt_nano_long_filt |>
# #   group_by(sample_id) |>
# #   dplyr::summarize(sum = sum(relative_abundance))
# 
# asv_tab_all_bloo |> 
#   dplyr::filter(abundance_type == 'relative_abundance') %$%
#   abundance_value |>
#   range() #max is 0.83
# 
# ## Recover z-score for each ASV ----
# ### I want to highlight anomalies for each ASV to do so I recover z-scores for those ASVs that that have high z-scores
# ### at some point of the dataset. Easy to observe if those ASVs are having random anomalies or all of them happen at the same time
# sample_id_pico <- asv_tab_10y_pico_zclr_bloo |>
#   ungroup() |>
#   dplyr::select(sample_id) |>
#   distinct(sample_id)
# 
# m_pico <- euk_all_pico |>
#   ungroup() |>
#   dplyr::select(sample_id, date, year, month, day, fraction) |>
#   distinct(sample_id, date, year, month, day, fraction) |>
#   right_join(sample_id_pico) |> #there's one missing sample, probably because of the transformation to CLR 
#   dplyr::mutate(sample_id_num = str_c(1:nrow(sample_id_pico)))
# 
# x <- 120*0.75
# z_scores_pico <- asv_tab_10y_pico_zclr |>
#   left_join(euk_all_pico) |>
#   group_by(seq) |>
#   dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
#   #dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 50% of the samples
#   #dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
#   group_by(seq) |>
#   dplyr::reframe(z_score_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
#                                             na_rm = TRUE, values = relative_abundance, 
#                                             plotting = FALSE)[c(2)]) |>
#   as_tibble() |>
#   unnest(cols = z_score_ra) |>
#   group_by(seq) |>
#   dplyr::mutate(sample_id_num = str_c(1:nrow(sample_id_pico))) |>
#   left_join(m_pico, by = 'sample_id_num')
# 
# asv_tab_10y_02_pseudo_zclr
# 
# x <- 117*0.75 #number of 0s that we accept in the dataset to filter it by to calculate anomalies
# 
# sample_id_nano <- asv_tab_10y_nano_zclr_bloo |>
#   ungroup() |>
#   dplyr::select(sample_id) |>
#   distinct(sample_id)
# 
# m_nano <- euk_all_nano |>
#   ungroup() |>
#   dplyr::select(sample_id, date, year, month, day, fraction) |>
#   distinct(sample_id, date, year, month, day, fraction) |>
#   right_join(sample_id_nano) |> #there's one missing sample, probably because of the transformation to CLR 
#   dplyr::mutate(sample_id_num = str_c(1:nrow(sample_id_nano)))
# 
# z_scores_nano <- asv_tab_10y_nano_zclr |>
#   left_join(euk_all_nano) |>
#   group_by(seq) |>
#   dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
#   #dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 50% of the samples
#   #dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
#   group_by(seq) |>
#   dplyr::reframe(z_score_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
#                                             na_rm = TRUE, values = relative_abundance, 
#                                             plotting = FALSE)[c(2)]) |>
#   as_tibble() |>
#   unnest(cols = z_score_ra) |>
#   group_by(seq) |>
#   dplyr::mutate(sample_id_num = str_c(1:nrow(sample_id_nano))) |>
#   left_join(m_nano, by = 'sample_id_num') 
# 
# z_scores_all <- z_scores_pico |>
#   bind_rows(z_scores_nano)
# 
# ##Why do I have inf values when I assess z-scores?
# ###understand why do we have zscores near infinite (it is because we have enormous changes of relative abundance coming from 0 values)
# z_scores_pico |>
#   dplyr::filter(z_score_ra ==  'Inf') |>
#   dplyr::filter(z_score_ra >= 1.96) #check that filtering by z_score is detecting infinitive values as a number
# 
# z_score_infinite <- asv_tab_all_bloo_z_tax |>
#   dplyr::select(z_score_ra, seq, sample_id, abundance_type, abundance_value) |>
#   #dplyr::filter(abundance_type == 'relative_abundance') |>
#   dplyr::filter(z_score_ra == is.infinite(z_score_ra)) 
# 
# z_scores_all |>
#   colnames()
# 
# # Dataset with all information ----
# asv_tab_all_bloo_z_tax <- asv_tab_all_bloo |>
#   left_join(z_scores_all)
# # |> 
#   # left_join(euk_all_tax, by = 'seq') I don't add the taxonomy because it is already present
# 
# 
# # write.csv2(asv_tab_all_bloo_z_tax, 'data/18S/euk_asv_tab_all_bloo_z_tax.csv')
# 
# asv_tab_all_bloo_z_tax |>
#   distinct(asv_num.x) |>
#   dim()
# 
# asv_tab_all_bloo_z_tax |>
#   colnames()
# 
# ##Upload euk bloomers data-----
# euk_asv_tab_all_bloo_z_tax <- read.csv2('data/18S/euk_asv_tab_all_bloo_z_tax.csv') |>
#   as_tibble() |>
#   dplyr::select(sample_id, seq, date, year, month, day, fraction, abundance_value, z_score_ra, abundance_type) |>
#   left_join(tax_euk_def, by = 'seq')  #just the most actualized taxonomy
# 
# euk_asv_tab_all_bloo_z_tax |>
#   colnames()
# 
# #labels----
# labs_fraction_euk <- as_labeller(c('pico' = 'Picoeukaryotes (0.2-3 um)',
#                                'nano' = 'Nanoeukaryotes (3-20 um)'))
# 
# #reorder taxonomy as factor
# euk_asv_tab_all_bloo_z_tax <- euk_asv_tab_all_bloo_z_tax |>
#   dplyr::mutate(supergroup_f = as_factor(supergroup),
#                 group_f = as_factor(group),
#                 genus_f = as_factor(genus))
# 
# euk_asv_tab_all_bloo_z_tax$group_f <-  factor(euk_asv_tab_all_bloo_z_tax$group_f, 
#                                           levels=unique(euk_asv_tab_all_bloo_z_tax$group_f[order(euk_asv_tab_all_bloo_z_tax$supergroup_f)]), 
#                                           ordered=TRUE)
# 
# euk_asv_tab_all_bloo_z_tax$genus_f <-  factor(euk_asv_tab_all_bloo_z_tax$genus_f, 
#                                           levels=unique(euk_asv_tab_all_bloo_z_tax$genus_f[order(euk_asv_tab_all_bloo_z_tax$supergroup_f,
#                                                                                              euk_asv_tab_all_bloo_z_tax$group_f)]), 
#                                           ordered=TRUE)
# ##palette
# ### tinc un problema amb la taxonomia perquè sembla que la raw i la definitive són de diferents bases de dades
# euk_asv_tab_all_bloo_z_tax |>
#   distinct(group_f, supergroup_f) |>
#   View()
#   as_vector()
# 
# palette_supergroup_assigned <- c('Alveolata' = "#fcca46",
#                                  "Archaeplastida" = '#0051BF', 
#                                  'Rhizaria' = '#B0413E',
#                               'Opisthokonta' = "#009e73",
#                                'Stramenopiles'= '#BE8DCB',
#                               'Amoebozoa' ="#fb9a99",
#                               'Cryptista' = "#87878b")
# 
# palette_group_assigned <- c("Ciliophora"  = "#fcca46",
#                             "MALV-I"    =   '#FFA737',
#                             "Ulvophyceae" = "#0051BF" ,
#                             "Acantharea"  = "#b0413e",
#                             "Basidiomycota" =    "#009e73",        
#                             "Diatomea" = "#69267e", 
#                             "Picozoa"  = '#002671', 
#                             "Colpodellida"    = '#ce7800', 
#                             "Chytridiomycota"   = '#00733C',    
#                             "Ichthyosporea"    = "#003029",         
#                             "MAST-3" = '#e3a6ce',
#                             "Choanoflagellata"  = '#4dbaa9',
#                             "Gracilipodida" =  "#fb9a99",  
#                             "Dinoflagellata"   = '#fbed5c',   
#                             "MALV-II"   = '#BE8DCB',   
#                             "Ascomycota" = '#2E5A51',
#                             "Chlorodendrophyceae" = '#3D518E',
#                             "Mamiellophyceae" =  '#1F78B4',
#                             "Cryptomonadales" = "#87878b",         
#                             "Cercozoa"  =  '#c55e5c',
#                             "Bicosoecida" = '#654584')
#   
# ## Plot the euk bloomers and observe what do we see-----
# euk_asv_tab_all_bloo_z_tax |>
#   colnames()
# 
# ##color by supergroup----
# euk_asv_tab_all_bloo_z_tax$date
# 
# euk_asv_tab_all_bloo_z_tax$fraction <- factor(euk_asv_tab_all_bloo_z_tax$fraction, levels = c('pico', 'nano'))
# 
# euk_asv_tab_all_bloo_z_tax |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   group_by(date, fraction) |>
#   dplyr::mutate(max_abund = sum(abundance_value)) |>
#   ungroup() |>
#   group_by(date, fraction, supergroup) |>
#   dplyr::mutate(abund_order = sum(abundance_value)) |>
#   ungroup() |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   ggplot(aes(date, max_abund))+
#   #geom_line(aes(date, max_abund))+
#   #geom_segment(aes(x = '2005-01-01', y = 0, xend = '2005-01-02', yend =0.57),color="black")+
#   scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
#                    #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
#                    #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
#   )+
#   geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
#                                                     ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
#   
#   #geom_stream(aes(fill = class_f, group = class_f), type = "ridge", bw=1)+
#   geom_area(aes(date, abund_order, fill = supergroup_f, group = supergroup_f), alpha = 0.8,  position='stack')+
#   #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
#   # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
#   # geom_point(data = community_eveness_all_m |>
#   #              dplyr::filter(anomaly_color == '#9F0011'),  
#   #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
#    scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
#   #,
#   #                    sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
#   scale_color_identity()+
#   scale_fill_manual(values = 
#                       palette_supergroup_assigned, na.value = "#000000")+
#   labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Supergroup')+
#   facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_euk)+
#   #facet_wrap(fraction~phylum_f, dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
#   guides(fill = guide_legend(ncol = 6, size = 10,
#                              override.aes = aes(label = '')),
#          alpha = 'none')+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(), strip.text = element_text(size = 7),
#         legend.position = 'bottom', axis.text.y = element_text(size = 8),
#         axis.title = element_text(size = 8), strip.background = element_blank(), 
#         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside') 
# 
# ##color by group----
# euk_asv_tab_all_bloo_z_tax |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   group_by(date, fraction) |>
#   dplyr::mutate(max_abund = sum(abundance_value)) |>
#   ungroup() |>
#   group_by(date, fraction, group) |>
#   dplyr::mutate(abund_order = sum(abundance_value)) |>
#   ungroup() |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   ggplot(aes(date, max_abund))+
#   #geom_line(aes(date, max_abund))+
#   #geom_segment(aes(x = '2005-01-01', y = 0, xend = '2005-01-02', yend =0.57),color="black")+
#   scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
#                    #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
#                    #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
#   )+
#   geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
#                                                     ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
#   
#   #geom_stream(aes(fill = class_f, group = class_f), type = "ridge", bw=1)+
#   geom_area(aes(date, abund_order, fill = group_f, group = group_f), alpha = 0.8,  position='stack')+
#   #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
#   # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
#   # geom_point(data = community_eveness_all_m |>
#   #              dplyr::filter(anomaly_color == '#9F0011'),  
#   #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
#   scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
#   #,
#   #                    sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
#   scale_color_identity()+
#   scale_fill_manual(values = 
#                       palette_group_assigned, na.value = "#000000")+
#   labs(x = 'Time', y = 'Relative abundance (%)', fill = 'group')+
#   facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_euk)+
#   #facet_wrap(fraction~phylum_f, dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
#   guides(fill = guide_legend(ncol = 6, size = 10,
#                              override.aes = aes(label = '')),
#          alpha = 'none')+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(), strip.text = element_text(size = 7),
#         legend.position = 'bottom', axis.text.y = element_text(size = 8),
#         axis.title = element_text(size = 8), strip.background = element_blank(), 
#         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside') 

### I redo the pipeline with the clr dataset because I'm trying to  recover more sample from the nanoplankton fraction----
## This function is better at least for the prokaryotes dataset so I would continue with it.
## I try the deconstant function with pseudoconunt and observe what do we get ----
euk_all_w_nano |>
  colnames()

asv_tab_10y_nano_clr <- decostand(euk_all_w_nano, method = 'rclr') |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with(c('A','C','G','T')), names_to = 'seq', values_to = 'clr')

asv_tab_10y_pico_clr <- decostand(euk_all_w_pico, method = 'rclr')|>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with(c('A','C','G','T')), names_to = 'seq', values_to = 'clr')

## dimensions are different from relative and pseudo dataset and clr understand why----
## with this approximation I get more AVSs

asv_tab_10y_nano_clr |>
  group_by(seq) |>
  summarize(n = n()) |>
  dim() #11837

asv_tab_10y_pico_clr  |>
  group_by(seq) |>
  summarize(n = n()) |>
  dim() #11837

# Discover anomalies----
## For each ASVs based on relative abundances and pseudoabundances -----
### those ASVs that are present > 50% of the samples
asv_tab_10y_nano_clr %$%
  sample_id |>
  unique() |>
  summary() #116 samples

asv_tab_10y_nano_clr |>
  colnames()

### I filter by relative abundance > 10% at some point so that the computer can compute this easily

#x <- 120*0.75 ##percentage of ASV present in the dataset that we want to subset by (occurrence)

z_nano <- asv_tab_10y_nano_clr |>
  inner_join(euk_all_nano, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_clr vull afegir el clr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(seq) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  #group_by(asv_num) |>
  dplyr::reframe(anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = clr, plotting = FALSE)[c(1,2,3)])

asv_tab_10y_nano_clr %$%
  sample_id |>
  unique() |>
  summary() #116

z_pico <- asv_tab_10y_pico_clr |>
  inner_join(euk_all_pico, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_clr vull afegir el clr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(seq) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  #group_by(asv_num) |>
  dplyr::reframe(anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = clr, plotting = FALSE)[c(1,2,3)])

asv_tab_10y_pico_clr %$%
  sample_id |>
  unique() |>
  summary() #120

##check if I'm finding the correct number of potential bloomers ASVs----
n_bloomers_pico <-  asv_tab_10y_pico_clr |>
  inner_join(euk_all_pico, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_clr vull afegir el clr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(seq) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(clr >= 1.96) |>
  dplyr::distinct(seq) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(sum = sum(n))

bloo_pico <- asv_tab_10y_pico_clr |>
  inner_join(euk_all_pico, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_clr vull afegir el clr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(seq) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(clr >= 1.96) |>
  dplyr::distinct(seq) |>
  as_vector()

#write_csv2(as_tibble(bloo_pico), 'data/bloo_02_euk.csv')

n_bloomers_nano <-  asv_tab_10y_nano_clr |>
  inner_join(euk_all_nano, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_clr vull afegir el clr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(seq) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(clr >= 1.96) |>
  dplyr::distinct(seq) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(sum = sum(n))

bloo_nano <- asv_tab_10y_nano_clr |>
  inner_join(euk_all_nano, by = c('sample_id', 'seq')) |> #asv_tab_10y_02_clr vull afegir el clr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(seq) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(clr >= 1.96) |>
  dplyr::distinct(seq) |>
  as_vector()

#write_csv2(as_tibble(bloo_nano), 'data/bloo_3_euk.csv')

# Filter the ASV_tab by only those ASVs that have an anomaly at some point of the dataset ----
##function
# find_asv_with_anomalies <- function(anomalies_result, anomaly_in1, 
#                                     anomaly_in2, anomaly_in3 = NULL, 
#                                     logic1 = TRUE, 
#                                     logic2 = NULL,
#                                     logic3 = NULL, 
#                                     asv_col = seq) {
#   # if(is.list(anomalies_result) == FALSE){
#   #   stop("Function stopped: anomalies_result needs to be a list form the get_anomalies function")
#   # }
#   # if(is.logical({{anomaly_in1}}) == FALSE){
#   #   stop("Function stopped: anomaly_in1 needs to be logical (TRUE/FALSE)")
#   # }
#   
#   asv_potential_bloomers <-
#     anomalies_result |>
#     dplyr::filter(if (!is.null(logic1)) {{anomaly_in1}} %in% logic1 else TRUE) |>
#     dplyr::filter(if (!is.null(logic2)) {{anomaly_in2}} %in% logic2 else TRUE) |>
#     dplyr::filter(if (!is.null(logic3)) {{anomaly_in3}} %in% logic3 else TRUE) |>
#     dplyr::select({{asv_col}}) |>
#     as_vector()
#   
#   return(asv_potential_bloomers)
# }

## I only filter for those anomalies in relative abundance because pseudoabundance I can only use it for fl not for PA and clr has a problem with 
## dealing with many 0.
asv_anom_pico <- find_asv_with_anomalies(anomalies_result = z_pico, anomaly_in1 = anomalies_ra, 
                                         anomaly_in2 = NULL,
                                         anomaly_in3 = anomalies_clr, 
                                         logic1 = 'TRUE', logic2 = NULL, 
                                         logic3 = NULL,
                                         asv_col = seq)
##me'n surten 37

asv_anom_nano <- find_asv_with_anomalies(anomalies_result = z_nano, anomaly_in1 = anomalies_ra, 
                                         anomaly_in2 = anomalies_clr, 
                                         anomaly_in3 = NULL, #anomalies_ps
                                         logic1 = 'TRUE', logic2 = NULL, 
                                         logic3 = NULL, ##per 3 no és representatiu la pseudoabund
                                         asv_col = seq)

##68 ASVS
## número d'ASVs comuns i unics a cada fracció
asv_anom_nano_tb <- asv_anom_nano |>
  as_tibble()

asv_anom_pico_tb <- asv_anom_pico |>
  as_tibble()

common_bloomers_tax_euk <- asv_anom_nano_tb |>
  bind_rows(asv_anom_pico_tb) |>
  unique() |> ##86 totals >50% del dataset
  #left_join(tax_euk_def, by = c('value' = 'seq')) |>
  left_join(tax_euk_all, by = c('value' = 'seq'))

asv_anom_nano_tb |>
  anti_join(asv_anom_pico_tb) #49 ASV only in 3

asv_anom_pico_tb |>
  anti_join(asv_anom_nano_tb) #18 ASVs only in 0.2

asv_tab_10y_pico_clr_bloo <-  asv_tab_10y_pico_clr |>
  left_join(tax_euk_all, by = c('seq')) |> ##add taxonomy
  left_join(euk_all_pico, by = c('sample_id', 'seq')) |>
  #group_by(seq) |>
  #dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(clr, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  dplyr::filter(seq %in% asv_anom_pico |
                  seq %in% asv_anom_nano) ##recover ASVs that presented anomalies in 02 or 3

asv_tab_10y_pico_clr_bloo |>
  group_by(seq) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n())

asv_tab_10y_nano_clr_bloo <-  asv_tab_10y_nano_clr |>
  left_join(tax_euk_all, by = c('seq')) |> ##add taxonomy
  left_join(euk_all_nano, by = c('sample_id', 'seq')) |>
  group_by(seq) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(clr, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  dplyr::filter(seq %in% asv_anom_nano |
                  seq %in% asv_anom_pico) ##recover ASVs that presented anomalies in 02 or 3

asv_tab_10y_nano_clr_bloo |>
  group_by(seq) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n())

asv_tab_all_bloo <- asv_tab_10y_pico_clr_bloo |>
  bind_rows(asv_tab_10y_nano_clr_bloo)

asv_tab_all_bloo |>
  colnames()

asv_tab_all_bloo |>
  group_by(seq) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n()) #17 ASVs

# asv_tab_all_perc_filt_nano_long_filt |>
#   group_by(sample_id) |>
#   dplyr::summarize(sum = sum(relative_abundance))

asv_tab_all_bloo |> 
  dplyr::filter(abundance_type == 'relative_abundance') %$%
  abundance_value |>
  range() #max is 0.83

## Recover z-score for each ASV ----
### I want to highlight anomalies for each ASV to do so I recover z-scores for those ASVs that that have high z-scores
### at some point of the dataset. Easy to observe if those ASVs are having random anomalies or all of them happen at the same time
sample_id_pico <- asv_tab_10y_pico_clr_bloo |>
  ungroup() |>
  dplyr::select(sample_id) |>
  distinct(sample_id)

m_pico <- euk_all_pico |>
  ungroup() |>
  dplyr::select(sample_id, date, year, month, day, fraction) |>
  distinct(sample_id, date, year, month, day, fraction) |>
  right_join(sample_id_pico) |> #there's one missing sample, probably because of the transformation to CLR 
  dplyr::mutate(sample_id_num = str_c(1:nrow(sample_id_pico)))

z_scores_pico <- asv_tab_10y_pico_clr |>
  left_join(euk_all_pico) |>
  group_by(seq) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 50% of the samples
  #dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
  group_by(seq) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  group_by(seq) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(sample_id_pico))) |>
  left_join(m_pico, by = 'sample_id_num')

#x <- 117*0.75 #number of 0s that we accept in the dataset to filter it by to calculate anomalies

sample_id_nano <- asv_tab_10y_nano_clr_bloo |>
  ungroup() |>
  dplyr::select(sample_id) |>
  distinct(sample_id)

m_nano <- euk_all_nano |>
  ungroup() |>
  dplyr::select(sample_id, date, year, month, day, fraction) |>
  distinct(sample_id, date, year, month, day, fraction) |>
  right_join(sample_id_nano) |> #there's one missing sample, probably because of the transformation to CLR 
  dplyr::mutate(sample_id_num = str_c(1:nrow(sample_id_nano)))

z_scores_nano <- asv_tab_10y_nano_clr |>
  left_join(euk_all_nano) |>
  group_by(seq) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 50% of the samples
  #dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
  group_by(seq) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  group_by(seq) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(sample_id_nano))) |>
  left_join(m_nano, by = 'sample_id_num') 

z_scores_all <- z_scores_pico |>
  bind_rows(z_scores_nano)

##Why do I have inf values when I assess z-scores?
###understand why do we have zscores near infinite (it is because we have enormous changes of relative abundance coming from 0 values)
z_scores_pico |>
  dplyr::filter(z_score_ra ==  'Inf') |>
  dplyr::filter(z_score_ra >= 1.96) #check that filtering by z_score is detecting infinitive values as a number

z_score_infinite <- asv_tab_all_bloo_z_tax |>
  dplyr::select(z_score_ra, seq, sample_id, abundance_type, abundance_value) |>
  #dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(z_score_ra == is.infinite(z_score_ra)) 

z_scores_all |>
  colnames()

# Dataset with all information ----
euk_asv_tab_all_bloo_z_tax <- asv_tab_all_bloo |>
  left_join(z_scores_all)
# |> 
# left_join(euk_all_tax, by = 'seq') I don't add the taxonomy because it is already present

#write.csv2(asv_tab_all_bloo_z_tax, 'data/18S/euk_asv_tab_all_bloo_z_tax_ed.csv')

asv_tab_all_bloo_z_tax |>
  distinct(asv_num.x) |>
  dim()

asv_tab_all_bloo_z_tax |>
  colnames()

##Upload euk bloomers data-----

### hi ha algún error perque la fració es NA per algún motiu REVISAR EL CODI!-----
euk_asv_tab_all_bloo_z_tax <- read.csv2('data/18S/euk_asv_tab_all_bloo_z_tax.csv') |>
   as_tibble() |>
  dplyr::select(sample_id, seq, date, year, month, day, fraction, abundance_value, z_score_ra, abundance_type) |>
   left_join(tax_euk_def, by = 'seq')  #just the most actualized taxonomy

euk_asv_tab_all_bloo_z_tax |>
  colnames()

euk_asv_tab_all_bloo_z_tax |>
  dplyr::filter(is.na(fraction)) ## crec que podria ser normal degut a les dades que tenim

#labels----
labs_fraction_euk <- as_labeller(c('pico' = 'Picoeukaryotes (0.2-3 um)',
                                   'nano' = 'Nanoeukaryotes (3-20 um)'))

#reorder taxonomy as factor
euk_asv_tab_all_bloo_z_tax <- euk_asv_tab_all_bloo_z_tax |>
  dplyr::mutate(supergroup_f = as_factor(supergroup),
                group_f = as_factor(group),
                genus_f = as_factor(genus))

euk_asv_tab_all_bloo_z_tax$group_f <-  factor(euk_asv_tab_all_bloo_z_tax$group_f, 
                                              levels=unique(euk_asv_tab_all_bloo_z_tax$group_f[order(euk_asv_tab_all_bloo_z_tax$supergroup_f)]), 
                                              ordered=TRUE)

euk_asv_tab_all_bloo_z_tax$genus_f <-  factor(euk_asv_tab_all_bloo_z_tax$genus_f, 
                                              levels=unique(euk_asv_tab_all_bloo_z_tax$genus_f[order(euk_asv_tab_all_bloo_z_tax$supergroup_f,
                                                                                                     euk_asv_tab_all_bloo_z_tax$group_f)]), 
                                              ordered=TRUE)
##palette
### tinc un problema amb la taxonomia perquè sembla que la raw i la definitive són de diferents bases de dades
euk_asv_tab_all_bloo_z_tax |>
  distinct(group_f, supergroup_f) |>
  View()
as_vector()

palette_supergroup_assigned <- c('Alveolata' = "#fcca46",
                                 "Archaeplastida" = '#0051BF', 
                                 'Rhizaria' = '#B0413E',
                                 'Opisthokonta' = "#009e73",
                                 'Stramenopiles'= '#BE8DCB',
                                 'Amoebozoa' ="#fb9a99",
                                 'Cryptista' = "#87878b")

palette_group_assigned <- c("Ciliophora"  = "#fcca46",
                            "MALV-I"    =   '#FFA737',
                            "Ulvophyceae" = "#0051BF" ,
                            "Acantharea"  = "#b0413e",
                            "Basidiomycota" =    "#009e73",        
                            "Diatomea" = "#69267e", 
                            "Picozoa"  = '#002671', 
                            "Colpodellida"    = '#ce7800', 
                            "Chytridiomycota"   = '#00733C',    
                            "Ichthyosporea"    = "#003029",         
                            "MAST-3" = '#e3a6ce',
                            "Choanoflagellata"  = '#4dbaa9',
                            "Gracilipodida" =  "#fb9a99",  
                            "Dinoflagellata"   = '#fbed5c',   
                            "MALV-II"   = '#BE8DCB',   
                            "Ascomycota" = '#2E5A51',
                            "Chlorodendrophyceae" = '#3D518E',
                            "Mamiellophyceae" =  '#1F78B4',
                            "Cryptomonadales" = "#87878b",         
                            "Cercozoa"  =  '#c55e5c',
                            "Bicosoecida" = '#654584')

## Plot the euk bloomers and observe what do we see-----
euk_asv_tab_all_bloo_z_tax |>
  colnames()

##color by supergroup----
euk_asv_tab_all_bloo_z_tax$date

euk_asv_tab_all_bloo_z_tax$fraction <- factor(euk_asv_tab_all_bloo_z_tax$fraction, levels = c('pico', 'nano'))

euk_asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, supergroup) |>
  dplyr::mutate(abund_order = sum(abundance_value)) |>
  ungroup() |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(!is.na(fraction)) |>
  ggplot(aes(date, max_abund))+
  #geom_line(aes(date, max_abund))+
  #geom_segment(aes(x = '2005-01-01', y = 0, xend = '2005-01-02', yend =0.57),color="black")+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  
  #geom_stream(aes(fill = class_f, group = class_f), type = "ridge", bw=1)+
  geom_area(aes(date, abund_order, fill = supergroup_f, group = supergroup_f), alpha = 0.9,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),  
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
  #,
  #                    sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
  scale_color_identity()+
  scale_fill_manual(values = 
                      palette_supergroup_assigned, na.value = "#000000")+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Supergroup')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_euk)+
  #facet_wrap(fraction~phylum_f, dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside') 

##color by group----
euk_asv_tab_all_bloo_z_tax |>
  dplyr::filter(!is.na(fraction)) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, group) |>
  dplyr::mutate(abund_order = sum(abundance_value)) |>
  ungroup() |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, max_abund))+
  #geom_line(aes(date, max_abund))+
  #geom_segment(aes(x = '2005-01-01', y = 0, xend = '2005-01-02', yend =0.57),color="black")+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  
  #geom_stream(aes(fill = class_f, group = class_f), type = "ridge", bw=1)+
  geom_area(aes(date, abund_order, fill = group_f, group = group_f), alpha = 0.9,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),  
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
  #,
  #                    sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
  scale_color_identity()+
  scale_fill_manual(values = 
                      palette_group_assigned, na.value = "#000000")+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'group')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_euk)+
  #facet_wrap(fraction~phylum_f, dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside') 

## qui són els que están durant les obres del port -----

blooms_port <- euk_asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(date, fraction)
  
  
  
  blooms_port |>  
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, supergroup) |>
  dplyr::mutate(abund_order = sum(abundance_value)) |>
  ungroup() |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(!is.na(fraction)) |>
  dplyr::filter(is.na(supergroup)) |>
  ggplot(aes(date, max_abund))+
  #geom_line(aes(date, max_abund))+
  #geom_segment(aes(x = '2005-01-01', y = 0, xend = '2005-01-02', yend =0.57),color="black")+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  
  #geom_stream(aes(fill = class_f, group = class_f), type = "ridge", bw=1)+
  geom_area(aes(date, abund_order, fill = supergroup_f, group = supergroup_f), alpha = 0.9,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),  
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
  #,
  #                    sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
  scale_color_identity()+
  scale_fill_manual(values = 
                      palette_supergroup_assigned, na.value = "#000000")+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Supergroup')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_euk)+
  #facet_wrap(fraction~phylum_f, dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  ## the most abundant euk just after the harbor restoration is Helkesimastix marina n. sp. (Cercozoa: Sainouroidea superfam. n.) a Gliding Zooflagellate of Novel Ultrastructure and Unusual Ciliary Behaviour 
  
  ### observe these "potential blooming euk" and their changes over the years ------
  euk_asv_tab_all_bloo_z_tax |>
    dplyr::filter(!is.na(fraction)) |>
    dplyr::filter(abundance_type == 'clr') |>
    # group_by(date, fraction) |>
    # #dplyr::mutate(max_abund = sum(abundance_value)) |>
    # ungroup() |>
    # group_by(date, fraction, group) |>
    # dplyr::mutate(abund_order = sum(abundance_value)) |>
    ungroup() |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
    ggplot(aes(date, abundance_value))+
    facet_grid(supergroup~fraction)+
    geom_line(aes(color = supergroup, group = asv_num))+
    scale_color_manual(values = palette_supergroup_assigned)+
    labs(x = 'Time', y = 'rCLR')+
    theme_bw()+
    theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), strip.text = element_text(size = 7),
          legend.position = 'bottom', axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 8), strip.background = element_blank(), 
          legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside') 
  

  