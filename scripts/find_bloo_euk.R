##upload data
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

## see if we want to recover some samples from these part
euk_raw <- read.table('data/18S/BLN_nanoASV_counts_v3.txt', header = T) |>
  as_tibble() ##in these table we have different replicates from the different samples.

tax_euk_raw <- read.csv2('data/18S/BLN_nanoASV_tax_v3.csv') |>
  as_tibble()

euk_raw |>
  colnames()

euk_raw |>
  dim()

euk_raw %>%
  filter(rowSums(select(., starts_with('BL'))) > 0)

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
euk_def

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
euk_def_sam_data$Original_name

samples_raw_18

euk_raw 

## calculate relative abundances for the raw samples and observe how different they are between replicates
euk_raw_rel_tax <- euk_raw |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'reads') |>
  calculate_rel_abund(group_cols = sample_id) |>
  dplyr::filter(str_detect(sample_id, '_nano_')) |> # we want to compare only the nano fraction
  left_join(tax_euk_raw, by = 'ASV') |>
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
  ggplot(aes(replicate, as.numeric(relative_abundance), fill = Supergroup))+
  geom_col()+
  facet_wrap(vars(date)) ##we observe if hte replicates are different between them or similar. Is seems that they are mainly equal. 
##I'll keep r3 in case there are two of them. 

euk_raw_rel_tax_rep <- euk_raw_rel_tax |>
  group_by(date) |>
  dplyr::filter(case_when(
    'r3' %in% replicate ~ replicate == 'r3',
    'r2' %in% replicate ~ replicate == 'r2',
    'r1' %in% replicate ~ replicate == 'r1',
    TRUE ~ n() == 1  # Keep if there's only one replicate
  ))

# 
# samples_raw_18_nano <- samples_raw_18 |>
#   dplyr::filter(fraction == 'nano')

# I'll add these samples to the general dataset and observe what do we see in the eukaryotic commmunity. Do we trust these samples?
## I need the same columns to add the samples that were eliminated to the dataset
euk_raw_rel_tax_rep |>
  colnames()

euk_raw_rel_tax_rep %$%
  sample_id

euk_def_l_rel %$%
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
  dplyr::select(date, year, month, day, relative_abundance, asv_num, group, supergroup, sample_id)
  
euk_raw_rel_tax_rep |>
  dplyr::select()
  