
# packages----
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

#palette seasons---
palette_seasons_4 <- c("winter" = "#002562", 'spring' = "#519741", 'summer' = "#ffb900",'autumn' =  "#96220a")
# palettes taxonomy assigned----
palette_phylums_assigned <- c('Proteobacteria' = "#fcca46","Bacteroidota" = "#669bbc" , 'Actinobacteriota' = "#b0413e",
                              'Cyanobacteria' = "#009e73",'Crenarchaeota' = "#ffa737", 'Verrucomicrobiota' = "#cccccc",
                              'Planctomycetota' = "#69267e", 'Acidobacteriota' = "#1f78b4",'Bdellovibrionota' = "#8c789d",
                              'Firmicutes' = "#637b88", 'Myxococcota'= "#003029", 'Nitrospinota'= "#e3a6ce",
                              'Campilobacterota' = "#002960", 'Deinococcota'= "#ba9864",'Fusobacteriota' ="#fb9a99",
                              'Desulfobacterota' = "#005c69", 'Methylomirabilota' = "#000000" ,
                              'Gemmatimonadota' = "#c55e5c", 'Chloroflexi' = "#00d198", 'Spirochaetota' = "#5cb1d6",
                              'Calditrichota' = "#8a007a", 'Halobacterota' = "#b79c64", 'Nitrospirota' = "#41815f",
                              'Dependentiae' = "#5b95e5", 'Patescibacteria' = "#33af9c",'Cloacimonadota' = "#fbed5c",
                              'Synergistota' = "#ce7800", 'Abditibacteriota' = "#87878b", 'Deferribacterota' = "#4dbaa9")

palette_class_assigned <- c('Gammaproteobacteria' = '#fcca46', 'Alphaproteobacteria' = '#d95726', 'Zetaproteobacteria' = '#EBCD92',
                            'Bacteroidia' = '#669BBC','Rhodothermia' =  '#00425E',
                            'Ignavibacteria' = '#CAFFFF', 'Chlorobia' = '#5BD1FF', 'Kryptonia' = '#0071A8',
                            'Nitrososphaeria' = '#FFA737',
                            'Cyanobacteriia' = '#009E73', 'Vampirivibrionia' = '#00733C',
                            'Acidimicrobiia' = '#B0413E','Actinobacteria' = '#C55E5C',
                            'Coriobacteriia' = '#B44240', 'Thermoleophilia' = '#AB3F3D',
                            'Verrucomicrobiae' = '#CCCCCC', 'Kiritimatiellae' = '#6D6D6D',
                            'Lentisphaeria' = '#424242', 'Omnitrophia' = '#6D6D6D','Chlamydiae' = '#9B9B9B',
                            'Planctomycetes' = '#69267E', 'Phycisphaerae' = '#BE8DCB','Blastocatellia' = '#00ADFF',
                            'Holophagae' = '#86A4C6', 'Vicinamibacteria' = '#002671',
                            'Acidobacteriae' = '#1F78B4', 'Thermoanaerobaculia' = '#0051BF',
                            'Bdellovibrionia' = '#8C789D','Oligoflexia' = '#654584',
                            'Bacilli' = '#637B88',  'Clostridia' = '#384F5B',
                            'Negativicutes' = '#91AAB8', 'Syntrophomonadia' = '#C2DCEA',
                            'Desulfitobacteriia' =  '#0E2732', 'Thermoanaerobacteria' = '#005474',
                            'Myxococcia' = '#2E5A51','Polyangia' = '#003029',
                            'Nitrospinia' = '#e3a6ce', 'Campylobacteria' = '#002960',
                            'Deinococci' = '#ba9864','Fusobacteriia' = '#fb9a99',
                            'Desulfovibrionia' = '#005c69', 'Desulfobulbia' = '#3D518E',
                            'Desulfuromonadia' = '#6D7DBE', 'Desulfobacteria' = '#000036', 'Methylomirabilia' = '#000000',
                            'Gemmatimonadetes' = '#c55e5c',
                            'Anaerolineae' = '#00d198', 'Chloroflexia' = '#009F6A',
                            'Spirochaetia' = '#5CB1D6', 'Leptospirae' = '#005576', 
                            'Calditrichia'  = '#8a007a', 'Halobacteria' = '#b79c64', 'Nitrospiria' = '#41815f',
                            'Babeliae' = '#5b95e5', 'Saccharimonadia' = '#33af9c', 'Cloacimonadia' = '#fbed5c',
                            'Synergistia' = '#ce7800', 'Abditibacteria' = '#87878b', 'Deferribacteres' = '#4dbaa9')


# functions----
#source('../Bloomers/R/community_evenness.R')
#source('../Bloomers/R/get_anomalies.R') #la poso així perquè l'he actualizat i com que no he compilat el paquet potser no funciona actualitzada 
#source('../Bloomers/R/calculate_relative_abundance.R')
#source('../Bloomers/R/blooming_event_summary.R')
source('../../Bloomers/R/find_asv_with_anomalies.R')
source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R')

# Prepare data----
## TREBALLO AMB LA TAULA DE L'ADRIÀ 10 YEARS BASE DE DADES COMPLETA
## Faig una taula per als filtres de 0.2 i de 3 amb totes les mostres sense filtrat asv>5% of all samples

setwd("~/Documentos/Doctorat/BBMO/Taula_BBMO_Adria_sense_filtre/")
bbmo_10y <-readRDS("blphy10years.rds") ##8052 asv amb totes les mostres, no está en % aquesta taula
#str(bbmo_10y)
#bbmo_10y <- prune_samples(sample_sums(bbmo_10y) > 10000, bbmo_10y) in case we want to filter samples by number of reads

bbmo_10y <-
  prune_taxa(taxa_sums(bbmo_10y@otu_table) >0, ##filtrar per les asv que son 0 tot el dataset
             bbmo_10y)

bbmo_10y |>
  nsamples() #237 samples 

bbmo_10y |>
  ntaxa() #8052 but if we filter for those that are 0 during the whole datsaeet then we got 7,849 ASVs

## separate datasets by ASV_tab, taxonomy and metadata
asv_tab_bbmo_10y_l <- bbmo_10y@otu_table |>
  as_tibble()

tax_bbmo_10y <- bbmo_10y@tax_table |> 
  #mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(bbmo_10y@otu_table))) |> #ja té uns nº d'OTUs
  as_tibble()

m_bbmo_10y <- bbmo_10y@sam_data |>
  as_tibble()

## tidy colnames----
asv_tab_bbmo_10y_l |>
  colnames()

tax_bbmo_10y |>
  colnames()

colnames(asv_tab_bbmo_10y_l) <- c('asv_num', "sample_id", 'reads')

colnames(tax_bbmo_10y) <- c("asv_num", "kingdom", "phylum", "class", "order", "family", "genus",
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
                          "year", "month", "day", "season",            
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

##Pivot into long format
# asv_tab_all_perc_filt_02_long <- asv_tab_all_perc_filt_02 |>
#   pivot_longer(cols = starts_with('asv'), values_to = 'relative_abundance', names_to = 'asv_num')
# 
# colnames(asv_tab_all_perc_filt_02_long) <- c('sample_id', 'asv_num','relative_abundance')
# 
# asv_tab_all_perc_filt_3_long <- asv_tab_all_perc_filt_3 |>
#   pivot_longer(cols = starts_with('asv'), values_to = 'relative_abundance', names_to = 'asv_num')
# 
# colnames(asv_tab_all_perc_filt_3_long) <- c('sample_id', 'asv_num','relative_abundance')
  
# Calculate pseudoabundances ----
## This type of normalization makes more sense with the 0.2 fraction since cytometry accounts mainly for this fraction of the biomass
# m_bbmo_10y |>
#   colnames()

asv_tab_10y_3_pseudo <- asv_tab_10y_3_rel |>
calculate_pseudoabund(abund_data = m_bbmo_10y, rel_abund = relative_abundance, 
                      total_abund = bacteria_joint, 
                      by_ = 'sample_id')

asv_tab_10y_02_pseudo <- asv_tab_10y_02_rel |>
  calculate_pseudoabund(abund_data = m_bbmo_10y, 
                        rel_abund = relative_abundance, 
                        total_abund = bacteria_joint, 
                        by_ = 'sample_id')

# Transform data to CLR we use the zCompositions package which helps us to deal with 0 before applying this transformation.----
## transform asv_tab into wider format to go into vegan package
asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_l |>
  pivot_wider(names_from = 'asv_num', values_from = 'reads', values_fill = 0) |>
  as.data.frame()

rownames(asv_tab_bbmo_10y_w) <- asv_tab_bbmo_10y_w$sample_id

# asv_tab_bbmo_10y_w |>
#   dim()

asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_w[,-1]

#geometric mean
gm <- function(x){
  exp(mean(log(x[x>0])))
}

##with this transformation I'm losing samples (due to too much 0 in some samples, z.warning set up t0 0.99 to keep all samples)
### at 0.8 (default) I lose 30 samples which belonged to the years corresponding to harbour remodelation 
zclr_df <- cmultRepl(asv_tab_bbmo_10y_w, method = 'CZM', output = 'p-count', z.warning = 0.99
                     #adjust = 0.2,    t = 237, s = 7849
                  ) |>
  as_tibble(rownames = "sample_id") %>%
  pivot_longer(-sample_id) %>%
  group_by(sample_id) %>%
  dplyr::mutate(zclr = log(value/gm(value))) %>%
  ungroup() %>%
  dplyr::select(-value) %>%
  pivot_wider(names_from = name, values_from = zclr, values_fill = 0) %>%
  column_to_rownames("sample_id")

##we create two datasets one for FL and one for PA
asv_tab_10y_02_zclr <- zclr_df |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'zclr') |>
  dplyr::filter(str_detect(sample_id, '_0.2_'))

asv_tab_10y_3_zclr <- zclr_df |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'zclr') |>
  dplyr::filter(str_detect(sample_id, '_3_'))

##dimesions are different from relative and pseudo dataset and zclr understand why
### i have less ASVs in the zclr dataset specifically 5282 less
#7849-2567

asv_tab_10y_02_pseudo |>
  group_by(asv_num) |>
  summarize(n = n()) |>
  dim() #7849

asv_tab_10y_02_zclr  |>
  group_by(asv_num) |>
  summarize(n = n()) |>
  dim() #2567

asv_tab_10y_3_pseudo |>
  group_by(asv_num) |>
  summarize(n = n()) |>
  dim() #7849

asv_tab_10y_3_zclr  |>
  group_by(asv_num) |>
  summarize(n = n()) |>
  dim() #2567

##creation of a complete dataset with all the normalitzations performed (relative_abundances, pseudoabundances and zclr)----
asv_tab_10y_02_pseudo_zclr <- asv_tab_10y_02_pseudo |>
  left_join(asv_tab_10y_02_zclr, by = c('sample_id', 'asv_num')) |>
  dplyr::mutate(zclr = case_when(is.na(zclr)~ 0,
                                  !is.na(zclr) ~ zclr))

asv_tab_10y_3_pseudo_zclr <- asv_tab_10y_3_pseudo |>
  left_join(asv_tab_10y_3_zclr, by = c('sample_id', 'asv_num')) |>
  dplyr::mutate(zclr = case_when(is.na(zclr)~ 0,
                                 !is.na(zclr) ~ zclr))

nrow(asv_tab_10y_02_pseudo_zclr) == nrow( asv_tab_10y_02_pseudo)
nrow(asv_tab_10y_3_pseudo_zclr) == nrow(asv_tab_10y_3_pseudo)
#if FALSE something is WRONG!

#Do I have ASVs that only have one positive value?
# Check columns with only 1 non-zero in the given data set.
# checkNumZerosCol <- apply(asv_tab_bbmo_10y_w,2,function(x) sum(x==0))
# cases <- which(checkNumZerosCol == (nrow(asv_tab_bbmo_10y_w) - 1)) 
# length(cases) # 1797 columns are all zeros but one positive value
# zcomp <- cmultRepl(asv_tab_bbmo_10y_w[,-cases]) # GBM imputation without them works
## no és aquest el problema perquè em segueixen faltant mostres.

#look for the lost samples why do they disappear? (to much 0)----
# zclr_default <- cmultRepl(asv_tab_bbmo_10y_w, method = 'CZM', output = 'p-count'#, z.warning = 0.99
#                      #adjust = 0.2,    t = 237, s = 7849
# )
# ## missing samples 3
# samples <- asv_tab_10y_3_zclr$sample_id |>
#   unique() |>
#   as_tibble()
# 
# y <- asv_tab_bbmo_10y_w |>
#   row.names() |>
#   as_tibble() |>
#   dplyr::filter(str_detect(value, '_3_'))
# 
# missing_samples_3 <- y |>
#   dplyr::filter(!(value %in% samples$value))
# 
# ##missing samples 02 
# y |>
#   as_tibble() |>
#   dplyr::filter(!(y %in% samples$value))
# 
# samples <- asv_tab_10y_02_zclr$sample_id |>
#   unique() |>
#   as_tibble()
# 
# y <- asv_tab_bbmo_10y_w |>
#   row.names() |>
#   as_tibble() |>
#   dplyr::filter(str_detect(value, '_0.2_'))
# 
# missing_samples_02 <- y |>
#   dplyr::filter(!(value %in% samples$value)) 
# 
# all_missing_samples <- missing_samples_02  |>
#   bind_rows(missing_samples_3)
# 
# asv_tab_missing_samples <- asv_tab_bbmo_10y_w |>
#   rownames_to_column( var= 'sample_id') |>
#   dplyr::filter(sample_id %in% all_missing_samples$value) |>
#   dplyr::select(-sample_id) |>
#   rowSums() |>
#   as_tibble()

## number samples that disappeared when we apply the function cumltRepl with the default parameters
# asv_tab_10y_02_zclr <- zclr_df |>
#   rownames_to_column(var = 'sample_id') |>
#   pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'zclr') |>
#   dplyr::filter(str_detect(sample_id, '_0.2_'))
# 
# asv_tab_10y_02_zclr |>
#   group_by(sample_id) |>
#   dplyr::summarize(n = n()) |>
#   count() #115 
# 
# asv_tab_10y_3_zclr <- zclr_df |>
#   rownames_to_column(var = 'sample_id') |>
#   pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'zclr') |>
#   dplyr::filter(str_detect(sample_id, '_3_'))
# 
# asv_tab_10y_3_zclr |>
#   group_by(sample_id) |>
#   dplyr::summarize(n = n()) |>
#   count() #87


# Calculate diversity parameters ----
## We calculate different diversity parameters to check if the blooming events detected have an effect on the community structure
## Community Evenness----
### We need to apply rarefaction because we have uneven sequencing effort, and this will affect community diversity analysis
### Rarefaction (process that repeats the subsampling many times and tries to overcome the uneven sequencing effort bias)----


#### transform asv_tab into wider format to go into vegan package
asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_l |>
  pivot_wider(names_from = 'asv_num', values_from = 'reads', values_fill = 0) |>
  as.data.frame()

rownames(asv_tab_bbmo_10y_w) <- asv_tab_bbmo_10y_w$sample_id

asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_w[,-1]

#### calculate the minimum read size of all samples to rarefy at that number
min_n_seqs <- asv_tab_bbmo_10y_l |>
  group_by(sample_id) |>
  dplyr::summarize(n_seqs = sum(reads)) |>
  dplyr::summarize(min = min(n_seqs)) |>
  pull(min)
 
##this function gives us a randomly rarefied community data
# rrarefy(asv_tab_bbmo_10y_w, sample = min_n_seqs) |>
#   as_tibble(rownames = 'sample_id') #just rarefying (one simgle subsampling)

## we use rarefaction (which repeats the subsampling step many times)
## perform this a 1000 times to get an empirical diversity values calculating the mean value for each timepoint.
asv_tab_bbmo_10y_w_rar <- rrarefy.perm(asv_tab_bbmo_10y_w, 
                                       sample = min_n_seqs, 
                                       n = 1000, 
                                       round.out = T)

### without rarefaction---
# asv_tab_bbmo_10y_l |> ### we need original reads
#   str()
# 
# community_eveness_02 <- asv_tab_bbmo_10y_l |>
#   mutate(reads = as.numeric(reads)) |>
#   dplyr::filter(str_detect(sample_id, '_0.2_')) |>
#   #dplyr::select(sample_id, reads, asv_num) |>
#   as_tibble() |>
#   group_by(sample_id) |>
#   #ungroup() |>
#   dplyr::reframe(community_eveness_result = community_evenness(abundances = reads, index = "Pielou"))
# 
# community_eveness_3 <- asv_tab_bbmo_10y_l |>
#   mutate(reads = as.numeric(reads)) |>
#   dplyr::filter(str_detect(sample_id, '_3_')) |>
#   #dplyr::select(sample_id, reads, asv_num) |>
#   as_tibble() |>
#   group_by(sample_id) |>
#   #ungroup() |>
#   dplyr::reframe(community_eveness_result = community_evenness(abundances = reads, index = "Pielou"))
# 
# asv_tab_10y_filt_3_pseudo <- asv_tab_10y_filt_3_pseudo  |>
#   dplyr::mutate(relative_abundance = as.numeric(relative_abundance))
# 
# asv_tab_10y_filt_02_pseudo <- asv_tab_10y_filt_02_pseudo  |>
#   dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

## Rarefied dataset to calculate Community eveness----
source('../../Bloomers/R/community_evenness.R')
community_eveness_02 <- asv_tab_bbmo_10y_w_rar |>
  as.data.frame() |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou'))

community_eveness_3 <- asv_tab_bbmo_10y_w_rar |>
  as.data.frame() |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou'))

###plot community evenness
community_eveness_02 |>
  bind_rows(community_eveness_3) |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, community_eveness_rar))+
  geom_point()+
  facet_grid(vars(fraction))+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Community Eveness')+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank())

## Bray Curtis dissimilarity----
### 0 means the two sites have the same composition (that is they share all the species), and 1 means the two sites 
### do not share any species.
source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R') #we need to upload it since we updated the function but I didn't compile the package
# 
# bray_curtis_02 <- dissimilarity_matrix(data = asv_tab_10y_02_rel, 
#                                            sample_id_col = sample_id,
#                                            values_cols_prefix = 'BL')
# 
# bray_curtis_3 <- dissimilarity_matrix(data = asv_tab_10y_3_rel, 
#                                           sample_id_col = sample_id,
#                                           values_cols_prefix = 'BL')

### We need the rarefied table transformed to relative abundances
asv_tab_bbmo_10y_l_rel_rar <- asv_tab_bbmo_10y_w_rar |>
  as.data.frame() |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), values_to = 'reads', names_to = 'asv_num') |>
  dplyr::mutate(reads = as.numeric(reads)) |>
  calculate_rel_abund(group_cols = sample_id) 

asv_tab_10y_3_rel_rar <- asv_tab_bbmo_10y_l_rel_rar |>
    dplyr::filter(sample_id %in% m_3$sample_id)

asv_tab_10y_filt_02_rel_rar <- asv_tab_bbmo_10y_l_rel_rar %>%
  dplyr::filter(sample_id %in% m_02$sample_id)

bray_curtis_02_rar <- dissimilarity_matrix(data = asv_tab_10y_filt_02_rel_rar, 
                                    sample_id_col = sample_id,
                                    values_cols_prefix = 'BL')

bray_curtis_3_rar <- dissimilarity_matrix(data = asv_tab_10y_3_rel_rar, 
                                       sample_id_col = sample_id,
                                       values_cols_prefix = 'BL')

###plot bray curtis dissimilarity
#### compare rarified dataset with non-rarified dataset (NO CHANGES, WHY?)
# bray_curtis_02 |>
#   bind_rows(bray_curtis_3) |>
#   dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
#   left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   ggplot(aes(date, bray_curtis_result))+
#   geom_point()+
#   facet_grid(vars(fraction))+
#   scale_x_datetime()+
#   labs(x = 'Time', y = 'Bray Curtis Dissimilarity')+
#   geom_line()+
#   theme_bw()+
#   theme(panel.grid = element_blank(), strip.background = element_blank())

bray_curtis_02_rar |>
  bind_rows(bray_curtis_3_rar) |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, bray_curtis_result))+
  geom_point()+
  facet_grid(vars(fraction))+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity')+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank())


# Discover anomalies----
## For each ASVs based on relative abundances and pseudoabundances-----
### those ASVs that are present > 50% of the sampless
asv_tab_10y_02_pseudo %$%
  sample_id |>
  unique() |>
  summary() #60 is half of the dataset

asv_tab_10y_02_zclr |>
  colnames()

##percentage of ASV present in the dataset that we want to subset by
x <- 120*0.75

z_02 <- asv_tab_10y_02_pseudo |>
  inner_join(asv_tab_10y_02_zclr, by = c('sample_id', 'asv_num')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  group_by(asv_num) |>
  dplyr::reframe(anomalies_ab = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = zclr, plotting = FALSE)[c(1,2,3)])

asv_tab_10y_3_pseudo %$%
  sample_id |>
  unique() |>
  summary() #

x <- 117*0.75
z_3 <- asv_tab_10y_3_pseudo |>
  inner_join(asv_tab_10y_3_zclr, by = c('sample_id', 'asv_num')) |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::filter(num_0 <= x) |>  #filter those ASVs that are 0 more than 50% of the dataset
  #as_tibble() |>
  group_by(asv_num) |>
  dplyr::reframe(anomalies_ab = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
    anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
    anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = zclr, plotting = FALSE)[c(1,2,3)])

## At the level of community, we use the Eveness result and Bray Curtis dissimilarity ----
z_diversity <- bray_curtis_02_rar |>
  dplyr::right_join(community_eveness_02, by = join_by("samples" == "sample_id")) |> 
  #ungroup() %>%
  #group_by(sample_id) %>%
  dplyr::reframe(anomalies_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = TRUE)[c(1,2,3)],# ),
                 anomalies_eveness = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = TRUE)[c(1,2,3)])

z_diversity %>%
  str()

# Filter the ASV_tab by only those ASVs that have an anomaly at some point of the dataset ----
##function
find_asv_with_anomalies <- function(anomalies_result, anomaly_in1, 
                                    anomaly_in2, anomaly_in3 = NULL, 
                                    logic1 = TRUE, 
                                    logic2 = TRUE,
                                    logic3 = NULL, 
                                    asv_col = asv_num) {
  # if(is.list(anomalies_result) == FALSE){
  #   stop("Function stopped: anomalies_result needs to be a list form the get_anomalies function")
  # }
  # if(is.logical({{anomaly_in1}}) == FALSE){
  #   stop("Function stopped: anomaly_in1 needs to be logical (TRUE/FALSE)")
  # }
  
  asv_potential_bloomers <-
    anomalies_result |>
    dplyr::filter(if (!is.null(logic1)) {{anomaly_in1}} %in% logic1 else TRUE) |>
    dplyr::filter(if (!is.null(logic2)) {{anomaly_in2}} %in% logic2 else TRUE) |>
    dplyr::filter(if (!is.null(logic3)) {{anomaly_in3}} %in% logic3 else TRUE) |>
    dplyr::select({{asv_col}}) |>
    as_vector()
  
  return(asv_potential_bloomers)
}

asv_anom_02 <- find_asv_with_anomalies(anomalies_result = z_02, anomaly_in1 = anomalies_ra, 
                                       anomaly_in2 = anomalies_ab,
                                       anomaly_in3 = anomalies_clr, 
                                    logic1 = 'TRUE', logic2 = 'TRUE', 
                                    logic3 = 'TRUE',
                                    asv_col = asv_num)
##86 ASVs cumplint les 3 condicions i sense CLR també

##no acabo d'entendre perquè torna un vector amb numero consecutiu d'ASVs...

asv_anom_3 <- find_asv_with_anomalies(anomalies_result = z_3, anomaly_in1 = anomalies_ra, 
                                      anomaly_in2 = anomalies_clr, 
                                      anomaly_in3 = NULL, #anomalies_ab
                                      logic1 = 'TRUE', logic2 = 'TRUE', 
                                      logic3 = NULL, ##per 3 no és representatiu la pseudoabund
                                      asv_col = asv_num)

##46 ASVS
## número d'ASVs comuns i unics a cada fracció
asv_anom_3_tb <- asv_anom_3 |>
  as_tibble()

asv_anom_02_tb <- asv_anom_02 |>
  as_tibble()

asv_anom_3_tb |>
  bind_rows(asv_anom_02_tb) |>
  unique() ##96 totals >50% del dataset

asv_anom_3_tb |>
  anti_join(asv_anom_02_tb) #10 ASV only in 3

asv_anom_02_tb |>
  anti_join(asv_anom_3_tb) #40 ASVs only in 0.2

asv_tab_10y_02_pseudo_zclr_bloo <-  asv_tab_10y_02_pseudo_zclr |>
  group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(zclr, relative_abundance, pseudoabundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  dplyr::filter(asv_num %in% asv_anom_02 |
           asv_num %in% asv_anom_3) ##recover ASVs that presented anomalies in 02 or 3

asv_tab_10y_02_pseudo_zclr_bloo |>
  group_by(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n())

asv_tab_10y_3_pseudo_zclr_bloo <-  asv_tab_10y_3_pseudo_zclr |>
  group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(zclr, relative_abundance, pseudoabundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  dplyr::filter(asv_num %in% asv_anom_3 |
           asv_num %in% asv_anom_02) ##recover ASVs that presented anomalies in 02 or 3

asv_tab_10y_3_pseudo_zclr_bloo |>
  group_by(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n())

asv_tab_all_bloo <- asv_tab_10y_02_pseudo_zclr_bloo |>
  bind_rows(asv_tab_10y_3_pseudo_zclr_bloo)

asv_tab_all_bloo |>
  group_by(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n()) #17 ASVs

# asv_tab_all_perc_filt_3_long_filt |>
#   group_by(sample_id) |>
#   dplyr::summarize(sum = sum(relative_abundance))

asv_tab_all_bloo |> 
  dplyr::filter(abundance_type == 'relative_abundance') %$%
  abundance_value |>
  range() #max is 0.52

## Recover z-score for each ASV ----
### I want to highlight anomalies for each ASV to do so I recover z-scores for those ASVs that that have high z-scores
### at some point of the dataset. Easy to observe if those ASVs are having random anomalies or all of them happen at the same time
x <- 120*0.75
z_scores_02 <- asv_tab_10y_02_pseudo_zclr |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 50% of the samples
  #dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
  group_by(asv_num) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')


x <- 117*0.75
z_scores_3 <- asv_tab_10y_3_pseudo_zclr |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 < 58) |> ##only anomalies for ASVs that are present in > 50% of the samples
  dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 25% of the samples
  group_by(asv_num) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = 3, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_all <- z_scores_02 |>
  bind_rows(z_scores_3)

##Dataset with all information ----
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') 
  
### since pseudoabundance is only useful for the 0.2 fraction we now only plot the relative abundance changes

##highligh harbour restoration period (NO funciona encara)
harbour_restoration <-  asv_tab_all_bloo |> 
  group_by(sample_id) |>
  subset(date = between(date, '2010-03-01','2012-04-01')) |>
  dplyr::summarize(xmin=min(date),xmax=max(date))

## color by class
asv_tab_all_bloo |>
  #left_join(m_3, by = 'sample_id') |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'zclr')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  # geom_rect(aes(xmin = '2010-03-01 01::00:00', xmax = '2012-04-01 01::00:00', ymin = 0, ymax = Inf),
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  # geom_rect(data = asv_tab_all_perc_filt_3_long_filt, mapping=aes(xmin = date, xmax = date, x=NULL, y=NULL,
  #                                     ymin = -Inf, ymax = Inf, fill = '#D0CFC8'), alpha = 0.4)+
  geom_line(aes(group = asv_num, color = class))+
  geom_point(aes(color = class))+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Class')+
  scale_color_manual(values = palette_class_assigned)+
  facet_grid(fraction~abundance_type, scales = 'free')+
  #facet_grid(fraction~class)+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

asv_tab_all_bloo |>
  #left_join(m_3, by = 'sample_id') |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'relative_abundance')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  # geom_rect(aes(xmin = '2010-03-01 01::00:00', xmax = '2012-04-01 01::00:00', ymin = 0, ymax = Inf),
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  # geom_rect(data = asv_tab_all_perc_filt_3_long_filt, mapping=aes(xmin = date, xmax = date, x=NULL, y=NULL,
  #                                     ymin = -Inf, ymax = Inf, fill = '#D0CFC8'), alpha = 0.4)+
  geom_line(aes(group = asv_num, color = class))+
  geom_point(aes(color = class))+
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'z CLR', color = 'Class')+
  scale_color_manual(values = palette_class_assigned)+
  facet_grid(fraction~abundance_type~asv_num, scales = 'free')+
  #facet_grid(fraction~class)+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

##shape anomaly
asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type != 'pseudoabundance') |>
  dplyr::filter(abundance_type != 'zclr') |>
  #ungroup() |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value, shape = ifelse(z_score_ra >= 1.96, '8', '19')))+ #, color = 'Class' 
  scale_x_datetime()+
  facet_wrap(vars(class), scales = 'free')+
  geom_line(aes(group = fraction, color = class))+
  geom_point(aes(color = class))+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Class')+
  scale_color_manual(values = palette_class_assigned)+
  facet_grid(vars(asv_num), scales = 'free')+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

##plot zclr scores
asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type != 'pseudoabundance') |>
  dplyr::filter(abundance_type != 'relative_abundance') |>
  #ungroup() |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value, shape = ifelse(z_score_ra >= 1.96, '8', '19')))+ #, color = 'Class' 
  scale_x_datetime()+
  facet_wrap(vars(class), scales = 'free')+
  geom_line(aes(group = fraction, color = class))+
  geom_point(aes(color = class))+
  ##geom_vline(xintercept = c('2006-01-01', '2005-01-01'))+ no va!!
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'CLR', color = 'Class')+
  scale_color_manual(values = palette_class_assigned)+
  facet_grid(vars(asv_num), scales = 'free')+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

##indentify those that are seasonal from those that are not seasonal
asv_tab_all_bloo |>
  colnames()

asv_tab_all_bloo_z_tax$season <- asv_tab_all_bloo_z_tax$season |>
  factor(levels = c('winter', 'spring', 'summer', 'autumn'))

asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type != 'pseudoabundance') |>
  dplyr::filter(abundance_type != 'zclr') |>
  ggplot(aes(day_of_year, abundance_value, color = season))+
  geom_line(aes(group = year))+
  geom_point(aes(color = season))+
  scale_color_manual(values = palette_seasons_4)+
  facet_grid(asv_num~fraction, scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance(%)', color = 'Season')+
  scale_y_continuous(labels = percent_format())+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())
  
##color anomaly----
### I would like to recover information from 0.2 and 3 for those ASVs that presented anomalies along the dataset.
x <- asv_anom_3 |>
  as_tibble()

y <- asv_anom_02 |>
  as_tibble()

asv_anom_all <- x |>
  bind_rows(y) |>
  as_tibble() |>
  unique() |>
  as_vector()

##plot those by anomaly presence colored in red
asv_tab_all_bloo |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
                                   z_score_ra == 'NaN' ~ 0,
                                   z_score_ra == Inf ~ 0,
                                   TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly_color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  #ungroup() |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  # ggplot(aes(date, abundance_value, color = ifelse(is.na(z_score_ra), "#080808", 
  #            if_else(z_score_ra >= 1.96, '#9F0011', '#080808'))))+ #, color = 'Class' , shape = class 
  #ggplot(aes(date, abundance_value, color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')))+
  ggplot(aes(date, abundance_value, color = anomaly_color))+  
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  geom_line(aes(group = asv_num), color = '#080808')+ #, color = '#3D3B3B'
  geom_point(alpha = 1)+ #aes(shape = class)
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Anomaly')+ #, shpae = 'Class'
      scale_color_identity()+
  #scale_color_manual(values = if_else(asv_tab_z_scores_all$z_score_ra >= 1.96,  '#9F0011', '#080808', missing = '#080808'))+
  facet_wrap(fraction~abundance_type, scales = 'free')+
  #facet_grid(fraction~class, scales = 'free')+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

##plot all month together (seasonal anomalies?)
asv_tab_all_bloo |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
                                   z_score_ra == 'NaN' ~ 0,
                                   z_score_ra == Inf ~ 0,
                                   TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly_color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  dplyr::filter(z_score_ra_ed >= 1.96) |>
  #ungroup() |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  # ggplot(aes(date, abundance_value, color = ifelse(is.na(z_score_ra), "#080808", 
  #            if_else(z_score_ra >= 1.96, '#9F0011', '#080808'))))+ #, color = 'Class' , shape = class 
  #ggplot(aes(date, abundance_value, color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')))+
  ggplot(aes(month, abundance_value))+  
  geom_violin(aes(month, abundance_value))+
  #scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  #geom_smooth(aes(month, abundance_value))+
  #geom_line(aes(group = asv_num), color = '#080808')+ #, color = '#3D3B3B'
  geom_point(aes(color = class), alpha = 1)+ #aes(shape = class)
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Month', y = 'Relative abundance (%)', color = 'Class')+ #, shpae = 'Class'
  #scale_color_identity()+
  scale_color_manual(values = palette_class_assigned)+
  #scale_color_manual(values = if_else(asv_tab_z_scores_all$z_score_ra >= 1.96,  '#9F0011', '#080808', missing = '#080808'))+
  facet_wrap(fraction~abundance_type, scales = 'free')+
  #facet_grid(fraction~class, scales = 'free')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

##plot all seasons together (seasonal anomalies?)
asv_tab_all_bloo |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
                                   z_score_ra == 'NaN' ~ 0,
                                   z_score_ra == Inf ~ 0,
                                   TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly_color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  dplyr::filter(z_score_ra_ed >= 1.96) |>
  #ungroup() |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  # ggplot(aes(date, abundance_value, color = ifelse(is.na(z_score_ra), "#080808", 
  #            if_else(z_score_ra >= 1.96, '#9F0011', '#080808'))))+ #, color = 'Class' , shape = class 
  #ggplot(aes(date, abundance_value, color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')))+
  ggplot(aes(season, abundance_value))+  
  geom_violin(aes(season, abundance_value), draw_quantiles = T)+
  #scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  #geom_smooth(aes(month, abundance_value))+
  #geom_line(aes(group = asv_num), color = '#080808')+ #, color = '#3D3B3B'
  geom_point(aes(color = class), alpha = 1)+ #aes(shape = class)
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Season', y = 'Relative abundance (%)', color = 'Class')+ #, shpae = 'Class'
  #scale_color_identity()+
  scale_color_manual(values = palette_class_assigned)+
  #scale_color_manual(values = if_else(asv_tab_z_scores_all$z_score_ra >= 1.96,  '#9F0011', '#080808', missing = '#080808'))+
  facet_wrap(fraction~abundance_type, scales = 'free')+
  #facet_grid(fraction~class, scales = 'free')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())
  
#pivot_wider(id_cols = sample_id_num, values_from = anomalies_ra, names_from = asv_num)

##ocurrence of this AVSs vs frequency of blooming events and magnitude----
###calculate occurency
nsamples_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3') %$%
  sample_id |>
  unique() |>
  length()

occurence_perc_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  dplyr::filter(abundance_value > 0) |>
  group_by(asv_num, sample_id) |>
  summarize(n = n()) |>
  group_by(asv_num) |>
  summarize(occurence = sum(n)) |>
  mutate(occurence_perc = occurence/nsamples_3,
         fraction = '3')

nsamples_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') %$%
  sample_id |>
  unique() |>
  length()

occurence_perc_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  dplyr::filter(abundance_value > 0) |>
  group_by(asv_num, sample_id) |>
  summarize(n = n()) |>
  group_by(asv_num) |>
  summarize(occurence = sum(n)) |>
  mutate(occurence_perc = occurence/nsamples_02,
         fraction = '0.2')

occurence_perc <- occurence_perc_02 |>
  bind_rows(occurence_perc_3)

###calculate number of anomalies > 0.1% in the dataset
#### we use number of total samples or number of presence of that ASV?
anom_perc_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  dplyr::filter(abundance_value >= 0.1) |>
  dplyr::mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
                                   z_score_ra == 'NaN' ~ 0,
                                   z_score_ra == Inf ~ 0,
                                   TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly = case_when(z_score_ra >= 1.96 ~ 1,
                                          z_score_ra < 1.96 ~ 0)) |>
  group_by(asv_num) |>
  dplyr::summarize(n_anom = sum(anomaly)) |>
  dplyr::mutate(anom_perc = n_anom/nsamples_3,
                fraction = '3')

anom_perc_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  dplyr::filter(abundance_value >= 0.1) |>
  dplyr::mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
                                          z_score_ra == 'NaN' ~ 0,
                                          z_score_ra == Inf ~ 0,
                                          TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly = case_when(z_score_ra >= 1.96 ~ 1,
                                    z_score_ra < 1.96 ~ 0)) |>
  group_by(asv_num) |>
  dplyr::summarize(n_anom = sum(anomaly)) |>
  dplyr::mutate(anom_perc = n_anom/nsamples_02,
                fraction = '0.2')

anom_perc <- anom_perc_3 |>
  bind_rows(anom_perc_02)

###max relative abundance (canviar a max in a blooming envent???)
max_rel_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  group_by(asv_num) |>
  dplyr::summarize(max_rel = max(abundance_value)) |>
  dplyr::mutate(fraction = '3')

max_rel_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  group_by(asv_num) |>
  dplyr::summarize(max_rel = max(abundance_value))|>
  dplyr::mutate(fraction = '0.2')

max_rel <- max_rel_3 |>
  bind_rows(max_rel_02)
  
###plot Occurence vs 
bloom_occurrence_tax <- occurence_perc |>
  left_join(anom_perc, by = c('asv_num','fraction')) |>
  left_join(max_rel, by = c('asv_num','fraction')) |>
  left_join(tax_bbmo_10y, by = 'asv_num') 

occurence_perc |>
 left_join(anom_perc, by = c('asv_num','fraction')) |>
  left_join(max_rel, by = c('asv_num','fraction')) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  ggplot(aes(occurence_perc, anom_perc, color = class, size = max_rel*10, shape = fraction))+
  geom_point()+
  scale_color_manual(values = palette_class_assigned)+
  labs(x = 'ASV Occurrence (%)', y = 'Anomaly (%)', color = 'Class')+
  scale_x_continuous(labels = percent_format())+
  scale_y_continuous(labels = percent_format())+
  #facet_wrap(vars(asv_num))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'right', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())
  

# Analysis of FL and PA blooming events ----
## Do blooming events synchronize in PA and FL? 
## Some ASVs prefer to bloom in one fraction than another? Is there a preference?
## plot by fractions
asv_tab_all_bloo$fraction <- asv_tab_all_bloo$fraction  |> 
  factor(levels =  (c('0.2', '3')))

asv_tab_all_bloo |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'zclr')) |>
  #dplyr::filter(Class != is.na(Class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value, color = 'fraction', group = 'fraction'))+
  scale_x_datetime()+
  #facet_grid(vars(Class), scales = 'free')+
  geom_line(aes(group = fraction, color = fraction))+
  geom_point(aes(color = fraction), size = 0.7)+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Abundance', color = 'Fraction')+
  #scale_color_manual(values = palette_class_assigned)+
  facet_wrap(vars(asv_num), scales = 'free')+
  guides(fill = 'none')+
  scale_color_manual(values = c('#9056a9' , '#6FA956'))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'right', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

asv_tab_all_bloo |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'relative_abundance')) |>
  #dplyr::filter(Class != is.na(Class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value, color = 'fraction', group = 'fraction'))+
  scale_x_datetime()+
  #facet_grid(vars(Class), scales = 'free')+
  geom_line(aes(group = fraction, color = fraction))+
  geom_point(aes(color = fraction), size = 0.7)+
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'CLR', color = 'Fraction')+
  #scale_color_manual(values = palette_class_assigned)+
  facet_wrap(vars(asv_num), scales = 'free')+
  guides(fill = 'none')+
  scale_color_manual(values = c('#9056a9' , '#6FA956'))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'right', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

# 
#   ggsave('bbmo_10y_bloomers_10perc.pdf', bbmo_10y_bloomers_10perc,
#        path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/",
#        width = 230,
#        height = 180,
#        units = 'mm')


# Analysis with years sequenced with Parada primers (2014-2015) -----

## import data
BBMO_parada <- read_xlsx('../BBMO_bloomers/data/2014-2015_Primers_Parada/ASVtab_BBMO_ParadaPrimers_2014to2016_CORRECTED_v1.xlsx') |>
  as_tibble()

BBMO_parada |>
  dim()

BBMO_parada <- BBMO_parada |>
  dplyr::mutate(asv_num = str_c( 'asv' , 1:n())) |> #create a new column with asv num
  #dplyr::mutate(asv_num = str_c( 'asv' , 1:row(BBMO_parada))) |> 
  dplyr::select(-contains('mock'))
  
BBMO_parada  %$%
  asv_num |>
  unique()

BBMO_parada |>
  colnames()

BBMO_parada_tax <- BBMO_parada |>
  dplyr::select(seq, Kingdom, Phylum, Class, Order, Family, Genus, asv_num)

BBMO_parada_filt <- BBMO_parada |>
  dplyr::select(-c(seq, Kingdom, Phylum, Class, Order, Family, Genus, reads)) ##to be done check even sequencing for all samples

BBMO_parada_filt |>
  colnames()

BBMO_parada_rel <- BBMO_parada_filt |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'reads') |>
  calculate_rel_abund(group_cols = sample_id) 

BBMO_parada_rel_02 <- BBMO_parada_rel |>
  dplyr::filter(str_detect(sample_id, '-02_'))

BBMO_parada_rel_3 <- BBMO_parada_rel |>
  dplyr::filter(str_detect(sample_id, '-3_'))

### creatae metadata from sample codes
metadata_parada <- BBMO_parada_rel |>
  ungroup() |>
  dplyr::select(sample_id) |>
  dplyr::distinct(sample_id) |>
  tidyr::separate(col = sample_id, into = c('station', 'year_month', 'mesh', 'sample'), 
                  remove = FALSE)

##detect anomalies   
  z_02 <- BBMO_parada_rel_02 |>
  #group_by(asv_num) |>
  # dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  # dplyr::filter(num_0 < 100) |>
  #as_tibble() |>
  group_by(asv_num) |>
  dplyr::reframe(#anomalies_ab = get_anomalies(time_lag = 3, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_ra = get_anomalies(time_lag = 3, negative = FALSE, 
                                              cutoff = 1.96, na_rm = TRUE, 
                                              values = relative_abundance, 
                                              plotting = FALSE)[c(1,2,3)])
  
  z_3 <- BBMO_parada_rel_3 |>
    #group_by(asv_num) |>
    # dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
    # dplyr::filter(num_0 < 100) |>
    #as_tibble() |>
    group_by(asv_num) |>
    dplyr::reframe(#anomalies_ab = get_anomalies(time_lag = 3, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
      anomalies_ra = get_anomalies(time_lag = 3, negative = FALSE, 
                                   cutoff = 1.96, na_rm = TRUE, 
                                   values = relative_abundance, 
                                   plotting = FALSE)[c(1,2,3)])
  
## Filter the ASV_tab by only those ASVs that have an anomaly at some point of the dataset----

  
  asv_anom_02 <- find_asv_with_anomalies(anomalies_result = z_02, anomaly_in1 = anomalies_ra, anomaly_in2 = NULL, 
                                         logic1 = 'TRUE', logic2 = NULL, 
                                         asv_col = asv_num)
  
  asv_tab_all_perc_filt_02_long_filt <-  BBMO_parada_rel |>
    group_by(asv_num) |>
    dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format percentatge!
    #pivot_longer(cols = c(relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
    dplyr::filter(asv_num %in% asv_anom_02)
  
  asv_tab_all_perc_filt_02_long_filt %$%
    asv_num |>
    unique()
  
  asv_tab_all_perc_filt_02_long_filt %$%
    sample_id |>
    unique()
  
  asv_tab_all_perc_filt_02_long_filt |>
    left_join(metadata_parada, by = 'sample_id') |>
    left_join(BBMO_parada_tax, by = 'asv_num') |>
    ggplot(aes(year_month, relative_abundance, color = Class))+
    geom_point()+
    geom_line(aes(group = asv_num))+
    facet_grid(vars(mesh))+
    scale_y_continuous(labels = percent_format())+
    scale_color_manual(values = palette_class_assigned)+
    labs(x = 'Year-month', y = 'Relative abundance (%)')+
    theme_bw()+
    theme(panel.grid = element_blank(), strip.background = element_blank())
  
  asv_anom_3 <- find_asv_with_anomalies(anomalies_result = z_3, anomaly_in1 = anomalies_ra, anomaly_in2 = NULL, 
                                         logic1 = 'TRUE', logic2 = NULL, 
                                         asv_col = asv_num)
  
  asv_tab_all_perc_filt_3_long_filt <-  BBMO_parada_rel |>
    group_by(asv_num) |>
    dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format percentatge!
    #pivot_longer(cols = c(relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
    dplyr::filter(asv_num %in% asv_anom_3)
  
  asv_tab_all_perc_filt_3_long_filt %$%
    asv_num |>
    unique()
  
  asv_tab_all_perc_filt_3_long_filt %$%
    sample_id |>
    unique()
  
  asv_tab_all_perc_filt_3_long_filt |>
    left_join(metadata_parada, by = 'sample_id') |>
    left_join(BBMO_parada_tax, by = 'asv_num') |>
    ggplot(aes(year_month, relative_abundance, color = Class))+
    geom_point()+
    geom_line(aes(group = asv_num))+
    facet_grid(vars(mesh))+
    scale_y_continuous(labels = percent_format())+
    scale_color_manual(values = palette_class_assigned)+
    labs(x = 'Year-month', y = 'Relative abundance (%)')+
    theme_bw()+
    theme(panel.grid = element_blank(), strip.background = element_blank())
  
  asv_tab_all_perc_filt_3_long_filt %$%
    asv_num |>
    unique()
  
  asv_tab_all_perc_filt_02_long_filt %$%
    asv_num |>
    unique()
  

  
  
## Testing CLR transformation and other transformations to calculate distances----

group_count <- asv_tab_bbmo_10y_l |>
    group_by(sample_id) |>
    dplyr::summarize(n_seqs = sum(reads))
  
group_count$n_seqs |>
  range()

asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_l |>
  pivot_wider(names_from = asv_num, values_from = reads, values_fill = 0) |>
  as.data.frame()

rownames(asv_tab_bbmo_10y_w) <- asv_tab_bbmo_10y_w$sample_id

asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_w[,-1]

norare_dist <- vegdist(asv_tab_bbmo_10y_w, 
                       method = 'euclidean')

rare_dist <- avgdist(asv_tab_bbmo_10y_w, 
                     dmethod = 'euclidean', 
                     sample = min(group_count$n_seqs))
#geometric mean
gm <- function(x){
  exp(mean(log(x[x>0])))
}
x <- c(3,8,9,0)

gm(x)

#what happens with 0? Robust CLR removes 0 from the transformations
asv_tab_bbmo_10y_rclr <- 
  asv_tab_bbmo_10y_l |>
  dplyr::mutate(reads = as.numeric(reads)) |>
   group_by(sample_id) |>
  #dplyr::mutate(gm = gm(reads))
   dplyr::mutate(rclr = log(reads/gm(reads))) |>
   ungroup() |>
   select(-reads) |>
   pivot_wider(names_from = asv_num, values_from = rclr, values_fill = 0) |>
   as.data.frame()
   
rownames(asv_tab_bbmo_10y_rclr) <- asv_tab_bbmo_10y_rclr$sample_id
asv_tab_bbmo_10y_rclr <- asv_tab_bbmo_10y_rclr[,-1]

rclr_dist <- vegdist(asv_tab_bbmo_10y_rclr, method = 'euclidean')

 #install.packages("zCompositions") 
 
library(zCompositions) #helps to deal with 0 in the dataset

 zclr_df <- cmultRepl(asv_tab_bbmo_10y_w, method = 'CZM', output = 'p-count') |>
   as_tibble(rownames = "sample_id") %>%
   pivot_longer(-samples_id) %>%
   group_by(samples_id) %>%
   mutate(zclr = log(value/gm(value))) %>%
   ungroup() %>%
   select(-value) %>%
   pivot_wider(names_from=name, values_from=zclr, values_fill=0) %>%
   column_to_rownames("samples")

zclr_dist <- vegdist(zclr_df, method = 'euclidean')

nonrare_tbl <- norare_dist |>
  as.matrix() |>
  as_tibble(rownames = 'sample_id') |>
  pivot_longer(cols = -sample_id) |>
  filter(name < sample_id)

rare_tbl <- rare_dist |>
  as.matrix() |>
  as_tibble(rownames = 'sample_id') |>
  pivot_longer(cols = -sample_id) |>
  filter(name < sample_id)

#rclr_dist
zclr_tbl <- zclr_dist |>
  as.matrix() |>
  as_tibble(rownames = 'sample_id') |>
  pivot_longer(cols = -sample_id) |>
  filter(name < sample_id)


 inner_join(nonrare_tbl, rare_tbl, by = c('sample_id', 'name')) |>
   inner_join(zclr_tbl, by = c('sample_id', 'name')) |>
   inner_join(group_count, by = 'sample_id') |>
   inner_join(group_count, by = c('name' = 'sample_id')) |>
   dplyr::mutate(diffs = abs(n_seqs.x - n_seqs.y)) |>
   dplyr::select(sample_id, name, norare = value.x, rare=value.y, zclr=value, diffs) |>
   dplyr::filter(str_detect(sample_id, '_0.2_') &
                   str_detect(name, '_0.2_')) |>
   pivot_longer(cols = c(norare, rare, zclr), names_to = 'method', values_to = 'dist') |>
   ggplot(aes(x=diffs, y = dist))+
   geom_point()+
   facet_wrap(~method, scales = 'free_y')+
   geom_smooth()

 
 inner_join(nonrare_tbl, rare_tbl, by = c('sample_id', 'name')) |>
   inner_join(zclr_tbl, by = c('sample_id', 'name')) |>
   inner_join(group_count, by = 'sample_id') |>
   inner_join(group_count, by = c('name' = 'sample_id')) |>
   dplyr::mutate(diffs = abs(n_seqs.x - n_seqs.y)) |>
   dplyr::select(sample_id, name, norare = value.x, rare=value.y, zclr=value, diffs) |>
   dplyr::filter(str_detect(sample_id, '_3_') &
                   str_detect(name, '_3_')) |>
   pivot_longer(cols = c(norare, rare, zclr), names_to = 'method', values_to = 'dist') |>
   ggplot(aes(x=diffs, y = dist))+
   geom_point()+
   facet_wrap(~method, scales = 'free_y')+
   geom_smooth()

 