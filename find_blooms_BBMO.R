
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
library(ggplot2)
library(scales)
library(zCompositions)

#labels----
labs_fraction <- as_labeller(c('0.2' = 'Free living (0.2 um)',
                               '3' = 'Particle attached (3 um)'))

labs_diversity <- as_labeller(c('community_eveness_rar' = 'Community Eveness', 
                                'bray_curtis_result' = 'Bray-Curtis dissimilarity'))
# palettes----
palette_seasons_4 <- c("winter" = "#002562", 'spring' = "#519741", 'summer' = "#ffb900",'autumn' =  "#96220a")

palette_fraction <- c('0.2' = '#00808F', '3' = '#454545')

## palettes taxonomy assigned----
palette_phylums_assigned <- c('Proteobacteria' = "#fcca46","Bacteroidota" = "#669bbc" , 'Actinobacteriota' = "#b0413e",
                              'Cyanobacteria' = "#009e73",'Crenarchaeota' = "#ffa737", 'Verrucomicrobiota' = "#0051BF",
                              'Planctomycetota' = "#69267e", 'Acidobacteriota' = "#1f78b4",'Bdellovibrionota' = "#8c789d",
                              'Firmicutes' = "#637b88", 'Myxococcota'= "#003029", 'Nitrospinota'= "#e3a6ce",
                              'Campilobacterota' = "#002960", 'Deinococcota'= "#ba9864",'Fusobacteriota' ="#fb9a99",
                              'Desulfobacterota' = "#005c69", 'Methylomirabilota' = "#000000" ,
                              'Gemmatimonadota' = "#c55e5c", 'Chloroflexi' = "#00d198", 'Spirochaetota' = "#5cb1d6",
                              'Calditrichota' = "#8a007a", 'Halobacterota' = "#b79c64", 'Nitrospirota' = "#41815f",
                              'Dependentiae' = "#5b95e5", 'Patescibacteria' = "#33af9c",'Cloacimonadota' = "#fbed5c",
                              'Synergistota' = "#ce7800", 'Abditibacteriota' = "#87878b", 'Deferribacterota' = "#4dbaa9")

palette_phylums_assigned_bloo <- c('Proteobacteria' = "#fcca46", 'Cyanobacteria' = "#009e73", 'Verrucomicrobiota' = "#cccccc",
                                   'Verrucomicrobiota' = "#0051BF", "Bacteroidetes" = "#669bbc")

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


palette_class_assigned_bloo <- c('Gammaproteobacteria' = '#FFA737', 'Alphaproteobacteria' = '#B0413E', 
                                'Verrucomicrobiae' = '#005c69', 'Opitutae' =   '#fcca46',  'Phycisphaerae' =  '#e3a6ce',   
                                 'Flavobacteriia' =  '#0051BF',  'Sphingobacteriia'  = '#69267E' ,
                                'Cyanobacteria'  =  '#009F6A',  'Deltaproteobacteria' = '#000000')

bloomers_tax |>
  dplyr::filter(taxonomy_rank == 'sum_f') |>
  left_join(tax_bbmo_10y, by = c('taxonomy' = 'family')) |>
  distinct(taxonomy, order) |>

palette_order_assigned_bloo <- c('Thiotrichales' = '#FFA737', 'SAR11_clade' = '#B0413E', 'Rhodobacterales' = '#C55E5C',
                                   'Sphingomonadales' = '#8C000A', 'Rickettsiales' = '#5A0000',
                                 'Rhizobiales' = '#682624', 'Rhodospirillales' = '#FFA197',
                                 'Verrucomicrobiales' = '#005c69', 'Puniceicoccales' =   '#fcca46',  'Phycisphaerales' =  '#e3a6ce',   
                                 ' Flavobacteriales' =  '#0051BF',  'Sphingobacteriales'  = '#69267E' ,
                                 'SubsectionI'  =  '#009F6A',  'Bdellovibrionales' = '#000000',
                                 'Oceanospirillales' = '#DA6D00',
                                   'Alteromonadales' = '#A63B00',
                                   'Vibrionales'= '#F2AC5D', 
                                   'Enterobacteriales' = '#FFA200',
                                   'Cellvibrionales' = '#F35900', 
                                   'Pseudomonadales' = '#FF8E00')

palette_family_assigned_bloo <- c('Thiotrichaceae' = '#FFA737', 'Surface_2' = '#B0413E', 'Surface_1'= '#CD7F78', 
                                  'Rhodobacteraceae' = '#C55E5C', 'Sphingomonadaceae' = '#8C000A', 'Erythrobacteraceae' = '#690000',
                                  'SAR116_clade' = '#5A0000',
                                  'Rhodobiaceae' = '#682624',
                                  'Rhodospirillaceae'= '#FFA197',
                                  'Verrucomicrobiaceae' = '#005c69',
                                  )

palette_genus_assigned_bloo <- c()

all_tax_bloomers <- c()
# 
# palette_bloomers_tax_order <- c('Thiotrichales' =      'SAR11_clade'=        'Rhodobacterales' =
#                                   'Sphingomonadales' =   'Rickettsiales'=
#  'Rhizobiales' =       'Rhodospirillales' =  'Verrucomicrobiales' = 'Puniceicoccales' =
#    'Phycisphaerales' =
#  'Flavobacteriales' =   'Sphingobacteriales' = 'SubsectionI' =        'Bdellovibrionales' =
#    'Oceanospirillales' =
# 'Alteromonadales' =    'Vibrionales' =       'Enterobacteriales' = 'Cellvibrionales' =
#   'Pseudomonadales' = )
# 
# palette_bloomers_tax_family <- c('Thiotrichaceae',  'Surface_2 ',       'Surface_1',
#                                  'Rhodobacteraceae' ,   'Sphingomonadaceae',
#                                  'Surface_2', 'Erythrobacteraceae',  'SAR116_clade'     ,   'Rhodobiaceae'   ,     'Rhodospirillaceae',
#                                  'Verrucomicrobiaceae',
#                                  'DEV007', 'Puniceicoccaceae',
#                                  'Phycisphaeraceae' ,   'NS7_marine_group' ,   'NS9_marine_group' ,
#                                  'Cryomorphaceae'  ,    'Saprospiraceae'   ,   'Flavobacteriaceae',
#                                  'FamilyI'  ,           'Bacteriovoracaceae' ,
#                                  'Alcanivoracaceae',    'Alteromonadaceae',    'Vibrionaceae',
#                                  'Enterobacteriaceae',  'Halieaceae',
#                                  'Moraxellaceae',   'SAR86_clade')
# 
#   bloo_taxonomy %$%
#   family_f |>
#   unique()
# 

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

### plot community evenness----
community_eveness_02 |>
  bind_rows(community_eveness_3) |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, community_eveness_rar))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#94969E')+
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

### plot bray curtis dissimilarity----
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

m_bbmo_10y$fraction <- m_bbmo_10y$fraction |>
  factor(levels = c('0.2', '3'))

asv_tab_all_bloo_z_tax$season <- asv_tab_all_bloo_z_tax$season |>
  factor(levels = c('winter', 'spring', 'summer', 'autumn'))

bray_curtis_02_rar |>
  bind_rows(bray_curtis_3_rar) |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, bray_curtis_result))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_point(aes(shape = fraction, color = fraction))+
  geom_line(aes(date, bray_curtis_result, group = fraction, color = fraction))+
  #facet_grid(vars(fraction))+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity', color = 'Fraction')+
  guides(shape = 'none')+
  scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

### plot Bray-Curtis dissimilarity and Community Eveness together----
community_eveness_all <- community_eveness_02 |>
  bind_rows(community_eveness_3) 

bray_curtis_rar_all <- bray_curtis_02_rar |> ##one sample less, the first one can't be compared with the previous
  bind_rows(bray_curtis_3_rar)
  
community_eveness_all |> 
  left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, value))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_point(aes(shape = fraction, color = fraction))+
  geom_line(aes(date, value, group = fraction, color = fraction))+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  facet_grid(vars(diversity_index), labeller = labs_diversity)+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Time', y = 'Community diversity', color = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

###axis x day of the year
m_bbmo_10y |>
  colnames()

community_eveness_all |> 
  left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(day_of_year, value))+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_point(aes(shape = fraction, color = fraction))+
  #geom_line(aes(day, value, group = fraction, color = fraction))+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  facet_grid(diversity_index~fraction)+
  scale_color_manual(values= palette_fraction)+
  geom_violin(aes(group = month), alpha = 0.5)+
  #scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  #geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Time', y = 'Community diversity', color = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

  
# Discover anomalies----
## For each ASVs based on relative abundances and pseudoabundances-----
### those ASVs that are present > 50% of the sampless
asv_tab_10y_02_pseudo %$%
  sample_id |>
  unique() |>
  summary() #60 is half of the dataset

asv_tab_10y_02_zclr |>
  colnames()

###I filter by relativeabundance > 10% at some point so that the computer can compute this easily

x <- 120*0.75 ##percentage of ASV present in the dataset that we want to subset by (occurrence)

z_02 <- asv_tab_10y_02_pseudo |>
  inner_join(asv_tab_10y_02_zclr, by = c('sample_id', 'asv_num')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  #group_by(asv_num) |>
  dplyr::reframe(anomalies_ab = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = zclr, plotting = FALSE)[c(1,2,3)])

asv_tab_10y_3_pseudo %$%
  sample_id |>
  unique() |>
  summary() #

x <- 117*0.75  ##percentage of ASV present in the dataset that we want to subset by (occurrence)
z_3 <- asv_tab_10y_3_pseudo |>
  inner_join(asv_tab_10y_3_zclr, by = c('sample_id', 'asv_num')) |>
  group_by(asv_num) |>
  # dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  # dplyr::filter(num_0 <= x) |>  #filter those ASVs that are 0 more than 50% of the dataset
  #as_tibble() |>
  #group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
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

common_bloomers_tax <- asv_anom_3_tb |>
  bind_rows(asv_anom_02_tb) |>
  unique() |> ##96 totals >50% del dataset
  left_join(tax_bbmo_10y, by = c('value' = 'asv_num'))

asv_anom_3_tb |>
  anti_join(asv_anom_02_tb) #10 ASV only in 3

asv_anom_02_tb |>
  anti_join(asv_anom_3_tb) #40 ASVs only in 0.2

asv_tab_10y_02_pseudo_zclr_bloo <-  asv_tab_10y_02_pseudo_zclr |>
  #group_by(asv_num) |>
  #dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
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
  range() #max is 0.56

## Recover z_scores of diversity----
z_scores_div_02_bray <- bray_curtis_02_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_div_02_ev <- community_eveness_02 |>
  dplyr::reframe(z_score_ev = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ev) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_div_3_bray <- bray_curtis_3_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_div_3_ev <- community_eveness_3 |>
  dplyr::reframe(z_score_ev = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ev) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_bray <- z_scores_div_02_bray |>
  bind_rows(z_scores_div_3_bray)

z_scores_ev <- z_scores_div_02_ev |>
  bind_rows(z_scores_div_3_ev)

## Recover z-score for each ASV ----
### I want to highlight anomalies for each ASV to do so I recover z-scores for those ASVs that that have high z-scores
### at some point of the dataset. Easy to observe if those ASVs are having random anomalies or all of them happen at the same time
x <- 120*0.75
z_scores_02 <- asv_tab_10y_02_pseudo_zclr |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 50% of the samples
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

x <- 117*0.75 #number of 0s that we accept in the dataset to filter it by to calculate anomalies
z_scores_3 <- asv_tab_10y_3_pseudo_zclr |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 < 58) |> ##only anomalies for ASVs that are present in > 50% of the samples
  #dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 25% of the samples
  group_by(asv_num) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = 3, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num') 

##Why do I have inf values when I assess z-scores?
###understand why do we have zscores near infinite (it is because we have enormous changes of relative abundance comming from 0 values)
 z_scores_02 |>
  dplyr::filter(z_score_ra ==  'Inf') |>
  dplyr::filter(z_score_ra >= 1.96) #check that filtering by z_score is detecting infinitive values as a number

z_score_infinite <- asv_tab_all_bloo_z_tax |>
  dplyr::select(z_score_ra, asv_num, sample_id, abundance_type, abundance_value) |>
  #dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(z_score_ra == is.infinite(z_score_ra)) 

z_scores_all |>
  colnames()

asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num == 'asv38') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(z_score_ra, date, asv_num, abundance_value, location, fraction) |>
  dplyr::mutate(infinit_color = if_else(z_score_ra == 'Inf',  '#9F0011', '#080808', missing = '#080808')) |>
  ggplot(aes(date, abundance_value, color = infinit_color))+
  scale_color_identity()+
  geom_line(aes(date, abundance_value, group = location))+
  geom_point(aes(color = infinit_color))+
  facet_wrap(vars(fraction))+
  theme_bw()+
  theme(legend.position = 'rigth')

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num == 'asv1') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(z_score_ra, date, asv_num, abundance_value, location, fraction) |>
  dplyr::mutate(infinit_color = if_else(z_score_ra == 'Inf',  '#9F0011', '#080808', missing = '#080808')) |>
  ggplot(aes(date, abundance_value, color = infinit_color))+
  scale_color_identity()+
  geom_line(aes(date, abundance_value, group = location))+
  geom_point(aes(color = infinit_color))+
  facet_wrap(vars(fraction))+
  theme_bw()+
  theme(legend.position = 'rigth')
  
# Dataset with all information ----
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') 

### Bloomers taxonomy exploration----

#### try to create a hierarchical piechart 

 <- asv_tab_all_bloo_z_tax |>
  dplyr::select(asv_num, phylum, class, order, family, genus) |>
  distinct() |>
  group_by(phylum, class, order, family, genus) |>
  dplyr::mutate(genus = case_when(genus = is.na(genus) ~ 'unclassified',
                                  genus != is.na(genus) ~ genus)) |>
  dplyr::summarize(n_bloom = n()) |>
  ungroup() |>
  dplyr::mutate(bloom_perc_genus = n_bloom/sum(n_bloom)) |>
  group_by(phylum, class, order, family, genus) |>
  mutate(sum_f = sum(bloom_perc_genus)) |>
  
  group_by(phylum, class, order, family) |>
  dplyr::mutate(n_bloom = n()) |>
  ungroup() |>
  dplyr::mutate(bloom_perc_family = n_bloom/sum(n_bloom)) |>
  group_by(phylum, class, order, family) |>
  mutate(sum_f = sum(bloom_perc_family)) |>
  
  group_by(phylum, class, order) |>
  dplyr::mutate(n_bloom = n()) |>
  ungroup() |>
  dplyr::mutate(bloom_perc_order = n_bloom/sum(n_bloom)) |>
  group_by(phylum, class, order) |>
  mutate(sum_o = sum(bloom_perc_order)) |>
  
  group_by(phylum, class) |>
  dplyr::mutate(n_bloom = n()) |>
  ungroup() |>
  dplyr::mutate(bloom_perc_class = n_bloom/sum(n_bloom)) |>
  group_by(phylum, class) |>
  mutate(sum_c = sum(bloom_perc_class)) |>
  
  group_by(phylum) |>
  dplyr::mutate(n_bloom = n()) |>
  ungroup() |>
  dplyr::mutate(bloom_perc_phylum = n_bloom/sum(n_bloom)) |>
  group_by(phylum) |>
  mutate(sum_p = sum(bloom_perc_phylum)) |>
  pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage')

bloomers_tax$taxonomy_rank <- bloomers_tax$taxonomy_rank |> 
  factor(levels = c('bloom_perc_phylum', 'bloom_perc_class', 'bloom_perc_order', 'bloom_perc_family', 'bloom_perc_genus'))

bloomers_tax$taxonomy_rank <- bloomers_tax$taxonomy_rank |> 
  factor(levels = c('sum_p', 'sum_c', 'sum_o', 'sum_f', 'sum_g'))

bloomers_tax_rank <- function(data){
    data_ed <- data |>
      dplyr::select(asv_num, phylum, class, order, family, genus) |>
      distinct() |>
      dplyr::mutate(genus = case_when(genus = is.na(genus) ~ 'unclassified',
                                      genus != is.na(genus) ~ genus))
    
    bloom_g <- data_ed |>
    group_by(phylum, class, order, family, genus) |>
      dplyr::summarize(n_bloom = n()) |>
      ungroup() |>
      dplyr::mutate(bloom_perc_genus = n_bloom/sum(n_bloom)) |>
      group_by(phylum, class, order, family, genus) |>
      mutate(sum_g = sum(bloom_perc_genus)) |>
      pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage') |>
      dplyr::select(genus, taxonomy_rank, percentage) |>
      ungroup() |>
      distinct(genus, family, taxonomy_rank, percentage) |>
      rename( taxonomy = genus) |>
      dplyr::select(-family)
    
    bloom_f <- data_ed |>
      group_by(phylum, class, order, family) |>
      dplyr::mutate(n_bloom = n()) |>
      ungroup() |>
      dplyr::mutate(bloom_perc_family = n_bloom/sum(n_bloom)) |>
      group_by(phylum, class, order, family) |>
      mutate(sum_f = sum(bloom_perc_family)) |>
      pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage') |>
      ungroup() |>
      dplyr::select(family, taxonomy_rank, percentage) |>
      distinct(family, taxonomy_rank, percentage) |>
      rename( taxonomy = family)
    
   bloom_o <-data_ed |>
      group_by(phylum, class, order) |>
      dplyr::mutate(n_bloom = n()) |>
      ungroup() |>
      dplyr::mutate(bloom_perc_order = n_bloom/sum(n_bloom)) |>
      group_by(phylum, class, order) |>
      mutate(sum_o = sum(bloom_perc_order)) |>
     pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage') |>
     ungroup() |>
     dplyr::select(order, taxonomy_rank, percentage) |>
     distinct(order, taxonomy_rank, percentage) |>
     rename( taxonomy = order)
   
   bloom_c <- data_ed |>
      group_by(phylum, class) |>
      dplyr::mutate(n_bloom = n()) |>
      ungroup() |>
      dplyr::mutate(bloom_perc_class = n_bloom/sum(n_bloom)) |>
      group_by(phylum, class) |>
      mutate(sum_c = sum(bloom_perc_class)) |>
     pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage') |>
     ungroup() |>
     dplyr::select(class, taxonomy_rank, percentage) |>
     distinct(class, taxonomy_rank, percentage)|>
     rename( taxonomy = class)
      
     bloom_p <- data_ed |>
      group_by(phylum) |>
      dplyr::mutate(n_bloom = n()) |>
      ungroup() |>
      dplyr::mutate(bloom_perc_phylum = n_bloom/sum(n_bloom)) |>
      group_by(phylum) |>
       mutate(sum_p = sum(bloom_perc_phylum)) |>
       pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage') |>
       ungroup() |>
       dplyr::select(phylum, taxonomy_rank, percentage) |>
       distinct(phylum, taxonomy_rank, percentage) |>
       rename( taxonomy = phylum)

   bloomers_tax_data <- bloom_p |>
     bind_rows(bloom_c) |>
     bind_rows(bloom_o) |>
     bind_rows(bloom_f) |>
     bind_rows(bloom_g)
   
   return(bloomers_tax_data)

}

bloomers_tax <- asv_tab_all_bloo_z_tax |>
  bloomers_tax_rank()

data_ed <- asv_tab_all_bloo_z_tax |>
  dplyr::select(asv_num, phylum, class, order, family, genus) |>
  distinct() |>
  dplyr::mutate(genus = case_when(genus = is.na(genus) ~ 'unclassified',
                                  genus != is.na(genus) ~ genus))


bloom_g <- data_ed  |>
  group_by(phylum, class, order, family, genus) |>
  dplyr::summarize(n_bloom = n()) |>
  ungroup() |>
  dplyr::mutate(bloom_perc_genus = n_bloom/sum(n_bloom)) |>
  group_by(phylum, class, order, family, genus) |>
  mutate(sum_g = sum(bloom_perc_genus)) |>
  pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage') |>
  dplyr::select(genus, taxonomy_rank, percentage) |>
  ungroup() |>
  distinct(genus, family, taxonomy_rank, percentage) |>
  rename( taxonomy = genus) |>
  dplyr::select(-family)

bloom_g |>
  group_by(taxonomy_rank) |>
  summarize(n = sum(percentage))

bloom_c |>
  bind_rows(bloom_p)

tax_bbmo_10y |>
  head()

tax_bbmo_10y |>
  dplyr::select(-c(asv_num, otu_corr, seq)) |>
  pivot_longer(cols = )

bloomers_tax |>
  dplyr::filter(taxonomy_rank == 'sum_p')

bloomers_tax |>
  #dplyr::filter(taxonomy_rank == 'bloom_perc_phylum') |>
  # group_by(phylum) |>
  # dplyr::mutate(sum = if_else(taxonomy_rank == 'bloom_perc_phylum',  sum(percentage), '0', missing = NULL)) |>
  #distinct() |>
 # dplyr::select(-n_bloom) |>
  #dplyr::mutate(tax_combined = paste(phylum, class, order, family, genus)) |> 
  #distinct(tax_combined,  percentage, taxonomy_rank) |>
  ggplot(aes(y = taxonomy_rank, x = percentage, fill = taxonomy))+
 # scale_fill_manual(values = paired_12_12)+ #values = palette_phylums_assigned
  labs(x='% of potential bloomers', 'Taxonomy rank')+
  coord_polar()+
  geom_col()+
  #scale_y_continuous(expand=c(0, 0)) +
  theme_bw()

paired_12_12 <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",
                  "#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")

### since pseudoabundance is only useful for the 0.2 fraction we now only plot the relative abundance changes

##highligh harbour restoration period ----
### The remodelation of the Blanes harbour strated on 24th March 2010 and finished on the 9th of june 2012
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

## Plot Eveness and Bray Curtis anomalies----
### crec que no té sentit perquè ja ens indica un canvi en la comunitat almenys la Bray-Curtis dissimilarity
bray_curtis_02_rar |>
  bind_rows(bray_curtis_3_rar) |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, bray_curtis_result))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_point(aes(shape = fraction, color = fraction))+
  geom_line(aes(date, bray_curtis_result, group = fraction, color = fraction))+
  #facet_grid(vars(fraction))+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity', color = 'Fraction')+
  guides(shape = 'none')+
  scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

z_scores_bray <- z_scores_div_02_bray |>
  bind_rows(z_scores_div_3_bray)

z_scores_ev <- z_scores_div_02_ev |>
  bind_rows(z_scores_div_3_ev)

m_bbmo_10y$fraction <- m_bbmo_10y$fraction |>
  factor(levels = c('0.2', '3'))

asv_tab_all_bloo_z_tax$season <- asv_tab_all_bloo_z_tax$season |>
  factor(levels = c('winter', 'spring', 'summer', 'autumn'))

 bray_curtis_02_rar |>
  bind_rows(bray_curtis_3_rar) |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  left_join(z_scores_bray) |>
  dplyr::mutate(z_scores_bray = case_when(is.na(z_score_bray) ~ 0,
                                   z_score_bray == 'NaN' ~ 0,
                                   z_score_bray == Inf ~ 10000,
                                   TRUE ~ z_score_bray)) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_bray) >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, bray_curtis_result))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_line(aes(date, bray_curtis_result, group = fraction))+
   geom_point(aes(shape = fraction, color = anomaly_color))+
  facet_grid(vars(fraction))+
  scale_color_identity()+
  #scale_shape_manual(labels = labs_fraction)+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity', color = 'Fraction')+
  guides(shape = 'none')+
  scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

## Plot ASV color by class----
asv_tab_all_bloo_z_tax %$%
  class |>
  unique()

asv_tab_all_bloo_z_tax %$%
  order |>
  unique()

bloo_taxonomy <- asv_tab_all_bloo_z_tax %>%
  dplyr::select(phylum_f, class_f, order_f, family_f, asv_num_f) |>
  unique()

asv_tab_all_bloo_z_tax |>
  colnames()

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
library(ggforce)
bloomers_bbmo <- asv_tab_all_bloo_z_tax |>
  #left_join(m_3, by = 'sample_id') |>
  #left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")),
                family_asv_num = as.factor(paste0(family_f, '.', asv_num_f))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'zclr')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  #geom_vline(xintercept = c('2004-03-01'))+
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  # geom_rect(aes(xmin = '2010-03-01 01::00:00', xmax = '2012-04-01 01::00:00', ymin = 0, ymax = Inf),
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                      ymin = -Inf, ymax = Inf), fill = '#C7C7C7')+ #, alpha = 0.4
  geom_line(aes(group = fraction, color = class))+
  geom_point(aes(color = class, shape = fraction), size = 0.9)+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Class', shape = 'Fraction')+
  scale_color_manual(values = palette_bloomers_tax_class)+
  #facet_grid(fraction~abundance_type, scales = 'free')+
  #facet_grid(asv_num~fraction)+
  #facet_wrap(vars(family_asv_num), scales = 'free', ncol = 5)+
  facet_wrap_paginate(vars(family_asv_num), ncol = 2, nrow = 3, scales = 'free_y')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 10),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10), strip.background = element_blank())
# 
# ggsave( plot = bloomers_bbmo, filename = 'blo_bbmo_10y.pdf',
#         path = 'Results/Figures/',
#         width = 180, height = 230, units = 'mm')

n <- n_pages(bloomers_bbmo)

pdf('bloomers_bbmo.pdf', paper= 'A4', w= 210/25.4, 297/25.4)
for(i in 1:n){
  print(bloomers_bbmo + facet_wrap_paginate(vars(family_asv_num), ncol = 2, nrow = 3, page = i))
}
dev.off()

##plots example for PhD committee----
asv11_bbmo <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv11')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")),
                family_asv_num = as.factor(paste0(family_f, '.', asv_num_f))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'zclr')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  #geom_vline(xintercept = c('2004-03-01'))+
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  # geom_rect(aes(xmin = '2010-03-01 01::00:00', xmax = '2012-04-01 01::00:00', ymin = 0, ymax = Inf),
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#C7C7C7')+ #, alpha = 0.4
  geom_line(aes(group = fraction, color = fraction), linewidth = 1)+
  geom_point(aes(color = fraction, shape = fraction), size = 3)+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Fraction', shape = 'Fraction')+
  scale_color_manual(values = palette_fraction)+
  #facet_grid(fraction~abundance_type, scales = 'free')+
  #facet_grid(asv_num~fraction)+
  facet_grid(vars(fraction), scales = 'fixed')+
  #facet_wrap_paginate(vars(family_asv_num), ncol = 2, nrow = 3, scales = 'free_y')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 28), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 28),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28), strip.background = element_blank(), strip.text.y = element_blank(), 
        legend.text = element_text(size = 28), legend.title = element_text(size = 28))
  
ggsave( plot = asv11_bbmo, filename = 'asv11_bbmo.pdf',
        path = '../BBMO_bloomers/Results/Figures/',
        width = 380, height = 300, units = 'mm')

asv1_bbmo <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv1')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")),
                family_asv_num = as.factor(paste0(family_f, '.', asv_num_f))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'zclr')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  #geom_vline(xintercept = c('2004-03-01'))+
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  # geom_rect(aes(xmin = '2010-03-01 01::00:00', xmax = '2012-04-01 01::00:00', ymin = 0, ymax = Inf),
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#C7C7C7')+ #, alpha = 0.4
  geom_line(aes(group = fraction, color = season), linewidth = 1)+
  geom_point(aes(color = season, shape = fraction), size = 3)+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Season', shape = 'Fraction')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_grid(fraction~abundance_type, scales = 'free')+
  #facet_grid(asv_num~fraction)+
  facet_grid(vars(fraction), scales = 'fixed')+
  #facet_wrap_paginate(vars(family_asv_num), ncol = 2, nrow = 3, scales = 'free_y')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 28), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 28),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28), strip.background = element_blank(), strip.text.y = element_blank(), 
        legend.text = element_text(size =18), legend.title = element_text(size = 28))

ggsave( plot = asv1_bbmo, filename = 'asv1_bbmo.pdf',
        path = '../BBMO_bloomers/Results/Figures/',
        width = 380, height = 300, units = 'mm')


asv17_bbmo <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv17')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")),
                family_asv_num = as.factor(paste0(family_f, '.', asv_num_f))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'zclr')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  #geom_vline(xintercept = c('2004-03-01'))+
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7')+ #, alpha = 0.4
  geom_line(aes(group = fraction, color = fraction), linewidth = 1)+
  geom_point(aes(color = fraction, shape = fraction), size = 3)+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Fraction', shape = 'Fraction')+
  scale_color_manual(values = palette_fraction)+
  #facet_grid(fraction~abundance_type, scales = 'free')+
  #facet_grid(asv_num~fraction)+
  facet_grid(vars(fraction), scales = 'fixed')+
  #facet_wrap_paginate(vars(family_asv_num), ncol = 2, nrow = 3, scales = 'free_y')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 28), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 28),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28), strip.background = element_blank(), strip.text.y = element_blank(), 
        legend.text = element_text(size = 28), legend.title = element_text(size = 28))
  

ggsave( plot = asv17_bbmo, filename = 'asv17_bbmo.pdf',
        path = '../BBMO_bloomers/Results/Figures/',
        width = 380, height = 300, units = 'mm')

asv192_bbmo <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv192')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")),
                family_asv_num = as.factor(paste0(family_f, '.', asv_num_f))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'zclr')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  #geom_vline(xintercept = c('2004-03-01'))+
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7')+ #, alpha = 0.4
  geom_line(aes(group = fraction, color = fraction), linewidth = 1)+
  geom_point(aes(color = fraction, shape = fraction), size = 3)+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Fraction', shape = 'Fraction')+
  scale_color_manual(values = palette_fraction)+
  #facet_grid(fraction~abundance_type, scales = 'free')+
  #facet_grid(asv_num~fraction)+
  facet_grid(vars(fraction), scales = 'fixed')+
  #facet_wrap_paginate(vars(family_asv_num), ncol = 2, nrow = 3, scales = 'free_y')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 28), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 28),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28), strip.background = element_blank(), strip.text.y = element_blank(), 
        legend.text = element_text(size = 28), legend.title = element_text(size = 28))


ggsave( plot = asv192_bbmo, filename = 'asv192_bbmo.pdf',
        path = '../BBMO_bloomers/Results/Figures/',
        width = 380, height = 300, units = 'mm')


asv282_bbmo <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv282')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")),
                family_asv_num = as.factor(paste0(family_f, '.', asv_num_f))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'zclr')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  #geom_vline(xintercept = c('2004-03-01'))+
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#C7C7C7')+ #, alpha = 0.4
  geom_line(aes(group = fraction, color = fraction), linewidth = 1)+
  geom_point(aes(color = fraction, shape = fraction), size = 3)+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Fraction', shape = 'Fraction')+
  scale_color_manual(values = palette_fraction)+
  #facet_grid(fraction~abundance_type, scales = 'free')+
  #facet_grid(asv_num~fraction)+
  facet_grid(vars(fraction), scales = 'fixed')+
  #facet_wrap_paginate(vars(family_asv_num), ncol = 2, nrow = 3, scales = 'free_y')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 28), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 28),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28), strip.background = element_blank(), strip.text.y = element_blank(), 
        legend.text = element_text(size = 28), legend.title = element_text(size = 28))


ggsave( plot = asv282_bbmo, filename = 'asv282_bbmo.pdf',
        path = '../BBMO_bloomers/Results/Figures/',
        width = 380, height = 300, units = 'mm')

asv282_bbmo <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv282')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")),
                family_asv_num = as.factor(paste0(family_f, '.', asv_num_f))) |>
  dplyr::filter(!abundance_type %in% c('pseudoabundance', 'relative_abundance')) |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  #geom_vline(xintercept = c('2004-03-01'))+
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#C7C7C7')+ #, alpha = 0.4
  geom_line(aes(group = fraction, color = fraction), linewidth = 1)+
  geom_point(aes(color = fraction, shape = fraction), size = 3)+
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'zCLR', color = 'Fraction', shape = 'Fraction')+
  scale_color_manual(values = palette_fraction)+
  #facet_grid(fraction~abundance_type, scales = 'free')+
  #facet_grid(asv_num~fraction)+
  facet_grid(vars(fraction), scales = 'fixed')+
  #facet_wrap_paginate(vars(family_asv_num), ncol = 2, nrow = 3, scales = 'free_y')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 28), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 28),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28), strip.background = element_blank(), strip.text.y = element_blank(), 
        legend.text = element_text(size = 28), legend.title = element_text(size = 28))


ggsave( plot = asv282_bbmo, filename = 'asv282_bbmo_zclr.pdf',
        path = '../BBMO_bloomers/Results/Figures/',
        width = 380, height = 300, units = 'mm')


##abundance changes over the dataset----
asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, bacteria_joint))+
  scale_x_datetime()+
  scale_y_continuous(labels = scales::scientific)+#, breaks = c(5e+00, 1e+01), limits = c(2e+00, 1.25e+01)
  geom_point()+
  geom_line()+
  labs(x = 'Time', y = 'Prokaryotic Abundance (cells/mL)', color = 'Class', shape = 'Fraction')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 5),
        axis.ticks = element_blank(), legend.position = 'right', axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 5), strip.background = element_blank())

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

###shape anomaly----
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

###plot zclr scores----
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

### Identify those that are seasonal from those that are not seasonal----
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

###plot those by anomaly presence colored in red----
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

##Ocurrence of this AVSs vs frequency of blooming events and magnitude of the events --------
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
                                   z_score_ra = is.na(z_score_ra) ~ 0,
                                   z_score_ra = is.infinite(z_score_ra) ~ 10000,
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
                                          z_score_ra = is.na(z_score_ra) ~ 0,
                                          z_score_ra = is.infinite(z_score_ra) ~ 10000,
                                          TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly = case_when(z_score_ra >= 1.96 ~ 1,
                                    z_score_ra < 1.96 ~ 0)) |>
  group_by(asv_num) |>
  dplyr::summarize(n_anom = sum(anomaly)) |>
  dplyr::mutate(anom_perc = n_anom/nsamples_02,
                fraction = '0.2')

anom_perc <- anom_perc_3 |>
  bind_rows(anom_perc_02)

##per year (explore the % of blooming events per tax)----
anom_perc_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  dplyr::filter(abundance_value >= 0.1) |>
  dplyr::mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
                                          z_score_ra == 'NaN' ~ 0,
                                          z_score_ra = is.na(z_score_ra) ~ 0,
                                          z_score_ra = is.infinite(z_score_ra) ~ 10000,
                                          TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly = case_when(z_score_ra >= 1.96 ~ 1,
                                    z_score_ra < 1.96 ~ 0)) |>
  group_by(year, month) |>
  dplyr::summarize(n_anom = sum(anomaly)) |>
  dplyr::mutate(#anom_perc = n_anom/12,
                fraction = '0.2')

anom_perc_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  dplyr::filter(abundance_value >= 0.1) |>
  dplyr::mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
                                          z_score_ra == 'NaN' ~ 0,
                                          z_score_ra = is.na(z_score_ra) ~ 0,
                                          z_score_ra = is.infinite(z_score_ra) ~ 10000,
                                          TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly = case_when(z_score_ra >= 1.96 ~ 1,
                                    z_score_ra < 1.96 ~ 0)) |>
  group_by(year, month) |>
  dplyr::summarize(n_anom = sum(anomaly)) |>
  dplyr::mutate(#anom_perc = n_anom/12,
                fraction = '3')

anom_perc_02 |>
  bind_rows(anom_perc_3) |>
  ggplot(aes(month, n_anom, shape = fraction))+
 # geom_point()+
  geom_col()+
  #geom_line()+
  geom_smooth(method= 'loess')+
  facet_grid(fraction~year)+
  #scale_y_continuous(labels = percent_format())+
  theme_bw()

###max relative abundance (in a blooming event)
asv_tab_all_bloo_z_tax |>
  colnames()

max_rel_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(z_score_ra > 1.96) |>
  dplyr::filter(fraction == '3' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  group_by(asv_num) |>
  dplyr::summarize(max_rel = max(abundance_value)) |>
  dplyr::mutate(fraction = '3')

max_rel_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(z_score_ra > 1.96) |>
  dplyr::filter(fraction == '0.2' &
                  !abundance_type %in% c('pseudoabundance', 'zclr')) |>
  group_by(asv_num) |>
  dplyr::summarize(max_rel = max(abundance_value))|>
  dplyr::mutate(fraction = '0.2')

max_rel <- max_rel_3 |>
  bind_rows(max_rel_02)
  
###plot Occurence vs frequency of blooming events
bloom_occurrence_tax <- occurence_perc |>
  left_join(anom_perc, by = c('asv_num','fraction')) |>
  left_join(max_rel, by = c('asv_num','fraction')) |>
  left_join(tax_bbmo_10y, by = 'asv_num') 

occurrence_perc_tax <-occurence_perc  |>
 left_join(anom_perc, by = c('asv_num','fraction')) |>
  left_join(max_rel, by = c('asv_num','fraction')) |>
  left_join(tax_bbmo_10y, by = 'asv_num')

occurence_perc_tax |>
  colnames()

occurence_perc_tax %$%
  unique(genus)

occurence_perc_tax %$%
  unique(family)

occurrence_perc_tax |>
  ggplot(aes(occurence_perc, anom_perc, color = class, size = max_rel*100, shape = fraction))+
  geom_point()+
  scale_color_manual(values = palette_class_assigned)+
  labs(x = 'ASV occurrence (%)', y = 'Anomaly (%)', color = 'Class', size = 'Max relative\nabundance (%)\nin a blooming event')+
  scale_x_continuous(labels = percent_format())+
  scale_y_continuous(labels = percent_format())+
  facet_wrap(vars(fraction))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'right', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

##distribution of zcores in the dataset----

asv_tab_all_bloo_z_tax |>
  colnames()

library(ggridges)

asv_tab_all_bloo_z_tax |>
  dplyr::mutate(sampling = '10Y_BBMO') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(z_score_ra != is.na(z_score_ra) &
                  z_score_ra != is.infinite(z_score_ra)) |>
  ggplot(aes(x = as.numeric(z_score_ra)))+
  geom_density_line()+
  scale_x_continuous(limits = c(-10, 10))
  # geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 1,
  #                     jittered_points = TRUE,
  #                     point_shape = 21, point_size = 0.2, point_alpha = 0.0,
  #                     quantile_lines = TRUE,
  #                     quantile_fun = mean)

# ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f,  group = phylum_f, label = counts))+
#   geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 1,
#                       jittered_points = TRUE,
#                       point_shape = 21, point_size = 0.2, point_alpha = 0.0,
#                       quantile_lines = TRUE,
#                       quantile_fun = mean

##number of events/ year




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

 