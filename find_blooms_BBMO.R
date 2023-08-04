
# packages----
library(readxl)
library(tidyverse)
library(janitor)
library(Bloomers)
library(magrittr)
library(scales)
library(phyloseq)
library(speedyseq)

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
  ntaxa() #8052 but if we filter for those that are 0 during the whole datsaeet then we got 7849 ASVs

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
  dplyr::filter(fraction == 0.2) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02)))

m_3 <- m_bbmo_10y |>
  dplyr::filter(fraction == 3.0) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3)))

# Calculate relative abundance----
asv_tab_bbmo_10y_l_rel <- asv_tab_bbmo_10y_l |>
  calculate_rel_abund(group_cols = sample_id)

asv_tab_10y_filt_3_rel <- asv_tab_bbmo_10y_l_rel %>%
  dplyr::filter(sample_id %in% m_3$sample_id)

asv_tab_10y_filt_02_rel <- asv_tab_bbmo_10y_l_rel %>%
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
## This type of normalization makes more sense with the 0.2 fraction since citometry accounts mainly for this fraction of the biomass
m_bbmo_10y |>
  colnames()

asv_tab_10y_filt_3_pseudo <- asv_tab_10y_filt_3_rel |>
calculate_pseudoabund(abund_data = m_bbmo_10y, rel_abund = relative_abundance, 
                      total_abund = bacteria_joint, 
                      by_ = 'sample_id')

asv_tab_10y_filt_02_pseudo <- asv_tab_10y_filt_02_rel |>
  calculate_pseudoabund(abund_data = m_bbmo_10y, 
                        rel_abund = relative_abundance, 
                        total_abund = bacteria_joint, 
                        by_ = 'sample_id')


# Calculate diversity parameters ----
## We calculate different diversity parameters to check if the blooming events detected have an effect on the community structure
## Community Evenness----
### Maybe we need to rarefy or rarification because we have uneven sequencing effort.
### we need original reads
asv_tab_bbmo_10y_l |>
  str()

##I obtain NaN values... why? Too much 0?
community_eveness_02 <- asv_tab_bbmo_10y_l |>
  mutate(reads = as.numeric(reads)) |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_result = community_evenness(abundances = reads, index = "Pielou"))

community_eveness_3 <- asv_tab_bbmo_10y_l |>
  mutate(reads = as.numeric(reads)) |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_result = community_evenness(abundances = reads, index = "Pielou"))

asv_tab_10y_filt_3_pseudo <- asv_tab_10y_filt_3_pseudo  |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

asv_tab_10y_filt_02_pseudo <- asv_tab_10y_filt_02_pseudo  |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

## Bray Curtis dissimilarity----
### 0 means the two sites have the same composition (that is they share all the species), and 1 means the two sites 
### do not share any species.

source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R') #we need to upload it since we I updated the function but I didn't compile the package

bray_curtis_02 <- dissimilarity_matrix(data = asv_tab_10y_filt_02_rel, 
                                    sample_id_col = sample_id,
                                    values_cols_prefix = 'BL')


bray_curtis_3 <- dissimilarity_matrix(data = asv_tab_10y_filt_3_rel, 
                                       sample_id_col = sample_id,
                                       values_cols_prefix = 'BL')

###plot bray curtis dissimilarity
bray_curtis_02 |>
  bind_rows(bray_curtis_3) |>
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
asv_tab_10y_filt_02_pseudo %$%
  sample_id |>
  unique() |>
  summary() #60 is half of the dataset

z_02 <- asv_tab_10y_filt_02_pseudo |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::filter(num_0 < 60) |>
  #as_tibble() |>
  group_by(asv_num) |>
  dplyr::reframe(anomalies_ab = get_anomalies(time_lag = 3, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_ra = get_anomalies(time_lag = 3, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)])

z_3 <- asv_tab_10y_filt_3_pseudo |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::filter(num_0 < 60) |>  #filter those ASVs that are 0 more than 50% of the dataset
  #as_tibble() |>
  group_by(asv_num) |>
  dplyr::reframe(anomalies_ab = get_anomalies(time_lag = 3, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
    anomalies_ra = get_anomalies(time_lag = 3, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)])

## At the level of community, we use the eveness result and bray curtis dissimilarity ----
### sine I have NaNs values in eveness I don't calculate the anomalies for this parameter (to be done!)
z_diversity <- bray_curtis_02 |>
  dplyr::right_join(community_eveness, by = join_by("samples" == "sample_id")) |> 
  #ungroup() %>%
  #group_by(sample_id) %>%
  dplyr::reframe(anomalies_bray = get_anomalies(time_lag = 4, values = bray_curtis_result, plotting = TRUE)[c(1,2,3)])# ,
                 #anomalies_eveness = get_anomalies(time_lag = 4, values = community_eveness_result, plotting = TRUE)[c(1,2,3)])

z_diversity %>%
  str()

# Filter the ASV_tab by only those ASVs that have an anomaly at some point of the dataset----
#find_asv_with_anomalies <- function(anomalies_result, anomaly_in1, anomaly_in2, logic1 = TRUE, logic2 = TRUE, asv_col = asv_num) {
  # if(is.list(anomalies_result) == FALSE){
  #   stop("Function stopped: anomalies_result needs to be a list form the get_anomalies function")
  # }
  # if(is.logical({{anomaly_in1}}) == FALSE){
  #   stop("Function stopped: anomaly_in1 needs to be logical (TRUE/FALSE)")
  # }
  
  asv_potential_bloomers <-
    anomalies_result |>
    dplyr::filter(if (!is.null(logic1)) {{anomaly_in1}} %in% logic1 
                  else TRUE) |>
    dplyr::filter(if (!is.null(logic2)) {{anomaly_in2}} %in% logic2 
                  else TRUE) |>
    dplyr::select({{asv_col}}) |>
    as_vector()
  
  return(asv_potential_bloomers)
#}

asv_anom_02 <- find_asv_with_anomalies(anomalies_result = z_02, anomaly_in1 = anomalies_ab, anomaly_in2 = anomalies_ra, 
                                    logic1 = 'TRUE', logic2 = 'TRUE', 
                                    asv_col = asv_num)

asv_tab_all_perc_filt_02_long_filt <-  asv_tab_10y_filt_02_pseudo |>
  group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(pseudoabundance, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  filter(asv_num %in% asv_anom_02)

asv_anom_3 <- find_asv_with_anomalies(anomalies_result = z_3, anomaly_in1 = anomalies_ab, anomaly_in2 = anomalies_ra, 
                                       logic1 = 'TRUE', logic2 = 'TRUE', 
                                       asv_col = asv_num)

asv_tab_all_perc_filt_3_long_filt <-  asv_tab_10y_filt_3_pseudo |>
  group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(pseudoabundance, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  filter(asv_num %in% asv_anom_3)

# asv_tab_all_perc_filt_3_long_filt |>
#   group_by(sample_id) |>
#   dplyr::summarize(sum = sum(relative_abundance))

asv_tab_all_perc_filt_3_long$relative_abundance |> 
  range() #max is 0.56
  
asv_tab_all_perc_filt_3_long_filt |>
  colnames()

### since pseudoabundance is only useful for the 0.2 fraction we now only plot the relative abundance changes
asv_tab_all_perc_filt_02_long_filt |>
  bind_rows(asv_tab_all_perc_filt_3_long_filt) %$%
  asv_num |>
  unique() 

asv_tab_all_perc_filt_02_long_filt |>
  colnames()

##highligh harbour restoration period
harbour_restoration <-  asv_tab_all_perc_filt_3_long_filt |> 
  subset(date = between(date, '2010-03-01','2012-04-01')) |>
  dplyr::summarize(xmin=min(date),xmax=max(date))

## color by class
asv_tab_all_perc_filt_3_long_filt |>
  bind_rows(asv_tab_all_perc_filt_02_long_filt) |>
  #left_join(m_3, by = 'sample_id') |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type != 'pseudoabundance') |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value))+ #, color = 'Class' 
  scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  geom_rect(aes(xmin = '2010-03-01 01::00:00', xmax = '2012-04-01 01::00:00', ymin = 0, ymax = Inf),
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_line(aes(group = asv_num, color = class))+
  geom_point(aes(color = class))+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Class')+
  scale_color_manual(values = palette_class_assigned)+
  facet_grid(vars(fraction), scales = 'free')+
  #facet_grid(fraction~class)+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())


## Recover z-score for each ASV ----
### I want to highlight anomlies for each ASV to do so I recover z-scores for those ASVs that that have high z-scores
### at some point of the dataset. Easy to observe if those ASVs are having random anomalies or all of them happen at the same time

z_scores_02 <- asv_tab_10y_filt_02_pseudo |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 < 60) |> ##only anomalies for ASVs that are present in > 50% of the samples
  dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
  group_by(asv_num) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = 3, negative = FALSE, cutoff = 1.96, 
                                              na_rm = TRUE, values = relative_abundance, 
                                              plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_3 <- asv_tab_10y_filt_3_pseudo |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 < 60) |> ##only anomalies for ASVs that are present in > 50% of the samples
  dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
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

##shape anomaly
asv_tab_all_perc_filt_all_long_filt |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type != 'pseudoabundance') |>
  #ungroup() |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  ggplot(aes(date, abundance_value, shape = ifelse(z_score_ra > 1.96, '8', '19')))+ #, color = 'Class' 
  scale_x_datetime()+
  facet_wrap(vars(class), scales = 'free')+
  geom_line(aes(group = asv_num, color = class))+
  geom_point(aes(color = class))+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Class')+
  scale_color_manual(values = palette_class_assigned)+
  facet_grid(vars(fraction), scales = 'free')+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

##color anomaly
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

##recover all ASVs that presentend and anomaly in PA or FL
asv_tab_all_perc_filt_3_long_filt <-  asv_tab_10y_filt_3_pseudo |>
  group_by(asv_num) |>
  #dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(pseudoabundance, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  filter(asv_num %in% asv_anom_all)

asv_tab_all_perc_filt_02_long_filt %$%
  asv_num |>
  unique()

asv_tab_all_perc_filt_02_long_filt <- asv_tab_10y_filt_02_pseudo |>
  group_by(asv_num) |>
  #dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(pseudoabundance, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  filter(asv_num %in% asv_anom_all)

asv_tab_z_scores_all <- asv_tab_all_perc_filt_3_long_filt |>
  bind_rows(asv_tab_all_perc_filt_02_long_filt) |>
  left_join(z_scores_all) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type != 'pseudoabundance') |>
  group_by(asv_num) |>
  dplyr::filter(any(abundance_value >=  0.05)) #only blooms that arribe at least at 5% of the community

asv_tab_z_scores_all |>
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
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Anomaly')+ #, shpae = 'Class'
      scale_color_identity()+
  #scale_color_manual(values = if_else(asv_tab_z_scores_all$z_score_ra >= 1.96,  '#9F0011', '#080808', missing = '#080808'))+
  facet_grid(vars(fraction), scales = 'free')+
  #facet_grid(fraction~class, scales = 'free')+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())
  
#pivot_wider(id_cols = sample_id_num, values_from = anomalies_ra, names_from = asv_num)
  
# Analysis of FL and PA blooming events ----
## Do blooming events synchronize in PA and FL? 
## Some ASVs prefer to bloom in one fraction than another? Is there a preference?
## plot by fractions
asv_tab_all_perc_filt_all_long_filt |>
  colnames()

asv_tab_all_perc_filt_all_long_filt <-  asv_tab_all_perc_filt_02_long_filt |>
  bind_rows(asv_tab_all_perc_filt_3_long_filt) |>
  as_tibble()

asv_tab_all_perc_filt_all_long_filt$fraction <- asv_tab_all_perc_filt_all_long_filt$fraction  |> 
  factor(levels =  (c('0.2', '3')))

asv_tab_all_perc_filt_all_long_filt |>
  #left_join(, by = c('sample_id' = 'Code_ed')) |>
  left_join(tax_bbmo_10y, by = 'asv_num') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type != 'pseudoabundance') |>
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
  

  