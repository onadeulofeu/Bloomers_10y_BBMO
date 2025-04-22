# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                                                             ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                         All figures code                    ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# packages ----
library(tidyverse)
library(ggplot2)
library(vegan)
library(purrr)
library(broom)
library(ggtree)
library(ape)
library(cowplot)
library(ggpubr)
library(stats)
library(Bloomers)
library(speedyseq)
library(car)  ## homocedasticity 
library(dunn.test) ## significant groups in Kruskal Test
library(ggpubr) ## ggqqplot

source('../../Bloomers/R/find_asv_with_anomalies.R')
source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R')

# upload data ----
asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked_rclr.csv') |> ##using dada2 classifier assign tax with silva 138.1 and correctly identifying bloomers
  as_tibble() |>
  dplyr::select(-X) |>
  pivot_longer(cols = c('relative_abundance', 'rclr'), values_to = 'abundance_value', names_to = 'abundance_type')

bloo_all_types_summary_tax <- read.csv('results/tables/bloo_all_types_summary_tb_tax_v2.csv') |>
  dplyr::select(-X)

occurrence_bloo_bbmo <- read.delim2('data/occurrence_bloo_bbmo.csv', sep = ',') |>
  dplyr::mutate(occurrence_category = ifelse(occurrence_perc > 2/3, 'broad',
                                             ifelse(occurrence_perc < 1/3, 'narrow',
                                                    'intermediate')))

harbour_restoration_dec <- tibble(xmin = '2010.24', xmax = '2010.52') |>
  dplyr::mutate(date_min = as.numeric(xmin),
                date_max = as.numeric(xmax))

bloo_02 <- read.csv('data/detect_bloo/bloo_02.csv') |>
  as_tibble() |>
  dplyr::filter(!value %in% c('asv2', 'asv3', 'asv5', 'asv8'))

bloo_3 <- read.csv('data/detect_bloo/bloo_3.csv') |>
  as_tibble()

asv_tab_bbmo_10y_w_rar <- read.csv2('data/asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()|> 
  rename('sample_id' = 'X') 

bbmo_10y <-readRDS("data/blphy10years.rds") ## 8052 all samples, no filtering

bbmo_10y <-
  prune_taxa(taxa_sums(bbmo_10y@otu_table) >0, ##filter ASVs that are 0 in the whole dataset
             bbmo_10y)

## separate datasets by ASV_tab, taxonomy and metadata
asv_tab_bbmo_10y_l <- bbmo_10y@otu_table |>
  as_tibble()

m_bbmo_10y <- bbmo_10y@sam_data |>  
  as_tibble()

colnames(asv_tab_bbmo_10y_l) <- c('asv_num', "sample_id", 'reads')

colnames(m_bbmo_10y) <- c('sample_id', "project", "location", "code",             
                          "type", "samname", "fraction", "run",               
                          "date", "basics", "julian_day", "day_of_year",       
                          "decimal_date", "position", "sampling_time", "day_length",        
                          "temperature", "secchi", "salinity", "chla_total",     
                          "chla_3um", "PO4" ,"NH4", "NO2" ,              
                          "NO3",  "Si", "BP_FC1.55", "PNF_Micro",         
                          "PNF2_5um_Micro", "PNF_5um_Micro", "cryptomonas", "micromonas",        
                          "HNF_Micro", "HNF2_5um_Micro", "HNF_5um_Micro", "LNA",               
                          "HNA", "prochlorococcus_FC", "Peuk1",  "Peuk2",          
                          "year", "month", "day", "season",            
                          "bacteria_joint", "synechococcus", "depth", "name_complete")

m_02 <- m_bbmo_10y  |>
  dplyr::filter(fraction == 0.2) 

m_02 <- m_02 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02)))

m_3 <- m_bbmo_10y |>
  dplyr::filter(fraction == 3.0)

m_3 <- m_3 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3)))

tree_complete <- ggtree::read.tree('data/raxml/complete_tree_cesga/bbmo.raxml.support')

## add blooming events or not 
bloom_event <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(bloom_event = case_when(abundance_type == 'relative_abundance' &
                                          abundance_value >= 0.1 &
                                          z_score_ra > cut_off_value_ra &
                                          z_score_rclr > cut_off_value_rclr ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  filter(bloom_event != 0) |>
  dplyr::select(date, asv_num, bloom_event, fraction) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ungroup() |>
  distinct()

### we work with two different datasets one for FL and the other for PA
wavelet_3_df <- read.csv2('data/wavelet_3_df_deconstand.csv', sep = ',') |>
  as_tibble() |>
  dplyr::select(-X)

wavelet_02_df <- read.csv2('data/wavelet_02_df_deconstand.csv') |>
  as_tibble() |>
  dplyr::select(-X)

### interpolated environmental variables
env_data_interpolated_values_all <- read.csv2('data/env_data/env_data_interpolated_values_all.csv') |>
  rename(sample_id_num = X)

## correlations datasets genes and community
bray_unifrac_eucl_tb_02 <- read.csv('data/bray_unifrac_eucl_tb_02.csv')
bray_unifrac_eucl_tb <- read.csv('data/bray_unifrac_eucl_tb.csv')
corr_bray_genes_community_tb <- read.csv('data/corr_bray_genes_community_tb_scg.csv')

## MDR - S map results
mdr_tb <- read.csv2('../EDM_carmen/MDR/Bl_nin120_cvunit0.025_aenet_jcof_Nmvx_Rallx_demo_v2.csv',
                    header = T ) |>
  as_tibble() ## new data

harbour_restoration_dec <- tibble(xmin = '2010.24', xmax = '2010.52') |>
  dplyr::mutate(date_min = as.numeric(xmin),
                date_max = as.numeric(xmax))

# palettes ----
palette_occurrence <- c(narrow = "#AE659B",
                        intermediate = "#3e3e3e",
                        broad = "#57a9a8")

palette_diversity <- c('genes' = "#270000",
                       "wunifrac_distance" = "#898989",
                       "euclidean_distance" = "#008241",
                       "bray_curtis_community" = "#a21e66", 
                       'bray_curtis_kmers' = '#FFCAFB')

palette_fraction_env <- c("env" = "#008241", 
                          'bio' =  '#008241',
                          'phch' = '#ffa737',
                          '0.2' = '#00808F', 
                          '3' = "#454545")

palette_phylums_assigned <- c('Proteobacteria' = "#ffb900","Bacteroidota" = "#0051BF" , 'Actinobacteriota' = "#b0413e",
                              'Cyanobacteria' = "#009e73",'Crenarchaeota' = "#ffa737", 'Verrucomicrobiota' = '#005c69',
                              'Planctomycetota' = "#69267e", 'Acidobacteriota' = "#1f78b4",'Bdellovibrionota' = "#8c789d",
                              'Firmicutes' = "#637b88", 'Myxococcota'= "#003029", 'Nitrospinota'= "#e3a6ce",
                              'Campilobacterota' = "#002960", 'Deinococcota'= "#ba9864",'Fusobacteriota' ="#fb9a99",
                              'Desulfobacterota' = "#005c69", 'Methylomirabilota' = "#000000" ,
                              'Gemmatimonadota' = "#c55e5c", 'Chloroflexi' = "#00d198", 'Spirochaetota' = "#5cb1d6",
                              'Calditrichota' = "#8a007a", 'Halobacterota' = "#b79c64", 'Nitrospirota' = "#41815f",
                              'Dependentiae' = "#5b95e5", 'Patescibacteria' = "#33af9c",'Cloacimonadota' = "#fbed5c",
                              'Synergistota' = "#ce7800", 'Abditibacteriota' = "#87878b", 'Deferribacterota' = "#4dbaa9") #,NA ==  "#000000"

palette_phylums_assigned_bloo <- c('Proteobacteria' = "#fcca46", 'Cyanobacteria' = "#009e73",
                                   "Bacteroidota" = "#0051BF",  'Verrucomicrobiota' = '#005c69', 'Planctomycetota' = "#69267e",
                                   'Bdellovibrionota' = "#8c789d") #, NA == "#000000"

palette_class_assigned <- c('Gammaproteobacteria' = '#fcca46', 'Alphaproteobacteria' = '#8C000A', 
                            'Zetaproteobacteria' = '#EBCD92',
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
                            'Desulfuromonadia' = '#6D7DBE', 'Desulfobacteria' = '#000036', 'Methylomirabilia' = '#ffff00',
                            'Gemmatimonadetes' = '#c55e5c',
                            'Anaerolineae' = '#00d198', 'Chloroflexia' = '#009F6A',
                            'Spirochaetia' = '#5CB1D6', 'Leptospirae' = '#005576', 
                            'Calditrichia'  = '#8a007a', 'Halobacteria' = '#b79c64', 
                            'Nitrospiria' = '#41815f',
                            'Babeliae' = '#5b95e5', 'Saccharimonadia' = '#33af9c', 
                            'Cloacimonadia' = '#fbed5c',
                            'Synergistia' = '#ce7800', 'Abditibacteria' = '#87878b', 
                            'Deferribacteres' = '#4dbaa9',
                            'env' = "#4A785C") #, NA ==  "#000000"

palette_class_assigned_bloo <- c('Gammaproteobacteria' = '#FFA737', 'Alphaproteobacteria' = '#B0413E', 
                                 'Cyanobacteriia'  =  '#009F6A', 'Bacteroidia' = '#0051BF',  
                                 'Verrucomicrobiae' = '#005c69', 'Phycisphaerae' =  '#e3a6ce',   
                                 'Bdellovibrionia' = '#8C789D') #, NA == "#000000"

palette_order_assigned_all <-  c("SAR11 clade" =      "#ca6094",
                                 "Rhodospirillales"  = '#FFA180', 
                                 "Sphingomonadales"  = '#8C000A', 
                                 'Balneolales' = '#00425E',
                                 'Bacillales' = '#00598D',
                                 'Bradymonadales' = '#B89EF1',
                                 "Puniceispirillales" = '#79AD5A',
                                 'Rhizobiales' = '#B31722',    
                                 "Rhodobacterales" = '#2d373b',
                                 "Verrucomicrobiales"= '#005c69',
                                 "Opitutales"   =   '#74B9C8',
                                 "Phycisphaerales"  = '#e3a6ce', 
                                 "Flavobacteriales"   =  '#0051BF', 
                                 'Chitinophagales' = '#92ABFF', 
                                 "Synechococcales"  = '#009F6A', 
                                 "Bacteriovoracales" = '#8C789D',
                                 "Pseudomonadales"  = '#FF8E00', 
                                 "Enterobacterales" = "#f1c510",
                                 'Thiotrichales' =  "#000000",
                                 'Nitrososphaerales' = '#E64D49',
                                 'Pelagibacterales' = '#AB7E81',
                                 'PCC-6307' =  '#7EAB9C',
                                 'Poseidoniales' = '#EBB473',
                                 'PS1' = '#60A65E',
                                 'SAR86' = '#EBE773', 
                                 'env' = "#4A785C",
                                 'Cellvibrionales' = '#6B5C57') 

palette_order_assigned_bloo <-  c("SAR11 clade" =      "#ca6094",
                                  "Rhodospirillales"  = '#FFA180', 
                                  "Sphingomonadales"  = '#8C000A', 
                                  "Puniceispirillales" = '#4cb50f',
                                  'Rhizobiales' = '#B31722',    
                                  "Rhodobacterales" = '#2d373b',
                                  "Verrucomicrobiales"= '#005c69',
                                  "Opitutales"   =   '#74B9C8',
                                  "Phycisphaerales"  = '#e3a6ce', 
                                  "Flavobacteriales"   =  '#0051BF', 
                                  'Chitinophagales' = '#92ABFF', 
                                  "Synechococcales"  = '#009F6A', 
                                  "Bacteriovoracales" = '#8C789D',
                                  "Pseudomonadales"  = '#FF8E00', 
                                  "Enterobacterales" = "#f1c510",
                                  'Thiotrichales' =  "#000000",
                                  'env' = "#4A785C") 

palette_family_assigned_bloo <- c(##gammaproteobacteria
  "Thiotrichaceae" = "#000000",  
  "Marinobacteraceae" =  '#B5A9A4',  
  "Alcanivoracaceae1"  =  '#FCE005',  
  "Moraxellaceae"  =  '#FF8E00', #pseudomonadales
  "Halieaceae" = '#FF7B38', "SAR86 clade" = '#FFA200', #pseudomonadales
  "Vibrionaceae"      = '#1C1B1F',   "Yersiniaceae"  = '#89584D',    #enterobacterales    
  "Alteromonadaceae"   =  '#5D5B59',    ##alphaproteobacteria
  "Clade II" = '#B0413E',  "Clade I" = '#cd7f78',#"#ca6094", 
  "Rhodobacteraceae" = '#C55E5C',
  "Sphingomonadaceae"   = '#8C000A',  
  "SAR116 clade"  = '#4cb50f', 
  "Stappiaceae" = '#B31722',
  "AEGEAN-169 marine group" =  '#3D0000', 
  "Cyanobiaceae"  = '#009F6A', 
  "NS7 marine group" =  '#92ABFF',
  "NS9 marine group"  =  '#3B52A3', "Cryomorphaceae" = '#002A8E',
  "Saprospiraceae" = '#5F7CCB',
  "Flavobacteriaceae"   =  '#0051BF',  
  "Rubritaleaceae"  =  '#005c69',  "DEV007" = '#74B9C8', 
  
  "Puniceicoccaceae"   = '#29335C',     
  "Phycisphaeraceae"   = '#e3a6ce',    
  "Bacteriovoracaceae" =  '#8C789D' 
)  # NA == "#000000" ## add genus color 

palette_genus_assigned_bloo <- c('unclassified' = '#534F4A', 
                                 "Marixanthomonas" =  '#0051BF',  "NS4 marine group" = '#46ACC2',                
                                 "NS5 marine group" = '#2B4162', "Peredibacter"  = '#8C789D',
                                 "Prochlorococcus MIT9313" =  '#009F6A',       
                                 "Synechococcus CC9902" = '#004501',
                                 "Urania-1B-19 marine sediment group" = '#e3a6ce',
                                 "Candidatus Puniceispirillum" =  '#D3BF27',    
                                 "Amylibacter" =  '#C55E5C', "Ascidiaceihabitans" = '#B0413E',
                                 "HIMB11" =   '#960200',                      
                                 "Limimaricola" = '#721817',
                                 "Clade Ia"  =   '#CD7F78',  "Clade Ib"   = '#B84A62',                     
                                 "Erythrobacter"  = '#8C000A',  "Sphingobium  " = '#70161E',          
                                 "Glaciecola" =  '#A63B00', 
                                 "Vibrio"    = '#F2AC5D',   "Serratia"  = '#FFA200',  
                                 "Alcanivorax"    =  '#A05C00',                    
                                 "OM60(NOR5) clade"  = '#F35900', 
                                 "Marinobacter" =  '#DE6931', 
                                 "Acinetobacter" = '#FF8E00', "Psychrobacter"  = '#E55934', 
                                 "Lentimonas"  = '#fcca46',   
                                 "Roseibacillus" =  '#005c69')

# labels -----
labs_occurrence <- as_labeller(c(narrow = "Narrow\n(<1/3)",
                                 intermediate = "Intermediate\n(1/3 < x < 2/3)",
                                 broad = "Broad\n(>2/3)"))

labs_fraction_rec_freq <-  as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                                         '3' = 'Particle attached (3-20 um)',
                                         no = 'Recurrent',
                                         yes = 'Non-Recurrent',
                                         "non-recurrent" = 'Chaotic',
                                         recurrent = 'Seasonal',
                                         seasonal = 'Seasonal',
                                         stochastic = 'Chaotic'))

labs_perturbation <- as_labeller(c( 'pre_perturbation' = 'Pre' ,
                                    'perturbation' = 'During',
                                    'post_perturbation' = 'Post'))

labs_season <- as_labeller(c( 'winter' = 'Winter' ,
                              'spring' = 'Spring',
                              'summer' = 'Summer',
                              'autumn' = 'Autumn')) 

labs_fraction <- as_labeller(c('0.2' = 'Free living\n(0.2-3 um)',
                                   '3' = 'Particle attached\n(3-20 um)',
                               'bloomer' = 'Bloomer',
                               'no-bloomer' = 'No Bloomer'))

labs_diversity <- as_labeller(c('genes' = 'Bray Curtis Genes',
                                "wunifrac_distance" = 'Weigthed\nUNIFRAC',
                                "euclidean_distance" = 'Euclidean Environmental\nDistance',
                                "bray_curtis_community" = 'Bray Curtis\nCommunity Composition',
                                'bray_curtis_kmers' = 'Bray Curtis k-mers'))

labs_fraction_env <- as_labeller(c('0.2' = 'Free living\n(0.2-3 um)',
                                   '3' = 'Particle attached\n(3-20 um)',
                                   'env' = 'Environmental\nvariables',
                                   'bio' =  'Biological variables',
                                   'phch' = 'Pysicochemical variables'
                                   ))

### The remodelation of the Blanes harbour strated on 24th March 2010 and finished on the 9th of june 2012
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

# ---------------------- METHODS ----------------------  ########## ------
# ------ ########## Figure define bloomers threshold ------  ########## ------
## Calculate relative abundance ----
asv_tab_10y_l_rel <- asv_tab_bbmo_10y_l |>
  calculate_rel_abund(group_cols = sample_id)

asv_tab_10y_3_rel <- asv_tab_10y_l_rel %>%
  dplyr::filter(sample_id %in% m_3$sample_id)

asv_tab_10y_02_rel <- asv_tab_10y_l_rel %>%
  dplyr::filter(sample_id %in% m_02$sample_id)

asv_tab_10y_3_rel |>
  dplyr::group_by(sample_id) |>
  dplyr::reframe(sum_rel = sum(relative_abundance)) |>
  dplyr::filter(!sum_rel == 1) ### check that the sum of the relative abundances is 1

asv_tab_10y_02_rel |>
  dplyr::group_by(sample_id) |>
  dplyr::reframe(sum_rel = sum(relative_abundance)) |>
  dplyr::filter(!sum_rel == 1) ### check that the sum of the relative abundances is 1

asv_tab_10y_3_pseudo <- asv_tab_10y_3_rel |>
  calculate_pseudoabund(abund_data = m_bbmo_10y, rel_abund = relative_abundance, 
                        total_abund = bacteria_joint, 
                        by_ = 'sample_id')

asv_tab_10y_02_pseudo <- asv_tab_10y_02_rel |>
  calculate_pseudoabund(abund_data = m_bbmo_10y, 
                        rel_abund = relative_abundance, 
                        total_abund = bacteria_joint, 
                        by_ = 'sample_id')

source('src/count_number_potential_bloomers_threshold.R')

# I apply a loop for all the thresholds that I'm interested in
# Define a vector of threshold values
threshold_values <- c(0, 0.0001, 0.00015, 0.00025, 0.0005, 0.00075, 
                      0.001, 0.0015, 0.0025, 0.005, 0.0075, 
                      0.01, 0.015, 0.0, 0.025, 0.05, 0.075,
                      0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6)

# Create an empty list to store the results
datasets <- list()

# Iterate over each threshold value and apply the function
for (threshold in threshold_values) {
  dataset_name <- paste0("n_0.2_", threshold * 100)  # Create dataset name
  dataset <- count_num_bloomers(threshold, 0.2, asv_tab_10y_02_pseudo, z_scores_tb = z_scores_02_red)
  datasets[[dataset_name]] <- dataset  # Store the dataset in the list
}

# Combine all datasets into a single dataframe
result_dataset_02 <- bind_rows(datasets)

# Create an empty list to store the results
datasets <- list()

# Iterate over each threshold value and apply the function
for (threshold in threshold_values) {
  dataset_name <- paste0("n_3_", threshold * 100)  # Create dataset name
  dataset <- count_num_bloomers(threshold, 3, asv_tab_10y_3_pseudo, z_scores_tb = z_scores_3_red)
  datasets[[dataset_name]] <- dataset  # Store the dataset in the list
}

result_dataset_3 <- bind_rows(datasets)

### plot the results
blooming_threshold <- bind_rows(result_dataset_02,
                                result_dataset_3) |>
  ggplot(aes(threshold, num))+
  geom_point(size = 1)+
  scale_x_continuous(#expand = c(0,0),
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = percent_format())+
  #scale_y_continuous(expand = c(0,0))+
  #scale_y_log10()+
  facet_grid(vars(fraction), labeller = labs_fraction)+
  geom_vline(xintercept = 0.1, linetype = 'dashed')+
  labs(x = 'Relative abundance (%) threshold of the potential blooming event', y = 'Number of potential bloomers detected')+
  geom_line()+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 14),
        axis.ticks = element_blank())

blooming_threshold

# ggsave(blooming_threshold,  filename = 'blooming_threshold_ed2.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 188, units = 'mm')

# ---------------------- RESULTS ----------------------  ########## ------
# ------ ########## Figure recurrency in Blanes vs. changes over the years ------ #########
bray_unifrac_eucl_tb$bray_curtis_type |>
  unique()

data <- bray_unifrac_eucl_tb |>
  dplyr::filter(!bray_curtis_type %in% c('bray_curtis_kmers', "genes")) 

data$bray_curtis_type <- factor(data$bray_curtis_type, levels = c('bray_curtis_community', 'wunifrac_distance', 'euclidean_distance'))
data$fraction |>
  unique()

data$fraction <- factor(data$fraction, levels = c('0.2', '3', 'env'))

data <- data |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

# legend 
legend_plot <- data |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  ggplot(aes(date, bray_curtis_result))+
  facet_wrap(fraction ~ bray_curtis_type, scales = "free_x", ncol = 1, 
             labeller = labeller(fraction = function(x) ifelse(x == "fraction", "", x), 
                                 bray_curtis_type = labs_diversity), 
             strip.position = 'left') + 
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
  scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(data = data |>
              dplyr::filter(fraction == '3'), aes(date, bray_curtis_result), alpha = 0.3)+
  labs(x = 'Date', y = '', color = '', linetype = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    panel.grid.major.y = element_blank(),
    legend.text = element_text(size = 7),
    panel.border = element_blank(),
    plot.margin = margin(t = 5, b = 5))
legend_plot 

legend_plot <- get_legend(legend_plot)

data$bray_curtis_type |>
  unique()

bray_unifrac_eucl_plot <- data |>
  dplyr::filter(fraction != c('env')) |>
  dplyr::filter(bray_curtis_type !=  'wunifrac_distance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  ggplot(aes(date, bray_curtis_result))+
  facet_wrap(fraction ~ bray_curtis_type, scales = "fixed", ncol = 1, 
             labeller = labeller(fraction = function(x) ifelse(x %in% c("0.2", '3'), "", x), 
                                 bray_curtis_type = labs_diversity), 
             strip.position = 'left') + 
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
  scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(data = data |>
              dplyr::filter(fraction == '3' &
                              bray_curtis_type !=  'wunifrac_distance'), aes(date, bray_curtis_result), alpha = 0.3)+
  labs(x = 'Date', y = '', color = '', linetype = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside',
    axis.ticks.length.y =  unit(0.2, "mm"),
    plot.margin = margin(l = 0, t = 5, b = 5, r = 20))

bray_unifrac_eucl_plot

# eucl_plot <- data |>
#   dplyr::filter(fraction == 'env',
#                 bray_curtis_type == 'euclidean_distance') |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   dplyr::filter(!str_detect(date, '2003-')) |>
#   ggplot(aes(date, bray_curtis_result))+
#   facet_wrap(fraction ~ bray_curtis_type, scales = "free_x", ncol = 1, 
#              labeller = labeller(fraction = function(x) ifelse(x == "env", "", x), 
#                                  bray_curtis_type = labs_diversity), 
#              strip.position = 'left') + 
#   geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
#   scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
#   scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
#   scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
#   labs(x = 'Date', y = '', color = '', linetype = '')+
#   scale_shape_discrete(labels = labs_fraction)+
#   guides(shape = 'none',
#          color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
#   theme_bw()+
#   theme(#panel.grid = element_blank(), 
#     strip.background = element_blank(), legend.position = 'none',
#     panel.grid.minor = element_blank(),
#     axis.title  = element_text(size = 7),
#     strip.text = element_text(size = 5),
#     axis.text = element_text(size = 5),
#     # axis.text.x = element_text(size = 7), 
#     panel.grid.major.y = element_blank(),
#     panel.border = element_blank(),
#     strip.placement = 'outside', 
#     plot.margin = margin(t = 5, l = 10))
# 
# eucl_plot

eucl_plot <- euclidean_distance_tb_bio_phch |>
  ggplot(aes(date, euclidean_distance))+
  geom_line(data = euclidean_distance_tb_bio_phch |>
              dplyr::filter(type == 'phch'), aes(date, group = type, color = type), linewidth = 0.75, alpha = 0.3)+
  geom_line(aes(date, group = type, color = type, linetype = type), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
  #scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Date', y = 'Euclidean Distance', color = '', linetype = '')+
  #scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 7),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.placement = 'outside', 
    plot.margin = margin(t = 20, l = 46, r = 20))

eucl_plot

bray_curtis_lag_plot <- bray_curtis_lag_all |>
  ggplot(aes(lag, bray_curtis_result))+
  scale_x_continuous(breaks = c(12, 24, 36, 48, 12*5, 12*6, 12*7, 12*8, 12*9, 120), labels = c(1, 2, 3,4, 5, 6, 7, 8,  9, 10))+
  geom_point(aes(shape = fraction), alpha = 0.02)+
  scale_shape_discrete(labels = labs_fraction)+
  labs(y = 'Bray Curtis Dissimilarity', x = 'Lag between samples', shape = '', color = '', fill = '', linetype = '')+
  geom_line(data = bray_curtis_lag_all |>
              group_by(fraction, lag) |>
              dplyr::reframe(mean = mean(bray_curtis_result)), aes(lag, mean, group = fraction, color = fraction, 
                                                                   linetype = fraction), linewidth = 1)+
  scale_linetype_discrete( labels = labs_fraction)+
  scale_color_manual(values = palette_fraction_env, labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction_env, labels = labs_fraction)+
  guides(shape = guide_legend(ncol = 2))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        strip.background = element_rect(fill = 'transparent'),
        legend.position = 'none',
        plot.margin = margin(l = 30, t = 5, b = 5, r = 20))

bray_curtis_lag_plot 

# Now arrange the full layout, with bray_unifrac_eucl_plot occupying the top row
BBMO_community_diversity_presentation_plot <- plot_grid(
  bray_curtis_lag_plot,
  bray_unifrac_eucl_plot,
  eucl_plot,
  legend_plot,
  ncol = 1,                # One column layout for the main grid
  rel_heights = c(1.2, 1.5, 1, 0.25),
  labels = c('A', 'B', 'c')
)

# Print the final plot
print(BBMO_community_diversity_presentation_plot)

# ggsave( plot = BBMO_community_diversity_presentation_plot,
#         filename = 'BBMO_community_diversity_presentation_plot_v5.pdf',
#         path = 'results/figures/',
#         width = 180, height = 150, units = 'mm')

# ------ ########## Figure bloomers community timeseries ------ ########## ----------

## Rarefied dataset to calculate Community Evenness----
source('../../Bloomers/R/community_evenness.R')

community_eveness_02 <- asv_tab_bbmo_10y_w_rar |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou'))

community_eveness_3 <- asv_tab_bbmo_10y_w_rar |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou'))

### do not share any species.
source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R') #we need to upload it since we updated the function but I didn't compile the package

### We need the rarefied table transformed to relative abundances
asv_tab_bbmo_10y_l_rel_rar <- asv_tab_bbmo_10y_w_rar |>
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

community_eveness_all <- community_eveness_02 |>
  bind_rows(community_eveness_3) 

bray_curtis_rar_all <- bray_curtis_02_rar |> ##one sample less, the first one can't be compared with the previous
  bind_rows(bray_curtis_3_rar)

community_eveness_all_m <- community_eveness_all |>
  left_join(m_bbmo_10y, by = c('sample_id')) 

##color by order
## I remove the 4 ASVs that clustered together in the seasonality analysis but not others that could be potential blooms
community_eveness_all_m <- community_eveness_all_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

## reorder taxonomy as factors ----
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

## FIG 2
bbmo_bloo_ev_order_plot <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num %in% bloo_02_tb$asv_num & fraction == '0.2' |
                  asv_num %in% bloo_3_tb$asv_num & fraction == '3' ) |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, order_f) |>
  dplyr::mutate(abund_order = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
  )+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_area(aes(date, abund_order, fill = order_f, group = order_f), alpha = 1,  position='stack')+
  geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1),
                     sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
  scale_color_identity()+
  scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_env)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

bbmo_bloo_ev_order_plot

order_legend <- get_legend(bbmo_bloo_ev_order_plot)

bbmo_bloo_ev_order_plot <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(order_f = factor(order_f, levels = rev(unique(order_f)))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num %in% bloo_02_tb$asv_num & fraction == '0.2' |
                  asv_num %in% bloo_3_tb$asv_num & fraction == '3' ) |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, order_f) |>
  dplyr::mutate(abund_order = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
  )+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_area(aes(date, abund_order, fill = fct_rev(order_f), group = order_f), alpha = 1,
            position='stack')+
  geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1),
                     sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
  scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'none', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

bbmo_bloo_ev_order_plot

# figure B
bloo_all_types_summary_tax$recurrency <- factor(bloo_all_types_summary_tax$recurrency , 
                                                   levels = c('recurrent', 'non-recurrent'))
data_text <- bloo_all_types_summary_tax |>
  group_by(recurrency, fraction, asv_num) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, recurrency) |>
  dplyr::reframe(n_general = paste0( 'n = ', n())) 

tax_seasonality_plot <- bloo_all_types_summary_tax |>
  group_by(recurrency, fraction, order) |>
  dplyr::reframe(n = n()) |>
  group_by( fraction, recurrency) |>
  dplyr::mutate(n_total = sum(n)) |>
  dplyr::mutate(perc = n/n_total) |>
  left_join(data_text) |>
  ggplot(aes(perc,as.factor(fraction), fill = order))+
  scale_y_discrete(labels = labs_fraction)+
  scale_x_continuous(labels = percent_format())+
  geom_col()+
  geom_text(aes(label = n_general), check_overlap = TRUE, size = 2, vjust = 2, hjust = 0.65) +
  labs(x='', y = '')+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(recurrency), labeller = labs_fraction_rec_freq)+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 8),
        plot.margin = margin(2,20,0,5))

tax_seasonality_plot 

bbmo_bloo_ev_order_no_sar_cluster_seas_plot  <- plot_grid(bbmo_bloo_ev_order_plot, 
                                                             tax_seasonality_plot,
                                                             order_legend,
             ncol = 1, rel_heights = c(1, 0.35, 0.5))

bbmo_bloo_ev_order_no_sar_cluster_seas_plot

# ggsave('bbmo_bloo_ev_order_no_sar_cluster_seas_plot_v2.pdf', bbmo_bloo_ev_order_no_sar_cluster_seas_plot,
#        path = "results/figures/",
#        width = 180,
#        height = 160,
#        units = 'mm')

## BC and Community Evenness during bloom events -----
m_bbmo_10y <- m_bbmo_10y |>
  dplyr::mutate(date = as.character(as.Date(date)))

bloom_event_highlight <- bloom_events_tb |>  
  dplyr::mutate(fraction = as.character(fraction)) |>
  dplyr::mutate(fraction = as.numeric(fraction)) |>
  left_join(bloo_all_types_summary_tax, by = c('asv_num', 'fraction')) |>
  dplyr::mutate(bloom_event = 'bloom') |>
  dplyr::select(date, fraction, bloom_event) |>
  distinct() |>
  dplyr::mutate(fraction = as.character(fraction)) |>
  left_join(m_bbmo_10y) |>
  dplyr::select(decimal_date, fraction, bloom_event)

data <- bray_curtis_rar_all |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  rename(sample_id = samples) |>
  dplyr::select(sample_id, fraction, bray_curtis_result, decimal_date) |>
  left_join(community_eveness_all_m) |>
  dplyr::select(fraction, community_eveness_rar, bray_curtis_result, decimal_date) |>
  left_join(bloom_event_highlight)

## stats 
## check normality 
shapiro.test(as.numeric(data |>
                          dplyr::filter(bloom_event == 'bloom') %$%
                          bray_curtis_result)) # => p-value = 0.08601 (NORMALITY)
shapiro.test(as.numeric(data |>
                          dplyr::filter(is.na(bloom_event)) %$%
                          bray_curtis_result)) # => p-value = 5.796e-05 (NO NORMALITY)

shapiro.test(as.numeric(data |>
                          dplyr::filter(bloom_event == 'bloom') %$%
                          community_eveness_rar)) # => p-value = 3.386e-05 ( NO NORMALITY)

ggqqplot(as.numeric(data |>
                      dplyr::filter(bloom_event == 'bloom') %$%
                      bray_curtis_result))

data <- data |>
  dplyr::mutate(bloom_event = case_when(is.na(bloom_event) ~ 'No Bloom',
                                        TRUE ~ bloom_event)) 

## check homocedasticity 
leveneTest(as.numeric(bray_curtis_result) ~ as.factor(bloom_event), data = data) #p-value > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met
leveneTest(as.numeric(community_eveness_rar) ~ as.factor(bloom_event), data = data) # p-value > 0.05

# Perform one-way ANOVA to compare the means of Rho across the three groups
## anova_result <- aov(rho ~ fraction_causal, data = data_ed3_f)  # Replace `bloom` with your actual group variable
kruskal_result <- kruskal.test(as.numeric(bray_curtis_result) ~ as.factor(bloom_event), data = data) # significant differences
kruskal_result <- kruskal.test(as.numeric(community_eveness_rar) ~ as.factor(bloom_event), data = data) # significant differences

boxplot_bc_bloom <- data |>
  dplyr::mutate(bloom_event = case_when(is.na(bloom_event) ~ 'No Bloom',
                                        TRUE ~ bloom_event)) |>
  ggplot(aes(bloom_event, bray_curtis_result))+
  geom_point(alpha = 0.6, position = position_jitter(width = 0.1))+
  geom_boxplot(notch = T)+
  scale_shape_identity()+
  scale_x_discrete(labels = c( 'bloom' = 'Bloom', 'No Bloom' = 'No Bloom'))+
  labs(y = 'Bray Curtis Dissimiliarity', x = '')+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()+
  theme(aspect.ratio = 4/4,
        axis.ticks.length.y =  unit(0.2, "mm")
  )

boxplot_bc_bloom 

boxplot_ev_bloom <- data |>
  dplyr::mutate(bloom_event = case_when(is.na(bloom_event) ~ 'No Bloom',
                                        TRUE ~ bloom_event)) |>
  ggplot(aes(bloom_event, community_eveness_rar))+
  geom_point(alpha = 0.6, position = position_jitter(width = 0.1))+
  geom_boxplot(notch = T)+
  scale_shape_identity()+
  scale_x_discrete(labels = c( 'bloom' = 'Bloom', 'No Bloom' = 'No Bloom'))+
  labs(y = 'Community Evenness', x = '')+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()+
  theme(aspect.ratio = 4/4,
        axis.ticks.length.y =  unit(0.2, "mm")
  )

boxplot_ev_bloom

boxplot_bloom_community <- plot_grid(boxplot_bc_bloom, boxplot_ev_bloom,
                                     ncol = 1)

# ggsave('boxplot_bc_ev_bloom.pdf', boxplot_bc_bloom,
#        path = "results/figures/",
#        width = 80,
#        height = 120,
#        units = 'mm')

bbmo_bloo_ev_order_seas_bloom_events_plot <- plot_grid(bbmo_bloo_ev_order_no_sar_cluster_seas_plot, 
          boxplot_bloom_community, ncol = 2,
          rel_widths = c(1, 0.3))

# ggsave('bbmo_bloo_ev_order_seas_bloom_events_plot.pdf', bbmo_bloo_ev_order_seas_bloom_events_plot,
#        path = "results/figures/",
#        width = 180,
#        height = 160,
#        units = 'mm')

## -------- ########## Supplementary figure: BBMO bloomers time series separated by the different types of bloomers ------ ########## ---------
##try to explain that seasonal blooms are not always lead by the same ASV
labs_fraction_rec_freq <-  as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                                         '3' = 'Particle attached (3-20 um)',
                                         no = 'Recurrent',
                                         yes = 'Non-recurrent',
                                         seasonal = 'Seasonal',
                                         stochastic = 'Non-seasonal'))

summary_types_of_blooms <- bloo_all_types_summary_tb_tax_v2 |>
  dplyr::mutate(fraction = as.character(fraction))

bloo_all_types_summary_tb <- bloo_all_types_summary_tb_tax_v2 |>
  dplyr::mutate(fraction = as.character(fraction))

bloo_all_types_summary_tb$recurrency <- factor(bloo_all_types_summary_tb$recurrency, levels = c('recurrent', 'non-recurrent'))

summary_types_of_blooms |>
  dplyr::filter(recurrency == 'recurrent' &
                  type_of_bloomer == 'Seasonal')

booming_events_seas_BBMO10Y <- asv_tab_all_bloo_z_tax |>
  left_join(bloo_all_types_summary_tb, by = c('asv_num', 'fraction')) |>
  dplyr::filter(asv_num %in% bloo_02$value & fraction == '0.2' |
                  asv_num %in% bloo_3$value & fraction == '3' ) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(date, fraction, recurrency) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, order_f, recurrency) |>
  dplyr::mutate(abund_class = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
  )+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), 
            fill = '#C7C7C7', alpha = 0.5)+
  

  geom_area(aes(date, abund_class, fill = fct_rev(order_f), group = fct_rev(order_f)), alpha = 1,  position='stack')+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
  scale_color_identity()+
  scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Order')+
  facet_grid(fraction~recurrency,  scales = 'free_y',  labeller =  labs_fraction_rec_freq)+
  guides(fill = guide_legend(ncol = 5, size = 7,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text.x = element_text(size = 14),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        strip.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), strip.background = element_blank(), 
        legend.text = element_text(size = 8), legend.title = element_text(size = 10), 
        strip.placement = 'outside', 
        panel.spacing.y = unit(1, "lines"), panel.border = element_blank()) 

booming_events_seas_BBMO10Y

# ggsave('booming_events_BBMO10Y_new_tax_seasonal_ed4.pdf', booming_events_seas_BBMO10Y,
#        path = "results/figures/",
#        width = 180,
#        height = 160,
#        units = 'mm')

# --------- ########## Figure example of different types of bloomers identified in our dataset ------ ########## ------------------------
bloo_all_types_summary_tb_tax_v2 |>
  dplyr::filter(recurrency == 'recurrent')

bloo_all_types_summary_tb_tax_v2 |>
  dplyr::filter(recurrency != 'recurrent') |>
  dplyr::filter(asv_num %in% c('asv11', 'asv555'))

bloo_all_types_summary_tb_tax_v2 |>
  dplyr::filter((asv_num %in% c('asv17', 'asv11',  'asv72') & fraction == '3') |
                  (asv_num %in% c('asv7',  'asv555','asv62') & fraction == '0.2') |
                  (asv_num == 'asv11' & fraction == '0.2'))

asv_tab_all_bloo_z_tax_examples <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(fraction = as.double(fraction)) |>
  dplyr::select(-order,-family, -phylum, -class) |>
  dplyr::left_join(bloo_all_types_summary_tb_tax, by = c( 'asv_num','fraction')) |>
  dplyr::filter(order != 'SAR11 clade') |>
  #dplyr::filter(asv_num_f %in% c('asv17', 'asv23', 'asv11', 'asv84', 'asv62', 'asv58', 'asv555', 'asv559')) |>
  dplyr::filter((asv_num %in% c('asv17', 'asv11',  'asv72') & fraction == '3') |
                  (asv_num %in% c('asv7',  'asv555','asv62') & fraction == '0.2') |
                  (asv_num == 'asv11' & fraction == '0.2')) |>
  #dplyr::left_join(occurrence_bloo_bbmo) |>
  dplyr::mutate(bloom = case_when (z_score_ra >= 1.96 &
                                     abundance_value >= 0.1~ 'Bloom',
                                   TRUE ~ 'No-bloom')) |>
  dplyr::mutate(bloom = as.factor(bloom)) |>
  dplyr::mutate(asv_num_f = asv_num)

asv_summary_bloom <- asv_tab_all_bloo_z_tax_examples |>
  distinct(asv_num, occurrence_category, recurrency, frequency) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num'))

##plot1 
plot1 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv7')) |>
  ggplot(aes(date, abundance_value))+
  geom_hline(yintercept = 0.1, linetype = 'dashed', linewidth = 0.3)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.3))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),  position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv7'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+ #, title = 'Recurrents'
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 0),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 7), strip.background = element_rect(fill = "#57a9a8",color = "#57a9a8" ), #strip.text.x = element_blank(),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6, title = element_text(size = 8),
        panel.border = element_blank())

##plot2
plot2 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv62')) |>
  ggplot(aes(date, abundance_value))+
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.2))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv62'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 0),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 7), strip.background = element_rect(fill = "#3e3e3e",color = "#3e3e3e"), #strip.text.x = element_blank(),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6, title = element_text(size = 8),
        panel.border = element_blank())
# plot 3
plot3 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv72')) |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0),  limits = c(0,0.3))+ #
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv72'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 0),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 7), strip.background = element_rect(fill = "#AE659B", color = "#AE659B"), #strip.text.x = element_blank(),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6, title = element_text(size = 8),
        panel.border = element_blank())

# plot 4
plot4 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv17')) |>
  ggplot(aes(date, abundance_value))+
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.35))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv17'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+ #, title = 'Non-Recurrents'
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#57a9a8", color = "#57a9a8" ),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6,
        panel.border = element_blank())

# plot 7
plot7 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv11')) |>
  dplyr::filter(fraction == '0.2') |>
  ggplot(aes(date, abundance_value))+
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.5))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv11',
                             fraction == '0.2'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+ #, title = 'Ephemeral'
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#3e3e3e", color = "#3e3e3e"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6,
        panel.border = element_blank())

# plot 8
plot8 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv555')) |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.2))+ #, limits = c(0,0.25)
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  geom_area( aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv555'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#AE659B", color = "#AE659B"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6,
        panel.border = element_blank())

# legend
legend_plot <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(bloom = case_when (z_score_ra >= 1.96 &
                                     abundance_value >= 0.1~ 'Bloom',
                                   TRUE ~ 'No-bloom')) |>
  dplyr::mutate(bloom = as.factor(bloom)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, abundance_value)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, 0.2)) +
  geom_area(aes(group = asv_num_f, fill = occurrence_category), position = 'identity') +
  geom_point(aes(color = bloom), size = 1) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_fill_manual(values = palette_occurrence, labels = labs_occurrence) +
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Occurrence category', color = 'Bloom')+
  guides(fill = guide_legend(ncol = 3, size = 6),
         color = guide_legend(ncol = 2, size = 6),  # Hide points in the color legend
         alpha = 'none') +
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(4, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#AE659B", color = "transparent"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'))

#### I decide to keep only 6 blooming examples ------
legend_only  <- cowplot::get_legend(legend_plot)

#library(cowplot)
#library(gridExtra)

# Define the layout matrix
layout_mat <- rbind(c(1, 2, 3),
                    c(4, 5, 6),
                    c(8,8,8),
                    c(7, 7, 7),# This specifies that legend_only occupies all columns in the last row
                    c(9,9,9))  

# Define heights for each row
heights <- c(0.95, 0.95,1.95, 0.15,  0.30)  # Adjust the height of the last row to be shorter

# Arrange the plots using the layout matrix and heights
examples_bloomers_types_titles_observed <- grid.arrange(plot1, plot2, plot3,
                                               plot4, plot7, plot8,
                                               types_of_bloomers_summary_plot, 
                                               legend_only,
                                               legend_order,
                                               layout_matrix = layout_mat,
                                               heights = heights)
# Add text to each row
#grid.text("Recurrent", x = unit(0.02, "npc"), y = unit(1, "npc") - unit(0.5 * heights[1], "cm"), just = "left", gp = gpar(fontsize = 8))
#grid.text("Non-recurrent", x = unit(0.02, "npc"), y = unit(1, "npc") - unit(heights[1] + 3 * heights[2], "cm"), just = "left", gp = gpar(fontsize = 8))
# Add grid line between rows 1 and 2 (still tpe be improved)
#grid.lines(x = unit(c(0, 1), "npc"), y = unit(1, "npc") - unit(heights[1], "npc"), gp = gpar(col = "black", lwd = 1))

#Save the combined plot (trobar la manera de guardar-ho amb el text!!)
ggsave(filename = 'examples_bloomers_types_titles_observed.pdf', plot = examples_bloomers_types_titles_observed,
       path = 'results/figures/',
       width = 188, height = 230, units = 'mm')

# ggsave(filename = 'examples_bloomers_types_red.svg', 
#        plot = examples_bloomers_types_titles,
#        path = 'results/figures/poster_svg_format/',
#        width = 188, height = 120, units = 'mm')

## same plot without colors ----
types_of_bloomers_summary_plot <- wavelet_results_all  |>
  dplyr::filter(!is.na(type_of_bloomer)) |>
  left_join(text_data) |>
  dplyr::filter(type_of_bloomer %in% c('Seasonal', 'No-significant Periodicity')) |>
  dplyr::mutate(asv_fraction = paste0(asv_num, '.', fraction)) |>
  #dplyr::filter(type_of_bloomer == 'Year-to-Year') |>
  ggplot(aes(period, power_average, group = asv_num))+
  geom_vline(xintercept = 4,  alpha = 0.1)+
  geom_vline(xintercept = 8,  alpha = 0.1)+
  geom_vline(xintercept = 16,  alpha = 0.1)+
  geom_vline(xintercept = 32,  alpha = 0.1)+
  geom_text( aes(label =  n, x = 30, y = 0.75), color = 'black', 
             check_overlap = TRUE, size = 3) + 
  geom_point(data = wavelet_results_all |>
               dplyr::filter(pvalue_power_average < 0.05 &
                               power_average > mean_power_sig$mean &
                               type_of_bloomer %in% c('Seasonal', 'No-significant Periodicity')))+
  # alpha = ifelse(pvalue_power_average < 0.05 &
  #                               power_average > mean_power_sig$mean, 1, 0), color = order))+
  geom_line(aes(group = asv_fraction, linetype = fraction), alpha = 1)+
  #scale_linetype_discrete(values = c('dashed' = '0.2', 'solid' = '3'))+
  # geom_line(aes( group = fraction,  linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  scale_y_continuous(limits = c(0, 0.95))+
  coord_flip()+
  #facet_wrap(vars(asv_num))+
  facet_grid(type_of_bloomer~occurrence_category, switch = 'y')+
  scale_x_continuous(limits = c(2, 33),
                     breaks = c(3, 6, 12 , 24 , 33),
                     labels = c("Short\n(<4 months)", "Half-yearly\n(4-8 months)", "Seasonal\n(8-16 months)", 
                                "Year-to-year\n(>16 months)", 'Inter-Annual(>32)'))+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y = 'Wavelet Power Average', x = 'Period (Months)', color = 'Order', linetype = 'Fraction')+
  theme_bw()+
  theme(legend.position = 'bottom', panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect('transparent'),
        strip.text = element_text(size = 6),
        text = element_text(size = 6),
        panel.spacing.x = unit(2, "lines"),  # adjust the spacing between panels in the x direction
        panel.spacing.y = unit(1, "lines"))+  # adjust the spacing between panels in the y direction)+
  guides(alpha = 'none',
         color = 'none',
         linetype = 'none')

types_of_bloomers_summary_plot

legend_order <- wavelet_results_all  |>
  dplyr::filter(!is.na(type_of_bloomer)) |>
  left_join(text_data) |>
  dplyr::filter(type_of_bloomer %in% c('Seasonal', 'No-significant Periodicity')) |>
  dplyr::mutate(asv_fraction = paste0(asv_num, '.', fraction)) |>
  #dplyr::filter(type_of_bloomer == 'Year-to-Year') |>
  ggplot(aes(period, power_average, group = asv_num, color = order))+
  geom_vline(xintercept = 4,  alpha = 0.1)+
  geom_vline(xintercept = 8,  alpha = 0.1)+
  geom_vline(xintercept = 16,  alpha = 0.1)+
  geom_vline(xintercept = 32,  alpha = 0.1)+
  geom_text( aes(label =  n, x = 30, y = 0.75), color = 'black', 
             check_overlap = TRUE, size = 3) + 
  geom_point(data = wavelet_results_all |>
               dplyr::filter(pvalue_power_average < 0.05 &
                               power_average > mean_power_sig$mean &
                               type_of_bloomer %in% c('Seasonal', 'No-significant Periodicity')))+
  # alpha = ifelse(pvalue_power_average < 0.05 &
  #                               power_average > mean_power_sig$mean, 1, 0), color = order))+
  geom_line(aes(group = asv_fraction, linetype = fraction), alpha = 1)+
  #scale_linetype_discrete(values = c('dashed' = '0.2', 'solid' = '3'))+
  # geom_line(aes( group = fraction,  linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  scale_y_continuous(limits = c(0, 0.95))+
  coord_flip()+
  #facet_wrap(vars(asv_num))+
  facet_grid(type_of_bloomer~occurrence_category, switch = 'y')+
  scale_x_continuous(limits = c(2, 33),
                     breaks = c(3, 6, 12 , 24 , 33),
                     labels = c("Short\n(<4 months)", "Half-yearly\n(4-8 months)", "Seasonal\n(8-16 months)", 
                                "Year-to-year\n(>16 months)", 'Inter-Annual(>32)'))+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y = 'Wavelet Power Average', x = 'Period (Months)', color = 'Order', linetype = 'Fraction')+
  theme_bw()+
  theme(legend.position = 'bottom', panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect('transparent'),
        strip.text = element_text(size = 6),
        text = element_text(size = 6),
        panel.spacing.x = unit(2, "lines"),  # adjust the spacing between panels in the x direction
        panel.spacing.y = unit(1, "lines"))+  # adjust the spacing between panels in the y direction)+
  guides(alpha = 'none',
         linetype = guide_legend(ncol = 1),
         color = 'none')

legend_order <- get_legend( legend_order )

# Define the layout matrix
layout_mat <- rbind(c(1, 2, 3),
                    c(4, 5, 6),
                    c(7,7,7),
                    c(8, 8, 9))  

# Define heights for each row
heights <- c(0.95, 0.95,1.95, 0.55)  # Adjust the height of the last row to be shorter

# Arrange the plots using the layout matrix and heights
examples_bloomers_types_titles_observed <- grid.arrange(plot1, plot2, plot3,
                                                        plot4, plot7, plot8,
                                                        types_of_bloomers_summary_plot, 
                                                        legend_only,
                                                        legend_order,
                                                        layout_matrix = layout_mat,
                                                        heights = heights)

# ggsave(filename = 'examples_bloomers_types_titles_observed_ed1.pdf', plot = examples_bloomers_types_titles_observed,
#        path = 'results/figures/',
#        width = 188, height = 230, units = 'mm')



#  ------ ########## Figure seasonal bloomers do not bloom every year ------ ########## ------ 

bloo_all_types_summary_tb <- bloo_all_types_summary_tb_tax_v2 |>
  dplyr::mutate(fraction = as.factor(fraction))

data <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(fraction = as.factor(fraction)) |>
  dplyr::left_join(bloo_all_types_summary_tb, by = c('asv_num_f' = 'asv_num','fraction', 'order', 'class', 'phylum', 'family')) |> 
  dplyr::filter(type_of_bloomer == 'Seasonal' )

data <- data |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num))

data$class_f <-  factor(data$class_f, 
                        levels=unique(data$class_f[order(data$phylum_f)]), 
                        ordered=TRUE)

data$order_f <-  factor(data$order_f, 
                        levels=unique(data$order_f[order(data$phylum_f,
                                                         data$class_f)]), 
                        ordered=TRUE)

data$family_f <-  factor(data$family_f, 
                         levels=unique(data$family_f[order(data$phylum_f,
                                                           data$class_f,
                                                           data$order_f)]), 
                         ordered=TRUE)

data$asv_num_f <-  factor(data$asv_num_f, 
                          levels=unique(data$asv_num_f[order(data$phylum_f,
                                                             data$class_f,
                                                             data$order_f,
                                                             data$family_f)]), 
                          ordered=TRUE)

data <- data |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) 

data$facet_labels <- factor(data$facet_labels, 
                            levels = unique(data$facet_labels[order(
                              data$family_f,
                              data$asv_num_f)]), 
                            ordered=TRUE)

seasonal_bloo_y <- data |>
  ggplot(aes(day_of_year, abundance_value, color = family_f))+
  geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_point(aes(color = family_f), size = 0.2)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = fraction, linetype = fraction), linewidth = 1)+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype(labels = labs_fraction)+
  labs(linetype = 'Fraction', color = 'Family')+
  facet_grid(facet_labels ~ year, scales = 'free_y')+
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance' )+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text.x = element_text(size = 7, margin = margin(5, 5, 5, 5)),
        strip.text.y = element_text(size = 5, margin = margin(1, 1, 1, 1)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 10), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))+
  guides(linetype = guide_legend(ncol = 1))

seasonal_bloo_y

# ggsave(filename = 'results/figures/seasonal_bloo_y_fraction_ed2.pdf', plot = seasonal_bloo_y,
#        # path = 'results/figures/',
#        width = 180, height = 200, units = 'mm')
# ------ Figure harbor restoration ------ ########## ------
## are there statistically significant differences? -----
bloo_all_types_summary_tax_3 <- bloo_all_types_summary_tax |>
  dplyr::filter(fraction == '3') |>
  dplyr::mutate(fraction = as.character(fraction)) 

positive_negative_effect <- asv_tab_all_bloo_z_tax |> 
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d")) |> 
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  )) |> 
  dplyr::group_by(harbor_restoration, abundance_type, asv_num) |> 
  dplyr::filter(fraction == '3') |> 
  dplyr::filter(abundance_type == 'rclr') |> 
  dplyr::reframe(mean_abund = mean(abundance_value)) |> 
  pivot_wider(names_from = 'harbor_restoration', values_from = starts_with('mean')) |> 
  dplyr::mutate(harbor_effect_rclr = case_when(
    pre_perturbation > perturbation ~ 'negative',
    pre_perturbation < perturbation ~ 'positive',
    after_perturbation > perturbation |
      after_perturbation > pre_perturbation ~ 'opportunistic',
    TRUE ~ 'neutral'
  )) |> 
  distinct(asv_num, harbor_effect_rclr)

data <- asv_tab_all_bloo_z_tax |> 
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d")) |> 
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  )) |> 
  dplyr::filter(fraction == '3') |> 
  dplyr::filter(abundance_type == 'rclr') 

positive_negative_effect_stats <- data |> 
  dplyr::group_by(asv_num) |>  # Group by ASV
  tidyr::nest() |>  # Nest the data for each ASV
  dplyr::mutate(
    # For each ASV, extract the nested data and calculate pre, perturbation, and after means
    pre_values = map(data, ~ .x$abundance_value[.x$harbor_restoration == "pre_perturbation"]),
    pert_values = map(data, ~ .x$abundance_value[.x$harbor_restoration == 'perturbation']),
    after_values = map(data, ~ .x$abundance_value[.x$harbor_restoration == 'after_perturbation']),
    
    # Perform the global test comparing pre, perturbation, and after
    test_result = pmap(list(pre_values, pert_values, after_values), ~ {
      pre_vals <- ..1
      pert_vals <- ..2
      after_vals <- ..3
      
      # Skip if any of the groups has fewer than 2 values
      if (length(pre_vals) > 1 & length(pert_vals) > 1 & length(after_vals) > 1) {
        # Check if all values in any group are identical (constant values)
        pre_identical = length(unique(pre_vals)) == 1
        pert_identical = length(unique(pert_vals)) == 1
        after_identical = length(unique(after_vals)) == 1
        
        # Perform normality test only if the values in the group are not constant
        pre_normal = if (!pre_identical) shapiro.test(pre_vals)$p.value > 0.05 else FALSE
        pert_normal = if (!pert_identical) shapiro.test(pert_vals)$p.value > 0.05 else FALSE
        after_normal = if (!after_identical) shapiro.test(after_vals)$p.value > 0.05 else FALSE
        
        # Choose between ANOVA or Kruskal-Wallis based on normality
        if (pre_normal & pert_normal & after_normal) {
          # Prepare data for ANOVA
          df = tibble(
            abundance_value = c(pre_vals, pert_vals, after_vals),
            harbor_restoration = rep(c("pre_perturbation", "perturbation", "after_perturbation"),
                                     times = c(length(pre_vals), length(pert_vals), length(after_vals)))
          )
          # Perform one-way ANOVA
          aov_res = aov(abundance_value ~ harbor_restoration, data = df)
          tidy(aov_res) |> dplyr::filter(term == "harbor_restoration") |> pull(p.value)
        } else {
          # Prepare data for Kruskal-Wallis test
          df = tibble(
            abundance_value = c(pre_vals, pert_vals, after_vals),
            harbor_restoration = rep(c("pre_perturbation", "perturbation", "after_perturbation"),
                                     times = c(length(pre_vals), length(pert_vals), length(after_vals)))
          )
          # Perform Kruskal-Wallis test
          kruskal_res = kruskal.test(abundance_value ~ harbor_restoration, data = df)
          kruskal_res$p.value
        }
      } else {
        NA_real_  # Skip if there is insufficient data
      }
    }),
    
    # Compute mean values for classification
    mean_pre = map_dbl(pre_values, ~ mean(.x, na.rm = TRUE)),
    mean_pert = map_dbl(pert_values, ~ mean(.x, na.rm = TRUE)),
    mean_after = map_dbl(after_values, ~ mean(.x, na.rm = TRUE)),
    
    # Classify harbor effect based on the mean comparisons
    harbor_effect = case_when(
      mean_pre > mean_pert ~ 'negative',
      mean_pre < mean_pert ~ 'positive',
      mean_after > mean_pert | mean_after > mean_pre ~ 'opportunistic',
      TRUE ~ 'neutral'
    ),
    
    # Mark if the global test is significant
    significant = ifelse(!is.na(test_result) & test_result < 0.05, TRUE, FALSE)
  ) |> 
  dplyr::select(asv_num, harbor_effect, test_result, significant) |>
  as_tibble()

positive_negative_effect_stats  |> 
  dplyr::filter(significant == 'TRUE') #32

positive_negative_effect |>
  left_join(positive_negative_effect_stats) |>
  dplyr::filter(harbor_effect_rclr != harbor_effect) ## none! 

clear_effect_asvs <- positive_negative_effect |>
  left_join(positive_negative_effect_stats) |>
  dplyr::filter(harbor_effect_rclr == harbor_effect) |>
  left_join(bloo_all_types_summary_tax_3)

clear_effect_asvs |>
  dplyr::filter(significant == 'TRUE') # 32 ASVs

## plot to check that stats are correct ----
data <- asv_tab_all_bloo_z_tax |> 
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d")) |> 
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  )) |> 
  dplyr::group_by(harbor_restoration, abundance_type, asv_num) |> 
  dplyr::filter(fraction == '3') |> 
  dplyr::filter(abundance_type == 'rclr') |>
  left_join(positive_negative_effect_stats) |>
  dplyr::mutate(asv_num = reorder(asv_num, significant == "Significant", FUN = max))

data$significant <- factor(data$significant, levels = c('TRUE', 'FALSE'), labels = c('Significant', 'Non-significant'))
data$harbor_restoration <- factor(data$harbor_restoration, levels = c('pre_perturbation', 'perturbation', 'after_perturbation'))

differential_abundance_harbor_plot <- data |>
  ggplot(aes(harbor_restoration, abundance_value))+
  geom_point(aes(color = significant, shape = significant), position = position_jitter(width = 0.25))+
  geom_boxplot(alpha = 0.2)+
  scale_color_manual(values = c('#4cb76a', '#2466BF'))+
  scale_x_discrete(labels = c('pre_perturbation' = 'Pre', 'perturbation' = 'Perturbation', 'after_perturbation' = 'Post'))+
  labs(x = '', y = 'rCLR', color = '', shape = '')+
  facet_wrap(vars(interaction(order_f, asv_num_f)), scales = 'free', ncol = 5)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

differential_abundance_harbor_plot

# ggsave(filename = 'differential_abundance_harbor_plot.pdf',
#        plot = differential_abundance_harbor_plot,
#        path = 'Results/Figures/',
#        width = 180, height = 230, units = 'mm')

asv_tab_all_bloo_z_tax_sig <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(fraction == '3') |>
  left_join(clear_effect_asvs) |>
  dplyr::mutate(type_of_bloomer = case_when(type_of_bloomer == 'Inter-Annual' ~ 'Chaotic',
                                            type_of_bloomer == "Year-to-Year" ~ 'Chaotic',
                                            type_of_bloomer ==  "No-significant Periodicity" ~ 'Chaotic',
                                            type_of_bloomer == 'Recurrent' ~ 'Seasonal',
                                            TRUE ~ type_of_bloomer))
### Plot with harbor restoration images ----
asv_tab_all_bloo_z_tax_sig$harbor_effect <- factor(asv_tab_all_bloo_z_tax_sig$harbor_effect, levels = c('negative', 'positive'))

asv_tab_all_bloo_z_tax_sig$type_of_bloomer <- factor(asv_tab_all_bloo_z_tax_sig$type_of_bloomer, levels = c('Seasonal', 'Chaotic'))

harbor_pa_stats_plot <- asv_tab_all_bloo_z_tax_sig |>
  dplyr::filter(significant == 'TRUE') |>
  dplyr::mutate(harbor_effect = case_when(harbor_effect == 'positive' ~ 'Positive',
                                          harbor_effect == 'negative' ~ 'Negative')) |>
  dplyr::group_by(date, fraction, order_f, family_f, abundance_type, type_of_bloomer, harbor_effect) |>
  dplyr::reframe(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abund_max))+
  geom_vline(xintercept = as.POSIXct('2010-03-24', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_vline(xintercept = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_vline(xintercept = as.POSIXct('2012-06-09', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), 
                                                        xmax = as.POSIXct('2010-08-31', format = "%Y-%m-%d"), x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2011-07-01', format = "%Y-%m-%d"), 
                                                        xmax = as.POSIXct('2011-08-31', format = "%Y-%m-%d"), x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(y = 'rCLR', x = 'Time', fill = 'Family')+
  geom_area(aes(group = fct_rev(family_f), fill = fct_rev(family_f)), position = 'stack')+
  facet_wrap(vars(interaction(type_of_bloomer, harbor_effect)),  scales = 'free_y', ncol = 1)+
  guides(fill = guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        panel.grid = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 0.5, 5.25, 0.5), "cm"),
        legend.key.size = unit(4, 'mm'))

composition <- plot_grid(harbor_pa_stats_plot, align = "hv", nrow = 1) + 
  draw_image("data/env_data/harbour_restoration_data/blanes_harbor_icgc/Blanes_2010.jpg", x = 0.05, y = 0, width = 0.25, height = 0.25) +
  draw_image("data/env_data/harbour_restoration_data/blanes_harbor_icgc/Blanes_2011.jpg", x = 0.33, y = 0, width = 0.25, height = 0.25) +
  draw_image("data/env_data/harbour_restoration_data/blanes_harbor_icgc/Blanes_2013.jpg", x = 0.61, y = 0, width = 0.25, height = 0.25)

# Add labels to each part of the composition
composition_with_labels <- composition +
  annotate("text", x = 0.03, y = 0.97, label = "A", size = 3) +
  annotate("text", x = 0.07, y = 0.17, label = "B", size = 3,  colour = "white") +
  annotate("text", x = 0.34, y = 0.17, label = "C", size = 3, colour = "white") +
  annotate("text", x = 0.62, y = 0.17, label = "D", size = 3, colour = "white")

# Display the composition with labels
print(composition_with_labels)

# Save the composition as a PDF
# ggsave(filename = 'harbor_restoration_plot_photo_rclr_stats_v5.pdf',
#        plot = composition_with_labels,
#        path = 'Results/Figures/',
#        width = 180, height = 200, units = 'mm')

### relative abundance (supplementary plot) ----
asv_tab_all_bloo_z_tax_sig <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '3') |>
  left_join(clear_effect_asvs) |>
  dplyr::mutate(type_of_bloomer = case_when(type_of_bloomer == 'Inter-Annual' ~ 'Chaotic',
                                            type_of_bloomer == "Year-to-Year" ~ 'Chaotic',
                                            type_of_bloomer ==  "No-significant Periodicity" ~ 'Chaotic',
                                            type_of_bloomer == 'Recurrent' ~ 'Seasonal',
                                            TRUE ~ type_of_bloomer))

asv_tab_all_bloo_z_tax_sig$harbor_effect <- factor(asv_tab_all_bloo_z_tax_sig$harbor_effect, levels = c('negative', 'positive'))

asv_tab_all_bloo_z_tax_sig$type_of_bloomer <- factor(asv_tab_all_bloo_z_tax_sig$type_of_bloomer, levels = c('Seasonal', 'Chaotic'))

harbor_pa_stats_plot <- asv_tab_all_bloo_z_tax_sig |>
  dplyr::filter(significant == 'TRUE') |>
  dplyr::mutate(harbor_effect = case_when(harbor_effect == 'positive' ~ 'Positive',
                                          harbor_effect == 'negative' ~ 'Negative')) |>
  dplyr::group_by(date, fraction, order_f, family_f, abundance_type, type_of_bloomer, harbor_effect) |>
  dplyr::reframe(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abund_max))+
  geom_vline(xintercept = as.POSIXct('2010-03-24', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_vline(xintercept = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_vline(xintercept = as.POSIXct('2012-06-09', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), 
                                                        xmax = as.POSIXct('2010-08-31', format = "%Y-%m-%d"), x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2011-07-01', format = "%Y-%m-%d"), 
                                                        xmax = as.POSIXct('2011-08-31', format = "%Y-%m-%d"), x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(y = 'Relative Abundance', x = 'Time', fill = 'Family')+
  geom_area(aes(group = fct_rev(family_f), fill = fct_rev(family_f)), position = 'stack')+
  facet_wrap(vars(interaction(type_of_bloomer, harbor_effect)),  scales = 'free_y', ncol = 1)+
  guides(fill = guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        panel.grid = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        legend.key.size = unit(4, 'mm'))

harbor_pa_stats_plot
#
# ggsave(filename = 'harbor_restoration_plot_photo_rel_abund_stats_v3.pdf',
#        plot = harbor_pa_stats_plot,
#        path = 'Results/Figures/',
#        width = 180, height = 200, units = 'mm')

## how was the diversity affected by the harbor restoration ----
asv_tab_bbmo_10y_w_rar_t <- asv_tab_bbmo_10y_w_rar |>
  t() |>
  as_tibble(rownames = NULL)

asv_tab_bbmo_10y_w_rar_t |>
  colnames() <- asv_tab_bbmo_10y_w_rar_t[1,]

asv_tab_bbmo_10y_w_rar_t <- asv_tab_bbmo_10y_w_rar_t[-1,]

shannon_rar <- asv_tab_bbmo_10y_w_rar_t |>
  dplyr:: mutate(across(everything(), as.numeric)) |>
  microbiome::alpha('shannon') |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble()

shannon_rar_m <- shannon_rar |>
  left_join(m_bbmo_10y) |>
  dplyr::mutate(harbor_restoration = case_when(  date < '2010-03-24' ~ 'pre_perturbation',
                                                 date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
                                                 date >= '2012-06-09' ~ 'post_perturbation'))

shannon_rar_m$harbor_restoration <- factor(shannon_rar_m$harbor_restoration, 
                                           levels = c('post_perturbation',
                                                      'perturbation',
                                                      'pre_perturbation'))

shannon_harbor_rar_plot <- shannon_rar_m |>
  ggplot(aes(fraction, diversity_shannon))+
  geom_point( position = position_jitter(width = 0.15), aes(shape = fraction))+#aes(color = season),
  scale_x_discrete(labels = labs_fraction)+
  labs(y = 'Shannon diversity', x = 'Harbor restoration', color = 'Season', shape = 'Fraction')+
  facet_wrap(vars(fct_rev(harbor_restoration)), labeller = labs_perturbation)+
  scale_shape_discrete(labels = labs_fraction)+
  geom_violin(aes(group = fraction), alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  theme_bw()+
  theme(text = element_text(size = 20),
        strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

shannon_harbor_rar_plot 

# ------ Figure Convergent Cross Mapping (CCM) results -------- # ------
# ------ ####### ------ EXPLORE EDM RESULTS ------ ####### ------
## CCM convergent cross mapping results computed in MARBITS (it's a bit computationally expensive) ----
### 1. DATASET: We removed seasonality & data is in the format of rCLR -------
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

# Set a threshold for the maximum number of allowed NA values
na_threshold <- 100 ## to filter by it i need to have a values lower than 85 which is the dimension of my matrix. I don't want to filter in this case.
## 

# Identify columns where the number of NA values is less than or equal to the threshold
cols_below_threshold <- apply(ccm_rho_filt , 2, function(col) sum(is.na(col)) <= na_threshold)

# Remove columns where the number of NA values exceeds the threshold
matrix_clean <- ccm_rho_filt [, cols_below_threshold]
matrix_clean <- ccm_rho_filt [cols_below_threshold, ]

#### ---- add stats for each group ----
### stats for environmental variables separated in physicochemical and biological (MAIN) --------
data_ed2 |>
  colnames()

select_bloomers <- data_ed2 |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::distinct(asv_num_predicted, frequency_predicted, fraction_predicted) |>
  dplyr::filter(frequency_predicted != 'No Bloomer')

### anova to compare mean of the different groups 
data_ed3 <- data_ed2 |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted)~ 'no_bloomer',
                                  asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::mutate(fraction_bloom = paste0(fraction_predicted,bloom))

# I need to perform an anova for each group in case that each group has a normal distribution and that they have homocedasticity
#### group I: 0.2 no bloomers  -----
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
# ------------------------------- #
#   3 |  -0.144353
# |     1.0000
# |
#   env |   15.66140   13.97643
# |    0.0000*    0.0000*

#### group II:  0.2 bloomers  ------
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

#### group III: 3 no bloomers----
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

#### group IV: 3 bloomers -----
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

# plot ccm results summarized by groups -----
## without color (non informative in this case)
data_ed2 <- data_ed |>
  dplyr::filter(frequency_predicted != 'env') |>
  dplyr::mutate(bloom = case_when(!asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'no-bloomer',
                                  asv_num_predicted %in%  unique(select_bloomers$asv_num_predicted) ~ 'bloomer')) |>
  dplyr::mutate(fraction_causal = case_when((fraction_causal == 'env' & asv_num_causal %in% c( "Bacterial Abundance" , "Synechococcus" ,         
                                                                                               "Chl-a"  ,      "Bacterial Production",    "PNF" ,              
                                                                                               "PNF 2-5 um" , "PNF 5 um" ,  "Cryptomonas"  ,       
                                                                                               "Micromonas", "HNF" , "HFN 2-5 um"   ,       
                                                                                               "HNF 5 um"  )) ~ 'bio',
                                            (fraction_causal == 'env' & 
                                               asv_num_causal %in% c( "Temperature" , "Day length",    "PO4"  ,  "NH4" , "NO2" , "NO3" , 
                                                                      "Si"  )) ~ 'phch',
                TRUE ~ fraction_causal))

data_ed2$asv_num_causal |>
  unique()

causal_effects_on_community_plot <- data_ed2 |>
  ggplot(aes(x = fraction_causal, y = rho))+
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
        strip.text = element_text(margin = margin(2, 5, 2, 5), size = 8),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        plot.margin = margin(l = 10, t = 5, b = 5, r = 10),
        legend.key.size = unit(3, 'mm'))+
  guides(alpha = 'none',
         color = guide_legend(ncol = 4))

causal_effects_on_community_plot
# 
# ggsave(plot = causal_effects_on_community_plot, filename = 'causal_effects_on_community_plot_v3.pdf',
#                 path = 'results/figures/',
#                 width = 180, height = 180, units = 'mm')


# ------ Figure ASV7 and ASV15 co-ocurrence and interaction ------ ########## ------
## ------ A panel -----
asv7_asv15_cluster <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv7', 'asv15')) |>
  dplyr::filter(abundance_type == 'rclr') |>
  ungroup() |>
  ggplot(aes(date, abundance_value))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  scale_y_continuous( expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 1,  position='stack')+
  scale_color_manual(values = palette_family_assigned_bloo)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(fraction), labeller = labs_fraction)+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(5,15,5,15)) 

asv7_asv15_cluster
# ggsave(filename = 'asv7_asv15_cluster_rclr_v1.pdf', plot = asv7_asv15_cluster,
#        path = 'results/figures/',
#        width = 188, height = 80, units = 'mm')

## ------ B panel -----
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

## ------ Panel construction A and B ----
asv7_asv15_cluster_interaction_plot <- plot_grid( asv7_asv15_cluster,
                                                  asv15_affectedby_asv7_plot,
                                                  labels = c('A', 'B'),
                                                  label_size = 10,
                                                  label_fontface = 'plain',
                                                  ncol = 1)

ggsave(
  # plot = asv7_asv15_cluster_interaction_plot, 
  # filename = 'asv7_asv15_cluster_interaction_plot.pdf',
  # path = 'results/figures/main_df3/',
  # width = 180, 
  # height = 150, 
  # units = 'mm'
)

###------ ASV7 is affected by ASV15 supplementary --------
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

# ggsave(
#   plot = asv7_affectedby_asv15_plot , 
#   filename = 'asv7_affectedby_asv15_plot.pdf',
#   path = 'results/figures/MDR/v2/',
#   width = 180, 
#   height = 100, 
#   units = 'mm'
# )

## Relationship between Bray Curtis Genes and Bray Curtis of the Community ----
#### first correlation: Bray-Curtis genes vs Bray-Curtis community, and Bray-Curtis genes vs. weighted UNIFRAC distance ----
corr_bray_02_plot <- bray_unifrac_eucl_tb_02 |>
  ggplot(aes(bray_curtis_community, genes))+
  scale_x_continuous(limits = c(0.05, 0.85))+
  scale_y_continuous(limits = c(0.05, 0.85))+
  #geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 1, alpha = 0.5)+
  geom_point(color = "#00808F", alpha = 0.8)+
  stat_cor(aes( #color = 'black', 
    
    label =   paste(..p.label..)), label.x = 0.35,
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#00808F"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.35, label.y = 0.2, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#00808F")+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(method = 'lm', color = "#00808F", fill = "#00808F")+
  
  labs(x = 'Bray-Curtis Community-Based', y = 'Bray-Curtis Genes-Based')+
  theme_bw()+
  theme(#axis.text.x = element_text(size = 6), 
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
    legend.position = 'bottom', axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10), strip.background = element_blank(), 
    legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
    #panel.border = element_blank(),
    strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_02_plot

corr_unifrac_02_plot <- bray_unifrac_eucl_tb_02 |>
  ggplot(aes(bray_curtis_community, genes))+
  #geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 1, alpha = 0.5)+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_point(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes), #shape = 2, 
             alpha = 0.8,  color = "#00808F")+
  scale_x_continuous(limits = c(0.05, 0.45))+
  scale_y_continuous(limits = c(0.15, 0.85))+
  stat_cor(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes, 
                                               label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.3,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#00808F"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes, 
                                               label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.25,  color = "#00808F" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes), color = "#00808F", fill ="#00808F")+
  labs(x = 'wUNIFRAC Community-Based', y = 'Bray-Curtis Genes-Based')+
  theme_bw()+
  theme(#axis.text.x = element_text(size = 6), 
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
    legend.position = 'bottom', axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10), strip.background = element_blank(), 
    legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
    panel.border = element_blank(),
    strip.placement = 'outside', aspect.ratio = 12/12)

corr_unifrac_02_plot

## I add KOs and genes in the same correlation plot. KOs will inform us if there is a change in the functionality -----
corr_bray_02_kos_plot <- corr_bray_genes_community_tb |>
  ggplot(aes(bray_curtis_community, KEGG_sgc))+
  scale_x_continuous(limits = c(0.05, 0.85))+
  scale_y_continuous(limits = c(0.05, 0.25))+
  #geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 1, alpha = 0.5)+
  geom_point( #shape = 2, 
    alpha = 0.8,  color = "#00808F")+
  stat_cor(aes( #color = 'black', 
    label =   paste(..p.label..)), label.x = 0.6,
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#00808F"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.6, label.y = 0.2, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#00808F")+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(method = 'lm', color = "#00808F", fill = "#00808F")+
  labs(x = 'Bray-Curtis Community-Based', y = 'Bray-Curtis KOs-Based')+
  theme_bw()+
  theme(#axis.text.x = element_text(size = 6), 
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
    legend.position = 'bottom', axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10), strip.background = element_blank(), 
    legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
    #panel.border = element_blank(),
    strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_02_kos_plot

corr_unifrac_02_ko_plot <- corr_bray_genes_community_tb |>
  ggplot(aes(wunifrac_distance,KEGG_sgc))+
  #geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 1, alpha = 0.5)+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_point(data = corr_bray_genes_community_tb, aes(wunifrac_distance ,KEGG), #shape = 2, 
             alpha = 0.8,  color = "#00808F")+
  scale_x_continuous(limits = c(0.05, 0.45))+
  scale_y_continuous(limits = c(0.05, 0.25))+
  stat_cor(data = corr_bray_genes_community_tb, aes(wunifrac_distance, KEGG, 
                                                    label =   paste(..p.label..)), 
           label.x = 0.3,  
           label.y = 0.25,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#00808F"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = corr_bray_genes_community_tb, aes(wunifrac_distance, KEGG, 
                                                    label =   paste(..r.label..)),
           label.x = 0.3, label.y = 0.2,  color = "#00808F" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  geom_smooth(method = 'lm', data = corr_bray_genes_community_tb, aes(wunifrac_distance, KEGG), color = "#00808F", fill ="#00808F")+
  labs(x = 'wUNIFRAC Community-Based', y = 'Bray-Curtis KOs-Based')+
  theme_bw()+
  theme(#axis.text.x = element_text(size = 6), 
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
    legend.position = 'bottom', axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10), strip.background = element_blank(), 
    legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
    panel.border = element_blank(),
    strip.placement = 'outside', aspect.ratio = 12/12)

corr_unifrac_02_ko_plot

# Now arrange the full layout, with bray_unifrac_eucl_plot occupying the top row
corr_bray_unifrac_kos_plot <- plot_grid(
  corr_bray_02_plot, #corr_unifrac_02_plot,
  corr_bray_02_kos_plot, #corr_unifrac_02_ko_plot,
  # Top plot (spanning both columns)          # Second row with two plots
  ncol = 1,                # One column layout for the main grid
  rel_heights = c(1, 1),   # First plot 3 times the height of the second row
  labels = c('A', 'B'), #'C', 'D',
  label_fontface = 'plain'
)

# Print the final plot
print(corr_bray_unifrac_kos_plot)

# ggsave( plot = corr_bray_unifrac_kos_plot,
#         filename = 'corr_bray_unifrac_plot_v7.pdf',
#         path = 'results/figures/',
#         width = 88, height = 120, units = 'mm')




