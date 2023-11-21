
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

##highligh harbour restoration period ----
### The remodelation of the Blanes harbour strated on 24th March 2010 and finished on the 9th of june 2012
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

# palettes----
palette_seasons_4 <- c("winter" = "#002562", 'spring' = "#519741", 'summer' = "#ffb900",'autumn' =  "#96220a")

palette_fraction <- c('0.2' = '#00808F', '3' = '#454545')

## palettes taxonomy assigned----
### they need to be updated to the new taxonomy

# asv_tab_all_bloo_z_tax %$%
#   domain |>
#   unique()

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
asv_tab_all_bloo_z_tax %$%
  phylum |>
  unique()

asv_tab_all_bloo_z_tax |> 
  dplyr::filter(is.na(phylum))

asv_tab_all_bloo_z_tax |> 
  dplyr::filter(is.na(domain)) |>
  distinct(asv_num)
palette_phylums_assigned_bloo <- c('Proteobacteria' = "#fcca46", 'Cyanobacteria' = "#009e73",
                                   "Bacteroidota" = "#0051BF",  'Verrucomicrobiota' = '#005c69', 'Planctomycetota' = "#69267e",
                                   'Bdellovibrionota' = "#8c789d") #, NA == "#000000"

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
# asv_tab_all_bloo_z_tax %$%
#   class |>
#   unique()

palette_class_assigned_bloo <- c('Gammaproteobacteria' = '#FFA737', 'Alphaproteobacteria' = '#B0413E', 
                                 'Cyanobacteriia'  =  '#009F6A', 'Bacteroidia' = '#0051BF',  
                                'Verrucomicrobiae' = '#005c69', 'Phycisphaerae' =  '#e3a6ce',   
                                'Bdellovibrionia' = '#8C789D') #, NA == "#000000"

# bloomers_tax |>
#   dplyr::filter(taxonomy_rank == 'sum_f') |>
#   left_join(tax_bbmo_10y, by = c('taxonomy' = 'family')) |>
#   distinct(taxonomy, order)

# asv_tab_all_bloo_z_tax$order |>
#   unique()

 # asv_tab_all_bloo_z_tax |>
 #  dplyr::filter(order == 'Synechococcales') |>
 #  dplyr::select(class) |>
 #  distinct()

asv_tab_all_bloo_z_tax |>
  dplyr::select(class_f, order_f) |>
  dplyr::filter(order_f %in% c('Oceanospirillales', 'Opitutales')) |>
  distinct()

palette_order_assigned_bloo <-  c('Thiotrichales' = '#BB4430', "Alteromonadales" =  '#A63B00',  
                                  "Vibrionales" = '#F2AC5D', "Enterobacterales" = '#FFA200',
                                  "Cellvibrionales"   = '#F35900',  "Pseudomonadales"  = '#FF8E00', 
                                  #"SAR86 clade" = '#FFBF45',
                                  "SAR11 clade" =  '#B0413E',  "Rhodobacterales" = '#C55E5C',
                                  "Sphingomonadales"  = '#8C000A', "Puniceispirillales" = '#931F1D',  'Rhizobiales' = '#B31722',
                                  "Rhodospirillales"  = '#FFA197', 
                                  "Synechococcales"  = '#009F6A', 
                                  "Flavobacteriales"   =  '#0051BF',  'Chitinophagales' = '#92ABFF', 
                                  "Verrucomicrobiales"= '#005c69', "Opitutales"   =   '#74B9C8',
                                  "Phycisphaerales"  = '#e3a6ce',  
                                  "Bacteriovoracales" = '#8C789D'
                                  #"Oceanospirillales" =  '#A05C00'
) # NA ==  "#000000",

# asv_tab_all_bloo_z_tax$family_f |>
# unique()

# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(family == "Puniceicoccaceae") |>
#   dplyr::select(order) |>
#   distinct()

palette_family_assigned_bloo <- c("Thiotrichaceae" = '#BB4430',  
                                  "Marinobacteraceae" =  '#DE6931',  "Alcanivoracaceae1"  =  '#A05C00',  "Moraxellaceae"  =  '#FF8E00', #pseudomonadales
                                  "Halieaceae" = '#F35900', "SAR86 clade" = '#FFBA00', #pseudomonadales
                                  "Vibrionaceae"      = '#F2AC5D',   "Yersiniaceae"  = '#FFA200',    #enterobacterales    
                                  "Alteromonadaceae"   =  '#A63B00',     
                                  "Clade II" = '#B0413E',  "Clade I" = '#CD7F78', "Rhodobacteraceae" = '#C55E5C',
                                  "Sphingomonadaceae"   = '#8C000A',  "SAR116 clade"  ='#931F1D', "Stappiaceae" = '#B31722',
                                  "AEGEAN-169 marine group" =  '#690000', 
                                  "Cyanobiaceae"  = '#009F6A', 
                                  "NS7 marine group" =  '#92ABFF',
                                  "NS9 marine group"  =  '#3B52A3',   "Cryomorphaceae" = '#002A8E',
                                  "Saprospiraceae" = '#5F7CCB',
                                  "Flavobacteriaceae"   =  '#0051BF',  
                                  "Rubritaleaceae"  =  '#005c69',       
                                  "DEV007" = '#74B9C8',  "Puniceicoccaceae"   = '#29335C',     
                                  "Phycisphaeraceae"   = '#e3a6ce',    
                                  "Bacteriovoracaceae" =  '#8C789D' 
                                 )  # NA == "#000000" 
##add genus color 
asv_tab_all_bloo_z_tax |>
  dplyr::filter(genus == "Roseibacillus") |>
  dplyr::select(family) |>
  distinct()

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

palette_all_tax_rank <- c('Proteobacteria' = "#fcca46", 'Cyanobacteria' = "#009e73",
                          "Bacteroidota" = "#0051BF",  'Verrucomicrobiota' = '#005c69', 'Planctomycetota' = "#69267e",
                          'Bdellovibrionota' = "#8c789d",
                          ##class
                          'Gammaproteobacteria' = '#FFA737', 'Alphaproteobacteria' = '#B0413E', 
                          'Cyanobacteriia'  =  '#009F6A', 'Bacteroidia' = '#0051BF',  
                          'Verrucomicrobiae' = '#005c69', 'Phycisphaerae' =  '#e3a6ce',   
                          'Bdellovibrionia' = '#8C789D',
                          ##order
                          'Thiotrichales' = '#BB4430', "Alteromonadales" =  '#A63B00',  
                          "Vibrionales" = '#F2AC5D', "Enterobacterales" = '#FFA200',
                          "Cellvibrionales"   = '#F35900',  "Pseudomonadales"  = '#FF8E00', 
                          #"SAR86 clade" = '#FFBF45',
                          "SAR11 clade" =  '#B0413E',  "Rhodobacterales" = '#C55E5C',
                          "Sphingomonadales"  = '#8C000A', "Puniceispirillales" = '#931F1D',  'Rhizobiales' = '#B31722',
                          "Rhodospirillales"  = '#FFA197', 
                          "Synechococcales"  = '#009F6A', 
                          "Flavobacteriales"   =  '#0051BF',  'Chitinophagales' = '#92ABFF', 
                          "Verrucomicrobiales"= '#005c69', "Opitutales"   =   '#74B9C8',
                          "Phycisphaerales"  = '#e3a6ce',  
                          "Bacteriovoracales" = '#8C789D',
                          ##family
                          "Thiotrichaceae" = '#BB4430',  
                          "Marinobacteraceae" =  '#DE6931',  "Alcanivoracaceae1"  =  '#A05C00',  "Moraxellaceae"  =  '#FF8E00', #pseudomonadales
                          "Halieaceae" = '#F35900', "SAR86 clade" = '#FFBA00', #pseudomonadales
                          "Vibrionaceae"      = '#F2AC5D',   "Yersiniaceae"  = '#FFA200',    #enterobacterales    
                          "Alteromonadaceae"   =  '#A63B00',     
                          "Clade II" = '#B0413E',  "Clade I" = '#CD7F78', "Rhodobacteraceae" = '#C55E5C',
                          "Sphingomonadaceae"   = '#8C000A',  "SAR116 clade"  ='#931F1D', "Stappiaceae" = '#B31722',
                          "AEGEAN-169 marine group" =  '#690000', 
                          "Cyanobiaceae"  = '#009F6A', 
                          "NS7 marine group" =  '#92ABFF',
                          "NS9 marine group"  =  '#3B52A3',   "Cryomorphaceae" = '#002A8E',
                          "Saprospiraceae" = '#5F7CCB',
                          "Flavobacteriaceae"   =  '#0051BF',  
                          "Rubritaleaceae"  =  '#005c69',       
                          "DEV007" = '#74B9C8',  "Puniceicoccaceae"   = '#29335C',     
                          "Phycisphaeraceae"   = '#e3a6ce',    
                          "Bacteriovoracaceae" =  '#8C789D',
                          ##genus
                          'unclassified' = '#534F4A', 
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
                          "Erythrobacter"  = '#8C000A',  "Sphingobium" = '#70161E',          
                          "Glaciecola" =  '#A63B00', 
                          "Vibrio"    = '#F2AC5D',   "Serratia"  = '#FFA200',  
                          "Alcanivorax"    =  '#A05C00',                    
                          "OM60(NOR5) clade"  = '#F35900', 
                          "Marinobacter" =  '#DE6931', 
                          "Acinetobacter" = '#FF8E00', "Psychrobacter"  = '#E55934', 
                          "Lentimonas"  = '#fcca46',   
                          "Roseibacillus" =  '#005c69'
                          
)     
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

setwd("~/Documentos/Doctorat/BBMO/BBMO_bloomers/")
bbmo_10y <-readRDS("data/blphy10years.rds") ##8052 asv amb totes les mostres, no está en % aquesta taula
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

tax_bbmo_10y_old <- bbmo_10y@tax_table |> 
  #mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(bbmo_10y@otu_table))) |> #ja té uns nº d'OTUs
  as_tibble()

m_bbmo_10y <- bbmo_10y@sam_data |>
  as_tibble()

m_bbmo_10y |>
  colnames()

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
                          "year", "month", "day", "season",            
                          "bacteria_joint", "synechococcus", "depth", "name_complete")

##updated tax with SILVA 138 and DECIPHER (it has many NAs, even at Phylum taxonomic level)
# new_tax <-  readRDS('data/03_tax_assignation/devotes_all_tax_assignation.rds') |>
#   as_tibble(rownames = 'sequence')

#new taxonomy created with the database SILVA 138.1 using Assign tax at 50 (default)
new_tax <-  readRDS('data/03_tax_assignation/devotes_all_assign_tax_assignation_v2.rds') |>
  as_tibble(rownames = 'sequence')

tax_bbmo_10y_new <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax, by = c('seq' = 'sequence'))

##create a fasta file with ASV sequences to update the taxonomy----
# tax_bbmo_10y |>
#   mutate(fasta_format = paste0(">",asv_num,"\n",as.character(seq))) |>
#   pull(fasta_format) |>
#   cat(file = 'asv_seqs_bbmo10y.fasta', sep = "\n")


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

ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = class_f), scale = 1, alpha = 0.7,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)

###small plot of the distribution of the relative abundances over all the dataset----
# library(ggridges)
# asv_tab_10y_l_rel |>
#   dplyr::mutate(project = 'BBMO10Y') |>
#   left_join(tax_bbmo_10y_new_assign, by = 'asv_num') |>
#   dplyr::filter(relative_abundance > 0.00) |>
#   ggplot(aes( x = relative_abundance, fill = project))+
#   geom_density(aes(group = project, fill = project))+
#  # geom_histogram(bins = 80)+
#   scale_y_continuous(expand = c(0,0))+
#   scale_x_continuous(labels = percent_format(), limits = c(0, 0.15))+ #, limits = c(0, 0.02)
#   #theme_ridges()+
#   labs(y = 'Density', x = 'Relative abundance (%)')+
#   theme_bw()+
#   theme(legend.position = 'none')

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

##with this transformation I'm losing samples (due to too much 0 in some samples, z.warning set up to 0.99 to keep all samples)
### at 0.8 (default) I lose 30 samples which belonged to the years corresponding to harbour remodelation 
### I don't lose samples but I lose ASVs.
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

##dimensions are different from relative and pseudo dataset and zclr understand why----
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

nrow(asv_tab_10y_02_pseudo_zclr) == nrow(asv_tab_10y_02_pseudo)
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

#this part is commented since is computationally slow so i don't want to repeat this step many times. 
##this function gives us a randomly rarefied community data
# rrarefy(asv_tab_bbmo_10y_w, sample = min_n_seqs) |>
#   as_tibble(rownames = 'sample_id') #just rarefying (one simgle subsampling)

## we use rarefaction (which repeats the subsampling step many times)
## perform this a 1000 times to get an empirical diversity values calculating the mean value for each timepoint.
# asv_tab_bbmo_10y_w_rar <- rrarefy.perm(asv_tab_bbmo_10y_w, 
#                                        sample = min_n_seqs, 
#                                        n = 1000, 
#                                        round.out = T)
#write.csv2(asv_tab_bbmo_10y_w_rar, file = 'asv_tab_bbmo_10y_w_rar.csv')

asv_tab_bbmo_10y_w_rar <- read.csv2('data/asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()|> 
  rename('sample_id' = 'X') 

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
source('../../../Bloomers/R/community_evenness.R')

community_eveness_02 <- asv_tab_bbmo_10y_w_rar |>
  #as.data.frame() |>
  #tibble::rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou'))

community_eveness_3 <- asv_tab_bbmo_10y_w_rar |>
  #as.data.frame() |>
  #rownames_to_column(var = 'sample_id') |>
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
source('../../../Bloomers/R/compute_bray_curtis_dissimilariy.R') #we need to upload it since we updated the function but I didn't compile the package
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
  # as.data.frame() |>
  # rownames_to_column(var = 'sample_id') |>
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
## For each ASVs based on relative abundances and pseudoabundances -----
### those ASVs that are present > 50% of the sampless
asv_tab_10y_02_pseudo %$%
  sample_id |>
  unique() |>
  summary() #60 is half of the dataset

asv_tab_10y_02_zclr |>
  colnames()

### I filter by relative abundance > 10% at some point so that the computer can compute this easily

x <- 120*0.75 ##percentage of ASV present in the dataset that we want to subset by (occurrence)

z_02 <- asv_tab_10y_02_pseudo |>
  inner_join(asv_tab_10y_02_zclr, by = c('sample_id', 'asv_num')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  #group_by(asv_num) |>
  dplyr::reframe(anomalies_ps = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
                 anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = zclr, plotting = FALSE)[c(1,2,3)])

asv_tab_10y_3_pseudo %$%
  sample_id |>
  unique() |>
  summary() 

##check if I'm finding the correct number of potential bloomers ASVs----
 n_bloomers_02 <-  asv_tab_10y_02_pseudo |>
  inner_join(asv_tab_10y_02_zclr, by = c('sample_id', 'asv_num')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(zclr >= 1.96) |>
  dplyr::distinct(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(sum = sum(n))

bloo_02 <- asv_tab_10y_02_pseudo |>
  inner_join(asv_tab_10y_02_zclr, by = c('sample_id', 'asv_num')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(zclr >= 1.96) |>
  dplyr::distinct(asv_num) |>
  as_vector()

n_bloomers_3 <-  asv_tab_10y_3_pseudo |>
  inner_join(asv_tab_10y_3_zclr, by = c('sample_id', 'asv_num')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(zclr >= 1.96) |>
  dplyr::distinct(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(sum = sum(n))

bloo_3 <- asv_tab_10y_3_pseudo |>
  inner_join(asv_tab_10y_3_zclr, by = c('sample_id', 'asv_num')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(zclr >= 1.96) |>
  dplyr::distinct(asv_num) |>
  as_vector()

x <- 117*0.75  ##percentage of ASV present in the dataset that we want to subset by (occurrence)
z_3 <- asv_tab_10y_3_pseudo |>
  inner_join(asv_tab_10y_3_zclr, by = c('sample_id', 'asv_num')) |>
  group_by(asv_num) |>
  # dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  # dplyr::filter(num_0 <= x) |>  #filter those ASVs that are 0 more than 50% of the dataset
  #as_tibble() |>
  #group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  dplyr::reframe(anomalies_ps = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
    anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
    anomalies_clr = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = zclr, plotting = FALSE)[c(1,2,3)])

asv_tab_10y_3_pseudo |>
  inner_join(asv_tab_10y_3_zclr, by = c('sample_id', 'asv_num')) |> #asv_tab_10y_02_zclr vull afegir el zclr per calcular també les seves anomalies i veure si veiem el mateix
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.1)) |>
  dplyr::filter(zclr >= 1.96) |>
  dplyr::distinct(asv_num) 

## At the level of community, we use the Evenness result and Bray Curtis dissimilarity ----
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
                                    logic2 = NULL,
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

## I only filter for those anomalies in relative abundance because pseudoabundance I can only use it for fl not for PA and zclr has a problem with 
## dealing with many 0.

asv_anom_02 <- find_asv_with_anomalies(anomalies_result = z_02, anomaly_in1 = anomalies_ra, 
                                       anomaly_in2 = anomalies_ps,
                                       anomaly_in3 = anomalies_clr, 
                                    logic1 = 'TRUE', logic2 = NULL, 
                                    logic3 = NULL,
                                    asv_col = asv_num)
##me'n surten 19 ara mateix perque?
##86 ASVs cumplint les 3 condicions i sense CLR també

##no acabo d'entendre perquè torna un vector amb numero consecutiu d'ASVs...

asv_anom_3 <- find_asv_with_anomalies(anomalies_result = z_3, anomaly_in1 = anomalies_ra, 
                                      anomaly_in2 = anomalies_clr, 
                                      anomaly_in3 = NULL, #anomalies_ps
                                      logic1 = 'TRUE', logic2 = NULL, 
                                      logic3 = NULL, ##per 3 no és representatiu la pseudoabund
                                      asv_col = asv_num)

##45 ASVS
## número d'ASVs comuns i unics a cada fracció
asv_anom_3_tb <- asv_anom_3 |>
  as_tibble()

asv_anom_02_tb <- asv_anom_02 |>
  as_tibble()

common_bloomers_tax <- asv_anom_3_tb |>
  bind_rows(asv_anom_02_tb) |>
  unique() |> ##96 totals >50% del dataset
  left_join(tax_bbmo_10y_new_assign, by = c('value' = 'asv_num'))

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
  colnames()

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

z_scores_all <- z_scores_02 |>
  bind_rows(z_scores_3)
  
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

# Dataset with all information ----
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo |>
  left_join(z_scores_all) |> 
  left_join(tax_bbmo_10y_new, by = 'asv_num') 

asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  rename(domain = Kingdom, phylum = Phylum, class = Class, order = Order, family = Family, genus = Genus)

#write.csv2(asv_tab_all_bloo_z_tax, 'data/asv_tab_all_bloo_z_tax_new_assign.csv')
asv_tab_all_bloo_z_tax |>
  distinct(asv_num) |>
  dim()

asv_tab_all_bloo_z_tax |>
  colnames()


##UPLOAD BLOOMERS DATA-----
asv_tab_all_bloo_z_tax_old <- read.csv2('data/asv_tab_all_bloo_z_tax.csv') ## using silva 132 database
asv_tab_all_bloo_z_tax <- read.csv2('data/
                                    ') ## using decipher with silva 138 
asv_tab_all_bloo_z_tax <- read.csv2('data/asv_tab_all_bloo_z_tax_new_assign.csv') ##using dada2 classifier assign tax with silva 138.1

##upload diversity data----
## Rarefied dataset to calculate Community Evenness----
source('~/Documentos/Doctorat/Bloomers/R/community_evenness.R')

community_eveness_02 <- asv_tab_bbmo_10y_w_rar |>
  #as.data.frame() |>
  #tibble::rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou'))

community_eveness_3 <- asv_tab_bbmo_10y_w_rar |>
  #as.data.frame() |>
  #rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'reads_rar') |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  #dplyr::select(sample_id, reads, asv_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou'))

### We need the rarefied table transformed to relative abundances
asv_tab_bbmo_10y_l_rel_rar <- asv_tab_bbmo_10y_w_rar |>
  # as.data.frame() |>
  # rownames_to_column(var = 'sample_id') |>
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

### plot Bray-Curtis dissimilarity and Community Eveness together----
community_eveness_all <- community_eveness_02 |>
  bind_rows(community_eveness_3)

bray_curtis_rar_all <- bray_curtis_02_rar |> ##one sample less, the first one can't be compared with the previous
  bind_rows(bray_curtis_3_rar)

## At the level of community, we use the Evenness result and Bray Curtis dissimilarity ----
z_diversity <- bray_curtis_02_rar |>
  dplyr::right_join(community_eveness_02, by = join_by("samples" == "sample_id")) |> 
  #ungroup() %>%
  #group_by(sample_id) %>%
  dplyr::reframe(anomalies_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = TRUE)[c(1,2,3)],# ),
                 anomalies_eveness = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = TRUE)[c(1,2,3)])

z_diversity %>%
  str()

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

##reorder taxonomy as factors ----
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

asv_tab_all_bloo_z_tax$season <- asv_tab_all_bloo_z_tax$season |>
  factor(levels = c('winter', 'spring', 'summer', 'autumn'))


# BLOOMERS TAXONOMY EXPLORATION ----
#### We create a hierarchical piechart ----
#### Separate between those that bloom in FL from those that bloom in PA

# bloomers_tax <- asv_tab_all_bloo_z_tax |>
#   dplyr::select(asv_num_f, phylum_f, class_f, order_f, family_f, genus) |>
#   distinct() |>
#   group_by(phylum_f, class_f, order_f, family_f, genus) |>
#   dplyr::mutate(genus = case_when(genus = is.na(genus) ~ 'unclassified',
#                                   genus != is.na(genus) ~ genus)) |>
#   dplyr::summarize(n_bloom = n()) |>
#   ungroup() |>
#   dplyr::mutate(bloom_perc_genus = n_bloom/sum(n_bloom)) |>
#   group_by(phylum_f, class_f, order_f, family_f, genus) |>
#   mutate(sum_f = sum(bloom_perc_genus)) |>
#   
#   group_by(phylum_f, class_f, order_f, family_f) |>
#   dplyr::mutate(n_bloom = n()) |>
#   ungroup() |>
#   dplyr::mutate(bloom_perc_family = n_bloom/sum(n_bloom)) |>
#   group_by(phylum_f, class_f, order_f, family_f) |>
#   mutate(sum_f = sum(bloom_perc_family)) |>
#   
#   group_by(phylum_f, class_f, order_f) |>
#   dplyr::mutate(n_bloom = n()) |>
#   ungroup() |>
#   dplyr::mutate(bloom_perc_order = n_bloom/sum(n_bloom)) |>
#   group_by(phylum_f, class_f, order_f) |>
#   mutate(sum_o = sum(bloom_perc_order)) |>
#   
#   group_by(phylum_f, class_f) |>
#   dplyr::mutate(n_bloom = n()) |>
#   ungroup() |>
#   dplyr::mutate(bloom_perc_class = n_bloom/sum(n_bloom)) |>
#   group_by(phylum_f, class_f) |>
#   mutate(sum_c = sum(bloom_perc_class)) |>
#   
#   group_by(phylum_f) |>
#   dplyr::mutate(n_bloom = n()) |>
#   ungroup() |>
#   dplyr::mutate(bloom_perc_phylum = n_bloom/sum(n_bloom)) |>
#   group_by(phylum_f) |>
#   mutate(sum_p = sum(bloom_perc_phylum)) |>
#   pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage')
# 
# bloomers_tax$taxonomy_rank <- bloomers_tax$taxonomy_rank |> 
#   factor(levels = c('bloom_perc_phylum', 'bloom_perc_class', 'bloom_perc_order', 'bloom_perc_family', 'bloom_perc_genus'))

### we create a function that will allow us to calculate the % of each bloomer for the different taxonomic ranks
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

# bloomers_tax$taxonomy_rank <- bloomers_tax$taxonomy_rank |> 
#   factor(levels = c('sum_p', 'sum_c', 'sum_o', 'sum_f', 'sum_g'))

bloomers_tax_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% bloo_02) |>
  bloomers_tax_rank() |>
  dplyr::mutate(community = 'Free living (0.2-3 um)',
                n_bloo = '19')

bloomers_tax_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% bloo_3) |>
  bloomers_tax_rank() |>
  dplyr::mutate(community = 'Particle attached (3-20 um)',
                n_bloo = '45') #number of ASVs that bloom in PA

# asv_tab_all_bloo_z_tax
# 
# data_ed <- asv_tab_all_bloo_z_tax |>
#   dplyr::select(asv_num_f, phylum_f, class_f, order_f, family_f, genus) |>
#   distinct() |>
#   dplyr::mutate(genus = case_when(genus = is.na(genus) ~ 'unclassified',
#                                   genus != is.na(genus) ~ genus))

# bloom_g <- data_ed  |>
#   group_by(phylum_f, class_f, order_f, family_f, genus) |>
#   dplyr::summarize(n_bloom = n()) |>
#   ungroup() |>
#   dplyr::mutate(bloom_perc_genus = n_bloom/sum(n_bloom)) |>
#   group_by(phylum_f, class_f, order_f, family_f, genus) |>
#   mutate(sum_g = sum(bloom_perc_genus)) |>
#   pivot_longer(cols = starts_with('sum'), names_to = 'taxonomy_rank', values_to = 'percentage') |>
#   dplyr::select(genus, taxonomy_rank, percentage) |>
#   ungroup() |>
#   distinct(genus, family_f, taxonomy_rank, percentage) |>
#   rename( taxonomy = genus) |>
#   dplyr::select(-family_f)
# 
# bloom_g |>
#   group_by(taxonomy_rank) |>
#   summarize(n = sum(percentage))
# 
# bloom_c |>
#   bind_rows(bloom_p)
# 
# tax_bbmo_10y |>
#   head()
# 
# tax_bbmo_10y |>
#   dplyr::select(-c(asv_num, otu_corr, seq)) |>
#   pivot_longer(cols = )
# 
# bloomers_tax |>
#   dplyr::filter(taxonomy_rank == 'sum_p')

bloomers_tax_02$taxonomy_rank <- factor(bloomers_tax_02$taxonomy_rank, 
                                     levels = c('sum_p', 'sum_c', 'sum_o', 'sum_f', 'sum_g'))

bloomers_tax_3$taxonomy_rank <- factor(bloomers_tax_3$taxonomy_rank, 
                                        levels = c('sum_p', 'sum_c', 'sum_o', 'sum_f', 'sum_g')) 

bloomers_tax <- bloomers_tax_02 |>
  bind_rows(bloomers_tax_3)

bloomers_tax$taxonomy_rank <- factor(bloomers_tax$taxonomy_rank, 
                                       levels = c('sum_p', 'sum_c', 'sum_o', 'sum_f', 'sum_g')) 

labs_taxonomic_rank <- as_labeller(c("sum_p" = 'Phylum' ,
                          "sum_c" = 'Class',
                          "sum_o"  = 'Order',    
                          "sum_f" = 'Family',
                          "sum_g" = 'Genus'))

bloomers_tax$taxonomy <-  factor(bloomers_tax$taxonomy, levels = c('Proteobacteria', 'Cyanobacteria',
                                                                   "Bacteroidota",  'Verrucomicrobiota' , 'Planctomycetota' ,
                                                                   'Bdellovibrionota',
                                                                   ##class
                                                                   'Gammaproteobacteria', 'Alphaproteobacteria' , 
                                                                   'Cyanobacteriia' , 'Bacteroidia',  
                                                                   'Verrucomicrobiae' , 'Phycisphaerae',   
                                                                   'Bdellovibrionia' ,
                                                                   ##order
                                                                   'Thiotrichales', "Alteromonadales" ,  
                                                                   "Vibrionales", "Enterobacterales" ,
                                                                   "Cellvibrionales",  "Pseudomonadales", 
                                                                   #"SAR86 clade" = '#FFBF45',
                                                                   "SAR11 clade",  "Rhodobacterales" ,
                                                                   "Sphingomonadales", "Puniceispirillales",  'Rhizobiales',
                                                                   "Rhodospirillales" , 
                                                                   "Synechococcales", 
                                                                   "Flavobacteriales" ,  'Chitinophagales', 
                                                                   "Verrucomicrobiales", "Opitutales",
                                                                   "Phycisphaerales" ,  
                                                                   "Bacteriovoracales",
                                                                   ##family
                                                                   "Thiotrichaceae",  
                                                                   "Marinobacteraceae",  "Alcanivoracaceae1",  "Moraxellaceae", #pseudomonadales
                                                                   "Halieaceae" , "SAR86 clade", #pseudomonadales
                                                                   "Vibrionaceae" ,   "Yersiniaceae" ,    #enterobacterales    
                                                                   "Alteromonadaceae" ,     
                                                                   "Clade II" ,  "Clade I", "Rhodobacteraceae",
                                                                   "Sphingomonadaceae",  "SAR116 clade"  , "Stappiaceae",
                                                                   "AEGEAN-169 marine group", 
                                                                   "Cyanobiaceae" , 
                                                                   "NS7 marine group" ,
                                                                   "NS9 marine group" ,   "Cryomorphaceae",
                                                                   "Saprospiraceae",
                                                                   "Flavobacteriaceae",  
                                                                   "Rubritaleaceae",       
                                                                   "DEV007",  "Puniceicoccaceae",     
                                                                   "Phycisphaeraceae" ,    
                                                                   "Bacteriovoracaceae",
                                                                   ##genus
                                                                   'unclassified', 
                                                                   "Marixanthomonas" ,  "NS4 marine group",                
                                                                   "NS5 marine group" , "Peredibacter"  ,
                                                                   "Prochlorococcus MIT9313" ,       
                                                                   "Synechococcus CC9902",
                                                                   "Urania-1B-19 marine sediment group",
                                                                   "Candidatus Puniceispirillum" ,    
                                                                   "Amylibacter" , "Ascidiaceihabitans" ,
                                                                   "HIMB11",                      
                                                                   "Limimaricola" ,
                                                                   "Clade Ia" ,  "Clade Ib" ,                     
                                                                   "Erythrobacter",  "Sphingobium" ,          
                                                                   "Glaciecola" , 
                                                                   "Vibrio" ,   "Serratia",  
                                                                   "Alcanivorax"  ,                    
                                                                   "OM60(NOR5) clade", 
                                                                   "Marinobacter", 
                                                                   "Acinetobacter", "Psychrobacter", 
                                                                   "Lentimonas",   
                                                                   "Roseibacillus"))


bloomers_tax |>
  #dplyr::filter(taxonomy_rank == 'bloom_perc_phylum') |>
  # group_by(phylum) |>
  # dplyr::mutate(sum = if_else(taxonomy_rank == 'bloom_perc_phylum',  sum(percentage), '0', missing = NULL)) |>
  #distinct() |>
 # dplyr::select(-n_bloom) |>
  #dplyr::mutate(tax_combined = paste(phylum, class, order, family, genus)) |> 
  #distinct(tax_combined,  percentage, taxonomy_rank) |>
  ggplot(aes(y =taxonomy_rank, x = percentage, fill = taxonomy))+
 # scale_fill_manual(values = paired_12_12)+ #values = palette_phylums_assigned
  labs(x='% of potential bloomers', y = 'Taxonomy rank', fill = 'Taxonomy')+
  scale_fill_manual(values = palette_all_tax_rank)+
  scale_y_discrete(labels = labs_taxonomic_rank)+
  coord_polar()+
  geom_col()+
  facet_wrap(vars(community))+
  geom_text(mapping = aes(x = 0.9, y = 5.55, label = paste0('n = ', n_bloo)),
            check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 0,
            color = 'black', size = 2.5)+
  guides(fill =guide_legend(ncol = 3, size = 4))+
  # annotate('text', x = 0, y = c('sum_c' ,'sum_f', 'sum_g', 'sum_o', 'sum_p'),
  #          label = c('Class', 'Family', 'Order', 'Genus', 'Phylum'))+
  #scale_y_continuous(expand=c(0, 0)) +
  theme_bw()+
  theme(legend.position = 'right', panel.grid = element_blank(), panel.border = element_blank(),
        strip.background = element_blank(), text = element_text(size = 5), legend.key.size = unit(4, 'mm'))

## Plot Evenness and Bray Curtis anomalies----
### crec que no té sentit perquè ja ens indica un canvi en la comunitat al menys la Bray-Curtis dissimilarity
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


library(ggforce)
library(scales)
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
  scale_color_manual(values = palette_class_assigned_bloo)+
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

### Plot blooming events with geom_area during the whole timeseries----
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 

bray_curtis_rar_all_m <- bray_curtis_rar_all |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result))

community_eveness_all_m <- community_eveness_all |>
  #left_join(m_bbmo_10y, by = c('sample_id' = 'sample_id')) |>
  dplyr::filter(community_eveness_rar != is.na(community_eveness_rar)) |>
  left_join(z_scores_ev) |>
  dplyr::mutate(z_scores_ev = case_when(is.na(z_score_ev) ~ 0,
                                        z_score_ev == 'NaN' ~ 0,
                                        z_score_ev == Inf ~ 10000,
                                        TRUE ~ z_score_ev)) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_ev) >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

community_eveness_all_m %$%
  max(community_eveness_rar)

community_eveness_all_m |>
  colnames()

community_eveness_all_m %$%
  anomaly_color |>
  unique()

asv_tab_all_bloo_z_tax|>
  dplyr::filter(abundance_type == 'relative_abundance') %$%
  max(abundance_value)

community_eveness_all_m %$%
  min(date)

asv_tab_all_bloo_z_tax  |>
  dplyr::filter(abundance_type == 'relative_abundance') %$%
  class_f |>
  unique()

# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(class == is.na(class))

 #booming_events_ev_BBMO10Y <-
  # asv_tab_all_bloo_z_tax |>
  # dplyr::filter(abundance_type == 'relative_abundance') |>
  #   # group_by(date) |>
  #   # dplyr::mutate(max_abund = sum(abundance_value)) |>
  #   # ungroup() |>
  # #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  # ggplot(aes(date, abundance_value))+
  #   #geom_line(aes(date, max_abund))+
  # #geom_segment(aes(x = '2005-01-01', y = 0, xend = '2005-01-02', yend =0.57),color="black")+
  # scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
  #                  #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
  #                             #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  #                  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  # 
  # #geom_stream(aes(fill = class_f, group = class_f), type = "ridge", bw=1)+
  # geom_area(aes(fill = class_f, group = class_f), alpha = 0.8,  position='stack')+
  # #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),  
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  # scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.57), 
  #                    sec.axis = sec_axis(~.* 1.6 , name = 'Community Evenness'))+
  # scale_color_identity()+
  # scale_fill_manual(values = palette_class_assigned_bloo, na.value = "#000000")+
  # labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Family')+
  # facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  #   #facet_wrap(fraction~phylum_f, dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  # guides(fill = guide_legend(ncol = 6, size = 10,
  #                             override.aes = aes(label = '')),
  #        alpha = 'none')+
  # theme_bw()+
  # theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
  #       panel.grid.major = element_blank(), strip.text = element_text(size = 7),
  #       legend.position = 'bottom', axis.text.y = element_text(size = 8),
  #       axis.title = element_text(size = 8), strip.background = element_blank(), 
  #       legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')

# ggsave('booming_events_ev_BBMO10Y_new_tax.pdf', booming_events_ev_BBMO10Y,
#        path = "../results/figures/",
#        width = 180,
#        height = 160,
#        units = 'mm')
#### plot a non overlaping area ----
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 

bray_curtis_rar_all_m <- bray_curtis_rar_all |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result))

community_eveness_all_m <- community_eveness_all |>
  #left_join(m_bbmo_10y, by = c('sample_id' = 'sample_id')) |>
  dplyr::filter(community_eveness_rar != is.na(community_eveness_rar)) |>
  left_join(z_scores_ev) |>
  dplyr::mutate(z_scores_ev = case_when(is.na(z_score_ev) ~ 0,
                                        z_score_ev == 'NaN' ~ 0,
                                        z_score_ev == Inf ~ 10000,
                                        TRUE ~ z_score_ev)) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_ev) >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

asv_tab_all_bloo_z_tax$class_f |>
  unique()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, class_f) |>
  dplyr::mutate(abund_class = sum(abundance_value)) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
  geom_area(aes(date, abund_class, fill = class_f, group = class_f), alpha = 0.8,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  geom_point(data = community_eveness_all_m |>
               dplyr::filter(anomaly_color == '#9F0011'),  
             aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1),
                     sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
  scale_color_identity()+
  scale_fill_manual(values = palette_class_assigned_bloo, na.value = "#000000")+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
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

##color by order
asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, order_f) |>
  dplyr::mutate(abund_order = sum(abundance_value)) |>
  ungroup() |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
  geom_area(aes(date, abund_order, fill = order_f, group = order_f), alpha = 0.8,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  geom_point(data = community_eveness_all_m |>
               dplyr::filter(anomaly_color == '#9F0011'),  
             aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1),
                     sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
  scale_color_identity()+
  scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
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

#### Another idea with the same time of plot is to plot only blooming events just to observe the events only not their dynamics over the years----
##color by order
asv_tab_all_bloo_z_tax |>
  colnames()

## I select only abundances for the anomalies of each ASV NOT the entire dataset for each ASVs that have an anomaly over the hole timeline.
asv_tab_all_bloo_z_tax %$%
  abundance_type |>
  unique()

asv_tab_all_bloo_z_tax_relabund <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance')

 asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'zclr' & 
                 abundance_value >= 1.96) |>
  group_by(asv_num, date) |>
  distinct(asv_num, date) |>
   left_join(asv_tab_all_bloo_z_tax_relabund, by = c('asv_num', 'date'))|>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  group_by(date, fraction, order_f) |>
  dplyr::mutate(abund_order = sum(abundance_value)) |>
  ungroup() |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
  geom_area(aes(date, abund_order, fill = order_f, group = order_f), alpha = 0.8,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  geom_point(data = community_eveness_all_m |>
               dplyr::filter(anomaly_color == '#9F0011'),  
             aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1),
                     sec.axis = sec_axis(~.* 1 , name = 'Community Evenness'))+
  scale_color_identity()+
  scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Family')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
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
# asv_tab_all_bloo |>
#   left_join(z_scores_all) |>
#   left_join(tax_bbmo_10y, by = 'asv_num') |>

asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  # mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
  #                                  z_score_ra == 'NaN' ~ 0,
  #                                  z_score_ra == Inf ~ 0,
  #                                  TRUE ~ z_score_ra)) |>
  # dplyr::mutate(anomaly_color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  #dplyr::filter(z_score >= 1.96) |>
  #ungroup() |>
  #left_join(m_bbmo_10y, by = 'sample_id') |>
  #left_join(tax, by = c('asv_num' = 'asv')) |>
  #dplyr::filter(class != is.na(class)) |> ##Raro tenir NAs a Class i que no estiguin filtrats?
  # ggplot(aes(date, abundance_value, color = ifelse(is.na(z_score_ra), "#080808", 
  #            if_else(z_score_ra >= 1.96, '#9F0011', '#080808'))))+ #, color = 'Class' , shape = class 
  #ggplot(aes(date, abundance_value, color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')))+
  dplyr::filter(abundance_type == 'relative_abundance') |>
  ggplot(aes(season, abundance_value))+  
  geom_violin(aes(season, abundance_value, group = season), position_jitter(width = 0,3) )+
  #scale_x_datetime()+
  #facet_wrap(vars(class), scales = 'free')+
  #geom_smooth(aes(month, abundance_value))+
  #geom_line(aes(group = asv_num), color = '#080808')+ #, color = '#3D3B3B'
  geom_point(aes(color = class), alpha = 1)+ #aes(shape = class)
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Month', y = 'Relative abundance (%)', color = 'Class')+ #, shpae = 'Class'
  #scale_color_identity()+
  scale_color_manual(values = palette_class_assigned_bloo)+
  #scale_color_manual(values = if_else(asv_tab_z_scores_all$z_score_ra >= 1.96,  '#9F0011', '#080808', missing = '#080808'))+
  #facet_wrap(fraction~abundance_type, scales = 'free')+
  #facet_grid(fraction~class, scales = 'free')+
  #guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), legend.position = 'bottom', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), strip.background = element_blank())

##plot all seasons together (seasonal anomalies?)
asv_tab_all_bloo_z_tax |>
  # left_join(z_scores_all) |>
  # left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  # dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  # mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
  #                                  z_score_ra == 'NaN' ~ 0,
  #                                  z_score_ra == Inf ~ 0,
  #                                  TRUE ~ z_score_ra)) |>
  # dplyr::mutate(anomaly_color = if_else(z_score_ra_ed >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  #dplyr::filter(z_score_ra_ed >= 1.96) |>
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

## Occurrence of this AVSs vs frequency of blooming events and magnitude of the events --------
###calculate occurrency
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

anom_perc <- anom_perc |>
  left_join(tax_bbmo_10y_old) |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num))

anom_perc$class_f <-  factor(anom_perc$class_f, 
                                          levels=unique(anom_perc$class_f[order(anom_perc$phylum_f)]), 
                                          ordered=TRUE)

anom_perc$order_f <-  factor(anom_perc$order_f, 
                                          levels=unique(anom_perc$order_f[order(anom_perc$phylum_f,
                                                                                             anom_perc$class_f)]), 
                                          ordered=TRUE)

anom_perc$family_f <-  factor(anom_perc$family_f, 
                                           levels=unique(anom_perc$family_f[order(anom_perc$phylum_f,
                                                                                               anom_perc$class_f,
                                                                                               anom_perc$order_f)]), 
                                           ordered=TRUE)


anom_perc$asv_num_f <-  factor(anom_perc$asv_num_f, 
                                            levels=unique(anom_perc$asv_num_f[order(anom_perc$phylum_f,
                                                                                                 anom_perc$class_f,
                                                                                                 anom_perc$order_f,
                                                                                                 anom_perc$family_f)]), 
                                            ordered=TRUE)

anom_perc |>
  ggplot(aes(interaction(family_f, asv_num_f), anom_perc, fill = order_f))+
  geom_col()+
  facet_wrap(vars(fraction), nrow = 1)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))



# EXPLORATION OF THE RELATIONSHIP BETWEEN BLOOMING EVENTS AND COMMUNITY ALTERATION----
## when blooming events happen we need the community evenness to be lower than 50%? Also high beta diversity?
### We explore this parammeters and see if we observe a pattern.
community_eveness_all_m |>
  dplyr::filter(community_eveness_rar < 0.5)

bray_curtis_rar_all_m |>
  dplyr::filter(bray_curtis_result > 0.9)

bray_curtis_rar_all_m_z <- bray_curtis_rar_all_m |>
  right_join(z_scores_bray) #, by = c('samples' = 'sample_id' )

bray_curtis_rar_all_m_Z$z_score_bray

bray_curtis_rar_all_m |>
  dim()

community_eveness_all_m |>
  left_join(bray_curtis_rar_all_m_Z) |>
  ggplot(aes(bray_curtis_result, community_eveness_rar))+
  geom_point(shape = ifelse(community_eveness_all_m$z_score_ev >= 1.96, 16, 17), 
             fill = ifelse(bray_curtis_rar_all_m_z$z_score_bray >= 1.69, '#B31722', '#000000'))+
  scale_color_identity()+
  theme_bw()

## Which is the mean evenness for all the dataset?
###I want to highlight the dates where there is a potential blooming event.
### What we do here is, when there's an anomaly in the relative abundance of an ASV, if we sum all
### ASVs relative abundance's that present this anomaly, do we observe a decrease in the community 
### Evenness?
asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  distinct(sample_id)

samples_id <- asv_tab_all_bloo_z_tax |>
  dplyr::select(sample_id) |>
  unique()

bloom_events_asvs <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(sample_id, asv_num, abundance_type, abundance_value) |>
  group_by(sample_id) |>
  dplyr::mutate(abund_anom = sum(abundance_value)) |>
  group_by(sample_id) |>
  dplyr::summarize(n_asv_bloom = n(), abund_anom = unique(abund_anom)) |>
  right_join(samples_id, by = 'sample_id') |>
  dplyr::mutate(n_asv_bloom = case_when(is.na(n_asv_bloom) ~0,
                                        !is.na(n_asv_bloom) ~ n_asv_bloom), 
                abund_anom = case_when(is.na(abund_anom) ~ 0,
                                       !is.na(abund_anom) ~ abund_anom))

bloom_events_asvs |>
  slice_max(order_by = n_asv_bloom, n = 5)

bloom_events_asvs |>
  summarize(max(n_asv_bloom), min(n_asv_bloom))

bloom_events_asvs$n_asv_bloom

# community_eveness_all_m |>
#   colnames()

palete_gradient_cb <- c(#"#240023",
  
                        "#4db2a2",
                        "#005a47" = 1,
                        na.value = '#000000') 

community_eveness_all_m |>
  left_join(bray_curtis_rar_all_m) |>
  left_join(bloom_events_asvs) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity', values_to = 'values') |>
  dplyr::filter(diversity != 'bray_curtis_result') |>
  ggplot(aes(y = diversity, x = values, color = as.numeric(abund_anom)))+
  geom_point(aes(size = as.numeric(n_asv_bloom)), position = position_jitter(0.03))+
  geom_violin(aes(group = diversity), alpha= 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_x_continuous(limits = c(0.3, 1.0), expand = c(0,0))+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  scale_color_gradientn(colors = palete_gradient_cb)+
  #scale_size( range = c(1,14))+
  scale_size_continuous(
                        breaks = seq(min(bloom_events_asvs$n_asv_bloom), max(bloom_events_asvs$n_asv_bloom), length.out = 3),
                        labels = round(seq(min(bloom_events_asvs$n_asv_bloom), max(bloom_events_asvs$n_asv_bloom), length.out = 3), 0))+ #, labels = c("Min Value", "Max Value"
guides(shape = guide_legend(ncol = 3, size = 1),
       size = guide_legend(ncol = 2,
                           override.aes = aes(label = '')))+
  labs(x = 'Community Evenness', 
       color = 'Sum of the relative abundance\nof all ASVs potentially blooming\nin that sample',
       y = '', size = 'Number of ASVs\npresenting a potential blooming event')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), legend.position = 'bottom', 
        text = element_text(size = 7),
        aspect.ratio = 3/9) #axis.text.x = element_blank()

## here the tendency is less clear, probably because I 
bray_curtis_rar_all_m |>
  left_join(community_eveness_all_m) |>
  left_join(bloom_events_asvs) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity', values_to = 'values') |>
  dplyr::filter(diversity != 'community_eveness_rar') |>
  ggplot(aes(y = diversity, x = values, color = as.numeric(abund_anom)))+
  geom_point(aes(size = as.numeric(n_asv_bloom)), position = position_jitter(0.03))+
  geom_violin(aes(group = diversity), alpha= 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_x_continuous(limits = c(0.3, 1.0), expand = c(0,0))+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  scale_color_gradientn(colors = palete_gradient_cb)+
  #scale_size( range = c(1,14))+
  scale_size_continuous(
    breaks = seq(min(bloom_events_asvs$n_asv_bloom), max(bloom_events_asvs$n_asv_bloom), length.out = 3),
    labels = round(seq(min(bloom_events_asvs$n_asv_bloom), max(bloom_events_asvs$n_asv_bloom), length.out = 3), 0))+ #, labels = c("Min Value", "Max Value"
  guides(shape = guide_legend(ncol = 3, size = 1),
         size = guide_legend(ncol = 2,
                             override.aes = aes(label = '')))+
  labs(x = 'Bray Curtis dissimilarity', 
       color = 'Sum of the relative abundance\nof all ASVs potentially blooming\nin that sample',
       y = '', size = 'Number of ASVs\npresenting a potential blooming event')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), legend.position = 'bottom', 
        text = element_text(size = 7),
        aspect.ratio = 3/9)

## NMDS, clustering. Do blooming events cluster together?-----
#### which samples present a potential blooming event?
samples_with_bloom_event <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(sample_id, asv_num, abundance_type, abundance_value) |>
  group_by(sample_id) |>
  dplyr::mutate(abund_anom = sum(abundance_value)) |>
  group_by(sample_id) |>
  dplyr::summarize(n_asv_bloom = n(), abund_anom = unique(abund_anom)) |>
  right_join(samples_id, by = 'sample_id') |>
  group_by(sample_id) |>
  slice_max(order_by = abund_anom, n = 1) |>
  ungroup() |>
  dplyr::filter(!is.na(n_asv_bloom)) 

#### which is the maximal relative abundance of the potential bloomer/ event
max_abund_asv_bloom_events <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(sample_id, asv_num, abundance_type, abundance_value) |>
  group_by(sample_id) |>
  slice_max(order_by = abundance_value, n = 3)
  
max_abund_asv_bloom_events |>
  dim() == samples_with_bloom_event |>
  dim() ##ROW should be TRUE (row - cols)

samples_with_bloom_event_asv_abund <- samples_with_bloom_event |>
  right_join(max_abund_asv_bloom_events)

## I would like to know more about these blooming events.
### Are they more present in PA or FL? Which is their taxonomy?----
new_tax <-  readRDS('~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/03_tax_assignation/devotes_all_assign_tax_assignation_v2.rds') |>
  as_tibble(rownames = 'sequence')

tax_bbmo_10y_new <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax, by = c('seq' = 'sequence')) |>
  rename(domain = Kingdom, phylum = Phylum, class = Class, order = Order, family = Family, genus = Genus) |>
  dplyr::mutate(phylum_f = as_factor(phylum),
              family_f = as_factor(family),
              order_f = as_factor(order),
              class_f = as_factor(class),
              asv_num_f = as_factor(asv_num))

tax_bbmo_10y_new$class_f <-  factor(tax_bbmo_10y_new$class_f, 
                                          levels=unique(tax_bbmo_10y_new$class_f[order(tax_bbmo_10y_new$phylum_f)]), 
                                          ordered=TRUE)

tax_bbmo_10y_new$order_f <-  factor(tax_bbmo_10y_new$order_f, 
                                          levels=unique(tax_bbmo_10y_new$order_f[order(tax_bbmo_10y_new$phylum_f,
                                                                                             tax_bbmo_10y_new$class_f)]), 
                                          ordered=TRUE)

tax_bbmo_10y_new$family_f <-  factor(tax_bbmo_10y_new$family_f, 
                                           levels=unique(tax_bbmo_10y_new$family_f[order(tax_bbmo_10y_new$phylum_f,
                                                                                               tax_bbmo_10y_new$class_f,
                                                                                               tax_bbmo_10y_new$order_f)]), 
                                           ordered=TRUE)


tax_bbmo_10y_new$asv_num_f <-  factor(tax_bbmo_10y_new$asv_num_f, 
                                            levels=unique(tax_bbmo_10y_new$asv_num_f[order(tax_bbmo_10y_new$phylum_f,
                                                                                                 tax_bbmo_10y_new$class_f,
                                                                                                 tax_bbmo_10y_new$order_f,
                                                                                                 tax_bbmo_10y_new$family_f)]), 
                                            ordered=TRUE)



blooming_events_fraction <- samples_with_bloom_event_asv_abund |>
  dplyr::mutate(fraction = case_when(str_detect(sample_id, '0.2_') ~ 0.2,
                                     str_detect(sample_id, '3_') ~ 3)) |>
  group_by(fraction, asv_num) |>
  dplyr::mutate(fraction_blooms_asv = n()) |>
  group_by(fraction) |>
  dplyr::mutate(fraction_blooms = n()) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 

blooming_events_fraction  |>
  colnames()

## create a hierarchical clustering to add to the heat_map----

## heat_map
blooming_events_fraction |>
  ggplot(aes( interaction(asv_num_f, family_f), date))+
  geom_tile(aes(fill = abundance_value))+
  scale_y_datetime(expand = c(0,0))+
  facet_wrap(vars(fraction.y), labeller = labs_fraction)+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  coord_flip()+
  labs(y = 'Date', x = 'ASV number')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0), text = element_text(size = 7))

## number of blooming events/fraction and their taxonomy----

unique(blooming_events_fraction$fraction_blooms)

unique(blooming_events_fraction$family_f)

### plot in columns colored by family
blooming_events_fraction |>
  dplyr::select(family_f, fraction_blooms_asv, fraction_blooms, asv_num_f, fraction.y, class_f) |>
  dplyr::distinct(fraction_blooms_asv, family_f, fraction_blooms, asv_num_f, fraction.y, class_f) |>
  group_by(family_f, fraction_blooms_asv) |>
  dplyr::mutate(sum_blooms_fam = sum(fraction_blooms_asv)) |>
  distinct(family_f, fraction_blooms, sum_blooms_fam, fraction.y, class_f) |>
  #mutate(abund_anom = sum(abundance_value)/fraction_blooms) |>
  #distinct(family_f, abund_anom, fraction_blooms, fraction.y) |>
  ggplot(aes(fraction.y, sum_blooms_fam, fill = family_f, label = fraction_blooms))+
  geom_col()+
  #geom_bar(stat = 'identity')+
  scale_x_discrete(labels = labs_fraction)+
  scale_y_continuous(expand = c(0,0), limits = c(0,93))+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  # geom_text(aes(fraction.y, sum_blooms_fam), unique(blooming_events_fraction$fraction_blooms, blooming_events_fraction$fraction.y), 
  #           color = "black") +  # Add text annotation for max value
  geom_text(nudge_y = 91, check_overlap = TRUE, size = 3)+
  scale_x_discrete(labels = labs_fraction) +
  labs(fill = "Family", y = 'Number of\npotential blooming events')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 7/3)

###plot by families 
blooming_events_fraction |>
  dplyr::select(family_f, fraction_blooms_asv, fraction_blooms, asv_num_f, fraction.y, class_f) |>
  dplyr::distinct(fraction_blooms_asv, family_f, fraction_blooms, asv_num_f, fraction.y, class_f) |>
  group_by(family_f, fraction_blooms_asv) |>
  dplyr::mutate(sum_blooms_fam = sum(fraction_blooms_asv)) |>
  distinct(family_f, fraction_blooms, sum_blooms_fam, fraction.y, class_f) |>
  #mutate(abund_anom = sum(abundance_value)/fraction_blooms) |>
  #distinct(family_f, abund_anom, fraction_blooms, fraction.y) |>
  ggplot(aes( sum_blooms_fam, interaction(family_f, class_f), fill = family_f, label = fraction_blooms))+
  geom_col()+
  #geom_bar(stat = 'identity')+
  #scale_x_discrete(labels = labs_fraction)+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_grid(cols = vars(fraction.y),  labeller = labs_fraction)+
  # geom_text(aes(fraction.y, sum_blooms_fam), unique(blooming_events_fraction$fraction_blooms, blooming_events_fraction$fraction.y), 
  #           color = "black") +  # Add text annotation for max value
  #geom_text(nudge_y = 91, check_overlap = TRUE, size = 3)+
  labs(fill = "Family", y = 'Number of\npotential blooming events')+
  guides(fill = guide_legend(ncol = 1, size = 3, keywidth = unit(2, 'mm'), keyheight = unit(2, 'mm')))+
  #coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank())

###plot by ASVs but colored by families 
blooming_events_fraction |>
  dplyr::select(family_f, fraction_blooms_asv, fraction_blooms, asv_num_f, fraction.y, class_f) |>
  dplyr::distinct(fraction_blooms_asv, family_f, fraction_blooms, asv_num_f, fraction.y, class_f) |>
  distinct(family_f, fraction_blooms_asv,  fraction.y, class_f, asv_num_f) |>
  #mutate(abund_anom = sum(abundance_value)/fraction_blooms) |>
  #distinct(family_f, abund_anom, fraction_blooms, fraction.y) |>
  ggplot(aes( fraction_blooms_asv, interaction(asv_num_f, family_f), fill = family_f, label = fraction_blooms_asv))+
  geom_col()+
  #geom_bar(stat = 'identity')+
  #scale_x_discrete(labels = labs_fraction)+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_grid(cols = vars(fraction.y),  labeller = labs_fraction)+
  # geom_text(aes(fraction.y, sum_blooms_fam), unique(blooming_events_fraction$fraction_blooms, blooming_events_fraction$fraction.y), 
  #           color = "black") +  # Add text annotation for max value
  #geom_text(nudge_y = 91, check_overlap = TRUE, size = 3)+
  labs(fill = "Family", y = 'Number of potential blooming events')+
  guides(fill = guide_legend(ncol = 1, size = 3, keywidth = unit(2, 'mm'), keyheight = unit(2, 'mm')))+
  #coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank())

### NMDS----
bbmo_10y@otu_table |>
  class()

row.names(asv_tab_bbmo_10y_w_rar) <- asv_tab_bbmo_10y_w_rar[,1]  

asv_tab_bbmo_10y_w_rar_ed <- asv_tab_bbmo_10y_w_rar[,-1]

data.hel <- asv_tab_bbmo_10y_w_rar_ed |>
  decostand(method="hellinger"); str(data.hel)

data.dist <- vegdist(data.hel, method="bray")
head(data.dist)
data.nmds<-metaMDS(data.dist)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.137 (per sota de 20; és acceptable)
data.nmds.points<-data.frame(data.nmds$points, Cluster = )  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)
data.nmds.points |>
  colnames()


# Create a data frame with NMDS coordinates and cluster information
nmds_data <- data.frame(nmds_result$points, Cluster = clusters)


nmds_bbmo_10y <- data.nmds.points |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble() |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  left_join(community_eveness_all, by = 'sample_id')

nmds_bbmo_10y |>
  colnames()

nmds_bbmo_10y |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar))+ # shape = fraction,
  geom_point(aes(color = season), alpha=3, shape = ifelse(nmds_bbmo_10y$community_eveness_rar <= 0.5, 16, 17))+
  #facet_grid(vars(year))+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  scale_color_manual(values = palette_seasons_4)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

#color year
palette_years <- c("#efd9ce",
                   "#d7b7c9",
                   "#be95c4",
                   "#af8ec2",
                   "#9f86c0",
                   "#7f6da7",
                   "#5e548e",
                   "#231942",
                   "#20173c",
                   "#1d1537")
nmds_bbmo_10y |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction))+ # shape = fraction,
  geom_point(aes(color = as.factor(year)), alpha=3)+
  #facet_grid(vars(year))+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

nmds_bbmo_10y |>
  ggplot(aes(MDS1, MDS2, shape = fraction, color = season))+
  geom_point(aes(color = season), alpha=3)+
  facet_grid(vars(fraction))+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  scale_color_manual(values = palette_seasons_4)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

nmds_bbmo_10y_02 <- nmds_bbmo_10y |>
  dplyr::filter(fraction == '0.2')

nmds_bbmo_10y |>
  dplyr::filter(fraction == '0.2') |>
  ggplot(aes(MDS1, MDS2))+ # shape = fraction,
  geom_point(aes(), alpha=3, shape = ifelse(nmds_bbmo_10y_02$community_eveness_rar <= 0.6, 16, 17))+#color = season
  facet_grid(vars(year))+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

## NMDS with data transformed to CLR (maybe we need to recalculate the transformation (some ASVs are removed)----
asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax %$%
  abundance_type 
  
## we need to transform the dataframe to wide format (this in only data from ASVs maybe I need the whole database transformed)
asv_tab_all_bloo_z_tax_w <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'zclr') |>
  dplyr::select(sample_id, asv_num, abundance_value) |>
  pivot_wider(id_cols = sample_id, names_from = asv_num, values_from = abundance_value, values_fill = 0) |>
  as.data.frame()

row.names(asv_tab_all_bloo_z_tax_w) <- asv_tab_all_bloo_z_tax_w[,1]  

asv_tab_all_bloo_z_tax_ed <- asv_tab_all_bloo_z_tax_w[,-1]

data.hel <- asv_tab_all_bloo_z_tax_ed |>
  decostand(method="hellinger"); str(data.hel)

data.dist <- vegdist(data.hel, method="euclidean", na.rm = TRUE)
head(data.dist)
data.nmds<-metaMDS(data.dist)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.0187 (per sota de 20; és acceptable)
data.nmds.points<-data.frame(data.nmds$points)  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)
data.nmds.points |>
  colnames()

nmds_bbmo_10y_bloo_zclr <- data.nmds.points |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble() |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  left_join(community_eveness_all, by = 'sample_id')

nmds_bbmo_10y_bloo_zclr |>
  colnames()

nmds_bbmo_10y_bloo_zclr |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction))+ # shape = fraction,
  geom_point(aes(color = as.factor(year)), alpha=3)+
  facet_wrap(vars(fraction), scales = 'free')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

nmds_bbmo_10y_bloo_zclr |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction))+ # shape = fraction,
  geom_point(aes(color = season), alpha=3)+
  facet_wrap(vars(fraction), scales = 'free')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  scale_color_manual(values = palette_seasons_4)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

nmds_bbmo_10y_bloo_zclr |>
  ggplot(aes(MDS1, MDS2, color = temperature, size = community_eveness_rar, shape = fraction))+ # shape = fraction,
  geom_point(aes(color = temperature), alpha=3)+
  facet_wrap(vars(fraction), scales = 'free')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  scale_color_gradientn(colors  = palete_gradient_cb)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

### NMDS with the whole community transformed to zCLR-----
asv_tab_10y_3_zclr |>
  colnames()

asv_tab_10y_02_zclr |>
  colnames()

#### We separate PA and FL because otherwise there are many NAs values and it can't calculate the distances.
##### 3
asv_tab_10y_3_zclr_w <- asv_tab_10y_3_zclr |>
  #bind_rows(asv_tab_10y_02_zclr) |>
  pivot_wider(id_cols = sample_id, names_from = asv_num, values_from = zclr, values_fill = 0) |>
  as.data.frame()

row.names(asv_tab_10y_3_zclr_w) <- asv_tab_10y_3_zclr_w[,1]  

asv_tab_10y_3_zclr_w_ed <- asv_tab_10y_3_zclr_w[,-1]

# my data is already transformed I don't need the transformation
# data.hel <- asv_tab_10y_3_zclr_w_ed |>
#   decostand(method="hellinger"); str(data.hel)

data.dist <- vegdist(asv_tab_10y_3_zclr_w_ed, method="euclidean", na.rm = TRUE)
head(data.dist)
data.nmds<-metaMDS(data.dist)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.16 (per sota de 20; és acceptable)
data.nmds.points<-data.frame(data.nmds$points)  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)
data.nmds.points |>
  colnames()

nmds_bbmo_10y_zclr_3 <- data.nmds.points |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble() |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  left_join(community_eveness_all, by = 'sample_id')

nmds_bbmo_10y_zclr_3 |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction))+ # shape = fraction,
  geom_point(aes(color = as.factor(year)), alpha=3)+
  facet_wrap(vars(fraction), scales = 'free')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

##### 0.2
asv_tab_10y_02_zclr_w <- asv_tab_10y_02_zclr |>
  #bind_rows(asv_tab_10y_02_zclr) |>
  pivot_wider(id_cols = sample_id, names_from = asv_num, values_from = zclr, values_fill = 0) |>
  as.data.frame()

row.names(asv_tab_10y_02_zclr_w) <- asv_tab_10y_02_zclr_w[,1]  

asv_tab_10y_02_zclr_w_ed <- asv_tab_10y_02_zclr_w[,-1]
# My data is already transformed to zCLR so I don't need the transformation 
# data.hel <- asv_tab_10y_02_zclr_w_ed |>
#   decostand(method=""); str(data.hel)

data.dist <- vegdist(asv_tab_10y_02_zclr_w_ed, method="euclidean", na.rm = TRUE)
head(data.dist)
data.nmds<-metaMDS(data.dist)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.13 (per sota de 20; és acceptable)
data.nmds.points<-data.frame(data.nmds$points)  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)
data.nmds.points |>
  colnames()

nmds_bbmo_10y_zclr_02 <- data.nmds.points |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble() |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  left_join(community_eveness_all, by = 'sample_id')

nmds_bbmo_10y_zclr_02 |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction))+ # shape = fraction,
  geom_point(aes(color = as.factor(year)), alpha=02)+
  facet_wrap(vars(fraction), scales = 'free')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

## When I transformed my ASVs data to zCLR some ASVs disappear, but using deconstand function they are kept (transformation rclr). To answer the question:----
## Do blooming events cluster together? 
### after the previous observations I observe that since PA and FL communities are well differentiated then I do the NMDS separated

asv_tab_bbmo_10y_l |> # upload the original table
  head()

asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_l |> #transform to wider fromat
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

data.dist <- vegdist(data.hel, method="euclidean", na.rm = TRUE)
head(data.dist)
data.nmds<-metaMDS(data.dist)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.15 (per sota de 20; és acceptable)
data.nmds.points<-data.frame(data.nmds$points)  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)
data.nmds.points |>
  colnames()

nmds_bbmo_10y_zclr_02_all <- data.nmds.points |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble() |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  left_join(community_eveness_all, by = 'sample_id')

nmds_bbmo_10y_zclr_02_all |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction))+ # shape = fraction,
  geom_point(aes(color = as.factor(year)), alpha=02)+
  facet_wrap(vars(fraction), scales = 'free')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

nmds_bbmo_10y_zclr_02_all |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction))+ # shape = fraction,
  geom_point(aes(color = as.factor(season)), alpha=02)+
  facet_wrap(vars(fraction), scales = 'free', labeller = labs_fraction)+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  scale_color_manual(values = palette_seasons_4)+
  #scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw()

### 3
data.hel <- asv_tab_bbmo_10y_w_3_ed |>
  decostand(method="rclr"); str(data.hel)

data.dist <- vegdist(data.hel, method="euclidean", na.rm = TRUE)
head(data.dist)
data.nmds<-metaMDS(data.dist)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.20 (per sota de 20; és acceptable)
data.nmds.points<-data.frame(data.nmds$points)  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)
data.nmds.points |>
  colnames()

nmds_bbmo_10y_zclr_3_all <- data.nmds.points |>
  rownames_to_column(var = 'sample_id') |>
  as_tibble() |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  left_join(community_eveness_all, by = 'sample_id')

nmds_bbmo_10y_zclr_3_all |>
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction, label = date))+ # shape = fraction,
  geom_point(aes(color = as.factor(year)), alpha=02)+
  facet_wrap(vars(fraction), scales = 'free')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  #scale_color_manual(values = palette_seasons_4)+
  scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw() 

### outlier BL091222_3_4022
nmds_bbmo_10y_zclr_3_all |>
  dplyr::filter(sample_id != 'BL091222_3_4022') |> #we remove the outlier to better observe the patterns
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = fraction, label = date))+ # shape = fraction,
  geom_point(aes(color = as.factor(season)), alpha=02)+
  facet_wrap(vars(fraction), scales = 'free', labeller = labs_fraction)+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  scale_color_manual(values = palette_seasons_4)+
  #scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw() 

### Conclusion: the community is basically structured by seasons

##Idea plot and highlight the blooming events on the NMDS.
nmds_bbmo_10y_zclr_3_all |>
  colnames()

samples_with_bloom_event ##list of ASVs that have a blooming event and the relative abundance of the anomaly and number of ASVs implicated

samples_with_bloom_event %$%
  n_asv_bloom |>
  range()

nmds_bbmo_10y_zclr_3_all |>
  dim()

nmds_bbmo_10y_zclr_3_all |>
left_join(samples_with_bloom_event, by = 'sample_id') |>
  dplyr::mutate(n_asv_bloom_ed = case_when(is.na(n_asv_bloom) ~ 0,
                                           n_asv_bloom == '1' ~ 1,
                                           n_asv_bloom == '2' ~ 1,
                                           n_asv_bloom == '3' ~ 1)) |>
  dplyr::filter(sample_id != 'BL091222_3_4022') |> #we remove the outlier to better observe the patterns
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = as.factor(n_asv_bloom_ed), label = date))+ # shape = fraction,
  geom_point(aes(color = as.factor(season)), alpha=02)+
  facet_wrap(vars(fraction), scales = 'free', labeller = labs_fraction)+
  labs(shape = 'Blooming event', size = 'Community\nevenness')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  scale_color_manual(values = palette_seasons_4)+
  #scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw() 

nmds_bbmo_10y_zclr_3_all |>
  colnames()

samples_with_bloom_event ##list of ASVs that have a blooming event and the relative abundance of the anomaly and number of ASVs implicated

samples_with_bloom_event %$%
  n_asv_bloom |>
  range()

nmds_bbmo_10y_zclr_3_all |>
  dim()

nmds_bbmo_10y_zclr_02_all |>
  left_join(samples_with_bloom_event, by = 'sample_id') |>
  dplyr::mutate(n_asv_bloom_ed = case_when(is.na(n_asv_bloom) ~ 0,
                                           n_asv_bloom == '1' ~ 1,
                                           n_asv_bloom == '2' ~ 1,
                                           n_asv_bloom == '3' ~ 1)) |>
  dplyr::filter(sample_id != 'BL091222_3_4022') |> #we remove the outlier to better observe the patterns
  ggplot(aes(MDS1, MDS2, color = season, size = community_eveness_rar, shape = as.factor(n_asv_bloom_ed), label = date))+ # shape = fraction,
  geom_point(aes(color = as.factor(season)), alpha=02)+
  facet_wrap(vars(fraction), scales = 'free', labeller = labs_fraction)+
  labs(shape = 'Blooming event', size = 'Community\nevenness')+
  #geom_text(aes(label=`Sample-Name`), check_overlap = TRUE, nudge_x = 0.02, nudge_y = 0.005)+
  scale_color_manual(values = palette_seasons_4)+
  #scale_color_manual(values = palette_years)+
  #scale_color_manual(values=palette_large)+
  theme_bw() 

####----

## per year (explore the % of blooming events per tax)----
asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax$month
  
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
  group_by(month, year, fraction) |>
  dplyr::summarize(n_anom_month = sum(n_anom)) |>
  ggplot(aes(month, n_anom_month, shape = fraction, group = month))+
 geom_point()+
  #geom_violin()+
  geom_line(aes(group = fraction))+
  scale_x_continuous(limits = c(0,12), expand = c(0,0))+
  #geom_line()+
  #geom_line(group = 'fraction')+
  #facet_grid(fraction~year)+
  #facet_grid(fraction~year)+
  facet_grid(vars(year))+
  #scale_y_continuous(labels = percent_format())+
  theme_bw()

anom_perc_02 |>
  bind_rows(anom_perc_3) |>
  group_by(month, year, fraction) |>
  dplyr::summarize(n_anom_month = sum(n_anom)) |>
  ggplot(aes(month, n_anom_month, shape = fraction, group = month))+
  geom_col()+
  #geom_violin()+
  #geom_line(aes(group = fraction))+
  scale_x_continuous(limits = c(0,12), expand = c(0,0))+
  #geom_boxplot(aes(group = month))+
  facet_grid(vars(fraction))+
  geom_smooth(aes(group = fraction), se = F)+
  #geom_line()+
  #geom_line(group = 'fraction')+
  #facet_grid(fraction~year)+
  facet_grid(fraction~year)+
  #facet_grid(vars(year))+
  coord_flip()+
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
  left_join(tax_bbmo_10y_new, by = 'asv_num') 

occurrence_perc_tax <-occurence_perc  |>
 left_join(anom_perc, by = c('asv_num','fraction')) |>
  left_join(max_rel, by = c('asv_num','fraction')) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num')

occurrence_perc_tax |>
  colnames()

occurrence_perc_tax %$%
  unique(genus)

occurrence_perc_tax %$%
  unique(family)

occurrence_perc_tax |>
  ggplot(aes(occurence_perc, anom_perc, color = class, size = max_rel*100, shape = fraction))+
  geom_point()+
  scale_color_manual(values = palette_class_assigned_bloo)+
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
  dplyr::mutate(project = '10Y_BBMO') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(z_score_ra != is.na(z_score_ra) &
                  z_score_ra != is.infinite(z_score_ra)) |>
  ggplot(aes(x = as.numeric(z_score_ra)))+
  geom_density(aes(group = project, fill = project))+
  scale_x_continuous(limits = c(-15, 50), expand = c(0,0))+
  labs(y = 'Density', x = 'z-scores')+
  geom_vline(xintercept =  1.96)+
  theme_bw()+
  theme(legend.position = 'none')
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

##I would like to plot the same but in this case for the whole dataset, so all z-scores calculated for all ASVs.

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

 
 
 