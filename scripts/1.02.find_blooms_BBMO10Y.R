# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# upload packages----
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
library(grateful) ## To cite the packages used.

# labels ----
labs_fraction <- as_labeller(c('0.2' = '0.2-3 µm',
                               '3' = '3-20 µm'))

labs_diversity <- as_labeller(c('community_eveness_rar' = 'Community Eveness', 
                                'bray_curtis_result' = 'Bray-Curtis dissimilarity'))

# highlight harbor restoration period ----
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |> ### The remodelation of the Blanes harbour strated on 24th March 2010 and finished on the 9th of june 2012
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

# palettes ----
palette_seasons_4 <- c("winter" = "#002562", 'spring' = "#519741", 'summer' = "#ffb900",'autumn' =  "#96220a")

palette_fraction <- c('0.2' = '#00808F', '3' = '#454545')


palette_gradient <- c('#DBDBDB', "#BDBEBE","#545454",
                      "#070607")
palette_gradient <- colorRampPalette(c('#DBDBDB', '#BDBEBE', '#545454', '#070607'))

## palettes taxonomy assigned ----
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
  "Clade II" = '#B0413E',  "Clade I" = "#ca6094", 
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

# functions----
source('../../Bloomers/R/get_anomalies.R') 
source('../../Bloomers/R/find_asv_with_anomalies.R')
source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R')

# Upload data and prepare it for the analysis ----
setwd("~/Documentos/Doctorat/BBMO/BBMO_bloomers/")

bbmo_10y <-readRDS("data/blphy10years.rds") ## 8052 all samples, no filtering

bbmo_10y <-
  prune_taxa(taxa_sums(bbmo_10y@otu_table) >0, ##filter ASVs that are 0 in the whole dataset
             bbmo_10y)

bbmo_10y |>
  nsamples() #237 samples 

bbmo_10y |>
  ntaxa() #8052 but if we filter for those that are 0 during the whole dataset then we got 7,849 ASVs

## separate datasets by ASV_tab, taxonomy and metadata
asv_tab_bbmo_10y_l <- bbmo_10y@otu_table |>
  as_tibble()

tax_bbmo_10y_old <- bbmo_10y@tax_table |> 
  as_tibble()

m_bbmo_10y <- bbmo_10y@sam_data |>
  as_tibble()

m_bbmo_10y |>
  colnames()

#write.csv(m_bbmo_10y, 'data/env_data/m_bbmo_10y.csv') 

## tidy colnames ----
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
                           "PNF2_5um_Micro", "PNF_5um_Micro", "cryptomonas", "micromonas",        
                           "HNF_Micro", "HNF2_5um_Micro", "HNF_5um_Micro", "LNA",               
                           "HNA", "prochlorococcus_FC", "Peuk1",  "Peuk2",          
                           "year", "month", "day", "season",            
                           "bacteria_joint", "synechococcus", "depth", "name_complete", 'sample_id')

# update taxonomy: new taxonomy created with the database SILVA 138.1 using Assign tax at 50 (default)
new_tax <-  readRDS('data/03_tax_assignation/devotes_all_assign_tax_assignation_v2.rds') |>
  as_tibble(rownames = 'sequence')

tax_bbmo_10y_new <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax, by = c('seq' = 'sequence'))

tax_bbmo_10y_new <- tax_bbmo_10y_new |>
  dplyr::select(domain = Kingdom, phylum = Phylum, class = Class,
                order = Order, family = Family, genus = Genus, asv_num, seq)

## Divide metadata into FL and PA ----
m_02 <- m_bbmo_10y  |>
  dplyr::filter(fraction == '0.2') 

m_02 <- m_02 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02)))

m_3 <- m_bbmo_10y |>
  dplyr::filter(fraction == 3.0)

m_3 <- m_3 |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3)))

## we lack of 3 samples in PA fraction which are [1] "2004-02-23" "2004-05-25" "2005-03-09" MISSING SAMPLES 

## --- Calculate relative abundance ----
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

### unique and exclusive ASVs in PA or FL fraction ----
asv_tab_10y_02_rel |>
  dplyr::filter(relative_abundance > 0 ) |>
  ungroup() |>
  distinct(asv_num) |>
  summary(n()) #3,390

asv_tab_10y_3_rel |>
  dplyr::filter(relative_abundance > 0 ) |>
  ungroup() |>
  distinct(asv_num) |>
  summary(n()) #6,486

### exclusive ASVs in PA and FL fraction ----
asvs_02 <- asv_tab_10y_02_rel |>
  dplyr::filter(relative_abundance > 0 ) |>
  ungroup() |>
  distinct(asv_num) 

asvs_3 <- asv_tab_10y_3_rel |>
  dplyr::filter(relative_abundance > 0 ) |>
  ungroup() |>
  distinct(asv_num) 

asvs_02 |>
  anti_join(asvs_3) #1,363

asvs_3 |>
  anti_join(asvs_02) #4,459

asvs_3 |>
  bind_rows(asvs_02) |>
  dplyr::filter(duplicated(asv_num)) #2,027

### general plot: the whole community along the years ----
asv_tab_10y_l_rel |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  ungroup() |>
  dplyr::select(class) |>
  distinct() |>
  as.vector()

m_bbmo_10y <- m_bbmo_10y |>
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d"))

bbmo_community_phylums <- asv_tab_10y_l_rel |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  dplyr::filter(!is.na(phylum))   |>
  group_by(phylum, sample_id, fraction, date) |>
  dplyr::summarize(sum_rel = sum(relative_abundance)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, sum_rel))+
  geom_area(aes(fill = phylum, group = phylum), alpha = 0.9,  position='stack')+
  scale_fill_manual(values = palette_phylums_assigned)+
  scale_x_datetime(expand = c(0,0), 
                   breaks = seq(min(m_bbmo_10y$date), 
                                max(m_bbmo_10y$date), by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  

bbmo_community_phylums

# ggsave('bbmo_community_phylums.pdf', bbmo_community_phylums,
#        path = '~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/',
#        width = 230,
#                height = 180,
#                units = 'mm')

## --- Transform data to rCLR ----
## transform asv_tab into wider format to go into vegan package
asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_l |>
  pivot_wider(names_from = 'asv_num', values_from = 'reads', values_fill = 0) |>
  as.data.frame()

rownames(asv_tab_bbmo_10y_w) <- asv_tab_bbmo_10y_w$sample_id

asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_w[,-1]

# --------
# In Coenen 2020 they say that
# adding a pseudocount disproportionately affects rare taxa, where the magnitude of differences between samples may 
# be similar to the magnitude of the added pseudocount and therefore obscured.
# However, when using the rCLR which is similar to the CLR it allows data that contains zeroes. This method does not use
# pseudocounts, unlike the standard CLR. Robust clr divides the vales by geometric mean o the observed features; zero values 
# are kept as zeroes, and not taken into account. In high dimensional data, the geometric mean of rclr is a good approximation
# of the true geometric mean (from deconstand documentation)
# ---------

### after checking the results from the different approximations we keep the rCLR transformation.

rclr_df <- decostand(asv_tab_bbmo_10y_w, method = 'rclr' ) # The decostand function will standardize the rows (samples) by default

##we create two datasets one for FL and one for PA
asv_tab_10y_02_rclr <- rclr_df |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'rclr') |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) 

asv_tab_10y_3_rclr <- rclr_df |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'rclr') |>
  dplyr::filter(str_detect(sample_id, '_3_'))

## creation of a complete dataset with all the normalizations performed (relative_abundances, pseudoabundances and rclr)----
asv_tab_10y_02_pseudo_rclr <- asv_tab_10y_02_rel|>
  left_join(asv_tab_10y_02_rclr, by = c('sample_id', 'asv_num')) |>
  dplyr::mutate(rclr = case_when(is.na(rclr)~ 0,
                                 !is.na(rclr) ~ rclr))

asv_tab_10y_3_pseudo_rclr <- asv_tab_10y_3_rel |>
  left_join(asv_tab_10y_3_rclr, by = c('sample_id', 'asv_num')) |>
  dplyr::mutate(rclr = case_when(is.na(rclr)~ 0,
                                 !is.na(rclr) ~ rclr))

nrow(asv_tab_10y_02_pseudo_rclr) == nrow(asv_tab_10y_02_rclr)
nrow(asv_tab_10y_3_pseudo_rclr) == nrow(asv_tab_10y_3_rclr)
#if FALSE something is WRONG!

#write.csv2(asv_tab_10y_02_pseudo_rclr, 'data/asv_tab_10y_02_pseudo_rclr.csv')
#write.csv2(asv_tab_10y_3_pseudo_rclr, 'data/asv_tab_10y_3_pseudo_rclr.csv')

# ---- Calculate diversity parameters ----
## We calculate different diversity parameters to check if the blooming events detected have an effect on the community structure
## Community Evenness ----
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

## Rarefied dataset to calculate Community eveness ----
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

### plot community evenness ----
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

## Bray Curtis dissimilarity ----
### 0 means the two sites have the same composition (that is they share all the species), and 1 means the two sites 
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

### plot bray curtis dissimilarity----
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
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity', color = 'Fraction')+
  guides(shape = 'none')+
  scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

### plot Bray-Curtis dissimilarity and Community Evenness together ----
community_eveness_all <- community_eveness_02 |>
  bind_rows(community_eveness_3) 

bray_curtis_rar_all <- bray_curtis_02_rar |> ##one sample less, the first one can't be compared with the previous
  bind_rows(bray_curtis_3_rar)

community_eveness_all_m <- community_eveness_all |>
  left_join(m_bbmo_10y, by = c('sample_id')) 

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
  facet_grid(vars(diversity_index), labeller = labs_diversity)+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Time', y = 'Community diversity', color = 'Fraction')+
  guides(shape = 'none')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

# ---- ## ------------- Discover anomalies ------------- ## ------------------------------------------------------------------------ 
## For each ASVs based on relative abundances and rCLR -----
time_lag_value <- 3 #define time lag number of samples before the bloom event considered
cut_off_value_ra <- 2.14 #define z-score cutoff 
cut_off_value_rclr <- 1 #define z-score cutoff 

z_02 <- asv_tab_10y_02_pseudo_rclr |>
  dplyr::group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> 
  dplyr::reframe(
    anomalies_ra = get_anomalies(time_lag = time_lag_value, negative = FALSE, 
                                 cutoff = cut_off_value_ra, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
    anomalies_clr = get_anomalies(time_lag = time_lag_value, negative = FALSE, 
                                  cutoff = cut_off_value_rclr, na_rm = TRUE, values = rclr, plotting = FALSE)[c(1,2,3)])

z_3 <- asv_tab_10y_3_pseudo_rclr |>
  dplyr::group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> 
  dplyr::reframe(
    anomalies_ra = get_anomalies(time_lag = time_lag_value, negative = FALSE, 
                                 cutoff = cut_off_value_ra, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)],
    anomalies_clr = get_anomalies(time_lag = time_lag_value, negative = FALSE, 
                                  cutoff = cut_off_value_rclr, na_rm = TRUE, values = rclr, plotting = FALSE)[c(1,2,3)])

# Filter the ASV_tab by only those ASVs that have an anomaly at some point of the dataset (IN THIS CASE WE ARE NOT FILTERING BY ANY RELATIVE ABUNDANCE IN THE COMMUITY JUST IF THEY PRESENT AN ANOMALY AT SOME POINT!)----
asv_anom_02 <- find_asv_with_anomalies(anomalies_result = z_02, 
                                       anomaly_in1 = anomalies_ra, 
                                       anomaly_in2 = NULL,  #anomalies_ps,
                                       anomaly_in3 = anomalies_clr, #anomalies_clr, 
                                       logic1 = 'TRUE',
                                       logic2 = NULL, 
                                       logic3 = 'TRUE',
                                       asv_col = asv_num)

## 21 apparent bloomers

asv_anom_3 <- find_asv_with_anomalies(anomalies_result = z_3, 
                                      anomaly_in1 = anomalies_ra, 
                                      anomaly_in2 = NULL, #anomalies_clr, 
                                      anomaly_in3 = anomalies_clr, #anomalies_ps
                                      logic1 = 'TRUE', 
                                      logic2 = NULL, 
                                      logic3 = 'TRUE', 
                                      asv_col = asv_num)
## 47 bloomers 

## Recover the z-score for each ASV at each time point of the time series  ----
### I want to highlight anomalies for each ASV to do so I recover z-scores for those ASVs that that have high z-scores
### at some point of the dataset. Easy to observe if those ASVs are having random anomalies or all of them happen at the same time
z_scores_02 <- asv_tab_10y_02_pseudo_rclr |>
  dplyr::group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::group_by(asv_num) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = time_lag_value, negative = FALSE, 
                                            cutoff = cut_off_value_ra, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  dplyr::group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num') 

z_scores_3 <- asv_tab_10y_3_pseudo_rclr |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::group_by(asv_num) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = time_lag_value, negative = FALSE, 
                                            cutoff = cut_off_value_ra, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  dplyr::group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num') 

z_scores_all <- z_scores_02 |>
  bind_rows(z_scores_3)

## anomalies based on rCLR abundances transfromed 
z_scores_02_rclr <- asv_tab_10y_02_pseudo_rclr |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::group_by(asv_num) |>
  dplyr::reframe(z_score_rclr = get_anomalies(time_lag = time_lag_value, negative = FALSE, 
                                              cutoff = cut_off_value_rclr, 
                                              na_rm = TRUE, values = rclr, 
                                              plotting = F)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_rclr) |>
  dplyr::group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_3_rclr <- asv_tab_10y_3_pseudo_rclr |>
  group_by(asv_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  dplyr::group_by(asv_num) |>
  dplyr::reframe(z_score_rclr = get_anomalies(time_lag = time_lag_value, negative = FALSE, cutoff = cut_off_value_rclr, 
                                              na_rm = TRUE, values = rclr, 
                                              plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_rclr) |>
  dplyr::group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num') 

z_scores_all_rclr <- z_scores_02_rclr |>
  bind_rows(z_scores_3_rclr)

## Now find potential bloomers ASVs that represent DURING the blooming event 10% of the community and a z-score of 1.96 ----
z_scores_02_red <- z_scores_02 |>
  dplyr::select(asv_num, z_score_ra, sample_id)

z_scores_02_red_rclr <- z_scores_02_rclr |>
  dplyr::select(asv_num, z_score_rclr, sample_id)

z_scores_02_red_all <- z_scores_02_red_rclr |>
  left_join(z_scores_02_red)

n_bloomers_02 <-  asv_tab_10y_02_pseudo_rclr |>
  dplyr::left_join(z_scores_02_red_all, by= c('sample_id', 'asv_num')) |>
  group_by(asv_num) |>
  dplyr::filter(relative_abundance >=  0.1 &
                  z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr) |>
  dplyr::distinct(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(sum = sum(n))

bloo_02 <- asv_tab_10y_02_pseudo_rclr |>
  dplyr::left_join(z_scores_02_red_all, by= c('sample_id', 'asv_num')) |>
  dplyr::select(z_score_ra, z_score_rclr, asv_num, relative_abundance) |>
  dplyr::filter((relative_abundance >=  0.1) &
                  (z_score_ra > cut_off_value_ra) &
                  (z_score_rclr > cut_off_value_rclr)) |>
  dplyr::ungroup() |>
  dplyr::distinct(asv_num) |>
  as_vector()

bloo_02 <- as_tibble_col(bloo_02, column_name = 'asv_num')

#write_csv2(as_tibble(bloo_02), 'data/detect_bloo/bloo_02.csv')

z_scores_3_red <- z_scores_3 |>
  dplyr::select(asv_num, z_score_ra, sample_id)

z_scores_3_red_rclr <- z_scores_3_rclr |>
  dplyr::select(asv_num, z_score_rclr, sample_id)

z_scores_3_red_all <- z_scores_3_red_rclr |>
  left_join(z_scores_3_red)

n_bloomers_3 <-  asv_tab_10y_3_pseudo_rclr |>
  dplyr::left_join(z_scores_3_red_all, by= c('sample_id', 'asv_num')) |>
  dplyr::filter(relative_abundance >=  0.1 &
                  z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr) |>
  dplyr::distinct(asv_num) |>
  dplyr::ungroup() |>
  dplyr::distinct(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(sum = sum(n))

bloo_3 <- asv_tab_10y_3_pseudo_rclr |>
  inner_join(asv_tab_10y_3_rclr, by = c('sample_id', 'asv_num')) |> 
  dplyr::left_join(z_scores_3_red_all, by= c('sample_id', 'asv_num')) |>
  dplyr::filter(relative_abundance >=  0.1 &
                  z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr) |>
  dplyr::ungroup() |>
  dplyr::distinct(asv_num) |>
  as_vector()

bloo_3 <- as_tibble_col(bloo_3, column_name = 'asv_num')

#write_csv2(as_tibble(bloo_3), 'data/detect_bloo/bloo_3.csv')

## CV needs to be greater than 1 ----
cv_02 <- asv_tab_10y_02_pseudo_rclr |>
  dplyr::filter(asv_num %in% bloo_02$asv_num) |>
  group_by(asv_num) |>
  dplyr::reframe(mean = mean(relative_abundance),
                 sd = sd(relative_abundance),
                 cv = sd(relative_abundance)/mean(relative_abundance)) |>
  dplyr::filter(cv < 1)

cv_3 <- asv_tab_10y_3_pseudo_rclr |>
  dplyr::filter(asv_num %in% bloo_3$asv_num) |>
  group_by(asv_num) |>
  dplyr::reframe(mean = mean(relative_abundance),
                 sd = sd(relative_abundance),
                 cv = sd(relative_abundance)/mean(relative_abundance)) |>
  dplyr::filter(cv < 1)

## bloomers tax table ---
bloo_02_tb <- bloo_02 |>
  dplyr::filter(!asv_num %in% c(cv_02$asv_num)) |>
  dplyr::mutate(fraction = '0.2')

bloo_3_tb <- bloo_3|>
  dplyr::mutate(fraction = '3')

bloo_tb_tax <- bloo_02_tb |>
  bind_rows(bloo_3_tb) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::select(-fraction) |>
  distinct() 

## bloom events table ----
bloom_events_3 <- asv_tab_10y_3_pseudo_rclr |>
  dplyr::filter(asv_num %in% bloo_3_tb$asv_num) |>
  left_join(z_scores_3_red_all) |>
  dplyr::filter(relative_abundance >=  0.1 &
                  z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr) |>
  distinct(sample_id, asv_num) |>
  left_join(m_bbmo_10y) |>
  dplyr::ungroup() |>
  dplyr::select(date, asv_num, fraction)

bloom_events_02 <- asv_tab_10y_02_pseudo_rclr |>
  dplyr::filter(asv_num %in% bloo_02_tb$asv_num) |>
  left_join(z_scores_02_red_all) |>
  dplyr::filter(relative_abundance >=  0.1 &
                  z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr) |>
  distinct(sample_id, asv_num) |>
  left_join(m_bbmo_10y) |>
  dplyr::ungroup() |>
  dplyr::select(date, asv_num, fraction)

bloom_events_tb <- bloom_events_3 |>
  bind_rows(bloom_events_02) |>
  dplyr::mutate(date = as.character(as.Date(date)))

#write.table(bloom_events_tb, 'results/tables/bloom_events_revised.txt', sep = '\t')

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

## Shared and exclusive bloomers for each fraction ----
asv_anom_3_tb <- bloo_3 |>
  as_tibble()

asv_anom_02_tb <- bloo_02_tb |>
  as_tibble()

common_bloomers_tax <- asv_anom_3_tb |>
  bind_rows(asv_anom_02_tb) |>
  unique() |> ##61 
  left_join(tax_bbmo_10y_new, by = c('asv_num' = 'asv_num'))

asv_anom_3_tb |>
  anti_join(asv_anom_02_tb) #41 ASV only in 3

asv_anom_02_tb |>
  anti_join(asv_anom_3_tb) #11 ASVs only in 0.2

# shared bloomers between fractions
asv_anom_3_tb |>
  bind_rows(asv_anom_02_tb) |>
  filter(duplicated(asv_num)) |> #6
  left_join(tax_bbmo_10y_new, by = c( 'asv_num')) |>
  dplyr::select(-seq)

asv_tab_10y_02_pseudo_rclr_bloo <-  asv_tab_10y_02_pseudo_rclr |>
  pivot_longer(cols = c(rclr, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  dplyr::filter(asv_num %in% bloo_02_tb$asv_num |
                  asv_num %in% bloo_3_tb$asv_num) ##recover ASVs that presented anomalies in 02 or 3

asv_tab_10y_02_pseudo_rclr_bloo |>
  group_by(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n())

asv_tab_10y_3_pseudo_rclr_bloo <- asv_tab_10y_3_pseudo_rclr |>
  group_by(asv_num) |>
  dplyr::filter(any(relative_abundance >=  0.1)) |> #estan en format 0-1
  pivot_longer(cols = c(rclr, relative_abundance), values_to = 'abundance_value', names_to = 'abundance_type') |>
  dplyr::filter(asv_num %in% asv_anom_02_tb$asv_num |
                  asv_num %in% asv_anom_3_tb$asv_num) ##recover ASVs that presented anomalies in 02 or 3

asv_tab_10y_3_pseudo_rclr_bloo |>
  group_by(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n())

asv_tab_all_bloo <- asv_tab_10y_02_pseudo_rclr_bloo |>
  bind_rows(asv_tab_10y_3_pseudo_rclr_bloo)

asv_tab_all_bloo |>
  colnames()

asv_tab_all_bloo |>
  group_by(asv_num) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n_num = n()) #58 ASVs

asv_tab_all_bloo |> 
  dplyr::filter(abundance_type == 'relative_abundance') %$%
  abundance_value |>
  range() #max is 0.56

## ----- Define relative abundance threshold, sensitivity test (defined 10%) -----
# Create a function to generate datasets based on the number of ASVs that fulfill the criteria of potential bloomers at different relative abundance thresholds.
source('src/count_number_potential_bloomers_threshold.R')

# I apply a loop for all the thresholds that I'm interested in
# Define a vector of threshold values
threshold_values <- c(0, 
                      0.00001, 0.000015, 0.000025, 0.00005, 0.000075, 
                      0.0001, 0.00015, 0.00025, 0.0005, 0.00075, 
                      0.001, 0.0015, 0.0025, 0.005, 0.0075, 
                      0.01, 0.015, 0.025, 0.05, 0.075,
                      0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                      0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75)

# Create an empty list to store the results
datasets <- list()

# Iterate over each threshold value and apply the function
for (threshold in threshold_values) {
  dataset_name <- paste0("n_0.2_", threshold * 100)  # Create dataset name
  dataset <- count_num_bloomers(threshold, 0.2, asv_tab_10y_02_pseudo_rclr, z_scores_tb = z_scores_02_red)
  datasets[[dataset_name]] <- dataset  # Store the dataset in the list
}

# Combine all datasets into a single dataframe
result_dataset_02 <- bind_rows(datasets)

# Create an empty list to store the results
datasets <- list()

# Iterate over each threshold value and apply the function
for (threshold in threshold_values) {
  dataset_name <- paste0("n_3_", threshold * 100)  # Create dataset name
  dataset <- count_num_bloomers(threshold, 3, asv_tab_10y_3_pseudo_rclr, z_scores_tb = z_scores_3_red)
  datasets[[dataset_name]] <- dataset  # Store the dataset in the list
}

result_dataset_3 <- bind_rows(datasets)

## supplementary figure 
blooming_threshold <- bind_rows(result_dataset_02,
                                result_dataset_3) |>
  ggplot(aes(threshold, num))+
  geom_point(size = 1)+
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = percent_format())+
  scale_y_log10()+
  facet_grid(vars(fraction), labeller = labs_fraction)+
  geom_vline(xintercept = 0.1, linetype = 'dashed')+
  labs(x = 'Relative abundance (%) threshold of the potential blooming event', y = 'Number of potential bloomers detected')+
  geom_line()+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        axis.ticks = element_blank())

blooming_threshold

# ggsave(blooming_threshold,  filename = 'blooming_threshold_log_ed3.pdf',
#        path = 'Results/Figures/',
#        width = 88, height = 88, units = 'mm')

## At the level of community, we use the Evenness result and Bray Curtis dissimilarity ----
z_diversity <- bray_curtis_02_rar |>
  dplyr::right_join(community_eveness_02, by = join_by("samples" == "sample_id")) |> 
  dplyr::reframe(anomalies_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = TRUE)[c(1,2,3)],# ),
                 anomalies_eveness = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = TRUE)[c(1,2,3)])

z_diversity %>%
  str()

## Recover z_scores of diversity----
z_scores_div_02_bray <- bray_curtis_02_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_div_02_ev <- community_eveness_02 |>
  dplyr::reframe(z_score_ev = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_ev) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_div_3_bray <- bray_curtis_3_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_div_3_ev <- community_eveness_3 |>
  dplyr::reframe(z_score_ev = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_ev) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_bray <- z_scores_div_02_bray |>
  bind_rows(z_scores_div_3_bray)

z_scores_ev <- z_scores_div_02_ev |>
  bind_rows(z_scores_div_3_ev)

# Dataset with all information ----
z_scores_all_red <- z_scores_all |>
  dplyr::ungroup() |>
  dplyr::select(asv_num, fraction, decimal_date, z_score_ra) |>
  dplyr::filter(asv_num %in% bloo_tb_tax$asv_num)

z_scores_all_rclr_red <- z_scores_all_rclr |>
  dplyr::ungroup() |>
  dplyr::select(asv_num, fraction, decimal_date, z_score_rclr)|>
  dplyr::filter(asv_num %in% bloo_tb_tax$asv_num)

asv_tab_all_bloo_z_tax <- asv_tab_all_bloo |>
  pivot_wider(values_from = abundance_value, names_from = abundance_type ) |>
  left_join(m_bbmo_10y, by =c('sample_id')) |>
  dplyr::select(asv_num, sample_id, reads, total_reads, relative_abundance, rclr, fraction, date, decimal_date, season) |>
  dplyr::ungroup() |>
  left_join(z_scores_all_red, by =c('asv_num', 'fraction', 'decimal_date')) |> 
  left_join(z_scores_all_rclr_red) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num')

#write.csv2(asv_tab_all_bloo_z_tax, 'data/asv_tab_all_bloo_z_tax_new_assign_zrclr.csv')
#write.csv2(asv_tab_all_bloo_z_tax, 'data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked_rclr.csv')

asv_tab_all_bloo_z_tax |>
  distinct(asv_num) |>
  dim()

asv_tab_all_bloo_z_tax |>
  colnames()

# ---- EXPLORE BLOOMING-LIKE BEHAVIOR IN OUR TIME SERIES ------------- ## -----------------------------------------------------------------------------
asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked_rclr.csv') |> ##using dada2 classifier assign tax with silva 138.1 and correctly identifying bloomers
  as_tibble() |>
  dplyr::select(-X)

## upload diversity data ----
## Rarefied dataset to calculate Community Evenness----
source('~/Documentos/Doctorat/Bloomers/R/community_evenness.R')

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

### plot Bray-Curtis dissimilarity and Community Evenness together----
community_eveness_all <- community_eveness_02 |>
  bind_rows(community_eveness_3)

bray_curtis_rar_all <- bray_curtis_02_rar |> ##one sample less, the first one can't be compared with the previous
  bind_rows(bray_curtis_3_rar)

## At the level of community, we use the Evenness result and Bray Curtis dissimilarity ----
z_diversity <- bray_curtis_02_rar |>
  dplyr::right_join(community_eveness_02, by = join_by("samples" == "sample_id")) |> 
  dplyr::reframe(anomalies_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = TRUE)[c(1,2,3)],# ),
                 anomalies_eveness = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = TRUE)[c(1,2,3)])

## Recover z_scores of diversity----
z_scores_div_02_bray <- bray_curtis_02_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_div_02_ev <- community_eveness_02 |>
  dplyr::reframe(z_score_ev = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_ev) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_div_3_bray <- bray_curtis_3_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_div_3_ev <- community_eveness_3 |>
  dplyr::reframe(z_score_ev = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = FALSE)[c(3)]) |>
  as_tibble() |>
  unnest(cols = z_score_ev) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_bray <- z_scores_div_02_bray |>
  bind_rows(z_scores_div_3_bray)

z_scores_ev <- z_scores_div_02_ev |>
  bind_rows(z_scores_div_3_ev)

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

asv_tab_all_bloo_z_tax$season <- asv_tab_all_bloo_z_tax$season |>
  factor(levels = c('winter', 'spring', 'summer', 'autumn'))

## Plot Evenness and Bray Curtis anomalies----
bray_curtis_02_rar |>
  bind_rows(bray_curtis_3_rar) |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, bray_curtis_result))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_point(aes(shape = fraction, color = fraction))+
  geom_line(aes(date, bray_curtis_result, group = fraction, color = fraction))+
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
  scale_x_datetime()+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity', color = 'Fraction')+
  guides(shape = 'none')+
  scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

### Plot blooming events with geom_area during the whole timeseries----
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  pivot_longer(cols = c('relative_abundance', 'rclr'), values_to = 'abundance_value', names_to = 'abundance_type')

bray_curtis_rar_all_m <- bray_curtis_rar_all |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result))

community_eveness_all_m <- community_eveness_all |>
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

#### plot a non overlaping area ----
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 

bray_curtis_rar_all_m <- bray_curtis_rar_all |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result))

community_eveness_all_m <- community_eveness_all |>
  dplyr::filter(community_eveness_rar != is.na(community_eveness_rar)) |>
  left_join(z_scores_ev) |>
  dplyr::mutate(z_scores_ev = case_when(is.na(z_score_ev) ~ 0,
                                        z_score_ev == 'NaN' ~ 0,
                                        z_score_ev == Inf ~ 10000,
                                        TRUE ~ z_score_ev)) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_ev) >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))


## try to explain that seasonal blooms are not always lead by the same ASV
labs_fraction_rec_freq <-  as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                                         '3' = 'Particle attached (3-20 um)',
                                         no = 'Recurrent',
                                         yes = 'Non-recurrent',
                                         seasonal = 'Seasonal',
                                         stochastic = 'Non-seasonal'))

summary_types_of_blooms <- bloo_all_types_summary_tax |>
  dplyr::mutate(fraction = as.character(fraction))

bloo_all_types_summary_tb <- bloo_all_types_summary_tax |>
  dplyr::mutate(fraction = as.character(fraction))

bloo_all_types_summary_tb$recurrency <- factor(bloo_all_types_summary_tb$recurrency, levels = c('recurrent', 'non-recurrent'))

summary_types_of_blooms |>
  dplyr::filter(recurrency == 'recurrent' &
                  type_of_bloomer == 'Seasonal')

booming_events_seas_BBMO10Y <- asv_tab_all_bloo_z_tax |>
  left_join(bloo_all_types_summary_tb, by = c('asv_num', 'fraction')) |>
  dplyr::filter(asv_num %in% bloo_02_tb$asv_num& fraction == '0.2' |
                  asv_num %in% bloo_3_tb$asv_num & fraction == '3' ) |>
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
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  
  geom_area(aes(date, abund_class, fill = order_f, group = order_f), alpha = 1,  position='stack')+
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
        panel.spacing.y = unit(1, "lines")) 

booming_events_seas_BBMO10Y

# ggsave('booming_events_BBMO10Y_new_tax_seasonal_ed3.pdf', booming_events_seas_BBMO10Y,
#        path = "results/figures/",
#        width = 180,
#        height = 160,
#        units = 'mm')

#bloo_all_types_summary_tb <- read.csv('results/tables/summary_types_of_blooms.csv')

## color by order
community_eveness_all_m <- community_eveness_all_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

bbmo_bloo_ev_order_no_sar_cluster <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num %in% bloo_02_tb$asv_num & fraction == '0.2' |
                  asv_num %in% bloo_3_tb$asv_num & fraction == '3' ) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
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
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

bbmo_bloo_ev_order_no_sar_cluster
# 
# ggsave('bbmo_bloo_ev_order_no_sar_cluster_ed2.pdf', bbmo_bloo_ev_order_no_sar_cluster,
#        path = '~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/',
#        width = 230,
#        height = 180,
#        units = 'mm')

### variability of bacterial abundances ----
m_bbmo_10y |>
  dplyr::select(bacteria_joint, decimal_date, type) |>
  distinct() |>
  dplyr::group_by(type) |>
  dplyr::reframe(mean = mean(bacteria_joint),
                 sd = sd(bacteria_joint),
                 se = sd(bacteria_joint) / sqrt(n()),
                 se_percent = (sd(bacteria_joint) / sqrt(n())) / mean(bacteria_joint) * 100)

## n ocurrences ---
summary_types_of_blooms |>
  distinct(asv_num, fraction, occurrence_category) |>
  dplyr::group_by(fraction, occurrence_category) |>
  dplyr::reframe(n = n())

### Identify those that are seasonal from those that are not seasonal----
asv_tab_all_bloo_z_tax$season <- asv_tab_all_bloo_z_tax$season |>
  factor(levels = c('winter', 'spring', 'summer', 'autumn'))

## I divide the different abundance values in separated plots and I filter only the potential bloomers in FL or PA-----
bloo_02 <- read.csv('data/detect_bloo/bloo_02.csv') |>
  as_tibble()

bloo_3 <-read.csv('data/detect_bloo/bloo_3.csv') |>
  as_tibble()

## preparation of the data for the plot
z_scores_all_red <- z_scores_all |>
  dplyr::select(asv_num, z_score_ra, sample_id, fraction, date, year)

## here we plot all potential bloomers for PA or FL and in which sampling points do they present an anomaly
blooming_events_rel_abund <- asv_tab_all_bloo_z_tax   |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  ggplot(aes(date, abundance_value, color = anomaly_color))+  
  scale_x_datetime(expand = c(0,0), 
                   breaks = seq(min(asv_tab_all_bloo_z_tax$date), 
                                max(asv_tab_all_bloo_z_tax$date), by = '1 year'),
                   date_labels = "%Y")+
  geom_hline(yintercept = 0.1, color = '#8B8989')+
  geom_line(aes(group = asv_num), color = '#5B5A5A')+ #, color = '#3D3B3B'
  geom_point(data = asv_tab_all_bloo_z_tax |>
               dplyr::filter(abundance_type == 'relative_abundance' &
                               z_score_ra >= cut_off_value_ra &
                               z_score_rclr >= cut_off_value_rclr &
                               abundance_value >= 0.1),  
             aes(date, abundance_value, color =  '#9F0011', alpha = 1))+ #size = z_score_ra
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  labs(x = 'Time', y = 'Relative abundance (%)', color = 'Anomaly')+ #, shpae = 'Class'
  scale_color_identity()+
  facet_wrap(vars(fraction), scales = 'fixed', ncol = 1, labeller = labs_fraction)+
  guides(fill = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size = 7),
        axis.ticks = element_blank(), legend.position = 'Bottom', axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7), strip.background = element_blank(),
        panel.border = element_blank())

# ggsave('blooming_events_rel_abund_v2.pdf', blooming_events_rel_abund,
#        path = '~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/',
#        width = 250,
#        height = 200,
#        units = 'mm')

## I plot rclr instead of relative abundances----
## preparation of the data for the plot
asv_tab_bloo_clr_z <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(z_score_ra_ed = case_when(is.na(z_score_ra) ~ 0,
                                          z_score_ra == 'NaN' ~ 0,
                                          z_score_ra == Inf ~ 0,
                                          TRUE ~ z_score_ra)) |>
  dplyr::mutate(anomaly_color = if_else((abundance_type == 'relative_abundance' &
                                                          z_score_ra_ed >= cut_off_value_ra &
                                                          z_score_rclr >= cut_off_value_rclr &
                                                          abundance_value >= 0.1),  '#9F0011', '#080808', missing = '#080808')) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(fraction == '0.2' & asv_num %in% bloo_02$asv_num |
                  fraction == '3' & asv_num %in% bloo_3$asv_num)

##blooming events per year colored by family
m_bbmo_10y <- m_bbmo_10y |>
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d"))
  
bloo_ev_year <- asv_tab_all_bloo_z_tax |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '0.2' & asv_num %in% bloo_02$asv_num |
                  fraction == '3' & asv_num %in% bloo_3$asv_num) |>
  dplyr::filter(z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |>
  dplyr::group_by(fraction, year) |>
  dplyr::summarize(n_year = n())

### mean events / year -----
asv_tab_all_bloo_z_tax |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '0.2' & asv_num %in% bloo_02_tb$asv_num |
                  fraction == '3' & asv_num %in% bloo_3_tb$asv_num) |>
  dplyr::filter(z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr &
                  abundance_value >= 0.1) |>
  distinct(fraction, year, date) |>
  dplyr::group_by(fraction, year) |>
  dplyr::summarize(n_year = n()) |>
  ungroup() |>
  dplyr::reframe(mean = mean(n_year), ## Bloom events / year
                 sd = sd(n_year))

m_bbmo_10y <- m_bbmo_10y |>
  dplyr::select(-date)

asv_tab_all_bloo_z_tax |>
  left_join(m_bbmo_10y) |>
  dplyr::mutate(bloom_event = case_when(abundance_type == 'relative_abundance' &
                                          abundance_value >= 0.1 &
                                          z_score_ra > cut_off_value_ra &
                                          z_score_rclr > cut_off_value_rclr ~ 1,
                                        TRUE ~ 0)) |>
  dplyr::distinct(bloom_event, date, fraction, decimal_date, year) |>
  dplyr::group_by(fraction, year) |>
  dplyr::reframe(n = sum(bloom_event)) |>
  group_by(fraction) |>
  dplyr::reframe(total_blooming_events = sum(n)) 

asv_tab_all_bloo_z_tax |>
  left_join(m_bbmo_10y) |>
  dplyr::mutate(bloom_event = case_when(abundance_type == 'relative_abundance' &
                                          abundance_value >= 0.1 &
                                          z_score_ra > cut_off_value_ra &
                                          z_score_rclr > cut_off_value_rclr ~ 1,
                                        TRUE ~ 0)) |>
  dplyr::distinct(bloom_event, date, fraction, decimal_date, year) |>
  dplyr::group_by(fraction, year) |>
  dplyr::reframe(n = sum(bloom_event)) |>
  group_by( fraction) |>
  reframe(mean_blooming_events = mean(n), ## Bloom events 0.2 & 3
          sd_blooming_events = sd(n))

asv_tab_all_bloo_z_tax |>
  left_join(m_bbmo_10y) |>
  dplyr::mutate(bloom_event = case_when(abundance_type == 'relative_abundance' &
                                          abundance_value >= 0.1 &
                                          z_score_ra > cut_off_value_ra &
                                          z_score_rclr > cut_off_value_rclr ~ 1,
                                        TRUE ~ 0)) |>
  dplyr::distinct(bloom_event, date, fraction, decimal_date, year) |>
  dplyr::group_by(fraction, year) |>
  dplyr::reframe(n = sum(bloom_event)) |>
  reframe(mean_blooming_events = mean(n),
          sd_blooming_events = sd(n))

asv_tab_all_bloo_z_tax |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '0.2' & asv_num %in% bloo_02$asv_num |
                  fraction == '3' & asv_num %in% bloo_3$asv_num) |>
  dplyr::filter(z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |>
  dplyr::group_by(fraction, year, family_f) |>
  dplyr::summarize(n_year_fam = n()) |>
  left_join(bloo_ev_year) |>
  dplyr::mutate(rel_ev_f = n_year_fam/n_year) |>
  ggplot(aes(year, n_year_fam))+  
  geom_col(aes(year, n_year_fam, fill = family_f))+ 
  labs(x = 'Year', y = 'Events / year', fill = 'Family')+ #, shpae = 'Class'
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(fraction), scales = 'fixed', labeller = labs_fraction, ncol = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), axis.title.x = element_text(size = 7),
        axis.ticks = element_blank(), legend.position = 'right', axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7), strip.background = element_blank())

## per month----
bloo_ev_month <- asv_tab_all_bloo_z_tax |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '0.2' & asv_num %in% bloo_02$asv_num |
                  fraction == '3' & asv_num %in% bloo_3$asv_num) |>
  dplyr::filter(z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr &
                  abundance_value >= 0.1) |>
  dplyr::group_by(fraction, month) |>
  dplyr::summarize(n_month = n())

asv_tab_all_bloo_z_tax |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '0.2' & asv_num %in% bloo_02$asv_num |
                  fraction == '3' & asv_num %in% bloo_3$asv_num) |>
  dplyr::filter(z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr &
                  abundance_value >= 0.1) |>
  dplyr::group_by(fraction, month, family_f) |>
  dplyr::summarize(n_month_fam = n()) |>
  left_join(bloo_ev_month) |>
  dplyr::mutate(rel_ev_f = n_month_fam/n_month) |>
  ggplot(aes(month, n_month_fam))+  
  geom_col(aes(month, n_month_fam, fill = family_f))+ 
  labs(x = 'month', y = 'Events / month', fill = 'Family')+ #, shpae = 'Class'
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(fraction), scales = 'fixed', labeller = labs_fraction, ncol = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), axis.title.x = element_text(size = 7),
        axis.ticks = element_blank(), legend.position = 'right', axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7), strip.background = element_blank())

## per season ----
bloo_ev_season <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '0.2' & asv_num %in% bloo_02$asv_num |
                  fraction == '3' & asv_num %in% bloo_3$asv_num) |>
  dplyr::filter(z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr &
                  abundance_value >= 0.1) |>
  dplyr::group_by(fraction, season) |>
  dplyr::summarize(n_season = n())

asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '0.2' & asv_num %in% bloo_02$asv_num |
                  fraction == '3' & asv_num %in% bloo_3$asv_num) |>
  dplyr::filter(z_score_ra >= cut_off_value_ra &
                  z_score_rclr >= cut_off_value_rclr &
                  abundance_value >= 0.1) |>
  dplyr::group_by(fraction, season, family_f) |>
  dplyr::summarize(n_season_fam = n()) |>
  left_join(bloo_ev_season) |>
  dplyr::mutate(rel_ev_f = n_season_fam/n_season) |>
  ggplot(aes(season, n_season_fam))+  
  geom_col(aes(season, n_season_fam, fill = family_f))+ 
  labs(x = 'Season', y = 'Events / season', fill = 'Family')+ #, shpae = 'Class'
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(fraction), scales = 'fixed', labeller = labs_fraction, ncol = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), axis.title.x = element_text(size = 7),
        axis.ticks = element_blank(), legend.position = 'right', axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7), strip.background = element_blank())

# Occurrence of this ASVs vs frequency of blooming events and magnitude of the events --------
###calculate occurrency
#asv_tab_all_bloo_z_tax <- read.delim2('~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/asv_tab_all_bloo_z_tax_new_assign.csv', sep = ';') 
nsamples_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3') %$%
  sample_id |>
  unique() |>
  length()

occurence_perc_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  !abundance_type %in% c('pseudoabundance', 'rclr')) |>
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
                  !abundance_type %in% c('pseudoabundance', 'rclr')) |>
  dplyr::filter(abundance_value > 0) |>
  group_by(asv_num, sample_id) |>
  summarize(n = n()) |>
  group_by(asv_num) |>
  summarize(occurence = sum(n)) |>
  mutate(occurence_perc = occurence/nsamples_02,
         fraction = '0.2')

occurence_perc <- occurence_perc_02 |>
  bind_rows(occurence_perc_3)

#### we can use number of samples in which that ASV has been detected
pres_s_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  !abundance_type %in% c('pseudoabundance', 'rclr')
                & asv_num %in% bloo_02$asv_num) |>
  dplyr::filter(abundance_value > 0) |>
  group_by(asv_num) |>
  dplyr::summarize(presence_asv = n())|>
  dplyr::mutate(fraction = '0.2')

pres_s_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  !abundance_type %in% c('pseudoabundance', 'rclr')
                & asv_num %in% bloo_3$asv_num) |>
  dplyr::filter(abundance_value > 0) |>
  group_by(asv_num) |>
  dplyr::summarize(presence_asv = n()) |>
  dplyr::mutate(fraction = '3')

pres_s <- pres_s_02 |>
  bind_rows(pres_s_3)

## ----------  Blooming Residuals ---------------
#### Here we would like to plot the difference between the mean abundance of an ASV and it's abundance during a blooming event.
source('src/calculate_and_plot_residuals.R')

library(purrr)
library(EnvStats)

## calculate the distance of each timepoint relative abundance from the geometric mean ----
mean_abund_f_asv <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(abundance_value = as.numeric(abundance_value),  # Convert abundance_value to numeric
                abundance_value = ifelse(is.na(abundance_value) | abundance_value == 0, 1, abundance_value * 100)) |> # geometric mean does not work with < 1 because it usees square root and it needs values > 1 
  dplyr::group_by(asv_num, fraction) |>
  dplyr::reframe(mean_abund = geoMean(abundance_value, na.rm = T),
                 sd_abund = geoSD(abundance_value, na.rm = T))

timeseries_limits <- tibble(xmin = '2004-01-26', xmax = '2013-10-15') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d"))

## recover the abundance during a blooming event
anom_rel_abund_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  abundance_type == 'relative_abundance' &
                  asv_num %in% bloo_3_tb$asv_num) |>
  dplyr::mutate(bloom_event = case_when(z_score_ra >= cut_off_value_ra &
                                          z_score_rclr >= cut_off_value_rclr &
                                          abundance_value >= 0.1 ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(date, sample_id, asv_num, fraction, abundance_value, family, bloom_event)

anom_rel_abund_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  abundance_type == 'relative_abundance' &
                  asv_num %in% bloo_02_tb$asv_num) |>
  dplyr::mutate(bloom_event = case_when(z_score_ra >= cut_off_value_ra &
                                          z_score_rclr >= cut_off_value_rclr &
                                          abundance_value >= 0.1 ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  dplyr::select(date, sample_id, asv_num, fraction, abundance_value, bloom_event, family, z_score_ra) 

anom_rel_abund <- anom_rel_abund_3 |>
  bind_rows(anom_rel_abund_02) |>
  arrange(-abundance_value) 

occurrence_bloo_bbmo_red <- occurrence_bloo_bbmo |>
  dplyr::select(asv_num, fraction, occurrence_perc) |>
  dplyr::mutate(fraction = as.double(fraction)) |>
  distinct()

## little function to compute residuals for all ASVs in PA or FL
bloo_02_filt <- bloo_02_tb 

summary_types_of_blooms_02 <-
  summary_types_of_blooms |>
  dplyr::filter(fraction == '0.2')

bloo_02_filt <-  bloo_02_filt |> 
  dplyr::mutate(fraction = as.double(fraction)) |>
  left_join(bloo_all_types_summary_tax, by = c('asv_num' = 'asv_num', 'fraction')) |>
  left_join(occurrence_bloo_bbmo_red) |>
  arrange(desc(occurrence_perc))

create_plot_02 <- function(asv_num) {
  plot <- plot_residuals(
    data_anom_abund = anom_rel_abund,
    data_mean_abund_asv = mean_abund_f_asv,
    asv_num =  {{asv_num}}, 
    community_fraction = '0.2'
  )
}

plots_list_02 <- purrr::map(bloo_02_filt$asv_num, create_plot_02)

residuals_02 <- gridExtra::grid.arrange(grobs = plots_list_02, ncol = 3) 

# ggsave("results/figures/residuals_02_plot_geo_mean_v2.pdf",
#        plot = residuals_02, width = 188, height = 220, units = "mm")

create_plot_3 <- function(asv_num) {
  plot_residuals(
    data_anom_abund = anom_rel_abund,
    data_mean_abund_asv = mean_abund_f_asv,
    asv_num = {{asv_num}},
    community_fraction = '3'
  )
}

summary_types_of_blooms_3 <- summary_types_of_blooms |>
  dplyr::filter(fraction == '3')

bloo_3 <-  bloo_3_tb |>
  dplyr::mutate(fraction = as.double(fraction)) |>
  left_join(bloo_all_types_summary_tax, by = c('asv_num' = 'asv_num', 'fraction')) |>
  left_join(occurrence_bloo_bbmo_red) |>
  arrange(desc(occurrence_perc))

bloo_3_1 <-  bloo_3 |>
  as_tibble() |>
  slice_head(n= 22)

bloo_3_2 <-  bloo_3 |>
  as_tibble() |>
  slice_tail(n = 23)

#i divide the PA community because there's too much ASVs for one plot
plots_list_3_1 <- map(bloo_3_1$asv_num, create_plot_3)

residuals_3 <- grid.arrange(grobs = plots_list_3_1, ncol = 3) 

# Assuming you have already created the plot and stored it in `residuals_3`

# Save the plot to a file
#ggsave("results/figures/residuals_3.1_plot_geo_mean_v2.pdf", plot = residuals_3, width = 188, height = 240, units = "mm")

#i divide the PA community because there's too much ASVs for one plot

plots_list_3_2 <- map(unique(bloo_3_2$asv_num), create_plot_3)

residuals_3 <- grid.arrange(grobs = plots_list_3_2, ncol = 3)

#ggsave("results/figures/residuals_3.2_plot_geo_mean_v2.pdf", plot = residuals_3, width = 188, height = 240, units = "mm")

bloo_02_tax <- bloo_02 |>
  left_join(tax_bbmo_10y_new)

#write.csv(bloo_02_tax, '../../Thesis/data/bloo_02_10y.csv', row.names = F)