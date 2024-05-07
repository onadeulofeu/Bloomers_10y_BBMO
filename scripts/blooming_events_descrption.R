## packages
library(tidyverse)
library(ggplot2)
library(magrittr)
library(scales)
library(phyloseq)
library(speedyseq)
library(Bloomers)
library(gridExtra)
#library(zCompositions)

##palettes
palete_gradient_cb <- c(#"#240023",
  "#4db2a2",
  "#005a47" = 1,
  na.value = '#000000') 

## functions
source('~/Documentos/Doctorat/Bloomers/R/community_evenness.R')

#upload data
asv_tab_all_bloo_z_tax <- read.csv2('detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') ##using dada2 classifier assign tax with silva 138.1


asv_tab_bbmo_10y_w_rar <- read.csv2('asv_tab_bbmo_10y_w_rar.csv') |> ## rarefied dataset for diversity analysis
  as_tibble()|> 
  rename('sample_id' = 'X') 

bloo_02 <- read.csv('~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/detect_bloo/bloo_02.csv') |>
  as_tibble()

bloo_3 <-read.csv('~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/detect_bloo/bloo_3.csv') |>
  as_tibble()

# z_scores_all <- read.csv('~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/z_scores_all.csv') |>
#   as_tibble()

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

## metadata
bbmo_10y <-readRDS("~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/blphy10years.rds") ##8052 asv amb totes les mostres, no está en % aquesta taula
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

tax_bbmo_10y_old |>
  str()

m_bbmo_10y <- bbmo_10y@sam_data |>
  as_tibble()

## tidy colnames----
colnames(asv_tab_bbmo_10y_l) <- c('asv_num', "sample_id", 'reads')

colnames(tax_bbmo_10y_old) <- c("asv_num", "kingdom", "phylum", "class", "order", "family", "genus",
                                "species", "curated", "otu_corr","seq")

colnames(m_bbmo_10y) <- c('sample_id', "project", "location", "code",             
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

#new taxonomy created with the database SILVA 138.1 using Assign tax at 50 (default)
new_tax <-  readRDS('~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/03_tax_assignation/devotes_all_assign_tax_assignation_v2.rds') |>
  as_tibble(rownames = 'sequence')

tax_bbmo_10y_old |>
  colnames()

tax_bbmo_10y_new <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax, by = c('seq' = 'sequence')) |>
  rename(domain = Kingdom, phylum = Phylum, class = Class, order = Order, family = Family, genus = Genus)

#labels----
labs_fraction <- as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                               '3' = 'Particle attached (3-20 um)'))

labs_diversity <- as_labeller(c('community_eveness_rar' = 'Community Eveness', 
                                'bray_curtis_result' = 'Bray-Curtis dissimilarity'))

## Taxonomic differences between whole community and the blooming community-----







## Relationship blooming events and community evenness----
### Rarefied dataset to calculate Community Evenness----
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

### plot Bray-Curtis dissimilarity and Community Eveness together----
community_eveness_all <- community_eveness_02 |>
  bind_rows(community_eveness_3)

bray_curtis_rar_all <- bray_curtis_02_rar |> ##one sample less, the first one can't be compared with the previous
  bind_rows(bray_curtis_3_rar)

## At the level of community, we use the Eveness result and Bray Curtis dissimilarity ----
z_diversity <- bray_curtis_02_rar |>
  dplyr::right_join(community_eveness_02, by = join_by("samples" == "sample_id")) |> 
  #ungroup() %>%
  #group_by(sample_id) %>%
  dplyr::reframe(anomalies_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = TRUE)[c(1,2,3)],# ),
                 anomalies_eveness = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = TRUE)[c(1,2,3)])

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


# Super blooming events -----
source('src/calculate_and_plot_residuals.R')
# asv_tab_bbmo_10y_w <- read.csv( 'data/.csv', 
#                                 row.names = 'X')
library(purrr)
library(EnvStats)

## calculate the distance of each timepoint relative abundance from the geometric mean ----
mean_abund_f_asv <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
 # dplyr::filter(abundance_value > 0) |>
  dplyr::mutate(abundance_value = as.numeric(abundance_value),  # Convert abundance_value to numeric
                abundance_value = ifelse(is.na(abundance_value) | abundance_value == 0, 1, abundance_value * 100)) |> # geometric mean does not work with < 1 because it usees square root and it needs values > 1 
  dplyr::group_by(asv_num, fraction) |>
  dplyr::reframe(mean_abund = geoMean(abundance_value, na.rm = T),
                   #sd_abund = geoSD(abundance_value, na.rm = T),
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
                  asv_num %in% bloo_3$value) |>
  # dplyr::filter(abundance_value >= 0.1 &
  #                 z_score_ra >= 1.96) |>
  dplyr::mutate(bloom_event = case_when(abundance_value > 0.1 &
                                          z_score_ra > 1.96 ~ 'bloom',
                                        abundance_value < 0.1 ~ 'no-bloom' ,
                                        z_score_ra < 1.96 ~ 'no-bloom')) |>
  dplyr::select(date, sample_id, asv_num, fraction, abundance_value, family, bloom_event)

anom_rel_abund_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  abundance_type == 'relative_abundance' &
                  asv_num %in% bloo_02$value) |>
  # dplyr::filter(abundance_value >= 0.1 &
  #                 z_score_ra >= 1.96) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(bloom_event = case_when(abundance_value > 0.1 &
                                          z_score_ra > 1.96 ~ 'bloom',
                                        abundance_value < 0.1 ~ 'no-bloom' ,
                                          z_score_ra < 1.96 ~ 'no-bloom')) |>
  dplyr::select(date, sample_id, asv_num, fraction, abundance_value, bloom_event, family, z_score_ra) 
  
anom_rel_abund <- anom_rel_abund_3 |>
  bind_rows(anom_rel_abund_02) |>
  arrange(-abundance_value) |>
  dplyr::mutate(abundance_value = abundance_value*100)
  
## little function to compute residuals for all ASVs in PA or FL
bloo_02_filt <- bloo_02 |>
  dplyr::filter(!value %in% c('asv2', 'asv3', 'asv5', 'asv8')) 

summary_types_of_blooms_02 <-
  summary_types_of_blooms |>
  dplyr::filter(fraction == '0.2')

bloo_02_filt <-  bloo_02_filt |>
  left_join(summary_types_of_blooms_02, by = c('value' = 'asv_num')) |>
  arrange(-occurrence_perc)

create_plot_02 <- function(asv_num) {
  plot <- plot_residuals(
    data_anom_abund = anom_rel_abund,
    data_mean_abund_asv = mean_abund_f_asv,
    asv_num = {{asv_num}},
    community_fraction = '0.2'
  )
}

plots_list_02 <- purrr::map(bloo_02_filt$value, create_plot_02)

residuals_02 <- gridExtra::grid.arrange(grobs = plots_list_02, ncol = 3) 

ggsave("results/figures/residuals_02_plot_geo_mean.pdf", 
       plot = residuals_02, width = 188, height = 220, units = "mm")

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

bloo_3 <-  bloo_3 |>
  left_join(summary_types_of_blooms_3, by = c('value' = 'asv_num')) |>
  arrange(-occurrence_perc)

bloo_3_1 <-  bloo_3 |>
  as_tibble() |>
  slice_head(n= 22)

bloo_3_2 <-  bloo_3 |>
  as_tibble() |>
  slice_tail(n= 23)

#i divide the PA community because there's too much ASVs for one plot
plots_list_3_1 <- map(bloo_3_1$value, create_plot_3)

residuals_3 <- grid.arrange(grobs = plots_list_3_1, ncol = 3) 

# Assuming you have already created the plot and stored it in `residuals_3`

# Save the plot to a file
ggsave("results/figures/residuals_3.1_plot_geo_mean.pdf", plot = residuals_3, width = 188, height = 220, units = "mm")

#i divide the PA community because there's too much ASVs for one plot

plots_list_3_2 <- map(unique(bloo_3_2$value), create_plot_3)

residuals_3 <- grid.arrange(grobs = plots_list_3_2, ncol = 3)

ggsave("results/figures/residuals_3.2_plot_geo_mean.pdf", plot = residuals_3, width = 188, height = 220, units = "mm")

## I do the same with mean (non geometric) because the geometric mean does not deal with 0 and it creates noise for those occurrently narrow ASVs------
mean_abund_f_asv <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::group_by(asv_num, fraction) |>
  dplyr::group_by(asv_num, fraction) |>
  dplyr::mutate(abundance_value = abundance_value*100) |>
 # dplyr::filter(abundance_value > 0) |>
  dplyr::reframe(mean_abund = mean(abundance_value, na.rm = T),
                 #sd_abund = geoSD(abundance_value, na.rm = T),
                 sd_abund = sd(abundance_value, na.rm = T))

timeseries_limits <- tibble(xmin = '2004-01-26', xmax = '2013-10-15') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d"))

## recover the abundance during a blooming event
anom_rel_abund_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  abundance_type == 'relative_abundance' &
                  asv_num %in% bloo_3$value) |>
  # dplyr::filter(abundance_value >= 0.1 &
  #                 z_score_ra >= 1.96) |>
  dplyr::mutate(bloom_event = case_when(abundance_value > 0.1 &
                                          z_score_ra > 1.96 ~ 'bloom',
                                        abundance_value < 0.1 ~ 'no-bloom' ,
                                        z_score_ra < 1.96 ~ 'no-bloom')) |>
  dplyr::select(date, sample_id, asv_num, fraction, abundance_value, family, bloom_event)

anom_rel_abund_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  abundance_type == 'relative_abundance' &
                  asv_num %in% bloo_02$value) |>
  # dplyr::filter(abundance_value >= 0.1 &
  #                 z_score_ra >= 1.96) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(bloom_event = case_when(abundance_value > 0.1 &
                                          z_score_ra > 1.96 ~ 'bloom',
                                        abundance_value < 0.1 ~ 'no-bloom' ,
                                        z_score_ra < 1.96 ~ 'no-bloom')) |>
  dplyr::select(date, sample_id, asv_num, fraction, abundance_value, bloom_event, family, z_score_ra) 

anom_rel_abund <- anom_rel_abund_3 |>
  bind_rows(anom_rel_abund_02) |>
  arrange(-abundance_value) |>
  dplyr::mutate(abundance_value = abundance_value*100)

## little function to compute residuals for all ASVs in PA or FL
bloo_02_filt <- bloo_02 |>
  dplyr::filter(!value %in% c('asv2', 'asv3', 'asv5', 'asv8')) 

summary_types_of_blooms_02 <-
  summary_types_of_blooms |>
  dplyr::filter(fraction == '0.2')

bloo_02_filt <-  bloo_02_filt |>
  left_join(summary_types_of_blooms_02, by = c('value' = 'asv_num')) |>
  arrange(-occurrence_perc)

create_plot_02 <- function(asv_num) {
  plot <- plot_residuals(
    data_anom_abund = anom_rel_abund,
    data_mean_abund_asv = mean_abund_f_asv,
    asv_num = {{asv_num}},
    community_fraction = '0.2'
  )
}

plots_list_02 <- purrr::map(bloo_02_filt$value, create_plot_02)

residuals_02 <- gridExtra::grid.arrange(grobs = plots_list_02, ncol = 3) 

ggsave("results/figures/residuals_02_plot_mean.pdf", 
       plot = residuals_02, width = 188, height = 220, units = "mm")

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

bloo_3 <-  bloo_3 |>
  left_join(summary_types_of_blooms_3, by = c('value' = 'asv_num')) |>
  arrange(-occurrence_perc)

bloo_3_1 <-  bloo_3 |>
  as_tibble() |>
  slice_head(n= 22)

bloo_3_2 <-  bloo_3 |>
  as_tibble() |>
  slice_tail(n= 23)

#i divide the PA community because there's too much ASVs for one plot
plots_list_3_1 <- map(bloo_3_1$value, create_plot_3)

residuals_3 <- grid.arrange(grobs = plots_list_3_1, ncol = 3) 

# Assuming you have already created the plot and stored it in `residuals_3`

# Save the plot to a file
ggsave("results/figures/residuals_3.1_plot_mean.pdf", plot = residuals_3, width = 188, height = 220, units = "mm")

#i divide the PA community because there's too much ASVs for one plot

plots_list_3_2 <- map(unique(bloo_3_2$value), create_plot_3)

residuals_3 <- grid.arrange(grobs = plots_list_3_2, ncol = 3)

ggsave("results/figures/residuals_3.2_plot_mean.pdf", plot = residuals_3, width = 188, height = 220, units = "mm")

### I would like to label the super blooming events (those that are outside the sd)----


## I want to check ASVs that are present in PA and FL present similar behavior in both fractions----
### for those that have blooming behavior in both size fractions----
bloo_both_fract <- asv_anom_3_tb |>
  bind_rows(asv_anom_02_tb) |>
  filter(duplicated(value)) |> #6
  left_join(tax_bbmo_10y_new, by = c('value' = 'asv_num')) |>
  dplyr::select(-seq)

asv_tab_10y_02_pseudo_rclr_bloo |>
  bind_rows(asv_tab_10y_3_pseudo_rclr_bloo) |>
  dplyr::filter(asv_num %in% bloo_both_fract$value) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(fraction, abundance_value, decimal_date, asv_num) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  pivot_wider(id_cols = c(asv_num, decimal_date, order_f), values_from = abundance_value, names_from = fraction) |>
  ggplot(aes(as.numeric(`0.2` ), as.numeric(`3`)))+
  geom_point(aes(color = order_f))+
  labs(x = 'Free living (0.2-3um)', y = 'Particle attached (3-20um)', color = 'Order')+
  geom_smooth(method = 'loess', se = F, span = 1, aes(group = asv_num, color = order_f))+
  scale_color_manual(values = palette_order_assigned_bloo)+
  theme_bw()+
  theme(panel.grid = element_blank())

asv_tab_10y_02_pseudo_rclr_bloo |>
  bind_rows(asv_tab_10y_3_pseudo_rclr_bloo) |>
  dplyr::filter(asv_num %in% bloo_both_fract$value) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::select(fraction, abundance_value, decimal_date, asv_num) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  pivot_wider(id_cols = c(asv_num, decimal_date, order_f, family_f), values_from = abundance_value, names_from = fraction) |>
  ggplot(aes(as.numeric(`0.2` ), as.numeric(`3`), group = order_f))+
  geom_point(aes(color = order_f))+
  labs(x = 'Free living (0.2-3um)', y = 'Particle attached (3-20um)', color = 'Order')+
  scale_color_manual(values = palette_order_assigned_bloo)+
  geom_smooth(method = 'loess', se = F, span = 4, aes(group = asv_num, color = order_f))+
  facet_wrap(vars(asv_num), scales = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

### for those that are present in both size fractions but are not considered bloomers in both-----
 asv_tab_10y_02_pseudo_rclr |>
   bind_rows( asv_tab_10y_3_pseudo_rclr) |>
   dplyr::filter(asv_num %in% bloo_taxonomy$asv_num_f) |>
  # dplyr::filter(abundance_type == 'rclr') |>
   dplyr::select(fraction, relative_abundance, decimal_date, asv_num) |>
   left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
   pivot_wider(id_cols = c(asv_num, decimal_date, order_f, family_f), values_from = relative_abundance, names_from = fraction) |>
   ggplot(aes(as.numeric(`0.2` ), as.numeric(`3`), group = order_f))+
   geom_point(aes(color = order_f))+
   labs(x = 'Free living (0.2-3um)', y = 'Particle attached (3-20um)', color = 'Order')+
   scale_color_manual(values = palette_order_assigned_bloo)+
   facet_wrap(vars(asv_num), scales = 'free')+
   geom_smooth(method = 'loess', se = F, span = 1, aes(group = asv_num, color = order_f))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

 asv_tab_10y_02_pseudo_rclr |>
   bind_rows( asv_tab_10y_3_pseudo_rclr) |>
   dplyr::filter(asv_num %in% bloo_taxonomy$asv_num_f) |>
   # dplyr::filter(abundance_type == 'rclr') |>
   dplyr::select(fraction, rclr, decimal_date, asv_num) |>
   left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
   dplyr::filter(!class_f %in% c('Phycisphaerae', 'Bdellovibrionia')) |>
   pivot_wider(id_cols = c(asv_num, decimal_date, order_f, family_f, class_f), values_from = rclr, names_from = fraction) |>
   ggplot(aes(as.numeric(`0.2` ), as.numeric(`3`), group = order_f))+
   geom_point(aes(color = order_f))+
   labs(x = 'Free living (0.2-3um)', y = 'Particle attached (3-20um)', color = 'Order')+
   scale_color_manual(values = palette_order_assigned_bloo)+
   geom_smooth(method = 'loess', se = F, span = 4, aes(group = asv_num, color = order_f))+
   facet_wrap(vars(order_f))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')
 
 asv_tab_10y_02_pseudo_rclr |>
   bind_rows( asv_tab_10y_3_pseudo_rclr) |>
   dplyr::filter(asv_num %in% bloo_taxonomy$asv_num_f) |>
   dplyr::select(fraction, rclr, decimal_date, asv_num, date) |>
   left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
   ggplot(aes(date, rclr))+
   geom_point(aes(color = fraction), alpha = 0.8)+
   labs(x = 'Time', y = 'rCLR', color = 'Fraction')+
   scale_color_manual(values = palette_fraction)+
   #geom_smooth(method = 'loess', se = F, span = 4, aes(group = asv_num, color = order_f))+
   geom_line(aes(group = fraction, color = fraction), alpha = 0.7)+
   facet_wrap(vars(asv_num))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')
 
 
 
## Types of blooms and their taxonomy------
asv_tab_all_bloo_z_tax_summary_all <- read_csv( 'asv_tab_all_bloo_z_tax_summary.csv') |>
   dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
   dplyr::mutate(occurrence_category = ifelse(occurrence_perc > 2/3, 'broad',
                                              ifelse(occurrence_perc < 1/3, 'narrow',
                                                     'intermediate')))
 
 asv_tab_all_bloo_z_tax_summary_all |>
   colnames()
 
asv_tab_all_bloo_z_tax_summary_all |>
  group_by(fraction, recurrency, asv_num, family_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, recurrency, family_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(perc = n/total) |>
  ggplot(aes(as.factor(fraction), perc))+
  geom_col(aes(fill = family_f))+
  facet_wrap(vars(recurrency))+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax_summary_all |>
  group_by(fraction, recurrency, asv_num, order_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, recurrency, order_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(perc = n/total) |>
  group_by(recurrency, fraction) |>
  slice_max(order_by = perc, n = 1)

asv_tab_all_bloo_z_tax_summary_all |>
  group_by(fraction, recurrency, asv_num, order_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, recurrency, order_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(perc = n/total) |>
  ggplot(aes(as.factor(fraction), perc))+
  geom_col(aes(fill = order_f))+
  facet_wrap(vars(recurrency))+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax_summary_all |>
  group_by(fraction, frequency, asv_num, family_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, frequency, family_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(perc = n/total) |>
  ggplot(aes(as.factor(fraction), perc))+
  geom_col(aes(fill = family_f))+
  facet_wrap(vars(frequency))+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  theme_bw()

# asv_tab_all_bloo_z_tax_summary_all |>
#   group_by(fraction, type_of_bloom, asv_num, family_f) |>
#   dplyr::reframe(n = n()) |>
#   group_by(fraction, type_of_bloom, family_f) |>
#   dplyr::reframe(n = n()) |>
#   group_by(fraction) |>
#   dplyr::mutate(total = sum(n)) |>
#   dplyr::mutate(perc = n/total) |>
#   ggplot(aes(as.factor(fraction), perc))+
#   geom_col(aes(fill = family_f))+
#   facet_wrap(vars(type_of_bloom))+
#   scale_y_continuous(labels = percent_format())+
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   theme_bw()
 
asv_tab_all_bloo_z_tax_summary_all |>
  group_by(fraction, occurrence_category, asv_num, family_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, occurrence_category, family_f) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(perc = n/total) |>
  ggplot(aes(as.factor(fraction), perc))+
  geom_col(aes(fill = family_f))+
  facet_wrap(vars(occurrence_category))+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  group_by(fraction, frequency, occurrence_category, asv_num) |>
  dplyr::reframe(total = n()) |>
  group_by(fraction, frequency, occurrence_category) |>
  dplyr::reframe(total = n())  

##number of blooms on the timeseries ----
asv_tab_all_bloo_z_tax_summary_all |>
  colnames()

asv_tab_all_bloo_z_tax_summary_all |>
  str()

asv_tab_all_bloo_z_tax_summary_all |>
  distinct(asv_num) |>
  as_vector()

total_blooming_events <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::mutate(z_score_ra = as.numeric(z_score_ra), 
                abundance_value = as.numeric(abundance_value)) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction) |>
  dplyr::reframe(total_blooming_events = n())

asv_tab_all_bloo_z_tax_summary_all |>
  colnames()

total_blooming_events <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::mutate(z_score_ra = as.numeric(z_score_ra), 
                abundance_value = as.numeric(abundance_value)) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, frequency) |>
  distinct(date, fraction) |>
  ungroup() |>
  distinct(date, frequency) |>
  dplyr::reframe(n = n()) |>
  group_by( frequency) |>
  dplyr::reframe(total_blooming_events = n())

asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(occurrence_category == 'broad') |>
  dplyr::select(occurrence_perc) |>
  slice_min(occurrence_perc, n = 4)

asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  ggplot(aes(date, abundance_value))+
  geom_col(aes(fill = family_f))+
  geom_smooth(method = 'loess', se = F, color = 'darkgrey')+
  scale_y_sqrt()+
  scale_x_datetime()+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_grid(fraction~frequency, labeller = labs_fraction_rec_freq)+ #, labeller = labs_fraction
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 8, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4), strip.background = element_blank(), 
        legend.text = element_text(size = 3), legend.title = element_text(size = 4), strip.placement = 'outside',
        plot.margin = margin(2,5,0,0),
        legend.key.size = unit(3, 'mm'))

bloom_events <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  dplyr::select(asv_num, date, fraction, recurrency, occurrence_category) |>
  dplyr::mutate(fraction = as_factor(fraction))

unique(asv_tab_all_bloo_z_tax$abundance_type)
bloom_events |>
  str()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
 right_join(bloom_events) |>
  ggplot(aes(date, abundance_value))+
  geom_col(aes(fill = order_f))+
  geom_smooth(method = 'loess', se = F)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_grid(fraction~recurrency~occurrence_category)+ #, labeller = labs_fraction
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 8, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6), strip.background = element_blank(), 
        legend.text = element_text(size = 3), legend.title = element_text(size = 4), strip.placement = 'outside',
        plot.margin = margin(2,5,0,0),
        legend.key.size = unit(3, 'mm'))

## We had 48 blooming events in FL and 62 blooming events on PA
total_blooming_events_y <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, year) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, year) |>
  dplyr::reframe(total_blooming_events = n())

blooming_events_year <- total_blooming_events_y |>
  ggplot(aes(year, total_blooming_events, group = as.factor(fraction)))+
  geom_line(aes(color = as.factor(fraction), linetype = as.factor(fraction)))+
  scale_linetype_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  labs(color = 'Fraction', x = 'Year', y ='Nº of blooming events')+
  guides(linetype = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 8, margin = margin(2, 2, 4, 2)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4), strip.background = element_blank(), 
        legend.text = element_text(size = 3), legend.title = element_text(size = 4), strip.placement = 'outside',
        plot.margin = margin(2,5,2,2),
        legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'blooming_events_year.pdf', plot = blooming_events_year,
#        path = 'results/figures/',
#        width = 88, height = 88, units = 'mm')

## number of types of blooms per year-----
total_blooming_events_y_resume <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  dplyr::mutate(blooming_summary = paste0(recurrency, '_', frequency, '_', type_of_bloom, '_', occurrence_category)) |>
  group_by(date, fraction, year, blooming_summary) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, year, blooming_summary) |>
  dplyr::reframe(total_blooming_events = n())

total_blooming_events_y_resume_plot <- total_blooming_events_y_resume |>
  ggplot(aes(year, total_blooming_events, group = as.factor(fraction)))+
  geom_line(aes(color = as.factor(fraction), linetype = as.factor(fraction)))+
  scale_linetype_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  facet_wrap(vars(blooming_summary))+
  labs(color = 'Fraction', x = 'Year', y ='Nº of blooming events')+
  guides(linetype = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 8, margin = margin(2, 2, 4, 2)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4), strip.background = element_blank(), 
        legend.text = element_text(size = 3), legend.title = element_text(size = 4), strip.placement = 'outside',
        plot.margin = margin(2,5,2,2),
        legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'blooming_events_year_resume.pdf', plot = total_blooming_events_y_resume_plot,
#        path = 'results/figures/',
#        width = 188, height = 110, units = 'mm')

## percentage and category of blooming events -----

###DUBTE COMPTEM UN BLOOM EVENT DE 2 ASV COM A UN UNIC BLOOM EVENT O COM A 2? Com a 1-----

## hi ha error al codi pq no sumen el mateix
blooming_events_with_duplicates <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, asv_num) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction) |>
  dplyr::reframe(total_blooming_events = n())

total_blooming_events_perc <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  dplyr::mutate(blooming_summary = paste0(recurrency, '_', frequency, '_', type_of_bloom, '_', occurrence_category)) |>
  group_by(date, fraction, blooming_summary) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, blooming_summary) |>
  dplyr::reframe(total_blooming_events_category = n()) |>
  group_by(fraction) |>
  reframe(sum = sum(total_blooming_events_category))
  left_join(blooming_events_with_duplicates) |>
  dplyr::mutate(perc_category = total_blooming_events_category/total_blooming_events)

group_by(fraction) |>
  reframe(n = sum(n))

total_blooming_events_perc |>
  ggplot(aes(as.factor(fraction), perc_category))+
  geom_col(aes(fill = blooming_summary))+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 8, margin = margin(2, 2, 4, 2)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4), strip.background = element_blank(), 
        legend.text = element_text(size = 8), legend.title = element_text(size = 4), strip.placement = 'outside',
        plot.margin = margin(2,5,2,2),
        legend.key.size = unit(3, 'mm'))

## Boxplot at different frequencies separate stochastic than seasonal blooming events-----
### in this case I count a blooming event unique even though there is more than one taxa involved

total_blooming_events_y <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, year, frequency) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, year, frequency) |>
  dplyr::reframe(total_blooming_events = n()) |>
  dplyr::mutate(time = 'year') |>
  rename(variable = 'year') |>
  dplyr::mutate(variable = as.character(variable))

asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, year) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, year,) |>
  dplyr::reframe(total_blooming_events = n()) |>
  dplyr::mutate(time = 'year') |>
  rename(variable = 'year') |>
  dplyr::mutate(variable = as.character(variable)) |>
  group_by( fraction) |>
  reframe(mean_blooming_events = mean(total_blooming_events),
          sd_blooming_events = sd(total_blooming_events))

total_blooming_events_y |>
  group_by(fraction, recurrency) |>
  reframe(total_timseries = sum(total_blooming_events))
  
total_blooming_events_s <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, season, frequency) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, season, frequency) |>
  dplyr::reframe(total_blooming_events = n()) |>
  dplyr::mutate(time = 'season') |>
  rename(variable = 'season') |>
  dplyr::mutate(variable = as.character(variable))

total_blooming_events_m <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, month, frequency) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, month, frequency) |>
  dplyr::reframe(total_blooming_events = n()) |>
  dplyr::mutate(time = 'month') |>
  rename(variable = 'month')|>
  dplyr::mutate(variable = as.character(variable))

total_blooming_events_time <- total_blooming_events_y |>
  bind_rows(total_blooming_events_s) |>
  bind_rows(total_blooming_events_m)

total_blooming_events_time$variable <- total_blooming_events_time$variable |>
  factor(
    levels = c('winter', 'spring', 'summer', 'autumn',
               '2004', '2005', '2006', '2007', '2008', '2009',
               '2010', '2011', '2012', '2013', 
               '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
               '11', '12'), 
    labels = c( winter = "Winter", spring = "Spring", summer = "Summer", autumn = "Autumn",
                "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013",
                "1" = "January", "2" = "February", "3" = "March", "4" = "April", "5" = "May", "6" = "June",
                "7" = "July", "8" = "August", "9" = "September", "10" = "October", "11" = "November", "12" = "December"
    )
  )

labs_fraction_frequency <-  as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                                          '3' = 'Particle attached (3-20 um)',
                                          seasonal = 'Seasonal',
                                          stochastic = 'Stochastic'))

labs_time <-  as_labeller(c(season = 'Season', month = 'Month', year = 'Year'))

palette_seasons_years_months <- c("Winter" = "#002562", "Spring" = "#519741", "Summer" = "#ffb900","Autumn" =  "#96220a",
                                  '2004' = "#efd9ce",
                                  '2005' = "#d7b7c9",
                                  '2006' =  "#be95c4",
                                  '2007' =  "#af8ec2",
                                  '2008' =  "#9f86c0",
                                  '2009' =  "#7f6da7",
                                  '2010' = "#5e548e",
                                  '2011' = "#231942",
                                  '2012' ="#20173c",
                                  '2013' = "#1d1537", 
                                  "January" = '#002182',
                                  "February" = '#779FFF',
                                  "March" = '#CFD5FF',
                                  "April" = '#003F00',
                                  "May" = '#629956',
                                  "June" = '#D6F7CC',
                                  "July" = '#F59F00',
                                  "August" = '#FEED00',
                                  "September" = '#FFFF8E',
                                  "October" = '#922930',
                                  "November" = '#A60025',
                                  "December" = '#E39793')

total_blooming_events_time |>
  ggplot(aes(time, total_blooming_events))+
  geom_point(aes(group = frequency, color = variable))+
  #geom_boxplot( alpha = 0.1)+
  geom_violin(alpha = 0.1)+
  facet_wrap(fraction~frequency, labeller = labs_fraction_frequency)+
  scale_x_discrete(labels = labs_time)+
  labs(y = 'Total blooming events', x = '', color = 'Variable')+
  scale_y_continuous(limits = c(0, 15), expand = c(0,0))+
  scale_color_manual(values = palette_seasons_years_months)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'))

total_blooming_events_time_rel_freq <- total_blooming_events_time |>
  dplyr::mutate(total_blooming_events_rel = case_when(time == 'year' ~ total_blooming_events/10,
                                                      time == 'season' ~ total_blooming_events/(4*10),
                                                      time == 'month' ~ total_blooming_events/12)) |>
  ggplot(aes(time, total_blooming_events_rel))+
  geom_point(aes(group = frequency, color = variable), position = position_jitter(width = 0.25))+
  #geom_boxplot( alpha = 0.1)+
  geom_violin(alpha = 0.1)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  facet_wrap(fraction~frequency, labeller = labs_fraction_frequency)+
  scale_x_discrete(labels = labs_time)+
  labs(y = 'Total blooming events/\n Time unit', x = '', color = '')+
  scale_y_continuous( )+
  scale_color_manual(values = palette_seasons_years_months)+

  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank())

# ggsave(filename = 'total_blooming_events_time_rel_freq.pdf', plot = total_blooming_events_time_rel_freq,
#        path = 'results/figures/',
#        width = 188, height = 130, units = 'mm')

## recurrent non recurrent instead of frequency-----
total_blooming_events_y <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, year, recurrency) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, year, recurrency) |>
  dplyr::reframe(total_blooming_events = n()) |>
  dplyr::mutate(time = 'year') |>
  rename(variable = 'year') |>
  dplyr::mutate(variable = as.character(variable))

total_blooming_events_s <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, season, recurrency) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, season, recurrency) |>
  dplyr::reframe(total_blooming_events = n()) |>
  dplyr::mutate(time = 'season') |>
  rename(variable = 'season') |>
  dplyr::mutate(variable = as.character(variable))

total_blooming_events_m <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(abundance_value > 0.1 &
                  z_score_ra >= 1.96) |>
  group_by(date, fraction, month, recurrency) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, month, recurrency) |>
  dplyr::reframe(total_blooming_events = n()) |>
  dplyr::mutate(time = 'month') |>
  rename(variable = 'month')|>
  dplyr::mutate(variable = as.character(variable))
total_blooming_events_time <- total_blooming_events_y |>
  bind_rows(total_blooming_events_s) |>
  bind_rows(total_blooming_events_m)

labs_fraction_rec <-  as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                                          '3' = 'Particle attached (3-20 um)',
                                          no = 'Non-recurrent',
                                          yes = 'Recurrent'))

total_blooming_events_time$variable <- total_blooming_events_time$variable |>
  factor(
    levels = c('winter', 'spring', 'summer', 'autumn',
               '2004', '2005', '2006', '2007', '2008', '2009',
               '2010', '2011', '2012', '2013', 
               '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
               '11', '12'), 
    labels = c( winter = "Winter", spring = "Spring", summer = "Summer", autumn = "Autumn",
                "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013",
                "1" = "January", "2" = "February", "3" = "March", "4" = "April", "5" = "May", "6" = "June",
                "7" = "July", "8" = "August", "9" = "September", "10" = "October", "11" = "November", "12" = "December"
    )
  )


total_blooming_events_time |>
  ggplot(aes(time, total_blooming_events))+
  geom_point(aes(group = recurrency, color = variable))+
  #geom_boxplot( alpha = 0.1)+
  geom_violin(alpha = 0.1)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  facet_wrap(fraction~recurrency, labeller = labs_fraction_rec)+
  scale_x_discrete(labels = labs_time)+
  labs(y = 'Total blooming events', x = '', color = 'Variable')+
  scale_y_continuous(limits = c(0, 15), expand = c(0,0))+
  scale_color_manual(values = palette_seasons_years_months)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'))

total_blooming_events_time |>
  ggplot(aes(time, total_blooming_events))+
  geom_point(aes(group = recurrency, color = variable))+
  #geom_boxplot( alpha = 0.1)+
  geom_violin(alpha = 0.1)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  facet_wrap(vars(fraction), labeller = labs_fraction_rec)+
  scale_x_discrete(labels = labs_time)+
  labs(y = 'Total blooming events', x = '', color = 'Variable')+
  scale_y_continuous(limits = c(0, 15), expand = c(0,0))+
  scale_color_manual(values = palette_seasons_years_months)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'))

total_blooming_events_time_rel_rec <- total_blooming_events_time |>
  dplyr::mutate(total_blooming_events_rel = case_when(time == 'year' ~ total_blooming_events/10,
                                                      time == 'season' ~ total_blooming_events/(4*10),
                                                      time == 'month' ~ total_blooming_events/12)) |>
  ggplot(aes(time, total_blooming_events_rel))+
  geom_point(aes(group = recurrency, color = variable), position = position_jitter(width = 0.25))+
  #geom_boxplot( alpha = 0.1)+
  geom_violin(alpha = 0.1)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  facet_wrap(fraction~recurrency, labeller = labs_fraction_rec)+
  scale_x_discrete(labels = labs_time)+
  labs(y = 'Total blooming events/\n Time unit', x = '', color = '')+
  scale_y_continuous( )+
  scale_color_manual(values = palette_seasons_years_months)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid.minor = element_blank())

ggsave(filename = 'total_blooming_events_time_rel_rec.pdf', plot = total_blooming_events_time_rel_rec,
       path = 'results/figures/',
       width = 188, height = 130, units = 'mm')
  
## Number of reads that "bloomers" represent from the whole community ----

#### could i divide the number of reads / sampling effort?
bloom_events_ed <- bloom_events |>
  dplyr::mutate(fraction = as.double(fraction))

asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::mutate(fraction = '0.2') |>
  dplyr::mutate(fraction = as.double(fraction)) |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  left_join(bloo_all_types_summary_tb_tax, by = c('asv_num', 'fraction')) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  left_join(bloom_events_ed) |>
  dplyr::group_by(sample_id, frequency)|> 
  dplyr::reframe(total_reads = sum(reads)) |>
  left_join(m_02) |>
  ggplot(aes(date, total_reads))+
  geom_line(aes(group = frequency))+
  geom_point()+
  facet_wrap(vars(frequency))+
  labs(x = 'Date', y = 'Total reads')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_rect('transparent'))


bloo_reads_02 <- asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::reframe(total_reads = sum(reads))

asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::reframe(total_reads = sum(reads)) |>
  bind_rows(bloo_reads_02)

perc_bloo_reads_02 <- 891501/3930422*100

bloo_reads3 <- asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::reframe(total_reads = sum(reads))

asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  dplyr::reframe(total_reads = sum(reads)) |>
  bind_rows(bloo_reads3)

perc_bloo_reads3 <- 1691591/4607899*100

asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  dplyr::mutate(fraction = '3') |>
  dplyr::mutate(fraction = as.double(fraction)) |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  left_join(bloo_all_types_summary_tb_tax, by = c('asv_num', 'fraction')) |>
  left_join(bloom_events_ed) |>
  dplyr::group_by(sample_id, frequency)|> 
  dplyr::reframe(total_reads = sum(reads)) |>
  left_join(m_3) |>
  ggplot(aes(date, total_reads))+
  geom_line(aes(group = frequency))+
  geom_point()+
  facet_wrap(vars(frequency))+
  labs(x = 'Date', y = 'Total reads')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_rect('transparent'))

asv_tab_bbmo_10y_l |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  dplyr::group_by(sample_id, fraction, date)|> 
  dplyr::reframe(total_reads = sum(reads)) |>
  ggplot(aes(date, total_reads))+
  geom_line()+
  geom_point()+
  facet_wrap(vars(fraction))+
  labs(x = 'Date', y = 'Total reads')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_rect('transparent'))

##when we have a blooming event do we have an increase in sampling effort??  

##bacterial production
m_bbmo_10y |>
  colnames()

bloom_events |>
  colnames()

m_bbmo_10y |>
  left_join(bloom_events) |>
  ggplot(aes(date, BP_FC1.55))+
  geom_line()+
  geom_vline(xintercept = bloom_events$date, color = 'black')+
  labs(x = 'Date', y = 'Total reads')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_rect('transparent'))

asv_tab_bbmo_10y_l |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  dplyr::group_by(sample_id, fraction, date,  BP_FC1.55)|> 
  dplyr::reframe(total_reads = sum(reads)) |>
  ggplot(aes(total_reads,  BP_FC1.55))+
  geom_point()

asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  dplyr::mutate(fraction = '3') |>
  dplyr::mutate(fraction = as.double(fraction)) |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  left_join(bloo_all_types_summary_tb_tax, by = c('asv_num', 'fraction')) |>
  left_join(bloom_events_ed) |>
  dplyr::group_by(sample_id, frequency)|> 
  dplyr::reframe(total_reads = sum(reads)) |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  ggplot(aes(total_reads,  BP_FC1.55))+
  geom_point()+
  facet_wrap(vars(frequency))

asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::mutate(fraction = '0.2') |>
  dplyr::mutate(fraction = as.double(fraction)) |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  left_join(bloo_all_types_summary_tb_tax, by = c('asv_num', 'fraction')) |>
  left_join(bloom_events_ed) |>
  dplyr::group_by(sample_id, frequency)|> 
  dplyr::reframe(total_reads = sum(reads)) |>
  left_join(m_bbmo_10y, by = 'sample_id') |>
  ggplot(aes(total_reads,  BP_FC1.55))+
  geom_point()+
  facet_wrap(vars(frequency))


m_02 <- m_02 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

bloom_events |>
  left_join(m_02, by = c('date')) |>
  ggplot(aes(date, BP_FC1.55))+
  geom_line()




