# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                   data analysis pipeline                    ++++++++++++++++++++++
# +++++++++++++++++++++++                   BBMO Dateseries 10-Y data                 ++++++++++++++++++++++
# +++++++++++++++++++++++                   metabarcoding & metagenomic               ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## packages ----
library(cowplot)
library(tidyverse)
library(ggpmisc)

# ------------------ ## --------- Blooming events and their functional implications ------------------ ## ---
## upload palettes ----

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




# ------------------ ## --------- Community level ----------------------- ## ----------------------- ## -----
## upload data correlations datasets genes and community
bray_unifrac_eucl_tb_02 <- read.csv('data/bray_unifrac_eucl_tb_02.csv')
bray_unifrac_eucl_tb <- read.csv('data/bray_unifrac_eucl_tb.csv')
corr_bray_genes_community_tb <- read.csv('data/corr_bray_genes_community_tb_scg.csv')

## Relationship between Bray Curtis Genes and Bray Curtis of the Community ----
#### first correlation: Bray-Curtis genes vs Bray-Curtis community, and Bray-Curtis genes vs. weighted UNIFRAC distance ----
corr_bray_02_plot <- bray_unifrac_eucl_tb_02 |>
  ggplot(aes(bray_curtis_community, genes))+
  scale_x_continuous(limits = c(0.05, 0.85))+
  scale_y_continuous(limits = c(0.05, 0.85))+
  geom_point(color = "#00808F", alpha = 0.8)+
  stat_cor(aes(
    label =   paste(..p.label..)), label.x = 0.35,
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#00808F"#,
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.35, label.y = 0.2, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#00808F")+
  geom_smooth(method = 'lm', color = "#00808F", fill = "#00808F")+
  
  labs(x = 'Bray-Curtis Community-Based', y = 'Bray-Curtis Genes-Based')+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'bottom', axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10), strip.background = element_blank(), 
    legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
    strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_02_plot

## I add KOs and genes in the same correlation plot. KOs will inform us if there is a change in the functionality -----
corr_bray_02_kos_plot <- corr_bray_genes_community_tb |>
  ggplot(aes(bray_curtis_community, KEGG_sgc))+
  scale_x_continuous(limits = c(0.05, 0.85))+
  scale_y_continuous(limits = c(0.05, 0.25))+
  geom_point( 
    alpha = 0.8,  color = "#00808F")+
  stat_cor(aes( 
    label =   paste(..p.label..)), label.x = 0.6,
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#00808F"#,
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.6, label.y = 0.2, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#00808F")+
  geom_smooth(method = 'lm', color = "#00808F", fill = "#00808F")+
  labs(x = 'Bray-Curtis Community-Based', y = 'Bray-Curtis KOs-Based')+
  theme_bw()+
  theme( 
    panel.grid.minor = element_blank(),
    legend.position = 'bottom', axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10), strip.background = element_blank(), 
    legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
    strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_02_kos_plot

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

# ------------------ ## --------- Bloomers level ----------------------- ## ----------------------- ## -----
## upload metabarcoding data ----

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


asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') |> ##using dada2 classifier assign tax with silva 138.1 and correctly identifying bloomers
  as_tibble() |>
  dplyr::select(-X)

## reorder taxonomy as factors
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

bloom_event <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(bloom_event = case_when(abundance_type == 'relative_abundance' &
                                          abundance_value >= 0.1 &
                                          z_score_ra > 1.96 ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  filter(bloom_event != 0) |>
  dplyr::select(date, asv_num, bloom_event, fraction) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ungroup() |>
  distinct()



## -------- Are there consistent trends between heterotrophic blooms and genetic content from metagenomes? -
### Heterotrophic bloomers during 2009-2014 FL fraction metabarcoding -----
heterotrophic_blooms <- bloom_event |>
  dplyr::filter(fraction == '0.2' & bloom_event == 'bloom') |>
  dplyr::filter(asv_num %in% c('asv11', 'asv27', 'asv58')) |>
  dplyr::filter(str_detect(date, '2009|2010|2011|2012|2013')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

heterotrophic_blooms <- heterotrophic_blooms |>
  dplyr::mutate(tax = case_when(asv_num == 'asv11' ~ 'Glaciecola',
                                asv_num == 'asv58'  ~ 'NS4 Marine Group',
                                asv_num == 'asv27'  ~ 'Amylibacter'))

asv_tab_0914 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv27', 'asv58'))

blooms_0914_plot <- asv_tab_0914 |>
  group_by(asv_num, order_f, genus, date) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, max_abund))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_area(aes(y = abundance_value, group = asv_num, fill = order_f), position = 'identity', alpha = 0.2)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
  )+
  geom_line(aes(y = abundance_value, group = asv_num, color = order_f), linewidth = 1)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  labs(x = 'Date', y = 'Relative Abundance', fill = '')+
  theme_bw()+
  theme(legend.position = 'none', strip.background = element_rect(fill = NA),
        panel.grid.minor = element_blank())

blooms_0914_plot

### Total Genes related to bloomers timeseries ---------
# from gene annotation in MARBITS I select those genes that have the closest tax to my bloomers 
genes_tax_tb <- read.csv('data/genes_sergio/genes_interest_bloo.tbl', sep = '\t')

genes_tax_tb |>
  colnames() <- c('gene', 'counts', 'tax_resolution', 'tax', 'tax_cat' )

genes_tax_tb_f <- genes_tax_tb  |>
  dplyr::filter(!str_detect(gene, 'SO')) # we remove genes that belong to SOLA Date series.

genes_tax_tb_f |>
  dim() #564 418 

genes_tax_tb_f <- genes_tax_tb_f  |>
  dplyr::select(gene, tax) |>
  dplyr::mutate(tax_bloo_close = case_when(str_detect(tax, 'Amylibacter') ~ 'Amylibacter',
                                           str_detect(tax, 'Glaciecola')  ~ 'Glaciecola',
                                           str_detect(tax, '121220')  ~      'NS4 marine group'))
## I filter the genes table in MARBITS 

#### tables of genes related to bloomers abundances ----------
# upload table of those genes that have a relationship with my bloomers
# genes_bloo_tb <- read.csv('data/genes_sergio/BBMO-GC_250bp_gene.lengthNorm.SingleCopyGeneNorm.counts_bloo_genes_bloo.tbl',
#                           sep = '\t')
genes_bloo_tb <- read.csv('data/genes_sergio/BBMO-GC_250bp_gene.lengthNorm.metaGsizeGbNorm.counts_bloo_genes_bloo.tbl',
                          sep = '\t')

genes_bloo_tb_tax <- genes_bloo_tb |>
  left_join(genes_tax_tb_f, by = c('gene')) |>
  dplyr::filter(!is.na(tax_bloo_close))

m_bbmo_10y_ed_4metag <- m_bbmo_10y |>
  separate(sample_id, into = c('sample_id_sim', 'fraction', 'code'), sep = '_', remove = F) |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::select(sample_id_sim, year, day_of_year, decimal_date, month, day, date, sample_id) |>
  arrange(date) |>
  dplyr::mutate(sample_id_num = row_number())

m_bbmo_10y_ed_4metag_red <- m_bbmo_10y_ed_4metag |>
  dplyr::select(-sample_id)

m_bbmo_10y_ed_4metag_red2 <- m_bbmo_10y_ed_4metag_red |>
  dplyr::mutate(sample_id_ed = paste0('BL', 1:nrow(m_bbmo_10y_ed_4metag)))

genes_bloo_tb_tax_l <- genes_bloo_tb_tax |>
  dplyr::select(tax_bloo_close, starts_with('BL')) |>
  pivot_longer(starts_with('BL'), names_to = 'sample_id') |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_sim')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

genes_bloo_tb_tax_l_sum <- genes_bloo_tb_tax_l |>
  dplyr::select(sample_id_ed, value, tax_bloo_close, date, sample_id ) |>
  dplyr::group_by(tax_bloo_close, date, sample_id) |>
  dplyr::reframe(sum_genes = sum(value))

genes_bloo_tb_tax_l_sum$tax_bloo_close |>
  unique()

plot_kos_counts_bloo <- genes_bloo_tb_tax_l_sum |>
  ggplot(aes(date, sum_genes))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_line(aes(group = tax_bloo_close, color = tax_bloo_close), linewidth = 1)+
  geom_area(aes(fill =  tax_bloo_close), position = 'identity', alpha = 0.1)+
  scale_color_manual(values = c('NS4 marine group' = "#0051BF" , 'Amylibacter' = "#2d373b",  'Glaciecola' = "#f1c510")) +
  scale_fill_manual(values = c('NS4 marine group' = "#0051BF" , 'Amylibacter' = "#2d373b",  'Glaciecola' = "#f1c510")) +
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'TPM')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'none',
        panel.grid.minor = element_blank())

plot_kos_counts_bloo  

# total counts / sample to know the relative abundance of those genes to the total sample
# genes_scg_sum_tb <- read.csv('data/genes_sergio/BBMO-GC_250bp_gene.lengthNorm.SingleCopyGeneNorm.counts_sum.tbl',
#                           sep = '\t', header = F)

genes_scg_sum_tb <- read.csv('data/genes_sergio/BBMO-GC_250bp_gene.lengthNorm.metaGsizeGbNorm.counts_bloo_genes_sum.tbl',
                             sep = '\t', header = F)

genes_scg_sum_tb |>
  colnames() <- c('sum') 

date_code <- genes_bloo_tb |>
  colnames() |>
  as_tibble_col('sample_id')

date_code <- date_code[-1,]

genes_scg_sum_tb <- genes_scg_sum_tb |>
  bind_cols(date_code)

##### calculate the % of genes potentially related with my bloomers -----
# perc_genes_tb <- genes_bloo_tb_tax_l_sum |>
#   left_join(genes_scg_sum_tb) |>
#   dplyr::mutate(perc = sum_genes/sum)
# 
# plot_perc_genes_bloo <- perc_genes_tb |>
#   ggplot(aes(date, perc)+
#   geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
#   geom_line(aes(group = tax_bloo_close, color = tax_bloo_close), linewidth = 1)+
#   geom_area(aes(fill =  tax_bloo_close), position = 'identity', alpha = 0.1)+
#   scale_color_manual(values = c('NS4 marine group' = "#0051BF" , 'Amylibacter' = "#2d373b",  'Glaciecola' = "#f1c510")) +
#   scale_fill_manual(values = c('NS4 marine group' = "#0051BF" , 'Amylibacter' = "#2d373b",  'Glaciecola' = "#f1c510")) +
#   scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
#   labs(x = 'Date', y = 'Relatvie Genes Counts (TPM)')+
#   theme_bw()+
#   theme(axis.text = element_text(size = 4), legend.position = 'none')
# 
# plot_perc_genes_bloo
# 
# perc_genes_tb_tpm <- perc_genes_tb
# 
# perc_genes_tb 

## -------- Are bloomers sustaining some functions (KOs) during bloom events? --------
# compare KOs abundance vs gene related to bloomers tax abundance -----
## We will study two bloom events, January 2011 (Glaciecola) and April 2009 Amylibacter 
## KOs table
kos_tb <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.tbl')
#kos_tb <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.metaGsizeNorm.counts.txt')

kos_tb <- kos_tb[,1:61]

kos_tb_names <- kos_tb %$%
  V1

kos_tb_t <- t(kos_tb)

kos_tb_t <- kos_tb_t[-1,]

kos_tb_t |>
  colnames() <- kos_tb_names

genes_bloo_tb_tax_l |>
  colnames()

## I add to the genes tb the KO they are related to 
kos_genes_relation_tb <- read.table( 'data/genes_sergio/kos_tax_interest_bloo_bbmo.tbl', sep = '\t')

kos_genes_relation_tb |>
  head()

kos_genes_relation_tb |>
  colnames() <- c('gene', 'KO')


genes_bloo_tb <- read.csv('data/genes_sergio/BBMO-GC_250bp_gene.lengthNorm.SingleCopyGeneNorm.counts_bloo_genes_bloo.tbl',
                          sep = '\t')

genes_bloo_tb_tax <- genes_bloo_tb |>
  left_join(genes_tax_tb_f, by = c('gene')) |>
  dplyr::filter(!is.na(tax_bloo_close))

### Bloom Glaciecola 
top25_kos_gla <- genes_bloo_tb_tax |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id') |>
  left_join(kos_genes_relation_tb, by = c('gene')) |>
  dplyr::filter(tax_bloo_close == 'Glaciecola') |>
  dplyr::filter(str_detect(sample_id, 'BL1101')) |>
  group_by(KO, sample_id, tax_bloo_close) |>
  dplyr::reframe(contribution_to_ko = sum(value)) |>
  dplyr::filter(contribution_to_ko != 0) |>
  ungroup() |>
  slice_max(order_by = contribution_to_ko, n = 25) %$%
  KO |>
  as.vector()

genes_bloo_ko_1101_tb <- genes_bloo_tb_tax |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id') |>
  left_join(kos_genes_relation_tb, by = c('gene')) |>
  dplyr::filter(tax_bloo_close == 'Glaciecola') |>
  dplyr::filter(str_detect(sample_id, 'BL1101')) |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_sim')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  group_by(KO, sample_id, tax_bloo_close, date) |>
  dplyr::reframe(contribution_to_ko = sum(value)) |>
  dplyr::filter(contribution_to_ko != 0) |>
  slice_max(order_by = contribution_to_ko, n = 25) |>
  dplyr::select(tax_bloo_close, KO, sample_id, contribution_to_ko, date)

genes_bloo_ko_1101_tb %$%
  contribution_to_ko |>
  range()

### Bloom Amylibacter 
top25_kos_amy <- genes_bloo_tb_tax |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id') |>
  left_join(kos_genes_relation_tb, by = c('gene')) |>
  dplyr::filter(tax_bloo_close == 'Amylibacter') |>
  dplyr::filter(str_detect(sample_id, 'BL0904')) |>
  group_by(KO, sample_id, tax_bloo_close) |>
  dplyr::reframe(contribution_to_ko = sum(value)) |>
  dplyr::filter(contribution_to_ko != 0) |>
  ungroup() |>
  slice_max(order_by = contribution_to_ko, n = 25) %$%
  KO |>
  as.vector()

genes_bloo_ko_0904_tb <- genes_bloo_tb_tax |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id') |>
  left_join(kos_genes_relation_tb, by = c('gene')) |>
  dplyr::filter(tax_bloo_close == 'Amylibacter') |>
  dplyr::filter(str_detect(sample_id, 'BL0904')) |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_sim')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  group_by(KO, sample_id, tax_bloo_close, date) |>
  dplyr::reframe(contribution_to_ko = sum(value)) |>
  dplyr::filter(contribution_to_ko != 0) |>
  slice_max(order_by = contribution_to_ko, n = 25) |>
  dplyr::select(tax_bloo_close, KO, sample_id, contribution_to_ko, date)

genes_bloo_ko_0904_tb %$%
  contribution_to_ko |>
  range()

## upload information about these KOs 
top25_kos_amy |>
  cat(top25_kos_gla)

top_25_kos_names <- read.csv2('data/genes_sergio/top25_kos_modules_scg.csv')

top_25_kos_names |>
  colnames()

top_25_kos_names |>
  distinct(module)

## upload KEGGs counts SCG
kos_tb_t <- kos_tb_t |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(kos_tb_t))) |>
  dplyr::filter(str_detect(annot, 'BL'))

kos_tb_t_l_bl1101 <- kos_tb_t |>
  dplyr::select(-sum_not_annotated, -sample_id) |>
  pivot_longer(starts_with('K'), values_to = 'scg_kos_counts', names_to = 'asv_num') |>
  dplyr::filter(asv_num %in% top25_kos_gla) |>
  dplyr::filter(str_detect(annot, 'BL1101'))

kos_tb_t_l_bl0904 <- kos_tb_t |>
  dplyr::select(-sum_not_annotated, -sample_id) |>
  pivot_longer(starts_with('K'), values_to = 'scg_kos_counts', names_to = 'asv_num') |>
  dplyr::filter(asv_num %in% top25_kos_amy) |>
  dplyr::filter(str_detect(annot, 'BL0904'))

## Are this taxa sustaining functions in the ecosystem during bloom events? 
### calculate relative contribution of my bloomers tax to each specific KO ----
relative_contribution_gla <- genes_bloo_ko_1101_tb |>
  left_join(kos_tb_t_l_bl1101, by = c('KO' = 'asv_num', 'sample_id' = 'annot')) |>
  dplyr::mutate(relative_contribution = as.numeric(contribution_to_ko)/as.numeric(scg_kos_counts)) |>
  left_join(top_25_kos_names, by = c('KO' = 'KO'))

relative_contribution_gla_plot <- relative_contribution_gla |>
  dplyr::mutate(ko_gene = paste0(KO, ' ', gene, module)) |>
  dplyr::mutate(ko_gene = as.factor(ko_gene)) |>
  dplyr::mutate(ko_gene = fct_reorder(ko_gene, relative_contribution, .desc = F)) |>  # Reorder factor
  ggplot(aes(relative_contribution, fct_infreq(ko_gene)))+
  geom_col(aes(), fill = "#f1c510")+
  labs(x = 'Relative contribution', y = '')+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 4),
        axis.title.x = element_text(size = 5),
        axis.ticks.length = unit(0.2, "mm"))

relative_contribution_gla_plot

# counts_contribution_gla_plot
relative_contribution_amy <- genes_bloo_ko_0904_tb |>
  left_join(kos_tb_t_l_bl0904, by = c('KO' = 'asv_num', 'sample_id' = 'annot')) |>
  dplyr::mutate(relative_contribution = as.numeric(contribution_to_ko)/as.numeric(scg_kos_counts)) |>
  left_join(top_25_kos_names, by = c('KO' = 'KO'))

relative_contribution_amy_plot <- relative_contribution_amy |>
  dplyr::mutate(ko_gene = paste0(KO, ' ', gene, module)) |>
  dplyr::mutate(ko_gene = as.factor(ko_gene)) |>
  dplyr::mutate(ko_gene = fct_reorder(ko_gene, relative_contribution, .desc = F)) |>  # Reorder factor
  ggplot(aes(relative_contribution, fct_infreq(ko_gene)))+
  geom_col(aes(), fill = "#2d373b")+
  labs(x = 'Relative contribution', y = '')+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_text(size = 5),
        axis.ticks.length = unit(0.2, "mm"))

relative_contribution_amy_plot

plot_grid(relative_contribution_gla_plot, relative_contribution_amy_plot, 
          nrow = 1,
          rel_widths = c(1, 1))

### Composition plot genes and heterotrophic bloom events ----
timeseries_bloom0914_fl_plot <- plot_grid(
  blooms_0914_plot,
  plot_kos_counts_bloo,
  rel_heights = c(1,1, 1, 2),
  labels = c('C', 'D'),
  ncol = 1)

kos_fl_plot  <- plot_grid(relative_contribution_amy_plot, 
                          relative_contribution_gla_plot,
                          nrow = 1,
                          rel_widths = c(1,  1),
                          labels = c('E', '', 'F'))

kos_genes_plot <- plot_grid(timeseries_bloom0914_fl_plot,
                            kos_fl_plot,
                            rel_heights = c(1.15,1),
                            ncol = 1)

composition_plot <- plot_grid( corr_bray_unifrac_kos_plot,
                               kos_genes_plot,
                               ncol = 2,
                               rel_widths = c(0.75, 1))

composition_plot

ggsave( plot = composition_plot,
        filename = 'genes_kos_bloomers_v3.pdf',
        path = 'results/figures/',
        width = 180, height = 160, units = 'mm')
