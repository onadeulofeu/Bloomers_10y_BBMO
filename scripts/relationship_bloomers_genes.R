# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


##packages
#library(cowplot)
library(tidyverse)
library(gridExtra)

##palettes-----
palette_genes <- c("#ffe355",
                   "#ef8d00",
                   "#ae659b",
                   "#2b347a",
                   "#c23939",
                   "#4cb76a",
                   "#804c90",
                   "#006cb0",
                   "#518535",
                   "#2d373b",
                   "#6a6964",
                   "#57a9a8",
                   "#c9acb8",
                   "#8289ba",
                   "#003000",
                   "#f1c549",
                   "#9f0085",
                   "#bec735",
                   "#c5ebff",
                   "#ca6094",
                   "#000000",
                   "#c7c7c7",
                   "#670000",
                   "#ff9d91")

# Upload data -----
## There are differences in the codes for metagenomic data and metabarcoding data. 
## From now on in this part of the analysis I will numerate my sampples fomr 1 to 60 which are the number of samples that are shared between both datasets.
## edit metadata sample_id column so that it matches the sample_id code that we have in the genes data
m_bbmo_10y |>
  colnames()

m_bbmo_10y_ed_4metag <- m_bbmo_10y |>
    separate(sample_id, into = c('sample_id_sim', 'fraction', 'code'), sep = '_', remove = F) |>
    dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::select(sample_id_sim, year, day_of_year, decimal_date, month, day, date, sample_id) |>
  arrange(date) |>
  dplyr::mutate(sample_id_num = row_number())

#upload data from metagenomes we have the table with counts for the genes related to phoyostynthesis----
m_bbmo_10y_ed_4metag_red <- m_bbmo_10y_ed_4metag |>
  dplyr::select(-sample_id)

genes_photo <- read.csv2('../data/genes_sergio/K02692_genes_detail.csv') |>
  dplyr::select(annot, starts_with('BL')) |> #filter only Blanes samples not SOLA 
  as_tibble() |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'gene_counts') |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  group_by(sample_id)  |>
  arrange(sample_id) |>
  mutate(sample_id_num = cur_group_id()) |>
  left_join(m_bbmo_10y_ed_4metag_red, by = c('sample_id_num')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013'))

genes_photo |>
  dplyr::filter(sample_id != sample_id_sim) |>
  distinct(sample_id, sample_id_sim)

## Differences in the dates from one dataset to the other. Why?
# A tibble: 4 × 2
# Groups:   sample_id [4]
# sample_id sample_id_sim
# <chr>     <chr>        
#   1 BL100413  BL100414     
# 2 BL110704  BL110705     
# 3 BL120518  BL120511     
# 4 BL131204  BL131215  

## explore a bit this type of data----
genes_photo %$%
  gene_counts |>
  range() #0.00000 17.37556

genes_photo <- genes_photo|>
  dplyr::mutate(annot = as.factor(annot))

genes$annot <- factor(genes$annot, levels = fct_infreq(levels(genes$annot)))

genes_photo |>
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 1)) |> #filter by genes_photo that have a minimum count at some point
  ggplot(aes(date, gene_counts))+
  facet_wrap(vars(annot))+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  geom_line(aes(group = annot), color = '#5B5A5A')+ #, color = '#3D3B3B'
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        legend.position = 'bottom')
  
## we would like to detect anomalies in genes_photoand see if they correlate with anomalies in ASVs that are selected as potential bloomers
z_02_genes_photo<- genes_photo |>
  group_by(annot) |>
  arrange(date) |>
  group_by(annot) |>
  dplyr::reframe(anomalies_gene = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, 
                                                cutoff = 1,96, values = gene_counts, 
                                                plotting = FALSE)[c(1,2,3)])

z_scores_02_genes_photo <- genes_photo |>
  group_by(annot) |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  left_join(m_bbmo_10y_ed_4metag, by = c('sample_id_num')) |>
  #arrange(date) |>
  group_by(annot) |>
  dplyr::reframe(z_score_gen = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = gene_counts, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_gen) |>
  group_by(annot) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_bbmo_10y_ed_4metag))) |>
  dplyr::mutate(sample_id_num = as.integer(sample_id_num)) |>
  left_join(m_bbmo_10y_ed_4metag, by = 'sample_id_num') 

## filter by the years that we share after computing the z_scores so that I am able to calculate the last ones (crec que no cal perque els NAs els trobaré al final)

## plot the genes_photo with their anomalies an observe if they have correspondence with the blooming events----
genes_photo |>
  left_join(z_scores_02_genes_photo) %$%
  gene_counts |>
  range()

genes_photo |>
  left_join(z_scores_02_genes_photo) |>
  dplyr::filter(gene_counts > 17)

z_scores_02_genes_photo_f <- z_scores_02_genes_photo |>
  dplyr::mutate(z_score_gen = case_when(is.na(z_score_gen) ~ 0,
                                        z_score_gen == 'NaN' ~ 0,
                                        z_score_gen == Inf ~ 10000,
                                          TRUE ~ z_score_gen)) |>
  dplyr::mutate(z_score_gen = as.numeric(z_score_gen)) |>
  group_by(annot) |>
  dplyr::filter(any(z_score_gen > 1.96, na.rm = TRUE)) |> # only plot those genes_photo that present an anomaly at some point
  dplyr::select(annot, z_score_gen, sample_id_sim, date)

genes_anom <- z_scores_02_genes_photo_f %$%
  annot |>
  unique() 

genes_counts_z_scores_photo <- genes_photo |>
  filter(annot %in% genes_anom) |>
  dplyr::select(sample_id, gene_counts, annot, date) |>
  left_join(z_scores_02_genes_f, by = c('sample_id' = 'sample_id_sim', 'annot')) ##

genes_counts_z_scores_photo_anom <- genes_counts_z_scores_photo |>
  group_by(annot) |>
  dplyr::filter(z_score_gen > 1.96 &
                  any(gene_counts > 5))
  
genes_photo_plot <- genes_counts_z_scores_photo |>  
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 3)) |> #filter by genes_photothat have a minimum count at some point
  dplyr::filter(any(z_score_gen > 5050)) |>
  #dplyr::filter(any(gene_counts > 4)) |> #filter by genes_photothat have a minimum count at some point
  ungroup() |>
  #dplyr::mutate(date.x = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_gen) >= 1.96, '#9F0011', '#080808', missing = '#080808')) |>
  group_by(date.x) |>
  dplyr::mutate(max_counts = sum(gene_counts)) |>
  ggplot(aes(date.x, max_counts))+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0))+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  # geom_point(aes(color = anomaly_color))+
  geom_vline(xintercept = abund_asv_bloom_events_02$date, color = '#8B8989')+
  geom_area(aes(group = annot, fill = annot, y = gene_counts))+
  # geom_point(data = genes_counts_z_scores_anom,
  #   aes(date.x, gene_counts, color =  '#9F0011'), size = 1)+
  scale_fill_manual(values = palette_genes)+
  labs(y = 'Gene counts', x = 'Time (Y)', title = 'K02692, photosyntesis')+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'none')

genes_photo_plot
  
## potential bloomers dataset for the same dates as we have metagenomes 2009-2013------
##blooming events 0.2 between 2009-2013 comparison with Metagenomes
asv_tab_all_bloo_z_tax |>
  colnames()

abund_asv_bloom_events_02 |>
  distinct(asv_num)

abund_asv_bloom_events_02 <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  as.numeric(z_score_ra) >= 1.96 &
                  as.numeric(abundance_value) >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value, asv_num) |>
  dplyr::filter(year %in% c('2008', '2009', '2010', '2011', '2012', '2013')) |>
  left_join(tax_bbmo_10y_new)

#write.csv(abund_asv_bloom_events_02, '../data/genes_sergio/abund_asv_bloom_events_02.csv') # removed the SAR11 clade that was not blooming

## plot them to see if there's some kind of correlation with KEGGs
abund_asv_bloom_events_02_unique <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  left_join(tax_bbmo_10y_new) |>
  distinct(asv_num)

abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'rclr' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value, date) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  left_join(tax_bbmo_10y_new)

abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value, date) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  left_join(tax_bbmo_10y_new)

#write.csv(abund_asv_bloo_fl_09_13, 'data/genes_sergio/abund_asv_bloo_fl_09_13.csv')

abund_asv_bloo_fl_09_13 <- abund_asv_bloo_fl_09_13 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

abund_asv_bloom_events_02 |>
  colnames()

abund_asv_bloom_events_02 <- abund_asv_bloom_events_02 |>
  distinct(sample_id, asv_num) |>
  left_join(m_bbmo_10y_ed_4metag) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

bloo_02_09_13_plot <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2008','2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(asv_num %in% abund_asv_bloom_events_02_unique$asv_num) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  group_by(date) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund, group = asv_num))+
  geom_vline(xintercept = abund_asv_bloom_events_02$date)+
  #facet_wrap(vars(asv_num), ncol = 1)+
  scale_x_datetime(expand = c(0,0),
                   # breaks = seq(min(asv_tab_bloo_rel_abund_z$date), 
                   #              max(asv_tab_bloo_rel_abund_z$date), by = '1 year'),
                   date_labels = "%Y")+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  #geom_line(aes(group = asv_num), color = '#5B5A5A')+ #, color = '#3D3B3B'
  # geom_point(data = abund_asv_bloo_fl_09_13 |>
  #              dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013') &
  #                              asv_num %in% abund_asv_bloom_events_02_unique$asv_num),
  #            aes(date, abundance_value, color =  '#9F0011'), size = 1)+
  scale_y_continuous(labels = percent_format())+
  geom_area(aes(fill = family_f, group = asv_num_f, y = abundance_value), position = 'stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(colour = '', alpha = '', y = 'Relative abundance (%)', x = 'Date')+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
       axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'bottom')

bloo_02_09_13_plot 

# ggsave( plot = bloo_02_09_13_plot, filename = 'bloo_02_09_13_plot.pdf',
#         path = 'results/figures/relationship_genes_blooms/',
#         width = 180, height = 160, units = 'mm')

grid.arrange(genes_plot,
             genes_caymes_plot,
             bloo_02_09_13_plot, ncol = 1)

bloo_02_09_13_plot <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(asv_num %in% abund_asv_bloom_events_02_unique$asv_num) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  # group_by(date) |>
  # dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, abundance_value, group = asv_num))+
  #geom_vline(xintercept = abund_asv_bloom_events_02$date)+
  facet_wrap(vars(interaction(asv_num, family_f)), ncol = 2)+
  scale_x_datetime(expand = c(0,0),
                   breaks =  '1 year',
                   date_labels = "%Y")+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  #geom_line(aes(group = asv_num), color = '#5B5A5A')+ #, color = '#3D3B3B'
  # geom_point(data = abund_asv_bloo_fl_09_13 |>
  #              dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013') &
  #                              asv_num %in% abund_asv_bloom_events_02_unique$asv_num &
  #                              abundance_type == 'rclr'),
  #            aes(date, abundance_value, color =  '#9F0011'), size = 1)+
  #scale_y_continuous(labels = percent_format())+
  geom_area(#aes(fill = family_f, group = asv_num_f, y = abundance_value), 
            position = 'stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(colour = '', alpha = '', y = 'rCLR', x = 'Date')+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
         axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'bottom')

bloo_02_09_13_plot 

# ggsave( plot = bloo_02_09_13_plot, filename = 'bloo_02_09_13_plot_rclr.pdf',
#         path = 'results/figures/relationship_genes_blooms/',
#         width = 180, height = 160, units = 'mm')

### cazymes and blooming events -----
#upload data from metagenomes we have the table with counts for the genes
genes_cazy <- read_table('data/genes_sergio/BBMOSOLA-GC_250bp_CAZy.lengthNorm.SCGnorm.counts.tbl') |>
  dplyr::select(annot, starts_with('BL')) |> #filter only Blanes samples not SOLA 
  as_tibble() |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'gene_counts') |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  group_by(sample_id)  |>
  arrange(sample_id) |>
  mutate(sample_id_num = cur_group_id()) |>
  left_join(m_bbmo_10y_ed_4metag, by = c('sample_id_num')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013'))

##filtrar GH GT 
genes_sample_id <- genes_cazy %$%
  sample_id |>
  unique() #I have 60 samples in common with metaG data (remember that the days from one dataset and the other sometimes do not match!)

## explore a bit this type of data----
genes_cazy %$%
  gene_counts |>
  range() #0.00000 2314.334

genes_cazy <- genes_cazy |>
  dplyr::mutate(annot = as.factor(annot))

genes_cazy$annot <- factor(genes_cazy$annot, levels = fct_infreq(levels(genes_cazy$annot)))

genes_cazy |>
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 5)) |> #filter by genes that have a minimum count at some point
  ggplot(aes(date, gene_counts))+
  facet_wrap(vars(annot), scales = 'free_y')+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  geom_line(aes(group = annot), color = '#5B5A5A')+ #, color = '#3D3B3B'
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 3/7)

## we would like to detect anomalies in genes and see if they correlate with anomalies in ASVs that are selected as potential bloomers
z_02_genes <- genes_cazy |>
  group_by(annot) |>
  arrange(date) |>
  group_by(annot) |>
  dplyr::reframe(anomalies_gene = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = gene_counts, plotting = FALSE)[c(1,2,3)])

z_scores_02_genes <- genes_cazy |>
  group_by(annot) |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  left_join(m_bbmo_10y_ed_4metag, by = c('sample_id_num')) |>
  #arrange(date) |>
  group_by(annot) |>
  dplyr::reframe(z_score_gen = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                             na_rm = TRUE, values = gene_counts, 
                                             plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_gen) |>
  group_by(annot) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_bbmo_10y_ed_4metag))) |>
  dplyr::mutate(sample_id_num = as.integer(sample_id_num)) |>
  left_join(m_bbmo_10y_ed_4metag, by = 'sample_id_num') 

## filter by the years that we share after computing the z_scores so that I am able to calculate the last ones (crec que no cal perque els NAs els trobaré al final)

## plot the genes with their anomalies an observe if they have correspondence with the blooming events----
z_scores_02_genes |>
  colnames()

z_scores_02_genes |>
  dim()

genes_cazy |>
  left_join(z_scores_02_genes) %$%
  gene_counts |>
  range()

genes_cazy |>
  left_join(z_scores_02_genes) |>
  dplyr::filter(gene_counts > 17)

z_scores_02_genes_f <- z_scores_02_genes |>
  dplyr::mutate(z_score_gen = case_when(is.na(z_score_gen) ~ 0,
                                        z_score_gen == 'NaN' ~ 0,
                                        z_score_gen == Inf ~ 10000,
                                        TRUE ~ z_score_gen)) |>
  dplyr::mutate(z_score_gen = as.numeric(z_score_gen)) |>
  group_by(annot) |>
  dplyr::filter(any(z_score_gen > 1.96, na.rm = TRUE)) |> # only plot those genes that present an anomaly at some point
  dplyr::select(annot, z_score_gen, sample_id_sim, date)

genes_anom <- z_scores_02_genes_f %$%
  annot |>
  unique() 

genes |>
  colnames()

z_scores_02_genes_f |>
  colnames()

genes_cazy$sample_id.x
z_scores_02_genes_f$sample_id_sim

genes_counts_z_scores <- genes_cazy |>
  filter(annot %in% genes_anom) |>
  dplyr::select(sample_id.x, gene_counts, annot, date) |>
  left_join(z_scores_02_genes_f, by = c('sample_id.x' = 'sample_id_sim', 'annot')) ##

genes_counts_z_scores_anom <- genes_counts_z_scores |>
  group_by(annot) |>
  dplyr::filter(z_score_gen > 1.96 &
                  any(gene_counts > 5))

##cazymes plot ----
abund_asv_bloom_events_02_het <- abund_asv_bloom_events_02 |>
  dplyr::filter(asv_num == 'asv11')

genes_cazymes_plot <- genes_counts_z_scores |>  
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 3)) |> #filter by genes that have a minimum count at some point
  # dplyr::filter(any(z_score_gen > 150)) |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::filter(str_detect(annot, 'GH')) |>
  ungroup() |>
  #dplyr::mutate(date.x = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_gen) >= 1.96, '#9F0011', '#080808', missing = '#080808')) |>
  group_by(date.x) |>
  dplyr::mutate(max_counts = sum(gene_counts)) |>
  ggplot(aes(date.x, max_counts))+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0))+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  # geom_point(aes(color = anomaly_color))+
  geom_vline(xintercept = abund_asv_bloom_events_02_het$date, color = '#8B8989')+
  geom_area(aes(group = fct_infreq(annot), fill = annot, y = gene_counts))+
  # geom_point(data = genes_counts_z_scores_anom,
  #   aes(date.x, gene_counts, color =  '#9F0011'), size = 1)+
  scale_fill_manual(values = palette_genes)+
  labs(x = 'Time (Y)', y = 'Gene counts', title = 'CAZYMES GH')+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'none')

genes_cazymes_plot

genes_cazymes_plot <- genes_counts_z_scores |>  
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 3)) |> #filter by genes that have a minimum count at some point
  # dplyr::filter(any(z_score_gen > 150)) |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::filter(str_detect(annot, 'GT')) |>
  ungroup() |>
  dplyr::mutate(annot = as.factor(annot)) |>
  #dplyr::mutate(date.x = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_gen) >= 1.96, '#9F0011', '#080808', missing = '#080808')) |>
  group_by(date.x) |>
  dplyr::mutate(max_counts = sum(gene_counts)) |>
  ggplot(aes(date.x, max_counts, group = fct_infreq(annot)))+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0))+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  # geom_point(aes(color = anomaly_color))+
  geom_vline(xintercept = abund_asv_bloom_events_02_het$date, color = '#8B8989')+
  geom_area(aes(group = fct_infreq(annot), fill = annot, y = gene_counts))+
  # geom_point(data = genes_counts_z_scores_anom,
  #   aes(date.x, gene_counts, color =  '#9F0011'), size = 1)+
  scale_fill_manual(values = palette_genes)+
  labs(x = 'Time (Y)', y = 'Gene counts', title = 'CAZYMES GT')+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'bottom')

genes_cazymes_plot

genes_cazymes_plot <- genes_counts_z_scores |>  
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 1)) |> #filter by genes that have a minimum count at some point
  # dplyr::filter(any(z_score_gen > 150)) |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::filter(str_detect(annot, 'PL')) |>
  ungroup() |>
  #dplyr::mutate(date.x = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_gen) >= 1.96, '#9F0011', '#080808', missing = '#080808')) |>
  group_by(date.x) |>
  dplyr::mutate(max_counts = sum(gene_counts)) |>
  ggplot(aes(date.x, max_counts))+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0))+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  # geom_point(aes(color = anomaly_color))+
  geom_vline(xintercept = abund_asv_bloom_events_02_het$date, color = '#8B8989')+
  geom_area(aes(group = annot, fill = annot, y = gene_counts))+
  # geom_point(data = genes_counts_z_scores_anom,
  #   aes(date.x, gene_counts, color =  '#9F0011'), size = 1)+
  scale_fill_manual(values = palette_genes)+
  labs(x = 'Time (Y)', y = 'Gene counts', title = 'CAZYMES PL')+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'bottom')

genes_cazymes_plot

genes_cazymes_plot <- genes_counts_z_scores |>  
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 1)) |> #filter by genes that have a minimum count at some point
  # dplyr::filter(any(z_score_gen > 150)) |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::filter(str_detect(annot, 'AA')) |>
  ungroup() |>
  #dplyr::mutate(date.x = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_gen) >= 1.96, '#9F0011', '#080808', missing = '#080808')) |>
  group_by(date.x) |>
  dplyr::mutate(max_counts = sum(gene_counts)) |>
  ggplot(aes(date.x, max_counts))+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0))+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  # geom_point(aes(color = anomaly_color))+
  geom_vline(xintercept = abund_asv_bloom_events_02_het$date, color = '#8B8989')+
  geom_area(aes(group = annot, fill = annot, y = gene_counts))+
  # geom_point(data = genes_counts_z_scores_anom,
  #   aes(date.x, gene_counts, color =  '#9F0011'), size = 1)+
  scale_fill_manual(values = palette_genes)+
  labs(x = 'Time (Y)', y = 'Gene counts', title = 'CAZYMES AA')+
  #scale_color_identity()+
  # geom_line(aes(group = annot), color = '#5B5A5A')+ #, color = '#3D3B3B'
  #facet_wrap(vars(annot), ncol = 2)+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'none')

genes_cazymes_plot

##### 
max_abund_asv_bloom_events |>
  dim() == samples_with_bloom_event |>
  dim() ##ROW should be TRUE (row - cols)

samples_with_bloom_event_asv_abund <- samples_with_bloom_event |>
  right_join(max_abund_asv_bloom_events)

## relationship between genes and blooming ASVs----
### I need the to compare normalized genes and normalized ASVs abundances.
rclr_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(asv_num %in% abund_asv_bloom_events_02_unique$asv_num) |>
  dplyr::filter(abundance_type == 'rclr')

genes_counts_z_scores_f <- genes_counts_z_scores |>  
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 5)) |> ## I do it for the most abundant genes to simplify
dplyr::select(-date.y) |>
  dplyr::mutate(date.x = (as.POSIXct(date.x, format = "%Y-%m-%d")))
  
zscores_gen_ra_fl_09_13 <- rclr_bloo_fl_09_13 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(genes_counts_z_scores_f, by = c('date' = 'date.x'), relationship = 'many-to-many') 

zscores_gen_ra_fl_09_13 |>
  colnames() 

zscores_gen_ra_fl_09_13 |>
  #dplyr::filter(!is.na(annot)) |>
  ggplot(aes(z_score_ra, z_score_gen))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  facet_wrap(asv_num~annot, scales = 'free')+
  theme(strip.text = element_text(size = 5), axis.text = element_text(size = 3))

##another option of plotting the same see if it is more clear
zclr_bloo_fl_09_13_sim <- zclr_bloo_fl_09_13 |>
  dplyr::select(date, z_score_ra, asv_num) |>
  rename(z_score = z_score_ra, variable = asv_num)

genes_counts_z_scores_f_sim <- genes_counts_z_scores_f |>
  dplyr::select(date.x, z_score_gen, annot) |>
  rename(date = date.x, z_score = z_score_gen, variable = annot)

zclr_bloo_fl_09_13_sim |>
  bind_rows(genes_counts_z_scores_f_sim) |>
  ggplot(aes(date, z_score))+
  geom_point()+
  geom_line(aes(group = 'annot'))+
  facet_grid(vars(variable), scales = 'free')

#### bray-curtis from Ramiro's data------
bc_kmrs <- read_delim('data/genes_sergio/mat_abundance_braycurtis.csv.gz', delim = ';')

bc_kmrs_red <- bc_kmrs |>
  dplyr::select(1:61) |>
  dplyr::mutate(sample_num = row_number()) |>
  dplyr::filter(sample_num %in% 1:60) 

sample_id_num_2 <- bc_kmrs_red |>
  dplyr::select(sample_num, ...1) |>
  dplyr::mutate(sample_num_2 = as.character(sample_num)) |>
  dplyr::select(-sample_num)

bc_ramiro <- bc_kmrs_red |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id_1', values_to = 'bc') |>
  left_join(sample_id_num_2, by = c('sample_id_1' = '...1')) |>
  dplyr::filter(as.numeric(sample_num) == (as.numeric(sample_num_2)-1)) |>
  dplyr::mutate(sample_num_2 = as.numeric(sample_num_2)) |>
  left_join(m_bbmo_10y_ed_4metag, by = c( 'sample_num_2' = 'sample_id_num'))

abund_asv_bloom_events_02 |>
  distinct(sample_id)|>
  separate(sample_id, into = c('sample_id', 'fraction', 'id'), sep = '_') 

bc_ramiro |>
  str()

bc_ramiro |>
  ungroup() |>
  #right_join(dates_bbmo, by = c('sample_id_1' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, bc))+
  geom_line()+
  geom_point()+
  #scale_x_datetime()+
  theme_bw()

## plot bc Ramiro in the same plot as the bloomers abundances----
asv_tab_all_bloo_z_tax$sample_id |>
  colnames()

bc_blooms <- asv_tab_all_bloo_z_tax |>
  separate(sample_id, into = c('sample_id', 'fraction', 'id'), sep = '_') |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(asv_num %in% abund_asv_bloom_events_02_unique$asv_num) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(-month, -day) |>
  left_join(bc_ramiro, by = 'date')

bloomers_bc <- bc_blooms |>
  group_by(date) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund, group = asv_num))+
  #geom_vline(xintercept = abund_asv_bloom_events_02$date, color = '#8B8989' )+
  facet_wrap(vars(asv_num_f), ncol = 2)+
  scale_x_datetime(expand = c(0,0),
                   # breaks = seq(min(asv_tab_bloo_rel_abund_z$date), 
                   #              max(asv_tab_bloo_rel_abund_z$date), by = '1 year'),
                   date_labels = "%Y")+
  geom_area(aes(fill = family_f, group = asv_num_f, y = abundance_value), position = 'stack')+
  geom_point(data = bc_blooms |>
               dplyr::filter(z_score_ra > 1.96 &
                               abundance_value > 0.1), aes(date, abundance_value, color =  '#9F0011'), size = 1, alpha = 1)+
  scale_fill_manual(values = palette_family_assigned_bloo, guide = guide_legend(keywidth = unit(2, "mm"),keyheight = unit(2, "mm"),
                                                                                ncol = 2))+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(~.* 2 , name = 'Bray-Curtis dissimilarity'))+
  geom_line(aes(y = bc/2), color = '#5B5A5A')+ #, color = '#3D3B3B'
  labs(fill = 'Family', alpha = '', y = 'Relative abundance (%)', x = 'Date')+
  theme_bw()+
  guides(color = "none",
         alpha = "none")+
  theme(panel.grid = element_blank(), text = element_text(size = 7),
        #panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_rect(fill = 'transparent'), legend.position = 'bottom')
  # guides(guide_legend(colour = 'none',
  #                     fill = guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)))

bloomers_bc

# ggsave(bloomers_bc, path = 'results/figures/',
#        file = 'bloomers_09_13_FL.pdf',
#        width = 180, height = 120, units = 'mm')

genes_bloomers <- plot_grid(bloomers_bc,
                            plot_grid(genes_cazymes_plot, genes_photo_plot, ncol = 1),
                            ncol = 2,
                            labels = c('A', "B", "C"),  # Remove label for the first plot in the first column
                            hjust = -1)

# ggsave(genes_bloomers, path = '../results/figures/', 
#        file = 'genes_bloomers.pdf',
#        width = 150, height = 160, units = 'mm')

### Look only for photosynthetic KO's------
##### Photosystem II [PATH:map00195] [BR:ko00194]: K02703, K02706, K02705, K02704, K02707, K02708, K02689, K02690, K02691, K02692, K02693, K02694
##### Anoxygenic photosystem II [BR:ko00194]: K08928, K08929, K08940, K08941, K08942, K08943
##### NblS-NblR (photosynthesis) two-component regulatory system [PATH:map02020] [BR:ko02022]	K07769
##### NblS-NblR (photosynthesis) two-component regulatory system [PATH:map02020] [BR:ko02022]	K11332

kegs_table <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.metaGsizeNorm.counts.txt', header = T) |>
  dplyr::select(annot, starts_with('BL')) 

photosystem_ii <- kegs_table |>
  dplyr::filter(annot %in% c('K02703', 'K02706', 'K02705', 'K02704', 'K02707', 'K02708', 'K02689', 
                             'K02690', 'K02691', 
                             'K02692', 'K02693', 'K02694')) |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'gene_counts') |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  group_by(sample_id)  |>
  arrange(sample_id) |>
  mutate(sample_id_num = cur_group_id()) |>
  left_join(m_bbmo_10y_ed_4metag_red, by = c('sample_id_num')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013'))

abund_asv_bloom_events_02_photo <- abund_asv_bloom_events_02 |>
  dplyr::filter(asv_num %in% c('asv1', 'asv7'))
  
photosystem_ii_plot <- photosystem_ii |>
  group_by(date) |>
  dplyr::mutate(max_counts = sum(gene_counts)) |>
  ggplot(aes(date, max_counts))+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0))+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  # geom_point(aes(color = anomaly_color))+
  geom_vline(xintercept = abund_asv_bloom_events_02_photo$date, color = '#8B8989')+
  geom_area(aes(group = annot, fill = annot, y = gene_counts))+
  # geom_point(data = genes_counts_z_scores_anom,
  #   aes(date.x, gene_counts, color =  '#9F0011'), size = 1)+
  scale_fill_manual(values = palette_genes)+
  labs(y = 'Gene counts', x = 'Time (Y)', title = 'Photosystem II')+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 6), panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'none')

## potential bloomers that presenta photosystem II -----
bloomers_bc_pII <- bc_blooms |>
  dplyr::filter(asv_num %in% c('asv1', 'asv7')) |>
  group_by(date) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund, group = asv_num))+
  #geom_vline(xintercept = abund_asv_bloom_events_02$date, color = '#8B8989' )+
  #facet_wrap(vars(asv_num_f), ncol = 1)+
  scale_x_datetime(expand = c(0,0),
                   # breaks = seq(min(asv_tab_bloo_rel_abund_z$date), 
                   #              max(asv_tab_bloo_rel_abund_z$date), by = '1 year'),
                   date_labels = "%Y")+
  geom_area(aes(fill = asv_num_f, group = asv_num_f, y = abundance_value), position = 'stack')+
  # geom_point(data = bc_blooms |>
  #            dplyr::filter(asv_num %in% c('asv1', 'asv7')) |>
  #              dplyr::filter(z_score_ra > 1.96), aes(date, abundance_value, color =  '#9F0011'), size = 1, alpha = 1)+
  scale_fill_manual(values = c('#000000', "#57a9a8"), guide = guide_legend(keywidth = unit(2, "mm"),keyheight = unit(2, "mm"),
                                                                                ncol = 2))+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(expand = c(0,0),
                     #sec.axis = sec_axis(~.* 4 , name = 'Bray-Curtis dissimilarity')
                     )+
  #geom_line(aes(y = bc/4), color = '#5B5A5A')+ #, color = '#3D3B3B'
  labs(fill = 'Family', alpha = '', y = 'Relative abundance (%)', x = 'Date')+
  theme_bw()+
  guides(color = "none",
         alpha = "none")+
  theme(panel.grid = element_blank(), text = element_text(size = 7),
        #panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_rect(fill = 'transparent'), legend.position = 'none')
# guides(guide_legend(colour = 'none',
#                     fill = guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)))

bloomers_bc_pII

genes_photosyntesis_bloo <- grid.arrange(photosystem_ii_plot, bloomers_bc_pII, ncol = 1,
             heights = c(1,1))

ggsave(genes_photosyntesis_bloo, 
       path = 'results/figures/',
       file = 'genes_phtosynyesis_bloo_together.pdf',
       width = 180, height = 200, units = 'mm')

## whole timeseries taxa that could present a phoyosystem II ----
asv_tab_10y_filt_02_rel_rar |>
  colnames()

asv_tab_bbmo_cyano <- asv_tab_10y_filt_02_rel_rar |>
  left_join(m_bbmo_10y) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(class == 'Cyanobacteriia')

cyanos_bbmo_genes <- asv_tab_bbmo_cyano |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(any(relative_abundance > 0)) |>
  group_by(date) |>
  dplyr::mutate(max_abund = sum(relative_abundance)) |>
  ggplot(aes(date, max_abund, group = asv_num))+
  #geom_vline(xintercept = abund_asv_bloom_events_02$date, color = '#8B8989' )+
  #facet_wrap(vars(asv_num_f), ncol = 1)+
  scale_x_datetime(expand = c(0,0),
                   # breaks = seq(min(asv_tab_bloo_rel_abund_z$date), 
                   #              max(asv_tab_bloo_rel_abund_z$date), by = '1 year'),
                   date_labels = "%Y")+
  geom_area(aes(fill = family, group = asv_num, y = relative_abundance), position = 'stack')+
  # geom_point(data = bc_blooms |>
  #            dplyr::filter(asv_num %in% c('asv1', 'asv7')) |>
  #              dplyr::filter(z_score_ra > 1.96), aes(date, abundance_value, color =  '#9F0011'), size = 1, alpha = 1)+
  scale_fill_manual(values = palette_family_assigned_bloo, guide = guide_legend(keywidth = unit(2, "mm"),keyheight = unit(2, "mm"),
                                                                           ncol = 2))+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(expand = c(0,0),
                     #sec.axis = sec_axis(~.* 4 , name = 'Bray-Curtis dissimilarity')
  )+
  #geom_line(aes(y = bc/4), color = '#5B5A5A')+ #, color = '#3D3B3B'
  labs(fill = 'Family', alpha = '', y = 'Relative abundance (%)', x = 'Date')+
  theme_bw()+
  guides(color = "none",
         alpha = "none")+
  theme(panel.grid = element_blank(), text = element_text(size = 7),
        #panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_rect(fill = 'transparent'), legend.position = 'none')

genes_photosyntesis_whole_community <- grid.arrange(photosystem_ii_plot, 
                                                    cyanos_bbmo_genes, ncol = 1,
                                                    heights = c(1,1))

ggsave(genes_photosyntesis_whole_community, 
       path = 'results/figures/',
       file = 'genes_photosyntesis_whole_community.pdf',
       width = 180, height = 200, units = 'mm')

## Do blooming events have an impact on the available ecosystemic functions?-----

### Relationships between cazymes and heterotrophic bloomers----
cazy_GT <- genes_cazy |>  
  pivot_longer(cols = c(!annot), names_to = 'sample_id', values_to = 'gene_counts') |>
  left_join(m_bbmo_10y_ed_4metag_red, by = c('sample_id' = 'sample_id_sim')) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 1)) |> #filter by genes that have a minimum count at some point
  # dplyr::filter(any(z_score_gen > 150)) |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::filter(str_detect(annot, 'GT')) |>
  ungroup() |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  #dplyr::mutate(anomaly_color = if_else(abs(z_score_gen) >= 1.96, '#9F0011', '#080808', missing = '#080808')) |>
  group_by(date) |>
  dplyr::mutate(max_counts = sum(gene_counts)) |>
  ggplot(aes(date, max_counts))+
  scale_x_datetime(expand = c(0,0),
                   breaks = seq(min(bc_blooms$date),
                                max(bc_blooms$date), by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0))+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  # geom_point(aes(color = anomaly_color))+
  geom_vline(xintercept = abund_asv_bloom_events_02_het$date, color = '#8B8989')+
  geom_area(aes(group = annot, fill = annot, y = gene_counts))+
  # geom_point(data = genes_counts_z_scores_anom,
  #   aes(date.x, gene_counts, color =  '#9F0011'), size = 1)+
  #scale_fill_manual(values = palette_genes)+
  labs(x = 'Time (Y)', y = 'Gene counts', title = 'CAZYMES GT')+
  #scale_color_identity()+
  # geom_line(aes(group = annot), color = '#5B5A5A')+ #, color = '#3D3B3B'
  #facet_wrap(vars(annot), ncol = 2)+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'none')

cazy_GH <- genes_cazy |>  
  pivot_longer(cols = c(!annot), names_to = 'sample_id', values_to = 'gene_counts') |>
  left_join(m_bbmo_10y_ed_4metag_red, by = c('sample_id' = 'sample_id_sim')) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 1)) |> #filter by genes that have a minimum count at some point
  # dplyr::filter(any(z_score_gen > 150)) |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::filter(str_detect(annot, 'GH')) |>
  ungroup() |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  #dplyr::mutate(anomaly_color = if_else(abs(z_score_gen) >= 1.96, '#9F0011', '#080808', missing = '#080808')) |>
  group_by(date) |>
  dplyr::mutate(max_counts = sum(gene_counts)) |>
  ggplot(aes(date, max_counts))+
  scale_x_datetime(expand = c(0,0),
                   breaks = seq(min(bc_blooms$date),
                                max(bc_blooms$date), by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0))+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  # geom_point(aes(color = anomaly_color))+
  geom_vline(xintercept = abund_asv_bloom_events_02_het$date, color = '#8B8989')+
  geom_area(aes(group = annot, fill = annot, y = gene_counts))+
  # geom_point(data = genes_counts_z_scores_anom,
  #   aes(date.x, gene_counts, color =  '#9F0011'), size = 1)+
  #scale_fill_manual(values = palette_genes)+
  labs(x = 'Time (Y)', y = 'Gene counts', title = 'CAZYMES GH')+
  #scale_color_identity()+
  # geom_line(aes(group = annot), color = '#5B5A5A')+ #, color = '#3D3B3B'
  #facet_wrap(vars(annot), ncol = 2)+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'none')


## potential heterotrophic bloomers -----
bloomers_bc_asv11 <- bc_blooms |>
  dplyr::filter(asv_num %in% c('asv11')) |>
  group_by(date) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund, group = asv_num))+
  #geom_vline(xintercept = abund_asv_bloom_events_02$date, color = '#8B8989' )+
  #facet_wrap(vars(asv_num_f), ncol = 1)+
  scale_x_datetime(expand = c(0,0),
                   breaks = seq(min(bc_blooms$date),
                                max(bc_blooms$date), by = '1 year'),
                   date_labels = "%Y")+
  geom_area(aes(fill = asv_num_f, group = asv_num_f, y = abundance_value), position = 'stack')+
  # geom_point(data = bc_blooms |>
  #            dplyr::filter(asv_num %in% c('asv1', 'asv7')) |>
  #              dplyr::filter(z_score_ra > 1.96), aes(date, abundance_value, color =  '#9F0011'), size = 1, alpha = 1)+
  scale_fill_manual(values = c('#000000', "#57a9a8"), guide = guide_legend(keywidth = unit(2, "mm"),keyheight = unit(2, "mm"),
                                                                           ncol = 2))+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(expand = c(0,0),
                     #sec.axis = sec_axis(~.* 4 , name = 'Bray-Curtis dissimilarity')
  )+
  #geom_line(aes(y = bc/4), color = '#5B5A5A')+ #, color = '#3D3B3B'
  labs(fill = 'Family', alpha = '', y = 'Relative abundance (%)', x = 'Date')+
  theme_bw()+
  guides(color = "none",
         alpha = "none")+
  theme(panel.grid = element_blank(), text = element_text(size = 7),
        #panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_rect(fill = 'transparent'), legend.position = 'none')
# guides(guide_legend(colour = 'none',
#                     fill = guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)))

bloomers_bc_asv11

genes_cazy_plot <- grid.arrange(bloomers_bc_asv11, cazy_GT,
                                cazy_GH, ncol = 1,
                                heights = c(1,1,1))

## relationship genes and bloomers ----
photosystem_ii <- kegs_table |>
  dplyr::filter(annot %in% c('K02703', 'K02706', 'K02705', 'K02704', 'K02707', 'K02708', 'K02689', 
                             'K02690', 'K02691', 
                             'K02692', 'K02693', 'K02694')) |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'gene_counts') |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  group_by(sample_id)  |>
  arrange(sample_id) |>
  mutate(sample_id_num = cur_group_id()) |>
  left_join(m_bbmo_10y_ed_4metag_red, by = c('sample_id_num')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::select(annot, sample_id_num, year, month, day, date, gene_counts)
  
photo_bloo_tb <- bc_blooms |>
  dplyr::filter(asv_num %in% c('asv1', 'asv7')) |>
  group_by(date) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  dplyr::select(date, asv_num, abundance_value, sample_num_2) |>
  dplyr::mutate(sample_num_2 = case_when(is.na(sample_num_2) ~ 1, 
                                       !is.na(sample_num_2) ~ sample_num_2)) 

abund_asv_bloom_events_02 |>
  colnames()

bc_blooms |>
  colnames()

ko_photosystem_ii_abund <- photosystem_ii |>
  left_join(photo_bloo_tb, by = c( 'sample_id_num' = 'sample_num_2')) |>
  ggplot(aes(abundance_value, gene_counts))+
  geom_point(aes(color = asv_num))+
  facet_wrap(vars(annot), scales = 'free')+
  geom_smooth(method = 'lm', aes(group = asv_num, color = asv_num))+
  labs(fill = 'Family', alpha = '', x = 'Relative abundance', y = 'Gene counts',
       color = 'ASV number')+
  theme_bw()+
  guides(#color = "b",
         alpha = "none")+
  scale_color_manual(values = c('#000000', "#57a9a8"), guide = guide_legend(keywidth = unit(2, "mm"),keyheight = unit(2, "mm"),
                                                                           ncol = 2))+
  theme(panel.grid = element_blank(), text = element_text(size = 7),
        #panel.border = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 3/3,
        strip.background = element_rect(fill = 'transparent'), legend.position = 'bottom')
# guides(guide_legend(colour = 'none',
#                     fill = guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)))

ko_photosystem_ii_abund

# ggsave(ko_photosystem_ii_abund, 
#        path = 'results/figures/',
#        file = 'ko_photosystem_ii_abund.pdf',
#        width = 180, height = 200, units = 'mm')

## I will try the code localy (Error: vector memory exhausted (limit reached?)) ----
### This packages are useful for working with big tables, but in this case it's too much
library(vegan)
#library("data.table")
#library("parallelDist")
#library("bigmemory")

## I filter genes in MARBITS and then i move to local -----
### What I did is filter the genes table by different thresholds, to be able to see if the tendencies are the same even if i'm filtering some genes.

## threshold 0.001 -------
threshold_value <- 0.001
filtered_genes <- read.table('data/genes_sergio/filtered_0.1threshold_BBMO_genes_cleaned2.tbl') ## threshold 0.001

filtered_genes_names <- filtered_genes %$%
  V1

filtered_genes_t <- t(filtered_genes)

filtered_genes_t <- filtered_genes_t[-1,]

filtered_genes_t |>
  colnames() <- filtered_genes_names

filtered_genes_t <- filtered_genes_t |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

filtered_genes_t_l <- filtered_genes_t |>
  pivot_longer(starts_with('BBMO'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_genes_filtered <- dissimilarity_matrix(data = filtered_genes_t_l, 
                                          sample_id_col = sample_id,
                                          values_cols_prefix = 'BL')

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

bray_curtis_genes_filtered_v1 <-  bray_curtis_genes_filtered |>
  left_join(abund_asv_bloo_fl_09_13, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(genes_threshold = threshold_value)

bray_curtis_genes_filtered_v1 |>
  ggplot(aes(date, bray_curtis_result))+
  geom_point()+
  geom_line()

## threshold 0.0001 -------
threshold_value <- 0.0001
filtered_genes <- read.table('data/genes_sergio/filtered_0.0001threshold_BBMO_genes_cleaned.tbl') ## threshold 0.001

filtered_genes_names <- filtered_genes %$%
  V1

filtered_genes_t <- t(filtered_genes)

filtered_genes_t <- filtered_genes_t[-1,]

filtered_genes_t |>
  colnames() <- filtered_genes_names

filtered_genes_t <- filtered_genes_t |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

filtered_genes_t_l <- filtered_genes_t |>
  pivot_longer(starts_with('BBMO'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance)) |>
  dplyr::mutate(genes_threshold = threshold_value)

bray_curtis_genes_filtered <- dissimilarity_matrix(data = filtered_genes_t_l, 
                                                   sample_id_col = sample_id,
                                                   values_cols_prefix = 'BL')

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

bray_curtis_genes_filtered_v2 <-  bray_curtis_genes_filtered |>
  left_join(abund_asv_bloo_fl_09_13, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(genes_threshold = threshold_value)

bray_curtis_genes_filtered_v2 |>
  ggplot(aes(date, bray_curtis_result))+
  geom_point()+
  geom_line()

## threshold 0.00001 -------
threshold_value <- 0.00001
filtered_genes <- read.table('data/genes_sergio/filtered_0.00001threshold_BBMO_genes_cleaned.tbl') ## threshold 0.001

filtered_genes_names <- filtered_genes %$%
  V1

filtered_genes_t <- t(filtered_genes)

filtered_genes_t <- filtered_genes_t[-1,]

filtered_genes_t |>
  colnames() <- filtered_genes_names

filtered_genes_t <- filtered_genes_t |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

filtered_genes_t_l <- filtered_genes_t |>
  pivot_longer(starts_with('BBMO'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_genes_filtered <- dissimilarity_matrix(data = filtered_genes_t_l, 
                                                   sample_id_col = sample_id,
                                                   values_cols_prefix = 'BL')

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

bray_curtis_genes_filtered_v3 <-  bray_curtis_genes_filtered |>
  left_join(abund_asv_bloo_fl_09_13, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(genes_threshold = threshold_value)

bray_curtis_genes_filtered_v3 |>
  ggplot(aes(date, bray_curtis_result))+
  geom_point()+
  geom_line()

## threshold 0.000001 -------
threshold_value <- 0.000001
filtered_genes <- read.table('data/genes_sergio/filtered_0.000001threshold_BBMO_genes_cleaned.tbl') ## threshold 0.001

filtered_genes_names <- filtered_genes %$%
  V1

filtered_genes_t <- t(filtered_genes)

filtered_genes_t <- filtered_genes_t[-1,]

filtered_genes_t |>
  colnames() <- filtered_genes_names

filtered_genes_t <- filtered_genes_t |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

filtered_genes_t_l <- filtered_genes_t |>
  pivot_longer(starts_with('BBMO'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

filtered_genes_t_l |>
  distinct(sample_id) |>
  as_vector()

filtered_genes_t_l |>
  arrange(sample_id) |>
  distinct(sample_id) |>
  as_vector()

bray_curtis_genes_filtered <- dissimilarity_matrix(data = filtered_genes_t_l, 
                                                   sample_id_col = sample_id,
                                                   values_cols_prefix = 'BL')

bray_curtis_genes_filtered |>
  distinct(samples) |>
  as_vector()

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

bray_curtis_genes_filtered_v4 <-  bray_curtis_genes_filtered |>
  left_join(abund_asv_bloo_fl_09_13, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(genes_threshold = threshold_value)

bray_curtis_genes_filtered_v4 |>
  ggplot(aes(date, bray_curtis_result))+
  geom_point()+
  geom_line()

## v2 filter by at least 1 unique value higher than the threshold ----
threshold_value <- (2.05712/100000*2)
#bray_curtis_genes_filtered <- read.csv('data/genes_sergio/bray_curtis_genes_threshold0.001_v2.csv') ## threshold 0.001 calculat directament a MARBITS 
bray_curtis_genes_filtered <- read.csv('data/genes_sergio/bray_curtis_genes_threshold_meanx2_v2.csv') ## mean abund of all genes * 2 (well ordered!)

#source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R')
## he detectat que l'ordre de les files no está sempre igual!

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t)))

abund_asv_bloo_fl_09_13$sample_id == bray_curtis_genes_filtered$samples

bray_curtis_genes_filtered_v5 <-  bray_curtis_genes_filtered |>
  separate(samples, into = c('camapin', 'sample_num'), sep = 2 , remove = F) |>
  left_join(abund_asv_bloo_fl_09_13, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(genes_threshold = paste0(threshold_value, 'v2'))

bray_curtis_genes_filtered_v5 |>
  ggplot(aes(date, bray_curtis_result))+
  geom_point()+
  geom_line()

## comparison of thresholds data ----
bray_curtis_genes_all <- bray_curtis_genes_filtered_v1 |>
  bind_rows(bray_curtis_genes_filtered_v2) |>
  bind_rows(bray_curtis_genes_filtered_v3) |>
  bind_rows(bray_curtis_genes_filtered_v4) |>
  dplyr::mutate(genes_threshold = as.character(genes_threshold)) |>
  bind_rows(bray_curtis_genes_filtered_v5)

comparison_bray_thresholds <- bray_curtis_genes_all |>
  ggplot(aes(date, bray_curtis_result))+
  #geom_point(aes(group = genes_threshold, color = as.character(genes_threshold)))+
  geom_line(aes(group = genes_threshold, color = as.character(genes_threshold)))+
  scale_color_manual(values = palette_genes)+
  labs(x = 'Date', y = 'Genes Bray Curtis Dissimilarity', color = 'Threshold')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), #panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        panel.border = element_blank()) 

comparison_bray_thresholds

# ggsave(comparison_bray_thresholds,
#        filename = 'comparison_bray_thresholds_plot.pdf',
#        path = 'Results/Figures/',
#        width = 180, height = 160, units = 'mm'
# )

## Plot community bray curtis with genes bray curtis + euclidean distance + unifrac distance -----
palette_diversity <- c('genes' = "#270000",
                       "wunifrac_distance" = "#898989",
                       "euclidean_distance" = "#008241",
                       "bray_curtis_community" = "#a21e66", 
                       'bray_curtis_kmers' = '#FFCAFB')


labs_diversity <- as_labeller(c('genes' = 'Bray Curtis Genes',
                                "wunifrac_distance" = 'Weigthed UNIFRAC',
                                "euclidean_distance" = 'Euclidean Environmental Distance',
                                "bray_curtis_community" = 'Bray Curtis Community Composition',
                                'bray_curtis_kmers' = 'Bray Curtis k-mers'))

bray_unifrac_eucl_tb$bray_curtis_type |>
  unique()

### prepare data 
# bray_curtis_genes_all |>
#   colnames()
# 
# unifrac_tibble_m |>
#   colnames()
# 
# euclidean_distance_tb |>
#   colnames()
# 
# bray_curtis_m_tb |>
#   colnames()
# 
# bray_curtis_genes_all$genes_threshold |>
#   unique()

bray_curtis_genes_all_red <- bray_curtis_genes_all |>
  dplyr::select( date, bray_curtis_result, genes_threshold) |>
  dplyr::mutate(bray_curtis_type = 'genes') |>
  dplyr::mutate(fraction = '0.2') |>
  dplyr::filter(genes_threshold == '4.11424e-05v2') |> # at least one gene value higher than twice the mean of gens over the dataset
  dplyr::select(-genes_threshold)

unifrac_tibble_m_red <- unifrac_tibble_m |>
  dplyr::select(date, bray_curtis_result = wunifrac_distance, fraction = fraction.x) |>
  dplyr::mutate(bray_curtis_type = 'wunifrac_distance')

euclidean_distance_tb_red <- euclidean_distance_tb |>
  dplyr::select(date, bray_curtis_result = euclidean_distance) |>
  dplyr::mutate(bray_curtis_type = 'euclidean_distance',
                fraction = 'env')

bray_curtis_m_tb_red <- bray_curtis_m_tb |>
  dplyr::filter(diversity_index == 'bray_curtis_result') |>
  dplyr::select(date, bray_curtis_result = value, fraction) |>
  dplyr::mutate(bray_curtis_type = 'bray_curtis_community')

bc_ramiro_red <- bc_ramiro |>
  dplyr::select(date, bray_curtis_result =  bc) |>
  dplyr::mutate(bray_curtis_type = 'bray_curtis_kmers',
                fraction = '0.2')

# bray_curtis_genes_all_red <- bray_curtis_genes_all_red |>
#   dplyr::mutate(date = as.Date(date)) |>
#   dplyr::mutate(date = as.character(date))

bray_unifrac_eucl_tb <- bray_curtis_genes_all_red |>
  bind_rows(unifrac_tibble_m_red) |>
  bind_rows(euclidean_distance_tb_red) |>
  bind_rows(bray_curtis_m_tb_red) |>
  bind_rows(bc_ramiro_red)

# plot weighted unifrac distances -----

bray_unifrac_eucl_plot <- bray_unifrac_eucl_tb |>
  dplyr::filter(bray_curtis_type != 'bray_curtis_kmers') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  ggplot(aes(date, bray_curtis_result))+
  facet_wrap(vars(fraction), scales = 'free_y', ncol = 1, labeller = labs_fraction_env)+
  geom_line(aes(date, group = bray_curtis_type, color = bray_curtis_type), linewidth = 1, alpha = 0.8)+ #, linetype = bray_curtis_type
  #scale_linetype_discrete(labels = labs_diversity)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  scale_color_manual(values= palette_diversity, labels = labs_diversity)+ #, labels = labs_fraction
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = '', color = '', linetype = '')+
  guides(shape = 'none',
         color = guide_legend(ncol = 2))+
  #guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank())

bray_unifrac_eucl_plot

# ggsave(bray_unifrac_eucl_plot,
#        filename = 'bray_unifrac_eucl_plot_v5.pdf',
#               path = 'Results/Figures/',
#               width = 180, height = 130, units = 'mm'
#        )

bray_unifrac_eucl_plot <- bray_unifrac_eucl_tb |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  dplyr::filter(fraction == '0.2') |>
  ggplot(aes(date, bray_curtis_result))+
  facet_wrap(vars(bray_curtis_type), scales = 'free_y', ncol = 1, labeller = labs_diversity)+
  geom_line(aes(date, group = bray_curtis_type, color = bray_curtis_type), linewidth = 1, alpha = 0.8)+ #, linetype = bray_curtis_type
  #scale_linetype_discrete(labels = labs_diversity)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  scale_color_manual(values= palette_diversity, labels = labs_diversity)+ #, labels = labs_fraction
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = '', color = '', linetype = '')+
  guides(shape = 'none',
         color = guide_legend(ncol = 2))+
  #guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank())

bray_unifrac_eucl_plot

# ggsave(bray_unifrac_eucl_plot,
#        filename = 'bray_unifrac_facet_eucl_plot_v3.pdf',
#               path = 'Results/Figures/',
#               width = 180, height = 160, units = 'mm'
#        )

## relationship between parameters ----

bray_unifrac_eucl_tb_02 <-  bray_unifrac_eucl_tb |>
  dplyr::filter(fraction %in% c('env', '0.2')) |>
  dplyr::mutate(date = as.Date(date)) |>
  pivot_wider(id_cols = c('date'), names_from = 'bray_curtis_type', values_from = 'bray_curtis_result')

bray_unifrac_eucl_tb_02 |>
  #dplyr::mutate(euclidean_distance = scale(euclidean_distance)) |>
  ggplot(aes(euclidean_distance, bray_curtis_community))+
  geom_point()+
  geom_smooth(method = 'loess')

bray_unifrac_eucl_tb_02 |>
  dplyr::mutate(euclidean_distance = scale(euclidean_distance)) |>
  ggplot(aes(wunifrac_distance, bray_curtis_community))+
  geom_point()+
  geom_smooth(method = 'loess')

bray_unifrac_eucl_tb_02 |>
  dplyr::mutate(euclidean_distance = scale(euclidean_distance)) |>
  ggplot(aes(wunifrac_distance, genes))+
  geom_point()+
  geom_smooth(method = 'loess')

bray_unifrac_eucl_tb_02 |>
  dplyr::mutate(euclidean_distance = scale(euclidean_distance)) |>
  ggplot(aes(bray_curtis_community, genes))+
  geom_point()+
  geom_smooth(method = 'loess')

bray_unifrac_eucl_tb_3 <- bray_unifrac_eucl_tb |>
  dplyr::filter(fraction %in% c('3', 'env')) |>
  dplyr::mutate(date = as.character(date)) |>
  pivot_wider(id_cols = c('date'), names_from = 'bray_curtis_type', values_from = 'bray_curtis_result')

bray_unifrac_eucl_tb_3 |>
  #dplyr::mutate(euclidean_distance = scale(euclidean_distance)) |>
  ggplot(aes(euclidean_distance, bray_curtis_community))+
  geom_point()+
  geom_smooth(method = 'loess')

bray_unifrac_eucl_tb_3 |>
  dplyr::mutate(euclidean_distance = scale(euclidean_distance)) |>
  ggplot(aes(wunifrac_distance, bray_curtis_community))+
  geom_point()+
  geom_smooth(method = 'loess')+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    #panel.grid.major.y = element_blank(),
    panel.border = element_blank())

### both fractions at the same time ----
bray_unifrac_eucl_tb_02 <- bray_unifrac_eucl_tb_02 |>
  dplyr::mutate(fraction = '0.2') |>
  dplyr::mutate(date = as.character(date))

bray_unifrac_eucl_tb_w <- bray_unifrac_eucl_tb_3 |>
  dplyr::mutate(fraction = '3') |>
  bind_rows(bray_unifrac_eucl_tb_02)

palette_diversity <- c('genes' = "#270000",
                       "wunifrac_distance" = "#898989",
                       "euclidean_distance" = "#008241",
                       "bray_curtis_community" = "#a21e66", 
                       'bray_curtis_kmers' = '#FFCAFB')

##before correlation we check normality
shapiro.test(as.numeric(bray_unifrac_eucl_tb_02$genes)) # => p-value = 0.166 (NORMALITY)
ggqqplot(as.numeric(bray_unifrac_eucl_tb_02$genes))

shapiro.test(as.numeric(bray_unifrac_eucl_tb_w$wunifrac_distance)) # => p-value = 6.031e-08 (NO NORMALITY)
ggqqplot(as.numeric(bray_unifrac_eucl_tb_w$wunifrac_distance))

shapiro.test(as.numeric(bray_unifrac_eucl_tb_w$euclidean_distance)) # => p-value = 7.75e-10 (NO NORMALITY)
ggqqplot(as.numeric(bray_unifrac_eucl_tb_w$wunifrac_distance))

shapiro.test(as.numeric(bray_unifrac_eucl_tb_w$bray_curtis_community)) # => p-value = 7.144e-06 (NO NORMALITY)
ggqqplot(as.numeric(bray_unifrac_eucl_tb_w$bray_curtis_community))

#### first correlation: Bray-Curtis genes vs Bray-Curtis community, and Bray-Curtis genes vs. weighted UNIFRAC distance ----

corr_bray_02_plot <- bray_unifrac_eucl_tb_02 |>
  ggplot(aes(bray_curtis_community, genes))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed')+
  geom_point(color = "#a21e66", alpha = 0.8)+
  stat_cor(aes( #color = 'black', 
    
    label =   paste(..p.label..)), label.x = 0.35,
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#a21e66"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.35, label.y = 0.2, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#a21e66")+
  geom_point(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes), #shape = 2, 
             alpha = 0.8,  color = "#898989")+
  stat_cor(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes, 
    label =   paste(..p.label..)), label.x = 0.05,  
    label.y = 0.3,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#898989"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes, 
                                               label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.25,  color = "#898989" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(method = 'lm', color = "#a21e66", fill = "#a21e66")+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes), color = "#898989", fill = "#898989")+
  labs(x = 'wUNIFRAC and Bray-CurtisCommunity-Based', y = 'Bray-Curtis Genes Based')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 6), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_02_plot

## Environment Bray - Curtis community and Environment Weigthed UNIFRAC distance ---- 
corr_euclidean_bray_plot <- bray_unifrac_eucl_tb_w |>
  dplyr::mutate(euclidean_sc = scale(euclidean_distance)) |>
  ggplot(aes(euclidean_sc, bray_curtis_community))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed')+
  geom_point(color = "#a21e66", alpha = 0.8, aes(shape = fraction))+
  stat_cor(aes( group = fraction, 
    label =   paste(..p.label..)), label.x = 0.35,
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#a21e66"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes(group = fraction, label = paste0(..r.label..)),label.x = 0.35, label.y = 0.2, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#a21e66")+
  geom_point(data = bray_unifrac_eucl_tb_w |>
               dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_sc,  wunifrac_distance), #shape = 2, 
             alpha = 0.8,  color = "#898989")+
  stat_cor(data = bray_unifrac_eucl_tb_w |>
             dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_sc,  wunifrac_distance, 
                                               label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.3,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#898989"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_w |>
             dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_sc,  wunifrac_distance, 
                                                                          label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.25,  color = "#898989" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(aes(group = fraction), method = 'lm', color = "#a21e66", fill = "#a21e66")+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_w |>
                dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_sc,  wunifrac_distance, group = fraction), color = "#898989", fill = "#898989")+
  labs(x = 'Scaled Euclidean Distance', y = 'Bray-Curtis Communtiy Based & Weigthed UNIFRAC distance')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 10), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

corr_euclidean_bray_plot

corr_euclidean_bray_plot <- bray_unifrac_eucl_tb_w |>
  dplyr::mutate(euclidean_sc = scale(euclidean_distance)) |>
  ggplot(aes(euclidean_distance, bray_curtis_community))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed')+
  geom_point(color = "#a21e66", alpha = 0.8, aes(shape = fraction))+
  stat_cor(aes(group = fraction, 
               label =   paste(..p.label..)), 
           label.x = 7,
           label.y = 0.75,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#a21e66"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(aes(group = fraction, label = paste0(..r.label..)),
           label.y = 0.8, 
           label.x = 7, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#a21e66")+
  geom_point(data = bray_unifrac_eucl_tb_w |>
               dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance), #shape = 2, 
             alpha = 0.8,  color = "#898989")+
  stat_cor(data = bray_unifrac_eucl_tb_w |>
             dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
                                                                          label =   paste(..p.label..)), 
           label.x = 8,  
           label.y = 0.3,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#898989"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_w |>
             dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
                                                                          label =   paste(..r.label..)),
           label.x = 8, label.y = 0.35,  color = "#898989" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(aes(group = fraction), method = 'lm', color = "#a21e66", fill = "#a21e66")+
  scale_shape_manual(values = c(19,17), labels = labs_fraction)+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_w |>
                dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, group = fraction), color = "#898989", fill = "#898989")+
  labs(x = 'Euclidean Distance', y = 'Bray-Curtis Communtiy Based & wUNIFRAC distance',
       shape = 'Fraction')+
  guides(shape = guide_legend(ncol = 1))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 10), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

corr_euclidean_bray_plot

# ggsave( plot = corr_euclidean_bray_plot, filename = 'corr_euclidean_bray_plot.pdf',
#         path = 'results/figures/',
#         width = 180, height = 180, units = 'mm')

#### Correlation between Bray - Curtis dissimiliary and UNIFRAC weighted distance ----
corr_unifrac_bray_plot <- bray_unifrac_eucl_tb_w |>
  ggplot(aes(bray_curtis_community, wunifrac_distance))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed')+
  geom_point(color = "#a21e66", alpha = 0.8, aes(shape = fraction))+
  stat_cor(aes( group = fraction, 
                label =   paste(..p.label..)), label.x = 0.4,
           label.y = 0.35,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#a21e66"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(aes(label = paste0(..r.label..)),label.x = 0.4, label.y = 0.3, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#a21e66")+
  # geom_point(data = bray_unifrac_eucl_tb_w |>
  #              dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance), #shape = 2, 
  #            alpha = 0.8,  color = "#898989")+
  # stat_cor(data = bray_unifrac_eucl_tb_w |>
  #            dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
  #                                                                         label =   paste(..p.label..)), label.x = 0.05,  
  #          label.y = 0.3,
  #          p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#898989"#,
  #          #position = position_jitter(0.0)
  # )+
  # stat_cor(data = bray_unifrac_eucl_tb_w |>
  #            dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
  #                                                                         label =   paste(..r.label..)),
  #          label.x = 0.05, label.y = 0.25,  color = "#898989" ,
  #          p.digits = 0.01, digits = 2, 
  #          p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(aes(), method = 'lm', color = "#a21e66", fill = "#a21e66")+
  scale_shape_manual(values = c(19,17), labels = labs_fraction)+
  # geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_w |>
  #               dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, group = fraction), color = "#898989", fill = "#898989")+
  labs(x = 'Bray-Curtis Communtiy Based ', y = 'Weigthed UNIFRAC distance', shape = 'Fraction')+
  guides(shape = guide_legend(ncol = 1, title.position = 'top'))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

corr_unifrac_bray_plot

# ggsave( plot = corr_unifrac_bray_plot, filename = 'corr_unifrac_bray_plot.pdf',
#         path = 'results/figures/',
#         width = 180, height = 180, units = 'mm')

fraction_legend <- get_legend(corr_unifrac_bray_plot)

corr_unifrac_bray_plot <- bray_unifrac_eucl_tb_w |>
  ggplot(aes(bray_curtis_community, wunifrac_distance))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed')+
  geom_point(color = "#a21e66", alpha = 0.8, aes(shape = fraction))+
  stat_cor(aes( group = fraction, 
                label =   paste(..p.label..)), label.x = 0.4,
           label.y = 0.35,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#a21e66"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(aes(label = paste0(..r.label..)),label.x = 0.4, label.y = 0.3, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#a21e66")+
  # geom_point(data = bray_unifrac_eucl_tb_w |>
  #              dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance), #shape = 2, 
  #            alpha = 0.8,  color = "#898989")+
  # stat_cor(data = bray_unifrac_eucl_tb_w |>
  #            dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
  #                                                                         label =   paste(..p.label..)), label.x = 0.05,  
  #          label.y = 0.3,
  #          p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#898989"#,
  #          #position = position_jitter(0.0)
  # )+
  # stat_cor(data = bray_unifrac_eucl_tb_w |>
#            dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
#                                                                         label =   paste(..r.label..)),
#          label.x = 0.05, label.y = 0.25,  color = "#898989" ,
#          p.digits = 0.01, digits = 2, 
#          p.accuracy = 0.01, method = 'spearman')+
#geom_smooth(method = 'loess', color = 'grey')+
geom_smooth(aes(), method = 'lm', color = "#a21e66", fill = "#a21e66")+
  scale_shape_manual(values = c(19,17), labels = labs_fraction)+
  # geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_w |>
  #               dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, group = fraction), color = "#898989", fill = "#898989")+
  labs(x = 'Bray-Curtis Communtiy Based ', y = 'wUNIFRAC distance', shape = 'Fraction')+
  guides(shape = guide_legend(ncol = 1))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'none', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

corr_unifrac_bray_plot

## Create a combination of plots ----
# Arrange the second row plots first
second_row <- plot_grid(
  corr_bray_02_plot, corr_unifrac_bray_plot, fraction_legend,
  ncol = 3,  # Three columns in the second row
  rel_widths = c(1, 1, 0.15),
  labels = c('B', 'C') # Adjust relative widths, with the legend smaller
)

# Now arrange the full layout, with bray_unifrac_eucl_plot occupying the top row
relationship_community_genes_all <- plot_grid(
  bray_unifrac_eucl_plot,  # Top plot (spanning both columns)
  second_row,              # Second row with two plots
  ncol = 1,                # One column layout for the main grid
  rel_heights = c(2, 1),   # First plot 3 times the height of the second row
  labels = c('A')
)

# Print the final plot
print(relationship_community_genes_all)

# ggsave( plot = relationship_community_genes_all, 
#         filename = 'relationship_community_genes_all_v2.pdf',
#         path = 'results/figures/',
#         width = 180, height = 180, units = 'mm')


# Is there a change in Bray - Curtis of genes when there is a blooming event?  ----
  
  abund_asv_bloom_events_02_unique <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  left_join(tax_bbmo_10y_new) |>
  separate(sample_id, sep = '_', into = c('sample_id_red', 'fraction', 'seq_code'), remove = F) |>
  separate(sample_id_red, into = c('station', 'year', 'month', 'day'), sep = c(2, 4, 6), remove = F) |> 
  dplyr::mutate(date = paste0(year, '-', month, '-', day)) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3','asv5', 'asv8')) |>
  distinct(date) |>
  dplyr::mutate(date = as.Date(date, '%y-%m-%d')) |>
  dplyr::mutate(date = as.character(date),
                bloom_event = 'bloom')
  
corr_bray_02_plot <- bray_unifrac_eucl_tb_02 |>
  left_join(  abund_asv_bloom_events_02_unique) |>
  ggplot(aes(bray_curtis_community, genes))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed')+
  geom_point(color = "#a21e66", aes(alpha = ifelse(bloom_event == 'bloom', 1, 0.2)))+
  stat_cor(aes( #color = 'black', 
    label =   paste(..p.label..)), label.x = 0.35,
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#a21e66"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.35, label.y = 0.2, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#a21e66")+
  geom_point(data = bray_unifrac_eucl_tb_02 |>
               left_join(abund_asv_bloom_events_02_unique), 
             aes(wunifrac_distance, genes, alpha = ifelse(bloom_event == 'bloom', 1, 0.2)), #shape = 2, 
             color = "#898989")+
  stat_cor(data = bray_unifrac_eucl_tb_02 |>
             left_join(abund_asv_bloom_events_02_unique), aes(wunifrac_distance, genes, 
                                                                label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.3,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#898989"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes, 
                                               label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.25,  color = "#898989" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(method = 'lm', color = "#a21e66", fill = "#a21e66")+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes), color = "#898989", fill = "#898989")+
  labs(x = 'wUNIFRAC and Bray-CurtisCommunity-Based', y = 'Bray-Curtis Genes Based', 
       alpha = 'Bloom event')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 6), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_02_plot

