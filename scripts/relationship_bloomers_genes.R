# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO Dateseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##packages ----
library(cowplot)
library(tidyverse)
library(gridExtra)

# labels -----
labs_fraction_env <- as_labeller(c('0.2' = 'Free living\n(0.2-3 um)',
                                   '3' = 'Particle attached\n(3-20 um)',
                                   'env' = 'Environmental\nvariables'))

## palettes-----
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
bray_unifrac_eucl_tb <- read_csv('data/bray_unifrac_eucl_tb.csv') ## for plotting figure 7 
 
## There are differences in the codes for metagenomic data and metabarcoding data. 
## From now on in this part of the analysis I will numerate my sampples from 1 to 60 which are the number of samples that are shared between both datasets.
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

# upload data from metagenomes we have the table with counts for the genes related to phoyostynthesis----
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
  dplyr::reframe(anomalies_gene = get_anomalies(Date_lag = 2, negative = FALSE, na_rm = TRUE, 
                                                cutoff = 1,96, values = gene_counts, 
                                                plotting = FALSE)[c(1,2,3)])

z_scores_02_genes_photo <- genes_photo |>
  group_by(annot) |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  left_join(m_bbmo_10y_ed_4metag, by = c('sample_id_num')) |>
  #arrange(date) |>
  group_by(annot) |>
  dplyr::reframe(z_score_gen = get_anomalies(Date_lag = 2, negative = FALSE, cutoff = 1.96, 
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
  labs(y = 'Gene counts', x = 'Date (Y)', title = 'K02692, photosyntesis')+
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
  unique() #I have 60 samples in common with metaG data (remember that the days from one dataset and the other someDates do not match!)

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
  dplyr::reframe(anomalies_gene = get_anomalies(Date_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = gene_counts, plotting = FALSE)[c(1,2,3)])

z_scores_02_genes <- genes_cazy |>
  group_by(annot) |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  left_join(m_bbmo_10y_ed_4metag, by = c('sample_id_num')) |>
  #arrange(date) |>
  group_by(annot) |>
  dplyr::reframe(z_score_gen = get_anomalies(Date_lag = 2, negative = FALSE, cutoff = 1.96, 
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
  labs(x = 'Date (Y)', y = 'Gene counts', title = 'CAZYMES GH')+
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
  labs(x = 'Date (Y)', y = 'Gene counts', title = 'CAZYMES GT')+
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
  labs(x = 'Date (Y)', y = 'Gene counts', title = 'CAZYMES PL')+
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
  labs(x = 'Date (Y)', y = 'Gene counts', title = 'CAZYMES AA')+
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
  dplyr::select(annot, starts_with('BL'))  ## be sure that is the correcy normalization.

kegs_table |>
  dim()

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
  labs(y = 'Gene counts', x = 'Date (Y)', title = 'Photosystem II')+
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

# ggsave(genes_photosyntesis_bloo, 
#        path = 'results/figures/',
#        file = 'genes_phtosynyesis_bloo_together.pdf',
#        width = 180, height = 200, units = 'mm')

## whole Dateseries taxa that could present a phoyosystem II ----
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
# 
# ggsave(genes_photosyntesis_whole_community, 
#        path = 'results/figures/',
#        file = 'genes_photosyntesis_whole_community.pdf',
#        width = 180, height = 200, units = 'mm')

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
  labs(x = 'Date (Y)', y = 'Gene counts', title = 'CAZYMES GT')+
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
  labs(x = 'Date (Y)', y = 'Gene counts', title = 'CAZYMES GH')+
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

## kegs anomalies -----
kegs_table <- read_table('data/genes_sergio/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.tbl') |>
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

kegs_table <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.tbl', header = T) |>
  dplyr::select(annot, starts_with('BL'))  

kegs_anomalies <- kegs_table |>
  pivot_longer(cols = starts_with('BL')) |>
  group_by(annot) |>
  dplyr::filter(any(value > 2)) |>
  dplyr::reframe(anomalies_kos = get_anomalies(Date_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 20, values = value, plotting = FALSE)[c(1,2,3)])

kos_anom <- find_asv_with_anomalies(anomalies_result = kegs_anomalies, 
                                       anomaly_in1 = anomalies_kos, 
                                       anomaly_in2 = NULL,  #anomalies_ps,
                                       anomaly_in3 = NULL, #anomalies_clr, 
                                       logic1 = 'TRUE',
                                       logic2 = NULL, 
                                       logic3 = NULL,
                                       asv_col = annot)

kos_anom <- kos_anom |>
  as_tibble()

 kegs_anomalies_tb <- kegs_table |>
    dplyr::filter(annot %in% kos_anom$value) |>
    pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'gene_counts') |>
    dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
    dplyr::filter(sum(gene_counts) > 0) |>
    group_by(sample_id)  |>
    arrange(sample_id) |>
    mutate(sample_id_num = cur_group_id()) |>
    left_join(m_bbmo_10y_ed_4metag_red, by = c('sample_id_num')) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
   dplyr::filter(annot != 'sum_not_annotated')
 
 kegs_anomalies_tb$annot |>
   unique()
 
 plot_kos <- kegs_anomalies_tb |>
   ggplot(aes(date, gene_counts))+
   #geom_point()+
   geom_line(aes(group = annot))+
   labs(y = 'KOs', x = '')+
   theme_bw()+
   theme(panel.grid = element_blank(),
         axis.text = element_text(size = 5))
 
 plot_kos
 
## cazy anomalies ----
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
 
 genes_cazy |>
   group_by(annot) |>
   dplyr::reframe(mean = mean(gene_counts)) 
 
cazy_anomalies <-  genes_cazy |>
   group_by(annot) |>
   dplyr::filter(any(gene_counts > 5)) |>
   dplyr::reframe(anomalies_cazy = get_anomalies(Date_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 7, values =  gene_counts, plotting = FALSE)[c(1,2,3)])
 
cazy_anom <- find_asv_with_anomalies(anomalies_result = cazy_anomalies, 
                                     anomaly_in1 = anomalies_cazy, 
                                     anomaly_in2 = NULL,  #anomalies_ps,
                                     anomaly_in3 = NULL, #anomalies_clr, 
                                     logic1 = 'TRUE',
                                     logic2 = NULL, 
                                     logic3 = NULL,
                                     asv_col = annot)
 
cazy_anom  <- cazy_anom  |>
   as_tibble()
 
cazy_anomalies_tb <-  genes_cazy |>
   dplyr::filter(annot %in% cazy_anom$value)  |>
  dplyr::filter(annot != 'sum_not_annotated')

cazy_anomalies_tb$annot |>
  unique()
 
 plot_cazy <- cazy_anomalies_tb |>
   ggplot(aes(date, gene_counts))+
   #geom_point()+
   #scale_y_sqrt()+
   geom_line(aes(group = annot))+
   labs(y = 'CAZymes counts', x = '')+
   theme_bw()+
   theme(panel.grid = element_blank(),
         axis.text.x = element_text(size = 5))
 plot_cazy
 
## pfam anomalies -----
 genes_pfam <- read_table('data/genes_sergio/BBMOSOLA-GC_250bp_pfam.lengthNorm.metaGsizeNorm.counts.tbl') |>
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
 
 pfam_mean <- genes_pfam |>
   group_by(annot) |>
   dplyr::reframe(mean = mean(gene_counts)) 
 
 pfam_mean |>
   arrange(mean)
 
 pfam_anomalies <-  genes_pfam |>
   group_by(annot) |>
   dplyr::filter(any(gene_counts > 2000)) |>
   dplyr::reframe(anomalies_pfam = get_anomalies(Date_lag = 2, negative = FALSE, na_rm = TRUE,
                                                 cutoff = 500, values =  gene_counts, plotting = FALSE)[c(1,2,3)])
 
 pfam_anom <- find_asv_with_anomalies(anomalies_result = pfam_anomalies, 
                                      anomaly_in1 = anomalies_pfam, 
                                      anomaly_in2 = NULL,  #anomalies_ps,
                                      anomaly_in3 = NULL, #anomalies_clr, 
                                      logic1 = 'TRUE',
                                      logic2 = NULL, 
                                      logic3 = NULL,
                                      asv_col = annot)
 
 pfam_anom  <- pfam_anom  |>
   as_tibble()
 
 pfam_anomalies_tb <-  genes_pfam |>
   dplyr::filter(annot %in% pfam_anom$value)
 
 pfam_anomalies_tb$annot |>
   unique()
 
 plot_pfam <- pfam_anomalies_tb |>
   dplyr::filter(annot != 'sum_not_annotated') |>
   ggplot(aes(date, gene_counts))+
   #geom_point()+
   geom_line(aes(group = annot))+
   labs(y = 'pfam', x = '')+
   theme_bw()+
   theme(panel.grid = element_blank(),
         axis.text = element_text(size = 5))
 
 plot_pfam
 
## eggnog anomalies -----
 genes_eggnogg <- read_table('data/genes_sergio/BBMOSOLA-GC_250bp_eggNOG.lengthNorm.metaGsizeNorm.counts.tbl') |>
   dplyr::select(annot, starts_with('BL')) |> #filter only Blanes samples not SOLA 
   as_tibble() |>
   pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'gene_counts') |>
   dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
   dplyr::filter(sum(gene_counts) > 1000) |>
   group_by(sample_id)  |>
   arrange(sample_id) |>
   mutate(sample_id_num = cur_group_id()) |>
   left_join(m_bbmo_10y_ed_4metag, by = c('sample_id_num')) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
   dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013'))
 
eggnogg_mean <-  genes_eggnogg |>
   group_by(annot) |>
   dplyr::reframe(mean = mean(gene_counts)) 
 
eggnogg_mean |>
   arrange(mean)
 
 eggnogg_anomalies <-  genes_eggnogg |>
   group_by(annot) |>
   dplyr::filter(any(gene_counts > 2000)) |>
   dplyr::reframe(anomalies_eggnogg = get_anomalies(Date_lag = 2, negative = FALSE, na_rm = TRUE, 
                                                    cutoff = 500, values =  gene_counts, plotting = FALSE)[c(1,2,3)])
 
 eggnogg_anom <- find_asv_with_anomalies(anomalies_result = eggnogg_anomalies, 
                                      anomaly_in1 = anomalies_eggnogg, 
                                      anomaly_in2 = NULL,  #anomalies_ps,
                                      anomaly_in3 = NULL, #anomalies_clr, 
                                      logic1 = 'TRUE',
                                      logic2 = NULL, 
                                      logic3 = NULL,
                                      asv_col = annot)
 
 eggnogg_anom  <- eggnogg_anom  |>
   as_tibble()

 eggnogg_anomalies_tb <-  genes_eggnogg |>
   dplyr::filter(annot %in% eggnogg_anom$value)  
 
 eggnogg_anomalies_tb$annot |>
   unique()
 
 plot_eggnogg <- eggnogg_anomalies_tb |>
   dplyr::filter(annot != 'sum_not_annotated') |>
   ggplot(aes(date, gene_counts))+
   #geom_point()+
   geom_line(aes(group = annot))+
   labs(y = 'eggnogg', x = '')+
   theme_bw()+
   theme(panel.grid = element_blank(),
         axis.text = element_text(size = 5))
 
 plot_eggnogg
 
### 
  plot_asvs_genes <-  test |>
    ggplot(aes(date, abundance_value, color = order))+
    #geom_point()+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y = 'Relative Abundance')+
  geom_line(aes(group = asv_num))+
  theme_bw()+
  theme(legend.position = 'none', panel.grid = element_blank())

plot_asvs_genes

genes_asvs_blooms <- plot_grid(plot_asvs_genes, 
          plot_kos,
          plot_cazy,
          plot_pfam,
          plot_eggnogg,
          labels = c('C', 'D', 'E', 'F'),
          label_fontface = 'plain',
          label_fontfamily = 7,
          ncol = 1)

genes_asvs_blooms

## geom lile ------
plot_eggnogg <- eggnogg_anomalies_tb |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::mutate(gene_counts = scale(gene_counts)) |>
  ggplot(aes(as.character(date), annot))+
  scale_fill_gradient()+
  #geom_point()+
  geom_tile(aes(fill = gene_counts ))+
  labs(y = 'eggnogg', x = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(size = 0))

plot_eggnogg

plot_kegs <- kegs_anomalies_tb |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::mutate(gene_counts = scale(gene_counts)) |>
  ggplot(aes(as.character(date), annot))+
  scale_fill_gradient()+
  #geom_point()+
  geom_tile(aes(fill = gene_counts ))+
  labs(y = 'eggnogg', x = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(size = 0))

plot_kegs

plot_cazy <- cazy_anomalies_tb |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::mutate(gene_counts = scale(gene_counts)) |>
  ggplot(aes(as.character(date), annot))+
  scale_fill_gradient()+
  #geom_point()+
  geom_tile(aes(fill = gene_counts ))+
  labs(y = 'eggnogg', x = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(size = 0))

plot_cazy

plot_pfam <- pfam_anomalies_tb |>
  dplyr::filter(annot != 'sum_not_annotated') |>
  dplyr::mutate(gene_counts = scale(gene_counts)) |>
  ggplot(aes(as.character(date), annot))+
  scale_fill_gradient()+
  #geom_point()+
  geom_tile(aes(fill = gene_counts ))+
  labs(y = 'eggnogg', x = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(size = 0))

plot_pfam

## pheatmap clustered do blooming events cluster together?
# Example matrix
## highlight events of bloom ---
## plot blooms during 2009-2014 FL ------
asv_tab_0914 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv15', 'asv27', 'asv58'))

data_heatmap_annot_c <- bloom_event |>
  dplyr::filter(asv_num %in% c('asv11', 'asv15', 'asv27', 'asv58', 'asv1', 'asv7')) |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(str_detect(date, '2009|2010|2011|2012|2013')) |>
  dplyr::filter(bloom_event == 'bloom') |>
  # dplyr::group_by(date) |>
  # dplyr::reframe(n = n()) |>
  # dplyr::filter(n >4)
  pivot_wider(names_from = 'asv_num', values_from = 'bloom_event', values_fill = 'no-bloom') |>
  dplyr::select(-fraction) |>
  column_to_rownames(var = "date") 

# add CLR abundance of my bloomers 
asv_tab_0914 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv15', 'asv27', 'asv58', 'asv1', 'asv7')) |>
  dplyr::select(abundance_value, asv_num, date) |>
  dplyr::mutate(date = as.character(as.Date(date))) |>
  pivot_wider(names_from = 'asv_num', values_from = 'abundance_value')

annotation_colors <- list(
  bloom_event = c("bloom" = "red", "no-bloom" = "gray")
)

pfam_anomalies_tb_l |>
  row.names()

pfam_anomalies_tb_l <- pfam_anomalies_tb |>
  dplyr::select(date, gene_counts, annot) |>
  pivot_wider(id_cols = date, names_from = annot, values_from = gene_counts, values_fill = 0) |>
  as.data.frame()

pfam_anomalies_tb_l <- pfam_anomalies_tb_l |>
  dplyr::mutate(date = as.character(as.Date(date))) |>
  column_to_rownames(var = 'date') |>
  as.matrix() 

pfam_anomalies_tb_l_num <- apply(pfam_anomalies_tb_l, 2, as.numeric)
rownames(pfam_anomalies_tb_l_num) <- rownames(pfam_anomalies_tb_l)

annotation_colors <- list(
  asv15 = c("bloom" = "red", "no-bloom" = "white"),
  asv27 = c("bloom" = "red", "no-bloom" = "white"),
  asv58 = c("bloom" = "red", "no-bloom" = "white"),
  asv11 = c("bloom" = "red", "no-bloom" = "white"),
  asv1 =  c("bloom" = "red", "no-bloom" = "white"),
  asv7 = c("bloom" = "red", "no-bloom" = "white")
)

# Create heatmap with clustering
pheatmap(pfam_anomalies_tb_l, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "complete", display_numbers = F,
         annotation_row = data_heatmap_annot_c ,
         annotation_colors = annotation_colors,
         cluster_cols = T)
## keggs
kegs_anomalies_tb_l <- kegs_anomalies_tb |>
  dplyr::select(date, gene_counts, annot) |>
  dplyr::mutate(gene_counts = scale(gene_counts)) |>
  pivot_wider(id_cols = date, names_from = annot, values_from = gene_counts, values_fill = 0) |>
  as.data.frame()

kegs_anomalies_tb_l <- kegs_anomalies_tb_l |>
  dplyr::mutate(date = as.character(as.Date(date))) |>
  column_to_rownames(var = 'date') |>
  as.matrix() 

kegs_anomalies_tb_l_num <- apply(kegs_anomalies_tb_l, 2, as.numeric)
rownames(kegs_anomalies_tb_l_num) <- rownames(kegs_anomalies_tb_l)

annotation_colors <- list(
  asv15 = c("bloom" = "red", "no-bloom" = "white"),
  asv27 = c("bloom" = "red", "no-bloom" = "white"),
  asv58 = c("bloom" = "red", "no-bloom" = "white"),
  asv11 = c("bloom" = "red", "no-bloom" = "white"),
  asv1 =  c("bloom" = "red", "no-bloom" = "white"),
  asv7 = c("bloom" = "red", "no-bloom" = "white")
)

# Create heatmap with clustering
pheatmap(kegs_anomalies_tb_l, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "complete", display_numbers = F,
         annotation_row = data_heatmap_annot_c ,
         #annotation_colors = annotation_colors,
         cluster_cols = T)

enteros <- tax_bbmo_10y_new |>
  dplyr::filter(order == 'Enterobacterales') |>
  distinct(family, genus)

## cazy 
cazy_anomalies_tb_l <- cazy_anomalies_tb |>
  dplyr::select(date, gene_counts, annot) |>
  dplyr::mutate(gene_counts = sqrt(gene_counts)) |>
  pivot_wider(id_cols = date, names_from = annot, values_from = gene_counts, values_fill = 0) |>
  as.data.frame()

cazy_anomalies_tb_l <- cazy_anomalies_tb_l |>
  dplyr::mutate(date = as.character(as.Date(date))) |>
  column_to_rownames(var = 'date') |>
  as.matrix() 

cazy_anomalies_tb_l_num <- apply(cazy_anomalies_tb_l, 2, as.numeric)
rownames(cazy_anomalies_tb_l_num) <- rownames(cazy_anomalies_tb_l)

annotation_colors <- list(
  asv15 = c("bloom" = "red", "no-bloom" = "white"),
  asv27 = c("bloom" = "red", "no-bloom" = "white"),
  asv58 = c("bloom" = "red", "no-bloom" = "white"),
  asv11 = c("bloom" = "red", "no-bloom" = "white")
)

# Create heatmap with clustering
pheatmap(cazy_anomalies_tb_l, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "complete", display_numbers = F,
         annotation_row = data_heatmap_annot_c ,
         annotation_colors = annotation_colors,
         cluster_cols = T)


## eggnogg
eggnogg_anomalies_tb_l <- eggnogg_anomalies_tb |>
  dplyr::select(date, gene_counts, annot) |>
  dplyr::mutate(gene_counts = scale(gene_counts)) |>
  pivot_wider(id_cols = date, names_from = annot, values_from = gene_counts, values_fill = 0) |>
  as.data.frame()

eggnogg_anomalies_tb_l <- eggnogg_anomalies_tb_l |>
  dplyr::mutate(date = as.character(as.Date(date))) |>
  column_to_rownames(var = 'date') |>
  as.matrix() 

eggnogg_anomalies_tb_l_num <- apply(eggnogg_anomalies_tb_l, 2, as.numeric)
rownames(eggnogg_anomalies_tb_l_num) <- rownames(eggnogg_anomalies_tb_l)

annotation_colors <- list(
  asv15 = c("bloom" = "red", "no-bloom" = "white"),
  asv27 = c("bloom" = "red", "no-bloom" = "white"),
  asv58 = c("bloom" = "red", "no-bloom" = "white"),
  asv11 = c("bloom" = "red", "no-bloom" = "white")
)

# Create heatmap with clustering
pheatmap(eggnogg_anomalies_tb_l, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "complete", display_numbers = F,
         annotation_row = data_heatmap_annot_c ,
         annotation_colors = annotation_colors,
         cluster_cols = T)

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
                       "wunifrac_distance" = "'#FFBAA6'",
                       "euclidean_distance" = "#008241",
                       "bray_curtis_community" = "#a21e66", 
                       'bray_curtis_kmers' = '#FFCAFB')

labs_diversity <- as_labeller(c('genes' = 'Bray Curtis Genes',
                                "wunifrac_distance" = 'Weigthed\nUNIFRAC',
                                "euclidean_distance" = 'Euclidean Environmental\nDistance',
                                "bray_curtis_community" = 'Bray Curtis\nCommunity Composition',
                                'bray_curtis_kmers' = 'Bray Curtis k-mers'))

labs_fraction_env <- as_labeller(c('0.2' = 'Free living\n(0.2-3 um)',
                                   '3' = 'Particle attached\n(3-20 um)',
                                   'env' = 'Environmental\nvariables'))

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

bray_curtis_genes_all_red <- bray_curtis_genes_filtered_v5 |>
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
  bind_rows(bray_curtis_m_tb_red)# |>
  #bind_rows(bc_ramiro_red)

#write.csv(bray_unifrac_eucl_tb, 'data/bray_unifrac_eucl_tb.csv', row.names = F)

## Bray Curtis KOs SCG normalized (the correct dataset) -----
filtered_genes <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.tbl') ## threshold 0.001

filtered_genes_names <- filtered_genes %$%
  V1

filtered_genes_t <- t(filtered_genes)

filtered_genes_t <- filtered_genes_t[-1,]

filtered_genes_t |>
  colnames() <- filtered_genes_names

filtered_genes_t <- filtered_genes_t |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t))) |>
  dplyr::filter(str_detect(annot, 'BL'))

filtered_genes_t_l <- filtered_genes_t |>
  dplyr::select(-sum_not_annotated) |>
  pivot_longer(starts_with('K'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_genes_filtered <- dissimilarity_matrix(data = filtered_genes_t_l, 
                                                   sample_id_col = sample_id,
                                                   values_cols_prefix = 'BL')

bray_curtis_keggs <- bray_curtis_genes_filtered |>
  rename(keggs_scg = bray_curtis_result)

## EGGNOGG SCG ----
filtered_genes <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_pfam34.0.lengthNorm.SCGnorm.counts.tbl') ## threshold 0.001

filtered_genes_names <- filtered_genes %$%
  V1

filtered_genes_t <- t(filtered_genes)

filtered_genes_t <- filtered_genes_t[-1,]

filtered_genes_t |>
  colnames() <- filtered_genes_names

filtered_genes_t <- filtered_genes_t |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t))) |>
  dplyr::filter(str_detect(annot, 'BL'))

filtered_genes_t_l <- filtered_genes_t |>
  dplyr::select(-sum_not_annotated) |>
  pivot_longer(starts_with('P'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_genes_filtered <- dissimilarity_matrix(data = filtered_genes_t_l, 
                                                   sample_id_col = sample_id,
                                                   values_cols_prefix = 'BL')

bray_curtis_pfam <- bray_curtis_genes_filtered |>
  rename(pfam_scg = bray_curtis_result)

## eggnogg ----
filtered_genes <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_eggNOG.lengthNorm.SCGnorm.counts.tbl') ## threshold 0.001

filtered_genes |>
  dim()

filtered_genes_names <- filtered_genes %$%
  V1

filtered_genes_t <- t(filtered_genes)

filtered_genes_t <- filtered_genes_t[-1,]

filtered_genes_t |>
  colnames() <- filtered_genes_names

filtered_genes_t <- filtered_genes_t  |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t))) |>
  dplyr::filter(str_detect(annot, 'BL'))

filtered_genes_t_l <- filtered_genes_t |>
  dplyr::select(-sum_not_annotated, -starts_with('EN')) |>
  pivot_longer(starts_with('C'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_genes_filtered <- dissimilarity_matrix(data = filtered_genes_t_l, 
                                                   sample_id_col = sample_id,
                                                   values_cols_prefix = 'BL')

bray_curtis_eggnogg <- bray_curtis_genes_filtered |>
  rename(eggnogg_scg = bray_curtis_result)

## CAZY genes sgc (already correct) ---- 
filtered_genes <- read_table('data/genes_sergio/BBMOSOLA-GC_250bp_CAZy.lengthNorm.SCGnorm.counts.tbl')

filtered_genes_names <- filtered_genes %$%
  V1

filtered_genes_t <- t(filtered_genes)

filtered_genes_t |>
  colnames() <- filtered_genes_t[1,]

filtered_genes_t <- filtered_genes_t[-1,]

filtered_genes_t <- filtered_genes_t |>
  as.data.frame() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(filtered_genes_t))) |>
  rownames_to_column(var = 'annot') |>
  dplyr::filter(str_detect(annot, 'BL'))

filtered_genes_t |>
  colnames()

filtered_genes_t_l <- filtered_genes_t |>
  dplyr::select(-sum_not_annotated) |>
  pivot_longer(cols = !c('sample_id', 'annot'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_genes_filtered <- dissimilarity_matrix(data = filtered_genes_t_l, 
                                                   sample_id_col = sample_id,
                                                   values_cols_prefix = 'BL')

bray_curtis_cazy <- bray_curtis_genes_filtered |>
  rename(cazy_scg = bray_curtis_result)

## merge the tb with bc estimated with scg data ---- 
bc_scg <- bray_curtis_pfam |>
  left_join(bray_curtis_eggnogg) |>
  left_join(bray_curtis_keggs) |>
  left_join(bray_curtis_cazy)

m_bbmo_10y_ed_4metag_red2 <- m_bbmo_10y_ed_4metag_red |>
  dplyr::mutate(sample_id_ed = paste0('BL', 1:nrow(m_bbmo_10y_ed_4metag)))

bc_scg <- bc_scg |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('samples' = 'sample_id_ed' )) |>
  dplyr::select(date, pfam_scg, eggnogg_scg, keggs_scg, cazy_scg )

bc_scg$date

bc_scg <- bc_scg |>
  dplyr::mutate(date = as.character(as.Date(date)))

## 
bc_scg |>
  ggplot(aes(row_index_2, as.numeric(pfam_scg)))+
  geom_line()+
  geom_line(aes(row_index_2, as.numeric(eggnogg_scg)))+
  geom_line(aes(row_index_2, as.numeric(keggs_scg)))+
  geom_line(aes(row_index_2, as.numeric(cazy_scg)))+
  scale_x_continuous(breaks = c(11, 11+12, 11+24, 11+12*3, 11+12*4, 11+12*5, 11+12*6))

bray_unifrac_eucl_tb <- bray_unifrac_eucl_tb_02  |>
  full_join(©)

bray_unifrac_eucl_tb |>
  dim()

## fold change during bloom event, is there a shift at the level of functions? -----
keggs_fold_change <- kegs_table |>
  dplyr::select(annot, sample_id = sample_id.x, sample_id_num, gene_counts) |>
  dplyr::mutate(gene_counts = if_else(gene_counts == 0, 1, gene_counts)) |>
  pivot_wider(names_from = 'annot', values_from = 'gene_counts') |>
  dplyr::mutate(across(starts_with("K"), ~ . / lag(.), .names = "fold_change_{.col}")) |>
  dplyr::select( sample_id, starts_with('fold_change'))

keggs_fold_change |>
  dim()

keggs_fold_change_mean <- keggs_fold_change |>
  pivot_longer(cols = starts_with('fold_change')) |>
  group_by(sample_id) |>
  dplyr::filter(value > 1) |>
  dplyr::reframe(mean_fold = mean(value, na.rm = T)) |>
  dplyr::mutate(n = 2:nrow(keggs_fold_change))

keggs_fold_change |>
  pivot_longer(cols = starts_with('fold_change')) |>
  group_by(sample_id) |>
  dplyr::filter(sample_id %in% c('BL110112')) |>
  ungroup() |>
  slice_max(order_by = value, n = 25) |>
  separate(name, into = c('fold_change', 'KO'), sep = 'e_') |>
  dplyr::select(KO) %$%
  cat(KO)

kegs_fold_plot <- keggs_fold_change_mean |>
  ggplot(aes(n, mean_fold))+
  geom_line()

kegs_fold_plot


asv_tab_0914 |>
  colnames()

blooms_0914_plot <- asv_tab_0914 |>
 # dplyr::select(-asv1, -asv7) |>
  #dplyr::mutate(n = 1:nrow(keggs_fold_change)) |>
  #pivot_longer(cols = starts_with('asv'), values_to = 'value') |>
  #left_join(tax_bbmo_10y_new, by = c('name' = 'asv_num')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, abundance_value))+
  geom_vline(xintercept = heterotrophic_blooms$date)+
  geom_area(aes(group = asv_num, fill = order_f))+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(interaction(family_f,genus_f,asv_num)))+
  labs(x = 'Date', y = 'rCLR')+
  theme_bw()+
  theme(legend.position = 'none', strip.background = element_rect(fill = NA))

blooms_0914_plot

asv_tab_0914 |>
  colnames()

plot_grid(kegs_fold_plot,
          blooms_0914_plot,
          ncol = 1)

corr_bray_genes_community_tb |>
  ggplot(aes(KEGG, KEGG_sgc))+
  geom_point()+
  geom_smooth(method = 'lm')

corr_bray_genes_community_tb |>
  ggplot(aes(pfam,pfam_scg))+
  geom_point()+
  geom_smooth(method = 'lm')

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

### both fractions at the same Date ----
bray_unifrac_eucl_tb_02 <- bray_unifrac_eucl_tb_02 |>
  dplyr::mutate(fraction = '0.2') |>
  dplyr::mutate(date = as.character(date))

bray_unifrac_eucl_tb_w <- bray_unifrac_eucl_tb_3 |>
  dplyr::mutate(fraction = '3') |>
  bind_rows(bray_unifrac_eucl_tb_02)

palette_diversity <- c('genes' = "#270000",
                       "wunifrac_distance" = "'#FFBAA6'",
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
             alpha = 0.8,  color = "'#FFBAA6'")+
  stat_cor(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes, 
    label =   paste(..p.label..)), label.x = 0.05,  
    label.y = 0.3,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "'#FFBAA6'"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes, 
                                               label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.25,  color = "'#FFBAA6'" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(method = 'lm', color = "#a21e66", fill = "#a21e66")+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes), color = "'#FFBAA6'", fill = "'#FFBAA6'")+
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

## Construction of new correlation plot for main text ----
# write.csv(corr_bray_genes_community_tb, 'data/corr_bray_genes_community_tb_scg.csv', row.names = F) 
# write.csv(bray_unifrac_eucl_tb_02, 'data/bray_unifrac_eucl_tb_02.csv', row.names = F) 

bray_unifrac_eucl_tb_02 <- read.csv('data/bray_unifrac_eucl_tb_02.csv')
corr_bray_genes_community_tb <- read.csv('data/corr_bray_genes_community_tb_scg.csv')

# corr_bray_genes_community_tb <- corr_bray_genes_community_tb |>
#   left_join(bc_scg) |>
#   dplyr::filter(!is.na(date))

# bray_unifrac_eucl_tb_02_scg <- bray_unifrac_eucl_tb_02
# 
# bray_unifrac_eucl_tb_02 |>
#   colnames()
# 
# corr_bray_genes_community_tb |>
#   colnames()
# 
# bray_curtis_genes_filtered |>
#   dim()
# 
# bray_curtis_genes_filtered <- bray_curtis_genes_filtered[1:61,]
# bray_curtis_genes_filtered <- bray_curtis_genes_filtered |>
#   dplyr::rename(KEGG_sgc = bray_curtis_result)
# 
# corr_bray_genes_community_tb |>
#   dim()
# 
# corr_bray_genes_community_tb <- corr_bray_genes_community_tb |>
#   bind_cols(bray_curtis_genes_filtered)

#### first correlation: Bray-Curtis genes vs Bray-Curtis community, and Bray-Curtis genes vs. weighted UNIFRAC distance ----
corr_bray_02_plot <- bray_unifrac_eucl_tb_02 |>
  ggplot(aes(bray_curtis_community, genes))+
  scale_x_continuous(limits = c(0.05, 0.85))+
  scale_y_continuous(limits = c(0.05, 0.85))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 1, alpha = 0.5)+
  geom_point(color = "#00808F", alpha = 0.5)+
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
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_02_plot

corr_unifrac_02_plot <- bray_unifrac_eucl_tb_02 |>
  ggplot(aes(bray_curtis_community, genes))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 1, alpha = 0.5)+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_point(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes), #shape = 2, 
             alpha = 0.5,  color = "#00808F")+
  scale_x_continuous(limits = c(0.05, 0.85))+
  scale_y_continuous(limits = c(0.05, 0.85))+
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

# Now arrange the full layout, with bray_unifrac_eucl_plot occupying the top row
# corr_bray_unifrac_plot <- plot_grid(
#   corr_bray_02_plot, corr_unifrac_02_plot,
#   # Top plot (spanning both columns)          # Second row with two plots
#   ncol = 2,                # One column layout for the main grid
#   #rel_heights = c(1, 2.75, 1, 0.25),   # First plot 3 Dates the height of the second row
#   labels = c('A', 'B')
# )
# 
# # Print the final plot
# print(corr_bray_unifrac_plot)

# ggsave( plot = corr_bray_unifrac_plot,
#         filename = 'corr_bray_unifrac_plot.pdf',
#         path = 'results/figures/',
#         width = 180, height = 100, units = 'mm')

## I add KOs and genes in the same correlation plot. KOs will inform us if there is a change in the functionality -----
corr_bray_02_kos_plot <- corr_bray_genes_community_tb |>
  ggplot(aes(bray_curtis_community, KEGG_sgc))+
  scale_x_continuous(limits = c(0.05, 0.85))+
  scale_y_continuous(limits = c(0.05, 0.3))+
  #geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 1, alpha = 0.5)+
  geom_point( #shape = 2, 
             alpha = 0.5,  color = "#00808F")+
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
    panel.border = element_blank(),
    strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_02_kos_plot

corr_unifrac_02_ko_plot <- corr_bray_genes_community_tb |>
  ggplot(aes(wunifrac_distance,KEGG_sgc))+
  #geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 1, alpha = 0.5)+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_point(data = corr_bray_genes_community_tb, aes(wunifrac_distance ,KEGG), #shape = 2, 
             alpha = 0.5,  color = "#00808F")+
  scale_x_continuous(limits = c(0.05, 0.85))+
  scale_y_continuous(limits = c(0.05, 0.3))+
  stat_cor(data = corr_bray_genes_community_tb, aes(wunifrac_distance, KEGG, 
                                               label =   paste(..p.label..)), label.x = 0.6,  
           label.y = 0.25,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#00808F"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = corr_bray_genes_community_tb, aes(wunifrac_distance, KEGG, 
                                               label =   paste(..r.label..)),
           label.x = 0.6, label.y = 0.2,  color = "#00808F" ,
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

# corr_unifrac_02_plot_legend <- corr_bray_genes_community_tb |> 
#   pivot_longer(cols = !c('date', 'fraction')) |> 
#   ggplot(aes(x = date, y =value, shape = name)) + 
#   geom_point(color = "#00808F") + 
#   scale_shape_manual(values = c('genes' = 16, 'KEGG' = 17),
#                      labels = c('genes' = 'Genes', 'KEGG' = 'KOs - KEGGs')) + 
#   labs(shape = 'Bray-Curtis based') +
#   guides(shape = guide_legend(ncol = 2, title.position = 'left'))+
#   theme_minimal() 
#   
# legend_shape_genes_kos <- get_legend(corr_unifrac_02_plot_legend)
 
# Now arrange the full layout, with bray_unifrac_eucl_plot occupying the top row
corr_bray_unifrac_kos_plot <- plot_grid(
  corr_bray_02_plot, corr_unifrac_02_plot,
  corr_bray_02_kos_plot, corr_unifrac_02_ko_plot,
  # Top plot (spanning both columns)          # Second row with two plots
  ncol = 2,                # One column layout for the main grid
  rel_heights = c(1, 1),   # First plot 3 Dates the height of the second row
  labels = c('A', 'B', 'C', 'D'),
  label_fontface = 'plain'
)

# Print the final plot
print(corr_bray_unifrac_kos_plot)

# ggsave( plot = corr_bray_unifrac_kos_plot,
#         filename = 'corr_bray_unifrac_plot_v5.pdf',
#         path = 'results/figures/',
#         width = 180, height = 180, units = 'mm')

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
             alpha = 0.8,  color = "'#FFBAA6'")+
  stat_cor(data = bray_unifrac_eucl_tb_w |>
             dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_sc,  wunifrac_distance, 
                                               label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.3,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "'#FFBAA6'"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_w |>
             dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_sc,  wunifrac_distance, 
                                                                          label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.25,  color = "'#FFBAA6'" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(aes(group = fraction), method = 'lm', color = "#a21e66", fill = "#a21e66")+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_w |>
                dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_sc,  wunifrac_distance, group = fraction), color = "'#FFBAA6'", fill = "'#FFBAA6'")+
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
             alpha = 0.8,  color = "'#FFBAA6'")+
  stat_cor(data = bray_unifrac_eucl_tb_w |>
             dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
                                                                          label =   paste(..p.label..)), 
           label.x = 8,  
           label.y = 0.3,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "'#FFBAA6'"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_w |>
             dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
                                                                          label =   paste(..r.label..)),
           label.x = 8, label.y = 0.35,  color = "'#FFBAA6'" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(aes(group = fraction), method = 'lm', color = "#a21e66", fill = "#a21e66")+
  scale_shape_manual(values = c(19,17), labels = labs_fraction)+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_w |>
                dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, group = fraction), color = "'#FFBAA6'", fill = "'#FFBAA6'")+
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
  #            alpha = 0.8,  color = "'#FFBAA6'")+
  # stat_cor(data = bray_unifrac_eucl_tb_w |>
  #            dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
  #                                                                         label =   paste(..p.label..)), label.x = 0.05,  
  #          label.y = 0.3,
  #          p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "'#FFBAA6'"#,
  #          #position = position_jitter(0.0)
  # )+
  # stat_cor(data = bray_unifrac_eucl_tb_w |>
  #            dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
  #                                                                         label =   paste(..r.label..)),
  #          label.x = 0.05, label.y = 0.25,  color = "'#FFBAA6'" ,
  #          p.digits = 0.01, digits = 2, 
  #          p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(aes(), method = 'lm', color = "#a21e66", fill = "#a21e66")+
  scale_shape_manual(values = c(19,17), labels = labs_fraction)+
  # geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_w |>
  #               dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, group = fraction), color = "'#FFBAA6'", fill = "'#FFBAA6'")+
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
  #            alpha = 0.8,  color = "'#FFBAA6'")+
  # stat_cor(data = bray_unifrac_eucl_tb_w |>
  #            dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
  #                                                                         label =   paste(..p.label..)), label.x = 0.05,  
  #          label.y = 0.3,
  #          p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "'#FFBAA6'"#,
  #          #position = position_jitter(0.0)
  # )+
  # stat_cor(data = bray_unifrac_eucl_tb_w |>
#            dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, 
#                                                                         label =   paste(..r.label..)),
#          label.x = 0.05, label.y = 0.25,  color = "'#FFBAA6'" ,
#          p.digits = 0.01, digits = 2, 
#          p.accuracy = 0.01, method = 'spearman')+
#geom_smooth(method = 'loess', color = 'grey')+
geom_smooth(aes(), method = 'lm', color = "#a21e66", fill = "#a21e66")+
  scale_shape_manual(values = c(19,17), labels = labs_fraction)+
  # geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_w |>
  #               dplyr::mutate(euclidean_sc = scale(euclidean_distance)), aes(euclidean_distance,  wunifrac_distance, group = fraction), color = "'#FFBAA6'", fill = "'#FFBAA6'")+
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
  rel_heights = c(2, 1),   # First plot 3 Dates the height of the second row
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
             color = "'#FFBAA6'")+
  stat_cor(data = bray_unifrac_eucl_tb_02 |>
             left_join(abund_asv_bloom_events_02_unique), aes(wunifrac_distance, genes, 
                                                                label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.3,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "'#FFBAA6'"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes, 
                                               label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.25,  color = "'#FFBAA6'" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(method = 'lm', color = "#a21e66", fill = "#a21e66")+
  geom_smooth(method = 'lm', data = bray_unifrac_eucl_tb_02, aes(wunifrac_distance, genes), color = "'#FFBAA6'", fill = "'#FFBAA6'")+
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

## ---- BLOOMS and CARD-FISH (CLARA RUIZ GONZÁLEZ DATA) ---- ##
abund_asv_bloom_events_02_unique <- asv_tab_all_bloo_z_tax_summary_all |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  left_join(tax_bbmo_10y_new) |>
  distinct(asv_num)

abund_asv_bloo_fl_08_10 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'rclr' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value, date) |>
  dplyr::filter(year %in% c('2008','2009', '2010')) |>
  left_join(tax_bbmo_10y_new)

abund_asv_bloo_fl_08_10 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value, date) |>
  dplyr::filter(year %in% c('2008','2009', '2010')) |>
  left_join(tax_bbmo_10y_new)

asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num %in% bloo_02_filt$value) |>
  dplyr::filter(date %in% abund_asv_bloo_fl_08_10$date) |>
  dplyr::mutate(family_asv_num = paste0(family, asv_num)) |>
  dplyr::filter(fraction == '0.2') |>
  ggplot(aes(date, abundance_value))+
  geom_line(aes(group = asv_num_f))+
  facet_wrap(vars(family_asv_num))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'))


## New combination of plots: recurrence in Blanes, Bray Curtis and wUNIFRAC over the years + Env distance (IF WE KEEP IT THEN I MOVE IT TO THE SUMMARY SCRIPT) ----

palette_diversity <- c('genes' = "#270000",
                       "wunifrac_distance" = "'#FFBAA6'",
                       "euclidean_distance" = "#008241",
                       "bray_curtis_community" = "#a21e66", 
                       'bray_curtis_kmers' = '#FFCAFB')

palette_fraction_env <- c("env" = "#008241", 
                          '0.2' = '#00808F', 
                          '3' = "#454545")

# plot weighted unifrac distances -----
bray_unifrac_eucl_tb$bray_curtis_type |>
  unique()

data <- bray_unifrac_eucl_tb |>
  dplyr::filter(!bray_curtis_type %in% c('bray_curtis_kmers', "genes")) 

data$bray_curtis_type <- factor(data$bray_curtis_type, levels = c('bray_curtis_community', 'wunifrac_distance', 'euclidean_distance'))
data$fraction |>
  unique()

data$fraction <- factor(data$fraction, levels = c('0.2', '3', 'env'))

# legend 
legend_plot <- data |>
  #dplyr::filter(fraction != 'env') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  ggplot(aes(date, bray_curtis_result))+
  facet_wrap(fraction ~ bray_curtis_type, scales = "free_x", ncol = 1, 
             labeller = labeller(fraction = function(x) ifelse(x == "fraction", "", x), 
                                 bray_curtis_type = labs_diversity), 
             strip.position = 'left') + 
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  #scale_linetype_discrete(labels = labs_diversity)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
  scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(data = data |>
              dplyr::filter(fraction == '3'), aes(date, bray_curtis_result), alpha = 0.3)+
  labs(x = 'Date', y = '', color = '', linetype = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  #guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    strip.text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    legend.text = element_text(size = 7),
    panel.border = element_blank(),
    plot.margin = margin(t = 5, b = 5))
legend_plot 

legend_plot <- get_legend(legend_plot)

bray_unifrac_eucl_plot <- data |>
  dplyr::filter(fraction != 'env') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  ggplot(aes(date, bray_curtis_result))+
  facet_wrap(fraction ~ bray_curtis_type, scales = "fixed", ncol = 1, 
             labeller = labeller(fraction = function(x) ifelse(x %in% c("0.2", '3'), "", x), 
                                 bray_curtis_type = labs_diversity), 
             strip.position = 'left') + 
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  #scale_linetype_discrete(labels = labs_diversity)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
  scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(data = data |>
              dplyr::filter(fraction == '3'), aes(date, bray_curtis_result), alpha = 0.3)+
  labs(x = 'Date', y = '', color = '', linetype = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  #guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)
  #scale_y_continuous(limits = c(0.15, 1.0))+
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
    axis.ticks.length.y =  unit(0.2, "mm"))

bray_unifrac_eucl_plot

eucl_plot <- data |>
  dplyr::filter(fraction == 'env',
                bray_curtis_type == 'euclidean_distance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  ggplot(aes(date, bray_curtis_result))+
  facet_wrap(fraction ~ bray_curtis_type, scales = "free_x", ncol = 1, 
             labeller = labeller(fraction = function(x) ifelse(x == "env", "", x), 
                                 bray_curtis_type = labs_diversity), 
             strip.position = 'left') + 
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  #scale_linetype_discrete(labels = labs_diversity)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
  scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Date', y = '', color = '', linetype = '')+
  scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  #guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)
  #scale_y_continuous(limits = c(0.15, 1.0))+
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
    plot.margin = margin(t = 5, l = 10))

eucl_plot

bray_curtis_lag_plot <- bray_curtis_lag_all |>
  #dplyr::filter(as.numeric(lag) < 108) |>
  ggplot(aes(lag, bray_curtis_result))+
  scale_x_continuous(breaks = c(12, 24, 36, 48, 12*5, 12*6, 12*7, 12*8, 12*9, 120), labels = c(1, 2, 3,4, 5, 6, 7, 8,  9, 10))+
  geom_point(aes(shape = fraction), alpha = 0.02)+
  # geom_point(data = bray_curtis_lag_all |>
  #              group_by(fraction, lag) |>
  #              dplyr::reframe(mean = mean(bray_curtis_result)), aes(lag, mean, group = fraction, color = fraction, fill =  fraction))+
  #facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1)+
  scale_shape_discrete(labels = labs_fraction)+
  labs(y = 'Bray Curtis Dissimilarity', x = 'Lag between samples', shape = '', color = '', fill = '', linetype = '')+
  geom_line(data = bray_curtis_lag_all |>
              group_by(fraction, lag) |>
              dplyr::reframe(mean = mean(bray_curtis_result)), aes(lag, mean, group = fraction, color = fraction, 
                                                                   linetype = fraction), linewidth = 1)+

  #geom_smooth(aes(group = fraction, color = fraction, fill =  fraction), method = 'loess', span = 0.1)+
  #geom_smooth(aes(group = fraction,  color = fraction, fill =  fraction), method = 'lm')+
  scale_linetype_discrete( labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
  guides(shape = guide_legend(ncol = 2))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        strip.background = element_rect(fill = 'transparent'),
        legend.position = 'none',
        plot.margin = margin(l = 45, t = 5, b = 5))

bray_curtis_lag_plot 

# Now arrange the full layout, with bray_unifrac_eucl_plot occupying the top row
BBMO_community_diversity_presentation_plot <- plot_grid(
  bray_curtis_lag_plot,
  bray_unifrac_eucl_plot,
  eucl_plot,
  legend_plot,
  # Top plot (spanning both columns)          # Second row with two plots
  ncol = 1,                # One column layout for the main grid
  rel_heights = c(1, 2.75, 1, 0.25),   # First plot 3 Dates the height of the second row
  labels = c('A', 'B', 'c')
)

# Print the final plot
print(BBMO_community_diversity_presentation_plot)

# ggsave( plot = BBMO_community_diversity_presentation_plot,
#         filename = 'BBMO_community_diversity_presentation_plot.pdf',
#         path = 'results/figures/',
#         width = 180, height = 200, units = 'mm')

asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num == 'asv11') |>
  ggplot(aes(date, abundance_value))+
  geom_line(aes(group = fraction))+
  facet_wrap(vars(fraction))+
  theme_bw()

## I calculate the Bray Curtis for the KOs dataset and see if they also correlate with my Bray Curtis community based ----
### I explore the other genes tables (eggNOG, CAZY, KEGGs, and PFAM)
#### upload data
cazy_tb <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_CAZy.lengthNorm.SCGnorm.counts.tbl') 
egg_tb <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_eggNOG.lengthNorm.metaGsizeNorm.counts.tbl') 
kegg_tb <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.metaGsizeNorm.counts.txt') 
pfam_tb <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_pfam.lengthNorm.metaGsizeNorm.counts.tbl') 

### cazy_tb Bray Curtis  ----
colnames(cazy_tb) <- cazy_tb[1,]
cazy_tb <- cazy_tb[-1,]

cazy_tb <- cazy_tb |>
  as_tibble() |>
  dplyr::select(annot, matches('^BL09|^BL10|^BL11|^BL12|^BL13'))## only Blanes not SOLA

cazy_tb |>
  dim()

cazy_tb_l <- cazy_tb |>
  pivot_longer(starts_with('BL'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::select(asv_num = annot, sample_id = asv_num, relative_abundance) |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_cazy <- dissimilarity_matrix(data = cazy_tb_l, 
                                                   sample_id_col = sample_id,
                                                   values_cols_prefix = 'BL')

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(bray_curtis_cazy)))

bray_curtis_cazy <-  bray_curtis_cazy |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(bray_curtis_cazy))) |>
  left_join(abund_asv_bloo_fl_09_13, by = c( 'sample_id')) |>
  dplyr::mutate(bray_type = 'cazymes')

### egg_tb Bray Curtis ----
colnames(egg_tb) <- egg_tb[1,]
egg_tb <- egg_tb[-1,]

egg_tb <- egg_tb |>
  as_tibble() |>
  dplyr::select(annot, matches('^BL09|^BL10|^BL11|^BL12|^BL13'))## only Blanes not SOLA

egg_tb |>
  dim() #61 

egg_tb_l <- egg_tb |>
  pivot_longer(starts_with('BL'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::select(asv_num = annot, sample_id = asv_num, relative_abundance) |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_egg <- dissimilarity_matrix(data = egg_tb_l, 
                                         sample_id_col = sample_id,
                                         values_cols_prefix = 'BL')

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(bray_curtis_egg)))

bray_curtis_egg <-  bray_curtis_egg |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(bray_curtis_egg))) |>
  left_join(abund_asv_bloo_fl_09_13, by = c( 'sample_id')) |>
  dplyr::mutate(bray_type = 'eggNOG')

### kegg_tb Bray Curtis ----
colnames(kegg_tb) <- kegg_tb[1,]
kegg_tb <- kegg_tb[-1,]

kegg_tb <- kegg_tb |>
  as_tibble() |>
  dplyr::select(annot, matches('^BL09|^BL10|^BL11|^BL12|^BL13'))## only Blanes not SOLA

kegg_tb |>
  dim()

# kegg_tb |>
#   dplyr::mutate(sample_id = paste0('BL', 1:nrow(kegg_tb_t)))

kegg_tb_l <- kegg_tb |>
  pivot_longer(starts_with('BL'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::select(asv_num = annot, sample_id = asv_num, relative_abundance) |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_kegg <- dissimilarity_matrix(data = kegg_tb_l, 
                                         sample_id_col = sample_id,
                                         values_cols_prefix = 'BL')

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(bray_curtis_kegg)))

bray_curtis_kegg <-  bray_curtis_kegg |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(bray_curtis_kegg))) |>
  left_join(abund_asv_bloo_fl_09_13, by = c( 'sample_id')) |>
  dplyr::mutate(bray_type = 'KEGG')

### pfam_tb Bray Curtis 
colnames(pfam_tb) <- pfam_tb[1,]
pfam_tb <- pfam_tb[-1,]

pfam_tb <- pfam_tb |>
  as_tibble() |>
  dplyr::select(annot, matches('^BL09|^BL10|^BL11|^BL12|^BL13'))## only Blanes not SOLA

pfam_tb |>
  dim() #61

pfam_tb_l <- pfam_tb |>
  pivot_longer(starts_with('BL'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::select(asv_num = annot, sample_id = asv_num, relative_abundance) |>
  dplyr::mutate(relative_abundance = as.numeric(relative_abundance))

bray_curtis_pfam <- dissimilarity_matrix(data = pfam_tb_l, 
                                         sample_id_col = sample_id,
                                         values_cols_prefix = 'BL')

## add dates 
abund_asv_bloo_fl_09_13 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  distinct(year, date) |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(bray_curtis_pfam)))

bray_curtis_pfam <-  bray_curtis_pfam |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(bray_curtis_pfam))) |>
  left_join(abund_asv_bloo_fl_09_13, by = c( 'sample_id')) |>
  dplyr::mutate(bray_type = 'pfam')

## plot all the different genetic datasets -----
bray_curtis_genes_filtered_v5 <- bray_curtis_genes_filtered_v5 |>
  dplyr::mutate(bray_type = 'genes')

bray_curtis_genes_all_types_tb <- bray_curtis_cazy |>
  bind_rows(bray_curtis_egg) |>
  bind_rows(bray_curtis_kegg) |>
  bind_rows(bray_curtis_pfam) |>
  bind_rows(bray_curtis_genes_filtered_v5)

bray_curtis_genes_all_types_tb %$%
  unique(bray_type)

## nice palette for the different genes content 
palette_bray_types <- c('#A7FFE6',
                        '#FFBAA6',
                       # '#FFBAA6',
                        '#AA928B',
                        '#8BAA90')

bray_genes_all <- bray_curtis_genes_all_types_tb |>
  ggplot(aes(date, bray_curtis_result))+
  geom_line(aes(group = bray_type, color = bray_type), linewidth = 0.75)+
  scale_color_manual(values = palette_bray_types)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = 'Bray Curtis Dissimilarity', color = '', linetype = '')+
  guides(shape = 'none',
         color = guide_legend(ncol = 2))+
  #guide_legend(keywidth = unit(0.5, "cm"))+
  guides(color = guide_legend(ncol = 5))+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank())

bray_genes_all

bray_02_community <- bray_unifrac_eucl_tb_02 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, bray_curtis_community))+
  geom_line(aes( ), linewidth = 0.75, color = '#00808F')+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = 'Bray Curtis Dissimilarity', color = '', linetype = '')+
  #guide_legend(keywidth = unit(0.5, "cm"), ncol = 2)
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    panel.border = element_blank())

bray_02_community 

bray_community_vs_genes <- plot_grid(bray_genes_all,
          bray_02_community,rel_heights = c(2,1),
          ncol = 1,
          labels = c('A', 'B'), label_fontface = 'plain',
          label_size = 9,
          hjust = 0.1
          )

# ggsave( plot = bray_community_vs_genes, filename = 'bray_community_vs_genes.pdf',
#         path = 'results/figures/relationship_genes_blooms/',
#         width = 180, height = 140, units = 'mm')

## Correlations between the different genes correlations and the community ----
bray_curtis_genes_filtered_v5_ed <- bray_curtis_genes_filtered_v5 |>
  dplyr::select(samples, row_index_2, bray_curtis_result,  year, date, bray_type) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(year = as.double(year)) |>
  dplyr::mutate(sample_id = 'sample_id.x') |>
  dplyr::select(samples, row_index_2, bray_curtis_result, sample_id, year, date, bray_type)

bray_curtis_genes_all_types_tb <- bray_curtis_genes_filtered_v5_ed |>
  bind_rows(bray_curtis_pfam) |>
  bind_rows(bray_curtis_egg) |>
  bind_rows(bray_curtis_kegg) |>
  bind_rows(bray_curtis_cazy)

bray_unifrac_eucl_tb_02_ed <- bray_unifrac_eucl_tb_02 |> 
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(-genes) |>
  dplyr::mutate(date = as.character(date)) 
 
corr_bray_genes_community_tb <- bray_curtis_genes_all_types_tb |>
  dplyr::select(-sample_id, -year, -row_index_2, -samples) |>
  pivot_wider(names_from = 'bray_type', values_from = 'bray_curtis_result') |>
  dplyr::mutate(date = as.Date(date)) |>
  dplyr::mutate(date = as.character(date)) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(bray_unifrac_eucl_tb_02_ed, by = c('date'))

corr_bray_genes_community_tb  |>
  colnames()

##before correlation we check normality ----
shapiro.test(as.numeric(corr_bray_genes_community_tb$cazymes)) # => p-value = 1.537e-06 (NO NORMALITY)
ggqqplot(as.numeric(corr_bray_genes_community_tb$cazymes))

shapiro.test(as.numeric(corr_bray_genes_community_tb$eggNOG)) # => p-value = 7.583e-06(NO NORMALITY)
ggqqplot(as.numeric(corr_bray_genes_community_tb$eggNOG))

shapiro.test(as.numeric(corr_bray_genes_community_tb$KEGG)) # => p-value = 4.927e-06 (NO NORMALITY)
ggqqplot(as.numeric(corr_bray_genes_community_tb$KEGG))

shapiro.test(as.numeric(corr_bray_genes_community_tb$pfam)) # => p-value = 1.287e-07 (NO NORMALITY)
ggqqplot(as.numeric(corr_bray_genes_community_tb$pfam))

shapiro.test(as.numeric(corr_bray_genes_community_tb$bray_curtis_result)) # => p-value = 0.166 (NORMALITY)
ggqqplot(as.numeric(corr_bray_genes_community_tb$bray_curtis_result))

shapiro.test(as.numeric(corr_bray_genes_community_tb$genes)) # => p-value = 0.166 (NORMALITY)
ggqqplot(as.numeric(corr_bray_genes_community_tb$genes))

shapiro.test(as.numeric(corr_bray_genes_community_tb$bray_curtis_community)) # =>p-value = 0.1728 (NORMALITY)
ggqqplot(as.numeric(corr_bray_genes_community_tb$bray_curtis_community))

#### plot the correlations ----
palette_bray_types <- c('#A7FFE6',
                        '#2D515E',
                        '#FFBAA6',
                        '#AA928B',
                        '#8BAA90') #,'#2C362E'

corr_bray_genes_types_plot <- corr_bray_genes_community_tb |>
  ggplot(aes(bray_curtis_community, cazymes))+
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed')+
  geom_point(color = "#A7FFE6", alpha = 0.8)+
  stat_cor(aes( #color = 'black', 
    label =   paste(..p.label..)), label.x = 0.05,
    label.y = 0.22,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "#A7FFE6"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.05, label.y = 0.20, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman',
           color = "#A7FFE6")+
  geom_smooth(method = 'lm', color = "#A7FFE6", fill = "#A7FFE6")+
  
  geom_point(data = corr_bray_genes_community_tb, aes(bray_curtis_community, eggNOG), #shape = 2, 
             alpha = 0.8,  color = '#2D515E')+
  stat_cor(data = corr_bray_genes_community_tb, aes(bray_curtis_community, eggNOG, 
                                               label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.35,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = '#2D515E'
           #position = position_jitter(0.0)
  )+
  stat_cor(data = corr_bray_genes_community_tb, aes(bray_curtis_community, eggNOG,
                                               label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.32,  color = '#2D515E' ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  geom_smooth(method = 'lm', data = corr_bray_genes_community_tb, aes(bray_curtis_community, eggNOG), color = '#2D515E', fill = '#2D515E')+ 
  
  geom_point(data = corr_bray_genes_community_tb, aes(bray_curtis_community, KEGG), #shape = 2, 
             alpha = 0.8,  color = '#AA928B')+
  stat_cor(data = corr_bray_genes_community_tb, aes(bray_curtis_community, KEGG, 
                                                    label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.24,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = '#AA928B'
           #position = position_jitter(0.0)
  )+
  stat_cor(data = corr_bray_genes_community_tb, aes(bray_curtis_community, KEGG,
                                                    label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.26,  color = '#AA928B' ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  geom_smooth(method = 'lm', data = corr_bray_genes_community_tb, aes(bray_curtis_community, KEGG), color = '#AA928B', fill = '#AA928B')+ 

  geom_point(data = corr_bray_genes_community_tb, aes(bray_curtis_community, pfam), #shape = 2, 
             alpha = 0.8,  color = '#8BAA90')+
  stat_cor(data = corr_bray_genes_community_tb, aes(bray_curtis_community, pfam, 
                                                    label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.3,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = '#8BAA90'
           #position = position_jitter(0.0)
  )+
  stat_cor(data = corr_bray_genes_community_tb, aes(bray_curtis_community, pfam,
                                                    label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.28,  color = '#8BAA90' ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  geom_smooth(method = 'lm', data = corr_bray_genes_community_tb, aes(bray_curtis_community, pfam), color = '#8BAA90', fill = '#8BAA90')+ 

  geom_point(data = corr_bray_genes_community_tb, aes(bray_curtis_community, genes), #shape = 2, 
             alpha = 0.8,  color = '#FFBAA6')+
  stat_cor(data = corr_bray_genes_community_tb, aes(bray_curtis_community, genes, 
                                                    label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.4,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = '#FFBAA6'
           #position = position_jitter(0.0)
  )+
  stat_cor(data = corr_bray_genes_community_tb, aes(bray_curtis_community, genes,
                                                    label =   paste(..r.label..)),
           label.x = 0.05, label.y = 0.38,  color = '#FFBAA6' ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  geom_smooth(method = 'lm', data = corr_bray_genes_community_tb, aes(bray_curtis_community, genes), color = '#FFBAA6', fill = '#FFBAA6')+
  
  labs(x = 'Bray-Curtis Community-Based', y = 'Bray-Curtis Genes Based')+
  theme_bw()+
  theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

corr_bray_genes_types_plot

bray_genes_all_legend <- get_legend(bray_genes_all)

corr_bray_community_types_genes <- plot_grid(
  corr_bray_genes_types_plot,
  bray_genes_all_legend,
                        rel_heights = c(2.25,0.5),
                                     ncol = 1,
                                    # labels = c('A'), label_fontface = 'plain',
                                     label_size = 9,
                                     hjust = 0.1
)
# 
# ggsave( plot = corr_bray_community_types_genes, filename = 'corr_bray_community_types_genes.pdf',
#         path = 'results/figures/relationship_genes_blooms/',
#         width = 180, height = 180, units = 'mm')

## genes from Adrià ------
genes_amylibacter <- read.csv('data/genes_sergio/genes_asv27.csv')

genes_amylibacter <- genes_amylibacter |>
  dplyr::select(count, sample, annotation) |>
  left_join(m_02, by = c('sample' = 'code')) |>
  dplyr::select(date, count, sample, annotation) |>
  dplyr::mutate(date = as.character(as.Date(date))) 

genes_amylibacter |>
  left_join(asv27)

m_02_ed <- m_02 |>
  dplyr::mutate(date = as.character(as.Date(date))) |>
  dplyr::select(date,code ) |>
  dplyr::mutate(date = ())

asv27 <- asv27  |>
  left_join(m_02_ed, by = c('date'))

m_02_ed$date
asv27$date

setdiff(m_02_ed$date, asv27$date)
setdiff(asv27$date, tibble1$date) 

#days are different between both datasets 

genes_amylibacter$date <- gsub("-..$", "", genes_amylibacter$date)
asv27$date <- gsub("-..$", "", asv27$date)

genes_amylibacter$date

genes_amylibacter |>
  left_join(asv27) |>
  dplyr::filter(asv27 > 0) |>
  ggplot(aes(count, asv27))+
  geom_point()+
  labs(y = 'Amylibacter')+
  geom_smooth(method = 'lm')

genes_glaciecola <- read.csv('data/genes_sergio/genes_asv11.csv') |>
  dplyr::select(count, sample, annotation, species) |>
  left_join(m_02, by = c('sample' = 'code')) |>
  dplyr::select(date, count, sample, annotation, species) |>
  dplyr::mutate(date = as.character(as.Date(date))) 

asv11 <- asv_tab_0914 |>
  dplyr::select(date, asv11)

asv11$date <- gsub("-..$", "", asv11$date)
genes_glaciecola$date <- gsub("-..$", "", genes_glaciecola$date)

genes_glaciecola |>
  dplyr::filter(species == 'Glaciecola sp000155775') |>
  left_join(asv11) |>
  dplyr::filter(asv11 > 0) |>
  ggplot(aes(log10(count), asv11))+
  geom_point()+
  labs(y = 'Glaciecola sp000155775')+
  geom_smooth(aes(group = annotation), method = 'lm')+
  facet_wrap(annotation~species)

asv11_rel <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv11', 'asv15', 'asv27', 'asv58')) |>
  dplyr::filter(asv_num == 'asv11') |>
  dplyr::select(date, abundance_value) |>
  dplyr::mutate(date = as.character(as.Date(date))) 

asv11_rel$date <- gsub("-..$", "", asv11_rel$date)
asv11

genes_glaciecola |>
  dplyr::filter(species == 'Glaciecola sp000155775') |>
  left_join(asv11_rel) |>
  dplyr::filter(abundance_value > 0) |>
  ggplot(aes(log10(count), log10(abundance_value)))+
  geom_point()+
  labs(y = 'Glaciecola sp000155775')+
  geom_smooth(aes(group = annotation), method = 'lm')+
  facet_wrap(annotation~species)

read.csv('data/genes_sergio/genes_asv11.csv') |>
  dplyr::select( annotation, cogs, cycle, gene_name, gene.long.name)  |>
  distinct()

## -------- KOs Abundance - Occurrence (ubiquity % of samples) --------
kos_tb <- read.table('data/genes_sergio/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.tbl') ## threshold 0.001

kos_tb_names <- kos_tb %$%
  V1

kos_tb_t <- t(kos_tb)

kos_tb_t <- kos_tb_t[-1,]

kos_tb_t |>
  colnames() <- kos_tb_names

kos_tb_t <- kos_tb_t |>
  as_tibble() |>
  dplyr::mutate(sample_id = paste0('BL', 1:nrow(kos_tb_t))) |>
  dplyr::filter(str_detect(annot, 'BL'))

kos_tb_t_l <- kos_tb_t |>
  dplyr::select(-sum_not_annotated) |>
  pivot_longer(starts_with('K'), values_to = 'relative_abundance', names_to = 'asv_num') |>
  dplyr::mutate(reads = as.numeric(relative_abundance)) |>
  dplyr::select(-relative_abundance) |>
  calculate_rel_abund(group_cols = annot) |>
  dplyr::filter(str_detect(annot, 'BL09|BL10|BL11|BL12|BL13')) |> #filter by samples that match
  rename(KO = asv_num)

kos_tb_t_l |>
  group_by(annot) |>
  dplyr::reframe(sum = sum(relative_abundance)) |>
  dplyr::filter(sum != 1) # check that they all sum 1

kos_tb_t_l |>
  arrange(-relative_abundance)

n_ocurrence <- kos_tb_t_l |>
  dplyr::filter(relative_abundance > 0) |>
  group_by(KO, annot)  |>
  dplyr::reframe(n = n()) |>
  dplyr::group_by(KO) |>
  dplyr::reframe(occurrence = sum(n)/60)

n_ocurrence |>
  dplyr::filter(is.na(occurrence))

kos_ocurrence_abund <- kos_tb_t_l |>
  dplyr::group_by(KO) |>
  dplyr::reframe(mean_abund = mean(relative_abundance, na.rm = T),
                 sd_abund = sd(relative_abundance, na.rm = T),
                 max = max(relative_abundance)) |>
  left_join(n_ocurrence) |>
  dplyr::filter(!is.na(occurrence)) # their abundance is always 0

kos_ocurrence_abund |>
  dplyr::filter(is.na(occurrence)) |>
  dplyr::filter(mean_abund == 0 &
                  max == 0)

kos_ocurrence_abund |>
  arrange(-sd_abund)

kos_ocurrence_abund |>
  ggplot(aes(log10(mean_abund), log10(occurrence)))+
  geom_point(aes(alpha = if_else(sd_abund > 2, 1, 0.1)))+
  geom_smooth(method = 'loess', span = 0.5, color = 'black')+
  theme_bw()

## identify outliers
threshold_residual = 0.0001
threshold_occurrence = log10(0.9)
threshold_mean_abund = log10(0.00001)
threshold_max = 0.002

kos_ocurrence_abund |>
  dplyr::mutate(variables = 'KO') |>
  ggplot(aes(occurrence, variables))+
  geom_density_ridges()

kos_ocurrence_abund |>
  dplyr::mutate(variables = 'KO') |>
  ggplot(aes(max, variables))+
  geom_density_ridges()+
  geom_vline(xintercept = 0.002)+
  scale_x_sqrt()

kos_ocurrence_abund |>
  dplyr::mutate(variables = 'KO') |>
  ggplot(aes(sqrt(mean_abund), variables))+
  geom_density_ridges()+
  geom_vline(xintercept = 0.02)

kos_ocurrence_abund %$%
  range(mean_abund)

kos_tb_t_l %$%
  range(relative_abundance)

kos_ocurrence_abund <- kos_ocurrence_abund %>%
  dplyr::filter(!is.na(occurrence)) %>%  # Remove rows with NA residuals
  # Calculate residuals and filter based on conditions
  dplyr::mutate(residuals = log10(occurrence) - predict(loess(log10(occurrence) ~ log10(mean_abund), data = kos_ocurrence_abund))) 

kos_ocurrence_abund |>
  dplyr::mutate(variables = 'KO') |>
  ggplot(aes(sqrt(residuals), variables))+
  geom_density_ridges()+
  geom_vline(xintercept = sqrt(0.015))

threshold_residual = 0.08
threshold_occurrence = log10(1/3)
threshold_mean_abund = log10(0.01)
threshold_max = 0.0000000001

kos_ocurrence_abund |>
  arrange(-mean_abund)

kos_ocurrence_abund |>
  arrange(-residuals)

kos_ocurrence_abund_res <- kos_ocurrence_abund |>
  dplyr::filter(!is.na(occurrence)) %>%  # Remove rows with NA residuals
  dplyr::mutate(residuals = log10(occurrence) - predict(loess(log10(occurrence) ~ log10(mean_abund), data = kos_ocurrence_abund))) |>
  dplyr::mutate(residual_cat = case_when((residuals > threshold_residual &
                                            log10(mean_abund) < threshold_mean_abund &
                                            log10(occurrence) < threshold_occurrence &
                                            max > threshold_max) ~ 'outlier',
                                         TRUE ~ 'default')) |>
  dplyr::filter(residual_cat == 'outlier')

kos_ocurrence_abund_res  |>
  # Calculate residuals and filter based on conditions
  ggplot(aes(mean_abund, occurrence)) +
  #geom_abline(slope = 0.8, intercept = log10(0.000001)) +
  #geom_smooth(method = 'loess', span = 0.5, color = 'black') +  
  theme_bw() +  
  # Calculate residuals and highlight outliers
  geom_point(aes(color = residual_cat), alpha = 0.2,
    size = 3) +  # Size of points
  scale_color_manual(values = c('outlier' = 'orange', 'default' = 'blue')) +
  facet_wrap(vars(residual_cat))+
  theme(legend.position = "none")  # Hide legend if needed

kos_ocurrence_abund %>%
  # Calculate residuals and filter based on conditions
  dplyr::mutate(residuals = log10(occurrence) - predict(loess(log10(occurrence) ~ log10(mean_abund), data = kos_ocurrence_abund))) %>%
  ggplot(aes(log10(mean_abund), log10(occurrence))) +
  #geom_abline(slope = 0.8, intercept = log10(0.000001)) +
  geom_smooth(method = 'loess', span = 0.5, color = 'black') +  
  theme_bw() +  
  # Calculate residuals and highlight outliers
  geom_point(aes(color = if_else((
    -residuals > threshold_residual &
      log10(mean_abund) < threshold_mean_abund &
      log10(occurrence) < threshold_occurrence &
      max > threshold_max ), 
    'outlier', 'default')), 
    size = 3, alpha = 0.8) +  # Size of points
  scale_color_manual(values = c('outlier' = 'orange', 'default' = 'blue')) +
  theme(legend.position = "none")  # Hide legend if needed


bloomers_kos_top25_oc_tb <- kos_ocurrence_abund |>
  dplyr::filter(KO %in% top25_kos_blooming_events$KO) 

bloomers_kos_oc_tb <- kos_ocurrence_abund |>
  dplyr::filter(KO %in% exclusive_kos$KO) 

kos_anomalies_excl_bloo <- kos_ocurrence_abund |>
  dplyr::filter(KO %in% kegs_anomalies_tb$annot) |>
  dplyr::filter(KO %in% exclusive_kos$KO)   

kos_ocurrence_abund %>%
  # Calculate residuals and filter based on conditions
  dplyr::mutate(residuals = log10(occurrence) - predict(loess(log10(occurrence) ~ log10(max), data = kos_ocurrence_abund))) %>%
  ggplot(aes(log10(max), log10(occurrence))) +
  #geom_abline(slope = 0.8, intercept = log10(0.000001)) +
  geom_smooth(method = 'loess', span = 0.5, color = 'black') +  
  theme_bw() +  
  geom_point(aes(color = if_else((
    -residuals > threshold_residual &
      log10(mean_abund) < threshold_mean_abund &
      log10(occurrence) < threshold_occurrence &
      max > threshold_max ),
    'outlier', 'default')),
    size = 3, alpha = 0.2) +  # Size of points
  scale_color_manual(values = c('outlier' = 'orange', 'default' = 'blue')) +
  geom_point(data = kos_anomalies_excl_bloo , aes(log10(mean_abund), log10(occurrence)), color = 'black')+
  theme(legend.position = "none")  # Hide legend if needed

kos_rare <- kos_ocurrence_abund_res %>%
  dplyr::filter((residuals > threshold_residual &
                   log10(mean_abund) < threshold_mean_abund &
                   log10(occurrence) < threshold_occurrence &
                   max > threshold_max)) |>
  distinct(KO)
   
kos_rare_tb <- kos_tb_t_l |>
  dplyr::filter(KO %in% kos_rare$KO) |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_ed')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013'))
  
plot_rare_kos <- kos_rare_tb |>
  ggplot(aes(date, relative_abundance))+
  geom_point()+
  geom_line(aes(group = KO))+
  scale_y_log10()+
  #geom_boxplot(aes(group = date))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4))

plot_rare_kos

kos_abundant_tb <- kos_tb_t_l |>
  dplyr::filter(!KO %in% kos_rare$KO) |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_ed')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 

plot_abundant_kos <-  kos_abundant_tb |>
  dplyr::filter(relative_abundance != 0) |>
  ggplot(aes(date, relative_abundance))+
  geom_area(aes(fill = KO), position = 'stack')+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  #geom_line(aes(group = KO))+
  #scale_y_log10()+
  #geom_boxplot(aes(group = date))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'none')

plot_abundant_kos

kos_rare_tb$asv_num |>
  unique()

kos_rare_tb |>
  dplyr::filter(asv_num %in% c(unique(kegs_anomalies_tb$annot)))

plot_rare_kos <-   kos_rare_tb |>
  dplyr::filter(relative_abundance != 0) |>
  ggplot(aes(date, relative_abundance))+
  geom_point()+
  geom_vline(xintercept = heterotrophic_blooms$date)+
  #geom_line(aes(group = asv_num))+
  #scale_y_log10()+
  geom_boxplot(aes(group = date))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4))

plot_rare_kos

plot_rare_kos <-  kos_rare_tb |>
  dplyr::filter(relative_abundance != 0) |>
  ggplot(aes(date, relative_abundance))+
  geom_area(aes(fill = asv_num), position = 'stack')+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  #geom_line(aes(group = asv_num))+
  #scale_y_log10()+
  #geom_boxplot(aes(group = date))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'none')

plot_rare_kos


# ------------------ ## --------- Blooming events and their functional implications ------------------ ## ---------
## upload data
asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') |> ##using dada2 classifier assign tax with silva 138.1 and correctly identifying bloomers
  as_tibble() |>
  dplyr::select(-X)

## -------- Are there consistent trends between heterotrophic blooms and genetic content from metagenomes? --------
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

### Bloom Glaciecola 
genes_bloo_tb <- read.csv('data/genes_sergio/BBMO-GC_250bp_gene.lengthNorm.SingleCopyGeneNorm.counts_bloo_genes_bloo.tbl',
                          sep = '\t')

genes_bloo_tb_tax <- genes_bloo_tb |>
  left_join(genes_tax_tb_f, by = c('gene')) |>
  dplyr::filter(!is.na(tax_bloo_close))

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

genes_bloo_ko_1102_tb %$%
  contribution_to_ko |>
  range()

## upload information about these KOs 
top25_kos_amy |>
  cat(top25_kos_gla)

top_25_kos_names <- read.csv2('data/genes_sergio/top25_genes_description.csv')

kos_modules <- read.csv('data/genes_sergio/kos_modules.tbl', sep = '\t') |>
  dplyr::mutate(module = str_replace(module, 'md:', ''),
                ko = str_replace(ko, 'ko:', ''))
  
top_25_kos_name_module <- top_25_kos_names |>
  left_join(kos_modules, by = c('KO' ='ko'))

top_25_kos_name_module |>
  distinct(module)

#write.table(top_25_kos_name_module, 'data(genes_sergio/top25_kos_modules_scg.txt', sep = '\t')

top_25_kos_names |>
  distinct(module)

top_25_kos_names |>
  colnames()

top_25_kos_names <- read.csv2('data/genes_sergio/top25_genes_description_scg.csv', sep = '\t')

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
### calculate relative participation of my bloomers tax to each specific KO ----
relative_contribution_gla <- genes_bloo_ko_1101_tb |>
  left_join(kos_tb_t_l_bl1101, by = c('KO' = 'asv_num', 'sample_id' = 'annot')) |>
  dplyr::mutate(relative_contribution = as.numeric(contribution_to_ko)/as.numeric(scg_kos_counts)) |>
  left_join(top_25_kos_names, by = c('KO' = 'KO'))

relative_contribution_gla_plot <- relative_contribution_gla |>
  dplyr::mutate(ko_gene = paste0(KO, ' ', gene)) |>
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

# counts_contribution_gla_plot <- relative_contribution_gla |>
#   dplyr::mutate(ko_gene = paste0(KO, ' ', gene)) |>
#   dplyr::mutate(ko_gene = as.factor(ko_gene)) |>
#   dplyr::mutate(ko_gene = fct_reorder(ko_gene, relative_contribution, .desc = F)) |>  # Reorder factor
#   ggplot(aes(value, fct_infreq(ko_gene)))+
#   geom_col(aes(), fill = "#f1c510")+
#   labs(x = 'SGC counts', y = '')+
#   #scale_x_continuous( expand = c(0,0), limits = c(0, 0.45))+
#   theme_bw()+
#   theme(panel.grid.minor = element_blank(),
#         axis.text.y = element_text(size = 0), axis.ticks.y = element_blank(),
#         axis.text.x = element_text(size = 4),
#         axis.title.x = element_text(size = 5))
# 
# counts_contribution_gla_plot
relative_contribution_amy <- genes_bloo_ko_0904_tb |>
  left_join(kos_tb_t_l_bl0904, by = c('KO' = 'asv_num', 'sample_id' = 'annot')) |>
  dplyr::mutate(relative_contribution = as.numeric(contribution_to_ko)/as.numeric(scg_kos_counts)) |>
  left_join(top_25_kos_names, by = c('KO' = 'KO'))

relative_contribution_amy_plot <- relative_contribution_amy |>
  dplyr::mutate(ko_gene = paste0(KO, ' ', gene)) |>
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

# counts_contribution_amy_plot <- relative_contribution_amy |>
#   dplyr::mutate(ko_gene = paste0(KO, ' ', gene)) |>
#   dplyr::mutate(ko_gene = as.factor(ko_gene)) |>
#   dplyr::mutate(ko_gene = fct_reorder(ko_gene, relative_contribution, .desc = F)) |>  # Reorder factor
#   ggplot(aes(value, fct_rev(ko_gene)))+
#   geom_col(aes(), fill = "#2d373b")+
#   labs(x = 'SGC counts', y = '')+
#   scale_x_continuous( expand = c(0,0), limits = c(0, 0.4))+
#   theme_bw()+
#   theme(axis.text.y = element_text(size = 0), axis.ticks.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size = 4),
#         axis.title.x = element_text(size = 5))
# 
# counts_contribution_amy_plot

plot_grid(relative_contribution_gla_plot, relative_contribution_amy_plot, 
          nrow = 1,
          rel_widths = c(1, 1))

## identify KOs not shared between bloomers tax -----
exclusive_kos <- kos_bloomers_tax |>
  distinct(tax_bloo_close, KO) |>
  dplyr::filter(!is.na(tax_bloo_close)) |>
  group_by( KO) |>
  dplyr::reframe(n = n()) |>
  dplyr::filter(n == 1 )

kos_bloomers_tax |>
  distinct(KO) |>
  as.vector()

exclusive_kos |>
  distinct(KO)

kos_rare$asv_num 

kos_bloomers_tax |>
  dplyr::filter(KO %in% kos_rare$asv_num )  #0

kos_tb_t_l |>
  ungroup() |>
  distinct(asv_num) # 8510

kos_bloomers_tax |>
  distinct(KO) # 1617

genes_kos_tax  <-  genes_tax_tb |>
  left_join(kos_bloomers_tax) |>
   dplyr::filter(!is.na(KO))
 
 genes_kos_tax |>
   distinct(KO)
  
kos_bloomers_tb <- kos_tb_t_l |>
    dplyr::filter(asv_num %in% kos_bloomers_tax$KO) |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_ed')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(genes_kos_tax, by = c('asv_num' = 'KO'))

kos_bloomers_tb |>
  distinct(tax)
 
kos_bloomers_tb |>
  dplyr::filter(asv_num %in% c('K13661', 'K04081')) |> #GumC protein K13661
  ggplot(aes(date, relative_abundance))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_area(aes(fill = asv_num), position = 'stack')+
  #geom_line(aes(group = asv_num))+
  #scale_y_log10()+
  #geom_boxplot(aes(group = date))+
  facet_wrap(vars(asv_num))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'none')

plot_bloomers_kos <- kos_bloomers_tb |>
  dplyr::filter(str_detect(tax, 'Glaciecola')) |>
  ggplot(aes(date, relative_abundance))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_area(aes(fill = asv_num), position = 'stack')+
  #geom_line(aes(group = asv_num))+
  #scale_y_log10()+
  #geom_boxplot(aes(group = date))+
  facet_wrap(vars(tax))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'none')

plot_bloomers_kos
  
plot_bloomers_kos <- kos_bloomers_tb |>
  dplyr::filter(str_detect(tax, 'Amylibacter')) |>
  ggplot(aes(date, relative_abundance))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_area(aes(fill = asv_num), position = 'stack')+
  #geom_line(aes(group = asv_num))+
  #scale_y_log10()+
  #geom_boxplot(aes(group = date))+
  facet_wrap(vars(tax))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'none')

plot_bloomers_kos

plot_bloomers_kos <- kos_bloomers_tb |>
  dplyr::filter(str_detect(tax, '121220-bin8')) |>
  ggplot(aes(date, relative_abundance))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_area(aes(fill = asv_num), position = 'stack')+
  #geom_line(aes(group = asv_num))+
  #scale_y_log10()+
  #geom_boxplot(aes(group = date))+
  facet_wrap(vars(tax))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'none')

plot_bloomers_kos

 kos_bloomers_tb_exclusive_tb <- kos_bloomers_tb |>
  dplyr::filter(asv_num %in% exclusive_kos$KO) |>
  dplyr::mutate(tax_bloo_close = case_when(str_detect(tax, '121220-bin8') ~ 'NS4 marine group',
                                           str_detect(tax, 'Amylibacter') ~ 'Amylibacter',
                                           str_detect(tax, 'Glaciecola') ~ 'Glaciecola')) |>
  dplyr::filter(!is.na(tax_bloo_close))
 
 KOs_exclusive_bloo_tax_plot <- kos_bloomers_tb_exclusive_tb |> 
   ggplot(aes(date, relative_abundance))+
   geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
   geom_area(aes(fill = asv_num), position = 'stack')+
   #geom_line(aes(group = asv_num))+
   #scale_y_log10()+
   #geom_boxplot(aes(group = date))+
   facet_wrap(vars(tax_bloo_close), ncol = 1)+
   scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
   labs(x = 'Date', y = 'KOs counts SCG')+
   theme_bw()+
   theme(axis.text.y = element_text(size = 4), legend.position = 'none')
 
 KOs_exclusive_bloo_tax_plot

## -----
kos_modules <- read.csv('data/genes_sergio/kos_modules.tbl', sep = '\t')

kos_modules |>
  colnames() <- c('module', 'KO')

kos_modules <- kos_modules |>
  dplyr::mutate(module = str_replace(module, 'md:|md.', ''),
                KO = str_replace(KO, 'ko.|ko:', ''))

genes_names <- read.csv('data/genes_sergio/kos_names_bloo_tax.txt', sep = '\t')

modules_description <- fread('data/genes_sergio/modules_kos_bloo.txt', sep = '\t', fill = T) |>
  dplyr::mutate(module_description = paste0(V2, ' ', V3, ' ', V4, ' ', V5, ' ', V6, ' ', V7, ' ', V8, ' ', V9, ' ', V10, ' ', V11, ' ', V12, ' ', V13)) |>
  dplyr::select(module = V1, module_description)

genes_names |>
  head()

genes_names_module_tb <- genes_names |>
  distinct(KO, gene, gene_description) |>
  left_join(kos_modules, by = c('KO')) |>
  left_join(modules_description) |>
  distinct(KO, gene, gene_description, module, module_description)

kos_modules |>
  head()

genes_names_module_tb |>
  colnames()

genes_names_module_tb |>
  distinct(module)

kos_abund <- kos_bloomers_tb_exclusive_tb |>
  dplyr::filter(relative_abundance > 0.001) |>
  distinct(asv_num)

str(kos_abund$asv_num) # Check type of asv_num
str(kos_bloomers_tb_exclusive_tb_ed$KO) # Check type of KO
unique(kos_bloomers_tb_exclusive_tb_ed$KO_ed)

kos_bloomers_tb_exclusive_tb_ed <- kos_bloomers_tb_exclusive_tb |>
  dplyr::select(annot, KO = asv_num, relative_abundance, tax_bloo_close, date ) |>
  left_join(genes_names, by = c('KO')) |>
  dplyr::mutate(gene = case_when(is.na(gene) ~ 'unknown',
                                 TRUE ~ gene)) |>
  dplyr::mutate(KO_ed = case_when(KO %in% kos_abund$asv_num ~ KO, TRUE ~ 'other KO'))

kos_bloomers_tb_exclusive_tb_ed$KO_ed |>
  unique()

kos_bloomers_tb_exclusive_tb_ed |>
  distinct(KO_ed, annot, date, gene, tax_bloo_close, relative_abundance) |>
  dplyr::group_by(KO_ed, annot, date, gene, tax_bloo_close) |>
  dplyr::mutate(sum_ko = sum(relative_abundance)) |>
  dplyr::filter(KO_ed != 'other KO') |>
  ggplot(aes(date, sum_ko))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_area(aes(fill = KO_ed), position = 'stack')+
  facet_wrap(vars(tax_bloo_close), ncol = 1)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'bottom')

kos_bloomers_tb_exclusive_tb |>
  dplyr::select(annot, KO = asv_num, relative_abundance, tax_bloo_close, date ) |>
  distinct(KO, annot, date,  tax_bloo_close, relative_abundance) |>
  dplyr::filter(tax_bloo_close == 'Amylibacter') |>
  dplyr::mutate(sum_ko = sum(relative_abundance)) |>
  ggplot(aes(date, sum_ko))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_line()+
  #geom_area(aes(fill = KO_ed), position = 'stack')+
  facet_wrap(vars(tax_bloo_close), ncol = 1)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'bottom')

kos_bloomers_tb_exclusive_tb |>
  dplyr::select(annot, KO = asv_num, relative_abundance, tax_bloo_close, date ) |>
  distinct(KO, annot, date,  tax_bloo_close, relative_abundance) |>
  dplyr::filter(tax_bloo_close == 'NS4 marine group') |>
  dplyr::mutate(sum_ko = sum(relative_abundance)) |>
  ggplot(aes(date, sum_ko))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_line()+
  #geom_area(aes(fill = KO_ed), position = 'stack')+
  facet_wrap(vars(tax_bloo_close), ncol = 1)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'bottom')

kos_bloomers_tb_exclusive_tb$tax_bloo_close |>
  unique()

plot_kos_percentage_bloo <- kos_bloomers_tb_exclusive_tb |>
  dplyr::select(annot, KO = asv_num, relative_abundance, tax_bloo_close, date ) |>
  distinct(KO, annot, date,  tax_bloo_close, relative_abundance) |>
  #dplyr::filter(tax_bloo_close == 'NS4 marine group') |>
  dplyr::group_by(tax_bloo_close, date) |>
  dplyr::mutate(sum_ko = sum(relative_abundance)) |>
  ggplot(aes(date, sum_ko))+
  geom_vline(xintercept = heterotrophic_blooms$date, linetype = 'dashed')+
  geom_line(aes(group = tax_bloo_close, color = tax_bloo_close), linewidth = 1)+
  #geom_area(aes(fill = KO_ed), position = 'stack')+
  #facet_wrap(vars(tax_bloo_close), ncol = 1)+
  scale_color_manual(values = c('NS4 marine group' = "#0051BF" , 'Amylibacter' = "#2d373b",  'Glaciecola' = "#f1c510")) +
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'KOs counts SCG')+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4), legend.position = 'bottom')
  
## TOP10 contributing KOs of each bloomer ----
  top10_kos_glaciecola <- kos_bloomers_tb_exclusive_tb |>
    distinct(KO = asv_num, annot, date,  tax_bloo_close, relative_abundance, annot) |>
    dplyr::filter(str_detect(annot, 'BL1101')) |>
    dplyr::filter(tax_bloo_close == 'Glaciecola') |>
    slice_max(order_by = relative_abundance, n = 20) %$%
    as.vector(KO)
  
# K17680: PEO1; twinkle protein [EC:5.6.2.3]
# Function: The Twinkle protein (PEO1) is involved in DNA replication. It is a helicase that unwinds DNA, and is particularly important for mitochondrial DNA replication. The enzyme has an EC number of 5.6.2.3, which means it is an ATP-dependent helicase involved in the unwinding of DNA.
# K03406: mcp; methyl-accepting chemotaxis protein
# Function: MCPs (Methyl-accepting chemotaxis proteins) are involved in bacterial chemotaxis, where they sense environmental changes (such as the presence of certain chemicals) and relay this information to the signaling pathways inside the cell, helping the organism move toward or away from these signals.
# K07497: putative transposase
# Function: This is a putative (suggested) transposase, which is an enzyme that catalyzes the movement of genetic material (transposition) within and between DNA molecules, particularly involved in the movement of transposons (genetic elements that can change positions within a genome).
# K18138: acrB, mexB, adeJ, smeE, mtrD, cmeB; multidrug efflux pump
# Function: This is a multidrug efflux pump, responsible for pumping out toxic substances, including antibiotics, from cells. These pumps contribute to antibiotic resistance in bacteria by actively removing drugs from the bacterial cell, thus preventing the drugs from reaching their target sites.
# K00163: aceE; pyruvate dehydrogenase E1 component [EC:1.2.4.1]
# Function: The pyruvate dehydrogenase E1 component is involved in the conversion of pyruvate to acetyl-CoA, a critical step in cellular respiration and energy production. It has an EC number of 1.2.4.1, meaning it catalyzes the oxidative decarboxylation of pyruvate.
# K03582: recB; exodeoxyribonuclease V beta subunit [EC:3.1.11.5]
# Function: The RecB subunit of exodeoxyribonuclease V is involved in the repair of DNA damage, particularly double-strand breaks. It works as part of the RecBCD complex in bacteria and is crucial for recombination and repair processes. The EC number 3.1.11.5 indicates it is an exonuclease that degrades single-stranded DNA.
# K01682: acnB; aconitate hydratase 2 / 2-methylisocitrate dehydratase [EC:4.2.1.3 4.2.1.99]
# Function: Aconitate hydratase 2 is involved in the citric acid cycle (Krebs cycle), catalyzing the conversion of citrate to aconitate. It also has activity as a 2-methylisocitrate dehydratase, which is involved in a specific metabolic pathway in some bacteria. The EC numbers 4.2.1.3 and 4.2.1.99 describe this enzyme's hydratase and dehydratase functions, respectively.
# K03308: TC.NSS; neurotransmitter:Na+ symporter, NSS family
# Function: This is a neurotransmitter sodium symporter (NSS), which is a type of membrane protein that uses the gradient of sodium ions (Na+) to transport neurotransmitters across membranes, playing a key role in neurotransmitter reuptake in the brain.
# K08307: mltD, dniR; peptidoglycan lytic transglycosylase D [EC:4.2.2.29]
# Function: This enzyme is involved in peptidoglycan degradation, an important process in bacterial cell wall remodeling. It has a lytic transglycosylase activity, where it catalyzes the cleavage of glycosidic bonds in peptidoglycan. The EC number 4.2.2.29 indicates it is a glycosylhydrolase involved in the breakdown of peptidoglycan.
# K03722: dinG; ATP-dependent DNA helicase DinG [EC:5.6.2.3]
# Function: DinG is an ATP-dependent DNA helicase that unwinds DNA, particularly in the context of DNA repair. It is involved in the repair of stalled replication forks and other forms of DNA damage. The EC number 5.6.2.3 indicates it is a helicase.

  heterotrophic_blooms
  
  top10_kos_amilybact <- kos_bloomers_tb_exclusive_tb |>
    distinct(KO = asv_num, annot, date,  tax_bloo_close, relative_abundance, annot) |>
    dplyr::filter(str_detect(annot, 'BL0904')) |>
    dplyr::filter(tax_bloo_close == 'Amylibacter') |>
    slice_max(order_by = relative_abundance, n = 20) %$%
    as.vector(KO)
  
  # K00315  DMGDH; dimethylglycine dehydrogenase [EC:1.5.8.4]
  # K01998  livM; branched-chain amino acid transport system permease protein
  # K00303  soxB; sarcosine oxidase, subunit beta [EC:1.5.3.24 1.5.3.1]
  # K02001  proW; glycine betaine/proline transport system permease protein
  # K01997  livH; branched-chain amino acid transport system permease protein
  # K17675  SUPV3L1, SUV3; ATP-dependent RNA helicase SUPV3L1/SUV3 [EC:5.6.2.6]
  # K01006  ppdK; pyruvate, orthophosphate dikinase [EC:2.7.9.1]
  # K02032  ddpF; peptide/nickel transport system ATP-binding protein
  # K02057  ABC.SS.P; simple sugar transport system permease protein
  # K19191  mabO; 4-methylaminobutanoate oxidase (formaldehyde-forming) [EC:1.5.3.19]
  
top10_kos_ns4 <-  kos_bloomers_tb_exclusive_tb |>
    distinct(KO = asv_num, annot, date,  tax_bloo_close, relative_abundance, annot) |>
    dplyr::filter(str_detect(annot, 'BL1205')) |>
    dplyr::filter(tax_bloo_close == 'NS4 marine group') |>
    slice_max(order_by = relative_abundance, n = 20) %$%
    as.vector(KO)
  
  # K15987  hppA; K(+)-stimulated pyrophosphate-energized sodium pump [EC:7.2.3.1]
  # K00184  dmsB; dimethyl sulfoxide reductase iron-sulfur subunit
  # K00284  GLU, gltS; glutamate synthase (ferredoxin) [EC:1.4.7.1]
  # K01880  GARS, glyS1; glycyl-tRNA synthetase [EC:6.1.1.14]
  # K16052  ynaI, mscMJ; MscS family membrane protein
  # K03405  chlI, bchI; magnesium chelatase subunit I [EC:6.6.1.1]
  # K07456  mutS2; DNA mismatch repair protein MutS2
  # K03696  clpC; ATP-dependent Clp protease ATP-binding subunit ClpC
  # K01697  CBS; cystathionine beta-synthase [EC:4.2.1.22]
  # K12257  secDF; SecD/SecF fusion protein
  
##  ------- Top 25 KOs for each tax: -------
## recover general abundance of the community 
heterotrophic_blooms <- heterotrophic_blooms |>
  dplyr::mutate(tax = case_when(asv_num == 'asv11' ~ 'Glaciecola',
                                asv_num == 'asv58'  ~ 'NS4 Marine Group',
                                asv_num == 'asv27'  ~ 'Amylibacter'))

top10_kos_ns4
top10_kos_amilybact
top10_kos_glaciecola

kos_top10_ns4_tb <- kos_bloomers_tb_exclusive_tb |>
  distinct(KO = asv_num, annot, date,  tax_bloo_close, relative_abundance, annot) |>
  dplyr::filter(str_detect(annot, 'BL1205')) |>
  dplyr::filter(tax_bloo_close == 'NS4 marine group') |>
  dplyr::filter(KO %in% top10_kos_ns4) |>
  rename(relative_abundance_ns4 = relative_abundance)

kos_top10_amy_tb <- kos_bloomers_tb_exclusive_tb |>
  distinct(KO = asv_num, annot, date,  tax_bloo_close, relative_abundance, annot) |>
  dplyr::filter(str_detect(annot, 'BL0904')) |>
  dplyr::filter(tax_bloo_close == 'Amylibacter') |>
  dplyr::filter(KO %in% top10_kos_amilybact) |>
  rename(relative_abundance_amy = relative_abundance)

kos_top10_gla_tb <- kos_bloomers_tb_exclusive_tb |>
  distinct(KO = asv_num, annot, date,  tax_bloo_close, relative_abundance, annot) |>
  dplyr::filter(str_detect(annot, 'BL1101')) |>
  dplyr::filter(tax_bloo_close == 'Glaciecola') |>
  dplyr::filter(KO %in% top10_kos_glaciecola) |>
  rename(relative_abundance_gla = relative_abundance)

## general table of the same KOs -----
kos_top10_ns4_com_tb <- kos_tb_t_l |>
  rename(KO = asv_num) |>
  dplyr::filter(KO %in% top10_kos_ns4) |>
  dplyr::filter(str_detect(annot, 'BL1205')) |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_ed')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::ungroup() |>
  dplyr::select( KO, annot, relative_abundance_com = relative_abundance)

kos_top10_ns4_tb_all <- kos_top10_ns4_tb |>
  left_join(kos_top10_ns4_com_tb) |>
  left_join(genes_names_module_tb) |>
  dplyr::mutate(module = case_when(is.na(module) ~ 'undefined',
                                   TRUE ~ module)) |>
  rename(relative_abundance_bloo = relative_abundance_ns4) |>
  dplyr::mutate(relative_abundance_bloo == relative_abundance_com)
  
kos_top10_ns4_tb_all |>
  dplyr::mutate(ko_gen = paste0(KO, ' ', gene)) |>
  ggplot(aes(relative_abundance_bloo, interaction(module, ko_gen)))+
  geom_col()

kos_top10_amy_com_tb <- kos_tb_t_l |>
  rename(KO = asv_num) |>
  dplyr::filter(KO %in% top10_kos_amilybact) |>
  dplyr::filter(str_detect(annot, 'BL0904')) |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_ed')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::ungroup() |>
  dplyr::select( KO, annot, relative_abundance_com = relative_abundance)

kos_top10_amy_tb_all <- kos_top10_amy_tb |>
  left_join(kos_top10_amy_com_tb) |>
  left_join(genes_names_module_tb) |>
  rename(relative_abundance_bloo = relative_abundance_amy) |>
  dplyr::mutate(module = case_when(is.na(module) ~ 'undefined',
                                   TRUE ~ module)) 

kos_top10_amy_tb_all |>
  dplyr::mutate(ko_gen = paste0(KO, ' ', gene)) |>
  ggplot(aes(relative_abundance_bloo, interaction(module, ko_gen)))+
  geom_col()

kos_top10_gla_com_tb <- kos_tb_t_l |>
  rename(KO = asv_num) |>
  dplyr::filter(KO %in% top10_kos_glaciecola) |>
  dplyr::filter(str_detect(annot, 'BL1101')) |>
  left_join(m_bbmo_10y_ed_4metag_red2, by = c('sample_id' = 'sample_id_ed')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::ungroup() |>
  dplyr::select( KO, annot, relative_abundance_com = relative_abundance)

kos_top10_gla_tb_all <- kos_top10_gla_tb |>
  left_join(kos_top10_gla_com_tb) |>
  left_join(genes_names_module_tb) |>
  dplyr::mutate(module = case_when(is.na(module) ~ 'undefined',
                                   TRUE ~ module)) |>
  rename(relative_abundance_bloo = relative_abundance_gla) |>
  dplyr::mutate(relative_abundance_bloo == relative_abundance_com)

kos_top10_gla_tb_all |>
  dplyr::mutate(ko_gen = paste0(KO, ' ', gene)) |>
  ggplot(aes(relative_abundance_bloo, interaction(module, ko_gen)))+
  geom_col()

## plot the all together
genes_names_top25 <- read.csv('data/genes_sergio/kos_names_top25.txt', sep = '\t')

top25_kos_blooming_events <- kos_top10_gla_tb_all |>
  bind_rows(kos_top10_amy_tb_all) |>
  bind_rows(kos_top10_ns4_tb_all) |>
  dplyr::select(-gene, -gene_description) |>
  left_join(genes_names_top25, by = c('KO' = 'KO')) |>
  dplyr::mutate(module_gene = paste0(module, ' ', gene, KO)) 

top25_kos_blooming_events |>
  ggplot(aes(relative_abundance_bloo, module_gene))+
  geom_col()+
  facet_wrap(vars(tax_bloo_close), nrow = 1, scales = 'free_y')+
  theme_bw()

top25_kos_blooming_events$KO |>
  as.vector()

kos_top10_gla_tb_all$KO |>
  unique()

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

## ---- Staphylococcus crisis -----
 asv_tab_10y_02_pseudo_rclr <- read.csv2('data/asv_tab_10y_02_pseudo_rclr.csv', sep = ';')  ## community table

 asv_tab_10y_02_tb <- asv_tab_10y_02_pseudo_rclr |>
   dplyr::select(asv_num, date, relative_abundance, year, fraction) |>
   left_join(tax_bbmo_10y_new, by = c('asv_num')) 

 asv_tab_10y_02_tb %$%
   range(relative_abundance)

others <-  asv_tab_10y_02_tb |>
   group_by(date, fraction) |>
   dplyr::reframe(sum = sum(relative_abundance)) |>
   dplyr::filter(sum != 1.00) |>
   arrange(sum) |>
   dplyr::mutate(missing = 1-sum) |>
   dplyr::mutate(order_ed_f = 'others',
                 relative_abundance = missing) |>
  dplyr::select(-missing, -sum)

 asv_tab_10y_02_tb_f  <-  asv_tab_10y_02_tb |>
   dplyr::mutate(order_ed = case_when(order == 'Staphylococcales'~  order ,
                                       TRUE ~  'others' ))

palette_order_assigned_bloo <-  c("SAR11 clade" =      "#ca6094",
                                  "SAR11 clade asv15" =      "#ca6094",
                                  "Rhodospirillales"  = '#FFA180', 
                                  "Sphingomonadales"  = '#8C000A', 
                                  "Puniceispirillales" = '#4cb50f',
                                  'Rhizobiales' = '#B31722',    
                                  "Rhodobacterales" = '#2d373b',
                                  "Rhodobacterales asv27" = '#2d373b',
                                  "Verrucomicrobiales"= '#005c69',
                                  "Opitutales"   =   '#74B9C8',
                                  "Phycisphaerales"  = '#e3a6ce', 
                                  "Flavobacteriales"   =  '#0051BF', 
                                  "Flavobacteriales asv58"   =  '#0051BF', 
                                  'Chitinophagales' = '#92ABFF', 
                                  "Synechococcales"  = '#009F6A', 
                                  "Synechococcales asv1"  = '#009F6A', 
                                  "Synechococcales asv7"  = '#009F6A', 
                                  "Bacteriovoracales" = '#8C789D',
                                  "Pseudomonadales"  = '#FF8E00', 
                                  "Enterobacterales" = "#f1c510",
                                  "Enterobacterales asv11" = "#f1c510",
                                  'Thiotrichales' =  "#000000",
                                  'env' = "#4A785C", 
                                  'Staphylococcales' = '#B31722',
                                  'others' = 'black') 

bbmo_staphylo_plot <-  asv_tab_10y_02_tb_f  |>
  group_by(date, fraction, order_ed) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::reframe(max_abund = sum(relative_abundance)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   # limits = c(min(asv_tab_10y_02_tb$date), max(asv_tab_10y_02_tb$date),
                   # limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  geom_area(aes( fill = order_ed), alpha = 0.8,  position ='stack', color = 'black', linewidth = 0.1)+
  geom_vline(xintercept = heterotrophic_blooms$date, color = 'black', linetype = 'dashed', linewidth = 0.1)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  #scale_color_identity()+
  scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  guides(alpha = 'none')+
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Order')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_env)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 3), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom',
        axis.title = element_text(size = 4), strip.background = element_blank(),
        legend.text = element_text(size = 3), legend.title = element_text(size = 5), strip.placement = 'outside',
        plot.margin = margin(5,5,5,5), axis.text.y = element_text(size = 4),
        legend.key.size = unit(0.25,'lines'))

bbmo_staphylo_plot

ggsave(filename = 'bbmo_staphylo_plot_02.pdf', plot =  bbmo_staphylo_plot,
       path = 'results/figures/',
       width = 100, height = 80, units = 'mm')

asv_tab_10y_3_pseudo_rclr <- read.csv2('data/asv_tab_10y_3_pseudo_rclr.csv', sep = ';')  ## community table

asv_tab_10y_3_tb <- asv_tab_10y_3_pseudo_rclr |>
  dplyr::select(asv_num, date, relative_abundance, year, fraction) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) 

asv_tab_10y_3_tb %$%
  range(relative_abundance)

others <-  asv_tab_10y_3_tb |>
  group_by(date, fraction) |>
  dplyr::reframe(sum = sum(relative_abundance)) |>
  dplyr::filter(sum != 1.00) |>
  arrange(sum) |>
  dplyr::mutate(missing = 1-sum) |>
  dplyr::mutate(order_ed_f = 'others',
                relative_abundance = missing) |>
  dplyr::select(-missing, -sum)

asv_tab_10y_3_tb_f  <-  asv_tab_10y_3_tb |>
  dplyr::mutate(order_ed = case_when(order == 'Staphylococcales'~  order ,
                                     TRUE ~  'others' ))

bbmo_staphylo_plot <-  asv_tab_10y_3_tb_f  |>
  group_by(date, fraction, order_ed) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::reframe(max_abund = sum(relative_abundance)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   # limits = c(min(asv_tab_10y_3_tb$date), max(asv_tab_10y_3_tb$date),
                   # limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  geom_area(aes( fill = order_ed), alpha = 0.8,  position ='stack', color = 'black', linewidth = 0.1)+
  geom_vline(xintercept = heterotrophic_blooms$date, color = 'black', linetype = 'dashed', linewidth = 0.1)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  #scale_color_identity()+
  scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  guides(alpha = 'none')+
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Order')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_env)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 3), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom',
        axis.title = element_text(size = 4), strip.background = element_blank(),
        legend.text = element_text(size = 3), legend.title = element_text(size = 5), strip.placement = 'outside',
        plot.margin = margin(5,5,5,5), axis.text.y = element_text(size = 4),
        legend.key.size = unit(0.25,'lines'))

bbmo_staphylo_plot

ggsave(filename = 'bbmo_staphylo_plot_3.pdf', plot =  bbmo_staphylo_plot,
       path = 'results/figures/',
       width = 100, height = 80, units = 'mm')


## BACILLI
asv_tab_10y_3_pseudo_rclr <- read.csv2('data/asv_tab_10y_3_pseudo_rclr.csv', sep = ';')  ## community table

asv_tab_10y_3_tb <- asv_tab_10y_3_pseudo_rclr |>
  dplyr::select(asv_num, date, relative_abundance, year, fraction) |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) 

asv_tab_10y_3_tb %$%
  range(relative_abundance)

others <-  asv_tab_10y_3_tb |>
  group_by(date, fraction) |>
  dplyr::reframe(sum = sum(relative_abundance)) |>
  dplyr::filter(sum != 1.00) |>
  arrange(sum) |>
  dplyr::mutate(missing = 1-sum) |>
  dplyr::mutate(order_ed_f = 'others',
                relative_abundance = missing) |>
  dplyr::select(-missing, -sum)

asv_tab_10y_3_tb_f  <-  asv_tab_10y_3_tb |>
  dplyr::mutate(order_ed = case_when(class == 'Bacilli'~  class ,
                                     TRUE ~  'others' ))

bbmo_staphylo_plot <- asv_tab_10y_3_tb_f  |>
  group_by(date, fraction, order_ed) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::reframe(max_abund = sum(relative_abundance)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   # limits = c(min(asv_tab_10y_3_tb$date), max(asv_tab_10y_3_tb$date),
                   # limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  geom_area(aes( fill = order_ed), alpha = 0.8,  position ='stack', color = 'black', linewidth = 0.1)+
  geom_vline(xintercept = heterotrophic_blooms$date, color = 'black', linetype = 'dashed', linewidth = 0.1)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  #scale_color_identity()+
  #scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  guides(alpha = 'none')+
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Order')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction_env)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 3), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom',
        axis.title = element_text(size = 4), strip.background = element_blank(),
        legend.text = element_text(size = 3), legend.title = element_text(size = 5), strip.placement = 'outside',
        plot.margin = margin(5,5,5,5), axis.text.y = element_text(size = 4),
        legend.key.size = unit(0.25,'lines'))

bbmo_staphylo_plot

# ggsave(filename = 'bbmo_staphylo_plot_3.pdf', plot =  bbmo_staphylo_plot,
#        path = 'results/figures/',
#        width = 100, height = 80, units = 'mm')