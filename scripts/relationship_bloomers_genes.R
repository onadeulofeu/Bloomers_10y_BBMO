
## edit metadata sample_id column so that it matches the sample_id code that we have in the genes data
m_bbmo_10y |>
  colnames()

m_bbmo_10y_ed <- m_bbmo_10y |>
    separate(sample_id, into = c('sample_id_sim', 'fraction', 'code'), sep = '_') |>
    dplyr::filter(fraction == '0.2')

#upload data from metagenomes we have the table with counts for the genes
genes <- read.csv2('data/genes_sergio/K02692_genes_detail.csv') |>
  dplyr::select(annot, starts_with('BL')) |> #filter only Blanes samples not SOLA 
  as_tibble() |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'gene_counts') |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  left_join(m_bbmo_10y_ed, by = c('sample_id' = 'sample_id_sim')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013'))

genes_sample_id <- genes %$%
  sample_id |>
  unique() #I have 56 samples in common with metaG data

m_bbmo_10y_ed_red <- m_bbmo_10y_ed |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(sample_id_sim %in% genes_sample_id)

m_02_f <- m_bbmo_10y_ed |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(sample_id_sim %in% genes_sample_id) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_bbmo_10y_ed_red)))

m_02_f$sample_id_num

genes |>
  colnames()

m_02 |>
  colnames()

genes %$%
  year |>
  range()

## explore a bit this type of data----
genes %$%
  gene_counts |>
  range() #0.00000 17.37556

genes |>
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 1)) |> #filter by genes that have a minimum count at some point
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
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 3/7)
  
## we would like to detect anomalies in genes and see if they correlate with anomalies in ASVs that are selected as potential bloomers
z_02_genes <- genes |>
  group_by(annot) |>
  arrange(date) |>
  group_by(annot) |>
  dplyr::reframe(anomalies_gene = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = gene_counts, plotting = FALSE)[c(1,2,3)])

z_scores_02_genes <- genes |>
  group_by(annot) |>
  dplyr::mutate(gene_counts = as.numeric(gene_counts)) |>
  dplyr::filter(sum(gene_counts) > 0) |>
  left_join(m_bbmo_10y_ed, by = c('sample_id' = 'sample_id_sim')) |>
  #arrange(date) |>
  group_by(annot) |>
  dplyr::reframe(z_score_gen = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = gene_counts, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_gen) |>
  group_by(annot) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02_f))) |>
  left_join(m_02_f, by = 'sample_id_num') 

## filter by the years that we share after computing the z_scores so that I am able to calculate the last ones (crec que no cal perque els NAs els trobaré al final)

## plot the genes with their anomalies an observe if they have correspondence with the blooming events----
z_scores_02_genes |>
  colnames()

genes |>
  dim()

z_scores_02_genes |>
  dim()

genes |>
  left_join(z_scores_02_genes) %$%
  gene_counts |>
  range()

genes |>
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
# |>
#   length()

genes |>
  colnames()

z_scores_02_genes_f |>
  colnames()

genes$sample_id
z_scores_02_genes_f$sample_id_sim

genes_counts_z_scores <- genes |>
  filter(annot %in% genes_anom) |>
  dplyr::select(sample_id, gene_counts, annot, date) |>
  left_join(z_scores_02_genes_f, by = c('sample_id' = 'sample_id_sim', 'annot')) ##

genes$annot
z_scores_02_genes_f$annot
z_scores_02_genes_f$z_score_gen

genes$sample_id
z_scores_02_genes_f$sample_id_sim

genes |>
  dim()
  
z_scores_02_genes_f |>
  dim()

genes_counts_z_scores |>  
  group_by(annot) |>
  dplyr::filter(any(gene_counts > 5)) |> #filter by genes that have a minimum count at some point
  ungroup() |>
  #dplyr::mutate(date.x = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_gen) >= 1.96, '#9F0011', '#080808', missing = '#080808')) |>
  ggplot(aes(date.x, gene_counts))+
  scale_x_datetime(expand = c(0,0),
                   breaks = ( by = '1 year'),
                   date_labels = "%Y")+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  geom_point(aes(color = anomaly_color))+
  scale_color_identity()+
  geom_line(aes(group = annot), color = '#5B5A5A')+ #, color = '#3D3B3B'
  facet_wrap(vars(annot))+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 3/7)
  
## potential bloomers dataset for the same dates as we have metagenomes 2009-2013----
##blooming events 0.2 between 2009-2013 comparison with Metagenomes
asv_tab_all_bloo_z_tax |>
  colnames()

abund_asv_bloom_events_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  as.numeric(z_score_ra) >= 1.96 &
                  as.numeric(abundance_value) >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  left_join(tax_bbmo_10y_new)

write.csv(abund_asv_bloom_events_02, 'data/genes_sergio/abund_asv_bloom_events_02.csv')

## plot them to see if there's some kind of correlation with KEGGs
abund_asv_bloom_events_02_unique <- asv_tab_all_bloo_z_tax |>
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
  dplyr::filter(abundance_type == 'relative_abundance' &
                  z_score_ra >= 1.96 &
                  abundance_value >= 0.1) |> #we add this line in case we want potential bloomers not just anomalies in their rel abund
  dplyr::select(year, sample_id, asv_num, abundance_type, abundance_value, date) |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  left_join(tax_bbmo_10y_new)

#write.csv(abund_asv_bloo_fl_09_13, 'data/genes_sergio/abund_asv_bloo_fl_09_13.csv')

bloo_02_09_13_plot <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013')) |>
  dplyr::filter(asv_num %in% abund_asv_bloom_events_02_unique$asv_num) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  ggplot(aes(date, abundance_value))+
  facet_wrap(vars(asv_num))+
  scale_x_datetime(expand = c(0,0),
                   breaks = seq(min(asv_tab_bloo_rel_abund_z$date), 
                                max(asv_tab_bloo_rel_abund_z$date), by = '1 year'),
                   date_labels = "%Y")+
  #geom_hline(yintercept = 0.1, color = '#8B8989')+
  geom_line(aes(group = asv_num), color = '#5B5A5A')+ #, color = '#3D3B3B'
  geom_point(data = abund_asv_bloo_fl_09_13 |>
               dplyr::filter(year %in% c('2009', '2010', '2011', '2012', '2013') &
                               asv_num %in% abund_asv_bloom_events_02_unique$asv_num),
             aes(date, abundance_value, color =  '#9F0011'), size = 1)+
  scale_y_continuous(labels = percent_format())+
  geom_area()+
  labs(colour = '', alpha = '', y = 'Relative abundance (%)', x = 'Date')+
  theme_bw()+
  guides(guide_legend(colour = 'none',
                      alpha = 'none'))+
  theme(panel.grid = element_blank(), text = element_text(size = 7), panel.border = element_blank(),
       axis.ticks.x = element_blank(), aspect.ratio = 3/7,
        strip.background = element_blank(), legend.position = 'none')

ggsave( plot = bloo_02_09_13_plot, filename = 'bloo_02_09_13_plot.pdf',
        path = 'results/figures/relationship_genes_blooms/',
        width = 180, height = 160, units = 'mm')

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
dplyr::select(-date.y)
  
zscores_gen_ra_fl_09_13 <- rclr_bloo_fl_09_13 |>
  left_join(genes_counts_z_scores_f, by = c('date' = 'date.x'), relationship = 'many-to-many') 

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

