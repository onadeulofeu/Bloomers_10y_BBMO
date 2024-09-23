# packages ----


# upload data ---
asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') |> ##using dada2 classifier assign tax with silva 138.1 and correctly identifying bloomers
  as_tibble() |>
  dplyr::select(-X)

bloo_all_types_summary_tax <- read.csv('results/tables/bloo_all_types_summary_tb_tax_v2.csv')

# palettes ----
palette_occurrence <- c(narrow = "#AE659B",
                        intermediate = "#3e3e3e",
                        broad = "#57a9a8")

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

### The remodelation of the Blanes harbour strated on 24th March 2010 and finished on the 9th of june 2012
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

# ---------------------- METHODS ----------------------  ########## ------
# ------ ########## Figure define bloomers threshold ------  ########## ------
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

# ------ ########## Figure bloomers community timeseries ------ ########## ----------
## I remove the 4 ASVs that clustered together in the seasonality analysis but not others that could be potential blooms
##color by order
## I remove the 4 ASVs that clustered together in the seasonality analysis but not others that could be potential blooms
community_eveness_all_m <- community_eveness_all_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

bbmo_bloo_ev_order_no_sar_cluster <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num %in% bloo_02$value & fraction == '0.2' |
                  asv_num %in% bloo_3$value & fraction == '3' ) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
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
  geom_area(aes(date, abund_order, fill = order_f, group = order_f), alpha = 1,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),  
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
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
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

bbmo_bloo_ev_order_no_sar_cluster

order_legend <- get_legend(bbmo_bloo_ev_order_no_sar_cluster)

bbmo_bloo_ev_order_no_sar_cluster <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num %in% bloo_02$value & fraction == '0.2' |
                  asv_num %in% bloo_3$value & fraction == '3' ) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
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
  geom_area(aes(date, abund_order, fill = order_f, group = order_f), alpha = 1,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),  
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
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
        legend.position = 'none', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

bbmo_bloo_ev_order_no_sar_cluster

# figure B
bloo_all_types_summary_tb_tax <-  bloo_all_types_summary_tb_tax_v2
bloo_all_types_summary_tb_tax$recurrency <- factor(bloo_all_types_summary_tb_tax$recurrency , 
                                                   levels = c('recurrent', 'non-recurrent'))
data_text <- bloo_all_types_summary_tb_tax |>
  group_by(recurrency, fraction, asv_num) |>
  dplyr::reframe(n = n()) |>
  group_by(fraction, recurrency) |>
  dplyr::reframe(n_general = paste0( 'n = ', n())) 

tax_seasonality_plot <- bloo_all_types_summary_tb_tax |>
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
        text = element_text(size = 8))

tax_seasonality_plot 

bbmo_bloo_ev_order_no_sar_cluster_seas_plot  <- grid.arrange(bbmo_bloo_ev_order_no_sar_cluster, 
                                                             tax_seasonality_plot,
                                                             order_legend,
             ncol = 1, heights = c(1, 0.35, 0.5))

# ggsave('bbmo_bloo_ev_order_no_sar_cluster_seas_plot.pdf', bbmo_bloo_ev_order_no_sar_cluster_seas_plot,
#        path = "results/figures/",
#        width = 180,
#        height = 160,
#        units = 'mm')

## -------- ########## Supplementray figure: BBMO bloomers time series separated by the different types of bloomers ------ ########## ---------
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
  #left_join(summary_types_of_blooms, by = c('asv_num', 'fraction')) |>
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
  geom_area(aes(date, abund_class, fill = order_f, group = order_f), alpha = 1,  position='stack')+
  #geom_line(data = bray_curtis_rar_all_m, aes(date, bray_curtis_result))+
  # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
  # geom_point(data = community_eveness_all_m |>
  #              dplyr::filter(anomaly_color == '#9F0011'),  
  #            aes(date, community_eveness_rar/1.6, color = anomaly_color, alpha = 0.8))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
  scale_color_identity()+
  scale_fill_manual(values = palette_order_assigned_bloo, na.value = "#000000")+
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Order')+
  facet_grid(fraction~recurrency,  scales = 'free_y',  labeller =  labs_fraction_rec_freq)+
  #facet_wrap(fraction~phylum_f, dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
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

# ggsave('booming_events_BBMO10Y_new_tax_seasonal_ed3.pdf', booming_events_seas_BBMO10Y,
#        path = "results/figures/",
#        width = 180,
#        height = 160,
#        units = 'mm')

# ------ ########## Figure example of different types of bloomers identified in our dataset ------ ########## ------------------------
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
## how was the diverstiy affected by the harbor restoration ----
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
  #scale_color_manual(values = palette_seasons_4, labels = labs_season)+
  facet_wrap(vars(fct_rev(harbor_restoration)), labeller = labs_perturbation)+
  scale_shape_discrete(labels = labs_fraction)+
  #coord_flip()+
  geom_violin(aes(group = fraction), alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  theme_bw()+
  theme(text = element_text(size = 20),
        strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

shannon_harbor_rar_plot 



# ---------------------- DISCUSSION ----------------------  ########## ------
