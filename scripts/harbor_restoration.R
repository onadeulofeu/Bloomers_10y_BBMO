# pacakges ----
library(tidyverse)
library(magrittr)
library(ggplot2)
library(scales)
library(ggdendro)

# Here we will analyze how did harbor restoration period affected the microbial community at the BBMO
## Based on Roca et al., which analyized the harbor effects on the Posidonia oceanica meadows we know that 
## the most important works were performed during march-june 2010 and then on june 2012 they complete finished.

#upload data

# differences in abundances of bloomers previous - during - after the restoration---- 
### The remodelation of the Blanes harbour strated on 24th March 2010 and finished on the 9th of june 2012
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

asv_tab_all_bloo_z_tax_harbor <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(harbor_restoration = case_when(  date < '2010-03-24' ~ 'pre_perturbation',
                                                 date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
                                                 date >= '2012-06-09' ~ 'post_perturbation'))

asv_tab_all_bloo_z_tax_harbor <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) ##remove those taxa that are not real bloomers

asv_tab_all_bloo_z_tax_harbor$harbor_restoration <- factor(asv_tab_all_bloo_z_tax_harbor$harbor_restoration, 
                                                           levels = c('pre_perturbation',
                                                           'perturbation',
                                                           'post_perturbation'))

asv_tab_all_bloo_z_tax_harbor |>  
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '3') |> ## the community most afected by the harbor restoration period was the PA
  group_by(harbor_restoration, asv_num) |>
  dplyr::reframe(n = n()) |>
  distinct(harbor_restoration, n)

asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '3') |>
  ggplot(aes(harbor_restoration, abundance_value, group = harbor_restoration))+
  geom_point(aes(shape = fraction))+
  geom_boxplot(alpha = 0.5)+
  #facet_wrap(vars(fraction))+
  facet_wrap(vars(asv_num_f))+
  geom_smooth(aes(harbor_restoration, abundance_value), method = 'loess', se = F)+
  theme_bw()+
  theme(text = element_text(size = 6),
        strip.background = element_blank())

asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '0.2') |>
  ggplot(aes(harbor_restoration, abundance_value, group = harbor_restoration))+
  geom_point(aes(shape = fraction))+
  geom_boxplot(alpha = 0.5)+
  geom_smooth(method = 'loess', se = F)+
  #facet_wrap(vars(fraction))+
  facet_wrap(vars(asv_num_f), scales = 'free')+
  geom_smooth(aes(harbor_restoration, abundance_value), method = 'loess', se = F)+
  theme_bw()+
  theme(text = element_text(size = 6),
        strip.background = element_blank())

## We identified that Cyanobacteriia follow a different trend, not exclusive of the harbor restoration period. 
cianobacteria_harbor_restoration <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(class == 'Cyanobacteriia') |>
 # dplyr::filter(fraction == '0.2') |>
  ggplot(aes(harbor_restoration, abundance_value, group = harbor_restoration))+
  geom_point(aes(shape = fraction, color = fraction), position = position_jitter(width = 0.2), size = 1)+
  geom_boxplot(alpha = 0.1)+
  geom_smooth(method = 'loess', se = T, linewidth = 1)+
  scale_x_discrete(labels = c('Pre-disturbance', 'Disturbance', 'Post-disturbance'))+
  #facet_wrap(vars(fraction))+
  facet_wrap(vars(asv_num_f), scales = 'free_y')+
  scale_y_continuous(labels = percent_format())+
  scale_shape_discrete(labels = labs_fraction)+
  labs(y = 'Relative abundance (%)', shape = 'Fraction', color = 'Fraction', x = '')+
  guides(shape = guide_legend(ncol = 2), color = guide_legend(ncol = 2))+
  #facet_wrap(asv_num_f~fraction, scales = 'free')+
  geom_smooth(aes(group = fraction, color = fraction), method = 'loess', se = F)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  theme_bw()+
  theme(text = element_text(size = 7),
        strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        legend.position = 'bottom')

# ggsave(cianobacteria_harbor_restoration, file = 'results/figures/cianobacteria_harbor_restoration.pdf',
#        width = 180, height = 100, units = 'mm')

# Cyanobacteriia----
## We expected that Cyanobacteriia are highly affected by the harbor restoration period since during the harbor restoration period a lot of sediments were resuspended,
## which could lead to an increase of the turbulence and less light available for them ----

asv_tab_all_bloo_z_tax_harbor <- asv_tab_all_bloo_z_tax_harbor %>%
  mutate(date = as.POSIXct(date)) # Convert date column to POSIXct format if needed

asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(class == 'Cyanobacteriia') |>
  ggplot(aes(date, abundance_value, group = asv_num))+
  facet_grid(fraction~asv_num, scales = 'free_y')+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_line(aes())+
  geom_smooth(method = 'loess', color = 'black')+
 #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)')+
  annotate(geom = "text", x = (as.POSIXct('2006-04-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Predisturbance", size = 2)+
  annotate(geom = "text", x = (as.POSIXct('2010-10-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Disturbance", size = 2)+
  annotate(geom = "text", x = (as.POSIXct('2013-04-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Postdisturbance", size = 2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(), text = element_text(size = 7))

## Exploration of data from the harbor restoration period in the study lead by Oscar and Teresa Alcoverro----
### They observed the following: 
## 8 month after the most important sand movements the sediments have 
## return to normality, however, the physiological indicators of the marine 
## phanerogams have not return to normality until 15 month after the most important
## sand movements. 
## But plants density did not return to normality even 38 months post-disturbance (Fig. 6E)
## They observed 

## I plot a heatmap of my potential bloomers to observe what happened with bloomer and the harbor restortion period----
heatmap_data_l <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'rclr') |>
  #dplyr::filter(class == 'Cyanobacteriia') |>
  dplyr::group_by(asv_num, fraction) |>
  dplyr::mutate(sample_num = row_number()) |>
  dplyr::select(asv_num, decimal_date, sample_num, abundance_value,fraction)

 heatmap_data_l |>
  ungroup() |>
  distinct(decimal_date) |>
   as_vector()

heatmap_data_l |>
  ggplot(aes(sample_num, asv_num, fill = abundance_value))+
  geom_tile()+
  facet_wrap(vars(fraction))+
  geom_vline(xintercept = 74, linetype = 'dashed', color = 'grey')+
  geom_vline(xintercept = 79, linetype = 'dashed', color = 'grey')+
  geom_vline(xintercept = 102, linetype = 'dashed', color = 'grey')+
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, na.value = '#D7D6D3')

## For evaluation of differences between samples during the disturbance, before and after I need to explore the patterns they present over the years----
## and remove from the analysis those that have a trend despite the harbor restoration ---
asv_tab_all_bloo_z_tax_harbor |>
  colnames()
asv_tab_all_bloo_z_tax_harbor %$%
  asv_num |>
  unique()

bloo_types_summary %$%
  asv_num |>
  unique()
# 
# asv_tab_all_bloo_z_tax_harbor |>
#   left_join(bloo_types_summary) |>
#   dplyr::filter(frequency != 'seasonal') |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   ggplot(aes(harbor_restoration, abundance_value, group = harbor_restoration))+
#   geom_point(aes(shape = fraction, color = fraction), position = position_jitter(width = 0.2), size = 1)+
#   geom_boxplot(alpha = 0.1)+
#   geom_smooth(method = 'loess', se = T, linewidth = 1)+
#   scale_x_discrete(labels = c('Pre-disturbance', 'Disturbance', 'Post-disturbance'))+
#   #facet_wrap(vars(fraction))+
#   facet_wrap(vars(asv_num_f), scales = 'free_y')+
#   scale_y_continuous(labels = percent_format())+
#   scale_shape_discrete(labels = labs_fraction)+
#   labs(y = 'Relative abundance (%)', shape = 'Fraction', color = 'Fraction', x = '')+
#   guides(shape = guide_legend(ncol = 2), color = guide_legend(ncol = 2))+
#   #facet_wrap(asv_num_f~fraction, scales = 'free')+
#   geom_smooth(aes(group = fraction, color = fraction), method = 'loess', se = F)+
#   scale_color_manual(values = palette_fraction, labels = labs_fraction)+
#   theme_bw()+
#   theme(text = element_text(size = 7),
#         strip.background = element_rect(fill = 'transparent'),
#         panel.grid = element_blank(),
#         legend.position = 'bottom')

asv_tab_all_bloo_z_tax_harbor |>
  left_join(bloo_types_summary) |>
  #dplyr::filter(frequency != 'seasonal') |>
  dplyr::filter(abundance_type == 'rclr') |>
  ggplot(aes(date, abundance_value, group = asv_num))+
  facet_wrap(vars(asv_num), scales = 'free_y',  ncol = 5)+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_line(aes(group = fraction, color = fraction))+
  geom_smooth(aes(group = fraction, color = fraction), method = 'loess')+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)')+
  annotate(geom = "text", x = (as.POSIXct('2006-04-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Predisturbance", size = 2)+
  annotate(geom = "text", x = (as.POSIXct('2010-10-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Disturbance", size = 2)+
  annotate(geom = "text", x = (as.POSIXct('2013-04-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Postdisturbance", size = 2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom', 
        panel.grid = element_blank(), text = element_text(size = 7))

## asv that increase their abundance after the harbor restoration or coupling with the harbor restoration------
### we detect that there are underlying trends under these pattern it is not just the harbor restoration
increase_after_harbor_restoration <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::mutate(famil_asv = paste0(family,' ', asv_num)) |>
  left_join(bloo_types_summary) |>
  #dplyr::filter(frequency != 'seasonal') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(asv_num %in% c('asv1', 'asv15', 'asv17', 'asv23', 'asv31', 'asv4', 'asv7', 'asv80')) |>
  ggplot(aes(date, abundance_value, group = asv_num))+
  facet_wrap(vars(famil_asv), scales = 'free_y',  ncol = 3)+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_point(aes(group = fraction, color = fraction), size = 0.5)+
  geom_smooth(aes(group = fraction, color = fraction), method = 'loess')+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'rCLR', color = 'Fraction')+
  annotate(geom = "text", x = (as.POSIXct('2006-04-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Pre", size = 1)+
  annotate(geom = "text", x = (as.POSIXct('2010-12-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Disturbance", size = 1)+
  annotate(geom = "text", x = (as.POSIXct('2013-06-12', format = "%Y-%m-%d")), 
           y = 6, 
           label = "Post", size = 1)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom', 
        panel.grid = element_blank(),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 4),
        axis.title = element_text(size = 6)) #, text = element_text(size = 5)
# 
# ggsave(increase_after_harbor_restoration, filename = 'increase_after_harbor_restoration.pdf',
#        path = 'Results/Figures/',
#        width = 180, height = 130, units = 'mm')

# asv_tab_all_bloo_z_tax_harbor |>
#   dplyr::mutate(famil_asv = paste0(family,' ', asv_num)) |>
#   left_join(bloo_types_summary) |>
#   #dplyr::filter(frequency != 'seasonal') |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   dplyr::filter(asv_num %in% c('asv1', 'asv15', 'asv17', 'asv23', 'asv31', 'asv4', 'asv7', 'asv80')) |>
#   ggplot(aes(date, abundance_value, group = asv_num))+
#   facet_wrap(vars(famil_asv), scales = 'free',  ncol = 5)+
#   # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
#   #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
#   geom_point(aes(group = fraction, color = fraction), size = 1)+
#   geom_smooth(aes(group = fraction, color = fraction), method = 'loess')+
#   scale_color_manual(values = palette_fraction, labels = labs_fraction)+
#   scale_y_continuous(labels = percent_format())+
#   labs(x = 'Time', y = 'Relative abundance (%)')+
#   annotate(geom = "text", x = (as.POSIXct('2006-04-12', format = "%Y-%m-%d")), 
#            y = 0.2, 
#            label = "Predisturbance", size = 2)+
#   annotate(geom = "text", x = (as.POSIXct('2010-10-12', format = "%Y-%m-%d")), 
#            y = 0.2, 
#            label = "Disturbance", size = 2)+
#   annotate(geom = "text", x = (as.POSIXct('2013-04-12', format = "%Y-%m-%d")), 
#            y = 0.2, 
#            label = "Postdisturbance", size = 2)+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom', 
#         panel.grid = element_blank(), text = element_text(size = 7))

## Taxa that increased or decreased their abundances over the years of the timeseries (INTERANNUAL TRENDS) -----
### Detected with the wavelets transformations residuals. 
s4_fraction_smooth <- wavelets_result_02 |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  bind_rows(wavelets_result_3) |>
  dplyr::filter(wavelets_transformation == 's4') |>
  dplyr::mutate(s4_cluster = case_when(asv_num %in% c('asv62', 'asv38', 'asv11', 'asv27') & fraction == '0.2' ~ 'fl_cl1' ,
                                       asv_num %in% c('asv7', 'asv15')  & fraction == '0.2'~ 'fl_cl2',
                                       asv_num =='asv1'  & fraction == '0.2' ~ 'fl_cl3',
                                       asv_num %in% c('asv249', 'asv114', 'asv237', 'asv282', 'asv563', 'asv555', 'asv178', 'asv58', 'asv17') & fraction == '0.2'~ 'fl_cl4',
                                       asv_num %in% c('asv4', 'asv31', 'asv23') & fraction == '3' ~'pa_cl1',
                                       asv_num %in% c('asv7', 'asv15') & fraction == '3' ~ 'pa_cl2',
                                       asv_num =='asv1'  & fraction == '3' ~ 'pa_cl3',
                                       asv_num %in% c('asv42', 'asv126', 'asv118', 'asv72', 'asv100', 'asv85') & fraction == '3' ~ 'pa_cl4',
                                       asv_num %in% c('asv84', 'asv116', 'asv25', 'asv28', 'asv182', 'asv27', 'asv11') & fraction =='3' ~ 'pa_cl5' ,
                                       asv_num %in% c('asv69', 'asv153', 'asv559', 'asv194', 'asv276', 'asv223', 'asv264')& fraction == '3' ~ 'pa_cl6',
                                       asv_num %in% c('asv317', 'asv200', 'asv311', 'asv511', 'asv113', 'asv43', 'asv225', 
                                                      'asv752', 'asv471', 'asv385') & fraction == '3' ~ 'pa_cl7' ,
                                       asv_num %in% c('asv192', 'asv163', 'asv49', 'asv302', 'asv179', 'asv77',
                                                      'asv219', 'asv105') & fraction == '3' ~ 'pa_cl8',
                                       asv_num %in% c('asv80', 'asv22', 'asv17') & fraction == '3' ~ 'pa_cl9')) |>
  dplyr::mutate(family_asv = paste0(family, ' ', asv_num))

## reodrer asv_num 
s4_fraction_smooth <- s4_fraction_smooth |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num),
                famiy_asv_f = as_factor(family_asv))

s4_fraction_smooth$famiy_asv_f

s4_fraction_smooth$class_f <-  factor(s4_fraction_smooth$class_f, 
                                                                   levels=unique(s4_fraction_smooth$class_f[order(s4_fraction_smooth$phylum_f)]), 
                                                                   ordered=TRUE)

s4_fraction_smooth$order_f <-  factor(s4_fraction_smooth$order_f, 
                                                                   levels=unique(s4_fraction_smooth$order_f[order(s4_fraction_smooth$phylum_f,
                                                                                                                                               s4_fraction_smooth$class_f)]), 
                                                                   ordered=TRUE)

s4_fraction_smooth$family_f <-  factor(s4_fraction_smooth$family_f, 
                                                                    levels=unique(s4_fraction_smooth$family_f[order(s4_fraction_smooth$phylum_f,
                                                                                                                                                 s4_fraction_smooth$class_f,
                                                                                                                                                 s4_fraction_smooth$order_f)]), 
                                                                    ordered=TRUE)


s4_fraction_smooth$asv_num_f <-  factor(s4_fraction_smooth$asv_num_f, 
                                                                     levels=unique(s4_fraction_smooth$asv_num_f[order(s4_fraction_smooth$phylum_f,
                                                                                                                                                   s4_fraction_smooth$class_f,
                                                                                                                                                   s4_fraction_smooth$order_f,
                                                                                                                                                   s4_fraction_smooth$family_f)]), 
                                                                     ordered=TRUE)

s4_plot_fraction_smooth <- s4_fraction_smooth |>
  #dplyr::filter(s4_cluster %in% c('fl_cl2', 'fl_cl3', 'pa_cl1', 'pa_cl2', 'pa_cl3', 'pa_cl5', 'pa_cl6', 'pa_cl8', 'pa_cl9')) |>
  ggplot(aes(decimal_date, wavelets_result)) +
  # geom_vline(xintercept = 2010.24, linetype = 'dashed') +
  # geom_vline(xintercept = 2010.52, linetype = 'dashed') +
  # geom_vline(xintercept = 2012.48, linetype = 'dashed') +
  geom_line(aes(group = fraction, color = as.factor(fraction))) +
  geom_smooth(aes(group = fraction, color = as.factor(fraction)))+
  #facet_grid(vars(fraction), labeller = labs_fraction) +
  facet_wrap(vars(interaction(family_f, asv_num_f)), ncol = 5)+
  # geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                       ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6) +
  labs(x = 'Time', y = 'Wavelet coefficient: s4 (residuals)', color = 'Fraction') +
  scale_color_manual(values = palette_fraction, labels = labs_fraction) +
  #scale_linetype_manual(values = c("solid", "dashed")) + # Define linetypes
  theme_bw() +
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 6),
        strip.text = element_text(margin = margin(2, 2, 2, 2)))

s4_plot_fraction_smooth

# ggsave(s4_plot_fraction_smooth, filename = 's4_plot_fraction_smooth.pdf',
#        path = 'Results/Figures/',
#        width = 180, height = 240, units = 'mm')

### plot the abundances of those taxa that presented an increase over the years -----
blooming_asvs_increasing_over_the_years <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::mutate(famil_asv = paste0(family,' ', asv_num)) |>
  left_join(bloo_types_summary) |>
  #dplyr::filter(frequency != 'seasonal') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(asv_num %in% c('asv1', 'asv15', 'asv7', 'asv17', 'asv23', 'asv31', 'asv4')) |>
  ggplot(aes(date, abundance_value, group = asv_num))+
  facet_wrap(vars(interaction(family_f, asv_num_f)), scales = 'free',  ncol = 3)+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_line(aes(group = fraction, color = fraction), linewidth = 0.3, alpha = 0.6)+
  geom_smooth(aes(group = fraction, color = fraction), method = 'loess', linewidth = 0.8)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'rCLR', color = 'Fraction')+
  # annotate(geom = "text", x = (as.POSIXct('2006-04-12', format = "%Y-%m-%d")), 
  #          y = 0.2, 
  #          label = "Predisturbance", size = 2)+
  # annotate(geom = "text", x = (as.POSIXct('2010-10-12', format = "%Y-%m-%d")), 
  #          y = 0.2, 
  #          label = "Disturbance", size = 2)+
  # annotate(geom = "text", x = (as.POSIXct('2013-04-12', format = "%Y-%m-%d")), 
  #          y = 0.2, 
  #          label = "Postdisturbance", size = 2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom', 
        panel.grid = element_blank(), text = element_text(size = 7))

blooming_asvs_increasing_over_the_years

# ggsave(blooming_asvs_increasing_over_the_years, filename = 'blooming_asvs_increasing_over_the_years.pdf',
#        path = 'Results/Figures/main/',
#        width = 180, height = 180, units = 'mm')

blooming_asvs_decreasing_over_the_years <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::mutate(famil_asv = paste0(family,' ', asv_num)) |>
  left_join(bloo_types_summary) |>
  #dplyr::filter(frequency != 'seasonal') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(asv_num %in% c('asv28')) |>
  ggplot(aes(date, abundance_value, group = asv_num))+
  facet_wrap(vars(interaction(family_f, asv_num_f)), scales = 'free',  ncol = 3)+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_line(aes(group = fraction, color = fraction), linewidth = 0.3, alpha = 0.6)+
  geom_smooth(aes(group = fraction, color = fraction), method = 'loess', linewidth = 0.8)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  #scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'rCLR', color = 'Fraction')+
  # annotate(geom = "text", x = (as.POSIXct('2006-04-12', format = "%Y-%m-%d")), 
  #          y = 0.2, 
  #          label = "Predisturbance", size = 2)+
  # annotate(geom = "text", x = (as.POSIXct('2010-10-12', format = "%Y-%m-%d")), 
  #          y = 0.2, 
  #          label = "Disturbance", size = 2)+
  # annotate(geom = "text", x = (as.POSIXct('2013-04-12', format = "%Y-%m-%d")), 
  #          y = 0.2, 
  #          label = "Postdisturbance", size = 2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom', 
        panel.grid = element_blank(), text = element_text(size = 7))

blooming_asvs_decreasing_over_the_years

# ggsave(blooming_asvs_decreasing_over_the_years, filename = 'blooming_asvs_decreasing_over_the_years.pdf',
#        path = 'Results/Figures/main/',
#        width = 88, height = 60, units = 'mm')

# Compositional plot with the blooming community affected by the harbor restoration -----
## residuals from the wavelet transformation ----
## The residual component, also known as the error or irregular component, captures the random fluctuations or noise in the time series 
## data that cannot be explained by the trend or seasonality. Residuals are the unexplained variability within the data, 
## often reflecting the influence of unpredictable or external factors

harbour_restoration_dec <- tibble(xmin = '2010.24', xmax = '2010.52') |>
  dplyr::mutate(date_min = as.numeric(xmin),
                date_max = as.numeric(xmax))

## ASVs identified to be highly affected by the harbor restoration ,  'asv80', 'asv77' 
wavelets_result_02 |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  bind_rows(wavelets_result_3) |>
  dplyr::filter(wavelets_transformation == 's4') |>
  dplyr::mutate(s4_cluster = case_when(asv_num %in% c('asv62', 'asv38', 'asv11', 'asv27') & fraction == '0.2' ~ 'fl_cl1' ,
                                       asv_num %in% c('asv7', 'asv15')  & fraction == '0.2'~ 'fl_cl2',
                                       asv_num =='asv1'  & fraction == '0.2' ~ 'fl_cl3',
                                       asv_num %in% c('asv249', 'asv114', 'asv237', 'asv282', 'asv563', 'asv555', 'asv178', 'asv58', 'asv17') & fraction == '0.2'~ 'fl_cl4',
                                       asv_num %in% c('asv4', 'asv31', 'asv23') & fraction == '3' ~'pa_cl1',
                                       asv_num %in% c('asv7', 'asv15') & fraction == '3' ~ 'pa_cl2',
                                       asv_num =='asv1'  & fraction == '3' ~ 'pa_cl3',
                                       asv_num %in% c('asv42', 'asv126', 'asv118', 'asv72', 'asv100', 'asv85') & fraction == '3' ~ 'pa_cl4',
                                       asv_num %in% c('asv84', 'asv116', 'asv25', 'asv28', 'asv182', 'asv27', 'asv11') & fraction =='3' ~ 'pa_cl5' ,
                                       asv_num %in% c('asv69', 'asv153', 'asv559', 'asv194', 'asv276', 'asv223', 'asv264') & fraction == '3' ~ 'pa_cl6',
                                       asv_num %in% c('asv317', 'asv200', 'asv311', 'asv511', 'asv113', 'asv43', 'asv225', 
                                                      'asv752', 'asv471', 'asv385') & fraction == '3' ~ 'pa_cl7' ,
                                       asv_num %in% c('asv192', 'asv163', 'asv49', 'asv302', 'asv179', 'asv77',
                                                      'asv219', 'asv105') & fraction == '3' ~ 'pa_cl8',
                                       asv_num %in% c('asv80', 'asv22', 'asv17') & fraction == '3' ~ 'pa_cl9')) |>
dplyr::filter(s4_cluster %in% c('fl_cl2', 'fl_cl3', 'pa_cl1', 'pa_cl2', 'pa_cl3', 'pa_cl5', 'pa_cl6', 'pa_cl8', 'pa_cl9')) |>
  ggplot(aes(decimal_date, wavelets_result)) +
  geom_vline(xintercept = 2010.24, linetype = 'dashed') +
  geom_vline(xintercept = 2010.52, linetype = 'dashed') +
  geom_vline(xintercept = 2012.48, linetype = 'dashed') +
  geom_line(aes(group = asv_num, color = order)) +
  facet_grid(vars(s4_cluster)) +
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6) +
  labs(x = 'Time', y = 'Wavelet coefficient: s4 (residuals)', color = 'Order') +
  scale_color_manual(values = palette_order_assigned_bloo) +
  #scale_linetype_manual(values = c("solid", "dashed")) + # Define linetypes
  theme_bw() +
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 7))

wavelets_s4_data <- wavelets_result_02 |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  bind_rows(wavelets_result_3) |>
  dplyr::filter(wavelets_transformation == 's4') |>
  dplyr::mutate(s4_cluster = case_when(asv_num %in% c('asv62', 'asv38', 'asv11', 'asv27') & fraction == '0.2' ~ 'fl_cl1' ,
                                       asv_num %in% c('asv7', 'asv15')  & fraction == '0.2'~ 'fl_cl2',
                                       asv_num =='asv1'  & fraction == '0.2' ~ 'fl_cl3',
                                       asv_num %in% c('asv249', 'asv114', 'asv237', 'asv282', 'asv563', 'asv555', 'asv178', 'asv58', 'asv17') & fraction == '0.2'~ 'fl_cl4',
                                       asv_num %in% c('asv4', 'asv31', 'asv23') & fraction == '3' ~'pa_cl1',
                                       asv_num %in% c('asv7', 'asv15') & fraction == '3' ~ 'pa_cl2',
                                       asv_num =='asv1'  & fraction == '3' ~ 'pa_cl3',
                                       asv_num %in% c('asv42', 'asv126', 'asv118', 'asv72', 'asv100', 'asv85') & fraction == '3' ~ 'pa_cl4',
                                       asv_num %in% c('asv84', 'asv116', 'asv25', 'asv28', 'asv182', 'asv27', 'asv11') & fraction =='3' ~ 'pa_cl5' ,
                                       asv_num %in% c('asv69', 'asv153', 'asv559', 'asv194', 'asv276', 'asv223', 'asv264')& fraction == '3' ~ 'pa_cl6',
                                       asv_num %in% c('asv317', 'asv200', 'asv311', 'asv511', 'asv113', 'asv43', 'asv225', 
                                                      'asv752', 'asv471', 'asv385') & fraction == '3' ~ 'pa_cl7' ,
                                       asv_num %in% c('asv192', 'asv163', 'asv49', 'asv302', 'asv179', 'asv77',
                                                      'asv219', 'asv105') & fraction == '3' ~ 'pa_cl8',
                                       asv_num %in% c('asv80', 'asv22', 'asv17') & fraction == '3' ~ 'pa_cl9')) |>
  dplyr::filter(s4_cluster %in% c('fl_cl2', 'fl_cl3', 'pa_cl1', 'pa_cl2', 'pa_cl3', 'pa_cl5', 'pa_cl6', 'pa_cl8', 'pa_cl9')) 

wavelets_s4_data$s4_cluster <- factor(wavelets_s4_data$s4_cluster, levels = c('fl_cl2', 'fl_cl3', 'pa_cl1', 'pa_cl2', 'pa_cl3', 'pa_cl5', 'pa_cl6', 'pa_cl8', 'pa_cl9'))
# 
# wavelets_s4_data |>
#   dplyr::filter(s4_cluster %in% c( 'pa_cl5', 'pa_cl6', 'pa_cl8', 'pa_cl9')) |>
#   ggplot(aes(decimal_date, interaction(asv_num, family))) +
#   geom_vline(xintercept = 2010.24, linetype = 'dashed') +
#   geom_vline(xintercept = 2010.52, linetype = 'dashed') +
#   geom_vline(xintercept = 2012.48, linetype = 'dashed') +
#   geom_point(aes(color = order, size = wavelets_result, group = s4_cluster)) +
#   facet_grid(vars(fraction), scales = 'free_y', labeller = labs_fraction) +
#   scale_size_continuous(range = c(0,10))+
#   # geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
#   #                                                       ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6) +
#   labs(x = 'Time', y = 'Wavelet coefficient: s4 (residuals)', color = 'Order') +
#   scale_color_manual(values = palette_order_assigned_bloo) +
#   #scale_linetype_manual(values = c("solid", "dashed")) + # Define linetypes
#   theme_bw() +
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 7))

# wavelets_result_02 |>
#   bind_rows(wavelets_result_3) |>
#   dplyr::group_by(asv_num, fraction) |>
#   dplyr::mutate(sample_num = row_number()) |>
#   dplyr::filter(wavelets_transformation == 's4') |>
#   ggplot(aes(sample_num, asv_num, fill = wavelets_result)) +
#   # geom_vline(xintercept = 2010.24, linetype = 'dashed') +
#   # geom_vline(xintercept = 2010.52, linetype = 'dashed') +
#   # geom_vline(xintercept = 2012.48, linetype = 'dashed') +
#   geom_tile(aes()) +
#   facet_wrap(vars(fraction), ncol = 1, scales = 'free_y') +
#   # geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
#   #                                                       ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6) +
#   labs(x = 'Time', y = 'Wavelet coefficient: s4 (residuals)', color = 'Order') +
#   scale_color_manual(values = palette_order_assigned_bloo) +
#   scale_linetype_manual(values = c("solid", "dashed")) + # Define linetypes
#   theme_bw() +
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 7))

## 
harbor_group_data <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(!asv_num %in% s4_fraction_smooth) |> # these taxa presented a different interannual tendency which is different from that only related to the harbor restoration
  dplyr::filter(fraction == '3') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  # dplyr::mutate(harbor_group = case_when(asv_num %in% c('asv105', 'asv113', 'asv192',
  #                                                       'asv194', 'asv22', 'asv223', 'asv225', #'asv264',
  #                                                       'asv276', 'asv311', 'asv317', 'asv43', 'asv49',
  #                                                       'asv69', 'asv752') ~ 'II',
  #                                        asv_num %in% c('asv116', 'asv118', 'asv126',
  #                                                       'asv182', 'asv28', 'asv72', 'asv84', 'asv85') ~ 'I',
  #                                        asv_num %in% c('asv80', 'asv77', 'asv264') ~  'III')) |>
  
  dplyr::mutate(harbor_group = case_when(
    # asv_num %in% c('asv42', 'asv126', 'asv118', ##groups based on clustering analyisis of s4
    #                                                     'asv72', 'asv100', 'asv25', 'asv28', #'asv264',
    #                                                     'asv182', 'asv27', 'asv11') ~ 'I',
    #                                      asv_num %in% c('asv69', 'asv153', 'asv559',
    #                                                     'asv194', 'asv276', 'asv223', 'asv264', 'asv85') ~ 'III',
    #                                      asv_num %in% c('asv217', 'asv200', 'asv311',
    #                                                     'asv511', 'asv113', 'asv43', 'asv225',
    #                                                     'asv752', 'asv471', 'asv385') ~  'II',
    #                                      asv_num %in% c('asv192', 'asv163', 'asv49',
    #                                                     'asv302', 'asv179', 'asv77', 'asv219',
    #                                                     'asv105', 'asv80', 'asv22') ~  'IV')) |>
    asv_num %in% c('asv42', 'asv126', 'asv118', 'asv72', 'asv100', 'asv85') & fraction == '3' ~ 'I',
    asv_num %in% c('asv84', 'asv116', 'asv25', 'asv28', 'asv182', 'asv27', 'asv11') & fraction =='3' ~ 'I' ,
    asv_num %in% c('asv69', 'asv153', 'asv559', 'asv194', 'asv276', 'asv223', 'asv264')& fraction == '3' ~ 'III',
    asv_num %in% c('asv317', 'asv200', 'asv311', 'asv511', 'asv113', 'asv43', 'asv225', 
                   'asv752', 'asv471', 'asv385') & fraction == '3' ~ 'II' ,
    asv_num %in% c('asv192', 'asv163', 'asv49', 'asv302', 'asv179', 'asv77',
                   'asv219', 'asv105', 'asv22') & fraction == '3' ~ 'II',
    asv_num %in% c('asv80') & fraction == '3' ~ 'III')) |>
  dplyr::filter(!is.na(harbor_group)) 

harbor_group <-harbor_group_data |> 
  dplyr::group_by(asv_num, date, harbor_group) |>
  dplyr::mutate(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format())+
  geom_vline(xintercept = as.POSIXct('2010-03-24', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_vline(xintercept = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_vline(xintercept = as.POSIXct('2012-06-09', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), 
                                                        xmax = as.POSIXct('2010-09-1', format = "%Y-%m-%d"), x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2011-07-01', format = "%Y-%m-%d"), 
                                                        xmax = as.POSIXct('2011-09-1', format = "%Y-%m-%d"), x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(y = 'Relative abundance (%)', x = 'Time', fill = 'Family')+
  geom_area(aes(group = asv_num, fill = family), position = 'stack')+
  facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 4/11,
        panel.grid = element_blank(), text = element_text(size = 6),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
          legend.key.size = unit(3, 'mm'))

harbor_group

##try to filter by max coefficent in s4----
wavelets_result_3 |>
  dplyr::filter(wavelets_transformation == 's4') |>
  dplyr::filter(!asv_num %in% c('asv4', 'asv31', 'asv23', 'asv7', 'asv15', 'asv1', 'asv17'))|>
  dplyr::group_by(asv_num) |>
dplyr::reframe(mean = mean(wavelets_result)) |>
  ungroup() |>
  reframe(mean = mean(mean))

asv_num_harbor <- wavelets_result_3 |>
  dplyr::filter(wavelets_transformation == 's4') |>
  dplyr::filter(!asv_num %in% c('asv4', 'asv31', 'asv23', 'asv7', 'asv15', 'asv1', 'asv17'))|>
  dplyr::group_by(asv_num) |>
  dplyr:: filter(any(wavelets_result > 0.257)) |>
  dplyr::select(asv_num) |>
  distinct()

harbor_group_rclr <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(!asv_num %in% s4_fraction_smooth) |> # these taxa presented a different interannual tendency which is different from that only related to the harbor restoration
  dplyr::filter(fraction == '3') |>
  dplyr::filter(abundance_type == 'rclr') |>
  # dplyr::mutate(harbor_group = case_when(asv_num %in% c('asv105', 'asv113', 'asv192', 
  #                                                       'asv194', 'asv22', 'asv223', 'asv225', #'asv264',
  #                                                       'asv276', 'asv311', 'asv317', 'asv43', 'asv49',
  #                                                       'asv69', 'asv752') ~ 'II',
  #                                        asv_num %in% c('asv116', 'asv118', 'asv126',
  #                                                       'asv182', 'asv28', 'asv72', 'asv84', 'asv85') ~ 'I',
  #                                        asv_num %in% c('asv80', 'asv77', 'asv264') ~  'III')) |> 
  
  dplyr::mutate(harbor_group = case_when(
    # asv_num %in% c('asv42', 'asv126', 'asv118', ##groups based on clustering analyisis of s4
    #                                                     'asv72', 'asv100', 'asv25', 'asv28', #'asv264',
    #                                                     'asv182', 'asv27', 'asv11') ~ 'I',
    #                                      asv_num %in% c('asv69', 'asv153', 'asv559',
    #                                                     'asv194', 'asv276', 'asv223', 'asv264', 'asv85') ~ 'III',
    #                                      asv_num %in% c('asv217', 'asv200', 'asv311',
    #                                                     'asv511', 'asv113', 'asv43', 'asv225',
    #                                                     'asv752', 'asv471', 'asv385') ~  'II',
    #                                      asv_num %in% c('asv192', 'asv163', 'asv49',
    #                                                     'asv302', 'asv179', 'asv77', 'asv219',
    #                                                     'asv105', 'asv80', 'asv22') ~  'IV')) |>
    asv_num %in% c('asv42', 'asv126', 'asv118', 'asv72', 'asv100', 'asv85') & fraction == '3' ~ 'I', ##groups based on clustering analyisis of s4
    asv_num %in% c('asv84', 'asv116', 'asv25', 'asv28', 'asv182', 'asv27', 'asv11') & fraction =='3' ~ 'I' ,
    asv_num %in% c('asv69', 'asv153', 'asv559', 'asv194', 'asv276', 'asv223', 'asv264')& fraction == '3' ~ 'III', # late responsive 
    asv_num %in% c('asv317', 'asv200', 'asv311', 'asv511', 'asv113', 'asv43', 'asv225', 
                   'asv752', 'asv471', 'asv385') & fraction == '3' ~ 'II' , #first responsive
    asv_num %in% c('asv192', 'asv163', 'asv49', 'asv302', 'asv179', 'asv77',
                   'asv219', 'asv105', 'asv22') & fraction == '3' ~ 'II', #first responsives
    asv_num %in% c('asv80') & fraction == '3' ~ 'III')) |> ## it is similar to late responsives
  dplyr::filter(!is.na(harbor_group)) |>
  dplyr::filter(asv_num %in% asv_num_harbor$asv_num) |>
  dplyr::filter(!is.na(harbor_group)) |> 
  dplyr::group_by(asv_num, date, harbor_group) |>
  dplyr::mutate(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abundance_value))+
  #scale_y_continuous(labels = percent_format())+
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
  geom_area(aes(group = asv_num, fill = family), position = 'stack')+
  facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 4/11,
        panel.grid = element_blank(), text = element_text(size = 6),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

harbor_group_rclr

# asv_tab_all_bloo_z_tax_harbor |>
#   dplyr::filter(!asv_num %in% s4_fraction_smooth) |> # these taxa presented a different interannual tendency which is different from that only related to the harbor restoration
#   dplyr::filter(fraction == '3') |>
#   dplyr::filter(abundance_type == 'rclr') |>
#   dplyr::mutate(harbor_group = case_when(asv_num %in% c('asv105', 'asv113', 'asv192',
#                                                         'asv194', 'asv22', 'asv223', 'asv225', 'asv264',
#                                                         'asv276', 'asv311', 'asv317', 'asv43', 'asv49',
#                                                         'asv69', 'asv752') ~ 'II',
#                                          asv_num %in% c('asv116', 'asv118', 'asv126',
#                                                         'asv182', 'asv28', 'asv72', 'asv84', 'asv85') ~ 'I',
#                                          asv_num %in% c('asv80', 'asv77') ~  'III')) |>
#   dplyr::mutate(harbor_group = as.factor(harbor_group)) |>
#   dplyr::filter(!is.na(harbor_group)) |> 
#   # dplyr::group_by(asv_num, date, harbor_group) |>
#   # dplyr::mutate(abund_max = sum(abundance_value)) |>
#   ggplot(aes(date, abundance_value))+
#   #scale_y_continuous(labels = percent_format())+
#   geom_vline(xintercept = as.POSIXct('2010-03-24', format = "%Y-%m-%d"), linetype = 'dashed') +
#   geom_vline(xintercept = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), linetype = 'dashed') +
#   geom_vline(xintercept = as.POSIXct('2012-06-09', format = "%Y-%m-%d"), linetype = 'dashed') +
#   geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), 
#                                                         xmax = as.POSIXct('2010-08-31', format = "%Y-%m-%d"), x=NULL, y=NULL,
#                                                         ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
#   geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2011-07-01', format = "%Y-%m-%d"), 
#                                                         xmax = as.POSIXct('2011-08-31', format = "%Y-%m-%d"), x=NULL, y=NULL,
#                                                         ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   labs(y = 'rCLR', x = 'Time', color = 'Family')+
#   geom_line(aes(group = asv_num, color = family), position = 'stack')+
#   facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'right',
#         panel.grid = element_blank(), text = element_text(size = 6),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)))

## upload photos from the harbor restoration period ---


# Load and prepare your image
# #install.packages('jpeg')
# img_1 <- jpeg::readJPEG("data/env_data/harbour_restoration_data/photo_1.jpg")
# img_grob <- grid::rasterGrob(img_1, interpolate = TRUE)

# Arrange plots and image using cowplot
#library(cowplot) 
# Create the composition of plots
composition <- plot_grid(harbor_group, align = "hv", nrow = 1) + 
  draw_image("data/env_data/harbour_restoration_data/photo_1.jpg", x = 0.725, y = 0.66, width = 0.25, height = 0.25) +
  draw_image("data/env_data/harbour_restoration_data/photo_2.jpg", x = 0.725, y = 0.435, width = 0.25, height = 0.25) +
  draw_image("data/env_data/harbour_restoration_data/photo_3.jpg", x = 0.725, y = 0.2, width = 0.25, height = 0.25)

# Add labels to each part of the composition
composition_with_labels <- composition +
  annotate("text", x = 0.03, y = 0.88, label = "A", size = 3) +
  annotate("text", x = 0.74, y = 0.88, label = "B", size = 3,  colour = "white") +
  annotate("text", x = 0.74, y = 0.66, label = "C", size = 3, colour = "white") +
  annotate("text", x = 0.74, y = 0.41, label = "D", size = 3, colour = "white")

# Display the composition with labels
print(composition_with_labels)

# Save the composition as a PDF
ggsave(filename = 'harbor_restoration_plot_photo_v2.pdf',
       plot = composition_with_labels,
       path = 'Results/Figures/',
       width = 180, height = 200, units = 'mm')

# Create the composition of plots
composition <- plot_grid(harbor_group_rclr, align = "hv", nrow = 1) + 
  draw_image("data/env_data/harbour_restoration_data/photo_1.jpg", x = 0.725, y = 0.66, width = 0.25, height = 0.25) +
  draw_image("data/env_data/harbour_restoration_data/photo_2.jpg", x = 0.725, y = 0.435, width = 0.25, height = 0.25) +
  draw_image("data/env_data/harbour_restoration_data/photo_3.jpg", x = 0.725, y = 0.2, width = 0.25, height = 0.25)

# Display the composition
print(composition)

# Add labels to each part of the composition
composition_with_labels <- composition +
  annotate("text", x = 0.03, y = 0.88, label = "A", size = 3) +
  annotate("text", x = 0.74, y = 0.88, label = "B", size = 3,  colour = "white") +
  annotate("text", x = 0.74, y = 0.66, label = "C", size = 3, colour = "white") +
  annotate("text", x = 0.74, y = 0.41, label = "D", size = 3, colour = "white")

# Display the composition with labels
print(composition_with_labels)

# Save the composition as a PDF
ggsave(filename = 'harbor_restoration_plot_photo_rclr_v2.pdf',
       plot = composition_with_labels,
       path = 'Results/Figures/',
       width = 180, height = 200, units = 'mm')
