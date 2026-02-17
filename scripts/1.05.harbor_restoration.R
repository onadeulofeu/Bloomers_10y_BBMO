# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# packages ----
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

asv_tab_bbmo_10y_w_rar_harbor <- asv_tab_bbmo_10y_w_rar |>
  left_join(m_bbmo_10y) |>
  dplyr::mutate(harbor_restoration = case_when(  date < '2010-03-24' ~ 'pre_perturbation',
                                                 date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
                                                 date >= '2012-06-09' ~ 'post_perturbation')) 
  
asv_tab_all_bloo_z_tax_harbor <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) ##remove those taxa that are not real bloomers

asv_tab_all_bloo_z_tax_harbor$harbor_restoration <- factor(asv_tab_all_bloo_z_tax_harbor$harbor_restoration, 
                                                           levels = c('pre_perturbation',
                                                           'perturbation',
                                                           'post_perturbation'))

labs_perturbation <- as_labeller(c( 'pre_perturbation' = 'Pre' ,
                                     'perturbation' = 'During',
                                    'post_perturbation' = 'Post'))

labs_season <- as_labeller(c( 'winter' = 'Winter' ,
                                    'spring' = 'Spring',
                                    'summer' = 'Summer',
                              'autumn' = 'Autumn'))

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
shannon_rar_m |>
  group_by(fraction, harbor_restoration) |>
  dplyr::reframe(n = n())

## add stats
data_f <- shannon_rar_m |>
  dplyr::filter(fraction == '0.2')

## add stats
data_f_3 <- shannon_rar_m |>
  dplyr::filter(fraction == '3')

data_f |>
  colnames()

data_f |>
  distinct(harbor_restoration)

## check normality 
shapiro.test(as.numeric(shannon_rar_m |>
                          dplyr::filter(fraction == '0.2') %$%
                          diversity_shannon)) # => p-value = 0.0172 (NO NORMALITY)

shapiro.test(as.numeric(shannon_rar_m |>
                          dplyr::filter(fraction == '3') %$%
                          diversity_shannon)) # => p-value = 0.0172 (NO NORMALITY)

ggqqplot(as.numeric(shannon_rar_m |>
                      dplyr::filter(fraction == '0.2') %$%
                      diversity_shannon))

## check homocedasticity 
leveneTest(diversity_shannon ~ harbor_restoration, data = data_f) #p-value 0.844 > 0.05: Variances are homogeneous (equal), and the assumption for ANOVA is met

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value Pr(>F)
# group   2  0.1698  0.844
#       117                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

leveneTest(diversity_shannon ~ harbor_restoration, data = data_f_3)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value Pr(>F)
# group   2  0.8704 0.4215
#       114  

## Since we don't have normality but we have homocedasticity we still need to change to a non parametric test kruskal.test

# Perform one-way ANOVA to compare the means of Rho across the three groups
## anova_result <- aov(rho ~ fraction_causal, data = data_ed2_f)  # Replace with your actual group variable
kruskal_result <- kruskal.test(diversity_shannon ~ harbor_restoration, data = data_f) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_f$diversity_shannon,                            # The values you are comparing
  g = as.factor(data_f$harbor_restoration),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)  # no significant differences in the FL fraction

# Perform one-way ANOVA to compare the means of Rho across the three groups
## anova_result <- aov(rho ~ fraction_causal, data = data_ed2_f)  # Replace with your actual group variable
kruskal_result <- kruskal.test(diversity_shannon ~ harbor_restoration, data = data_f_3) 

# View the ANOVA summary
##summary(anova_result)
summary(kruskal_result)

# Perform Dunn's test correctly
dunn_result <- dunn.test(
  x = data_f_3$diversity_shannon,                            # The values you are comparing
  g = as.factor(data_f_3$harbor_restoration),    # The grouping variable
  method = 'bonferroni'                          # Bonferroni correction
)

# Print the result
print(dunn_result)  # no significant differences in the FL fraction

shannon_harbor_rar_plot <- shannon_rar_m |>
  ggplot(aes(harbor_restoration, diversity_shannon))+
  geom_point(aes(color = season), position = position_jitter(width = 0.15))+
  scale_x_discrete(labels = labs_perturbation)+
  labs(y = 'Shannon diversity', x = 'Harbor restoration', color = 'Season')+
  scale_color_manual(values = palette_seasons_4, labels = labs_season)+
  facet_wrap(vars(fraction), labeller = labs_fraction)+
  coord_flip()+
  geom_violin(aes(group = harbor_restoration), alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  theme_bw()+
  theme(text = element_text(size = 10),
        strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        legend.position = 'bottom')

shannon_harbor_rar_plot 
# 
# ggsave(shannon_harbor_rar_plot, file = 'results/figures/shannon_harbor_rar_plot.pdf',
#        width = 110, height = 110, units = 'mm')

shannon_harbor_rar_plot <- shannon_rar_m |>
  ggplot(aes(fraction, diversity_shannon))+
  geom_point( position = position_jitter(width = 0.15), aes(shape = fraction))+#aes(color = season),
  scale_x_discrete(labels = labs_fraction)+
  labs(y = 'Shannon diversity', x = 'Harbor restoration', color = 'Season')+
  #scale_color_manual(values = palette_seasons_4, labels = labs_season)+
  facet_wrap(vars(harbor_restoration), labeller = labs_perturbation)+
  scale_shape_discrete(labels = labs_fraction)+
  coord_flip()+
  geom_violin(aes(group = fraction), alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  theme_bw()+
  theme(text = element_text(size = 20),
        strip.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        legend.position = 'bottom')

shannon_harbor_rar_plot 
# 
# ggsave(shannon_harbor_rar_plot, 
#        file = 'results/figures/poster_svg_format/shannon_harbor_rar_plot.svg',
#        width = 260, height = 120, units = 'mm')

## exploring the community affected by the restoration -----
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

bloo_types_summary <- bloo_all_types_summary_tb_tax |>
  dplyr::mutate(asv_num_f = asv_num) |>
  dplyr::mutate(fraction = as.character(fraction))
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

## for the ISME conference -----
blooming_asvs_increasing_over_the_years <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::mutate(famil_asv = paste0(family,' ', asv_num)) |>
  left_join(bloo_types_summary) |>
  dplyr::mutate(fraction_asv_num = paste0(asv_num, '.', fraction)) |>
  #dplyr::filter(frequency != 'seasonal') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(asv_num %in% c('asv1', 'asv15', 'asv7', 'asv17', 'asv23', 'asv31', 'asv4')) |>
  ggplot(aes(date, abundance_value, group = asv_num))+
  facet_wrap(vars(fraction), scales = 'fixed')+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_line(aes(group = fraction_asv_num, color = family), linewidth = 0.2, alpha = 0.6)+
  geom_smooth(aes(group = fraction_asv_num, color = family, fill = family,
                  linetype = fraction), method = 'loess', linewidth = 0.8, alpha = 0.5)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #scale_y_continuous(labels = percent_format())+
  scale_linetype_discrete(labels = labs_fraction)+
  labs(x = 'Time', y = 'rCLR', color = 'Family', fill = 'Family', linetype = 'Fraction')+
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
  theme(#strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right', 
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_blank())

blooming_asvs_increasing_over_the_years

# ggsave(blooming_asvs_increasing_over_the_years, filename = 'blooming_asvs_increasing_over_the_years.svg',
#        path = 'Results/Figures/poster_svg_format/',
#        width = 300, height = 180, units = 'mm')

blooming_asvs_decreasing_over_the_years <- asv_tab_all_bloo_z_tax_harbor |>
  dplyr::mutate(famil_asv = paste0(family,' ', asv_num)) |>
  left_join(bloo_types_summary) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(asv_num %in% c('asv28')) |>
  ggplot(aes(date, abundance_value, group = asv_num))+
  facet_wrap(vars(interaction(family_f, asv_num_f)), scales = 'free',  ncol = 3)+
  geom_line(aes(group = fraction, color = fraction), linewidth = 0.3, alpha = 0.6)+
  geom_smooth(aes(group = fraction, color = fraction), method = 'loess', linewidth = 0.8)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  labs(x = 'Time', y = 'rCLR', color = 'Fraction')+
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

## upload data
wavelets_result_tibble_tax_3_biased <- read.csv('data/wavelets_analysis/wavelets_result_ed_tibble_tax_3_biased_red.csv')
wavelets_result_tibble_tax_02_biased <- read.csv('data/wavelets_analysis/wavelets_result_ed_tibble_tax_02_biased_red.csv')

wavelets_result_02 <- wavelets_result_tibble_tax_02_biased
wavelets_result_3 <- wavelets_result_tibble_tax_3_biased

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

harbor_group <- harbor_group_data |> 
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
        aspect.ratio = 4/9,
        panel.grid = element_blank(), text = element_text(size = 18),
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
  draw_image("data/env_data/harbour_restoration_data/blanes_harbor_icgc/Blanes_2010.jpg", x = 0.725, y = 0.66, width = 0.25, height = 0.25) +
  draw_image("data/env_data/harbour_restoration_data/blanes_harbor_icgc/Blanes_2011.jpg", x = 0.725, y = 0.435, width = 0.25, height = 0.25) +
  draw_image("data/env_data/harbour_restoration_data/blanes_harbor_icgc/Blanes_2013.jpg", x = 0.725, y = 0.2, width = 0.25, height = 0.25)

# Add labels to each part of the composition
composition_with_labels <- composition +
  annotate("text", x = 0.03, y = 0.88, label = "A", size = 3) +
  annotate("text", x = 0.74, y = 0.88, label = "B", size = 3,  colour = "white") +
  annotate("text", x = 0.74, y = 0.66, label = "C", size = 3, colour = "white") +
  annotate("text", x = 0.74, y = 0.41, label = "D", size = 3, colour = "white")

# Display the composition with labels
print(composition_with_labels)

# Save the composition as a PDF
ggsave(filename = 'harbor_restoration_plot_photo_v3.pdf',
       plot = composition_with_labels,
       path = 'Results/Figures/',
       width = 180, height = 200, units = 'mm')

# save composition as .svg for a poster
ggsave(filename = 'harbor_restoration_plot_photo_v3.svg',
       plot = harbor_group_rclr,
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
# ggsave(filename = 'harbor_restoration_plot_photo_rclr_v2.pdf',
#        plot = composition_with_labels,
#        path = 'Results/Figures/poster_svg_format/',
#        width = 180, height = 200, units = 'mm')

## Harbor restoration figure with the results obtained from the continuous wavelet transformation -----
asv_tab_all_bloo_z_tax |>
  left_join(bloo_all_types_summary_tb, by = c('asv_num', 'fraction')) |>
  #left_join(summary_types_of_blooms, by = c('asv_num', 'fraction')) |>
  dplyr::filter(asv_num %in% bloo_02$value & fraction == '0.2' |
                  asv_num %in% bloo_3$value & fraction == '3' ) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '3') |>
  dplyr::mutate(harbor_groups = case_when(type_of_bloomer == 'Seasonal' ~ 'Seasonal',
                                          type_of_bloomer %in% c('Year-to-Year', 'Inter-Annual') ~ 'Siginificant change',
                                          type_of_bloomer %in% c('No-significant Periodicity') ~ 'No-significant change')) |>
  #dplyr::filter(type_of_bloomer %in% c('Year-to-Year', 'Inter-Annual', 'No-significant Periodicity')) |> 
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
  geom_area(aes(group = asv_num_f, fill = family_f), position = 'stack')+
  facet_wrap(vars(harbor_groups), ncol = 1, scales = 'fixed')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 4/11,
        panel.grid = element_blank(), text = element_text(size = 6),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

## i try a heatmap
asv_tab_all_bloo_z_tax |>
  left_join(bloo_all_types_summary_tb, by = c('asv_num', 'fraction')) |>
  #left_join(summary_types_of_blooms, by = c('asv_num', 'fraction')) |>
  dplyr::filter(asv_num %in% bloo_02$value & fraction == '0.2' |
                  asv_num %in% bloo_3$value & fraction == '3' ) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(fraction == '3') |>
  dplyr::mutate(harbor_groups = case_when(type_of_bloomer == 'Seasonal' ~ 'Seasonal',
                                          type_of_bloomer %in% c('Year-to-Year', 'Inter-Annual') ~ 'Siginificant change',
                                          type_of_bloomer %in% c('No-significant Periodicity') ~ 'No-significant change')) |>
  #dplyr::filter(type_of_bloomer %in% c('Year-to-Year', 'Inter-Annual', 'No-significant Periodicity')) |> 
  ggplot(aes(date, y = asv_num, fill = family))+
  #scale_y_continuous(labels = percent_format())+
  geom_vline(xintercept = as.POSIXct('2010-03-24', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_vline(xintercept = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_vline(xintercept = as.POSIXct('2012-06-09', format = "%Y-%m-%d"), linetype = 'dashed') +
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2010-07-01', format = "%Y-%m-%d"), 
                                                        xmax = as.POSIXct('2010-09-1', format = "%Y-%m-%d"), x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = as.POSIXct('2011-07-01', format = "%Y-%m-%d"), 
                                                        xmax = as.POSIXct('2011-09-1', format = "%Y-%m-%d"), x=NULL, y=NULL,
                                                        ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  #scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(y = 'Relative abundance (%)', x = 'Time', fill = 'Family')+
  #geom_area(aes(group = asv_num_f, fill = family_f), position = 'stack')+
  geom_tile(aes(fill = as.numeric(abundance_value)))+
  facet_wrap(vars(harbor_groups), ncol = 1, scales = 'free_y')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 4/11,
        panel.grid = element_blank(), text = element_text(size = 6),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

### Different try for the harbor restoration plot ----
# #### I will keep any significant wavelet power 
# wavelet_results_all <- wavelet_results_tibble_3 |>
#   bind_rows(wavelet_results_tibble_02) |>
#   left_join(occurrence_bloo_bbmo_wav) |>
#   left_join(bloo_taxonomy)
# 
# mean_power_sig <- wavelet_results_all |>
#   dplyr::filter(pvalue_power_average < 0.05) |>
#   ungroup() |>
#   reframe(mean = mean(power_average))
# 
# wavelet_results_category <- wavelet_results_all |>
#   dplyr::filter(pvalue_power_average <= 0.05) |> # & power_average > mean_power_sig$mean
#   dplyr::group_by(asv_num, fraction) |>
#   #slice_max(order_by = power_average, n = 1 ) |>
#   dplyr::mutate(type_of_bloomer = case_when(pvalue_power_average < 0.05 #& power_average > 0.25
#                                             & period <= 4 ~ 'short',
#                                             pvalue_power_average < 0.05 & # power_average > 0.25 &
#                                               period > 4 &
#                                               period < 8 ~ 'half-yearly',
#                                             pvalue_power_average < 0.05 &# power_average > 0.25 &
#                                               period > 8  &
#                                               period == 8 |
#                                               period <= 16 ~ 'seasonal',
#                                             pvalue_power_average < 0.05 #& power_average > 0.25 
#                                             & period > 16 
#                                             & period < 32 ~ 'year-to-year', 
#                                             pvalue_power_average < 0.05 #& power_average > 0.25 
#                                             & period == 32 |
#                                               period >32 ~ 'Inter-annual'
#   )) |>
#   distinct(type_of_bloomer, fraction, asv_num) 
# 
# wavelet_results_all  <- wavelet_results_all |>
#   left_join(wavelet_results_category, relationship = 'many-to-many', by = c('asv_num', 'fraction')) |>
#   dplyr::mutate(type_of_bloomer = case_when(is.na(type_of_bloomer) ~ 'no-significant periodicity',
#                                             !is.na(type_of_bloomer) ~ type_of_bloomer)) 
# 
# wavelet_results_all$type_of_bloomer <- factor(wavelet_results_all$type_of_bloomer,
#                                               levels = c('short', 'half-yearly', 'seasonal', 'year-to-year',
#                                                          'Inter-annual', 'no-significant periodicity'),
#                                               labels = c('Fine Scale', 'Half-Yearly', 'Seasonal',
#                                                          'Year-to-Year', 'Inter-Annual', 'No-significant Periodicity'))
# 
# wavelet_results_all$occurrence_category <- factor(wavelet_results_all$occurrence_category,
#                                                   levels = c('broad', 'intermediate', 'narrow'),
#                                                   labels = c('Broad', 'Intermediate', 'Narrow'))
# 
# bloo_all_types_summary_tb <- bloo_all_types_summary_tb |>
#   dplyr::mutate(fraction = as.character(fraction))
# 
# wavelet_results_category_harbor_affected <- wavelet_results_category |>
#   dplyr::filter(fraction == '3') |>
#   distinct(asv_num)
# 
# wavelet_results_category |>
#   dplyr::filter(fraction == '3') |>
#   distinct(asv_num, type_of_bloomer) |>
#   distinct(type_of_bloomer)
# 
# asv_tab_all_bloo_z_tax |>
#   #dplyr::filter(asv_num %in% wavelet_results_category_harbor_affected$asv_num) |>
#   #left_join(summary_types_of_blooms, by = c('asv_num', 'fraction')) |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   dplyr::filter(abundance_type == 'rclr') |>
#   dplyr::filter(fraction == '3') |>
#   left_join(wavelet_results_category, by = c('asv_num', 'fraction'), relationship = "many-to-many") |>
#   dplyr::mutate(harbor_groups = case_when(type_of_bloomer == 'seasonal' ~ 'Seasonal',
#                                           type_of_bloomer %in% c('year-to-year', 'Inter-annual') ~ 'Siginificant change',
#                                           type_of_bloomer %in% c('No-significant Periodicity') ~ 'No-significant change')) |>
#   dplyr::filter(harbor_groups == 'Siginificant change') |>
#   dplyr::mutate(abund_max = sum(abundance_value)) |>
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
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   labs(y = 'rCLR', x = 'Time', fill = 'Family')+
#   geom_area(aes(group = asv_num, fill = family), position = 'stack')+
#   #facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
#   facet_wrap(vars(asv_num),  scales = 'free_y')+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         aspect.ratio = 4/11,
#         panel.grid = element_blank(), text = element_text(size = 6),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# asv_tab_all_bloo_z_tax |>
#   #dplyr::filter(asv_num %in% wavelet_results_category_harbor_affected$asv_num) |>
#   #left_join(summary_types_of_blooms, by = c('asv_num', 'fraction')) |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   dplyr::filter(abundance_type == 'rclr') |>
#   dplyr::filter(fraction == '3') |>
#   left_join(wavelet_results_category, by = c('asv_num', 'fraction'), relationship = "many-to-many") |>
#   dplyr::mutate(harbor_groups = case_when(type_of_bloomer == 'seasonal' ~ 'Seasonal',
#                                           type_of_bloomer %in% c('year-to-year', 'Inter-annual') ~ 'Siginificant change',
#                                           type_of_bloomer %in% c('No-significant Periodicity') ~ 'No-significant change')) |>
#   dplyr::filter(harbor_groups == 'Seasonal') |>
#   dplyr::mutate(abund_max = sum(abundance_value)) |>
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
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   labs(y = 'rCLR', x = 'Time', fill = 'Family')+
#   geom_area(aes(group = asv_num, fill = family), position = 'stack')+
#   #facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
#   facet_wrap(vars(asv_num),  scales = 'free_y', ncol = 2)+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         aspect.ratio = 4/11,
#         panel.grid = element_blank(), text = element_text(size = 6),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# ## plot by groups after checking in detail  -----
# wavelet_results_category_3 <- wavelet_results_category |>
#   dplyr::filter(fraction == '3') 
# 
# wavelet_results_category_3 |>
#   ungroup() |>
#   distinct(asv_num) # I have 38 ASVs with a significant tendency

## I separate those taxa that were positively affected by harbor restoration and negatively affected ----
bloo_all_types_summary_tax_3 <- bloo_all_types_summary_tax |>
  dplyr::filter(fraction == '3') |>
  dplyr::mutate(fraction = as.character(fraction))

# positive_negative_effect_rel <- asv_tab_all_bloo_z_tax |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   dplyr::mutate(harbor_restoration = case_when(  date < '2010-03-24' ~ 'pre_perturbation',
#                                                  date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
#                                                  date >= '2012-06-09' ~ 'after_perturbation')) |>
#   dplyr::group_by(harbor_restoration, abundance_type, asv_num) |>
#   dplyr::filter(fraction == '3') |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   dplyr::reframe(mean_abund = mean(abundance_value)) |>
#   pivot_wider(names_from = 'harbor_restoration', values_from = starts_with('mean')) |> 
#   dplyr::mutate(harbor_effect_rel = case_when(
#     pre_perturbation > perturbation ~ 'negative',
#     pre_perturbation < perturbation ~ 'positive',
#     after_perturbation > perturbation |
#       after_perturbation > pre_perturbation ~ 'opportunistic',
#     TRUE ~ 'neutral'
#   )) |> 
#   distinct(asv_num, harbor_effect_rel)

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

# positive_negative_effect_rel |>
#   bind_cols(positive_negative_effect) |>
#   dplyr::filter(harbor_effect_rel != harbor_effect_rclr) # there is only one taxa that has different effects when calculating them using rCLR or relative abundance
# 
# asv_tab_all_bloo_z_tax |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   dplyr::filter(asv_num == 'asv85') |>
#   dplyr::filter(fraction == '3') |>
#   dplyr::filter(abundance_type != 'pseudoabundance') |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   slice_max(order_by = abundance_value, n = 1) ## the bloom of this taxa happens when the harbor restoration stops during summer. I belive the rCLR results
# 
# asv_tab_all_bloo_z_tax |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   dplyr::filter(asv_num == 'asv85') |>
#   dplyr::filter(fraction == '3') |>
#   dplyr::filter(abundance_type != 'pseudoabundance') |>
#   ggplot(aes(date, abundance_value))+
#   geom_line()+
#   facet_wrap(vars(abundance_type), scales = 'free')
#   
# bloo_all_types_summary_tax_3|>
#   dplyr::filter(asv_num == 'asv85') 

## are there statistically significant differences? -----
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
  #scale_color_manual(values = c('#4cb76a', '#2466BF'))+
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

## viruses differential abundances during the harbor restoration ----
m_vir_harbor <- m_vir_tb |>
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d")) |> 
  dplyr::mutate(harbor_restoration = case_when(
    date < '2010-03-24' ~ 'pre_perturbation',
    date >= '2010-03-24' & date < '2012-06-09' ~ 'perturbation',
    date >= '2012-06-09' ~ 'after_perturbation'
  ))

m_vir_harbor$harbor_restoration <- factor(m_vir_harbor$harbor_restoration, 
                                          levels = c('pre_perturbation', 'perturbation', 'after_perturbation') )

library(dplyr)
library(rstatix)

# Pairwise Wilcoxon test
diff_abundance_wilcox <- m_vir_harbor |>
  dplyr::filter(!is.na(total_vlp)) |>
  pairwise_wilcox_test(total_vlp ~ year, p.adjust.method = "BH")

print(diff_abundance_wilcox)

plot_virus_abund_harbor <- m_vir_harbor |>
  ggplot(aes( as.character(year), total_vlp))+
  geom_point(position = position_jitter(width = 0.15), size = 0.25)+
  scale_x_discrete()+
  geom_boxplot(alpha = 0.2)+
  #scale_x_discrete(labels = c('pre_perturbation' = 'Pre', 'perturbation' = 'Perturbation', 'after_perturbation' = 'Post'))+
  labs(x = '', y = 'virus/mL', color = '', shape = '')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("2010", "2012"), c("2011", "2012"))) # Add significance comparisons to boxplot

plot_virus_abund_harbor

# ggsave(filename = 'plot_virus_abund_harbor.pdf', 
#               plot = plot_virus_abund_harbor,
#               path = 'Results/Figures/',
#               width = 100, height = 88, units = 'mm')

diff_abundance_wilcox <- m_vir_harbor |>
  dplyr::filter(!is.na(HNF_Micro)) |>
  pairwise_wilcox_test(HNF_Micro ~ year, p.adjust.method = "BH")

print(diff_abundance_wilcox)
## non-significant the HNF 

plot_HNF_abund_harbor <- m_vir_harbor |>
  ggplot(aes( as.character(year), HNF_Micro))+
  geom_point(position = position_jitter(width = 0.15))+
  geom_boxplot(alpha = 0.2, notch = F)+
  scale_x_discrete(labels = c('pre_perturbation' = 'Pre', 'perturbation' = 'Perturbation', 'after_perturbation' = 'Post'))+
  labs(x = '', y = 'cells/mL', color = '', shape = '')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("2010", "2012"), c("2011", "2012"))) # Add significance comparisons to boxplot

plot_HNF_abund_harbor

## new categories related to harbor restoration -----
## relative abundance ----
asv_tab_all_bloo_z_tax_sig <- asv_tab_all_bloo_z_tax |>
  #dplyr::filter(asv_num %in% wavelet_results_category_harbor_affected$asv_num) |>
  #left_join(summary_types_of_blooms, by = c('asv_num', 'fraction')) |>
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

## rCLR plots ---
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

asv_tab_all_bloo_z_tax_sig$type_of_bloomer |>
  unique()

asv_tab_all_bloo_z_tax_sig$harbor_effect <- factor(asv_tab_all_bloo_z_tax_sig$harbor_effect, 
                                                   levels = c('negative', 'positive'), labels = c('Harbor restoration negative effect', 'Harbor restoration positive effect'))

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
  geom_area(aes(group = family_f, fill = family_f), position = 'stack')+
  facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
  #facet_grid(type_of_bloomer~harbor_effect,  scales = 'free_y')+
  facet_wrap(vars(interaction(type_of_bloomer, harbor_effect)),  scales = 'free_y', ncol = 1)+
  guides(fill = guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        #aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(4, 'mm'))

harbor_pa_stats_plot
# 
#ggsave(filename = 'harbor_restoration_plot_photo_rclr_stats_v3.pdf',
       # plot = harbor_pa_stats_plot,
       # path = 'Results/Figures/',
       # width = 180, height = 200, units = 'mm')

harbor_pa_stats_plot <- asv_tab_all_bloo_z_tax_sig |>
  dplyr::filter(significant == 'TRUE') |>
  dplyr::mutate(harbor_effect = case_when(harbor_effect == 'positive' ~ 'Positive',
                                          harbor_effect == 'negative' ~ 'Negative')) |>
  dplyr::group_by(date, fraction, order_f,  abundance_type, type_of_bloomer, harbor_effect) |>
  dplyr::reframe(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abund_max))+
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
  scale_fill_manual(values = palette_order_assigned_bloo)+
  labs(y = 'rCLR', x = 'Time', fill = 'Order')+
  geom_area(aes(group = order_f, fill = order_f), position = 'stack')+
  #facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
  #facet_grid(type_of_bloomer~harbor_effect,  scales = 'free_y')+
  facet_wrap(vars(interaction(type_of_bloomer,  harbor_effect)),  scales = 'free_y', ncol = 1)+
  guides(fill = guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        #aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(4, 'mm'))

harbor_pa_stats_plot

# ggsave(filename = 'harbor_restoration_plot_photo_rclr_stats_v3.pdf',
#        plot = harbor_pa_stats_plot,
#        path = 'Results/Figures/',
#        width = 180, height = 200, units = 'mm')

asv_tab_all_bloo_z_tax_sig |>
  dplyr::filter(significant == 'FALSE') |>
  distinct(asv_num)

harbor_pa_stats_plot_unsignificant <- asv_tab_all_bloo_z_tax_sig |>
  dplyr::filter(significant == 'FALSE') |>
  dplyr::mutate(harbor_effect = case_when(harbor_effect == 'positive' ~ 'Positive',
                                          harbor_effect == 'negative' ~ 'Negative')) |>
  dplyr::group_by(date, fraction, order_f, family_f, abundance_type, type_of_bloomer, harbor_effect) |>
  dplyr::reframe(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abund_max))+
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
  geom_area(aes(group = family_f, fill = family_f), position = 'stack')+
  #facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
  #facet_grid(type_of_bloomer~harbor_effect,  scales = 'free_y')+
  facet_wrap(vars(interaction(type_of_bloomer, harbor_effect)),  scales = 'free_y', ncol = 1)+
  guides(fill = guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        #aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 10),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

harbor_pa_stats_plot_unsignificant

# ggsave(filename = 'harbor_restoration_plot_photo_rclr_stats_unsignificant.pdf',
#        plot = harbor_pa_stats_plot_unsignificant,
#        path = 'Results/Figures/',
#        width = 180, height = 100, units = 'mm')
  
## unfiltered by the significance in the statistics test -----
asv_tab_all_bloo_z_tax_sig |>
  dplyr::mutate(type_of_bloomer = case_when(type_of_bloomer == 'Inter-Annual' ~ 'Chaotic',
                                            type_of_bloomer == "Year-to-Year" ~ 'Chaotic',
                                            type_of_bloomer ==  "No-significant Periodicity" ~ 'Chaotic',
                                            type_of_bloomer == 'Recurrent' ~ 'Seasonal',
                                            TRUE ~ type_of_bloomer)) |>
  #dplyr::filter(significant == 'no') |>
  dplyr::group_by(date, fraction, order_f, family_f, abundance_type, type_of_bloomer, harbor_effect) |>
  dplyr::reframe(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abund_max))+
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
  geom_area(aes(group = family_f, fill = family_f), position = 'stack')+
  #facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
  facet_wrap(type_of_bloomer~harbor_effect,  scales = 'free_y', ncol = 1)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        aspect.ratio = 4/11,
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

## Create a composition with the images and the plots ---- 
### Add harbor restoration images ----
asv_tab_all_bloo_z_tax_sig$harbor_effect <- factor(asv_tab_all_bloo_z_tax_sig$harbor_effect, levels = c('negative', 'positive'))

asv_tab_all_bloo_z_tax_sig$type_of_bloomer <- factor(asv_tab_all_bloo_z_tax_sig$type_of_bloomer, levels = c('Seasonal', 'Chaotic'))

harbor_pa_stats_plot <- asv_tab_all_bloo_z_tax_sig |>
  dplyr::filter(significant == 'TRUE') |>
  dplyr::mutate(harbor_effect = case_when(harbor_effect == 'positive' ~ 'Positive',
                                          harbor_effect == 'negative' ~ 'Negative')) |>
  dplyr::group_by(date, fraction, order_f, family_f, abundance_type, type_of_bloomer, harbor_effect) |>
  dplyr::reframe(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abund_max))+
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
  geom_area(aes(group = family_f, fill = family_f), position = 'stack')+
  #facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
  #facet_grid(type_of_bloomer~harbor_effect,  scales = 'free_y')+
  facet_wrap(vars(interaction(type_of_bloomer, harbor_effect)),  scales = 'free_y', ncol = 1)+
  guides(fill = guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        #aspect.ratio = 6/11,
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
# ggsave(filename = 'harbor_restoration_plot_photo_rclr_stats_v4.pdf',
#        plot = composition_with_labels,
#        path = 'Results/Figures/',
#        width = 180, height = 200, units = 'mm')


## harbor restoration reads ----
asv_tab_10y_l_rel |>
  colnames()

asv_tab_10y_l_rel |>
  dplyr::group_by(sample_id) |>
  reframe(total_reads = sum(reads)) |>
  arrange(total_reads)
  
