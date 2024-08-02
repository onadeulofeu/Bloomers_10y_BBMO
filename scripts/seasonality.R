# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# packages
#library(fdrtool)
#library(GeneCycle)
library(patchwork)
library(lubridate)
library(ggthemes)
library(stringr)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(waveslim)
library(vegan) #deconstand

# Differentiate different types of bloomers studying their signals across time----
## upload data----
asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') |> ## in this dataset I have all the information from potential blooming ASVs
  as_tibble() |>
  dplyr::select(-X)
  
bloo_02 <- read.csv('data/detect_bloo/bloo_02.csv') |>
  as_tibble()

bloo_3 <-read.csv('data/detect_bloo/bloo_3.csv') |>
  as_tibble()

## Previous to the analysis we need to interpolate missing samples----
### we use the most straightforward way which is linear interpolation
#### for the FL faction I have already a complete dataset
asv_tab_all_bloo_z_tax_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2')

asv_tab_all_bloo_z_tax_02_sam <- asv_tab_all_bloo_z_tax_02 |>
  dplyr::select(date) |>
  distinct(date) |>
  as_tibble()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  group_by(sample_id) |>
  dplyr::reframe(n_samples = n())

#### for the particle attached fraction
asv_tab_all_bloo_z_tax_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3')

asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3') |>
  group_by(sample_id) |>
  dplyr::reframe(n_samples = n()) #117 samples

asv_tab_all_bloo_z_tax_3_sam <- asv_tab_all_bloo_z_tax_3 |>
  dplyr::select(date) |>
  distinct(date) |>
  as_tibble()
 
asv_tab_all_bloo_z_tax_02_sam |>
  anti_join(asv_tab_all_bloo_z_tax_3_sam)

#### Missing samples to be interpolated in the particle attached fraction-----
### 2004-03-22
### 2005-02-15
### 2005-05-10

## prepare the dataset previous to computing the interpolation----
asv_tab_all_bloo_3_interpol <- asv_tab_10y_3_rel |>
  dplyr::left_join(m_3) |>
  dplyr::select(date, asv_num, relative_abundance) |>
  # dplyr::filter(#asv_num == 'asv179' &
  #                 abundance_type == 'relative_abundance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

# Convert the target date to a Date object
target_date1 <- as.POSIXct("2004-03-22", format = "%Y-%m-%d")
target_date2 <- as.POSIXct("2005-02-15", format = "%Y-%m-%d")
target_date3 <- as.POSIXct("2005-05-10", format = "%Y-%m-%d")

# Perform linear interpolation
asv_tab_all_bloo_3_interpol <- asv_tab_all_bloo_3_interpol |>
  group_by(asv_num) |>
  dplyr::reframe(target_date1 = approx(date, relative_abundance, xout = target_date1)$y, #The approx function fits a linear spline to the data and estimates values at the specified target dates.
                 target_date2 = approx(date, relative_abundance, xout = target_date2)$y,
                 target_date3 = approx(date, relative_abundance, xout = target_date3)$y) |>
  ungroup() |>
  pivot_longer(cols = starts_with('target'), values_to = 'abundance_value', names_to = 'date') |> #prepare the dataset to add it to the general dataset
  dplyr::mutate(date = case_when(date == 'target_date1' ~ '2004-03-22',
                                 date == 'target_date2' ~ '2005-02-15',
                                 date == 'target_date3' ~ '2005-05-10')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

asv_tab_all_bloo_3_interpol_filt <- asv_tab_all_bloo_3_interpol |>
  dplyr::filter(asv_num %in% bloo_3$value)

## check that the interpolation is functioning well and mantaining the sum of total relative abundances at 1
 asv_tab_all_bloo_3_interpol |>
  group_by(date) |>
  dplyr::reframe(sum = sum(abundance_value))

##check how do the interpolations look like are they normal?----
 asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(date, abundance_value, asv_num, abundance_type) |>
  dplyr::filter(#asv_num == 'asv179' &
    abundance_type == 'relative_abundance') |>
  dplyr::select(-abundance_type) |>
  bind_rows(asv_tab_all_bloo_3_interpol_filt) |>
  dplyr::mutate(interpolation = if_else(date %in% c(as.POSIXct('2004-03-22', format = "%Y-%m-%d"),
                                                   as.POSIXct('2005-02-15', format = "%Y-%m-%d"),
                                                   as.POSIXct('2005-05-10', format = "%Y-%m-%d")),  '#9F0011', '#080808', missing = '#080808')) |>
  ggplot(aes(date, abundance_value))+
  facet_wrap(vars(asv_num))+
  scale_color_identity()+
  geom_vline(xintercept = c(as.POSIXct('2004-03-22', format = "%Y-%m-%d"),
                            as.POSIXct('2005-02-15', format = "%Y-%m-%d"),
                            as.POSIXct('2005-05-10', format = "%Y-%m-%d")))+
  geom_point(aes(color = interpolation), size = 1, alpha = 0.8)+
  geom_line()+
  theme_bw()+
  theme(text = element_text(size = 3))

##they all look fine so I will add them to the main dataset to proceed with the seasonality analysis

## add the interpolated samples to the dataset----

## I add the metadata and taxonomic data to them, and then I add it to the PA general dataset
asv_tab_all_bloo_z_tax_02_subset <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(date %in% c(as.POSIXct('2004-03-22', format = "%Y-%m-%d"),
                            as.POSIXct('2005-02-15', format = "%Y-%m-%d"),
                            as.POSIXct('2005-05-10', format = "%Y-%m-%d"))) |>
  dplyr::filter(#asv_num == 'asv179' &
    abundance_type == 'relative_abundance') |>
  dplyr::select(-seq, -asv_num, -family, - phylum, -class, -order, -genus, -abundance_type, -abundance_value, -z_score_ra, -sample_id_num, -domain, -sample_id, -reads) |>
  dplyr::distinct()

asv_tab_all_bloo_z_tax_3 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 

asv_tab_all_bloo_3_complete <- asv_tab_all_bloo_3_interpol_filt |>
  right_join(asv_tab_all_bloo_z_tax_02_subset, by = ('date')) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  bind_rows(asv_tab_all_bloo_z_tax_3)

asv_tab_all_bloo_3_complete %$%
  date |>
  unique() #120 we have all the data!!

asv_tab_all_bloo_z_tax_02 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2')

### Dataset to perform the seasonality analysis on----
asv_tab_all_bloo_3_complete #PA datset
asv_tab_all_bloo_z_tax_02 #FL dataset

####### MAKING A DECISION ON USING WAVELETS // FOURIER ANALYSIS // ARIMA FOR OUR DATASET --------

## The seasonal dynamics of persistent and intermittent individual bacterial taxa was analysed by Fourier time-series analysis. 
## Periodic components were extracted, and their significance was determined by Anderson and the Fisher G-test.

## Fourier transform works when no variation in time happens, when we are in the frequency-domain we lose this dependency.
 
# fisher.g.test calculates the p-value(s) according to Fisher's exact g test for one or more time series. 
## This test is useful to detect hidden periodicities of unknown frequency in a data set. For an application to microarray data see Wichert, Fokianos, and Strimmer (2004).

## Wavelet analysis enables investigation of time series characterized by different periodicities and is particularly suited for time series which are 
## not stationary, as applies to many biological systems. (nonstationarity is the status of a time series whose statistical properties are changing through time)

## We saw in the last section that we need a method to handle signals whose constituent frequencies vary over time (e.g. the ECG data). 
## We need a tool that has high resolution in the frequency domain and also in the time domain, that allows us to know at which frequencies the signal oscillates, and at which time these oscillations occur.
## The Wavelet transform fullfils these two conditions.
# 
# By means of wavelet transform a time series can be decomposed into
# a time dependent sum of frequency components. As a result we are able
# to capture seasonalities with time-varying period and intensity, which
# nourishes the belief that incorporating the wavelet transform in existing forecasting methods can improve their quality.

# A widely used approach is the Autoregressive Moving Average (ARIMA) model,
# which captures intertemporal linear dependence in the data itself as well as in the
# error term.
# 
# This is why Wong et al. (2003) use the wavelet transform, which is able to
# capture dynamics in period and intensity, to model both the trend and the seasonality. 
# By means of the wavelet transform we can decompose a time series into a
# linear combination of different frequencies

## Why do we choosed wavelet transform: traditional Fourier transform (decomposition using sine 
## and cosine) gives global average over entire signal, thus may obscure local information 
## (which we are very interested in)

## Continuous wavelet transform (CWT) we pick basicaly every possible scale in location
## Discrete wavelet transform (DWT) there's a discrete number of wavelets scales and locations

## Wavelet analysis is also well suited when signals exhibit sharp spikes or local
## discontinuities, features that are very poorly represented by Fourier analysis
## file:///Users/onadeulofeucapo/Downloads/engeland.pdf

## prior to the analysis we need to use CLR transformed data------
wavelet_02_df <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::select(abundance_value, asv_num, decimal_date)

# admissibility average arround 0-----
wavelet_02_df |>
  dplyr::group_by(asv_num) |>
  dplyr::filter(!asv_num %in%  c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::reframe(mean = mean(abundance_value))

# ## I run again the wavelets analysis with rCLR values from deconstant function not zcompositions, to use the same transformation for both PA and FL----
# asv_tab_bbmo_10y_w_02 <- asv_tab_bbmo_10y_l |>
#   dplyr::filter(str_detect(sample_id, '_0.2_')) |>
#   left_join(m_02) |>
#   dplyr::select(asv_num, reads, date) |>
#   pivot_wider(names_from = 'asv_num', values_from = 'reads', values_fill = 0) |>
#   as.data.frame()
# 
# rownames(asv_tab_bbmo_10y_w_02) <- asv_tab_bbmo_10y_w_02$date
# 
# # asv_tab_bbmo_10y_w_02_inter |>
# #   dim()
# 
# asv_tab_bbmo_10y_w_02 <- asv_tab_bbmo_10y_w_02[,-1]
# 
# #geometric mean
# gm <- function(x){
#   exp(mean(log(x[x>0])))
# }
# 
# ### i use the deconstant function 
# zclr_df_02 <- decostand(asv_tab_bbmo_10y_w_02, method = 'rclr') |>
#   as_tibble(rownames = "date") %>%
#   pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'zclr') 
# 
# ## filter it for my bloomers (the one's I want to perform the wavelets analysis on)
# wavelet_02_df_deconstand <-  zclr_df_02 |>
#   dplyr::filter(asv_num %in% bloo_02$value) |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   left_join(m_02, by = c('date')) |>
#   dplyr::select(decimal_date, asv_num, zclr)
# 
# asv_tab_10y_02_zclr_inter_bloo %$%
#   unique(decimal_date) ## one sample is missing needs to be solved.
# 
# ## observe them
# zclr_df_02 |>
#   dplyr::filter(asv_num %in% bloo_02$value) |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   ggplot(aes(date, zclr, group = asv_num))+
#   geom_line()+
#   geom_point()+
#   facet_wrap(vars(asv_num))+
#   theme_bw()
# 
# zclr_df_02  %$%
#   zclr |>
#   range()
# 
# ### comparison of the results using the zCLR transformation and the deconstand function (observe if there are important differences) ----
# wavelet_02_df |>
#   dplyr::filter(asv_num %in% bloo_02$value) |>
#   #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   ggplot(aes(decimal_date, abundance_value, group = asv_num))+
#   geom_line()+
#   geom_point()+
#   facet_wrap(vars(asv_num))+
#   theme_bw()

##compare the values obtained using one approximation or the other (we keep the analysis with the pseudcount approximation).-----
# x <- wavelet_02_df |>
#   dplyr::filter(asv_num == 'asv237') |>
#   rename(variable = abundance_value) |>
#   dplyr::mutate(type = 'zclr_zcomp')
# 
# y <- wavelet_02_df_deconstand |>
#   dplyr::filter(asv_num == 'asv237') |>
#   rename(variable = zclr) |>
#   dplyr::mutate(type = 'zclr_deconst')
# 
# reads <- asv_tab_bbmo_10y_l |>
#   dplyr::filter(str_detect(sample_id, '_0.2_')) |>
#   left_join(m_02) |>
#   dplyr::select(asv_num, reads, decimal_date) |>
#   dplyr::filter(asv_num == 'asv237') |>
#   rename(variable = reads)|>
#   dplyr::mutate(type = 'reads')
# 
# asv237_transfor <- bind_rows(reads, x, y) |>
#   pivot_wider(id_cols = c(asv_num, decimal_date), names_from = type, values_from = variable)

#write.csv(asv237_transfor, 'data/asv237_transfor.csv')

## I use deconstand in both fractions in the following analysis----
# wavelet_02_df <- wavelet_02_df_deconstand

##I save both datasets
#write.csv2(wavelet_02_df, 'data/wavelet_02_df_deconstand.csv') #using the deconstand function (the same we need to use for the PA wavelets)
#write.csv2(wavelet_02_df, 'data/wavelet_02_df_zclr.csv') #using the zcompositions function which deals with 0 in a more elegant way.

## for the PA fraction we need to calculate again the rCLR transformation after the interpolation of the missing samples ----
### i feel that i makes more sense to make the interpolation from relative_abundances therefore I need to convert it's relative abundance interpolated
### value to interpolated reads, for this I imagine that they have the total reads that we have in the mean dataset

mean_reads <- asv_tab_10y_3_rel |> #### total reads for the whole dataset not only for bloomers
  group_by(sample_id) |>
  dplyr::reframe(sum = sum(reads)) |>
  dplyr::mutate(mean_reads = mean(sum)) ##39383 reads / sample

##interpolated relative abundances from missing samples an transform them to reads
asv_tab_all_3_reads <-  asv_tab_all_bloo_3_interpol |>
  dplyr::select(asv_num, date, abundance_value) |> ## the abundances interpolated are in relative abundance
   dplyr::filter(date %in% c(as.POSIXct('2004-03-22', format = "%Y-%m-%d"),
                             as.POSIXct('2005-02-15', format = "%Y-%m-%d"),
                             as.POSIXct('2005-05-10', format = "%Y-%m-%d"))) |>
   #dplyr::group_by(date) |>
   #dplyr::reframe(sum = sum(abundance_values))
   dplyr::mutate(reads = abundance_value*39383) |>
  dplyr::select(asv_num, reads, date)

### add to the main dataset and calculate the rCLR transformation using the interpolated reads for the missing samples
asv_tab_bbmo_10y_w_3_inter <- asv_tab_bbmo_10y_l |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  left_join(m_3) |>
  dplyr::select(asv_num, reads, date) |>
  bind_rows(asv_tab_all_3_reads) |>
  pivot_wider(names_from = 'asv_num', values_from = 'reads', values_fill = 0) |>
  as.data.frame()

rownames(asv_tab_bbmo_10y_w_3_inter) <- asv_tab_bbmo_10y_w_3_inter$date

# asv_tab_bbmo_10y_w_3_inter |>
#   dim()

asv_tab_bbmo_10y_w_3_inter <- asv_tab_bbmo_10y_w_3_inter[,-1]

#geometric mean
gm <- function(x){
  exp(mean(log(x[x>0])))
}

## with this transformation I'm losing samples (due to too much 0 in some samples, z.warning set up to 0.99 to keep all samples)
### at 0.8 (default) I lose 30 samples which belonged to the years corresponding to harbour remodelation 
# ### I don't lose samples but I lose ASVs.
# zclr_df_inter <- cmultRepl(asv_tab_bbmo_10y_w_3_inter, method = 'CZM', output = 'p-count', z.warning = 0.99
#                      #adjust = 0.2,    t = 237, s = 7849
# ) |>
#   as_tibble(rownames = "sample_id") %>%
#   pivot_longer(-sample_id, names_to = 'asv_num') %>%
#   group_by(sample_id) %>%
#   dplyr::mutate(zclr = log(value/gm(value))) %>%
#   ungroup() %>%
#   dplyr::select(-value) %>%
#   pivot_wider(names_from = name, values_from = zclr, values_fill = 0) %>%
#   column_to_rownames("sample_id")
# 
# zclr_df_inter |>
#   dim()
# 
# zclr_df_inter_deconstand |>
#   dim()

### i use the deconstant function in this case because the cmultRepl was making me lose one sample
zclr_df_inter_deconstand <- decostand(asv_tab_bbmo_10y_w_3_inter, method = 'rclr') |>
  as_tibble(rownames = "date") %>%
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'rclr') 

# clr_df_inter_deconstand <- decostand(asv_tab_bbmo_10y_w_3_inter, method = 'clr', pseudocount = 1) |>
#   as_tibble(rownames = "date") %>%
#   pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'clr') 

## filter it for my bloomers (the one's I want to perform the wavelets analysis on)
wavelet_3_df <-  zclr_df_inter_deconstand |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(m_02, by = c('date')) |>
  dplyr::select(decimal_date, asv_num, rclr)

### admissibility average around 0-----
wavelet_3_df |>
  dplyr::group_by(asv_num) |>
  dplyr::mutate(rclr = as.numeric(rclr)) |>
  dplyr::reframe(mean = mean(rclr)) |>
  arrange(-mean)
# 
# asv_tab_10y_3_zclr_inter_bloo %$%
#   unique(decimal_date) ## one sample is missing needs to be solved.

## observe them
 zclr_df_inter |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, zclr, group = asv_num))+
  geom_line()+
  geom_point()+
  facet_wrap(vars(asv_num))+
  theme_bw()

zclr_df_inter  %$%
  zclr |>
  range()

#write.csv(wavelet_3_df, 'data/wavelet_3_df_deconstand.csv') #using the deconstand function (the same we need to use for the PA wavelets)

#### Notice that deconstant robust CLR transformation gives me more negative values than the cmultRepl step

## Check sampling frequency----
## In timeseries analysis whether if we perform Fourier analysis or Wavelets they assume constant sampling. We check it before we go any further.
###check normality before correlation
shapiro.test(as.numeric(sampling_freq$real_sampling)) # => 0.0004153 (NO-NORMALITY)
ggqqplot(as.numeric(sampling_freq$real_sampling))

shapiro.test(as.numeric(perfect_sampling_freq$ideal_sampling)) #p-value = 0.0004823 (NO-NORMALITY)
ggqqplot(as.numeric(perfect_sampling_freq$ideal_sampling))

check_sampling_constant <-  sampling_freq |>
  bind_cols(perfect_sampling_freq) |>
  ggplot(aes(real_sampling, ideal_sampling))+
  geom_smooth(method = 'lm', color = 'black', alpha = 0.3)+
  geom_point(alpha = 0.6)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  stat_cor(aes(label = paste(..p.label..)), label.x = 2006,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman')+
  labs(x = 'Real sampling', y = 'Ideal sampling')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 5),
        axis.ticks.length = unit(.25, "mm"))

# ggsave(check_sampling_constant,  filename = 'check_sampling_frequency.pdf',
#        path = 'Results/Figures/',
#        width = 88, height = 88, units = 'mm')

## Adria's functions
## https://gitlab.com/aauladell/AAP_time_series here we have the source.

# Takes a data.frame, reformulates the date and arranges rows in function of 
# this value 
reorder_chrono <- function(df) {
  dplyr::mutate(df, year = str_c("20", year),
         month = str_to_title(month)) %>%
    unite(dat, year, month, day, sep = "-", remove = F) %>%
    mutate(dat = lubridate::ymd(dat),
           decimaldat = decimal_date(dat)) %>%
    arrange(dat)
}

check_seasonality_fisher <- function(ls){
  ls %>% 
    map(~ts_conversor(.x)) %>% 
    map(~data.frame(per = periodogram(.)$freq[
      which.max(periodogram(.)$spec)],
      whole = periodogram(.),
      pval = fisher.g.test(.))) %>% 
    bind_rows( .id = "asv") 
}

ts_conversor <- function(df, s=c(2004, 1), f=12){
  ts(df$abundance_value, start = s,frequency = f) 
}

check_periodogram <- function(df,wrap="phylogroup"){
  
  left_join(df, BL.psmelt.raw, by = c("asv" = "OTU")) %>%
    distinct(asv,whole.freq, whole.spec, .keep_all = TRUE) %>% 
    ggplot() +
    geom_line(aes(whole.freq, whole.spec)) +
    facet_wrap(wrap, scales = "free_y") +
    lil.strip +
    ggtitle(paste0("Periodogram distribution"))
  
}

## Fourier time series analysis

## We need values transformed to CLR

four <- asv_tab_all_bloo_z_tax_02 |>
  dplyr::filter(abundance_type == 'zclr') |>
  dplyr::select(asv_num,abundance_value, month, year, day, family_f, decimal_date) |>
  arrange(decimal_date)

## Seasonality of each ASV that we believe can be a potential bloomer

## I divide the different abundance values in separated plots and I filter only the potential bloomers in FL or PA-----

bloo_02$value

# [1] "asv38"  "asv8"   "asv5"   "asv3"   "asv2"   "asv15"  "asv18"  "asv27"  "asv17"  "asv555"
# [11] "asv237" "asv563" "asv62"  "asv58"  "asv1"   "asv7"   "asv178" "asv282" "asv11" 

bloo_3$value

# [1] "asv179" "asv225" "asv264" "asv200" "asv15"  "asv385" "asv72"  "asv471" "asv27"  "asv153"
# [11] "asv17"  "asv77"  "asv43"  "asv192" "asv84"  "asv118" "asv311" "asv23"  "asv85"  "asv25" 
# [21] "asv163" "asv80"  "asv69"  "asv219" "asv116" "asv182" "asv126" "asv223" "asv105" "asv28" 
# [31] "asv1"   "asv7"   "asv4"   "asv31"  "asv276" "asv22"  "asv194" "asv559" "asv11"  "asv100"
# [41] "asv302" "asv113" "asv511" "asv42"  "asv49" 

#ADRIÀ CODE: FOURIER ANALYSIS Establish the numbers as a time series----
four |>
  colnames()

asv1.ts <- four %>%
  filter(asv_num == "asv8") %>%
  ts_conversor()

## I try to apply this part of the code to all potential bloomers ASVs-------
# Applying the code for different asv_num values
results_list <- map(bloo_02$value, ~four %>%
                      filter(asv_num == .x) %>%
                      ts_conversor())

# Combining the results into a single tibble
results_list |>
  class()

# Convert each time series to a data frame
results_list_df <- lapply(results_list, function(ts) as.data.frame(ts))

# Bind the rows
results_tibble <- bind_rows(results_list_df)

## results_tibble <- bind_rows(results_list)

# Is it seasonal?
## The fisher.g.test calculates the p-value(s) according to Fisher's exact g test for one or more time series. This test is
## useful to detect hidden periodicities of unknown frequency in a dataset. 
check_seasonality_fisher <- function(ls){
  ls %>% 
    map(~ts_conversor(.x)) %>% 
    map(~data.frame(per = periodogram(.)$freq[
      which.max(periodogram(.)$spec)],
      whole = periodogram(.),
      pval = fisher.g.test(.))) %>% 
    bind_rows( .id = "asv") 
}

fisher.g.test(asv1.ts)

descomposition <- stl(asv1.ts, s.window = "periodic")
plot(descomposition, main = "ASV1")
summary(descomposition)

# Does the periodogram presents a clear seasonality
ggplot(as.data.frame(periodogram(asv1.ts)), aes(freq, spec)) +
  geom_line() +
  geom_vline(xintercept = 0.083, color = "red")

# How does the ACF behaves?
acf(asv1.ts, lag.max = 120)


## Fisher G-test to determine the significance of periodic components
library(GeneCycle)

## The time series was decomposed in three components, seasonal periodicity 

## (oscillation inside each period), trend (evolution over time) and residuals, 
## through local regression by the stl function.


## Autocorrelogram was calculated using the acf function


### WAVELETS ANALYSIS-----
## need to pick a family of wavelets we want to work with.
# Here I try to apply the code that is in appendix A of the thesis -----

## MODWT (Maximal overlap discrete wavelet transforms):

### Create useful code for just one ASV----

# We need the abundance vector.
### example for asv11
abund11 <- asv11$abundance_value
modwtasv11 <- modwt(abund11, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
  phase.shift(wf = 'la8') # this correction is necessary because the filtering with wavelet filter induces a displacement of features relative to the original time series.
  
modwtasv11_biased <- modwt(abund11, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  #brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
  phase.shift(wf = 'la8') # this correction is necessary because the filtering with wavelet filter induces a displacement of features relative to the original time series.

### example for asv178
abund178 <- asv178$abundance_value
modwtasv178 <- modwt(abund178, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
  phase.shift(wf = 'la8') 

modwtasv178 |>
  str()

modwtasv178_biased <- modwt(abund178, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  #brick.wall(wf = 'la8') |>
  phase.shift(wf = 'la8') # this correction is necessary because the filtering with wavelet filter induces a displacement of features relative to the original time series.

### example for asv194----
abund194 <- asv194$abundance_value
modwtasv194 <- modwt(abund194, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
  phase.shift(wf = 'la8') 

modwtasv194_biased <- modwt(abund194, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  #brick.wall(wf = 'la8') |>
  phase.shift(wf = 'la8') # this correction is necessary because the filtering with wavelet filter induces a displacement of features relative to the original time series.


## I try to apply this part of the code for the whole dataset by grouping----

### We would like to perform the modwt on all ASVs and have them in a tibble with many lists inside.

# modwt.function <-function(abundance){
#  modwt_result <-  abundance |>
#    modwt( wf = 'la8', boundary = "periodic", n.levels = 4) |>
#     brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
#     phase.shift(wf = 'la8')
#   return(modwt_result)
#  }  
# 
# ## FL fraction
# # wavelet_02_df |>
# #   dplyr::group_by(asv_num) |>
# #   #dplyr::filter(asv_num == 'asv11') |>
# #   dplyr::select(abundance_value) %$%
# #   modwt.funtion(abundance = abundance_value)
# 
# modwt_results_02 <- wavelet_02_df |>
#   group_by(asv_num) %>%
#   dplyr::summarize(modwt_result = list(modwt.function(abundance_value)))
# 
# modwt_results_02 |>
#   str()
# 
# modwt_results_3 |>
#   str()
# 
# # # Continue with phase shift
# # shifted_result <- phase.shift(brick.wall(result, filter = db5_filter), filter = db5_filter)
# 
# # e-folding estimates (modified from Appendix A)
# x <- rep(0, 10001) 
# x[5000] <- 1 
# n.levels <- 4 
# len <- length(abund11) 
# temp <- phase.shift(modwt(x, n.levels = n.levels, wf = "la8"), wf = "la8")
# 
#  waveExtremes <- matrix(nrow = 3, ncol = n.levels + 1) 
#  colnames(waveExtremes) <- c(paste("d", 1:n.levels, sep = ""), paste("s", n.levels,  sep = "")) 
#  rownames(waveExtremes) <- c("left", "right", "top")
# 
#   for (i in 1:(n.levels + 1)) waveExtremes[, i] <- c(range(which(abs(temp[[i]]) 
#                                                                  >= max(abs(temp[[i]]))/(exp(2)))), which.max(abs(temp[[i]])))
# 
#  boundaries <- data.frame(end = len - (waveExtremes[3, ] - waveExtremes[1, ]), 
#                           start = waveExtremes[2, ] - waveExtremes[3, ])
#  
#  for (j in 1:(n.levels + 1)) { 
#    is.na(modwtasv11[[i]]) <- c(1:boundaries$start[i], boundaries$end[i]:length(modwtasv11[[i]])) 
#    is.na(modwtasv178[[i]]) <- c(1:boundaries$start[i], boundaries$end[i]:length(modwtasv178[[i]])) 
#  }
#  
#  ## i modify the function to apply it to the results I created before, which is the tibble with many lists inside-----
#  ### Elimination of the e-folding coefficients is accomplished using the following script where a constant zero valued series with a perturbation of 
#  ### value 1 in the middle (=x) is transformed.
#  x <- rep(0, 10001) 
#  x[5000] <- 1 
#  n.levels <- 4 
#  len <- length(abund11) 
#  temp <- phase.shift(modwt(x, n.levels = n.levels, wf = "la8"), wf = "la8")
#  
#  ## The positions to the left and to the right of the maximal influence of this spike are recorded in a matrix (left, right) together with the 
#  ## position of the maximum itself (top).
#  waveExtremes <- matrix(nrow = 3, ncol = n.levels + 1) 
#  colnames(waveExtremes) <- c(paste("d", 1:n.levels, sep = ""), paste("s", n.levels,  sep = "")) 
#  rownames(waveExtremes) <- c("left", "right", "top")
#  
#  ## The distance to the maximum from both sides of the influence is determined as 1/e2 times the maximum within a specific coefficient vector.
#  for (i in 1:(n.levels + 1)) waveExtremes[, i] <- c(range(which(abs(temp[[i]]) 
#                                                                 >= max(abs(temp[[i]]))/(exp(2)))), which.max(abs(temp[[i]])))
#  
#  ## The positions (waveExtremes) are used to calculate the distances to the left and to the right of the influence maximum. 
#  ## The distance to the left of the maximum is called "right" because it will serve to calculate the distance at the end of the series.
#  
#  boundaries <- data.frame(end = len - (waveExtremes[3, ] - waveExtremes[1, ]), 
#                           start = waveExtremes[2, ] - waveExtremes[3, ])
#  
#  ## The actual elimination of the coefficients within the margin defined by the e-folding boundaries: (this is what I want to apply to each ASV)----
 # modwt_results_02$asv_num
 # 
 # modwt_results_02$modwt_result
 # 
 # for (j in 1:(n.levels + 1)) { 
 #   is.na(modwt_results_02$ modwtasv11[[i]]) <- c(1:boundaries$start[i], boundaries$end[i]:length(modwtasv11[[i]])) 
 #   
 # }
 # 
 # ##another try
 # 
 # for (i in 1:length(modwt_results_02$asv_num)) {
 #   asv_num_index <- i
 #   modwt_results <- modwt_results_02$modwt_result[[asv_num_index]]
 #   
 #   for (j in 1:(n.levels + 1)) {
 #     # Check if boundaries$start[i] is not NA before proceeding
 #     if (!is.na(boundaries$start[i])) {
 #       is.na(modwt_results[[j]]) <- c(1:boundaries$start[i], boundaries$end[i]:length(modwt_results[[j]]))
 #     } else {
 #       # Handle the case where boundaries$start[i] is NA (modify as needed)
 #       warning(paste("Skipping index", i, "due to NA in boundaries$start."))
 #     }
 #   }
 # }
 # 
 # 
 # 
 ##PLOTS----
 
 #### Interpretation of the following plots----
 
 #### MODWTs of CLR transformation of the relative abundance from the potential blooming ASVs. The series d1-d4 are wavelet coefficient
 ### vectors; the series s4 is the scaling coefficient vector at the coarsest scale. 
 
 ### d1 is the finest scale and d4 is the coarsest. The s4 depicts the scaling coefficients as a function of time, representing the remainder
 ### of the variability after the finer scales (d1-d4) have been isolated. 
 
 ### The inner red, solid, vertical lines represent the limits of the boundary coefficients. The regions enclosed by either set of limits become
 ### progressively smaller as wavelet filters become wider. However, by allowing some bias (i.e as given by the region within the blue lines in each
 ### graph) event the coarser scales can be investigated over a relatively larger temporal scale.
 
 ### To observe which component dominates we need to observe the one that has the highest coefficients in absolute value.
 
### MODWT of the ASV11-----
 # Create an empty list to store your plots
 par(mfrow = c(n.levels + 1, 1), mar = c(2, 3, 0.5, 0.2))
 
 # Assuming n.levels is defined before this code snippet
 
 # Create an empty list to store your plots
 plotList <- vector("list", n.levels)
 
 for (i in 1:n.levels) { 
   plot(as.vector(time(abund178)), modwtasv11_biased[[i]], xlab = "", ylab = "", type = "h", 
        bty = "n", axes = F, ylim = range(modwtasv11[1:n.levels], na.rm = TRUE)) 
   axis(2, cex.axis = 1.5, las = 2) +
     abline(v = seq(start(abund178)[1], end(abund178)[1], by = 1), col = "gray", 
            lty = "dashed") + abline(h = 0) 
   abline(v = time(abund178)[range(which(!is.na(modwtasv11[[i]])))], lwd = 2, 
          col = "red") + 
     abline(v = time(abund178)[c(boundaries$start[i], boundaries$end[i])], col = "blue", 
            lwd = 3, lty = "dashed") +
     mtext(side = 4, line = -2, cex = 1.2, outer = FALSE, at = par("usr")[4], text = paste("d", i, sep = ""), las = 2) 
   
   # Save each plot to the list
   plotList[[i]] <- recordPlot()
 }
 
 # Set the title for the entire series of plots
 title(main = "ASV 178")
 
 # Replay and display the saved plots
 for (i in 1:n.levels) {
   replayPlot(plotList[[i]])
 }
 
 plot(as.vector(time(abund178)), modwtasv11_biased$s4, xlab = "", ylab = "", type = "h",
      bty = "n", axes = F) > axis(1, cex.axis = 1.5) 
 axis(2, cex.axis = 1.5, las = 2) 
 mtext(side = 4, line = -2, cex = 1.2, outer = F, at = par("usr")[4], text = paste("s4", 
                                                                                   sep = ""), las = 2) +
   abline(v = time(abund178)[range(which(!is.na(modwtasv11[[5]])))], lwd = 2, 
          col = "red") 
 title(main = "ASV 178")
 abline(v = time(abund178)[c(boundaries$start[5], boundaries$end[5])], col = "blue", 
        lwd = 3, lty = "dashed")

##MODWT of the ASV178
par(mfrow = c(n.levels + 1, 1), mar = c(2, 3, 0.5, 0.2))
# Assuming n.levels is defined before this code snippet

# Create an empty list to store your plots
plotList <- vector("list", n.levels)

for (i in 1:n.levels) { 
  plot(as.vector(time(abund178)), modwtasv178_biased[[i]], xlab = "", ylab = "", type = "h", 
       bty = "n", axes = F, ylim = range(modwtasv178[1:n.levels], na.rm = TRUE)) 
  axis(2, cex.axis = 1.5, las = 2) +
    abline(v = seq(start(abund178)[1], end(abund178)[1], by = 12), col = "gray", 
           lty = "dashed") + abline(h = 0) 
  abline(v = time(abund178)[range(which(!is.na(modwtasv178[[i]])))], lwd = 2, 
         col = "red") + 
    abline(v = time(abund178)[c(boundaries$start[i], boundaries$end[i])], col = "blue", 
           lwd = 3, lty = "dashed") +
    mtext(side = 4, line = -2, cex = 1.2, outer = FALSE, at = par("usr")[4], text = paste("d", i, sep = ""), las = 2) 
  
  # Save each plot to the list
  plotList[[i]] <- recordPlot()
}

# Set the title for the entire series of plots
title(main = "ASV 178")

# Replay and display the saved plots
for (i in 1:n.levels) {
  replayPlot(plotList[[i]])
}

plot(as.vector(time(abund178)), modwtasv178_biased$s4, xlab = "", ylab = "", type = "h",
     bty = "n", axes = F) > axis(1, cex.axis = 1.5) 
axis(2, cex.axis = 1.5, las = 2) 
mtext(side = 4, line = -2, cex = 1.2, outer = F, at = par("usr")[4], text = paste("s4", 
                                                                                  sep = ""), las = 2) +
abline(v = time(abund178)[range(which(!is.na(modwtasv178[[5]])))], lwd = 2, 
       col = "red") 

title(main = "ASV 178")
abline(v = time(abund178)[c(boundaries$start[5], boundaries$end[5])], col = "blue", 
       lwd = 3, lty = "dashed")

## MODWT for 194
par(mfrow = c(n.levels + 1, 1), mar = c(2, 3, 0.5, 0.2))
# Assuming n.levels is defined before this code snippet

# Create an empty list to store your plots
plotList <- vector("list", n.levels)

for (i in 1:n.levels) { 
  plot(as.vector(time(abund194)), modwtasv194_biased[[i]], xlab = "", ylab = "", type = "h", 
       bty = "n", axes = F, ylim = range(modwtasv194[1:n.levels], na.rm = TRUE)) 
  axis(2, cex.axis = 1.5, las = 2) +
    abline(v = seq(start(abund194)[1], end(abund194)[1], by = 12), col = "gray", 
           lty = "dashed") + abline(h = 0) 
  abline(v = time(abund194)[range(which(!is.na(modwtasv194[[i]])))], lwd = 2, 
         col = "red") + 
    abline(v = time(abund194)[c(boundaries$start[i], boundaries$end[i])], col = "blue", 
           lwd = 3, lty = "dashed") +
    mtext(side = 4, line = -2, cex = 1.2, outer = FALSE, at = par("usr")[4], text = paste("d", i, sep = ""), las = 2) 
  
  # Save each plot to the list
  plotList[[i]] <- recordPlot()
}

# Set the title for the entire series of plots
title(main = "ASV 194")

# Replay and display the saved plots
for (i in 1:n.levels) {
  replayPlot(plotList[[i]])
}

plot(as.vector(time(abund194)), modwtasv194_biased$s4, xlab = "", ylab = "", type = "h",
     bty = "n", axes = F) > axis(1, cex.axis = 1.5) 
axis(2, cex.axis = 1.5, las = 2) 
mtext(side = 4, line = -2, cex = 1.2, outer = F, at = par("usr")[4], text = paste("s4", 
                                                                                  sep = ""), las = 2) +
  abline(v = time(abund194)[range(which(!is.na(modwtasv194[[5]])))], lwd = 2, 
         col = "red") 

title(main = "ASV 194")
abline(v = time(abund194)[c(boundaries$start[5], boundaries$end[5])], col = "blue", 
       lwd = 3, lty = "dashed")

## Magnitude of the coefficients----- 
# Assuming you have computed the MODWT
modwt_result <- modwt(abund178, wf = 'la8', boundary = "periodic", n.levels = 4)

# Calculate the magnitude of coefficients for each level
coeff_magnitudes <- sapply(modwt_result, function(level) sqrt(sum(level^2)))

# Find the level with the maximum magnitude
selected_level <- which.max(coeff_magnitudes)

# Use the coefficients from the selected level for further analysis
selected_coefficients <- modwt_result[[selected_level]]

# Assuming you have computed the MODWT
modwt_result <- modwt(abund11, wf = 'la8', boundary = "periodic", n.levels = 4)

# Calculate the magnitude of coefficients for each level
coeff_magnitudes <- sapply(modwt_result, function(level) sqrt(sum(level^2)))

# Find the level with the maximum magnitude
selected_level <- which.max(coeff_magnitudes)

# Use the coefficients from the selected level for further analysis
selected_coefficients <- modwt_result[[selected_level]]

##Gain function----
#### No va la lliberia!!
# install.packages('wmsta', repos = "https://CRAN.R-project.org/package=wmtsa")
# library(wmtsa)
# 
# FreqResponsFunctionS8 <- wavGain(wavelet = "s8", n.levels = 4, n.fft = 1024, normalize = TRUE)
# SqrGainHigh <- t(matrix(FreqResponsFunctionS8$sqrgain$high, nrow = 4, byrow = F)) 
# colnames(SqrGainHigh) <- c("Level1", "Level2", "Level3", "Level4") 
# SqrGainLow <- t(matrix(FreqResponsFunctionS8$sqrgain$low, nrow = 4, byrow = F)) 
# colnames(SqrGainLow) <- c("Level1", "Level2", "Level3", "Level4") 
# par(mar = c(3, 4, 1, 0.1)) 
# plot(1:512, type = "n", axes = F, ylab = "Squared Gain Function", xlab = "", xlim = c(1, + 512), ylim = c(0, 1), font.lab = 1) 
# polygon(x = c(1/16, 1/8, 1/8, 1/16) * 1024, y = c(-0.1, -0.1, 1.1, 1.1), col = "lightgray", 
#         border = NA) 
# polygon(x = c(1/4, 1/2, 1/2, 1/4) * 1024, y = c(-0.1, -0.1, 1.1, 1.1), col = "lightgray", 
#         border = NA) 
# polygon(x = c(0, 1/32, 1/32, 0) * 1024, y = c(-0.1, -0.1, 1.1, 1.1), col = "lightgray", 
#         border = NA) 
# box(bty = "L", lwd = 1) 
# matplot(SqrGainHigh[1:512, ], type = "l", lty = "solid", axes = F, ylab = "Squared Gain", 
#         xlab = "Scales [month]", add = T) > lines(1:512, SqrGainLow[1:512, 4], col = "blue") 
# abline(v = c(1/24, 1/12, 1/6, 1/3) * 1024, lty = "solid", lwd = 1) 
# axis(2, font.axis = 1.3, las = 2) 
# mtext(side = 3, line = 0.25, outer = F, text = c(32, 16, 8, 4, 2), at = c(1/32, 1/16,  1/8, 1/4, 1/2) * 1024, font = 1, cex = 0.8) 
# text(x = c(1/72, 1/24, 1/12, 1/6, 1/3) * 1024, y = rep(0.5, 4), labels = c("s4", "d4", + "d3", "d2", "d1"), font = 1, cex = 0.75, pos = 4, offset = 0.1) 
# axis(1, font.axis = 1, at = c(1/24, 1/12, 1/6, 1/3) * 1024, labels = c(24, 12, 6, 3), 
#      cex.axis = 1) 
# mtext(side = 1, line = 1.6, outer = F, at = 256, text = "Periodicity (# obs.)", font = 1)



# WAVELETS ANALYSIS FOR ALL THE ASVs AT THE SAME TIME----
## PERFORM WAVELETS ANALYSIS APPLYING THE BRICK WALL FUNCTION WHICH REMOVES ALL VALUES AFFECTED BY MARGINS EFFECTS----

## upload data
### we work with two different datasets one for FL and the other for PA
wavelet_3_df <- read.csv2('data/wavelet_3_df_deconstand.csv', sep = ',') |>
  as_tibble() |>
  dplyr::select(-X)

wavelet_02_df <- read.csv2('data/wavelet_02_df_deconstand.csv') |>
  as_tibble() |>
  dplyr::select(-X)

  #### 4 steps APPLY BRICK WALL FUNCTION (REMOVE ALL SAMPLES AFFECTED BY THE MARGINS EFFECT)
### 1. modwt computation----
  ## FL fraction
  modwt_results_02 <- wavelet_02_df |>
  arrange(decimal_date) |>
    group_by(asv_num) %>%
    summarize(modwt_result = list(modwt.function(abundance_value)))
  
  ## PA fraction
  modwt_results_3 <- wavelet_3_df |>
    arrange(decimal_date) |>
    group_by(asv_num) %>%
    dplyr::summarize(modwt_result = list(modwt.function(zclr)))
  
### 2. e-folding-----
  ###### commmon for all 
  x <- rep(0, 10001) 
  x[5000] <- 1 
  n.levels <- 4 
  len <- length(m_02$sample_id) ##120 (length of my dataset)
  temp <- phase.shift(modwt(x, n.levels = n.levels, wf = "la8"), wf = "la8")
  
  ## The positions to the left and to the right of the maximal influence of this spike are recorded in a matrix (left, right) together with the 
  ## position of the maximum itself (top).
  waveExtremes <- matrix(nrow = 3, ncol = n.levels + 1) 
  colnames(waveExtremes) <- c(paste("d", 1:n.levels, sep = ""), paste("s", n.levels,  sep = "")) 
  rownames(waveExtremes) <- c("left", "right", "top")
  
  ## The distance to the maximum from both sides of the influence is determined as 1/e2 times the maximum within a specific coefficient vector.
  for (i in 1:(n.levels + 1)) waveExtremes[, i] <- c(range(which(abs(temp[[i]]) 
                                                                 >= max(abs(temp[[i]]))/(exp(2)))), which.max(abs(temp[[i]])))
  
  ## The positions (waveExtremes) are used to calculate the distances to the left and to the right of the influence maximum. 
  ## The distance to the left of the maximum is called "right" because it will serve to calculate the distance at the end of the series.
  
  boundaries <- data.frame(end = len - (waveExtremes[3, ] - waveExtremes[1, ]), 
                           start = waveExtremes[2, ] - waveExtremes[3, ])
  
  ##### specific for each ASV
  asv_num_index <- i
  modwt_results <- modwt_results_02$modwt_result[[asv_num_index]]
  
  i = 1 #i defined it because it was not working, but it should be fixed.
  for (j in 1:(n.levels + 1)) { 
    is.na(modwt_results[[i]]) <- c(1:boundaries$start[i], boundaries$end[i]:length(modwt_results[[i]])) 
    
  }
  
### 3. Visualize the results obtained from the modwt transfomation ------
#  relation_asvs_list <-  modwt_results_02$asv_num |>
#     as_tibble() |>
#     dplyr::mutate(num = row_number())
#   
#   ##create a dataframe with the results from the wavelet transfromation
#   
#   d1 <- modwt_results_02$modwt_result[[1]][[1]] |>
#     as_tibble_col(column_name = 'd1') |>
#     dplyr::mutate(sample_num = row_number())
#     
#  d2 <-  modwt_results_02$modwt_result[[1]][[2]] |>
#     as_tibble_col(column_name = 'd2')
#   
#   d3 <-modwt_results_02$modwt_result[[1]][[3]] |>
#     as_tibble_col(column_name = 'd3')
#   
#   d4 <- modwt_results_02$modwt_result[[1]][[4]] |>
#     as_tibble_col(column_name = 'd4')
#   
#   s4 <- modwt_results_02$modwt_result[[1]][[5]] |>
#     as_tibble_col(column_name = 's4')
#   
# asv_x <-  relation_asvs_list |>
#     dplyr::filter(num == 1) 
#   
#   # Define the number of replications
#   num_replications <- 120
#   
#   # Replicate the tibble
#   replicated_tibbles <- replicate(num_replications, asv_x, simplify = FALSE)
#   
#   # Combine the replicated tibbles into one
#   asv_num_tibble <- dplyr::bind_rows(replicated_tibbles) |>
#     dplyr::select(-num) |>
#     rename(asv_num = value)
#   
#   wavelets_result_tibble <- bind_cols(d1, d2, d3, d4, s4, asv_num_tibble)
#   
#  decimal_date_tibble <-  wavelet_02_df %$%
#     decimal_date |>
#     unique() |>
#     as_tibble_col(column_name = 'decimal_date') |>
#     dplyr::mutate(sample_num = row_number())
#   
#  tax <- asv_tab_all_bloo_z_tax |>
#    dplyr::select(asv_num, phylum, class, order, family, genus) |>
#    distinct()
#  
#  wavelets_result_tibble_tax  <- wavelets_result_tibble |>
#     pivot_longer(cols = !c(asv_num, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
#      left_join( decimal_date_tibble) |>
#      left_join(tax, by = asv_num)
#  
#  wavelets_result_tibble_tax  |>
#     ggplot(aes(decimal_date, wavelets_result))+
#     geom_col()+
#     ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_num), ' ', unique(wavelets_result_tibble_tax$family))) +
#     labs(x = 'Decimal date', y = 'Wavelets results')+
#     facet_grid(vars(wavelets_transformation))+
#      scale_x_continuous(expand = c(0,0))+
#      theme_bw()+
#      theme(panel.grid = element_blank(), strip.background = element_blank(),
#            aspect.ratio = 4/10,
#            text = element_text(size = 5))
#   
#  ##test for another ASV
#  d1 <- modwt_results_02$modwt_result[[2]][[1]] |>
#    as_tibble_col(column_name = 'd1') |>
#    dplyr::mutate(sample_num = row_number())
#  
#  d2 <-  modwt_results_02$modwt_result[[2]][[2]] |>
#    as_tibble_col(column_name = 'd2')
#  
#  d3 <-modwt_results_02$modwt_result[[2]][[3]] |>
#    as_tibble_col(column_name = 'd3')
#  
#  d4 <- modwt_results_02$modwt_result[[2]][[4]] |>
#    as_tibble_col(column_name = 'd4')
#  
#  s4 <- modwt_results_02$modwt_result[[2]][[5]] |>
#    as_tibble_col(column_name = 's4')
#  
#  asv_x <-  relation_asvs_list |>
#    dplyr::filter(num == 2) 
#  
#  # Define the number of replications
#  num_replications <- 120
#  
#  # Replicate the tibble
#  replicated_tibbles <- replicate(num_replications, asv_x, simplify = FALSE)
#  
#  # Combine the replicated tibbles into one
#  asv_num_tibble <- dplyr::bind_rows(replicated_tibbles) |>
#    dplyr::select(-num) |>
#    rename(asv_num = value) 
#  
#  wavelets_result_tibble <- bind_cols(d1, d2, d3, d4, s4, asv_num_tibble) 
#  
#  wavelets_result_tibble_tax  <- wavelets_result_tibble |>
#    pivot_longer(cols = !c(asv_num, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
#    left_join( decimal_date_tibble) |>
#    left_join(tax, by = asv_num)
#  
#  wavelets_result_tibble_tax  |>
#    ggplot(aes(decimal_date, wavelets_result))+
#    geom_col()+
#    ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_num), ' ', unique(wavelets_result_tibble_tax$family))) +
#    labs(x = 'Decimal date', y = 'Wavelets results')+
#    facet_grid(vars(wavelets_transformation))+
#    scale_x_continuous(expand = c(0,0))+
#    theme_bw()+
#    theme(panel.grid = element_blank(), strip.background = element_blank(),
#          aspect.ratio = 4/10,
#          text = element_text(size = 5))
 
 #### try to extract the results for all ASVs at the same time-----
 
 ###### This part works but still I need to do it for all ASVs.
 # Initialize lists to store tibbles
 # d_tibbles <- list()
 # 
 # # Loop through modwt_result[[2]]
 # for (i in 1:5) {
 #   d_tibbles[[i]] <- modwt_results_02$modwt_result[[2]][[i]] |>
 #     as_tibble_col(column_name = paste0("d", i))
 # }
 # 
 # # Combine the tibbles
 # wavelets_result_tibble <- bind_cols(d_tibbles) |>
 #   dplyr::mutate(sample_num = row_number())
 # 
 # # Replicate the asv_x tibble
 # num_replications <- 120
 # replicated_tibbles <- replicate(num_replications, asv_x, simplify = FALSE)
 # 
 # # Combine the replicated tibbles into one
 # asv_num_tibble <- dplyr::bind_rows(replicated_tibbles) |>
 #   dplyr::select(-num) |>
 #   rename(asv_num = value)
 # 
 # # Combine the wavelets_result_tibble and asv_num_tibble
 # final_tibble <- bind_cols(wavelets_result_tibble, asv_num_tibble)
 # 
 # wavelets_result_tibble_tax  <-  final_tibble |>
 #   pivot_longer(cols = !c(asv_num, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
 #   left_join( decimal_date_tibble) |>
 #   left_join(tax, by = asv_num)
 # 
 # wavelets_result_tibble_tax  |>
 #   ggplot(aes(decimal_date, wavelets_result))+
 #   geom_col()+
 #   ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_num), ' ', unique(wavelets_result_tibble_tax$family))) +
 #   labs(x = 'Decimal date', y = 'Wavelets results')+
 #   facet_grid(vars(wavelets_transformation))+
 #   scale_x_continuous(expand = c(0,0))+
 #   theme_bw()+
 #   theme(panel.grid = element_blank(), strip.background = element_blank(),
 #         aspect.ratio = 4/10,
 #         text = element_text(size = 5))
 
 ### I extract the wavelets results at the same time for ALL ASVs that I have in modwt_results 
 ## FL FRACTION----

 # Initialize a list to store tibbles
 all_tibbles <- list()
 
 # Loop through the rows of modwt_results_02
 for (i in seq_len(nrow(modwt_results_02))) {
   # Extract the current row
   current_row <- modwt_results_02[i, ]
   current_asv_row.number <- current_row |>
     dplyr::mutate(row_number_asv = row_number()) |>
     dplyr::select(row_number_asv) 
   
   current_asv_row.number <- current_asv_row.number$row_number_asv[1]
   
   # Extract asv_num and modwt_result list
   current_asv_num <- current_row$asv_num
   current_modwt_result <- current_row$modwt_result
   
   # Initialize a list to store tibbles for the current asv_num
   d_tibbles <- list()
   
   # Loop through modwt_result
   for (j in 1:5) {
     # Extract the asv_num tibble
     d_tibbles[[j]] <- current_modwt_result[[current_asv_row.number]][[j]] %>%
       as_tibble_col(column_name = paste0("d", j))
   }
   
   ##create a column with the asv_num
   current_row <- current_asv_num %>%
     as_tibble() |>
     mutate(asv_name = list(rep(value, each = 120))) |>
     unnest(asv_name) |>
     dplyr::select(-value)
   
   # Combine the tibbles for the current asv_num
   if (length(d_tibbles) > 0) {
     all_tibbles[[current_asv_num]] <- bind_cols(d_tibbles) %>%
       dplyr::mutate(sample_num = row_number()) %>%
       bind_cols(current_row)
   }
 }
 
 # Combine all the tibbles into one
 final_tibble <- bind_rows(all_tibbles)
 
  decimal_date_tibble <-  wavelet_02_df %$%
     decimal_date |>
     unique() |>
     as_tibble_col(column_name = 'decimal_date') |>
     dplyr::mutate(sample_num = row_number())
  
   tax <- asv_tab_all_bloo_z_tax |>
     dplyr::select(asv_num, phylum, class, order, family, genus) |>
     distinct()
 
 wavelets_result_tibble_tax_02  <-  final_tibble |>
   pivot_longer(cols = !c(asv_name, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
   left_join( decimal_date_tibble) |>
   left_join(tax, by = c('asv_name' = 'asv_num')) |>
   rename(asv_num = asv_name) |>
   dplyr::mutate(wavelets_transformation = str_replace(wavelets_transformation, 'd5', 's5'))
 
 ## PA fraction----
 # Initialize a list to store tibbles
 all_tibbles <- list()
 
 # Loop through the rows of modwt_results_3
 for (i in seq_len(nrow(modwt_results_3))) {
   # Extract the current row
   current_row <- modwt_results_3[i, ]
   current_asv_row.number <- current_row |>
     dplyr::mutate(row_number_asv = row_number()) |>
     dplyr::select(row_number_asv) 
   
   current_asv_row.number <- current_asv_row.number$row_number_asv[1]
   
   # Extract asv_num and modwt_result list
   current_asv_num <- current_row$asv_num
   current_modwt_result <- current_row$modwt_result
   
   # Initialize a list to store tibbles for the current asv_num
   d_tibbles <- list()
   
   # Loop through modwt_result
   for (j in 1:5) {
     # Extract the asv_num tibble
     d_tibbles[[j]] <- current_modwt_result[[current_asv_row.number]][[j]] %>%
       as_tibble_col(column_name = paste0("d", j))
   }
   
   ##create a column with the asv_num
   current_row <- current_asv_num %>%
     as_tibble() |>
     mutate(asv_name = list(rep(value, each = 120))) |>
     unnest(asv_name) |>
     dplyr::select(-value)
   
   # Combine the tibbles for the current asv_num
   if (length(d_tibbles) > 0) {
     all_tibbles[[current_asv_num]] <- bind_cols(d_tibbles) %>%
       dplyr::mutate(sample_num = row_number()) %>%
       bind_cols(current_row)
   }
 }
 
 # Combine all the tibbles into one
 final_tibble <- bind_rows(all_tibbles)
 
 decimal_date_tibble <-  wavelet_3_df %$%
   decimal_date |>
   unique() |>
   as_tibble_col(column_name = 'decimal_date') |>
   dplyr::mutate(sample_num = row_number())
 
 tax <- asv_tab_all_bloo_z_tax |>
   dplyr::select(asv_num, phylum, class, order, family, genus) |>
   distinct()
 
 wavelets_result_tibble_tax_3  <-  final_tibble |>
   pivot_longer(cols = !c(asv_name, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
   left_join( decimal_date_tibble) |>
   left_join(tax, by = c('asv_name' = 'asv_num')) |>
   rename(asv_num = asv_name) |>
   dplyr::mutate(wavelets_transformation = str_replace(wavelets_transformation, 'd5', 's5'))
 
### PLOT THE WAVELETS TRANSFROMATINS COMPUTED (visually inspect them, to be sure of the results)----
 ## FL FRACTION----
  ##i divide them so that i can see the plots better
 bloo_02
 
 wavelets_result_tibble_tax_02 |>
   dplyr::filter(asv_name %in%  (bloo_02 |>
                                   slice_head(n = 10) %$%
                                   value)) |>
   ggplot(aes(decimal_date, wavelets_result))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_name), ' ', unique( wavelets_result_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~asv_name)+
   scale_x_continuous(expand = c(0,0))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         #aspect.ratio = 4/10,
         text = element_text(size = 5))

 wavelets_result_tibble_tax_02 |>
   dplyr::filter(asv_num %in%  (bloo_02 |>
                                   slice_tail(n = 11) %$%
                                   value)) |>
   ggplot(aes(decimal_date, wavelets_result))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_name), ' ', unique(wavelets_result_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~asv_name)+
   scale_x_continuous(expand = c(0,0))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         #aspect.ratio = 4/10,
         text = element_text(size = 5))
 
 ## PA FRACTION---- 
 bloo_3$value #45 in total
 
 wavelets_result_tibble_tax_3 |>
   dplyr::filter(asv_name %in%  (bloo_3 |>
                                   slice_head(n = 20) %$%
                                   value)) |>
   ggplot(aes(decimal_date, wavelets_result))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_name), ' ', unique( wavelets_result_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~asv_name)+
   scale_x_continuous(expand = c(0,0))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         #aspect.ratio = 4/10,
         text = element_text(size = 5))
 
 wavelets_result_tibble_tax_3 |>
   dplyr::filter(asv_name %in%  (bloo_3 |>
                                   slice_tail(n = 20) %$%
                                   value)) |>
   ggplot(aes(decimal_date, wavelets_result))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_name), ' ', unique(wavelets_result_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~asv_name)+
   scale_x_continuous(expand = c(0,0))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         #aspect.ratio = 4/10,
         text = element_text(size = 5))
 

 ### create a loop to save all the wavelets transformations computed and visually inspect them
 # Create a folder to save the PDF files
 dir.create("results/figures/wavelets_plots", showWarnings = FALSE)
 
 ### FL----
 # Get unique asv_num values
 asv_nums <- unique(wavelets_result_tibble_tax_02$asv_num)
 
 for (asv_num in asv_nums) {
   # Filter data for the current asv_num
   plot_data <- wavelets_result_tibble_tax_02 %>%
     dplyr::filter(asv_num == !!asv_num)  # Use !! to unquote asv_num
   
   title_plot <- plot_data |>
     group_by(asv_num, family) |>
     dplyr::reframe(title_ed = paste0(asv_num,', ', family, ', ', order)) |>
     distinct() |>
     dplyr::select(title_ed) |>
     as.character()
   
   # Create the plot
   p <- ggplot(plot_data, aes(decimal_date, wavelets_result)) +
     geom_col() +
     labs(x = 'Decimal date', y = 'Wavelets results') +
     facet_grid(vars(wavelets_transformation)) +
     scale_x_continuous(expand = c(0,0)) +
     labs(title = title_plot)+
     theme_bw() +
     theme(panel.grid = element_blank(), strip.background = element_blank(),
           aspect.ratio = 4/10,
           text = element_text(size = 5))
   
   # Save the plot as a PDF file
   pdf_file <- paste0("results/figures/wavelets_plots/", asv_num, "_plot_02.pdf")
   ggsave(pdf_file, p, width = 6, height = 3, units = "in")
   
   # Print a message indicating the plot has been saved
   cat("Plot for", asv_num, "saved as", pdf_file, "\n")
 }
 
 
 ## PA----
 # Get unique asv_num values
 asv_nums <- unique(wavelets_result_tibble_tax_3$asv_num)
 
 for (asv_num in asv_nums) {
   # Filter data for the current asv_num
   plot_data <- wavelets_result_tibble_tax_3 %>%
     dplyr::filter(asv_num == !!asv_num)  # Use !! to unquote asv_num
   
   title_plot <- plot_data |>
     group_by(asv_num, family) |>
     dplyr::reframe(title_ed = paste0(asv_num,', ', family, ', ', order)) |>
     distinct() |>
     dplyr::select(title_ed) |>
     as.character()
   
   # Create the plot
   p <- ggplot(plot_data, aes(decimal_date, wavelets_result)) +
     geom_col() +
     labs(x = 'Decimal date', y = 'Wavelets results') +
     facet_grid(vars(wavelets_transformation)) +
     scale_x_continuous(expand = c(0,0)) +
     labs(title = title_plot)+
     theme_bw() +
     theme(panel.grid = element_blank(), strip.background = element_blank(),
           aspect.ratio = 4/10,
           text = element_text(size = 5))
   
   # Save the plot as a PDF file
   pdf_file <- paste0("results/figures/wavelets_plots/", asv_num, "_plot_3.pdf")
   ggsave(pdf_file, p, width = 6, height = 3, units = "in")
   
   # Print a message indicating the plot has been saved
   cat("Plot for", asv_num, "saved as", pdf_file, "\n")
 }
 
 
##### original code------
  # par(mfrow = c(n.levels + 1, 1), mar = c(2, 3, 0.5, 0.2))
  # # Assuming n.levels is defined before this code snippet
  # 
  # # Create an empty list to store your plots
  # plotList <- vector("list", n.levels)
  # 
  # for (i in 1:n.levels) { 
  #   plot(as.vector(time(abund11)), modwtasv11_biased[[i]], xlab = "", ylab = "", type = "h", 
  #        bty = "n", axes = F, ylim = range(modwtasv11[1:n.levels], na.rm = TRUE)) 
  #   axis(2, cex.axis = 1.5, las = 2) +
  #     abline(v = seq(start(abund194)[1], end(abund194)[1], by = 12), col = "gray", 
  #            lty = "dashed") + abline(h = 0) 
  #   abline(v = time(abund194)[range(which(!is.na(modwtasv11[[i]])))], lwd = 2, 
  #          col = "red") + 
  #     abline(v = time(abund11)[c(boundaries$start[i], boundaries$end[i])], col = "blue", 
  #            lwd = 3, lty = "dashed") +
  #     mtext(side = 4, line = -2, cex = 1.2, outer = FALSE, at = par("usr")[4], text = paste("d", i, sep = ""), las = 2) 
  #   
  #   # Save each plot to the list
  #   plotList[[i]] <- recordPlot()
  # }
  # 
  # # Set the title for the entire series of plots
  # title(main = "ASV 194")
  # 
  # # Replay and display the saved plots
  # for (i in 1:n.levels) {
  #   replayPlot(plotList[[i]])
  # }
  # 
  # plot(as.vector(time(abund11)), modwtasv11_biased$s4, xlab = "", ylab = "", type = "h",
  #      bty = "n", axes = F) > axis(1, cex.axis = 1.5) 
  # axis(2, cex.axis = 1.5, las = 2) 
  # mtext(side = 4, line = -2, cex = 1.2, outer = F, at = par("usr")[4], text = paste("s4", 
  #                                                                                   sep = ""), las = 2) +
  #   abline(v = time(abund11)[range(which(!is.na(modwtasv11[[5]])))], lwd = 2, 
  #          col = "red") 
  # 
  # title(main = "ASV 11")
  # abline(v = time(abund11)[c(boundaries$start[5], boundaries$end[5])], col = "blue", 
  #        lwd = 3, lty = "dashed")
  
 
### 4. Magnitude (coefficients) observe were they got the highest coefficients for the wavelet analysis--------

 #### Wavelet coefficient magnitude indicates how strongly the data are correlated with the mother wavelet at a given frequency and distance.
 
bloo_02_type <- wavelets_result_tibble_tax_02 %>%
   group_by(asv_num, wavelets_transformation) %>%
   dplyr::filter(!is.na(wavelets_result)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result^2))) |>  # Calculate the magnitude of coefficients for each level
   group_by(asv_num) |>
   top_n(1, wt = coefficients) |>  # Find the level with the maximum magnitude
   rename(max_coeff = wavelets_transformation) |>
   dplyr::mutate(bloomer_type = case_when(max_coeff == 's5' ~ 'inter-annual',
                                          max_coeff == 'd1' ~ 'fine-scale',
                                          max_coeff == 'd2' ~ 'half-yearly', 
                                          max_coeff == 'd3' ~ 'seasonal', 
                                          max_coeff == 'd4' ~ 'year-to-year')) |>
   left_join(tax, by = 'asv_num')
 
 bloo_02_type |>
   ungroup() |>
   dplyr::filter(bloomer_type == 'seasonal') |>
   dplyr::reframe(n = n())
 
 3/19 ## 15.58% of the free living blooms seem to have a seasonal pattern
 
 bloo_02_type_tax <-  bloo_02_type |>
   group_by(bloomer_type, family ) |>
   dplyr::reframe(n = n()) |>
   dplyr::mutate(fraction = '0.2')
 
 bloo_3_type <- wavelets_result_tibble_tax_3 %>%
   group_by(asv_num, wavelets_transformation) %>%
   dplyr::filter(!is.na(wavelets_result)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result^2))) |>  # Calculate the magnitude of coefficients for each level
   group_by(asv_num) |>
   top_n(1, wt = coefficients) |>  # Find the level with the maximum magnitude
   rename(max_coeff = wavelets_transformation) |>
   dplyr::mutate(bloomer_type = case_when(max_coeff == 's5' ~ 'inter-annual',
                                          max_coeff == 'd1' ~ 'fine-scale',
                                          max_coeff == 'd2' ~ 'half-yearly', 
                                          max_coeff == 'd3' ~ 'seasonal', 
                                          max_coeff == 'd4' ~ 'year-to-year')) |>
   left_join(tax, by = 'asv_num')
 
 bloo_3_type |>
   ungroup() |>
   dplyr::filter(bloomer_type == 'seasonal') |>
   dplyr::reframe(n = n())
 
 6/45 ## 13.33% of the particle attached seem to have a seasonal pattern
 
 bloo_3_type_tax <-  bloo_3_type |>
   group_by(bloomer_type, family) |>
   dplyr::reframe(n = n()) |>
   dplyr::mutate(fraction = '3')
 
 ## I create a table that summarizes the different type of bloomers that we have, observe their taxonomy too.
 bloo_02_type_tax |>
   bind_rows(bloo_3_type_tax) |>
   group_by( bloomer_type, fraction) |>
   dplyr::mutate(n_total = sum(n)) |>
   ggplot(aes(bloomer_type, n_total))+
     geom_col(aes(y = n, fill = family))+
     scale_fill_manual(values = palette_family_assigned_bloo)+
   facet_wrap(vars(fraction), labeller = labs_fraction)+
   scale_y_continuous(expand = c(0,0))+
   geom_text(aes(bloomer_type, 1, label = n_total),
             color = "black", nudge_y = 10) +
   theme_bw()+
   #coord_flip()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', axis.ticks.x = element_blank(),
         panel.border = element_blank(), text = element_text(size = 5))
 
## in relatives  
 bloo_all_type_tax <-  bloo_02_type_tax |>
   bind_rows(bloo_3_type_tax) |>
   group_by( bloomer_type, fraction) |>
   dplyr::mutate(n_total = sum(n)) |>
   dplyr::group_by(fraction) |>
   dplyr::mutate(n_blooms_fraction = sum(n)) |>
   dplyr::mutate(n_blooms_fraction_perc = n/n_blooms_fraction) |>
   dplyr::mutate(n_blooms_fraction_perc_type = n/n_total) 
 
 bloo_all_type_tax$bloomer_type <- factor( bloo_all_type_tax$bloomer_type,
                                           levels = c('fine-scale', 'seasonal', 'half-yearly', 'inter-annual'))
 
 ### by family
 bloo_all_type_tax |>
   ggplot(aes(bloomer_type, 1))+
   geom_col(aes(y = n_blooms_fraction_perc_type, fill = family))+
   scale_fill_manual(values = palette_family_assigned_bloo)+
   facet_wrap(vars(fraction), labeller = labs_fraction)+
   scale_y_continuous(expand = c(0,0), labels = percent_format())+
   labs(fill = 'Family', y ='Number of bloomers', x = 'Bloomer period')+
   geom_text(aes(bloomer_type, 1, label = n_total),
             color = "black", nudge_y = -0.3) +
   theme_bw()+
   #coord_flip()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', axis.ticks.x = element_blank(),
         panel.border = element_blank(), text = element_text(size = 5))
 
### by order ----
 bloo_02_type_tax_order <-  bloo_02_type |>
   group_by(bloomer_type, order ) |>
   dplyr::reframe(n = n()) |>
   dplyr::mutate(fraction = '0.2')
 
 
 bloo_3_type_tax_order <-  bloo_3_type |>
   group_by(bloomer_type, order) |>
   dplyr::reframe(n = n()) |>
   dplyr::mutate(fraction = '3')
 
 
 bloo_all_type_tax_order <-  bloo_02_type_tax_order |>
   bind_rows(bloo_3_type_tax_order) |>
   group_by( bloomer_type, fraction) |>
   dplyr::mutate(n_total = sum(n)) |>
   dplyr::group_by(fraction) |>
   dplyr::mutate(n_blooms_fraction = sum(n)) |>
   dplyr::mutate(n_blooms_fraction_perc = n/n_blooms_fraction) |>
   dplyr::mutate(n_blooms_fraction_perc_type = n/n_total) 
 
 bloo_all_type_tax_order$bloomer_type <- factor(bloo_all_type_tax_order$bloomer_type,
                                           levels = c('fine-scale', 'seasonal', 'half-yearly', 'inter-annual'))
 
 bloo_all_type_tax_order |>
   ggplot(aes(bloomer_type, 1))+
   geom_col(aes(y = n_blooms_fraction_perc_type, fill = order))+
   scale_fill_manual(values = palette_order_assigned_bloo)+
   facet_wrap(vars(fraction), labeller = labs_fraction)+
   scale_y_continuous(expand = c(0,0), labels = percent_format())+
   geom_text(aes(bloomer_type, 1, label = n_total),
             color = "black", nudge_y = -0.3) +
   labs(fill = 'Family', y ='Number of bloomers', x = 'Bloomer period')+
   theme_bw()+
   #coord_flip()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', axis.ticks.x = element_blank(),
         panel.border = element_blank(), text = element_text(size = 5))
 
 
 
## PERFORM WAVELETS ANALYSIS THIS TIME WITHOUT APPLYING THE BRICK WALL FUNCTION WHICH REMOVES VALUES AFFECTED BY MARGINS EFFECTS-----
 #### 4 steps 
 ### 1. modwt computation----
 
 ##n levels must be a number less or equal to log(length(x))
 log(120, base = 2) #6.9 max levels we could use.
 
 modwt.function.biased <- function(abundance){
   modwt_result <-  abundance |>
     modwt( wf = 'la8', boundary = "periodic", n.levels = 4) |> #If boundary=="periodic" the default TRUE, then the vector you decompose is assumed to be periodic on its defined interval,
     #brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
     phase.shift(wf = 'la8')
   return(modwt_result)
 }
 
 ## FL fraction
 modwt_results_02_biased <- wavelet_02_df |>
   arrange(decimal_date) |>
   group_by(asv_num) %>%
   summarize(modwt_result = list(modwt.function.biased(abundance_value)))
 
 ## PA fraction
 modwt_results_3_biased <- wavelet_3_df |>
   arrange(decimal_date) |>
   group_by(asv_num) %>%
   dplyr::summarize(modwt_result = list(modwt.function.biased(rclr)))
 
 ### 2. e-folding-----
 ###### commmon for all 
 x <- rep(0, 10001) 
 x[5000] <- 1 
 n.levels <- 4 
 len <- length(m_02$sample_id) ##120 (length of my dataset) we do not use 3 because they have 3 datapoints less for the metadata 
 temp <- phase.shift(modwt(x, n.levels = n.levels, wf = "la8"), wf = "la8")
 
 ## The positions to the left and to the right of the maximal influence of this spike are recorded in a matrix (left, right) together with the 
 ## position of the maximum itself (top).
 waveExtremes <- matrix(nrow = 3, ncol = n.levels + 1) 
 colnames(waveExtremes) <- c(paste("d", 1:n.levels, sep = ""), paste("s", n.levels,  sep = "")) 
 rownames(waveExtremes) <- c("left", "right", "top")
 
 ## The distance to the maximum from both sides of the influence is determined as 1/e2 times the maximum within a specific coefficient vector.
 for (i in 1:(n.levels + 1)) waveExtremes[, i] <- c(range(which(abs(temp[[i]]) 
                                                                >= max(abs(temp[[i]]))/(exp(2)))), which.max(abs(temp[[i]])))
 
 ## The positions (waveExtremes) are used to calculate the distances to the left and to the right of the influence maximum. 
 ## The distance to the left of the maximum is called "right" because it will serve to calculate the distance at the end of the series.
 
 boundaries <- data.frame(end = len - (waveExtremes[3, ] - waveExtremes[1, ]), 
                          start = waveExtremes[2, ] - waveExtremes[3, ])
 
 ##### General for the whole dataset (this part should be fixed but maybe I do not needed because I remove the boundaries after)
 # asv_num_index <- i
 # modwt_results <- modwt_results_02$modwt_result[[asv_num_index]]
 # 
 # i = 1 #i defined it because it was not working, but it should be fixed.
 # for (j in 1:(n.levels + 1)) { 
 #   is.na(modwt_results[[i]]) <- c(1:boundaries$start[i], boundaries$end[i]:length(modwt_results[[i]])) 
 #   
 # }
 
 ### 3. Visualize the results obtained from the modwt transfomation ------
 ### I extract the wavelets results at the same time for ALL ASVs that I have in modwt_results 
 
 ## FL FRACTION----
 
 # Initialize a list to store tibbles
 all_tibbles <- list()
 
 # Loop through the rows of modwt_results_02
 for (i in seq_len(nrow(modwt_results_02_biased))) {
   # Extract the current row
   current_row <- modwt_results_02_biased[i, ]
   current_asv_row.number <- current_row |>
     dplyr::mutate(row_number_asv = row_number()) |>
     dplyr::select(row_number_asv) 
   
   current_asv_row.number <- current_asv_row.number$row_number_asv[1]
   
   # Extract asv_num and modwt_result list
   current_asv_num <- current_row$asv_num
   current_modwt_result <- current_row$modwt_result
   
   # Initialize a list to store tibbles for the current asv_num
   d_tibbles <- list()
   
   # Loop through modwt_result
   for (j in 1:5) {
     # Extract the asv_num tibble
     d_tibbles[[j]] <- current_modwt_result[[current_asv_row.number]][[j]] %>%
       as_tibble_col(column_name = paste0("d", j))
   }
   
   ##create a column with the asv_num
   current_row <- current_asv_num %>%
     as_tibble() |>
     mutate(asv_name = list(rep(value, each = 120))) |>
     unnest(asv_name) |>
     dplyr::select(-value)
   
   # Combine the tibbles for the current asv_num
   if (length(d_tibbles) > 0) {
     all_tibbles[[current_asv_num]] <- bind_cols(d_tibbles) %>%
       dplyr::mutate(sample_num = row_number()) %>%
       bind_cols(current_row)
   }
 }
 
 # Combine all the tibbles into one
 final_tibble <- bind_rows(all_tibbles)
 
 decimal_date_tibble <-  wavelet_02_df %$%
   decimal_date |>
   unique() |>
   as_tibble_col(column_name = 'decimal_date') |>
   dplyr::mutate(sample_num = row_number())
 
 tax <- asv_tab_all_bloo_z_tax |>
   dplyr::select(asv_num, phylum, class, order, family, genus) |>
   distinct()
 
 wavelets_result_tibble_tax_02_biased  <-  final_tibble |>
   pivot_longer(cols = !c(asv_name, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
   left_join( decimal_date_tibble) |>
   left_join(tax, by = c('asv_name' = 'asv_num')) |>
   rename(asv_num = asv_name) |>
   dplyr::mutate(wavelets_transformation = str_replace(wavelets_transformation, 'd5', 's4'))
 
 ## I remove the most afected samples by the boundaries it is less biased but still more biased than when applying the brick wall function
 ## as we increase the signal the wavelet gets more affected by the margin effect
 boundaries 
 
 wavelets_result_tibble_tax_02_biased_red <-  wavelets_result_tibble_tax_02_biased |>
   dplyr::mutate(wavelets_result_ed = case_when(wavelets_transformation == 'd1' &
                                                  sample_num %in% c(1, 119, 120) ~ 'NA',
                                                wavelets_transformation == 'd2' &
                                                  sample_num %in% c(1,2,3,  117, 118, 119, 120) ~ 'NA',
                                                wavelets_transformation == 'd3' &
                                                  sample_num %in% c(1,2,3,4,5,6, 113,114,115,116,  117, 118, 119, 120) ~ 'NA',
                                                wavelets_transformation == 'd4' &
                                                  sample_num %in% c(1,2,3,4,5,6,7,8,9,10,11,12, 105,106,107,108,109,110,111,112,113,114,
                                                                    115,116,117, 118, 119, 120) ~ 'NA',
                                                wavelets_transformation == 's4' &
                                                  sample_num %in% c(1,2,3,4,5,6,7,8,9,10,11,12, 
                                                                    13,14,15,16,17,
                                                                    108,109,110,111,112,113,114,
                                                                    115,116,117, 118, 119, 120) ~ 'NA',
                                                TRUE ~ as.character(wavelets_result)))
 
 ## PA fraction----
 # Initialize a list to store tibbles
 all_tibbles <- list()
 
 # Loop through the rows of modwt_results_3
 for (i in seq_len(nrow(modwt_results_3_biased))) {
   # Extract the current row
   current_row <- modwt_results_3_biased[i, ]
   current_asv_row.number <- current_row |>
     dplyr::mutate(row_number_asv = row_number()) |>
     dplyr::select(row_number_asv) 
   
   current_asv_row.number <- current_asv_row.number$row_number_asv[1]
   
   # Extract asv_num and modwt_result list
   current_asv_num <- current_row$asv_num
   current_modwt_result <- current_row$modwt_result
   
   # Initialize a list to store tibbles for the current asv_num
   d_tibbles <- list()
   
   # Loop through modwt_result
   for (j in 1:5) {
     # Extract the asv_num tibble
     d_tibbles[[j]] <- current_modwt_result[[current_asv_row.number]][[j]] %>%
       as_tibble_col(column_name = paste0("d", j))
   }
   
   ##create a column with the asv_num
   current_row <- current_asv_num %>%
     as_tibble() |>
     mutate(asv_name = list(rep(value, each = 120))) |>
     unnest(asv_name) |>
     dplyr::select(-value)
   
   # Combine the tibbles for the current asv_num
   if (length(d_tibbles) > 0) {
     all_tibbles[[current_asv_num]] <- bind_cols(d_tibbles) %>%
       dplyr::mutate(sample_num = row_number()) %>%
       bind_cols(current_row)
   }
 }
 
 # Combine all the tibbles into one
 final_tibble <- bind_rows(all_tibbles)
 
 decimal_date_tibble <-  wavelet_3_df  |>
   arrange(as.numeric(decimal_date)) %$%
   decimal_date |>
   unique() |>
   as_tibble_col(column_name = 'decimal_date') |>
   dplyr::mutate(sample_num = row_number())
 
 tax <- asv_tab_all_bloo_z_tax |>
   dplyr::select(asv_num, phylum, class, order, family, genus) |>
   distinct()
 
 wavelets_result_tibble_tax_3_biased  <-  final_tibble |>
   pivot_longer(cols = !c(asv_name, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
   left_join( decimal_date_tibble) |>
   left_join(tax, by = c('asv_name' = 'asv_num')) |>
   rename(asv_num = asv_name) |>
   dplyr::mutate(wavelets_transformation = str_replace(wavelets_transformation, 'd5', 's4'))
 
 wavelets_result_tibble_tax_3_biased_red <-  wavelets_result_tibble_tax_3_biased |>
   dplyr::mutate(wavelets_result_ed = case_when(wavelets_transformation == 'd1' &
                                                  sample_num %in% c(1, 119, 120) ~ 'NA',
                                                wavelets_transformation == 'd2' &
                                                  sample_num %in% c(1,2,3,  117, 118, 119, 120) ~ 'NA',
                                                wavelets_transformation == 'd3' &
                                                  sample_num %in% c(1,2,3,4,5,6, 113,114,115,116,  117, 118, 119, 120) ~ 'NA',
                                                wavelets_transformation == 'd4' &
                                                  sample_num %in% c(1,2,3,4,5,6,7,8,9,10,11,12, 105,106,107,108,109,110,111,112,113,114,
                                                                    115,116,117, 118, 119, 120) ~ 'NA',
                                                wavelets_transformation == 's4' &
                                                  sample_num %in% c(1,2,3,4,5,6,7,8,9,10,11,12, 
                                                                    13,14,15,16,17,
                                                                    108,109,110,111,112,113,114,
                                                                    115,116,117, 118, 119, 120) ~ 'NA',
                                                TRUE ~ as.character(wavelets_result)))
 
 wavelets_result_ed_tibble_tax_02_biased_red <-  wavelets_result_tibble_tax_02_biased_red |>
   dplyr::mutate(wavelets_result_ed = as.numeric(wavelets_result_ed))
 
 wavelets_result_ed_tibble_tax_3_biased_red <-  wavelets_result_tibble_tax_3_biased_red |>
   dplyr::mutate(wavelets_result_ed = as.numeric(wavelets_result_ed))
 
 ### PLOT THE WAVELETS TRANSFROMATINS COMPUTED (visually inspect them, to be sure of the results)----
 ## FL FRACTION----
 ##i divide them so that i can see the plots better
 bloo_02
 
 wavelets_result_ed_tibble_tax_02_biased_red |>
   dplyr::filter(asv_num %in%  (bloo_02 |>
                                   slice_head(n = 10) %$%
                                   value)) |>
   ggplot(aes(decimal_date, wavelets_result_ed))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_ed_tibble_tax$asv_name), ' ', unique( wavelets_result_ed_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~asv_num)+
   scale_x_continuous(expand = c(0,0))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         #aspect.ratio = 4/10,
         text = element_text(size = 5))
 
 wavelets_result_ed_tibble_tax_02_biased_red |>
   dplyr::filter(asv_num %in%  (bloo_02 |>
                                  slice_tail(n = 11) %$%
                                  value)) |>
   ggplot(aes(decimal_date, wavelets_result_ed))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_ed_tibble_tax$asv_num), ' ', unique(wavelets_result_ed_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~asv_num)+
   scale_x_continuous(expand = c(0,0))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         #aspect.ratio = 4/10,
         text = element_text(size = 5))
 
 ## PA FRACTION---- 
 bloo_3$value #45 in total
 
 wavelets_result_ed_tibble_tax_3_biased_red |>
   dplyr::filter(asv_num %in%  (bloo_3 |>
                                   slice_head(n = 20) %$%
                                   value)) |>
   ggplot(aes(decimal_date, wavelets_result_ed))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_ed_tibble_tax$asv_num), ' ', unique( wavelets_result_ed_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~asv_num)+
   scale_x_continuous(expand = c(0,0))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         #aspect.ratio = 4/10,
         text = element_text(size = 5))
 
 wavelets_result_ed_tibble_tax_3_biased_red |>
   dplyr::filter(asv_num %in%  (bloo_3 |>
                                   slice_tail(n = 20) %$%
                                   value)) |>
   ggplot(aes(decimal_date, wavelets_result_ed))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_ed_tibble_tax$asv_num), ' ', unique(wavelets_result_ed_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~asv_num)+
   scale_x_continuous(expand = c(0,0))+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         #aspect.ratio = 4/10,
         text = element_text(size = 5))
 
 ### create a loop to save all the wavelets transformations computed and visually inspect them
 # Create a folder to save the PDF files
 dir.create("../results/figures/wavelets_plots/biased_red", showWarnings = FALSE)
 
 ### FL----
 # Get unique asv_num values
 asv_nums <- unique(wavelets_result_ed_tibble_tax_02_biased_red$asv_num)
 
 for (asv_num in asv_nums) {
   # Filter data for the current asv_num
   plot_data <- wavelets_result_ed_tibble_tax_02_biased_red %>%
     dplyr::filter(asv_num == !!asv_num)  # Use !! to unquote asv_num
   
   title_plot <- plot_data |>
     group_by(asv_num, family) |>
     dplyr::reframe(title_ed = paste0(asv_num,', ', family, ', ', order)) |>
     distinct() |>
     dplyr::select(title_ed) |>
     as.character()
   
   # Create the plot
   p <- ggplot(plot_data, aes(decimal_date, wavelets_result_ed)) +
     geom_col() +
     labs(x = 'Decimal date', y = 'Wavelets results') +
     facet_grid(vars(wavelets_transformation)) +
     scale_x_continuous(expand = c(0,0)) +
     labs(title = title_plot)+
     theme_bw() +
     theme(panel.grid = element_blank(), strip.background = element_blank(),
           aspect.ratio = 4/10,
           text = element_text(size = 5))
   
   # Save the plot as a PDF file
   pdf_file <- paste0("results/figures/wavelets_plots/biased_red_checked/", asv_num, "_plot_02_biased_red.pdf")
   ggsave(pdf_file, p, width = 6, height = 3, units = "in")
   
   # Print a message indicating the plot has been saved
   cat("Plot for", asv_num, "saved as", pdf_file, "\n")
 }
 
 
 ## PA----
 # Get unique asv_num values
 asv_nums <- unique(wavelets_result_ed_tibble_tax_3_biased_red$asv_num)
 
 for (asv_num in asv_nums) {
   # Filter data for the current asv_num
   plot_data <- wavelets_result_ed_tibble_tax_3_biased_red %>%
     dplyr::filter(asv_num == !!asv_num)  # Use !! to unquote asv_num
   
   title_plot <- plot_data |>
     group_by(asv_num, family) |>
     dplyr::reframe(title_ed = paste0(asv_num,', ', family, ', ', order)) |>
     distinct() |>
     dplyr::select(title_ed) |>
     as.character()
   
   # Create the plot
   p <- ggplot(plot_data, aes(decimal_date, wavelets_result_ed)) +
     geom_col() +
     labs(x = 'Decimal date', y = 'Wavelets results') +
     facet_grid(vars(wavelets_transformation)) +
     scale_x_continuous(expand = c(0,0)) +
     labs(title = title_plot)+
     theme_bw() +
     theme(panel.grid = element_blank(), strip.background = element_blank(),
           aspect.ratio = 4/10,
           text = element_text(size = 5))
   
   # Save the plot as a PDF file
   pdf_file <- paste0("../results/figures/wavelets_plots/biased_red_checked/", asv_num, "_plot_3_biased_red.pdf")
   ggsave(pdf_file, p, width = 6, height = 3, units = "in")
   
   # Print a message indicating the plot has been saved
   cat("Plot for", asv_num, "saved as", pdf_file, "\n")
 }
 
 ### 4. Magnitude (coefficients) observe were they got the highest coefficients for the wavelet analysis--------
 
 #### Wavelet coefficient magnitude indicates how strongly the data are correlated with the mother wavelet at a given frequency and distance.
 
bloo_02_type_biased_red <- wavelets_result_ed_tibble_tax_02_biased_red %>%
   group_by(asv_num, wavelets_transformation) %>%
   dplyr::filter(!is.na(wavelets_result_ed)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) |>  # Calculate the magnitude of coefficients for each level
   group_by(asv_num) |>
   top_n(1, wt = coefficients) |>  # Find the level with the maximum magnitude
   rename(max_coeff = wavelets_transformation) |>
   dplyr::mutate(bloomer_type = case_when(max_coeff == 's4' ~ 'inter-annual',
                                          max_coeff == 'd1' ~ 'fine-scale',
                                          max_coeff == 'd2' ~ 'half-yearly', 
                                          max_coeff == 'd3' ~ 'seasonal', 
                                          max_coeff == 'd4' ~ 'year-to-year')) |>
   left_join(tax, by = 'asv_num') |>
   dplyr::mutate(fraction = '0.2')
 
 bloo_02_type_biased_red |>
   ungroup() |>
   dplyr::filter(bloomer_type == 'seasonal') |>
   dplyr::reframe(n = n())
 
 1/20 ## 15.58% of the free living blooms seem to have a seasonal pattern
 
 bloo_02_type_tax_biased_red <-  bloo_02_type_biased_red |>
   group_by(bloomer_type, family ) |>
   dplyr::reframe(n = n()) |>
   dplyr::mutate(fraction = '0.2')
 
 bloo_3_type_biased_red <- wavelets_result_ed_tibble_tax_3_biased_red %>%
   group_by(asv_num, wavelets_transformation) %>%
   dplyr::filter(!is.na(wavelets_result_ed)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) |>  # Calculate the magnitude of coefficients for each level
   group_by(asv_num) |>
   top_n(1, wt = coefficients) |>  # Find the level with the maximum magnitude
   rename(max_coeff = wavelets_transformation) |>
   dplyr::mutate(bloomer_type = case_when(max_coeff == 's4' ~ 'inter-annual',
                                          max_coeff == 'd1' ~ 'fine-scale',
                                          max_coeff == 'd2' ~ 'half-yearly', 
                                          max_coeff == 'd3' ~ 'seasonal', 
                                          max_coeff == 'd4' ~ 'year-to-year')) |>
   left_join(tax, by = 'asv_num') |>
   dplyr::mutate(fraction = '3')
 
 ##save a dataframe with the type of bloomer according to the wavelets result highest coefficient
 bloo_type_biased_all <- bloo_3_type_biased_red |>
   bind_rows( bloo_02_type_biased_red)
 
 #write.csv(bloo_type_biased_all, '../data/bloo_type_biased_all_checked.csv')
 
 bloo_3_type_biased_red |>
   ungroup() |>
   dplyr::filter(bloomer_type == 'seasonal') |>
   dplyr::reframe(n = n())
 
 6/45 ## 13.33% of the particle attached seem to have a seasonal pattern
 
 bloo_3_type_tax_biased_red <-  bloo_3_type_biased_red |>
   group_by(bloomer_type, family) |>
   dplyr::reframe(n = n()) |>
   dplyr::mutate(fraction = '3')
 
 ## I create a table that summarizes the different type of bloomers that we have, observe their taxonomy too.
 bloo_02_type_tax_biased_red |>
   bind_rows(bloo_3_type_tax_biased_red) |>
   group_by( bloomer_type, fraction) |>
   dplyr::mutate(n_total = sum(n)) |>
   ggplot(aes(bloomer_type, n_total))+
   geom_col(aes(y = n, fill = family))+
   scale_fill_manual(values = palette_family_assigned_bloo)+
   facet_wrap(vars(fraction), labeller = labs_fraction)+
   scale_y_continuous(expand = c(0,0))+
   geom_text(aes(bloomer_type, 1, label = n_total),
             color = "black", nudge_y = 10) +
   theme_bw()+
   #coord_flip()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', axis.ticks.x = element_blank(),
         panel.border = element_blank(), text = element_text(size = 5))
 
 ## in relatives  
 bloo_all_type_tax_biased_red <-  bloo_02_type_tax_biased_red |>
   bind_rows(bloo_3_type_tax_biased_red) |>
   group_by( bloomer_type, fraction) |>
   dplyr::mutate(n_total = sum(n)) |>
   dplyr::group_by(fraction) |>
   dplyr::mutate(n_blooms_fraction = sum(n)) |>
   dplyr::mutate(n_blooms_fraction_perc = n/n_blooms_fraction) |>
   dplyr::mutate(n_blooms_fraction_perc_type = n/n_total) 
 
 bloo_all_type_tax_biased_red$bloomer_type <- factor( bloo_all_type_tax_biased_red$bloomer_type,
                                           levels = c('fine-scale', 'seasonal', 'half-yearly', 'inter-annual'))
 
 ### by family
 library(scales)
 bloo_all_type_tax_biased_red |>
   ggplot(aes(bloomer_type, 1))+
   geom_col(aes(y = n_blooms_fraction_perc_type, fill = family))+
   scale_fill_manual(values = palette_family_assigned_bloo)+
   facet_wrap(vars(fraction), labeller = labs_fraction)+
   scale_y_continuous(expand = c(0,0), labels = percent_format())+
   labs(fill = 'Family', y ='Number of bloomers', x = 'Bloomer period')+
   geom_text(aes(bloomer_type, 1, label = n_total),
             color = "black", nudge_y = -0.3) +
   theme_bw()+
   #coord_flip()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', axis.ticks.x = element_blank(),
         panel.border = element_blank(), text = element_text(size = 5))
 
 ##i create a table with all the coefficients (I do not decide which is the most important)-----
 wavelets_result_ed_tibble_tax_3_biased_red_coeff <- wavelets_result_ed_tibble_tax_3_biased_red %>%
   group_by(asv_num, wavelets_transformation) %>%
   dplyr::filter(!is.na(wavelets_result_ed)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) |>
   dplyr::mutate(fraction = '3') |>
   dplyr::mutate(wavelets_fraction = paste0(wavelets_transformation,'_', fraction)) |>
   dplyr::select(-fraction, -wavelets_transformation)

 wavelets_result_ed_tibble_tax_02_biased_red_coeff <- wavelets_result_ed_tibble_tax_02_biased_red %>%
   group_by(asv_num, wavelets_transformation) %>%
   dplyr::filter(!is.na(wavelets_result_ed)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) |>
   dplyr::mutate(fraction = '0.2') |>
   dplyr::mutate(wavelets_fraction = paste0(wavelets_transformation,'_', fraction)) |>
   dplyr::select(-fraction, -wavelets_transformation)
 
 wavelets_result_ed_tibble_biased_red_coeff_all <-  wavelets_result_ed_tibble_tax_3_biased_red_coeff |>
 bind_rows(wavelets_result_ed_tibble_tax_02_biased_red_coeff) |>
   pivot_wider(id_cols = asv_num, values_from = coefficients, names_from = wavelets_fraction, values_fill = 0)
 
 ### by order ----
 bloo_02_type_tax_order_biased_red <-  bloo_02_type_biased_red |>
   group_by(bloomer_type, order ) |>
   dplyr::reframe(n = n()) |>
   dplyr::mutate(fraction = '0.2')
 
 bloo_3_type_tax_order_biased_red <-  bloo_3_type_biased_red |>
   group_by(bloomer_type, order) |>
   dplyr::reframe(n = n()) |>
   dplyr::mutate(fraction = '3')
 
 bloo_all_type_tax_order_biased_red <-  bloo_02_type_tax_order_biased_red |>
   bind_rows(bloo_3_type_tax_order) |>
   group_by( bloomer_type, fraction) |>
   dplyr::mutate(n_total = sum(n)) |>
   dplyr::group_by(fraction) |>
   dplyr::mutate(n_blooms_fraction = sum(n)) |>
   dplyr::mutate(n_blooms_fraction_perc = n/n_blooms_fraction) |>
   dplyr::mutate(n_blooms_fraction_perc_type = n/n_total) 
 
 bloo_all_type_tax_order_biased_red$bloomer_type <- factor(bloo_all_type_tax_order_biased_red$bloomer_type,
                                                levels = c('fine-scale', 'seasonal', 'half-yearly', 'inter-annual'))
 
 bloo_all_type_tax_order_biased_red |>
   ggplot(aes(bloomer_type, 1))+
   geom_col(aes(y = n_blooms_fraction_perc_type, fill = order))+
   scale_fill_manual(values = palette_order_assigned_bloo)+
   facet_wrap(vars(fraction), labeller = labs_fraction)+
   scale_y_continuous(expand = c(0,0), labels = percent_format())+
   geom_text(aes(bloomer_type, 1, label = n_total),
             color = "black", nudge_y = -0.3) +
   labs(fill = 'Family', y ='Number of bloomers', x = 'Bloomer period')+
   theme_bw()+
   #coord_flip()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', axis.ticks.x = element_blank(),
         panel.border = element_blank(), text = element_text(size = 5))
 
 
 wavelets_result_ed_tibble_tax_3_biased_red <-  wavelets_result_ed_tibble_tax_3_biased_red |>
   dplyr::mutate(fraction = '3')

 wavelets_result_ed_tibble_tax_02_biased_red <-  wavelets_result_ed_tibble_tax_02_biased_red |>
   dplyr::mutate(fraction = '0.2')
 
 #write.csv( wavelets_result_ed_tibble_tax_3_biased_red, '../data/wavelets_analysis/wavelets_result_ed_tibble_tax_3_biased_red.csv')
 #write.csv( wavelets_result_ed_tibble_tax_02_biased_red, '../data/wavelets_analysis/wavelets_result_ed_tibble_tax_02_biased_red.csv')
 
 ## General plots for observing different patterns at different periods of blooming ASVs----
 ## labs fraction and wavelets 
 labs_wavelets_fraction <- as_labeller(c('d1' = 'Fine-scale',
                                'd2' = 'Half-yearly',
                                'd3' = 'Seasonal',
                                'd4' = 'Year-to-year',
                                's4' = 'Inter-annual',
                                '3' = 'Particle attached (3-20 (um)',
                                '0.2' = 'Free living (0.2 - 3 (um)'))
 
 wavelets_result_ed_tibble_tax_02_biased_red_all <-  wavelets_result_ed_tibble_tax_02_biased_red |>
   bind_rows(wavelets_result_ed_tibble_tax_3_biased_red)
 
 wavelets_result_ed_tibble_tax_02_biased_red_all <- wavelets_result_ed_tibble_tax_02_biased_red_all |>
   dplyr::mutate(phylum_f = as_factor(phylum),
                 family_f = as_factor(family),
                 order_f = as_factor(order),
                 class_f = as_factor(class),
                 asv_num_f = as_factor(asv_num))
 
 wavelets_result_ed_tibble_tax_02_biased_red_all$class_f <-  factor(wavelets_result_ed_tibble_tax_02_biased_red_all$class_f, 
                                           levels=unique(wavelets_result_ed_tibble_tax_02_biased_red_all$class_f[order(wavelets_result_ed_tibble_tax_02_biased_red_all$phylum_f)]), 
                                           ordered=TRUE)
 
 wavelets_result_ed_tibble_tax_02_biased_red_all$order_f <-  factor(wavelets_result_ed_tibble_tax_02_biased_red_all$order_f, 
                                           levels=unique(wavelets_result_ed_tibble_tax_02_biased_red_all$order_f[order(wavelets_result_ed_tibble_tax_02_biased_red_all$phylum_f,
                                                                                              wavelets_result_ed_tibble_tax_02_biased_red_all$class_f)]), 
                                           ordered=TRUE)
 
 wavelets_result_ed_tibble_tax_02_biased_red_all$family_f <-  factor(wavelets_result_ed_tibble_tax_02_biased_red_all$family_f, 
                                            levels=unique(wavelets_result_ed_tibble_tax_02_biased_red_all$family_f[order(wavelets_result_ed_tibble_tax_02_biased_red_all$phylum_f,
                                                                                                wavelets_result_ed_tibble_tax_02_biased_red_all$class_f,
                                                                                                wavelets_result_ed_tibble_tax_02_biased_red_all$order_f)]), 
                                            ordered=TRUE)
 
 
 wavelets_result_ed_tibble_tax_02_biased_red_all$asv_num_f <-  factor(wavelets_result_ed_tibble_tax_02_biased_red_all$asv_num_f, 
                                             levels=unique(wavelets_result_ed_tibble_tax_02_biased_red_all$asv_num_f[order(wavelets_result_ed_tibble_tax_02_biased_red_all$phylum_f,
                                                                                                  wavelets_result_ed_tibble_tax_02_biased_red_all$class_f,
                                                                                                  wavelets_result_ed_tibble_tax_02_biased_red_all$order_f,
                                                                                                  wavelets_result_ed_tibble_tax_02_biased_red_all$family_f)]), 
                                             ordered=TRUE)
 
 most_important_transformations <- wavelets_result_ed_tibble_tax_02_biased_red_all |>
   group_by(asv_num, wavelets_transformation, fraction) %>%
   dplyr::filter(!is.na(wavelets_result_ed)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) |>  # Calculate the magnitude of coefficients for each level
   group_by(asv_num, fraction) |>
   top_n(2, wt = coefficients) |>
   dplyr::filter(wavelets_transformation != 'd2') |>
   dplyr::mutate(asv_w_f = paste0(asv_num, wavelets_transformation, fraction))
 
 wavelets_transformation_summary_top2_nod2 <- wavelets_result_ed_tibble_tax_02_biased_red_all |>
   dplyr::mutate(asv_w_f = paste0(asv_num, wavelets_transformation, fraction)) |>
   dplyr::filter(asv_w_f %in%  most_important_transformations$asv_w_f) |>
   ggplot(aes(decimal_date, wavelets_result_ed))+
   labs(x = 'Date', y = 'Wavelet coefficient', color = 'Order')+
   facet_grid(wavelets_transformation~fraction, labeller = labs_wavelets_fraction)+
   geom_line(aes(group = asv_num, color = order))+
   scale_color_manual(values = palette_order_assigned_bloo)+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', text = element_text(size = 8),
         legend.key.width = unit(1, "lines"))+
   guides(color = guide_legend(override.aes = list(size = 4)))
 
 # ggsave(wavelets_transformation_summary_top2_nod2, filename = 'wavelets_transformation_summary_top2_nod2.pdf',
 #        path = 'results/figures/',
 #        width = 188, height = 230, units = 'mm')
 
 wavelets_transformation_summary_top2_nod2 <- wavelets_result_ed_tibble_tax_02_biased_red_all |>
   dplyr::mutate(asv_w_f = paste0(asv_num, wavelets_transformation, fraction)) |>
   dplyr::filter(asv_w_f %in%  most_important_transformations$asv_w_f) |>
   dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
   ggplot(aes(decimal_date, wavelets_result_ed))+
   labs(x = 'Date', y = 'Wavelet coefficient', color = 'Order')+
   facet_grid(wavelets_transformation~fraction, labeller = labs_wavelets_fraction)+
   geom_line(aes(group = asv_num))+
   #scale_color_manual(values = palette_order_assigned_bloo)+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', text = element_text(size = 12),
         legend.key.width = unit(1, "lines"))+
   guides(color = guide_legend(override.aes = list(size = 4)))
 
 # ggsave(wavelets_transformation_summary_top2_nod2, filename = 'wavelets_transformation_summary_top2_nod2_bw_nosar11_cluster.pdf',
        # path = 'results/figures/',
        # width = 188, height = 230, units = 'mm')
        # 
 
abundance_asv7 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'rclr') |>
   #dplyr::filter(asv_w_f %in%  most_important_transformations$asv_w_f) |>
   dplyr::filter(asv_num %in% c('asv7')) |>
   ggplot(aes(decimal_date, abundance_value))+
   labs(x = 'Date', y = 'rCLR', color = 'Order')+
   facet_wrap(vars(fraction), labeller = labs_wavelets_fraction)+
   geom_line(aes(group = asv_num))+
  scale_y_continuous()+
   #scale_color_manual(values = palette_order_assigned_bloo)+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', text = element_text(size = 12),
         legend.key.width = unit(1, "lines"))+
   guides(color = guide_legend(override.aes = list(size = 4)))

ggsave( abundance_asv7, filename = 'abundance_asv7_rclr.pdf',
       path = 'results/figures/',
       width = 188, height = 88, units = 'mm')
 
 wavelets_transformation_summary_top2_nod2_asv7 <- wavelets_result_ed_tibble_tax_02_biased_red_all |>
   dplyr::mutate(asv_w_f = paste0(asv_num, wavelets_transformation, fraction)) |>
   dplyr::filter(!wavelets_transformation %in% c('d2', 'd4')) |>
   #dplyr::filter(asv_w_f %in%  most_important_transformations$asv_w_f) |>
   dplyr::filter(asv_num %in% c('asv7')) |>
   ggplot(aes(decimal_date, wavelets_result_ed))+
   labs(x = 'Date', y = 'Wavelet coefficient', color = 'Order')+
   facet_grid(wavelets_transformation~fraction, labeller = labs_wavelets_fraction)+
   geom_line(aes(group = asv_num))+
   #scale_color_manual(values = palette_order_assigned_bloo)+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', text = element_text(size = 12),
         legend.key.width = unit(1, "lines"))+
   guides(color = guide_legend(override.aes = list(size = 4)))
 
 ggsave(wavelets_transformation_summary_top2_nod2_asv7, filename = 'wavelets_transformation_summary_top2_nod2_bw_asv7.pdf',
        path = 'results/figures/',
        width = 188, height = 230, units = 'mm')

 
### variance of the coefficients at different scales for each ASV-----
 wavelets_result_ed_tibble_tax_02_biased_red <- wavelets_result_ed_tibble_tax_02_biased_red  |>
   dplyr::mutate(fraction = '0.2')
 
 ### variance of the coefficients at different scales for each ASV-----
 labs_wavelets <- as_labeller(c('d1' = 'Fine-scale',
                                'd2' = 'Half-yearly',
                                'd3' = 'Seasonal',
                                'd4' = 'Year-to-year',
                                's4' = 'Inter-annual'))
 
wavelets_variance <-  wavelets_result_ed_tibble_tax_3_biased_red |>
   dplyr::mutate(fraction = '3') |>
   bind_rows( wavelets_result_ed_tibble_tax_02_biased_red) |>
   dplyr::filter(!is.na(wavelets_result_ed)) |>
   dplyr::group_by(asv_num, wavelets_transformation, class, fraction) |>
   dplyr::summarize(variance = var(wavelets_result_ed)) |>
   ungroup() |>
   ggplot(aes(wavelets_transformation, variance))+
   geom_line(aes(group = asv_num, color =  class), linewidth = 0.3 )+
   geom_point(aes(color =  class), size = 0.3)+
   scale_color_manual(values = palette_class_assigned_bloo)+
   geom_violin(alpha = 0.3)+
   facet_wrap(vars(fraction), labeller = labs_fraction, nrow = 2)+
  scale_x_discrete(labels = labs_wavelets)+
   theme_bw()+
   labs(x = 'Wavelet vectors', y = 'Variance', color = 'Class')+
   scale_y_continuous(expand = c(0,0))+
   theme(panel.grid = element_blank(), strip.background =  element_blank(),
         legend.position = 'bottom',
         text = element_text(size = 5),
         legend.key.size = unit(2, 'mm'),
         axis.ticks = element_line(unit(1, 'mm')))
 4459+1363+2027
 wavelets_variance
 
 # ggsave(wavelets_variance, filename = 'wavelets_variance.pdf',
 #        path = 'results/figures/',
 #        width = 88, height = 100, units = 'mm')
 
 #### GENERAL PLOTS FOR THOSE ASVs THAT ARE ACTUALLY SEASONAL FORM THOSE THAT ARE NOT-----
 ##reorder taxonomy as factors ----
 asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
   dplyr::mutate(phylum_f = as.factor(phylum),
                 family_f = as.factor(family),
                 order_f = as.factor(order),
                 class_f = as.factor(class),
                 asv_num_f = as.factor(asv_num))
 
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
 
 
 
 ##EXPLORE FL FRACTION-----
 asv_tab_all_bloo_z_tax |>
   dplyr::filter(fraction == '0.2') |>
   dplyr::filter(asv_num %in%  (bloo_02_type |>
                                  dplyr::filter(bloomer_type == 'seasonal') %$%
                   asv_num)) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
   # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
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
         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  
 
 asv_tab_all_bloo_z_tax |>
   dplyr::filter(fraction == '0.2') |>
   dplyr::filter(asv_num %in%  (bloo_02_type |>
                                  dplyr::filter(bloomer_type == 'random2') %$%
                                  asv_num)) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
   # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
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
         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  

 asv_tab_all_bloo_z_tax |>
   dplyr::filter(fraction == '0.2') |>
   dplyr::filter(asv_num %in%  (bloo_02_type |>
                                  dplyr::filter(bloomer_type == 'random') %$%
                                  asv_name)) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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

 asv_tab_all_bloo_z_tax |>
   dplyr::filter(fraction == '0.2') |>
   dplyr::filter(asv_num %in%  (bloo_02_type |>
                                  dplyr::filter(bloomer_type == 'residual') %$%
                                  asv_name)) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
 
 ## EXPLORE PA RESULTS ----
 asv_tab_all_bloo_z_tax |>
   dplyr::filter(fraction == '3') |>
   dplyr::filter(asv_num %in%  (bloo_3_type |>
                                  dplyr::filter(bloomer_type == 'seasonal') %$%
                                  asv_name)) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
   # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
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
         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  
 
 asv_tab_all_bloo_z_tax |>
   dplyr::filter(fraction == '3') |>
   dplyr::filter(asv_num %in%  (bloo_3_type |>
                                  dplyr::filter(bloomer_type == 'random2') %$%
                                  asv_name)) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
   # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
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
         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  
 
 asv_tab_all_bloo_z_tax |>
   dplyr::filter(fraction == '3') |>
   dplyr::filter(asv_num %in%  (bloo_3_type |>
                                  dplyr::filter(bloomer_type == 'random') %$%
                                  asv_name)) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
   # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
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
         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside') 
 
 asv_tab_all_bloo_z_tax |>
   dplyr::filter(fraction == '3') |>
   dplyr::filter(asv_num %in%  (bloo_3_type |>
                                  dplyr::filter(bloomer_type == 'residual') %$%
                                  asv_name)) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
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
   # geom_line(data = community_eveness_all_m, aes(date, community_eveness_rar/1.6), color = '#2D2A2B', alpha = 0.8)+
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
         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  
 
 
 ## DESIGN OF SEASONAL PLOTS -----
 
 asv_tab_all_bloo_z_tax |>
   dplyr::filter(asv_num %in%  (bloo_3_type |>
                                  dplyr::filter(bloomer_type == 'seasonal') %$%
                                  asv_name) & fraction == '3' |
                   asv_num %in%  (bloo_02_type |>
                                    dplyr::filter(bloomer_type == 'seasonal') %$%
                                    asv_name) & fraction == '0.2' ) |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
   dplyr::filter(abundance_type == 'zclr') |>
   ggplot(aes(day_of_year, abundance_value, color = as.factor(year), group = fraction))+
   geom_point(aes(group = year, shape = fraction), alpha = 0.6)+
   scale_color_manual(values = palette_years, na.value = "#000000")+
   labs(x = 'Day of the year', y = 'rCLR', color = 'Year', shape = 'Fraction')+
   facet_wrap(vars(asv_num_f))+
   geom_boxplot(aes(group = season), alpha = 0.1)+
   geom_smooth(method = 'loess', aes(group = asv_num), color = 'black')+ ##afegir grup per fraction asv27 té les dues
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

 
 ## Apply Fourier analysis to compare with wavelets results -----
 
 par(fig = c(0, 0.48, 0.5, 1), mar = c(3, 5, 0, 0)) 
 plot((abund11), bty = "n", cex.axis = 1.3, ylab = "", xlab = "", axes = F) 
 axis(1, cex.axis = 1.3) 
 axis(2, cex.axis = 1.3, las = 2) 
 box(bty = "L", lwd = 1.1) 
 mtext(side = 2, line = 3.5, outer = F, text = expression("[N" * O[3]^"-" * "] (" * + mu * "mol/l)"), cex = 1.3) 
 par(new = T, fig = c(0, 0.24, 0, 0.5), mar = c(4, 5, 0, 0)) 
 spec.pgram(abund11, detrend = T, spans = c(5, 5), bty = "n", main = "", sub = "", 
            xlab = "", ylab = "", demean = T, axes = F) 
 axis(1, cex.axis = 1.3) 
 axis(2, cex.axis = 1.3, las = 2) 
 box(bty = "L", lwd = 1.1) 
 mtext(side = 2, line = 3.5, outer = F, at = 0.15, text = "Periodogram", cex = 1.3) 
 mtext(side = 1, line = 2, at = 3, outer = F, text = expression("Freq. (year"^"-1" * + ")"), cex = 0.5) 
 abline(v = 1, col = "red") 
 par(new = T, fig = c(0.24, 0.48, 0, 0.5), mar = c(4, 5, 0, 0)) 
 acf(abund11, demean = T, ci.type = "white", bty = "n", main = "", sub = "", 
     xlab = "", ylab = "", axes = F) 
 axis(1, cex.axis = 1.3) 
 axis(2, cex.axis = 1.3, las = 2) 
 box(bty = "L", lwd = 1.1) 
 mtext(side = 2, line = 3.5, outer = F, at = 0.25, text = "Autocorrelation", cex = 1.3) 
 mtext(side = 1, line = 2, at = 0.75, outer = F, text = "Time lag (year)", cex = 1.3) 
 
 Spline <- smooth.spline(seq(1995, 2004 + 11/12, by = 1/12), abund178, spar = 0.75)$y 
 par(new = T, fig = c(0.52, 1, 0.5, 1), mar = c(3, 5, 0, 0)) 
 plot(abund178, bty = "n", cex.lab = 1, cex.axis = 1, ylab = "", xlab = "", axes = F) 
 axis(1, cex.axis = 1.3) 
 axis(2, cex.axis = 1.3, las = 2) > box(bty = "L", lwd = 1.1) 
 mtext(side = 2, line = 3.5, outer = F, text = expression("[DON] (" * mu * "mol/l)"),  cex = 1.3) 
 lines(seq(1995, 2004 + 11/12, by = 1/12), Spline, col = "red")
 par(new = T, fig = c(0.52, 0.76, 0, 0.5), mar = c(4, 5, 0, 0)) 
 spec.pgram(abund178, detrend = T, spans = c(5, 5), bty = "n", main = "", sub = "", 
            xlab = "", ylab = "", demean = T, axes = F) 
 axis(1, cex.axis = 1.3) 
 axis(2, cex.axis = 1.3, las = 2) 
 box(bty = "L", lwd = 1.1) 
 mtext(side = 2, line = 3.5, outer = F, at = 0.055, text = "Periodogram", cex = 1.3) 
 mtext(side = 1, line = 2, at = 3, outer = F, text = expression("Freq. (year"^"-1" * + ")"), cex = 1.3) 
 abline(v = 1, col = "red")
 
 par(new = T, fig = c(0.76, 1, 0, 0.5), mar = c(4, 5, 0, 0)) 
 acf(abund178 - Spline, demean = T, ci.type = "white", bty = "n", main = "", sub = "", 
     xlab = "", ylab = "", axes = F) 
 axis(1, cex.axis = 1.3) 
 axis(2, cex.axis = 1.3, las = 2) 
 box(bty = "L", lwd = 1.1) 
 mtext(side = 2, line = 3.5, outer = F, at = 0.4, text = "Autocorrelation", cex = 1.3) 
 mtext(side = 1, line = 2, at = 0.75, outer = F, text = "Time lag (year)", cex = 1.3)
 
 ## code it in a more friendly way
 library(ggplot2)
 library(tidyr)
 
 ##trying to solve the error 
 # Compute the periodogram
 periodogram <- spec.pgram(abund178, detrend = TRUE, spans = c(5, 5),
                           taper = 0.1, pad = 1)
 
 # Extract frequencies and corresponding spectral density values
 frequencies <- periodogram$freq
 spec_density <- periodogram$spec
 
 # Create a data frame for plotting
 periodogram_df <- data.frame(frequencies, spec_density)
 
 sampling_freq <- asv178$decimal_date |>
   as_tibble_col(column_name = "real_sampling") 
 
 perfect_sampling_freq <- seq(2004, 2013 + 11/12, by = 1/12)|>
   as_tibble_col(column_name = 'ideal_sampling') 
 
 
 # Data preparation
 asv11$decimal_date
 abund178_df <- data.frame(year = asv178$decimal_date, #seq(2004, 2013 + 11/12, by = 1/12),
                           abund178 = abund178,
                           Spline = Spline,
                           abund178_minus_Spline = abund178 - Spline)
 
 abund11_df <- data.frame(year = asv11$decimal_date, #seq(2004, 2013 + 11/12, by = 1/12),
                          abund11 = abund11,
                          Spline = Spline,
                          abund11_minus_Spline = abund11 - Spline)
 
 acf_df <- acf(abund11_df$abund11_minus_Spline, lag.max = 120, demean = TRUE)$acf |>
   as_tibble() |>
   dplyr::mutate(n = row_number())
 
 acf_df <- acf(abund178_df$abund178_minus_Spline, lag.max = 120, demean = TRUE)$acf |>
   as_tibble() |>
   dplyr::mutate(n = row_number())
 
 # Plotting using ggplot2
 plots <- list(
   ggplot(abund178_df, aes(x = year, y = abund178)) +
     geom_line(color = "black") +
     theme_minimal() +
     labs(x = NULL, y = expression("[N" * O[3]^"-" * "] (" * mu * "mol/l)")),
   
   # Plot the periodogram
   ggplot(periodogram_df, aes(x = frequencies, y = spec_density)) +
     geom_line(color = "black") +
     theme_minimal() +
     labs(x = "Frequency", y = "Spectral Density", title = "Spectral Analysis") +
     geom_vline(xintercept = 1, linetype = "dashed", color = "red"),
   
   ggplot(acf_df, aes(n, y = V1)) +
     geom_line(color = "black") +
     geom_col()+
     theme_minimal() +
     labs(x = "Time lag (year)", y = "Autocorrelation", title = "Autocorrelation Analysis") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "red")
 )
 
 # Arrange and display the plots
 do.call(gridExtra::grid.arrange, c(plots, ncol = 3))
 
 
 
 
 
 
 ## Strong s4 signal 
 most_important_transformations <- wavelets_result_ed_tibble_tax_02_biased_red_all |>
   group_by(asv_num, wavelets_transformation, fraction) %>%
   dplyr::filter(!is.na(wavelets_result_ed)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) |>  # Calculate the magnitude of coefficients for each level
   group_by(asv_num, fraction) |>
   top_n(2, wt = coefficients) |>
   dplyr::mutate(asv_w_f = paste0(asv_num, wavelets_transformation, fraction))
 
 wavelets_transformation_summary_top3_s4 <- wavelets_result_ed_tibble_tax_02_biased_red_all |>
   dplyr::filter(wavelets_transformation == 's4') |>
   left_join(bloo_all_types_summary, by = c('asv_num', 'fraction')) |>
   dplyr::mutate(asv_w_f = paste0(asv_num, wavelets_transformation, fraction)) |>
   dplyr::filter(asv_w_f %in%  most_important_transformations$asv_w_f) |>
   ggplot(aes(decimal_date, wavelets_result_ed))+
   labs(x = 'Date', y = 'Wavelet coefficient', color = 'Family')+
   facet_wrap(fraction~clustering_group)+ #, labeller = labs_clusters_pa_fl
   scale_x_continuous(expand = c(0,0))+
   #geom_line(aes(group = asv_num, color = order))+
   geom_smooth(aes(group = asv_num, color = family_f), span = 1)+
   scale_color_manual(values = palette_family_assigned_bloo)+
   theme_bw()+
   theme(panel.grid = element_blank(), strip.background = element_blank(),
         legend.position = 'bottom', text = element_text(size = 8),
         legend.key.width = unit(1, "lines"))+
   guides(color = guide_legend(override.aes = list(size = 4)))
 
 # ggsave( wavelets_transformation_summary_top3_s4, filename = 'wavelets_transformation_summary_top2_s4.pdf',
 #        path = 'results/figures/',
 #        width = 188, height = 150, units = 'mm')
 
 
 # APPLY CONTINUOUS WAVELET TRANSFORMATION #### ------
 library(WaveletComp)
 
 ## we need dates in date format not decimal date 
 m_02_dates <- m_02 |>
   dplyr::select(date, decimal_date) |>
   dplyr::mutate(date = as.POSIXct(date, "%Y-%M-%D"))

 wavelet_02_df_date <- wavelet_02_df |>
   left_join( m_02_dates) |>
   dplyr::select(date = date, abundance_value, asv_num )
 
 wavelet_3_df_date <- wavelet_3_df |>
   dplyr::mutate(decimal_date = as.numeric(decimal_date)) |> 
   left_join( m_02_dates) |>
   dplyr::select(date = date, abundance_value = rclr, asv_num ) |>
   dplyr::mutate(abundance_value = as.numeric(abundance_value))
 
 # asv38_test <- wavelet_02_df |>
 #   left_join( m_02_dates) |>
 #   dplyr::filter(asv_num == 'asv38') |>
 #   dplyr::select(-asv_num) |>
 #   dplyr::select(date = date, abundance_value )
 # 
 # asv11_test <- wavelet_02_df |>
 #   left_join( m_02_dates) |>
 #   dplyr::filter(asv_num == 'asv11') |>
 #   dplyr::select(-asv_num) |>
 #   dplyr::select(date = date, abundance_value )
 # 
 # asv27_test <- wavelet_02_df |>
 #   left_join( m_02_dates) |>
 #   dplyr::filter(asv_num == 'asv27') |>
 #   dplyr::select(-asv_num) |>
 #   dplyr::select(date = date, abundance_value )

#  # Compute the continuous wavelet transform
#  cwt_chirp <- analyze.wavelet(asv38_test, 
#                               my.series = 2, loess.span = 0, 
#                               dt = 1, dj = 1/20, 
#                               lowerPeriod = 2, 
#                               upperPeriod = 32, 
#                               make.pval = TRUE, method = "white.noise", params = NULL,
#                               n.sim = 100, 
#                               date.format = '%Y-%M-%d', date.tz = NULL, 
#                               verbose = TRUE)
#  
#  wt.image( cwt_chirp, color.key = "interval", main = "wavelet power spectrum",
#           legend.params = list(lab = "wavelet power levels"),
#           periodlab = "period (months)")
#  
#  wt.avg(cwt_chirp, siglvl = 0.05, sigcol = "red", 
#         periodlab = "period (months)")
#  
#  cwt_chirp <- analyze.wavelet(asv11_test, 
#                               my.series = 2, loess.span = 0, 
#                               dt = 1, dj = 1/20, 
#                               #lowerPeriod = 2*dt, 
#                               #upperPeriod = floor(nrow(my.data)/3)*dt, 
#                               make.pval = TRUE, method = "white.noise", params = NULL,
#                               n.sim = 100, 
#                               date.format = '%Y-%M-%d', date.tz = NULL, 
#                               verbose = TRUE)
#  
#  wt.image( cwt_chirp, color.key = "interval", main = "wavelet power spectrum",
#            legend.params = list(lab = "wavelet power levels"),
#            periodlab = "period (months)")
#  
#  wt.avg(cwt_chirp, siglvl = 0.05, sigcol = "red", 
#         periodlab = "period (months)")
#  
#  cwt_chirp <- analyze.wavelet(asv27_test, 
#                               my.series = 2, loess.span = 0, 
#                               dt = 1, dj = 1/20, 
#                               lowerPeriod = 2, 
#                               upperPeriod = 32, 
#                               make.pval = TRUE, method = "white.noise", params = NULL,
#                               n.sim = 100, 
#                               date.format = '%Y-%M-%d', date.tz = NULL, 
#                               verbose = TRUE)
#  
# wt.image( cwt_chirp, color.key = "interval", main = "wavelet power spectrum",
#            legend.params = list(lab = "wavelet power levels"),
#            periodlab = "period (months)")
#  
# wt.avg(cwt_chirp, siglvl = 0.05, sigcol = "red", 
#         periodlab = "period (months)")
 
 
 
 #### do it at the same time for all my asvs 
# Define a function to perform wavelet analysis and plot the results
 

wavelet_analysis <- function(asv_num, data, fraction) {
  # Filter the data for the current ASV
  asv_data <- data |>
    dplyr::filter(asv_num == !!asv_num) |>
    dplyr::select(date, abundance_value)
  
  # Compute the continuous wavelet transform
  cwt_chirp <- analyze.wavelet(asv_data, ## this function uses a Morlet wavelet.
                               my.series = 2, loess.span = 0, 
                               dt = 1, # number of observations per time unit
                               dj = 1/12, 
                               lowerPeriod = 2, 
                               upperPeriod = 32, 
                               make.pval = TRUE, method = "white.noise", params = NULL,
                               n.sim = 100, 
                               date.format = '%Y-%M-%d', date.tz = NULL, 
                               verbose = TRUE)
  
  # Plot the wavelet power spectrum
  wt.image(cwt_chirp, color.key = "interval", 
           main = paste0("Wavelet Power Spectrum - ASV", asv_num), 
           legend.params = list(lab = "wavelet power levels"),
           periodlab = "period (months)")
  
  # Plot the wavelet average
  wt.avg(cwt_chirp, siglvl = 0.05, sigcol = "red", 
         periodlab = "period (months)")
  
  # Create the directory if it doesn't exist
  dir.create(paste0("results/figures/wavelets_continuous_", fraction, "_ed/"), recursive = TRUE)
  
  # Save the wavelet power spectrum plot
  pdf(file = paste0("results/figures/wavelets_continuous_", fraction, "_ed/", asv_num, "_wavelet_power_spectrum.pdf"), 
      width = 40, height = 34
      # , units = "in"
  )
  wt.image(cwt_chirp, color.key = "interval", 
           main = paste0("Wavelet Power Spectrum - ASV", asv_num), 
           legend.params = list(lab = "wavelet power levels"),
           periodlab = "period (months)")
  dev.off()
  
  # Save the wavelet average plot
  pdf(file = paste0("results/figures/wavelets_continuous_", fraction, "_ed/", asv_num, "_wavelet_average.pdf"), 
      width = 40, height = 34#, units = "in"
  )
  wt.avg(cwt_chirp, siglvl = 0.05, sigcol = "red", 
         periodlab = "period (months)")
  dev.off()
}

# Apply the function to all unique ASV numbers
asv_nums <- unique(wavelet_02_df_date$asv_num)
for (asv_num in asv_nums) {
  wavelet_analysis(asv_num, wavelet_02_df_date, fraction = 0.2)
}
 
asv_nums <- unique(wavelet_3_df_date$asv_num)
for (asv_num in asv_nums) {
  wavelet_analysis(asv_num, wavelet_3_df_date, fraction = 3)
}

## apply the function to all env data ----
env_data_interpolated_values_all_z_score <- read.csv2('data/env_data/env_data_interpolated_values_all_z_score.csv', sep = ';')

env_data_interpolated_values_all_z_score  <- env_data_interpolated_values_all_z_score  |>
  dplyr::select(-'X')

env_data_interpolated_values_all <- env_data_interpolated_values_all_z_score |>
left_join(m_02_dates) |>
  dplyr::mutate(date = as.POSIXct(date, "%Y-%M-%D")) |>
  dplyr::select(date, abundance_value = z_score_environmental_variable, asv_num = environmental_variable) |>
  dplyr::mutate(abundance_value = as.numeric(abundance_value))


env_data_interpolated_values_all 

asv_nums <- unique(env_data_interpolated_values_all$asv_num)
for (asv_num in asv_nums) {
  wavelet_analysis(asv_num, env_data_interpolated_values_all)
}


 cwt(chirp, sj = 1:128, dj = 1/4, mother = "morlet") 
 
 ## Barplot with seasonality and the taxonomy ----
 
 ### recurrency
 labs_fraction_rec_freq <-  as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                                          '3' = 'Particle attached (3-20 um)',
                                          no = 'Recurrent',
                                          yes = 'Non-recurrent',
                                          seasonal = 'Seasonal',
                                          stochastic = 'Chaotic'))
 summary_types_of_blooms |>
   colnames()

 bloo_all_types_summary_tb_tax

 tax_seasonality_plot <- bloo_all_types_summary_tb_tax |>
   group_by(frequency, fraction, order) |>
   dplyr::reframe(n = n()) |>
   group_by(frequency, fraction) |>
   dplyr::mutate(n_total = sum(n)) |>
   dplyr::mutate(perc = n/n_total) |>
   ggplot(aes(perc,as.factor(fraction), fill = order))+
   scale_y_discrete(labels = labs_fraction)+
   scale_x_continuous(labels = percent_format())+
   geom_col()+
   labs(x='', y = '')+
   scale_fill_manual(values = palette_order_assigned_bloo)+
   facet_wrap(vars(frequency), labeller = labs_fraction_rec_freq)+
   theme_bw()+
   theme(legend.position = 'none',
         panel.grid = element_blank(),
         panel.border = element_blank(),
         strip.background = element_blank(),
         text = element_text(size = 14))

 # ggsave('tax_seasonality.svg',  tax_seasonality_plot,
 #        path = "results/figures/poster_svg_format/",
 #        width = 220,
 #        height = 120,
 #        units = 'mm')

 
 