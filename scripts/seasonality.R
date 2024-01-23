# packages
library(fdrtool)
library(GeneCycle)
library(patchwork)
library(lubridate)
library(ggthemes)
library(stringr)
library(tidyverse)
library(phyloseq)
library(ggpubr)

# Differentiate seasonal from non-seasonal bloomers to analyze them separately----

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

## prepare the dataset previous to computing the interpolation  
asv_tab_all_bloo_3_interpol <- asv_tab_all_bloo_z_tax_3 |>
  dplyr::select(date, abundance_value, asv_num, abundance_type) |>
  dplyr::filter(#asv_num == 'asv179' &
                  abundance_type == 'relative_abundance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

# Convert the target date to a Date object
target_date1 <- as.POSIXct("2004-03-22", format = "%Y-%m-%d")
target_date2 <- as.POSIXct("2005-02-15", format = "%Y-%m-%d")
target_date3 <- as.POSIXct("2005-05-10", format = "%Y-%m-%d")

# Perform linear interpolation
asv_tab_all_bloo_3_interpol <- asv_tab_all_bloo_3_interpol |>
  group_by(asv_num) |>
  dplyr::reframe(target_date1 = approx(date, abundance_value, xout = target_date1)$y, #The approx function fits a linear spline to the data and estimates values at the specified target dates.
                 target_date2 = approx(date, abundance_value, xout = target_date2)$y,
                 target_date3 = approx(date, abundance_value, xout = target_date3)$y) |>
  ungroup() |>
  pivot_longer(cols = starts_with('target'), values_to = 'abundance_value', names_to = 'date') |> #prepare the dataset to add it to the general dataset
  dplyr::mutate(date = case_when(date == 'target_date1' ~ '2004-03-22',
                                 date == 'target_date2' ~ '2005-02-15',
                                 date == 'target_date3' ~ '2005-05-10')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

##check how do the interpolations look like are they normal?----
 asv_tab_all_bloo_z_tax_3 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(date, abundance_value, asv_num, abundance_type) |>
  dplyr::filter(#asv_num == 'asv179' &
    abundance_type == 'relative_abundance') |>
  dplyr::select(-abundance_type) |>
  bind_rows(asv_tab_all_bloo_3_interpol) |>
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

## I add the metadata and taxonmical data to them, and then I add it to the PA general dataset
asv_tab_all_bloo_z_tax_02_subset <- asv_tab_all_bloo_z_tax_02 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(date %in% c(as.POSIXct('2004-03-22', format = "%Y-%m-%d"),
                            as.POSIXct('2005-02-15', format = "%Y-%m-%d"),
                            as.POSIXct('2005-05-10', format = "%Y-%m-%d"))) |>
  dplyr::filter(#asv_num == 'asv179' &
    abundance_type == 'relative_abundance') |>
  dplyr::select(-seq, -asv_num, -family, - phylum, -class, -order, -genus, -abundance_type, -abundance_value, -z_score_ra, -sample_id_num, -domain, -sample_id, -reads, -X) |>
  dplyr::distinct()

asv_tab_all_bloo_z_tax_3 <- asv_tab_all_bloo_z_tax_3 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 

asv_tab_all_bloo_3_complete <- asv_tab_all_bloo_3_interpol |>
  right_join(asv_tab_all_bloo_z_tax_02_subset, by = ('date')) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  bind_rows(asv_tab_all_bloo_z_tax_3)


asv_tab_all_bloo_3_complete %$%
  date |>
  unique() #120 we have all the data!!!


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

##Wavelet analysis is also well suited when signals exhibit sharp spikes or local
## discontinuities, features that are very poorly represented by Fourier analysis
## file:///Users/onadeulofeucapo/Downloads/engeland.pdf

## prior to the analysis we need to use CLR transformed data------


##(need to be re-calcualted for the PA fraction!!! )

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
bloo_02 <- read.csv('data/bloo_02.csv') |>
  as_tibble()

bloo_3 <-read.csv('data/bloo_3.csv') |>
  as_tibble()

bloo_02$value

# [1] "asv38"  "asv8"   "asv5"   "asv3"   "asv2"   "asv15"  "asv18"  "asv27"  "asv17"  "asv555"
# [11] "asv237" "asv563" "asv62"  "asv58"  "asv1"   "asv7"   "asv178" "asv282" "asv11" 

bloo_3$value

# [1] "asv179" "asv225" "asv264" "asv200" "asv15"  "asv385" "asv72"  "asv471" "asv27"  "asv153"
# [11] "asv17"  "asv77"  "asv43"  "asv192" "asv84"  "asv118" "asv311" "asv23"  "asv85"  "asv25" 
# [21] "asv163" "asv80"  "asv69"  "asv219" "asv116" "asv182" "asv126" "asv223" "asv105" "asv28" 
# [31] "asv1"   "asv7"   "asv4"   "asv31"  "asv276" "asv22"  "asv194" "asv559" "asv11"  "asv100"
# [41] "asv302" "asv113" "asv511" "asv42"  "asv49" 

#Establish the numbers as a time series
four |>
  colnames()

asv1.ts <- four %>%
  filter(asv_num == "asv8") %>%
  ts_conversor()

## I try to apply this part of the code to all potential bloomers ASVs
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
##I will try the biwavelet package
#install.packages('biwavelet')
library(biwavelet)

##we run the analysis for just one ASV
asv11 <- asv_tab_all_bloo_z_tax_02 |> #FL dataset 
  dplyr::filter(asv_num == 'asv11') |>
  dplyr::filter(abundance_type == 'relative_abundance') |> #see if we can perform the same analysis with zclr transformations
  dplyr::select(decimal_date, abundance_value)

asv_tab_all_bloo_z_tax_02 |>
  colnames()

asv11 |>
  colnames()

# Assuming 'your_data' is your microbiome time series data frame with 'date' and 'abundance_value'
# Convert 'date' to numeric if not already in a suitable format
asv11$date_numeric <- as.numeric(asv11$decimal_date)

# Perform wavelet analysis
wt_result <- wt(asv11)

# Plot the wavelet power spectrum
plot(wt_result)

# Assuming 'asv11' is your microbiome time series data frame with 'date' and 'abundance_value'
# Convert 'date' to Date class
asv11$date <- as.Date(asv11$date)

# Create a regular time series with constant time increment
regular_time_series <- zoo::zoo(asv11$abundance_value, asv11$date)
regular_time_series <- zoo::na.approx(regular_time_series)

# Perform wavelet analysis
wt_result <- wt(regular_time_series) #continuous wavelet transform

# Plot the wavelet power spectrum
plot(wt_result)


## I try the discrete wavelet transform-----

## need to pick a family of wavelets we want to work with.
### using the waveslim package 
library(waveslim)
asv11 <- asv_tab_all_bloo_z_tax_02 |> #FL dataset 
  dplyr::filter(asv_num == 'asv11') |>
  dplyr::filter(abundance_type == 'zclr') |> #see if we can perform the same analysis with zclr transformations
  dplyr::select(decimal_date, abundance_value)

asv178 <- asv_tab_all_bloo_z_tax_02 |> #FL dataset 
  dplyr::filter(asv_num == 'asv178') |>
  dplyr::filter(abundance_type == 'zclr') |> #see if we can perform the same analysis with zclr transformations
  dplyr::select(decimal_date, abundance_value)

## observe my data
asv11 |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line()+
  theme_bw()

asv178 |>
  ggplot(aes(decimal_date, abundance_value))+
  geom_line()+
  theme_bw()

# Extract the abundance values
abundance_values <- asv11$abundance_value
abundance_values <- asv178$abundance_value

ibmr.la8 <- modwt(abundance_values, n.levels = 4, wf = "la8")
names(ibmr.la8) <- c("w1", "w2", "w3", "w4", "v4")

# Access the wavelet coefficients for the first scale
coefficients_scale1 <- ibmr.la8$d4

# Plot the wavelet coefficients for the first scale
plot(coefficients_scale1, type = "l", col = "blue", xlab = "Decimal date", ylab = "Wavelet Coefficients", main = "Wavelet Coefficients - Scale 1")

# Reconstruct the signal from the wavelet coefficients and scaling coefficients
reconstructed_signal <- imodwt(modwt_result)

# Compare the original signal with the reconstructed signal
plot(asv11$decimal_date, asv11$abundance_value, type = "l", col = "blue", lty = 1, xlab = "Decimal Date", ylab = "Abundance Value", main = "Original vs Reconstructed Signal")
lines(asv11$decimal_date, reconstructed_signal, col = "red", lty = 6)
legend("topright", legend = c("Original Signal", "Reconstructed Signal"), col = c("blue", "red"), lty = 1:2)

## try help of the package waveslim
ibm.returns <- diff(abundance_values)

# Perform Maximal Overlap Discrete Wavelet Transform (MODWT)
ibmr.la8 <- modwt(ibm.returns, n.levels = 4, wf = "la8")
names(ibmr.la8) <- c("w1", "w2", "w3", "w4", "v4")

ibmr.la8 <- phase.shift(ibmr.la8, "la8")
## plot partial MODWT for IBM data
par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
plot.ts(ibm.returns, axes=FALSE, ylab="", main="(a)")
for(i in 1:5)
  plot.ts(ibmr.haar[[i]], axes=FALSE, ylab=names(ibmr.haar)[i])
axis(side=1, at=seq(0,368,by=23), 
     labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
plot.ts(ibm.returns, axes=FALSE, ylab="", main="(b)")
for(i in 1:5)
  plot.ts(ibmr.la8[[i]], axes=FALSE, ylab=names(ibmr.la8)[i])
axis(side=1, at=seq(0,368,by=23), 
     labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))

# Here I try to apply the code that is in apendix A of the thesis -----

### Wavelets assume constant a constant sampling interval
### How constant is my timeseries data?
## MODWT (Maximal overlap discrete wavelet transforms):

# abundance vector. 
abund11 <- asv11$abundance_value
modwtasv11 <- modwt(abund11, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
  phase.shift(wf = 'la8') # this correction is necessary because the filtering with wavelet filter induces a displacement of features relative to the original time series.
  
modwtasv11_biased <- modwt(abund11, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  #brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
  phase.shift(wf = 'la8') # this correction is necessary because the filtering with wavelet filter induces a displacement of features relative to the original time series.

abund178 <- asv178$abundance_value
modwtasv178 <- modwt(abund178, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
  phase.shift(wf = 'la8') 

modwtasv178_biased <- modwt(abund178, wf = 'la8', boundary = "periodic", n.levels = 4) |>
  #brick.wall(wf = 'la8') |> #elimination of all boundary coefficients is accomplished by the function 'brick.wall' prior to the phase shift correction
  phase.shift(wf = 'la8') # this correction is necessary because the filtering with wavelet filter induces a displacement of features relative to the original time series.

# Continue with phase shift
shifted_result <- phase.shift(brick.wall(result, filter = db5_filter), filter = db5_filter)

# e-folding estimates
x <- rep(0, 10001) 
x[5000] <- 1 
n.levels <- 4 
len <- length(abund) 
temp <- phase.shift(modwt(x, n.levels = n.levels, wf = "la8"), wf = "la8")

 waveExtremes <- matrix(nrow = 3, ncol = n.levels + 1) 
 colnames(waveExtremes) <- c(paste("d", 1:n.levels, sep = ""), paste("s", n.levels,  sep = "")) 
 rownames(waveExtremes) <- c("left", "right", "top")

  for (i in 1:(n.levels + 1)) waveExtremes[, i] <- c(range(which(abs(temp[[i]]) 
                                                                 >= max(abs(temp[[i]]))/(exp(2)))), which.max(abs(temp[[i]])))

 boundaries <- data.frame(end = len - (waveExtremes[3, ] - waveExtremes[1, ]), 
                          start = waveExtremes[2, ] - waveExtremes[3, ])
 
 for (j in 1:(n.levels + 1)) { 
   is.na(modwtasv11[[i]]) <- c(1:boundaries$start[i], boundaries$end[i]:length(modwtasv11[[i]])) 
   is.na(modwtasv178[[i]]) <- c(1:boundaries$start[i], boundaries$end[i]:length(modwtasv178[[i]])) 
   }
 
 ##PLOTS----
 #### Interpretation of the following plots
 #### MODWTs of CLR transformation of the relative abundance from the potential blooming ASVs. The series d1-d4 are wavelet coefficient
 ### vectors; the series s4 is the scaling coefficient vector at the coarsest scale. Coefficients not enclosed by the red solid lines are
 ### boundary coefficients (coefficients affected by circular filtering). The blue dashed lines are delineate the filter-specific e-folding
 ### distances from the beginning and end of the series (cone-of-influence in Torrence and Compo, 1998).
 
 ### d1 is the finest scale and d4 is the coarsest. The s4 depicts the scaling coefficients as a function of time, representing the remainder
 ### of the variability after the finder sacles (d1-d4) have been isolated. 
 
 ### The inner red, solid, vertical lines represent the limits of the boundary coefficients. The regions enclosed by either set of limits become
 ### progressively smaller as wavelet filters become wider. However, by allowing some bias (i.e as given by the region within the blue lines in each
 ### graph) event the coarser scales can be investigated over a relatively larger temporal scale.
 
 ### To observe which component dominates we need to observe the one that has the highest coefficients in absolute value.
 
 
### MODWT of the ASV11
 # Create an empty list to store your plots
 plotList <- vector("list", n.levels)
 
 for (i in 1:n.levels) { 
   plot(as.vector(time(abund11)), modwtasv11_biased[[i]], xlab = "", ylab = "", type = "h", 
        bty = "n", axes = F, ylim = range(modwtasv11[1:n.levels], na.rm = TRUE)) 
   axis(2, cex.axis = 1.5, las = 2) +
     abline(v = seq(start(abund11)[1], end(abund11)[1], by = 1), col = "gray", 
            lty = "dashed") + abline(h = 0) 
   abline(v = time(abund11)[range(which(!is.na(modwtasv11[[i]])))], lwd = 2, 
          col = "red") + 
     abline(v = time(abund11)[c(boundaries$start[i], boundaries$end[i])], col = "blue", 
            lwd = 3, lty = "dashed") +
     mtext(side = 4, line = -2, cex = 1.2, outer = FALSE, at = par("usr")[4], text = paste("d", i, sep = ""), las = 2) 
   
   # Save each plot to the list
   plotList[[i]] <- recordPlot()
 }
 
 # Set the title for the entire series of plots
 title(main = "ASV 11")
 
 # Replay and display the saved plots
 for (i in 1:n.levels) {
   replayPlot(plotList[[i]])
 }

plot(as.vector(time(abund11)), modwtasv11_biased$s4, xlab = "", ylab = "", type = "h",
     bty = "n", axes = F) > axis(1, cex.axis = 1.5) 
axis(2, cex.axis = 1.5, las = 2) +
mtext(side = 4, line = -2, cex = 1.2, outer = F, at = par("usr")[4], text = paste("s4", 
                                                                                  sep = ""), las = 2) 
title(main = "ASV 11")
abline(v = time(abund11)[range(which(!is.na(modwtasv11[[5]])))], lwd = 2, 
       col = "red") 
abline(v = time(abund11)[c(boundaries$start[5], boundaries$end[5])], col = "blue", 
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
    abline(v = seq(start(abund178)[1], end(abund178)[1], by = 1), col = "gray", 
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


## Apply Fourier analysis to compare-----
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
  
  perfect_sampling_freq<- seq(2004, 2013 + 11/12, by = 1/12)|>
    as_tibble_col(column_name = 'ideal_sampling') 
  

  
  # Data preparation
  asv11$decimal_date
  abund178_df <- data.frame(year = asv178$decimal_date, #seq(2004, 2013 + 11/12, by = 1/12),
                            abund178 = abund178,
                            Spline = Spline,
                            abund178_minus_Spline = abund178 - Spline)
  
  acf_df <- acf(abund178_df$abund178_minus_Spline, demean = TRUE)$acf |>
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
  
  

  