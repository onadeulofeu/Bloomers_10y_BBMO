# packages
library(fdrtool)
library(GeneCycle)
library(patchwork)
library(lubridate)
library(ggthemes)
library(stringr)
library(tidyverse)
library(phyloseq)

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

## maybe we should interpolate the missing samples...


## We need values transformed to CLR

four <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'zclr') |>
  dplyr::select(asv_num,abundance_value, month, year, day, family_f, decimal_date) |>
  arrange(decimal_date)

## Seasonality of each ASV that we believe can be a potential bloomer
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
  filter(asv_num == "asv1") %>%
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
