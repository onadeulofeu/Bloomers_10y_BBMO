# Differentiate seasonal from non-seasonal bloomers

## Adria's functions
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
  ts(df$Abundance, start = s,frequency = f) 
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


hammingdist.df <- as.data.frame(nwhamming) %>%
  rownames_to_column(var = "asv") %>%
  gather("otu", "hamming.dist", -asv) %>%
  group_by(asv) %>%
  filter(hamming.dist == min(hamming.dist))  %>%
  distinct(otu,asv, .keep_all = T)



## Fisher G-test to determine the significance of periodic components
library(GeneCycle)

## The time series was decomposed in three components, seasonal periodicity 
## (oscillation inside each period), trend (evolution over time) and residuals, 
## through local regression by the stl function.



## Autocorrelogram was calculated using the acf function
