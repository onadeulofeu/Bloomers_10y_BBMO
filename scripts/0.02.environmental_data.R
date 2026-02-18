# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                 Environmental data processing               ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## upload packages ----
library(vegan)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(multipanelfigure) #merge plots with different sizes
library(forcats) 
library(scales)
library(ggpmisc)
library(readxl) 

## upload functions ----
source('src/calculate_z_scores.R')

## environmental variables labs ----
labs_env <- as_labeller(c("day_length" = 'Day length' ,
                                "temperature" = 'Temperature',
                                "secchi"  = 'Turbididty\n(Secchi disck)',    
                                "salinity" = 'Salinity',
                                "chla_total" = 'Total chl-a',
                                "chla_3um" = 'Chl-a 3 um',
                                "PO4"  =   'Phosphate',       
                                "NH4"     = 'Ammonia',
                                "NO2"   =  'Nitrite',          
                                "NO3" = 'Nitrate',
                                "Si"  = 'Silicate',
                                "BP_FC1.55" = 'Bacterial production (µgC l-1 d-1)',     
                                "PNF_Micro" = 'Phototrophic nanoflagellates',
                                "PNF2_5um_Micro"  = 'Phototrophic nanoflagellates (2-5um)',
                                "PNF_5um_Micro"   = 'Phototrophic nanoflagellates (5um)',
                                "cryptomonas"   = 'Cryptomonas',
                                "micromonas"  = 'Micromonas',
                                "HNF_Micro"  =  'Heterotrophic nanoflagellates',   
                                "HNF2_5um_Micro" ='Heterotrophic nanoflagellates (2-5um)',
                                "HNF_5um_Micro"  ='Heterotrophic nanoflagellates (5um)',
                                "LNA" = 'LNA',
                                "HNA"  = 'HNA',
                                "prochlorococcus_FC" = 'Prochlorococcus',
                                "Peuk1"  = 'Picoeukaryotes population 1',       
                                "Peuk2"   = 'Picoeukaryotes population 2',
                                "bacteria_joint" = 'Bacterial abundance',
                                "synechococcus"= 'Synechococcus',
                          'low_vlp' = 'Low viruses',
                          'med_vlp' = 'Mid viruses',
                          'high_vlp' = 'High viruses',
                          'total_vlp' = 'Total viruses'))

##upload data ----
asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv')
asv_tab_rar <- read.csv2('data/asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()

##I upload the general metadata from the whole BBMO 20Y    
bbmo_20y <- read_xlsx('data/main_databaseMOSTREIGBLANES_March23_od.xlsx', skip = 0 ) |>
  as_tibble()

##I would like to add viruses to the analysis ----
bbmo_20y_v <- bbmo_20y |>
  dplyr::select(`NOM MOSTRA`, "Low VLP", "Med VLP" ,                           
                "High VLP" ,"total VLP")

bbmo_20y_v <- bbmo_20y_v |>
  rename(sample_id = `NOM MOSTRA`,
         low_vlp = "Low VLP",
         med_vlp = "Med VLP" ,  
         high_vlp = "High VLP" ,
         total_vlp = "total VLP") 
  
bbmo_20y_v_red <-  bbmo_20y_v |>
  dplyr::mutate(sample_id_ed = str_remove(sample_id, '_'),
                sample_id_ed2 = substr(sample_id_ed, 1, 8)) |>
  dplyr::filter(!is.na(sample_id_ed2)) |>
  dplyr::mutate(across(contains('_vlp'), as.numeric)) |>
  dplyr::filter(str_detect(sample_id_ed2, 'BL04|BL05|BL06|BL07|BL08|BL09|BL10|BL11|BL12|BL13|BL14')) |>
  dplyr::filter(!is.na(total_vlp)) |>
  dplyr::filter(!sample_id_ed2 == 'BL130709') |> #transecte DEVOTES
  dplyr::distinct(sample_id_ed2, total_vlp, low_vlp, med_vlp, high_vlp)

m_02_ed2 <- m_02 |>
  dplyr::select(sample_id, HNF_Micro, HNF2_5um_Micro, HNF_5um_Micro, decimal_date, date, year) |>
  separate(sample_id, sep = '_', into = c('sample_id', 'fraction', 'code'), remove = F)

m_vir_tb <- bbmo_20y_v_red |>
  right_join(m_02_ed2, by = c('sample_id_ed2' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

plot_virus_abund <- m_vir_tb |>
  ggplot(aes(date, total_vlp))+
  # geom_rect(data = harbour_restoration_dec, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_area(colour = '#182533')+
  labs(y = 'Virus/mL', x = 'Date')+
  geom_line(aes(date, low_vlp), color = '#69BFAE')+
  geom_line(aes(date, med_vlp), color = '#4E8AC7')+
  geom_line(aes(date, high_vlp), color = '#988B99')+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
  )+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

plot_virus_abund 

plot_HNF_abund <- m_vir_tb |>
  ggplot(aes(date, HNF_Micro))+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_line(colour = '#182533')+
  labs(y = 'cells/mL', x = 'Date')+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
  )+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

plot_HNF_abund

plot_grid(plot_virus_abund,
          plot_HNF_abund,
          cols = 1)

asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax$asv_num |>
  unique()

asv11_tb <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num == 'asv11') |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::select(abundance_value, sample_id) |>
  separate(sample_id, sep = '_', into = c('sample_id', 'fraction', 'code'), remove = F)

m_vir_tb |>
  left_join(asv11_tb, by = c('sample_id_ed2' = 'sample_id')) |>
  dplyr::mutate(lag_vir = lag(total_vlp)) |>
  ggplot(aes(lag_vir, abundance_value))+
  geom_point()

m_vir_tb |> 
  ggplot(aes(log10(total_vlp), log10(HNF_Micro))) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_poly_eq(
    aes(label = paste( ..rr.label.., ..p.value.label.., sep = "~~~")), 
    formula = y ~ x,
    method = "lm",
    parse = TRUE
  )

## we also have some data for the radiometer even though it has been broken for a period -----
### radiometer units llum (µE m-2 s-1)
m_bbmo_10y_ed <- m_bbmo_10y |>
  separate(sample_id, into = c('sample_id_sim', 'fraction', 'code'), sep = '_') |>
  dplyr::filter(fraction == '0.2') 

radiometer_data <- read_xlsx('data/raw/llum_Blanes_2005-14_ed.xlsx') |>
  mutate_all(~str_replace_all(., '/', 'NA')) |>
  pivot_longer(starts_with('BL'), names_to = 'sample_id', values_to = 'light') |>
  mutate_all(~str_replace_all(., '_', '')) |>
  right_join(m_bbmo_10y_ed, by = c('sample_id' = 'sample_id_sim'), relationship = "many-to-many") |>
  dplyr::mutate(light = as.numeric(light))

radiometer_data$sample_id

m_bbmo_10y_ed$samname

radiometer_data |>
  dplyr::filter(prof == '0') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, light))+
  geom_col(alpha = 0.6)+
  geom_line()+
  labs(x = 'Time', y = 'Light')+
  scale_x_datetime(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(), 
        axis.text.y = element_text(size = 5),
        panel.grid.major.y = element_blank())

radiometer_data |>
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d"),
                        prof = as.numeric(prof)) |>
  dplyr::filter(light != is.na(light) &
                  !is.na(prof)) |>
  ggplot(aes(month, prof))+
  geom_tile(aes(fill = light))+
  scale_y_reverse()+
  labs(x = 'Month', y = 'Depth (m)', fill = 'Light')+
  facet_wrap(vars(year), ncol = 1)+
  scale_fill_viridis_c(option = "magma")+
  theme_light()+
  theme(panel.grid.major.y = element_blank())

radiometer_data_mean <- radiometer_data |>
  dplyr::mutate(prof = as.numeric(prof)) |>
  dplyr::filter(prof >= 5) |>
  group_by(sample_id) |>
  dplyr::summarize(mean_surf_light = mean(light),
                sd_surf_light = sd(light)) |>
  left_join(m_bbmo_10y_ed, by = c('sample_id' = 'sample_id_sim')) |>
  dplyr::mutate(restoration = case_when(
    between(date, as.Date('2010-03-24', format = "%Y-%m-%d"), as.Date('2012-06-09', format = "%Y-%m-%d")) ~ 'restoration',
    date < as.Date('2010-03-24', format = "%Y-%m-%d") ~ 'no-restoration',
    date > as.Date('2012-06-09', format = "%Y-%m-%d") ~ 'no-restoration'
  )) |>
  dplyr::mutate(year = as.factor(year),
                restoration = as.factor(restoration))

##the same but by season
radiometer_data_mean_seas <- radiometer_data |>
  dplyr::mutate(prof = as.numeric(prof),
                      light = as.numeric(light)) |>
  dplyr::mutate(restoration = case_when(
    between(date, as.Date('2010-03-24', format = "%Y-%m-%d"), as.Date('2012-06-09', format = "%Y-%m-%d")) ~ 'restoration',
    date < as.Date('2010-03-24', format = "%Y-%m-%d") ~ 'no-restoration',
    date > as.Date('2012-06-09', format = "%Y-%m-%d") ~ 'no-restoration'
  )) |>
  dplyr::filter(prof >= 5) |>
  dplyr::filter(!is.na(light)) |>
  group_by(year, season, restoration) |>
  dplyr::summarize(mean_surf_light = mean(light),
                   sd_surf_light = sd(light)) |>
  dplyr::mutate(year = as.factor(year),
                restoration = as.factor(restoration))

radiometer_data_mean_seas$season <- radiometer_data_mean_seas$season |>
  factor(levels = c('winter', 'spring', 'summer', 'autumn'))

labs_season <- as_labeller(c('winter' = 'Winter',
                           'spring' = 'Spring', 
                           'summer' = 'Summer',
                           'autumn' = 'Autumn'))

radiometer_data_mean_seas |>
  ggplot(aes(x = season, y = mean_surf_light, shape = restoration)) +
  geom_boxplot(aes(x = season, y = mean_surf_light, shape = restoration), alpha = 0.1) +
  geom_point(aes(color = year, shape = restoration), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1)) +
  scale_color_manual(values = palette_years) +
  facet_wrap(vars(restoration))+
  labs(x = 'Season', color = 'Year', y = 'Light', shape = 'Harbour restoration')+
  scale_fill_manual(values = palette_years) +
  scale_x_discrete(labels = labs_season)+
  theme_light()

### restoration and no-restoration for other environmental data ----
m_bbmo_10y_ed_f <- m_bbmo_10y |>
  dplyr::select(date, year, season, month, sample_id) |>
  dplyr::mutate(restoration = case_when(
    between(date, as.Date('2010-03-24', format = "%Y-%m-%d"), as.Date('2012-06-09', format = "%Y-%m-%d")) ~ 'restoration',
    date < as.Date('2010-03-24', format = "%Y-%m-%d") ~ 'no-restoration',
    date > as.Date('2012-06-09', format = "%Y-%m-%d") ~ 'no-restoration'
  ))

bbmo_env_02_seas <- bbmo_env_02 |>
  pivot_longer(cols = (-sample_id), values_to = 'values', names_to = 'environmental_var') |>
  left_join(m_bbmo_10y_ed_f) |>
  dplyr::mutate(year = as.factor(year),
                restoration = as.factor(restoration))

bbmo_env_02_seas$season <- bbmo_env_02_seas$season |>
  factor(levels = c('winter', 'spring', 'summer', 'autumn'))

bbmo_env_02_seas$environmental_var |>
  unique()

bbmo_env_02_seas |>
  dplyr::filter(environmental_var %in% c("day_length",  "temperature", "secchi",  "salinity", "chla_total",  "chla_3um")) |>
  ggplot(aes(x = season, y = values, shape = restoration)) +
  geom_boxplot(aes(x = season, y = values, shape = restoration), alpha = 0.1) +
  geom_point(aes(color = year, shape = restoration), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1)) +
  scale_color_manual(values = palette_years) +
  facet_grid(environmental_var~restoration , scales = 'free_y', labeller = labs_env)+
  labs(x = 'Season', color = 'Year', y = '', shape = 'Harbour restoration')+
  scale_fill_manual(values = palette_years) +
  scale_x_discrete(labels = labs_season)+
  theme_light()+
  theme(text = element_text(size = 12), strip.text.x = element_blank())

bbmo_env_02_seas |>
  dplyr::filter(environmental_var %in% c("PO4", "NH4","NO2", "NO3","Si","BP_FC1.55")) |>
  ggplot(aes(x = season, y = values, shape = restoration)) +
  geom_boxplot(aes(x = season, y = values, shape = restoration), alpha = 0.1) +
  geom_point(aes(color = year, shape = restoration), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1)) +
  scale_color_manual(values = palette_years) +
  facet_grid(environmental_var~restoration , scales = 'free_y', labeller = labs_env)+
  labs(x = 'Season', color = 'Year', y = '', shape = 'Harbour restoration')+
  scale_fill_manual(values = palette_years) +
  scale_x_discrete(labels = labs_season)+
  theme_light()+
  theme(text = element_text(size = 12), strip.text.x = element_blank())

bbmo_env_02_seas |>
  dplyr::filter(environmental_var %in% c("PNF_Micro", "PNF2_5um_Micro" ,"PNF_5um_Micro", "dryptomonas", "micromonas" ,"HNF_Micro"  )) |>
  ggplot(aes(x = season, y = values, shape = restoration)) +
  geom_boxplot(aes(x = season, y = values, shape = restoration), alpha = 0.1) +
  geom_point(aes(color = year, shape = restoration), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1)) +
  scale_color_manual(values = palette_years) +
  facet_grid(environmental_var~restoration , scales = 'free_y', labeller = labs_env)+
  labs(x = 'Season', color = 'Year', y = '', shape = 'Harbour restoration')+
  scale_fill_manual(values = palette_years) +
  scale_x_discrete(labels = labs_season)+
  theme_light()+
  theme(text = element_text(size = 12), strip.text.x = element_blank())

bbmo_env_02_seas |>
  dplyr::filter(environmental_var %in% c("HNF2_5um_Micro", "HNF_5um_Micro", "LNA", "HNA", "prochlorococcus_FC", "Peuk1", "Peuk2" )) |>
  ggplot(aes(x = season, y = values, shape = restoration)) +
  geom_boxplot(aes(x = season, y = values, shape = restoration), alpha = 0.1) +
  geom_point(aes(color = year, shape = restoration), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1)) +
  scale_color_manual(values = palette_years) +
  facet_grid(environmental_var~restoration , scales = 'free_y', labeller = labs_env)+
  labs(x = 'Season', color = 'Year', y = '', shape = 'Harbour restoration')+
  scale_fill_manual(values = palette_years) +
  scale_x_discrete(labels = labs_season)+
  theme_light()+
  theme(text = element_text(size = 12), strip.text.x = element_blank())

bbmo_env_02_seas |>
  dplyr::filter(environmental_var %in% c( "bacteria_joint", "synechococcus", "low_vlp", "med_vlp", "high_vlp", "total_vlp" )) |>
  ggplot(aes(x = season, y = values, shape = restoration)) +
  geom_boxplot(aes(x = season, y = values, shape = restoration), alpha = 0.1) +
  geom_point(aes(color = year, shape = restoration), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1)) +
  scale_color_manual(values = palette_years) +
  facet_grid(environmental_var~restoration , scales = 'free_y', labeller = labs_env)+
  labs(x = 'Season', color = 'Year', y = '', shape = 'Harbour restoration')+
  scale_fill_manual(values = palette_years) +
  scale_x_discrete(labels = labs_season)+
  theme_light()+
  theme(text = element_text(size = 12), strip.text.x = element_blank())

### The restoration of the Blanes harbor strated on 24th March 2010 and finished on the 9th of june 2012----
harbour_restoration <- tibble(xmin = '2010-03-24', xmax = '2012-06-09') |>
  dplyr::mutate(date_min = as.POSIXct(xmin, format = "%Y-%m-%d"),
                date_max = (as.POSIXct(xmax, format = "%Y-%m-%d")))

sample_id <- asv_tab_all_bloo_z_tax$sample_id

bbmo_env <- asv_tab_all_bloo_z_tax |>
  dplyr::select(sample_id,
                day_length,
                #sampling_time
                temperature,
                secchi,
                salinity,
                chla_total,
                chla_3um,
                PO4,
                NH4, NO2, NO3,
                Si, BP_FC1.55,
                PNF_Micro, PNF2_5um_Micro,
                PNF_5um_Micro, 
                cryptomonas, micromonas,
                HNF_Micro, HNF2_5um_Micro,   
                HNF_5um_Micro ,  LNA,              
                HNA,   prochlorococcus_FC,
                Peuk1,   Peuk2,                 
                bacteria_joint, synechococcus) |>
  dplyr::select(-sample_id) |>
  dplyr::mutate_if( is.character, as.numeric) |>
  bind_cols(sample_id) |>
  rename(sample_id = '...28') |>
  distinct(sample_id,
           day_length,
           #sampling_time
           temperature,
           secchi,
           salinity,
           chla_total,
           chla_3um,
           PO4,
           NH4, NO2, NO3,
           Si, BP_FC1.55,
           PNF_Micro, PNF2_5um_Micro,
           PNF_5um_Micro, 
           cryptomonas, micromonas,
           HNF_Micro, HNF2_5um_Micro,   
           HNF_5um_Micro ,  LNA,              
           HNA,   prochlorococcus_FC,
           Peuk1,   Peuk2,                 
           bacteria_joint, synechococcus ) |>
  tidyr::separate(sample_id, into = c('sample_id_ed', 'filter', 'sequencing_num'), sep = '_', remove = FALSE) |>
  left_join(bbmo_20y_v_red, by = c('sample_id_ed' = 'sample_id_ed2'), relationship = "many-to-many") |>  ## add virus data to the environmental data 
  #rename('sample_id' = sample_id.x ) |>
  distinct(sample_id,
           day_length,
           #sampling_time
           temperature,
           secchi,
           salinity,
           chla_total,
           chla_3um,
           PO4,
           NH4, NO2, NO3,
           Si, BP_FC1.55,
           PNF_Micro, PNF2_5um_Micro,
           PNF_5um_Micro, 
           cryptomonas, micromonas,
           HNF_Micro, HNF2_5um_Micro,   
           HNF_5um_Micro ,  LNA,              
           HNA,   prochlorococcus_FC,
           Peuk1,   Peuk2,                 
           bacteria_joint, synechococcus,  low_vlp ,
           med_vlp ,  
           high_vlp ,
           total_vlp)


## Using the function created to normalize environmental data to z-scores----
bbmo_env_z <- 
  bbmo_env |>
  as_tibble() |>
  pivot_longer(cols = c( day_length,
                         #sampling_time
                         temperature,
                         secchi,
                         salinity,
                         chla_total,
                         chla_3um,
                         PO4,
                         NH4, NO2, NO3,
                         Si, BP_FC1.55,
                         PNF_Micro, PNF2_5um_Micro,
                         PNF_5um_Micro, 
                         cryptomonas, micromonas,
                         HNF_Micro, HNF2_5um_Micro,   
                         HNF_5um_Micro ,  LNA,              
                         HNA,   prochlorococcus_FC,
                         Peuk1,   Peuk2,                 
                         bacteria_joint, synechococcus, 
                         low_vlp ,
                         med_vlp ,  
                         high_vlp ,
                         total_vlp), values_to = 'env_values', names_to = 'environmental_variable') |>
  dplyr::filter(!is.na(env_values)) |>
  dplyr::mutate(env_values = as.numeric(env_values)) |>
  calculate_z_score(col = 'env_values', name = 'environmental_variable', group = 'environmental_variable') 
  
## reorder environmental factors ----
bbmo_env_z$environmental_variable <- factor(bbmo_env_z$environmental_variable, levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                                                                                          "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                                                                                          "bacteria_joint", "LNA", "HNA", "Peuk1",             
                                                                                          "Peuk2",
                                                                                          "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                                                                                          "HNF2_5um_Micro", "HNF_5um_Micro",
                                                                                          "low_vlp" ,
                                                                                          "med_vlp" ,  
                                                                                          "high_vlp" ,
                                                                                          "total_vlp"))
## plot environmental variables normalized by z-scores-----
bbmo_env_z |>
    left_join(m_bbmo_10y) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    ggplot(aes(date, z_score_environmental_variable))+
    geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                      ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
    geom_point(alpha = 0.6, size = 1/2)+
    geom_line()+
    labs(x = 'Time', y = 'z-scores')+
    facet_wrap(vars(environmental_variable), ncol = 2, scales = 'free_y', labeller = labs_env)+
    #facet_grid(vars(environmental_variable), scales = 'free_y', cols = vars(type_of_env), drop = T)+
    scale_y_continuous(labels = function(x) sprintf("%.1f", x),
                       expand = c(0,0))+
    scale_x_datetime(expand = c(0,0))+
    theme_bw()+
    theme(panel.grid.minor.y = element_blank(), 
          strip.background = element_blank(), 
          axis.text.y = element_text(size = 5),
          panel.grid.major.y = element_blank())

## I use a heat map so that maybe we can observe the changes better---
m_bbmo_10y_red <- m_bbmo_10y |>
  dplyr::select(date, decimal_date, sample_id, year)

bbmo_env_z_ed <- bbmo_env_z |>
  dplyr::filter(str_detect(sample_id, '0.2_')) |> #remove duplicates I had the same information for FL and PA
  dplyr::select(-env_values) |>
  pivot_wider( values_from = z_score_environmental_variable, values_fill = 0, names_from = environmental_variable) |>
  pivot_longer(cols = -sample_id, values_to = 'z_score_environmental_variable', names_to = 'environmental_variable') |>
left_join(m_bbmo_10y_red) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  arrange(decimal_date) |>
  group_by(environmental_variable, year) |>
  dplyr::mutate(consecutive_number = row_number()) |>
  mutate(type_of_env = case_when(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,
                                                               "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um") ~ 'physico_chem',
                                 environmental_variable %in% c( 
                                                                "bacteria_joint", "LNA", "HNA",
                                                                "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",
                                        "HNF2_5um_Micro", "HNF_5um_Micro", "prochlorococcus_FC", "Peuk1",
                                        "Peuk2", "synechococcus", "low_vlp" ,
                                        "med_vlp" ,  
                                        "high_vlp" ,
                                        "total_vlp") ~ 'biological'))
## reorder environmental factors ----
bbmo_env_z_ed$environmental_variable <- factor(bbmo_env_z_ed$environmental_variable, levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                                                                                          "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" , "BP_FC1.55", "bacteria_joint", "LNA", "HNA", "chla_total" ,  "chla_3um", 
                                                                                          "synechococcus", "prochlorococcus_FC", 
                                                                                           "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "HNF_Micro",         
                                                                                          "HNF2_5um_Micro", "HNF_5um_Micro", "Peuk1",             
                                                                                          "Peuk2",
                                                                                          "dryptomonas", "micromonas",
                                                                                          
                                                                                          "low_vlp" ,
                                                                                          "med_vlp" ,  
                                                                                          "high_vlp" ,
                                                                                          "total_vlp"))

# Define a custom color palette with fixed 0 color----
palete_gradient_cb5 <- c("#005300",
                         "#8ec68b",
                         "#FFFFFF", 
                         "#8011b8", 
                         "#4d009c")

bbmo_env_z_ed |>
  dplyr::filter(environmental_variable %in% 
                  c("day_length", "temperature" ,"secchi" , "salinity" ,
                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" )) %$% 
  z_score_environmental_variable |>
  range() ##which is the range of our palette, to have 0 at the middle

phyico_chemical <- bbmo_env_z_ed |>
  dplyr::filter(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,
                "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" )) |>
  ggplot(aes(consecutive_number, fct_rev(environmental_variable), fill = z_score_environmental_variable))+
  scale_y_discrete(labels = labs_env)+
  scale_fill_gradientn(colors = palete_gradient_cb5,
                       breaks=c(-5, 0, 5),
                       limits=c(-5.9,  5.9)) +
  facet_wrap(.~year, scales = 'free_x', nrow = 1, switch = 'x')+
  geom_tile(alpha = 1)+
  geom_vline(xintercept = seq(0.5, 12, by = 12), linetype = "dashed", color = "darkgrey") +
  scale_x_continuous(expand = c(0,0), labels = unique(bbmo_env_z_ed$year), breaks = unique(bbmo_env_z_ed$year))+
  labs(x = 'Year', y = '', fill = 'z-score')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 5), legend.position = 'bottom',
        panel.border = element_blank(), strip.background = element_blank(), plot.margin = unit(c(1,3,1,1), "mm"), 
        legend.key.size = unit(3, "mm"))

  # ggsave('physico_chem_z_scores_bbmo.pdf', phyico_chemical,
  #      path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/",
  #      width = 188,
  #      height = 60,
  #      units = 'mm')

bbmo_env_z_ed |>
    dplyr::filter(!(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,
                                                  "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" )))  %$% 
    z_score_environmental_variable |>
    range() ##which is the range of our palette, to have 0 at the middle
  
 biological <-
  bbmo_env_z_ed |>
  dplyr::filter(!(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,
                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ))) |>
  ggplot(aes(consecutive_number, fct_rev(environmental_variable), fill = z_score_environmental_variable))+
  scale_y_discrete(labels = labs_env)+
    scale_fill_gradientn(colors = palete_gradient_cb5,
                         breaks=c(-6, 0, 6),
                         limits=c(-6,  6), na.value="#4d009c")+ # I add this since I have one value which is outside the range, to define the color it should have
  facet_wrap(.~year, scales = 'free_x', nrow = 1, switch = 'x')+
  geom_tile(alpha = 1)+
  geom_vline(xintercept = seq(0.5, 12, by = 12), linetype = "dashed", color = "darkgrey") +
  scale_x_continuous(expand = c(0,0), labels = unique(bbmo_env_z_ed$year), breaks = unique(bbmo_env_z_ed$year))+
  labs(x = 'Year', y = '', fill = 'z-score')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 5), legend.position = 'bottom',
        panel.border = element_blank(), strip.background = element_blank(), plot.margin = unit(c(0,3,0,1), "mm"), 
        legend.key.size = unit(3, "mm"))

# ggsave('biological_bbmo.pdf', biological,
#        path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/",
#        width = 188,
#        height = 90,
#        units = 'mm')

## I divide them in different groups so that we can observe them better----
bbmo_env_z |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(environmental_variable %in% c('day_length', 'temperature', 'secchi', 'salinity', 'chla_total', 'chla_3um', 'PO4', 'NH4',
                                              'NO2', 'NO3', "Si", "BP_FC1.55" )) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, z_score_environmental_variable))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_point(alpha = 0.6, size = 1/2)+
  geom_line()+
  labs(x = 'Time', y = 'z-scores')+
  facet_wrap(vars(environmental_variable), ncol = 2, scales = 'free_y', labeller = labs_env)+
  scale_y_continuous(labels = function(x) sprintf("%.1f", x),
                     expand = c(0,0))+
  scale_x_datetime(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(), 
        axis.text.y = element_text(size = 5),
        panel.grid.major.y = element_blank())

bbmo_env_z |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(environmental_variable %in% c( "synechococcus", "prochlorococcus_FC", 
                                              "bacteria_joint", "LNA", "HNA", "Peuk1",             
                                              "Peuk2",
                                               "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                                              "HNF2_5um_Micro", "HNF_5um_Micro",
                                              "low_vlp" ,
                                              "med_vlp" ,  
                                              "high_vlp" ,
                                              "total_vlp")) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, z_score_environmental_variable))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_point(alpha = 0.6, size = 1/2)+
  geom_line()+
  labs(x = 'Time', y = 'z-scores')+
  facet_wrap(vars(environmental_variable), ncol = 2, scales = 'free_y', labeller = labs_env)+
  scale_y_continuous(labels = function(x) sprintf("%.1f", x),
                     expand = c(0,0))+
  scale_x_datetime(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(), 
        axis.text.y = element_text(size = 5),
        panel.grid.major.y = element_blank())

## plot the original data without the zscores normalization----
env_data_physico_chem <- bbmo_env_z |>
left_join(m_bbmo_10y) |>
  dplyr::filter(environmental_variable %in% c('day_length', 'temperature', 'secchi', 'salinity', 'chla_total', 'chla_3um', 'PO4', 'NH4',
                                              'NO2', 'NO3', "Si", "BP_FC1.55" )) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, env_values))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_point(alpha = 0.6, size = 1/2)+
  geom_line()+
  labs(x = 'Time', y = '')+
  facet_wrap(vars(environmental_variable), ncol = 2, scales = 'free_y', labeller = labs_env)+
  scale_y_continuous(labels = function(x) sprintf("%.1f", x),
                     expand = c(0,0))+
  scale_x_datetime(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(), 
        axis.text.y = element_text(size = 5),
        panel.grid.major.y = element_blank())

# ggsave('env_data_physico_chem.pdf', env_data_physico_chem,
#        path = '~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/',
#        width = 230,
#        height = 180,
#        units = 'mm')

env_data_cyto <- bbmo_env_z |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(environmental_variable %in% c( "synechococcus", "prochlorococcus_FC", 
                                               "bacteria_joint", "LNA", "HNA", "Peuk1",             
                                               "Peuk2",
                                               "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                                               "HNF2_5um_Micro", "HNF_5um_Micro",
                                               "low_vlp" ,
                                               "med_vlp" ,  
                                               "high_vlp" ,
                                               "total_vlp")) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, env_values))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_point(alpha = 0.6, size = 1/2)+
  geom_line()+
  labs(x = 'Time', y = 'Abundances (cells/mL) ')+
  facet_wrap(vars(environmental_variable), ncol = 2, scales = 'free_y', labeller = labs_env)+
  scale_y_continuous(labels = function(x) sprintf("%.1f", x),
                     expand = c(0,0))+
  scale_x_datetime(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(), 
        strip.background = element_blank(), 
        axis.text.y = element_text(size = 5),
        panel.grid.major.y = element_blank())

# ggsave('env_data_cyto.pdf', env_data_cyto,
#        path = '~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/',
#        width = 230,
#                height = 180,
#                units = 'mm')

# INTERPOLATE MISSING ENVIRONMENTAL VARIALBES TO HAVE A COMPLETE DATASET-------
 # data
 bbmo_env <- asv_tab_all_bloo_z_tax |>
   dplyr::select(decimal_date,
                 sample_id,
                 day_length,
                 #sampling_time
                 temperature,
                 secchi,
                 salinity,
                 chla_total,
                 chla_3um,
                 PO4,
                 NH4, NO2, NO3,
                 Si, BP_FC1.55,
                 PNF_Micro, PNF2_5um_Micro,
                 PNF_5um_Micro, 
                 cryptomonas, micromonas,
                 HNF_Micro, HNF2_5um_Micro,   
                 HNF_5um_Micro ,  LNA,              
                 HNA,   prochlorococcus_FC,
                 Peuk1,   Peuk2,                 
                 bacteria_joint, synechococcus) |>
   dplyr::select(-sample_id) |>
   dplyr::mutate_if( is.character, as.numeric) |>
   bind_cols(sample_id) |>
   rename(sample_id = '...28') |>
   distinct(decimal_date,
            sample_id,
            day_length,
            #sampling_time
            temperature,
            secchi,
            salinity,
            chla_total,
            chla_3um,
            PO4,
            NH4, NO2, NO3,
            Si, BP_FC1.55,
            PNF_Micro, PNF2_5um_Micro,
            PNF_5um_Micro, 
            cryptomonas, micromonas,
            HNF_Micro, HNF2_5um_Micro,   
            HNF_5um_Micro ,  LNA,              
            HNA,   prochlorococcus_FC,
            Peuk1,   Peuk2,                 
            bacteria_joint, synechococcus ) |>
   tidyr::separate(sample_id, into = c('sample_id_ed', 'filter', 'sequencing_num'), sep = '_', remove = FALSE)
 
 date_sample <- asv_tab_all_bloo_z_tax |>
   dplyr::filter(str_detect(sample_id, '_0.2_')) |>
   dplyr::filter(abundance_type == 'relative_abundance') |>
   dplyr::select(decimal_date, sample_id) |>
   group_by(decimal_date, sample_id) |>
   distinct()
 
 bbmo_env_ed <-  bbmo_env |>
   dplyr::filter(str_detect(sample_id, '_0.2_')) |>
   left_join(date_sample) |>
   dplyr::select(-sample_id)
 
## Interpolation for just one environmental variable-----
 #### we do it with the loess function
 # Select data without missing values
 env_to_fit <- bbmo_env_ed |>
   select(day_length, decimal_date) %>%
   filter(!is.na(day_length))
 
 # Select all data including missing values
 env_fitted <- bbmo_env_ed %>%
   select(day_length, decimal_date)
 
 # Fit LOESS model
 fit <- loess(day_length ~ decimal_date, data = env_to_fit, span = 0.1)
 
 # Predict missing values
 env_fitted$day_length_interpolated <- predict(fit, newdata = env_fitted)
 
 ## plot to observe the predicted values and the previous ones before the fit
env_fitted |>
  pivot_longer(cols = !decimal_date) |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  facet_grid(vars(name))+
  theme_bw()
 
## I do it for all the environmental variables with missing values ------
# Function to interpolate missing values for a single variable
source('src/interpolate_missing_values.R')
# interpolate_missing <- function(data, variable_to_inter, span_value) {
#   
#   env_to_fit <- data |>
#     dplyr::filter(variable == {{variable_to_inter}}) |>
#     dplyr::select(value, decimal_date) |>
#     dplyr::filter(!is.na(value))
# 
#   env_fitted <- data %>%
#     dplyr::filter(variable == {{variable_to_inter}}) |>
#     dplyr::select(value, decimal_date)
#   
#   fit <- loess(value ~ decimal_date, data = env_to_fit, span = span_value)
#   
#   env_fitted[[paste0(variable_to_inter, "_interpolated")]] <- predict(fit, newdata = env_fitted)
#   
#   env_fitted <- env_fitted |>
#     rename(!!variable_to_inter := value)
#   
#   return(env_fitted)
# }

# Pivot to long format
bbmo_env_long <- bbmo_env_ed |>
  pivot_longer(cols = -decimal_date, names_to = "variable", values_to = "value")

## I count the number of NAs that I have for each environmental variable 
na_counts <- bbmo_env_long |>
  group_by(variable) |>
  summarise(na_count = sum(is.na(value))) |>
  dplyr::filter(na_count >= 1)

no_na_counts <- bbmo_env_long |>
  group_by(variable) |>
  summarise(na_count = sum(is.na(value))) |>
  dplyr::filter(na_count == 0) # I don't need to interpolate data for synechococcus or bacterial abundance

## I will not interpolate data that have many missing values (>10% of the dataset)
remove <- na_counts |>
  dplyr::filter(na_count > 12/120*100) #basically the viruses dataset

## between 5%-10% of the dataset missing
pay_attention <- na_counts |>
  dplyr::filter(na_count > 6/120*100) |> #basically the viruses dataset
  dplyr::filter(!variable %in% remove$variable)

## no consecutive missing values ----
just_one <- bbmo_env_long |>
  group_by(variable) |>
  summarise(na_count = sum(is.na(value))) |>
  dplyr::filter(na_count == 1)

na_consecutive <- bbmo_env_long |>
  arrange(decimal_date) |>
  group_by(variable) |>
  mutate(n_sample = row_number()) |>
  ungroup() |>
  dplyr::filter(is.na(value)) |>
  dplyr::filter(!variable %in% remove$variable) |>
  dplyr::filter(!variable %in% just_one$variable)

na_consecutive |>
  distinct(variable)

consecutive_missing_values <- c('chla_3um', 'HNA', 'LNA', 'Peuk1', 'Peuk2', 'prochlorococcus_FC', 
                                        'salinity')

# Get unique variables (to interpolate)
variables <- unique(bbmo_env_long$variable) |>
  as_tibble_col(column_name = 'variable') |>
  dplyr::filter(!variable %in% remove$variable) |>
  dplyr::filter(!variable %in% no_na_counts$variable) |>
  dplyr::filter(!variable %in% consecutive_missing_values) |>
  dplyr::filter(variable != 'secchi') |> # i remove the secchi because i do not think it makes sense to interpolate the turbidity (no pattern)
  as_vector()

# Apply interpolation function to each variable ( i do it individually because the span value should be adapted to each environmental variable)------

### day length----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "day_length",
                    span_value = 0.1)

interpolated_data |>
  ggplot(aes(day_length, day_length_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(day_length) ~ day_length_interpolated,
                                     !is.na(day_length) ~ day_length)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == 'day_length') |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_dl <- interpolated_data |>
  dplyr::mutate(day_length_no_nas = case_when(is.na(day_length) ~ day_length_interpolated,
                                     !is.na(day_length) ~ day_length)) |>
  dplyr::select(decimal_date, day_length_no_nas)


### temperature----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "temperature",
                                         span_value = 0.09)

interpolated_data |>
  ggplot(aes(temperature, temperature_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(temperature) ~ temperature_interpolated,
                                     !is.na(temperature) ~ temperature)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == 'temperature') |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_t <- interpolated_data |>
  dplyr::mutate(temperature_no_nas = case_when(is.na(temperature) ~ temperature_interpolated,
                                              !is.na(temperature) ~ temperature)) |>
  dplyr::select(decimal_date, temperature_no_nas)

### secchi ----
# interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "secchi",
#                                          span_value = 0.1)
# 
# interpolated_data |>
#   ggplot(aes(secchi, secchi_interpolated))+
#   geom_point()+
#   geom_smooth(method = 'lm')
# 
# ## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
# interpolated_data |>
#   dplyr::mutate(new_data = case_when(is.na(secchi) ~ secchi_interpolated,
#                                      !is.na(secchi) ~ secchi)) |>
#   ggplot(aes(decimal_date, new_data))+
#   geom_point()+
#   geom_line()+
#   theme_bw()
# 
# new_data_with_inter_values_secc <- interpolated_data |>
#   dplyr::mutate(secchi_no_nas = case_when(is.na(secchi) ~ secchi_interpolated,
#                                               !is.na(secchi) ~ secchi)) |>
#   dplyr::select(decimal_date, secchi_no_nas)

### chl-a total ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "chla_total",
                                         span_value = 0.09)

interpolated_data |>
  ggplot(aes(chla_total, chla_total_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(chla_total) ~ chla_total_interpolated,
                                     !is.na(chla_total) ~ chla_total)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == 'chla_total') |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_chla_tot <- interpolated_data |>
  dplyr::mutate(chla_total_no_nas = case_when(is.na(chla_total) ~ chla_total_interpolated,
                                              !is.na(chla_total) ~ chla_total)) |>
  dplyr::select(decimal_date, chla_total_no_nas)

### PO4 ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "PO4",
                                         span_value = 0.085)

interpolated_data |>
  ggplot(aes(PO4, PO4_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(PO4) ~ PO4_interpolated,
                                     !is.na(PO4) ~ PO4)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == 'PO4') |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_po4 <- interpolated_data |>
  dplyr::mutate(PO4_no_nas = case_when(is.na(PO4) ~ PO4_interpolated,
                                              !is.na(PO4) ~ PO4)) |>
  dplyr::select(decimal_date, PO4_no_nas)

### NH4 ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "NH4",
                                         span_value = 1)

interpolated_data |>
  ggplot(aes(NH4, NH4_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(NH4) ~ NH4_interpolated,
                                     !is.na(NH4) ~ NH4)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == 'NH4') |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_nh4 <- interpolated_data |>
  dplyr::mutate(NH4_no_nas = case_when(is.na(NH4) ~ NH4_interpolated,
                                              !is.na(NH4) ~ NH4)) |>
  dplyr::select(decimal_date, NH4_no_nas)

### NO2 total ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "NO2",
                                         span_value = 2)

interpolated_data |>
  ggplot(aes(NO2, NO2_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(NO2) ~ NO2_interpolated,
                                     !is.na(NO2) ~ NO2)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == 'NO2') |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_no2 <- interpolated_data |>
  dplyr::mutate(NO2_no_nas = case_when(is.na(NO2) ~ NO2_interpolated,
                                              !is.na(NO2) ~ NO2)) |>
  dplyr::select(decimal_date, NO2_no_nas)

### NO3 total ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "NO3",
                                         span_value = 0.09)

interpolated_data |>
  ggplot(aes(NO3, NO3_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(NO3) ~ NO3_interpolated,
                                     !is.na(NO3) ~ NO3)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == 'NO3') |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_no3 <- interpolated_data |>
  dplyr::mutate(NO3_no_nas = case_when(is.na(NO3) ~ NO3_interpolated,
                                              !is.na(NO3) ~ NO3)) |>
  dplyr::select(decimal_date, NO3_no_nas)

### Si ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "Si",
                                         span_value = 1)

interpolated_data |>
  ggplot(aes(Si, Si_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(Si) ~ Si_interpolated,
                                     !is.na(Si) ~ Si)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == 'Si') |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_si <- interpolated_data |>
  dplyr::mutate(Si_no_nas = case_when(is.na(Si) ~ Si_interpolated,
                                              !is.na(Si) ~ Si)) |>
  dplyr::select(decimal_date, Si_no_nas)

### BP_FC1.55 ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "BP_FC1.55",
                                         span_value = 1)

interpolated_data |>
  ggplot(aes(BP_FC1.55, BP_FC1.55_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(BP_FC1.55) ~ BP_FC1.55_interpolated,
                                     !is.na(BP_FC1.55) ~ BP_FC1.55)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "BP_FC1.55") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_bp <- interpolated_data |>
  dplyr::mutate(BP_FC1.55_no_nas = case_when(is.na(BP_FC1.55) ~ BP_FC1.55_interpolated,
                                              !is.na(BP_FC1.55) ~ BP_FC1.55)) |>
  dplyr::select(decimal_date, BP_FC1.55_no_nas)

### PNF_Micro ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "PNF_Micro",
                                         span_value = 0.09)

interpolated_data |>
  ggplot(aes(PNF_Micro, PNF_Micro_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(PNF_Micro) ~ PNF_Micro_interpolated,
                                     !is.na(PNF_Micro) ~ PNF_Micro)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "PNF_Micro") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_PNF_Micro <- interpolated_data |>
  dplyr::mutate(PNF_Micro_no_nas = case_when(is.na(PNF_Micro) ~ PNF_Micro_interpolated,
                                              !is.na(PNF_Micro) ~ PNF_Micro)) |>
  dplyr::select(decimal_date, PNF_Micro_no_nas)

### PNF2_5um_Micro ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "PNF2_5um_Micro",
                                         span_value = 0.9)

interpolated_data |>
  ggplot(aes(PNF2_5um_Micro, PNF2_5um_Micro_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(PNF2_5um_Micro) ~ PNF2_5um_Micro_interpolated,
                                     !is.na(PNF2_5um_Micro) ~ PNF2_5um_Micro)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "PNF2_5um_Micro") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_PNF2_5um_Micro <- interpolated_data |>
  dplyr::mutate(PNF2_5um_Micro_no_nas = case_when(is.na(PNF2_5um_Micro) ~ PNF2_5um_Micro_interpolated,
                                              !is.na(PNF2_5um_Micro) ~ PNF2_5um_Micro)) |>
  dplyr::select(decimal_date, PNF2_5um_Micro_no_nas)

### PNF_5um_Micro ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "PNF_5um_Micro",
                                         span_value = 2)

interpolated_data |>
  ggplot(aes(PNF_5um_Micro, PNF_5um_Micro_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(PNF_5um_Micro) ~ PNF_5um_Micro_interpolated,
                                     !is.na(PNF_5um_Micro) ~ PNF_5um_Micro)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "PNF_5um_Micro") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_PNF_5um_Micro <- interpolated_data |>
  dplyr::mutate(PNF_5um_Micro_no_nas = case_when(is.na(PNF_5um_Micro) ~ PNF_5um_Micro_interpolated,
                                              !is.na(PNF_5um_Micro) ~ PNF_5um_Micro)) |>
  dplyr::select(decimal_date, PNF_5um_Micro_no_nas)


### cryptomonas ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "cryptomonas",
                                         span_value = 0.09)

interpolated_data |>
  ggplot(aes(cryptomonas, cryptomonas_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(cryptomonas) ~ cryptomonas_interpolated,
                                     !is.na(cryptomonas) ~ cryptomonas)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "cryptomonas") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_cryptomonas <- interpolated_data |>
  dplyr::mutate(cryptomonas_no_nas = case_when(is.na(cryptomonas) ~ cryptomonas_interpolated,
                                              !is.na(cryptomonas) ~ cryptomonas)) |>
  dplyr::select(decimal_date, cryptomonas_no_nas)

### micromonas ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "micromonas",
                                         span_value = 1)

interpolated_data |>
  ggplot(aes(micromonas, micromonas_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(micromonas) ~ micromonas_interpolated,
                                     !is.na(micromonas) ~ micromonas)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "micromonas") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_micromonas <- interpolated_data |>
  dplyr::mutate(micromonas_no_nas = case_when(is.na(micromonas) ~ micromonas_interpolated,
                                              !is.na(micromonas) ~ micromonas)) |>
  dplyr::select(decimal_date, micromonas_no_nas)

### HNF_Micro ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "HNF_Micro",
                                         span_value = 0.09)

interpolated_data |>
  ggplot(aes(HNF_Micro, HNF_Micro_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(HNF_Micro) ~ HNF_Micro_interpolated,
                                     !is.na(HNF_Micro) ~ HNF_Micro)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "HNF_Micro") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_HNF_Micro <- interpolated_data |>
  dplyr::mutate(HNF_Micro_no_nas = case_when(is.na(HNF_Micro) ~ HNF_Micro_interpolated,
                                              !is.na(HNF_Micro) ~ HNF_Micro)) |>
  dplyr::select(decimal_date, HNF_Micro_no_nas)

### HNF2_5um_Micro ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "HNF2_5um_Micro",
                                         span_value = 0.09)

interpolated_data |>
  ggplot(aes(HNF2_5um_Micro, HNF2_5um_Micro_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(HNF2_5um_Micro) ~ HNF2_5um_Micro_interpolated,
                                     !is.na(HNF2_5um_Micro) ~ HNF2_5um_Micro)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "HNF2_5um_Micro") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_HNF2_5um_Micro <- interpolated_data |>
  dplyr::mutate(HNF2_5um_Micro_no_nas = case_when(is.na(HNF2_5um_Micro) ~ HNF2_5um_Micro_interpolated,
                                              !is.na(HNF2_5um_Micro) ~ HNF2_5um_Micro)) |>
  dplyr::select(decimal_date, HNF2_5um_Micro_no_nas)

### HNF_5um_Micro ----
interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "HNF_5um_Micro",
                                         span_value = 0.09)

interpolated_data |>
  ggplot(aes(HNF_5um_Micro, HNF_5um_Micro_interpolated))+
  geom_point()+
  geom_smooth(method = 'lm')

## what I do is I substitute the NAs for the interpolated value, but I maintain the real values in the other situations
interpolated_data |>
  dplyr::mutate(new_data = case_when(is.na(HNF_5um_Micro) ~ HNF_5um_Micro_interpolated,
                                     !is.na(HNF_5um_Micro) ~ HNF_5um_Micro)) |>
  ggplot(aes(decimal_date, new_data))+
  geom_point()+
  geom_line()+
  theme_bw()

bbmo_env_long |>
  dplyr::filter(variable == "HNF_5um_Micro") |>
  ggplot(aes(decimal_date, value))+
  geom_point()+
  geom_line()+
  theme_bw()

new_data_with_inter_values_HNF_5um_Micro <- interpolated_data |>
  dplyr::mutate(HNF_5um_Micro_no_nas = case_when(is.na(HNF_5um_Micro) ~ HNF_5um_Micro_interpolated,
                                              !is.na(HNF_5um_Micro) ~ HNF_5um_Micro)) |>
  dplyr::select(decimal_date, HNF_5um_Micro_no_nas)


## new environmental dataset with no missing values-----
env_data_interpolated_values <- bind_cols(new_data_with_inter_values_t, new_data_with_inter_values_dl, new_data_with_inter_values_chla_tot,
          new_data_with_inter_values_po4, new_data_with_inter_values_nh4,
          new_data_with_inter_values_no2, new_data_with_inter_values_no3,
          new_data_with_inter_values_si, #it creates a peak which could be true or artefactual just one point.
          new_data_with_inter_values_bp,
          new_data_with_inter_values_PNF_Micro, new_data_with_inter_values_PNF2_5um_Micro, # nit creates a peak which could be true or artefactual just one point in PNF2-5
          new_data_with_inter_values_PNF_5um_Micro, 
          new_data_with_inter_values_cryptomonas, 
          new_data_with_inter_values_micromonas, # it creates a peak which could be true or artefactual just one point.
          new_data_with_inter_values_HNF_Micro, new_data_with_inter_values_HNF2_5um_Micro,
          new_data_with_inter_values_HNF_5um_Micro)

## add the columns for the environmental data with no missing values---
no_na_counts$variable

decimal_date <- env_data_interpolated_values$decimal_date...1 |>
  as_tibble_col(column_name = 'decimal_date')

env_data_interpolated_values_all <- bbmo_env_ed |>
  dplyr::select(no_na_counts$variable, decimal_date) |>
  bind_cols(env_data_interpolated_values) |>
  dplyr::select(-starts_with('decimal_date')) |>
  bind_cols(decimal_date)

#write.csv2(env_data_interpolated_values_all, 'data/env_data/env_data_interpolated_values_all.csv')

### calculate z-scores for the interpolated data----
env_data_interpolated_values_all_z_score <- env_data_interpolated_values_all |>
  pivot_longer(cols = -decimal_date, names_to = 'environmental_variable', values_to = 'env_values') |>
  dplyr::mutate(env_values = as.numeric(env_values)) |>
  calculate_z_score(col = 'env_values', name = 'environmental_variable', group = 'environmental_variable')

#write.csv2(env_data_interpolated_values_all_z_score, 'data/env_data/env_data_interpolated_values_all_z_score.csv')

## visually inspect the results and see if I'm adding weird values----
env_data_interpolated_values_all_z_score

year_date <- asv_tab_all_bloo_z_tax |>
  dplyr::select(decimal_date, year) |>
  distinct(decimal_date, year)

env_data_interpolated_values_all_z_score_ed <- env_data_interpolated_values_all_z_score |>
  left_join(year_date) |>
  arrange(decimal_date) |>
  group_by(environmental_variable, year) |>
  dplyr::mutate(consecutive_number = row_number()) |>
  dplyr::mutate(environmental_variable = str_replace(environmental_variable, "_no_nas", ""))

env_data_interpolated_values_all_z_score |>
  dplyr::filter(environmental_variable %in% 
                  c("day_length_no_nas", "temperature_no_nas" , 
                    "PO4_no_nas",  "NH4_no_nas" ,  "NO2_no_nas" , "NO3_no_nas" , "Si_no_nas" )) %$% 
  z_score_environmental_variable |>
  range() ##which is the range of our palette, to have 0 at the middle

## reorder environmental factors ----
env_data_interpolated_values_all_z_score_ed$environmental_variable <- factor(env_data_interpolated_values_all_z_score_ed$environmental_variable, levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                                                                                                "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" , "BP_FC1.55", "bacteria_joint", "LNA", "HNA", "chla_total" ,  "chla_3um", 
                                                                                                "synechococcus", "prochlorococcus_FC", 
                                                                                                "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "HNF_Micro",         
                                                                                                "HNF2_5um_Micro", "HNF_5um_Micro", "Peuk1",             
                                                                                                "Peuk2",
                                                                                                "cryptomonas", "micromonas",
                                                                                                
                                                                                                "low_vlp" ,
                                                                                                "med_vlp" ,  
                                                                                                "high_vlp" ,
                                                                                                "total_vlp"))

phyico_chemical <- env_data_interpolated_values_all_z_score_ed |>
  dplyr::filter(environmental_variable %in% c("day_length", "temperature" ,
                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" )) |>
  ggplot(aes(consecutive_number, fct_rev(environmental_variable), fill = z_score_environmental_variable))+
  scale_y_discrete(labels = labs_env)+
  scale_fill_gradientn(colors = palete_gradient_cb5,
                       breaks=c(-5, 0, 5),
                       limits=c(-5.9,  5.9)) +
  facet_wrap(.~year, scales = 'free_x', nrow = 1, switch = 'x')+
  geom_tile(alpha = 1)+
  geom_vline(xintercept = seq(0.5, 12, by = 12), linetype = "dashed", color = "darkgrey") +
  #scale_x_datetime(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), labels = unique(env_data_interpolated_values_all_z_score_ed$year), 
                     breaks = unique(env_data_interpolated_values_all_z_score_ed$year))+
  labs(x = 'Year', y = '', fill = 'z-score')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 5), legend.position = 'bottom',
        panel.border = element_blank(), strip.background = element_blank(), plot.margin = unit(c(1,3,1,1), "mm"), 
        legend.key.size = unit(3, "mm"))

# ggsave('physico_chem_z_scores_bbmo_interpolated.pdf', phyico_chemical,
#      path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/",
#      width = 188,
#      height = 60,
#      units = 'mm')

biological <-
  env_data_interpolated_values_all_z_score_ed |>
  dplyr::filter(!(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,
                                                "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ))) |>
  ggplot(aes(consecutive_number, fct_rev(environmental_variable), fill = z_score_environmental_variable))+
  scale_y_discrete(labels = labs_env)+
  scale_fill_gradientn(colors = palete_gradient_cb5,
                       breaks=c(-6, 0, 6),
                       limits=c(-6,  6), na.value="#4d009c")+ # I add this since I have one value which is outside the range, to define the color it should have
  facet_wrap(.~year, scales = 'free_x', nrow = 1, switch = 'x')+
  geom_tile(alpha = 1)+
  geom_vline(xintercept = seq(0.5, 12, by = 12), linetype = "dashed", color = "darkgrey") +
  #scale_x_datetime(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), labels = unique(env_data_interpolated_values_all_z_score_ed$year), 
                     breaks = unique(env_data_interpolated_values_all_z_score_ed$year))+
  labs(x = 'Year', y = '', fill = 'z-score')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 5), legend.position = 'bottom',
        panel.border = element_blank(), strip.background = element_blank(), plot.margin = unit(c(0,3,0,1), "mm"), 
        legend.key.size = unit(3, "mm"))

# ggsave('biological_bbmo_interpolated.pdf', biological,
#        path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/",
#        width = 188,
#        height = 90,
#        units = 'mm')

# When we have a blooming event do we observe a general increase in the bacterial production?-----
 sample_id <- asv_tab_all_bloo_z_tax$sample_id
 
 bbmo_env <- asv_tab_all_bloo_z_tax |>
   dplyr::select(decimal_date,
     sample_id,
                 day_length,
                 #sampling_time
                 temperature,
                 secchi,
                 salinity,
                 chla_total,
                 chla_3um,
                 PO4,
                 NH4, NO2, NO3,
                 Si, BP_FC1.55,
                 PNF_Micro, PNF2_5um_Micro,
                 PNF_5um_Micro, 
                 cryptomonas, micromonas,
                 HNF_Micro, HNF2_5um_Micro,   
                 HNF_5um_Micro ,  LNA,              
                 HNA,   prochlorococcus_FC,
                 Peuk1,   Peuk2,                 
                 bacteria_joint, synechococcus) |>
   dplyr::select(-sample_id) |>
   dplyr::mutate_if( is.character, as.numeric) |>
   bind_cols(sample_id) |>
   rename(sample_id = '...28') |>
   distinct(decimal_date,
     sample_id,
            day_length,
            #sampling_time
            temperature,
            secchi,
            salinity,
            chla_total,
            chla_3um,
            PO4,
            NH4, NO2, NO3,
            Si, BP_FC1.55,
            PNF_Micro, PNF2_5um_Micro,
            PNF_5um_Micro, 
            cryptomonas, micromonas,
            HNF_Micro, HNF2_5um_Micro,   
            HNF_5um_Micro ,  LNA,              
            HNA,   prochlorococcus_FC,
            Peuk1,   Peuk2,                 
            bacteria_joint, synechococcus ) |>
   tidyr::separate(sample_id, into = c('sample_id_ed', 'filter', 'sequencing_num'), sep = '_', remove = FALSE)
 
 bbmo_env |>
   dplyr::filter(str_detect(sample_id, '_0.2_')) |>
   ggplot(aes(decimal_date, BP_FC1.55))+
   geom_point()+
   scale_y_continuous(expand = c(0,0))+
   theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(decimal_date, BP_FC1.55, abundance_value) |>
  dplyr::summarize(bloom_event = case_when(any(z_score_ra >= 1.96 &
                                                abundance_value >= 0.1 ) ~ 'bloom',
                                           TRUE ~ 'no-bloo')) |>
  ggplot(aes(sqrt(abundance_value), sqrt(BP_FC1.55)))+
  geom_point(aes(color = bloom_event))+ 
  scale_y_continuous(expand = c(0,0))+
  theme_bw()
 
bloo_events <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(decimal_date, BP_FC1.55) |>
  dplyr::summarize(bloom_event = case_when(any(z_score_ra >= 1.96 &
                                                 abundance_value >= 0.1 ) ~ 'bloom',
                                           TRUE ~ 'no-bloo')) |>
  left_join(bbmo_env)

asv_tab_all_bloo_z_tax |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |> 
  dplyr::select(decimal_date, abundance_value) |>
  dplyr::group_by(decimal_date) |>
  dplyr::summarize(max_abundance = max(abundance_value, na.rm = TRUE)) |>
  left_join(bloo_events) |>
  ggplot(aes(log(max_abundance), log(BP_FC1.55)))+
  geom_point(aes(color = bloom_event))+ 
  scale_y_continuous(expand = c(0,0))+
  geom_smooth(method = 'loess')+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |> 
  dplyr::select(decimal_date, abundance_value) |>
  dplyr::group_by(decimal_date) |>
  dplyr::summarize(max_abundance = max(abundance_value, na.rm = TRUE)) |>
  left_join(bloo_events) |>
  dplyr::filter(!is.na(BP_FC1.55)) |>
  ungroup() |>
  ggplot(aes(bloom_event, sqrt(BP_FC1.55)))+
  geom_point(position = position_jitter(width = 0.1), aes(color = bloom_event))+
  geom_violin(aes(color = bloom_event, fill = bloom_event))+
  geom_boxplot(aes(color = bloom_event), notch = T)+ 
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  scale_y_continuous(expand = c(0,0))+
  labs(x = '', y = 'Bacterial production (µgC l-1 d-1)', size = 'Maximal\nrelative\nabundance (%)')+
  theme_bw()
  
### When we have a blooming event do we observe a general increase/decrease in some environmental variables?------
bloo_events <- bloo_events |>
  dplyr::mutate(HNA_perc = HNA/(LNA+HNA))

asv_tab_all_bloo_z_tax |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |> 
  dplyr::select(decimal_date, abundance_value) |>
  dplyr::group_by(decimal_date) |>
  dplyr::summarize(max_abundance = max(abundance_value, na.rm = TRUE)) |>
  left_join(bloo_events) |>
  ungroup() |>
  dplyr::select(decimal_date, bloom_event,  temperature, day_length, secchi, salinity, chla_total, chla_3um, PO4, NH4, NO2, NO3, Si, PNF_Micro,
                PNF2_5um_Micro, PNF_5um_Micro, cryptomonas, micromonas, HNF_Micro, HNF2_5um_Micro, HNF_5um_Micro, LNA, HNA, prochlorococcus_FC,
                Peuk1, Peuk2, bacteria_joint, synechococcus, max_abundance, HNA_perc) |>
  pivot_longer(cols = -c(decimal_date, bloom_event, max_abundance)) |>
  ggplot(aes(bloom_event, as.numeric(value)))+
  geom_point( position = position_jitter(width = 0.05))+ 
  scale_x_discrete(label = c('Bloom', 'No bloom'))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  facet_wrap(vars(name), scale = 'free')+
  geom_boxplot(alpha = 0.3, notch = T)+
  labs(x = '', y = '', size = 'Maximal\nrelative\nabundance (%)')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 6))

### environmental distance between samples ----
## Upload environmental data to model with interpolated missing variables----
env_data_interpolated_values_all <- read.csv2('data/env_data/env_data_interpolated_values_all.csv') |>
  rename(sample_id_num = X)

m_02_red2 <- m_02 |>
  dplyr::select(sample_id, date, fraction, sample_id_num) |>
  dplyr::select(-sample_id)

env_data_interpolated_values_all |>
  colnames()

env_data_interpolated_values_all |>
  dim()

env_data_interpolated_values_all_scaled <- env_data_interpolated_values_all |>
  arrange(decimal_date) |>
  dplyr::select(-decimal_date) |>
  dplyr::mutate_if(is.numeric, scale)

decimal_date_tb_unique <- env_data_interpolated_values_all |>
  arrange(decimal_date) |>
  dplyr::select(decimal_date, sample_id_num) 

# Dissimilarity matrix
distances_env <- env_data_interpolated_values_all_scaled  |>
  dplyr::select(-BP_FC1.55_no_nas, -sample_id_num) |>
  stats::dist( method = "euclidean")

euclidean_distance <- distances_env |>
  as.matrix() |>
  as_tibble() |>
  bind_cols(decimal_date_tb_unique) |>
  pivot_longer(cols = -c('decimal_date', 'sample_id_num'), values_to = 'euclidean_distance', names_to = 'sample_id_num_2')

decimal_date_tb_unique <- decimal_date_tb_unique |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

euclidean_distance_tb_g <- euclidean_distance |>
  dplyr::mutate(sample_distance = (as.numeric(sample_id_num_2) - as.numeric(sample_id_num))) |>
  dplyr::filter(sample_distance == 1) |>
  dplyr::select(-decimal_date) |>
  left_join(decimal_date_tb_unique, by = c('sample_id_num_2' = 'sample_id_num'))

euclidean_distance_tb_g |>
  ggplot(aes(decimal_date, euclidean_distance))+
  geom_line()

## less env variables -----
# Dissimilarity matrix (non bio parammeters nutrients temperature and day length)
env_data_interpolated_values_all_scaled |>
  colnames()

distances_env <- env_data_interpolated_values_all_scaled  |>
  dplyr::select("temperature_no_nas"  ,  "day_length_no_nas"   ,  "chla_total_no_nas"    ,
                "PO4_no_nas"    ,        "NH4_no_nas"       ,     "NO2_no_nas"    ,        "NO3_no_nas"     ,       "Si_no_nas" ) |>
  stats::dist( method = "euclidean")

euclidean_distance <- distances_env |>
  as.matrix() |>
  as_tibble() |>
  bind_cols(decimal_date_tb_unique) |>
  pivot_longer(cols = -c('decimal_date', 'sample_id_num'), values_to = 'euclidean_distance', names_to = 'sample_id_num_2')

decimal_date_tb_unique <- decimal_date_tb_unique |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

euclidean_distance_tb_phch <- euclidean_distance |>
  dplyr::mutate(sample_distance = (as.numeric(sample_id_num_2) - as.numeric(sample_id_num))) |>
  dplyr::filter(sample_distance == 1) |>
  dplyr::select(-decimal_date) |>
  left_join(decimal_date_tb_unique, by = c('sample_id_num_2' = 'sample_id_num'))

euclidean_distance_tb_phch |>
  ggplot(aes(decimal_date, euclidean_distance))+
  geom_line()

## biological env variables ----
env_data_interpolated_values_all_scaled |>
  colnames()

distances_env <- env_data_interpolated_values_all_scaled  |>
  dplyr::select(!c("temperature_no_nas"  ,  "day_length_no_nas"   ,  "chla_total_no_nas"    ,
                "PO4_no_nas"    ,        "NH4_no_nas"       ,     "NO2_no_nas"    ,        "NO3_no_nas"     ,       "Si_no_nas", sample_id_num)) |>
  stats::dist( method = "euclidean")

euclidean_distance <- distances_env |>
  as.matrix() |>
  as_tibble() |>
  bind_cols(decimal_date_tb_unique) |>
  pivot_longer(cols = -c('decimal_date', 'sample_id_num'), values_to = 'euclidean_distance', names_to = 'sample_id_num_2')

decimal_date_tb_unique <- decimal_date_tb_unique |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

euclidean_distance_tb_bio <- euclidean_distance |>
  dplyr::mutate(sample_distance = (as.numeric(sample_id_num_2) - as.numeric(sample_id_num))) |>
  dplyr::filter(sample_distance == 1) |>
  dplyr::select(-decimal_date) |>
  left_join(decimal_date_tb_unique, by = c('sample_id_num_2' = 'sample_id_num'))

euclidean_distance_tb_bio |>
  ggplot(aes(decimal_date, euclidean_distance))+
  geom_line()

# plot euclidean distance phyiscochemical and biological ----
euclidean_distance_tb_bio |>
  colnames()

euclidean_distance_tb_phch <- euclidean_distance_tb_phch |>
  dplyr::mutate(type = 'phch')

euclidean_distance_tb_g <- euclidean_distance_tb_g  |>
  dplyr::mutate(type = 'g') |>
  dplyr::mutate(sample_id_num = as.character(sample_id_num))

euclidean_distance_tb_bio |>
  dplyr::mutate(type = 'bio') |>
  bind_rows(euclidean_distance_tb_phch) |>
  ggplot(aes(decimal_date, euclidean_distance))+
  geom_line(aes(group = type))

euclidean_distance_tb_bio_phch <- euclidean_distance_tb_bio |>
  dplyr::mutate(type = 'bio') |>
  bind_rows(euclidean_distance_tb_phch) |>
  left_join(m_02_red2, by = c('sample_id_num_2' = 'sample_id_num')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(date, euclidean_distance, type)

eucl_plot <- euclidean_distance_tb_bio_phch |>
  ggplot(aes(date, euclidean_distance))+
  geom_line(data = euclidean_distance_tb_bio_phch |>
              dplyr::filter(type == 'phch'), aes(date, group = type, color = type), linewidth = 0.75, alpha = 0.3)+
  geom_line(aes(date, group = type, color = type, linetype = type), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Date', y = 'Euclidean Distance', color = '', linetype = '')+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 7),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    panel.grid.major.y = element_blank(),
    strip.placement = 'outside', 
    plot.margin = margin(t = 5, l = 20, r = 20))

eucl_plot
 
m_02 |>
  ggplot(aes(decimal_date, HNA))+
  geom_line()+
  geom_line(aes(y = LNA), linetype = 'dashed')

m_lna_hna <- m_02 |>
  dplyr::mutate(date = as.character.Date(date)) |>
  dplyr::select(date, LNA, HNA) |>
  pivot_longer(cols = !c('date')) |>
  add_row(date = '2003-12-10') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

plot_HNA_LNA <- m_lna_hna |>
  dplyr::filter(!is.na(name)) |>
  ggplot(aes(date, value))+
  scale_x_datetime(date_labels = '%Y', expand = c(0,0), date_breaks = '1 year', limits = c(min(m_lna_hna$date), max(m_lna_hna$date))
  )+
  geom_area(aes(date, value, fill = name, group = name), alpha = 1,  position='stack')+
  scale_y_continuous(labels = scientific_format(), expand = c(0,0))+
  scale_fill_manual(values = c('LNA' = '#9E364B', 'HNA' =  "#52343A"), na.value = "#000000")+
  scale_color_manual(values = c('LNA' = '#9E364B', 'HNA' =  "#52343A"), na.value = "#000000")+
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Family')+
  facet_grid(vars(name))+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         color = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'none', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

plot_HNA_LNA

# ggsave('plot_HNA_LNA.pdf', plot_HNA_LNA,
#        path = "results/figures/main_df3/supplementary/",
#        width = 180,
#        height = 120,
#        units = 'mm')

### some correlations (non significative)
## correlation PA FL
euclidean_distance_tb_bio_phch_w <- euclidean_distance_tb_bio_phch  |>
  pivot_wider(values_from = 'euclidean_distance', names_from = 'type')|>
  dplyr::mutate(date = as.Date.character(date)) |>
  dplyr::mutate(date = (as.character(date))) 

bray_unifrac_eucl_tb_02 |>
  left_join(euclidean_distance_tb_bio_phch_w, by = 'date') |>
  dplyr::filter(fraction == '0.2') |>
  ggplot(aes(bray_curtis_community, bio))+
  geom_point()+
  geom_smooth(method = 'lm')

bray_unifrac_eucl_tb_02 |>
  left_join(euclidean_distance_tb_bio_phch_w, by = 'date') |>
  dplyr::filter(fraction == '0.2') |>
  ggplot(aes(bray_curtis_community, phch))+
  geom_point()+
  geom_smooth(method = 'lm')

bray_unifrac_eucl_tb |>
  pivot_wider(values_from = 'bray_curtis_result', names_from = 'bray_curtis_type') |>
  dplyr::filter(fraction == '3') |>
  dplyr::mutate(date = as.Date.character(date)) |>
  dplyr::mutate(date = (as.character(date))) |>
  left_join(euclidean_distance_tb_bio_phch_w, by = 'date') |>
  dplyr::filter(fraction == '3') |>
  ggplot(aes(bray_curtis_community, bio))+
  geom_point()+
  geom_smooth(method = 'lm')

bray_unifrac_eucl_tb |>
  pivot_wider(values_from = 'bray_curtis_result', names_from = 'bray_curtis_type') |>
  dplyr::filter(fraction == '3') |>
  dplyr::mutate(date = as.Date.character(date)) |>
  dplyr::mutate(date = (as.character(date))) |>
  left_join(euclidean_distance_tb_bio_phch_w, by = 'date') |>
  dplyr::filter(fraction == '3') |>
  ggplot(aes(bray_curtis_community, phch))+
  geom_point()+
  geom_smooth(method = 'lm')

## NMDS env data
data.hel <- env_data_interpolated_values_all_scaled |>
  dplyr::select(-sample_id_num) |>
  dplyr::select("temperature_no_nas"  ,  "day_length_no_nas"   ,  "chla_total_no_nas"    ,
                "PO4_no_nas"    ,        "NH4_no_nas"       ,     "NO2_no_nas"    ,        "NO3_no_nas"     ,       "Si_no_nas" ) |>
  decostand(method="range"); str(data.hel)

data.dist <- vegdist(data.hel, method="euclidean")
head(data.dist)
data.nmds <- metaMDS(data.dist, k = 2)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.085 (per sota de 20; és acceptable)
data.nmds.points <- data.frame(data.nmds$points)  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)

data.nmds.points |>
  colnames()

data.nmds.points |>
  dim() ## 120 samples

decimal_date_tb_unique

nmds_10y_env <- data.nmds.points |>
  bind_cols(decimal_date_tb_unique) |>
  dplyr::select(-sample_id_num) |>
  left_join(m_02, by = c('decimal_date'))

nmds_10y_env$season <- factor(nmds_10y_env$season, levels = c('winter', 'spring', 'summer', 'autumn'))

nmds_10y_env |>
  ggplot(aes(MDS1, MDS2, color = season))+ # shape = fraction,
  geom_point(aes(shape = fraction))+
  scale_shape_discrete(labels = labs_fraction)+
  labs(shape = 'Fraction', color = 'Season')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = 'bottom', axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), strip.background = element_blank(), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
        strip.placement = 'outside')

nmds_10y_env |>
  ggplot(aes(MDS1, MDS2, color = as.factor(year)))+ # shape = fraction,
  geom_point(aes(shape = fraction))+
  scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_years)+
  labs(shape = 'Fraction', color = 'Season')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = 'bottom', axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), strip.background = element_blank(), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
        strip.placement = 'outside')

## with all my env variables ---
data.hel <- env_data_interpolated_values_all_scaled |>
  dplyr::select(-sample_id_num) |>
  decostand(method="range"); str(data.hel)

data.dist <- vegdist(data.hel, method="euclidean")
head(data.dist)
data.nmds <- metaMDS(data.dist, k = 2)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.085 (per sota de 20; és acceptable)
data.nmds.points <- data.frame(data.nmds$points)  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)

data.nmds.points |>
  colnames()

data.nmds.points |>
  dim() ## 120 samples

decimal_date_tb_unique

nmds_10y_env <- data.nmds.points |>
  bind_cols(decimal_date_tb_unique) |>
  dplyr::select(-sample_id_num) |>
  left_join(m_02, by = c('decimal_date'))

nmds_10y_env$season <- factor(nmds_10y_env$season, levels = c('winter', 'spring', 'summer', 'autumn'))

nmds_10y_env |>
  ggplot(aes(MDS1, MDS2, color = season))+ # shape = fraction,
  geom_point(aes(shape = fraction))+
  scale_shape_discrete(labels = labs_fraction)+
  labs(shape = 'Fraction', color = 'Season')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = 'bottom', axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), strip.background = element_blank(), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
        strip.placement = 'outside')

nmds_10y_env |>
  ggplot(aes(MDS1, MDS2, color = as.factor(year)))+ # shape = fraction,
  geom_point(aes(shape = fraction))+
  scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_years)+
  labs(shape = 'Fraction', color = 'Season')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = 'bottom', axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), strip.background = element_blank(), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
        strip.placement = 'outside')

#### Euclidean distance between samples ##### ----
### we need env data in wider format
env_data_interpolated_values_all |>
  colnames()

env_w <- env_data_interpolated_values_all |>
  dplyr::select(-sample_id_num, -BP_FC1.55_no_nas) |> 
  arrange(decimal_date) |>
  dplyr::mutate_if(is.numeric, scale)

env_w |>
  dim()

sample_num <- env_data_interpolated_values_all %$%
  sample_id_num |>
  as_tibble_col(column_name = 'sample_id_num')

# Dissimilarity matrix
distances_env <- env_w |>
  dplyr::select(-decimal_date) |>
  stats::dist( method = "euclidean")

distances_env_tb <- distances_env |>
  as.matrix() |>
  as_tibble() 

euclidean_distance <- distances_env_tb |>
  bind_cols(sample_num) |>
  pivot_longer(cols = -c('sample_id_num'), values_to = 'euclidean_distance', names_to = 'sample_num_2')

m <- m_02 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d/%m/%y")))

euclidean_distance_tb <- euclidean_distance |>
  dplyr::mutate(sample_distance = (as.numeric(sample_num_2) - as.numeric(sample_id_num))) |>
  dplyr::filter(sample_distance == 1) |>
  right_join(m, by = c( 'sample_num_2' = 'sample_id_num')) 

euclidean_distance_tb$season <- factor(euclidean_distance_tb$season, levels = c('winter', 'spring',
                                                                                'summer', 'autumn'))

euclidean_distance_plot <- euclidean_distance_tb |>  
  ggplot(aes(date, euclidean_distance))+
  geom_line()+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(y = 'Euclidean distance', x = 'Date')+
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7), 
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        aspect.ratio = 5/6.5)+
  guides(shape = 'none',
         color = guide_legend(ncol =1, size = 10,
                              override.aes = aes(label = ''))) 

euclidean_distance_plot

# ggsave(filename = 'euclidean_distance_plot.pdf',
#        plot = euclidean_distance_plot,
#        path = 'results/figures/',
#        width = 180, height = 80, units = 'mm')

## mean and sd prokaryotes abundances ---
m_bbmo_10y |>
  colnames()

m_bbmo_10y |>
  dplyr::select(date, bacteria_joint) |>
  distinct(date, bacteria_joint) |>
  dplyr::mutate(timeseries = 'BBMO') |>
  dplyr::group_by(timeseries) |>
  dplyr::reframe(mean = mean(bacteria_joint), sd = sd(bacteria_joint), se = sd(bacteria_joint)/sqrt(length(bacteria_joint))) |>
  dplyr::mutate(mean = scientific(mean), 
                sd = scientific(sd),
                se = scientific(se),
                cv = as.numeric(mean)/as.numeric(sd)*100)

### euclidean distances community -------
### we need env data in wider format
tb_w <- asv_tab_bbmo_10y_w_02_rclr  |>
  rownames_to_column(var = 'sample_id') |>
  dplyr::filter( str_detect(sample_id, '_0.2_')) 

tb_w |>
  dim()

sample_num <- tb_w  |>
  dplyr::mutate(sample_id_num = 1:nrow(tb_w))%$%
  sample_id_num |>
  as_tibble_col(column_name = 'sample_id_num')

# Dissimilarity matrix
distances_comm <- tb_w |>
  dplyr::select(-sample_id) |>
  stats::dist( method = "euclidean")

distances_comm_tb <- distances_comm |>
  as.matrix() |>
  as_tibble() 

euclidean_distance <- distances_comm_tb |>
  bind_cols(sample_num) |>
  pivot_longer(cols = -c('sample_id_num'), values_to = 'euclidean_distance', names_to = 'sample_num_2')

m <- m_02 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d/%m/%y")))

euclidean_distance_tb_02 <- euclidean_distance |>
  dplyr::mutate(sample_distance = (as.numeric(sample_num_2) - as.numeric(sample_id_num))) |>
  dplyr::filter(sample_distance == 1) |>
  right_join(m, by = c( 'sample_num_2' = 'sample_id_num'))   |>
  dplyr::mutate(fraction == '0.2')

euclidean_distance_tb$season <- factor(euclidean_distance_tb$season, levels = c('winter', 'spring',
                                                                                'summer', 'autumn'))

euclidean_distance_plot <- euclidean_distance_tb |>  
  ggplot(aes(date, euclidean_distance))+
  geom_line()+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(y = 'Euclidean distance', x = 'Date')+
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7), 
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        aspect.ratio = 5/6.5)+
  guides(shape = 'none',
         color = guide_legend(ncol =1, size = 10,
                              override.aes = aes(label = ''))) 

euclidean_distance_plot

#### 3-20 µm
tb_w <- asv_tab_bbmo_10y_w_3_rclr  |>
  rownames_to_column(var = 'sample_id') |>
  dplyr::filter( str_detect(sample_id, '_3_')) 

tb_w |>
  dim()

sample_num <- tb_w  |>
  dplyr::mutate(sample_id_num = 1:nrow(tb_w))%$%
  sample_id_num |>
  as_tibble_col(column_name = 'sample_id_num')

# Dissimilarity matrix
distances_comm <- tb_w |>
  dplyr::select(-sample_id) |>
  stats::dist( method = "euclidean")

distances_comm_tb <- distances_comm |>
  as.matrix() |>
  as_tibble() 

euclidean_distance <- distances_comm_tb |>
  bind_cols(sample_num) |>
  pivot_longer(cols = -c('sample_id_num'), values_to = 'euclidean_distance', names_to = 'sample_num_2')

m <- m_3 |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d/%m/%y")))

euclidean_distance_tb_all <- euclidean_distance |>
  dplyr::mutate(sample_distance = (as.numeric(sample_num_2) - as.numeric(sample_id_num))) |>
  dplyr::filter(sample_distance == 1) |>
  right_join(m, by = c( 'sample_num_2' = 'sample_id_num'))  |>
  dplyr::mutate(fraction == '3') |>
  bind_rows(euclidean_distance_tb_02)

euclidean_distance_tb_all$season <- factor(euclidean_distance_tb_all$season, levels = c('winter', 'spring',
                                                                                'summer', 'autumn'))

euclidean_distance_plot <- euclidean_distance_tb_all |>  
  ggplot(aes(date, euclidean_distance))+
  geom_line(aes(date, group = fraction, color = fraction, linetype = fraction), linewidth = 0.75, alpha = 1)+ #, linetype = bray_curtis_type
  scale_color_manual(values= palette_fraction_env, labels = labs_fraction_env)+ #, labels = labs_fraction
  scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(data = euclidean_distance_tb_all |>
              dplyr::filter(fraction == '3'), aes(date, euclidean_distance), alpha = 0.3)+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(y = 'Euclidean distance', x = 'Date')+
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7), 
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        aspect.ratio = 5/6.5)+
  guides(shape = 'none',
         color = guide_legend(ncol =1, size = 10,
                              override.aes = aes(label = ''))) 

euclidean_distance_plot

euclidean_distance_tb_bio_phch |>
  pivot_wider(values_from = euclidean_distance, names_from = type) |>
  dplyr::select(euclidean_distance_env_bio = bio, euclidean_distance_env_phch = phch, date) |>
  left_join(euclidean_distance_tb_all) |>
  ggplot(aes(euclidean_distance_env_bio, euclidean_distance))+
  geom_point(aes(shape = fraction))+
  facet_wrap(vars(fraction))+
  geom_smooth(method = 'lm')+
  stat_cor( aes(
                                               label =   paste(..p.label..)), 
            label.x = 5,  
           label.y = 15,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "black"#,
           #position = position_jitter(0.0)
  )+
  stat_cor(aes(
                                               label =   paste(..r.label..)),
           label.x = 5, label.y = 10,  color = "black" ,
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  theme_bw()

euclidean_distance_tb_bio_phch |>
  pivot_wider(values_from = euclidean_distance, names_from = type) |>
  dplyr::select(euclidean_distance_env_bio = bio, euclidean_distance_env_phch = phch, date) |>
  left_join(euclidean_distance_tb_all) |>
  dplyr::mutate(  euclidean_distance = as.numeric(scale(euclidean_distance)),
                  euclidean_distance_tb_bio = as.numeric(scale(euclidean_distance_env_bio)),
                  euclidean_distance_phch = as.numeric(scale(euclidean_distance_env_phch))) |>
  ggplot(aes(euclidean_distance_phch , euclidean_distance))+
  geom_point(aes(shape = fraction))+
  facet_wrap(vars(fraction))+
  geom_smooth(method = 'lm')+
  stat_cor( aes(
    label =   paste(..p.label..)), 
    label.x = 5,  
    label.y = 15,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = "black"#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes(
    label =   paste(..r.label..)),
    label.x = 5, label.y = 10,  color = "black" ,
    p.digits = 0.01, digits = 2, 
    p.accuracy = 0.01, method = 'spearman')+
  theme_bw()
  
