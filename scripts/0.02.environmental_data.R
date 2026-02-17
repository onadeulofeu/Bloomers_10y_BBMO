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

##upload data----
asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv')
asv_tab_rar <- read.csv2('data/asv_tab_bbmo_10y_w_rar.csv') |>
  as_tibble()

library(readxl) ##I upload the general metadata from the whole BBMO 20Y
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
    # mutate(type_of_env = case_when(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,      
    #                                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", 
    #                                                              "bacteria_joint", "LNA", "HNA") ~ 'physico_chem',
    #                                environmental_variable %in% c(  "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
    #                                       "HNF2_5um_Micro", "HNF_5um_Micro", "prochlorococcus_FC", "Peuk1",             
    #                                       "Peuk2", "synechococcus") ~ 'biological')) |>
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
# palete_gradient_cb <- c("#86a0df",
#                      #"#666585" = 0,
#                      '#bbbbbb' = 0,
#                      "#b81131") 

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
  #scale_fill_gradientn( colours = palete_gradient_cb2)+
  #scale_fill_gradient2(low = palete_gradient_cb2[1], high = palete_gradient_cb2[4], midpoint = 0) +
  scale_fill_gradientn(colors = palete_gradient_cb5,
                       breaks=c(-5, 0, 5),
                       limits=c(-5.9,  5.9)) +
  # scale_fill_gradient2(low = "#0049B7",
  #                      high = "#b81131", midpoint = 0)+
  facet_wrap(.~year, scales = 'free_x', nrow = 1, switch = 'x')+
  geom_tile(alpha = 1)+
  geom_vline(xintercept = seq(0.5, 12, by = 12), linetype = "dashed", color = "darkgrey") +
  #scale_x_datetime(expand = c(0,0))+
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
  #scale_fill_gradientn( colours = palete_gradient_cb)+
  # scale_fill_gradient2(low = "#0049B7",
  #                      high = "#b81131", midpoint = 0)+
  # scale_fill_gradient2(colours = palete_gradient_cb2,
  #                      midpoint = 0  # Set the midpoint value)+
    scale_fill_gradientn(colors = palete_gradient_cb5,
                         breaks=c(-6, 0, 6),
                         limits=c(-6,  6), na.value="#4d009c")+ # I add this since I have one value which is outside the range, to define the color it should have
  facet_wrap(.~year, scales = 'free_x', nrow = 1, switch = 'x')+
  geom_tile(alpha = 1)+
  geom_vline(xintercept = seq(0.5, 12, by = 12), linetype = "dashed", color = "darkgrey") +
  #scale_x_datetime(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), labels = unique(bbmo_env_z_ed$year), breaks = unique(bbmo_env_z_ed$year))+
  labs(x = 'Year', y = '', fill = 'z-score')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 5), legend.position = 'bottom',
        panel.border = element_blank(), strip.background = element_blank(), plot.margin = unit(c(0,3,0,1), "mm"), 
        legend.key.size = unit(3, "mm"))

ggsave('biological_bbmo.pdf', biological,
       path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/",
       width = 188,
       height = 90,
       units = 'mm')

## I divide them in different groups so that we can observe them better----
bbmo_env_z |>
  left_join(m_bbmo_10y) |>
  dplyr::filter(environmental_variable %in% c('day_length', 'temperature', 'secchi', 'salinity', 'chla_total', 'chla_3um', 'PO4', 'NH4',
                                              'NO2', 'NO3', "Si", "BP_FC1.55" )) |>
  # mutate(type_of_env = case_when(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,      
  #                                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", 
  #                                                              "bacteria_joint", "LNA", "HNA") ~ 'physico_chem',
  #                                environmental_variable %in% c(  "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
  #                                       "HNF2_5um_Micro", "HNF_5um_Micro", "prochlorococcus_FC", "Peuk1",             
  #                                       "Peuk2", "synechococcus") ~ 'biological')) |>
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
  # mutate(type_of_env = case_when(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,      
  #                                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", 
  #                                                              "bacteria_joint", "LNA", "HNA") ~ 'physico_chem',
  #                                environmental_variable %in% c(  "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
  #                                       "HNF2_5um_Micro", "HNF_5um_Micro", "prochlorococcus_FC", "Peuk1",             
  #                                       "Peuk2", "synechococcus") ~ 'biological')) |>
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

#plot the original data without the zscores normalization----
env_data_physico_chem <- bbmo_env_z |>
left_join(m_bbmo_10y) |>
  dplyr::filter(environmental_variable %in% c('day_length', 'temperature', 'secchi', 'salinity', 'chla_total', 'chla_3um', 'PO4', 'NH4',
                                              'NO2', 'NO3', "Si", "BP_FC1.55" )) |>
  # mutate(type_of_env = case_when(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,      
  #                                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", 
  #                                                              "bacteria_joint", "LNA", "HNA") ~ 'physico_chem',
  #                                environmental_variable %in% c(  "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
  #                                       "HNF2_5um_Micro", "HNF_5um_Micro", "prochlorococcus_FC", "Peuk1",             
  #                                       "Peuk2", "synechococcus") ~ 'biological')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, env_values))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_point(alpha = 0.6, size = 1/2)+
  geom_line()+
  labs(x = 'Time', y = '')+
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

ggsave('env_data_physico_chem.pdf', env_data_physico_chem,
       path = '~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/',
       width = 230,
       height = 180,
       units = 'mm')

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
  # mutate(type_of_env = case_when(environmental_variable %in% c("day_length", "temperature" ,"secchi" , "salinity" ,      
  #                                                              "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", 
  #                                                              "bacteria_joint", "LNA", "HNA") ~ 'physico_chem',
  #                                environmental_variable %in% c(  "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
  #                                       "HNF2_5um_Micro", "HNF_5um_Micro", "prochlorococcus_FC", "Peuk1",             
  #                                       "Peuk2", "synechococcus") ~ 'biological')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, env_values))+
  geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                    ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  geom_point(alpha = 0.6, size = 1/2)+
  geom_line()+
  labs(x = 'Time', y = 'Abundances (cells/mL) ')+
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

ggsave('env_data_cyto.pdf', env_data_cyto,
       path = '~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/',
       width = 230,
               height = 180,
               units = 'mm')

##calculate distance matrices----
asv_tab_bbmo_10y_l |> # upload the original table
  head()

####community data distance-----
## If we transform data to z-scores or rCLR we need to use Euclidean distance because we will have negative values, if
## we don't transform our data, then we can use Bray curtis distance

#####Rarefied dataset---
  abund_rar_02 <- asv_tab_rar |>
    dplyr::filter(str_detect(X, '_0.2_'))|>
    as.data.frame() #|>
    # dplyr::select(-X) |>
    # vegdist(method = 'bray')
  
  abund_rar_3 <- asv_tab_rar |>
    dplyr::filter(str_detect(X, '_3_'))|>
    as.data.frame() #|>
    # dplyr::select(-X) |>
    # vegdist(method = 'bray')
  
  # asv_tab_bbmo_10y_l |> # upload the original table
  #   head()
  
  asv_tab_bbmo_10y_w <- asv_tab_bbmo_10y_l |> #transform to wider format
    pivot_wider(id_cols = sample_id, names_from = asv_num, values_from = reads) |>
    as.data.frame()
  
  #####Transformed to rCLR----
  asv_tab_bbmo_10y_w_02 <- asv_tab_bbmo_10y_w |>
    dplyr::filter(str_detect(as.character(sample_id), '_0.2'))
  
  asv_tab_bbmo_10y_w_3 <- asv_tab_bbmo_10y_w |>
    dplyr::filter(str_detect(as.character(sample_id), '_3_'))
  
  row.names(abund_rar_02) <- abund_rar_02[,1] 
  row.names(abund_rar_3) <- abund_rar_3[,1] 

  row.names(asv_tab_bbmo_10y_w_02) <- asv_tab_bbmo_10y_w_02[,1]  
  row.names(asv_tab_bbmo_10y_w_3) <- asv_tab_bbmo_10y_w_3[,1] 
  
  abund_rar_02_ed <- abund_rar_02[,-1]
  abund_rar_3_ed <- abund_rar_3[,-1]
  
  asv_tab_bbmo_10y_w_02_ed <- asv_tab_bbmo_10y_w_02[,-1]
  asv_tab_bbmo_10y_w_3_ed <- asv_tab_bbmo_10y_w_3[,-1]
  
  abund_rar_02_rclr <-  abund_rar_02_ed |>
    decostand(method="rclr")
  abund_rar_3_rclr <-  abund_rar_3_ed |>
    decostand(method="rclr")
  
  asv_tab_bbmo_10y_w_02_rclr <-  asv_tab_bbmo_10y_w_02_ed |>
    decostand(method="rclr")
  asv_tab_bbmo_10y_w_3_rclr <-  asv_tab_bbmo_10y_w_3_ed |>
    decostand(method="rclr")
  
  # data.hel <- asv_tab_bbmo_10y_w_02_ed |>
  #   decostand(method="rclr"); str(data.hel)
  # 
  # data.hel |>
  #   dim()
  # 
  # data.dist_02 <- vegdist(data.hel, method="euclidean", na.rm = TRUE)
  # head(data.dist)
  # 
  # data.dist_02 |>
  #   class()
  
  ####env variables-----
# bbmo_env_02_dist <- bbmo_env |>
#   dplyr::filter(str_detect(sample_id, '0.2')) |>
#   dist(method = 'euclidean')
# 
# bbmo_env_3_dist <- bbmo_env |>
#   dplyr::filter(str_detect(...28, '3')) |>
#   dist(method = 'euclidean')

## here we use environmental data normalized using z_scores
# bbmo_envz_02_dist <- bbmo_env_zscore_w |>
#   dplyr::filter(str_detect(sample_id, '_0.2')) |>
#   dist(method = 'euclidean')
# 
# bbmo_envz_3_dist <- bbmo_env_zscore_w |>
#   dplyr::filter(str_detect(sample_id, '_3')) |>
#   dist(method = 'euclidean')

#Mantel test----
## The test statistic is the correlation coefficient. r falls in the range
## of -1 to +1, where being close to -1 indicates strong negative correlation
## and +1 indicates strong positive correlation. An r value of 0 indicates no
## correlation.
## 02
# abund_dist_02 |>
#   class()
# 
# abund_env_02 <- mantel(abund_dist_02, bbmo_env_02_dist,
#                        method = 'spearman', permutations = 9999,
#                        na.rm = TRUE)
# 
# #with env data normalized with z_scores
# abund_envz_02 <- mantel(abund_dist_02, bbmo_envz_02_dist,
#                        method = 'spearman', permutations = 9999,
#                        na.rm = TRUE)
# 
# ##3
# abund_env_3 <- mantel(abund_dist_3, bbmo_env_3_dist,
#                        method = 'spearman', 
#                       permutations = 9999,
#                        na.rm = TRUE)

## MANTEL TEST environmental variable (z-scores) vs community structure (transformed to zCLR)----

# mantel(data.dist_02, bbmo_envz_02_dist,
#        method = 'spearman', 
#        permutations = 9999,
#        na.rm = TRUE)

##one environmental variable at a time
  
##inputs FL (02) ----
  abund_rar_02_rclr
  bbmo_envz_02 <- bbmo_env_zscore_w |>
    dplyr::filter(str_detect(sample_id, '_0.2_')) |>
    as.data.frame()

  row.names(bbmo_envz_02) <-  bbmo_envz_02[,1] 
  bbmo_envz_02_ed <-  bbmo_envz_02[,-1]
# Extract common sample names
common_sample_names <- intersect(rownames( abund_rar_02_rclr), rownames(bbmo_envz_02_ed))

# Subset both matrices based on common sample names
abundance_data <- abund_rar_02_rclr[common_sample_names, ]
environmental_data <- bbmo_envz_02_ed[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_02_ed[, colnames(bbmo_envz_02_ed) != "other_variables"]  # Adjust as needed

# check that the names match
row.names(abund_rar_02_rclr) == row.names(bbmo_envz_02_ed)

# Initialize a list to store Mantel test results
mantel_results_02 <- list()

# Iterate over each environmental variable
for (variable in colnames(environmental_variables)) {
  
  # Extract the current environmental variable
  current_variable <- environmental_variables[, variable]
  
  # Calculate Euclidean distance
  euclidean_distance <- vegdist(current_variable, method = "euclidean", 
                                na.rm = TRUE)
  
  # Assuming you have another matrix or data frame for the second set of variables
  data.dist_02 <- vegdist(abund_rar_02_rclr, method="euclidean", 
                          na.rm = TRUE)
  # Replace 'your_second_set_of_variables' with your actual data
  # second_set_of_variables <-  asv_tab_bbmo_10y_w_02_ed |>
  #   dplyr::select(starts_with('asv'))

  # Calculate Euclidean distance for the second set of variables
  # euclidean_distance_second_set <- dist(second_set_of_variables, method = "euclidean")
  
  # Perform Mantel test
  mantel_test_result_02 <- mantel(euclidean_distance, data.dist_02, method = "pearson", permutations = 999,
                               na.rm = TRUE)
  
  # Store the Mantel test result in the list
  mantel_results_02[[variable]] <- mantel_test_result_02
}

# 'mantel_results_02' now contains Mantel test results for each environmental variable
# You can access the results using mantel_results$variable_name
# For example, mantel_results$"your_variable_name"$

# mantel_results_02$day_length
# mantel_results_02$temperature
# mantel_results_02$secchi
# mantel_results_02$salinity

# Create an empty tibble to store the summary results
results_mantel_02 <- tibble(
  variable_name = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each variable in 'mantel_results' and add the results to the tibble
for (variable_name in names(mantel_results_02)) {
  results <- mantel_results_02[[variable_name]]
  
    mantel_correlation <- results$statistic

  results_mantel_02 <- bind_rows(results_mantel_02, tibble(
    variable_name = variable_name,
    mantel_correlation = mantel_correlation,
    p_value = results$signif
  ))
}

# Print the summary tibble
print(results_mantel_02)

#write.csv2(results_mantel_02, 'results/mantel_community_02.csv')

#plot the results for the Mantel test between environmental variables and the community structure 
results_mantel_02$variable_name <- results_mantel_02$variable_name  |>
  factor(levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                    "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                    "bacteria_joint", "LNA", "HNA", "Peuk1",             
                    "Peuk2", 'total_vlp', 'low_vlp', 'med_vlp', 'high_vlp',
                    "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                    "HNF2_5um_Micro", "HNF_5um_Micro"))


plot_mantel_02_community <- results_mantel_02 |>
  mutate(community = 'Free living (0.2-3 um)') |>
  #left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  ggplot(aes(community, variable_name, fill = mantel_correlation))+
  geom_tile()+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  scale_y_discrete(labels = labs_env)+
  geom_text(aes(label = ifelse(p_value < 0.05, '*', '')))+
  labs(x = 'Community', y = 'Environmental variables', fill = 'Mantel correlation', title = 'Particle Attached')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), text = element_text(size = 5),
        axis.text.x = element_text(angle = 65, hjust = 1))

## inputs second PA (3) ----
abund_rar_3_rclr |>
  row.names()
bbmo_envz_3 <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  as.data.frame()

row.names(bbmo_envz_3) <-  bbmo_envz_3[,1] 
bbmo_envz_3_ed <-  bbmo_envz_3[,-1]
# Extract common sample names
common_sample_names <- intersect(rownames( abund_rar_3_rclr), rownames(bbmo_envz_3_ed))

# Subset both matrices based on common sample names
abundance_data <- abund_rar_3_rclr[common_sample_names, ]
environmental_data <- bbmo_envz_3_ed[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_3_ed[, colnames(bbmo_envz_3_ed) != "other_variables"]  # Adjust as needed

# check that the names match
row.names(abund_rar_3_rclr) == row.names(bbmo_envz_3_ed)

# Initialize a list to store Mantel test results
mantel_results_3 <- list()

# Iterate over each environmental variable
for (variable in colnames(environmental_variables)) {
  
  # Extract the current environmental variable
  current_variable <- environmental_variables[, variable]
  
  # Calculate Euclidean distance
  euclidean_distance <- vegdist(current_variable, method = "euclidean", 
                                na.rm = TRUE)
  
  # Assuming you have another matrix or data frame for the second set of variables
  data.dist_3 <- vegdist(abund_rar_3_rclr, method="euclidean", 
                          na.rm = TRUE)
  # Replace 'your_second_set_of_variables' with your actual data
  # second_set_of_variables <-  asv_tab_bbmo_10y_w_3_ed |>
  #   dplyr::select(starts_with('asv'))
  
  # Calculate Euclidean distance for the second set of variables
  # euclidean_distance_second_set <- dist(second_set_of_variables, method = "euclidean")
  
  # Perform Mantel test
  mantel_test_result_3 <- mantel(euclidean_distance, data.dist_3, method = "pearson", permutations = 999,
                                  na.rm = TRUE)
  
  # Store the Mantel test result in the list
  mantel_results_3[[variable]] <- mantel_test_result_3
}

# 'mantel_results_3' now contains Mantel test results for each environmental variable
# You can access the results using mantel_results$variable_name
# For example, mantel_results$"your_variable_name"$

# mantel_results_3$day_length
# mantel_results_3$temperature
# mantel_results_3$secchi
# mantel_results_3$salinity

# Create an empty tibble to store the summary results
results_mantel_3 <- tibble(
  variable_name = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each variable in 'mantel_results' and add the results to the tibble
for (variable_name in names(mantel_results_3)) {
  results <- mantel_results_3[[variable_name]]
  
  mantel_correlation <- results$statistic
  
  results_mantel_3 <- bind_rows(results_mantel_3, tibble(
    variable_name = variable_name,
    mantel_correlation = mantel_correlation,
    p_value = results$signif
  ))
}

# Print the summary tibble
print(results_mantel_3)

#plot the results for the Mantel test between environmental variables and the community structure 
results_mantel_3$variable_name <- results_mantel_3$variable_name  |>
  factor(levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                    "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                    "bacteria_joint", "LNA", "HNA", "Peuk1",             
                    "Peuk2", 'total_vlp', 'low_vlp', 'med_vlp', 'high_vlp',
                    "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                    "HNF2_5um_Micro", "HNF_5um_Micro"))

results_mantel_02 <-   results_mantel_02 |>
  dplyr::mutate(community = 'Free living (0.2-3 um)')

#write.csv2(results_mantel_3, 'results/mantel_community_3.csv')

##plot results from FL and PA together-----
plot_mantel_3_02_community <- results_mantel_3 |>
  dplyr::mutate(community = 'Particle attached (3-20 um)') |>
  bind_rows(results_mantel_02) |>
  dplyr::filter(p_value < 0.05 &
                  abs(mantel_correlation) > 0.1) |>
  #left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  ggplot(aes(community, variable_name, fill = mantel_correlation))+
  geom_tile()+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  scale_y_discrete(labels = labs_env)+
  #geom_text(aes(label = ifelse(p_value < 0.05, '*', ''), colour = '#757C76'))+
  labs(x = 'Community', y = 'Environmental variables', fill = 'Mantel correlation')+
  guides( colour = 'none')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), text = element_text(size = 5),
        axis.text.x = element_text(angle = 0), axis.ticks.x = element_blank(), panel.border = element_blank())

# ggsave('mantel_test_community_sig.pdf', plot_mantel_3_02_community, path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/results/figures/",
#        width = 88,
#        height = 88,
#        units = 'mm')

## I perform the Mantel test only with my potential bloomers vs each environmental variable----

### Prepare the inputs for the Mantel test
asv_anom_3_tb <- asv_anom_3 |>
  as_tibble()

asv_anom_02_tb <- asv_anom_02 |>
  as_tibble()

abundance_data_bloo02 <- asv_tab_bbmo_10y_w_02_ed |>
  dplyr::select(asv_anom_02_tb$value) # filter by potential bloomers in FL

abundance_data_bloo3 <- asv_tab_bbmo_10y_w_3_ed |>
  dplyr::select(asv_anom_3_tb$value) # filter by potential bloomers in PA

#abundance data transformed to rCLR: 
abundance_data_clr02 <- abundance_data_bloo02 |>
  decostand(method="clr", na.rm = TRUE, pseudocount = 1) 

abundance_data_clr02 |>
  dim()

abundance_data_clr3 <- abundance_data_bloo3 |>
  decostand(method="clr", na.rm = TRUE, pseudocount = 1) 

abundance_data_clr3 |>
  dim()

###FL -----
bbmo_envz_02 <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  as.data.frame()

row.names(bbmo_envz_02) <-  bbmo_envz_02[,1] 
bbmo_envz_02_ed <-  bbmo_envz_02[,-1]

# Extract common sample names
common_sample_names <- intersect(rownames(abundance_data_clr02), rownames(bbmo_envz_02_ed))

# Subset both matrices based on common sample names
abundance_data <- abundance_data_clr02[common_sample_names, ]
environmental_data <- abundance_data_clr02[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_02_ed[, colnames(bbmo_envz_02_ed) != "other_variables"]  # Adjust as needed

# check that the names match
row.names(abundance_data_clr02) == row.names(bbmo_envz_02_ed)

# Create an empty tibble to store the Mantel test results
mantel_results_02 <- tibble(
  taxon_name = character(),
  variable_name_environmental = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each taxon in the abundance matrix
for (taxon_name in colnames(abundance_data_clr02)) {
  
  # Extract the abundance data for the current taxon
  current_variable_abundance <- abundance_data_clr02[, taxon_name, drop = FALSE]
  
  # Calculate distance for the abundance variable
  distance_abundance <- vegdist(current_variable_abundance, method = "euclidean", na.rm = TRUE)
  
  # Iterate over each environmental variable
  for (variable_name_environmental in colnames(bbmo_envz_02_ed)) {
    
    # Extract the environmental data for the current variable
    current_variable_environmental <- bbmo_envz_02_ed[, variable_name_environmental, drop = FALSE]
    
    # Calculate distance for the environmental variable
    distance_environmental <- vegdist(current_variable_environmental, method = "euclidean", na.rm = TRUE)
    
    # Perform Mantel test
    mantel_test_result <- mantel(distance_abundance, distance_environmental, method = "pearson", permutations = 999, na.rm = TRUE)
    
    # Store the Mantel test result in the tibble
    mantel_results_02 <- bind_rows(mantel_results_02, tibble(
      taxon_name = taxon_name,
      variable_name_environmental = variable_name_environmental,
      mantel_correlation = mantel_test_result$statistic,
      p_value = mantel_test_result$signif
    ))
  }
}

# Print the summary tibble
print(mantel_results_02)

#write.csv2(mantel_results_02, 'results/mantel_results_02_bloo.csv')

## reorder environmental factors ----
mantel_results_02$variable_name_environmental <- mantel_results_02$variable_name_environmental  |>
  factor(levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                    "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                    "bacteria_joint", "LNA", "HNA", "Peuk1",             
                    "Peuk2", 'total_vlp', 'low_vlp', 'med_vlp', 'high_vlp',
                    "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                    "HNF2_5um_Micro", "HNF_5um_Micro"))

# plot the results
plot_mantel_02 <- mantel_results_02 |>
  left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  dplyr::mutate(community = 'Free living') |>
  ggplot(aes(interaction(taxon_name, family_f), variable_name_environmental, fill = mantel_correlation))+
  geom_tile()+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  scale_y_discrete(labels = labs_env)+
  facet_wrap(vars(community))+
  labs(x = 'Taxonomy', y = 'Environmental variables', fill = 'Mantel correlation')+ #, title = 'Free living fraction'
  geom_text(aes(label = ifelse(p_value < 0.05, '*', ''), colour = '#757C76'))+
  guides( colour = 'none')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), text = element_text(size = 5),
        axis.text.x = element_text(angle = 55, hjust = 1), strip.background = element_blank(),
        panel.border = element_blank())

####for the PA fraction------
bbmo_envz_3 <- bbmo_env_zscore_w |>
  dplyr::filter(str_detect(sample_id, '_3_')) |>
  as.data.frame()

row.names(bbmo_envz_3) <-  bbmo_envz_3[,1] 
bbmo_envz_3_ed <-  bbmo_envz_3[,-1]

# Extract common sample names
common_sample_names <- intersect(rownames(abundance_data_clr3), rownames(bbmo_envz_3_ed))

# Subset both matrices based on common sample names
abundance_data <- abundance_data_clr3[common_sample_names, ]
environmental_data <- abundance_data_clr3[common_sample_names, ]

# Extract the environmental variables (assuming they are in columns)
environmental_variables <- bbmo_envz_3_ed[, colnames(bbmo_envz_3_ed) != "other_variables"]  # Adjust as needed

# check that the names match
row.names(abundance_data_clr3) == row.names(bbmo_envz_3_ed)

# Create an empty tibble to store the Mantel test results
mantel_results_3 <- tibble(
  taxon_name = character(),
  variable_name_environmental = character(),
  mantel_correlation = double(),
  p_value = double()
)

# Iterate over each taxon in the abundance matrix
for (taxon_name in colnames(abundance_data_clr3)) {
  
  # Extract the abundance data for the current taxon
  current_variable_abundance <- abundance_data_clr3[, taxon_name, drop = FALSE]
  
  # Calculate distance for the abundance variable
  distance_abundance <- vegdist(current_variable_abundance, method = "euclidean", na.rm = TRUE)
  
  # Iterate over each environmental variable
  for (variable_name_environmental in colnames(bbmo_envz_3_ed)) {
    
    # Extract the environmental data for the current variable
    current_variable_environmental <- bbmo_envz_3_ed[, variable_name_environmental, drop = FALSE]
    
    # Calculate distance for the environmental variable
    distance_environmental <- vegdist(current_variable_environmental, method = "euclidean", na.rm = TRUE)
    
    # Perform Mantel test
    mantel_test_result <- mantel(distance_abundance, distance_environmental, method = "pearson", permutations = 999, na.rm = TRUE)
    
    # Store the Mantel test result in the tibble
    mantel_results_3 <- bind_rows(mantel_results_3, tibble(
      taxon_name = taxon_name,
      variable_name_environmental = variable_name_environmental,
      mantel_correlation = mantel_test_result$statistic,
      p_value = mantel_test_result$signif
    ))
  }
}

# Print the summary tibble
print(mantel_results_3)

#write.csv2(mantel_results_3, 'results/mantel_results_3_bloo.csv')

##reorder taxonomy as factors 
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num))

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

# Plot the results
tax_factors <- asv_tab_all_bloo_z_tax |>
  dplyr::select(asv_num_f, family_f, order_f, class_f, phylum_f) |>
  distinct(asv_num_f, family_f, order_f, class_f, phylum_f)

asv_tab_all_bloo_z_tax |>
  colnames()

palete_gradient_cb <- c(#"#240023",
  
  "#4db2a2" = 0,
  "#005a47",
  na.value = '#000000') 

## reorder environmental factors ----
mantel_results_3$variable_name_environmental <- mantel_results_3$variable_name_environmental  |>
  factor(levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                    "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                    "bacteria_joint", "LNA", "HNA", "Peuk1",             
                    "Peuk2", 'total_vlp', 'low_vlp', 'med_vlp', 'high_vlp',
                    "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "dryptomonas", "micromonas","HNF_Micro",         
                    "HNF2_5um_Micro", "HNF_5um_Micro"))

## add FL results to plot them together
mantel_results_02_tax <- mantel_results_02 |>
  left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  dplyr::mutate(community = 'Free living (0.2-3 um)')


plot_mantel_all <- mantel_results_3 |>
  left_join(tax_factors, by = c('taxon_name' = 'asv_num_f')) |>
  dplyr::mutate(community = 'Particle Attached (3-20 um)') |>
  bind_rows(mantel_results_02_tax) |>
  dplyr::filter(p_value < 0.05 &
                  abs(mantel_correlation) > 0.2) |>
  ggplot(aes(interaction(taxon_name, family_f), variable_name_environmental, fill = mantel_correlation))+
  geom_tile()+
  scale_fill_gradientn(colors = palete_gradient_cb)+
  scale_y_discrete(labels = labs_env)+
  #coord_flip()+
  facet_wrap(vars(community), nrow = 2)+
  labs(x = 'Taxonomy', y = 'Environmental variables', fill = 'Mantel\ncorrelation')+ #, title = 'Free living fraction'
  #geom_text(aes(label = ifelse(p_value < 0.05, '*', ''), colour = '#757C76'))+
  guides(colour = 'none')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), text = element_text(size = 5),
        axis.text.x = element_text(angle = 55, hjust = 1), strip.background = element_blank(),
        panel.border = element_blank(), legend.key.size = unit(4, 'mm'), plot.margin = margin(0,0,5,15),
        legend.margin = margin(0,15,0,5))

# mantel_test_all <-  multi_panel_figure(columns = 1, rows = 2, width = 180, height = 200, 
#                                                row_spacing = 0.2, unit = 'mm',
#                                                panel_label_type = 'upper-alpha')
# 
# mantel_test_all  %<>%
#   fill_panel(plot_mantel_02, column = 1, row = 1) %<>%
#   fill_panel(plot_mantel_3, column = 1, row = 2)
# 
# ggsave('mantel_test_all_bloo_sig.pdf', plot_mantel_all ,
#        path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/results/figures/",
#        width = 88,
#        height = 100,
#        units = 'mm')

## TEST: plot the results of the Mantel test-------
# correlation_coefficient <- abund_env_02$statistic
# p_value <- abund_env_02$signif
# 
# results_mantel_02 <- tibble(Mantel_correlation = correlation_coefficient,
#                             p_value = p_value)
# 
# results_mantel_3 |>
#   ggplot(aes(Mantel_correlation, p_value))+
#   geom_point()
# 
# 
# # Option, correlation matrix with bloomers and environmental data
# ## subset ASV tab by potential bloomers
# asv_bloo <- asv_tab_all_bloo_z_tax |>
#   dplyr::select(asv_num) |>
#   as_vector() |>
#   unique() 
# 
# asv_tab_all_bloo_z_tax |>
#   colnames()
# 
# abund_dist_02_bloo <- asv_tab_rar |>
#   dplyr::select(X, all_of(asv_bloo)) |>
#   dplyr::filter(str_detect(X, '0.2')) |>
#   dplyr::select(-X) |>
#   vegdist(method = 'bray')
# 
# abund_env_02 <- mantel(abund_dist_02_bloo, bbmo_env_02_dist,
#                        method = 'spearman', permutations = 9999,
#                        na.rm = TRUE)
# library(corrplot)
# qcorrplot(correlate(bbmo_env))
# 
# 
# ##visualization of correlation matrix
# corr_bloo_02 <- cor(abund_dist_02_bloo, bbmo_envz_02_dist)
# 
# ggcorrplot::ggcorrplot(as.matrix(corr_bloo_02))
# ggcorrplot::ggcorrplot(as.matrix(abund_dist_02))
# ggcorrplot::ggcorrplot(as.matrix(abund_dist_3))
# ggcorrplot::ggcorrplot(as.matrix(bbmo_env_02_dist))
# ggcorrplot::ggcorrplot(as.matrix(bbmo_env_3_dist))
# 
# ggcorrplot::ggcorrplot(as.matrix(bbmo_envz_02_dist))
# ggcorrplot::ggcorrplot(as.matrix(bbmo_envz_3_dist))
# 
# ##
# rda1 <- rda(abund_dist_02, bbmo_env3)
# 
# ## Mantel correlogram
# ## This function comptues a multivariate Mantel correlogram.
# mantel.correlog(abund_dist_02, bbmo_env_02_dist,
#                 nperm = 9999)
# 
# mantel.correlog(abund_dist_02, bbmo_envz_02_dist,
#                 nperm = 9999)
# 
# |>
#   print.mantel.correlog()
# 
# bbmo_env_02_dist |>
#   head()
# 
# bbmo_env_02_dist |>
#   dim()
# 
# bbmo_env |>
#   dim()

## other ideas
### perform a mantel test for each environmental variable and the community 


### perform a mantel test for each environmental variable and each potential bloomer


## Visualizing the data-----
asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  distinct(asv_num)

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
                fraction == 3 &
                  abundance_type == 'relative_abundance') |>
  ggplot(aes(temperature, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(temperature, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(secchi, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(salinity, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(chla_total, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(chla_3um, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(PO4, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(NH4, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(NO2, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(NO3, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  ggplot(aes(dryptomonas, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

## try the same plots but with z_scores
asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(dryptomonas_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(temperature_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(day_length_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(secchi_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(salinity_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(chla_total_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(#asv_num == 'asv179' &
    fraction == 3 &
      abundance_type == 'relative_abundance') |>
  left_join(bbmo_env_zscore_w, suffix = c('no_norm', '_z_norm'), by = c('sample_id')) |>
  ggplot(aes(PO4_z_norm, abundance_value, group = asv_num, color = family))+
  geom_point()+
  geom_smooth(aes(group = asv_num), se = F)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  theme_bw()

## Multipatt function:
### This function studies the association between species patterns and combinations of groups of sites.
library(indicspecies)

multipatt(asv_tab_rar)


## multicorrelation
library(psych)
pairs.panels(bbmo_env[2:13], method = 'pearson')

library(ggcorrplot)
cor_pmat(asv_tab_rar)

ggcorr(asv_tab_rar, method = c("everything", "pearson")) 


# Wavelets analysis for the environmental variables-------
library(waveslim)

## input data
### we need environmental variables z-scores----
sample_id_dec_date <- asv_tab_all_bloo_z_tax_02 |>
  dplyr::select(sample_id, decimal_date) |>
  distinct() 

env_z_wavelets_df <- bbmo_env_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  left_join(sample_id_dec_date) |>
  dplyr::select(-sample_id, -env_values) #|>
  #pivot_wider(id_cols = decimal_date, names_from = environmental_variable, values_from = z_score_environmental_variable)

##before imputing values I use the complete variables which are abundance and synechococcus
env_z_wavelets_df <- env_z_wavelets_df |>
  group_by(environmental_variable) |>
  dplyr::filter(n() >= 120)

#### 4 steps 
### 1. modwt computation----

  modwt_results_env <- env_z_wavelets_df  |>
    group_by(environmental_variable) %>%
    summarize(modwt_result = list(modwt.function.biased(z_score_environmental_variable)))

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

 ### I extract the wavelets results at the same time for ALL environmental variables that I have in modwt_results -----

 # Initialize a list to store tibbles
 all_tibbles <- list()
 
 # Loop through the rows of modwt_results_env
 for (i in seq_len(nrow(modwt_results_env))) {
   # Extract the current row
   current_row <- modwt_results_env[i, ]
   current_asv_row.number <- current_row |>
     dplyr::mutate(row_number_asv = row_number()) |>
     dplyr::select(row_number_asv) 
   
   current_asv_row.number <- current_asv_row.number$row_number_asv[1]
   
   # Extract environmental_variable and modwt_result list
   current_environmental_variable <- current_row$environmental_variable
   current_modwt_result <- current_row$modwt_result
   
   # Initialize a list to store tibbles for the current environmental_variable
   d_tibbles <- list()
   
   # Loop through modwt_result
   for (j in 1:5) {
     # Extract the environmental_variable tibble
     d_tibbles[[j]] <- current_modwt_result[[current_asv_row.number]][[j]] %>%
       as_tibble_col(column_name = paste0("d", j))
   }
   
   ##create a column with the environmental_variable
   current_row <- current_environmental_variable %>%
     as_tibble() |>
     mutate(asv_name = list(rep(value, each = 120))) |>
     unnest(asv_name) |>
     dplyr::select(-value)
   
   # Combine the tibbles for the current environmental_variable
   if (length(d_tibbles) > 0) {
     all_tibbles[[current_environmental_variable]] <- bind_cols(d_tibbles) %>%
       dplyr::mutate(sample_num = row_number()) %>%
       bind_cols(current_row)
   }
 }
 
 # Combine all the tibbles into one
 final_tibble <- bind_rows(all_tibbles)
 
  decimal_date_tibble <-  env_z_wavelets_df %$%
     decimal_date |>
     unique() |>
     as_tibble_col(column_name = 'decimal_date') |>
     dplyr::mutate(sample_num = row_number())
  
   # tax <- asv_tab_all_bloo_z_tax |>
   #   dplyr::select(environmental_variable, phylum, class, order, family, genus) |>
   #   distinct()
 
 wavelets_result_env_tibble  <-  final_tibble |>
   pivot_longer(cols = !c(asv_name, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
   left_join( decimal_date_tibble) |>
   #left_join(tax, by = c('asv_name' = 'environmental_variable')) |>
   rename(environmental_variable = asv_name) |>
   dplyr::mutate(wavelets_transformation = str_replace(wavelets_transformation, 'd5', 's4'))
 
 ## I remove the most afected samples by the boundaries it is less biased but still more biased than when applying the brick wall function
 ## as we increase the signal the wavelet gets more affected by the margin effect
 boundaries 
 
 wavelets_result_env_tibble_red <-   wavelets_result_env_tibble |>
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
 
### PLOT THE WAVELETS TRANSFROMATINS COMPUTED (visually inspect them, to be sure of the results)----
 wavelets_result_env_tibble_red |>
   str()
 
 wavelets_result_env_tibble_red |>
   ggplot(aes(as.numeric(decimal_date), as.numeric(wavelets_result_ed)))+
   geom_col()+
   #ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_name), ' ', unique( wavelets_result_tibble_tax$family))) +
   labs(x = 'Decimal date', y = 'Wavelets results')+
   facet_grid(wavelets_transformation~environmental_variable)+
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
   plot_data <- wavelets_result_env_tibble %>%
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
 
 
### 4. Magnitude (coefficients) observe were they got the highest coefficients for the wavelet analysis--------

 #### Wavelet coefficient magnitude indicates how strongly the data are correlated with the mother wavelet at a given frequency and distance.
 
env_type <-  wavelets_result_env_tibble_red %>%
   group_by(environmental_variable, wavelets_transformation) %>%
   dplyr::filter(!is.na(wavelets_result)) |>
   dplyr::summarize(coefficients = sqrt(sum(wavelets_result^2))) |>  # Calculate the magnitude of coefficients for each level
   group_by(environmental_variable) |>
   top_n(1, wt = coefficients) |>  # Find the level with the maximum magnitude
   rename(max_coeff = wavelets_transformation) |>
   dplyr::mutate(bloomer_type = case_when(max_coeff == 's4' ~ 'inter-annual',
                                          max_coeff == 'd1' ~ 'fine-scale',
                                          max_coeff == 'd2' ~ 'half-yearly', 
                                          max_coeff == 'd3' ~ 'seasonal', 
                                          max_coeff == 'd4' ~ 'year-to-year')) 
 
 
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
 
 date_sample <- asv_tab_bloo_rel_abund_z |>
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

year_date <- asv_tab_10y_02_pseudo_rclr_bloo |>
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
  #scale_fill_gradientn( colours = palete_gradient_cb2)+
  #scale_fill_gradient2(low = palete_gradient_cb2[1], high = palete_gradient_cb2[4], midpoint = 0) +
  scale_fill_gradientn(colors = palete_gradient_cb5,
                       breaks=c(-5, 0, 5),
                       limits=c(-5.9,  5.9)) +
  # scale_fill_gradient2(low = "#0049B7",
  #                      high = "#b81131", midpoint = 0)+
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
  #scale_fill_gradientn( colours = palete_gradient_cb)+
  # scale_fill_gradient2(low = "#0049B7",
  #                      high = "#b81131", midpoint = 0)+
  # scale_fill_gradient2(colours = palete_gradient_cb2,
  #                      midpoint = 0  # Set the midpoint value)+
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
 
bloo_events_list |>
   asv_tab_bloo_rel_abund_z |>
   dplyr::filter()

asv_tab_bloo_rel_abund_z |>
  colnames()

asv_tab_bloo_rel_abund_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>

  group_by(decimal_date, BP_FC1.55) |>
  dplyr::summarize(bloom_event = case_when(any(z_score_ra >= 1.96 &
                                                abundance_value >= 0.1 ) ~ 'bloom',
                                           TRUE ~ 'no-bloo')) |>
  ggplot(aes(decimal_date, BP_FC1.55))+
  geom_point(aes(color = bloom_event))+ 
  scale_y_continuous(expand = c(0,0))+
  theme_bw()
 
asv_tab_bloo_rel_abund_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  
  group_by(decimal_date, BP_FC1.55) |>
  dplyr::summarize(bloom_event = case_when(any(z_score_ra >= 1.96 &
                                                 abundance_value >= 0.1 ) ~ 'bloom',
                                           TRUE ~ 'no-bloo')) |>
  ggplot(aes(decimal_date, BP_FC1.55))+
  geom_point(aes(color = bloom_event))+ 
  scale_y_continuous(expand = c(0,0))+
  theme_bw()

bloo_events <- asv_tab_bloo_rel_abund_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(decimal_date, BP_FC1.55) |>
  dplyr::summarize(bloom_event = case_when(any(z_score_ra >= 1.96 &
                                                 abundance_value >= 0.1 ) ~ 'bloom',
                                           TRUE ~ 'no-bloo')) |>
  left_join(bbmo_env)

asv_tab_bloo_rel_abund_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |> 
  dplyr::select(decimal_date, abundance_value) |>
  dplyr::group_by(decimal_date) |>
  dplyr::summarize(max_abundance = max(abundance_value, na.rm = TRUE)) |>
  left_join(bloo_events) |>
  ggplot(aes(max_abundance, BP_FC1.55))+
  geom_point(aes(color = bloom_event))+ 
  scale_y_continuous(expand = c(0,0))+
  geom_smooth(method = 'loess')+
  theme_bw()
  
asv_tab_bloo_rel_abund_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |> 
  dplyr::select(decimal_date, abundance_value) |>
  dplyr::group_by(decimal_date) |>
  dplyr::summarize(max_abundance = max(abundance_value, na.rm = TRUE)) |>
  left_join(bloo_events) |>
  ungroup() |>
  dplyr::filter(!is.na(BP_FC1.55)) |>
  ggplot(aes(bloom_event, as.numeric(BP_FC1.55)))+
  geom_point(aes(size = max_abundance), position = position_jitter(width = 0.05))+ 
  #scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(label = c('Bloom', 'No bloom'))+
  geom_violin(alpha = 0.3, draw_quantiles = T)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  #geom_boxplot(alpha = 0.3)+
  labs(x = '', y = 'Bacterial production (µgC l-1 d-1)', size = 'Maximal\nrelative\nabundance (%)')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 6))

bloo_events |>
  colnames()

### When we have a blooming event do we observe a general increase/decrease in some environmental variables?------
asv_tab_bloo_rel_abund_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |> 
  dplyr::select(decimal_date, abundance_value) |>
  dplyr::group_by(decimal_date) |>
  dplyr::summarize(max_abundance = max(abundance_value, na.rm = TRUE)) |>
  left_join(bloo_events) |>
  ungroup() |>
  dplyr::select(decimal_date, bloom_event,  temperature, day_length, secchi, salinity, chla_total, chla_3um, PO4, NH4, NO2, NO3, Si, PNF_Micro,
                PNF2_5um_Micro, PNF_5um_Micro, cryptomonas, micromonas, HNF_Micro, HNF2_5um_Micro, HNF_5um_Micro, LNA, HNA, prochlorococcus_FC,
                Peuk1, Peuk2, bacteria_joint, synechococcus, max_abundance) |>
  pivot_longer(cols = -c(decimal_date, bloom_event, max_abundance)) |>
  #dplyr::filter(!is.na(BP_FC1.55)) |>
  ggplot(aes(bloom_event, as.numeric(value)))+
  geom_point(aes(size = max_abundance), position = position_jitter(width = 0.05))+ 
  #scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(label = c('Bloom', 'No bloom'))+
  #geom_violin(alpha = 0.3, draw_quantiles = T)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  facet_wrap(vars(name), scale = 'free')+
  geom_boxplot(alpha = 0.3)+
  labs(x = '', y = '', size = 'Maximal\nrelative\nabundance (%)')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 6))

##same plot but with z-scores
date_sample <- asv_tab_bloo_rel_abund_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(decimal_date, sample_id) |>
  group_by(decimal_date, sample_id) |>
  distinct()

bbmo_env_z_w <- bbmo_env_z |>
  dplyr::select(-env_values) |>
  pivot_wider(id_cols = sample_id, values_from = z_score_environmental_variable, names_from = environmental_variable ) |>
  left_join(date_sample) |>
  dplyr::select(-sample_id)

max_abund <- asv_tab_bloo_rel_abund_z |>
  dplyr::select(decimal_date, abundance_value) |>
  dplyr::group_by(decimal_date) |>
  dplyr::summarize(max_abundance = max(abundance_value, na.rm = TRUE))

bloom_events_env_variables <-  asv_tab_bloo_rel_abund_z |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(decimal_date, BP_FC1.55) |>
  dplyr::summarize(bloom_event = case_when(any(z_score_ra >= 1.96 &
                                                 abundance_value >= 0.1 ) ~ 'bloom',
                                           TRUE ~ 'no-bloo')) |>
  dplyr::select(-BP_FC1.55) |>
  left_join(bbmo_env_z_w) |>
  left_join(max_abund) |>
  ungroup() |>
  # dplyr::select(decimal_date, bloom_event,  temperature, day_length, secchi, salinity, chla_total, chla_3um, PO4, NH4, NO2, NO3, Si, PNF_Micro,
  #               PNF2_5um_Micro, PNF_5um_Micro, cryptomonas, micromonas, HNF_Micro, HNF2_5um_Micro, HNF_5um_Micro, LNA, HNA, prochlorococcus_FC,
  #               Peuk1, Peuk2, bacteria_joint, synechococcus, max_abundance) |>
  pivot_longer(cols = -c(decimal_date, bloom_event, max_abundance)) |>
  #dplyr::filter(!is.na(BP_FC1.55)) |>
  ggplot(aes(bloom_event, as.numeric(value)))+
  geom_point(aes(size = max_abundance), position = position_jitter(width = 0.05))+ 
  #scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(label = c('Bloom', 'No bloom'))+
  scale_size_continuous(range = c(0, 2)) + 
  geom_violin(alpha = 0.3, draw_quantiles = T)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  facet_wrap(vars(name),  labeller = labs_env, ncol = 3)+
  #geom_boxplot(alpha = 0.3)+
  labs(x = '', y = 'z-score', size = 'Maximum\nrelative\nabundance (%)')+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 6),
        strip.background = element_blank())

ggsave('bloom_events_env_variables.pdf', bloom_events_env_variables,
      path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/",
      width = 188,
      height = 280,
      units = 'mm')
 
## Do the same but separating no-bloom / bloom / super blooming events------


### other idea from Marta
##### y has probado haciendo distancia euclidea de vbles ambientales respecto al tiempo anterior y ver si hay relación con los blooms? aunque quizá eso ya te lo dice el análisis de Carmen?

## Wavelets analysis with environmental data (NO MISSING VALUES) and cross correlate ------
library(waveslim)

## dataset with al the environmental values interpolated
env_data_interpolated_values_all_z_score 

## input data
### we need environmental variables z-scores----
sample_id_dec_date <- asv_tab_all_bloo_z_tax_02 |>
  dplyr::select(sample_id, decimal_date) |>
  distinct() 

env_z_wavelets_df <- env_data_interpolated_values_all_z_score |>
  dplyr::select(-env_values) #|>
#pivot_wider(id_cols = decimal_date, names_from = environmental_variable, values_from = z_score_environmental_variable)

#### 4 steps 
### 1. modwt computation----

modwt_results_env <- env_z_wavelets_df  |>
  group_by(environmental_variable) %>%
  summarize(modwt_result = list(modwt.function.biased(z_score_environmental_variable)))

### 2. e-folding-----
###### commmon for all 
x <- rep(0, 10001) 
x[5000] <- 1 
n.levels <- 4 
len <- 120 ##120 (length of my dataset)
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

### 3. Visualize the results obtained from the modwt transfomation ------

### I extract the wavelets results at the same time for ALL environmental variables that I have in modwt_results -----

# Initialize a list to store tibbles
all_tibbles <- list()

# Loop through the rows of modwt_results_env
for (i in seq_len(nrow(modwt_results_env))) {
  # Extract the current row
  current_row <- modwt_results_env[i, ]
  current_asv_row.number <- current_row |>
    dplyr::mutate(row_number_asv = row_number()) |>
    dplyr::select(row_number_asv) 
  
  current_asv_row.number <- current_asv_row.number$row_number_asv[1]
  
  # Extract environmental_variable and modwt_result list
  current_environmental_variable <- current_row$environmental_variable
  current_modwt_result <- current_row$modwt_result
  
  # Initialize a list to store tibbles for the current environmental_variable
  d_tibbles <- list()
  
  # Loop through modwt_result
  for (j in 1:5) {
    # Extract the environmental_variable tibble
    d_tibbles[[j]] <- current_modwt_result[[current_asv_row.number]][[j]] %>%
      as_tibble_col(column_name = paste0("d", j))
  }
  
  ##create a column with the environmental_variable
  current_row <- current_environmental_variable %>%
    as_tibble() |>
    mutate(asv_name = list(rep(value, each = 120))) |>
    unnest(asv_name) |>
    dplyr::select(-value)
  
  # Combine the tibbles for the current environmental_variable
  if (length(d_tibbles) > 0) {
    all_tibbles[[current_environmental_variable]] <- bind_cols(d_tibbles) %>%
      dplyr::mutate(sample_num = row_number()) %>%
      bind_cols(current_row)
  }
}

# Combine all the tibbles into one
final_tibble <- bind_rows(all_tibbles)

decimal_date_tibble <-  env_z_wavelets_df %$%
  decimal_date |>
  unique() |>
  as_tibble_col(column_name = 'decimal_date') |>
  dplyr::mutate(sample_num = row_number())

# tax <- asv_tab_all_bloo_z_tax |>
#   dplyr::select(environmental_variable, phylum, class, order, family, genus) |>
#   distinct()

wavelets_result_env_tibble  <-  final_tibble |>
  pivot_longer(cols = !c(asv_name, sample_num), values_to = 'wavelets_result', names_to = 'wavelets_transformation') |>
  left_join( decimal_date_tibble) |>
  #left_join(tax, by = c('asv_name' = 'environmental_variable')) |>
  rename(environmental_variable = asv_name) |>
  dplyr::mutate(wavelets_transformation = str_replace(wavelets_transformation, 'd5', 's4'))

## I remove the most afected samples by the boundaries it is less biased but still more biased than when applying the brick wall function
## as we increase the signal the wavelet gets more affected by the margin effect
boundaries 

wavelets_result_env_tibble_red <-   wavelets_result_env_tibble |>
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

### PLOT THE WAVELETS TRANSFROMATINS COMPUTED (visually inspect them, to be sure of the results)----
wavelets_result_env_tibble_red |>
  str()

wavelets_result_env_tibble_red |>
  ggplot(aes(as.numeric(decimal_date), as.numeric(wavelets_result_ed)))+
  geom_col()+
  #ggtitle(paste0(unique(wavelets_result_tibble_tax$asv_name), ' ', unique( wavelets_result_tibble_tax$family))) +
  labs(x = 'Decimal date', y = 'Wavelets results')+
  facet_grid(wavelets_transformation~environmental_variable)+
  scale_x_continuous(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        #aspect.ratio = 4/10,
        text = element_text(size = 5))

### create a loop to save all the wavelets transformations computed and visually inspect them
# Create a folder to save the PDF files
dir.create("../results/figures/wavelets_plots_env/", showWarnings = T)

### ENVIRONMENTAL VARIABLES
# Get unique environmental_variable values
environmental_variables <- unique(wavelets_result_env_tibble_red$environmental_variable)

for (environmental_variable in environmental_variables) {
  # Filter data for the current environmental_variable
  plot_data <- wavelets_result_env_tibble %>%
    dplyr::filter(environmental_variable == !!environmental_variable)  # Use !! to unquote environmental_variable
  
  title_plot <- plot_data |>
    group_by(environmental_variable) |>
    dplyr::reframe(title_ed = paste0(environmental_variable)) |>
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
  pdf_file <- paste0("../results/figures/wavelets_plots_env/", environmental_variable, "_plot_env.pdf")
  ggsave(pdf_file, p, width = 6, height = 3, units = "in")
  
  # Print a message indicating the plot has been saved
  cat("Plot for", environmental_variable, "saved as", pdf_file, "\n")
}

dev.off()

### 4. Magnitude (coefficients) observe were they got the highest coefficients for the wavelet analysis--------

#### Wavelet coefficient magnitude indicates how strongly the data are correlated with the mother wavelet at a given frequency and distance.

env_type <-  wavelets_result_env_tibble_red %>%
  group_by(environmental_variable, wavelets_transformation) %>%
  dplyr::filter(!is.na(wavelets_result)) |>
  dplyr::summarize(coefficients = sqrt(sum(wavelets_result^2))) |>  # Calculate the magnitude of coefficients for each level
  group_by(environmental_variable) |>
  top_n(2, wt = coefficients) |>  # Find the level with the maximum magnitude
  rename(max_coeff = wavelets_transformation) |>
  dplyr::mutate(bloomer_type = case_when(max_coeff == 's4' ~ 'inter-annual',
                                         max_coeff == 'd1' ~ 'fine-scale',
                                         max_coeff == 'd2' ~ 'half-yearly', 
                                         max_coeff == 'd3' ~ 'seasonal', 
                                         max_coeff == 'd4' ~ 'year-to-year')) 

##i create a table with all the coefficients (I do not decide which is the most important)-----
wavelets_result_env_tibble_red_coeff <- wavelets_result_env_tibble_red %>%
  dplyr::mutate(wavelets_result_ed = as.numeric(wavelets_result_ed)) |>
  dplyr::filter(!is.na(wavelets_result_ed)) |>
  group_by(environmental_variable, wavelets_transformation) %>%
  dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) 

wavelets_result_env_tibble_red_coeff <- wavelets_result_env_tibble_red_coeff |>
  dplyr::mutate(environmental_variable = str_replace(environmental_variable,'_no_nas',''))

wavelets_result_env_tibble_red_coeff$environmental_variable <- factor(wavelets_result_env_tibble_red_coeff$environmental_variable, levels = c("day_length", "temperature" ,"secchi" , "salinity" ,      
                                                                                          "PO4",  "NH4" ,  "NO2" , "NO3" , "Si" ,  "chla_total" ,  "chla_3um", "synechococcus", "prochlorococcus_FC", 
                                                                                          "bacteria_joint", "LNA", "HNA", "Peuk1",             
                                                                                          "Peuk2",
                                                                                          "BP_FC1.55", "PNF_Micro" , "PNF2_5um_Micro",  "PNF_5um_Micro", "cryptomonas", "micromonas","HNF_Micro",         
                                                                                          "HNF2_5um_Micro", "HNF_5um_Micro",
                                                                                          "low_vlp" ,
                                                                                          "med_vlp" ,  
                                                                                          "high_vlp" ,
                                                                                          "total_vlp"))

labs_wavelets <- as_labeller(c('d1' = 'Fine-scale',
                               'd2' = 'Half-yearly',
                               'd3' = 'Seasonal',
                               'd4' = 'Year-to-year',
                               's4' = 'Inter-annual'))

wavelets_result_env_tibble_red_coeff_plot <- wavelets_result_env_tibble_red_coeff |>
  ggplot(aes( wavelets_transformation, environmental_variable, fill = coefficients))+
  geom_tile()+
  scale_y_discrete(labels = labs_env)+
  scale_fill_gradientn(colors = palette_gradient_bw)+
  scale_x_discrete(labels = labs_wavelets)+
  theme_bw()+
  labs(x = 'Wavelet vectors', y = 'Environmental variables', fill = 'Coefficient')+
  theme(panel.border = element_blank(), text = element_text(size = 4))
  
wavelets_result_env_tibble_red_coeff_plot

# ggsave('wavelets_result_env_tibble_red_coeff.pdf', wavelets_result_env_tibble_red_coeff_plot,
#      path = "~/Documentos/Doctorat/BBMO/BBMO_bloomers/Results/Figures/wavelets_plots/",
#      width = 88,
#      height = 100,
#      units = 'mm')

### environmental distance between samples ----
## Upload environmental to model with interpolated missing variables----
env_data_interpolated_values_all <- read.csv2('data/env_data/env_data_interpolated_values_all.csv') |>
  rename(sample_id_num = X)

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
  #bind_rows(euclidean_distance_tb_g) |>
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
  #scale_linetype_manual( labels = labs_fraction_env, values = c('0.2' = 1, '3' = 2, 'env' = 1))+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Date', y = 'Euclidean Distance', color = '', linetype = '')+
  #scale_shape_discrete(labels = labs_fraction)+
  guides(shape = 'none',
         color = guide_legend(ncol = 3, keywidth = unit(1, "cm")))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 7),
    strip.text = element_text(size = 5),
    axis.text = element_text(size = 5),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(),
    strip.placement = 'outside', 
    plot.margin = margin(t = 5, l = 20, r = 20))

eucl_plot
 
m_02 |>
  colnames()

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
  #geom_hline(yintercept = 0.1, linetype = 'dashed')+
  scale_y_continuous(labels = scientific_format(), expand = c(0,0))+
  #scale_color_identity()+
  #geom_smooth( aes(color = name, fill = name), span = 0.1)+
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
  scale_color_manual(values = palette_seasons_4)+
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
  scale_color_manual(values = palette_seasons_4)+
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
asv_tab_10y_02_pseudo_rclr |>
  colnames()

tb_w <- asv_tab_10y_02_pseudo_rclr  |>
  dplyr::ungroup() |>
  dplyr::select(rclr, asv_num, sample_id) |>
  dplyr::filter( str_detect(sample_id, '_0.2_')) |>
  pivot_wider(values_from = rclr, names_from = asv_num)

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
tb_w <- asv_tab_10y_3_pseudo_rclr  |>
  dplyr::ungroup() |>
  dplyr::select(rclr, asv_num, sample_id) |>
  dplyr::filter( str_detect(sample_id, '_3_')) |>
  pivot_wider(values_from = rclr, names_from = asv_num)

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
  ggplot(aes(euclidean_distance_env_phch, euclidean_distance))+
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
  
