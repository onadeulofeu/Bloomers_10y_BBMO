# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                         mTAGs & mOTUs analysis              ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                       reads from metagenome                 ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# upload packages ---- 
library(stringr)
library(magrittr)
library(tidyverse)
library(Bloomers)
library(scales)
library(readxl)
library(EcolUtils)
library(gridExtra)
library(cowplot)

# upload functions ----
source('../../Bloomers/R/find_asv_with_anomalies.R')
source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R')

# mTAGs come from metagenomes, are reads that map with the 16s rRNA SILVA database.

## palettes ---- 
palette_motus_asvs_mtags <-  c('mTAGs' = '#344848',
                               'ASVs' = '#85DFC9',
                               'mOTUs' = '#EB9073',
                               'Flow Citomentry' = '#6B5C57')

palette_long <- c( "#fcca46",  "#009e73",
                   "#0051BF",   '#005c69',  "#69267e",
                   "#8c789d",'#FFA737',  '#B0413E', '#009F6A',  '#BB4430',   
                   '#F2AC5D', '#FFA200','#F35900',  '#8C000A', 
                   '#FFA197', '#74B9C8','#e3a6ce',  '#534F4A', 
                   '#46ACC2', '#2B4162', '#ffffff','#bbbbbb','#000000','#004501','#e3a6ce', '#D3BF27',    
                   '#C55E5C',  '#B0413E','#960200',  '#721817', '#CD7F78',  '#B84A62',                     
                   '#8C000A',  '#70161E', '#A63B00', '#F2AC5D',  '#FFA200',  '#A05C00', '#F35900', 
                   '#DE6931', '#FF8E00', '#fcca46',"#fcca46",  "#009e73","#0051BF",   '#005c69',  "#69267e",
                   "#8c789d",'#FFA737',  '#B0413E', '#009F6A',  '#BB4430', '#F2AC5D', '#FFA200','#F35900',  '#8C000A', 
                   '#FFA197', '#74B9C8','#e3a6ce',  '#534F4A', 
                   '#46ACC2', '#2B4162', '#ffffff','#bbbbbb','#000000',
                   '#004501','#e3a6ce', '#D3BF27', '#C55E5C',  '#B0413E','#960200',  '#721817','#CD7F78',  '#B84A62',                     
                   '#8C000A',  '#70161E', '#A63B00','#F2AC5D',  '#FFA200',  '#A05C00', '#F35900', 
                   '#DE6931', '#FF8E00', '#fcca46',      "#0051BF",   '#005c69',  "#69267e",
                   "#8c789d",'#FFA737',  '#B0413E', '#009F6A',  '#BB4430', '#F2AC5D', '#FFA200',
                   '#F35900',  '#8C000A', '#FFA197', '#74B9C8','#e3a6ce',  '#534F4A', 
                   '#46ACC2', '#2B4162', '#ffffff','#bbbbbb','#000000',
                   '#004501','#e3a6ce', '#D3BF27',    '#C55E5C',  '#B0413E','#960200',  '#721817',
                   '#CD7F78',  '#B84A62', '#8C000A',  '#70161E', '#A63B00', 
                   '#F2AC5D',  '#FFA200',  '#A05C00', '#F35900', 
                   '#DE6931', '#FF8E00', '#fcca46',"#0051BF",   '#005c69',  "#69267e",
                   "#8c789d",'#FFA737',  '#B0413E', '#009F6A',  '#BB4430', '#F2AC5D', '#FFA200','#F35900',  '#8C000A', 
                   '#FFA197', '#74B9C8','#e3a6ce',  '#534F4A', 
                   '#46ACC2', '#2B4162', '#ffffff','#bbbbbb','#000000',
                   '#004501','#e3a6ce', '#D3BF27','#C55E5C',  '#B0413E',
                   '#960200',  '#721817','#CD7F78',  '#B84A62','#8C000A',  '#70161E', '#A63B00', 
                   '#F2AC5D',  '#FFA200',  '#A05C00', '#F35900', 
                   '#DE6931', '#FF8E00', '#fcca46') 
               
palette_class_assigned_mtags_motus <- c('Gammaproteobacteria' = '#fcca46', 'Alphaproteobacteria' = '#8C000A', 
                                'Zetaproteobacteria' = '#EBCD92',
                                'Bacteroidia' = "#0051BF",'Rhodothermia' =  '#00425E',
                                'Ignavibacteria' = '#CAFFFF', 'Chlorobia' = '#5BD1FF', 'Kryptonia' = '#0071A8',
                                'Nitrososphaeria' = '#FFA737',
                                'Cyanobacteriia' = '#009E73', 'Vampirivibrionia' = '#00733C',
                                'Acidimicrobiia' = '#B0413E','Actinobacteria' = '#C55E5C',
                                'Coriobacteriia' = '#B44240', 'Thermoleophilia' = '#AB3F3D',
                                'Verrucomicrobiae' = '#CCCCCC', #'Kiritimatiellae' = '#6D6D6D',
                                'Lentisphaeria' = '#424242', 'Omnitrophia' = '#6D6D6D','Chlamydiae' = '#9B9B9B',
                                'Planctomycetes' = '#69267E', 'Phycisphaerae' = '#BE8DCB','Blastocatellia' = '#00ADFF',
                                'Holophagae' = '#86A4C6', 'Vicinamibacteria' = '#002671',
                                'Acidobacteriae' = '#1F78B4', 'Thermoanaerobaculia' = '#0051BF',
                                'Bdellovibrionia' = '#8C789D','Oligoflexia' = '#654584',
                                'Bacilli' = '#637B88',  'Clostridia' = '#384F5B',
                                'Negativicutes' = '#91AAB8', 'Syntrophomonadia' = '#C2DCEA',
                                'Desulfitobacteriia' =  '#0E2732', 'Thermoanaerobacteria' = '#005474',
                                'Myxococcia' = '#2E5A51','Polyangia' = '#003029',
                                'Kiritimatiellae' = '#e3a6ce', 'Campylobacteria' = '#002960',
                                'Deinococci' = '#ba9864','Fusobacteriia' = '#fb9a99',
                                'Desulfovibrionia' = '#005c69', 'Desulfobulbia' = '#3D518E',
                                'Desulfuromonadia' = '#6D7DBE', 'Desulfobacteria' = '#000036', 'Methylomirabilia' = '#ffff00',
                                'Gemmatimonadetes' = '#c55e5c',
                                'Anaerolineae' = '#00d198', 'Chloroflexia' = '#009F6A',
                                'Spirochaetia' = '#5CB1D6', 'Leptospirae' = '#005576', 
                                'Calditrichia'  = '#8a007a', 'Halobacteria' = '#b79c64', 
                                'Nitrospiria' = '#41815f',
                                'Babeliae' = '#5b95e5', 'Saccharimonadia' = '#33af9c', 
                                'Cloacimonadia' = '#fbed5c',
                                'Synergistia' = '#ce7800', 'Abditibacteria' = '#87878b', 'Dehalococcoidia' =  '#005c69',                                      
                                'Thermoplasmata' =     '#5b95e5',
                                'Deferribacteres' = '#4dbaa9', 'BD2-11terrestrialgroup' = '#ba9864', 'Dadabacteriia' = '#1F78B4',
                                'NA' = "#ffffff",'others' = '#000000', 'unknown' = "#ffffff",
                                'Unassigned' = '#FBFBFB',
                                'Unaligned' = '#000036',
                                '<NA>' = '#fbed5c', 'Gammaproteobacteria' = '#fcca46', 'Alphaproteobacteria' = '#8C000A', 
                                  'Zetaproteobacteria' = '#EBCD92',
                                  'Bacteroidia' = "#0051BF",'Rhodothermia' =  '#00425E',
                                  'Ignavibacteria' = '#CAFFFF', 'Chlorobia' = '#5BD1FF', 'Kryptonia' = '#0071A8',
                                  'Nitrososphaeria' = '#FFA737',
                                  'Cyanobacteriia' = '#009E73', 'Vampirivibrionia' = '#00733C',
                                  'Acidimicrobiia' = '#B0413E','Actinobacteria' = '#C55E5C',
                                  'Coriobacteriia' = '#B44240', 'Thermoleophilia' = '#AB3F3D',
                                  'Verrucomicrobiae' = '#CCCCCC', #'Kiritimatiellae' = '#6D6D6D',
                                  'Lentisphaeria' = '#424242', 'Omnitrophia' = '#6D6D6D','Chlamydiae' = '#9B9B9B',
                                  'Planctomycetes' = '#69267E', 'Phycisphaerae' = '#BE8DCB','Blastocatellia' = '#00ADFF',
                                  'Holophagae' = '#86A4C6', 'Vicinamibacteria' = '#002671',
                                  'Acidobacteriae' = '#1F78B4', 'Thermoanaerobaculia' = '#0051BF',
                                  'Bdellovibrionia' = '#8C789D','Oligoflexia' = '#654584',
                                  'Bacilli' = '#637B88',  'Clostridia' = '#384F5B',
                                  'Negativicutes' = '#91AAB8', 'Syntrophomonadia' = '#C2DCEA',
                                  'Desulfitobacteriia' =  '#0E2732', 'Thermoanaerobacteria' = '#005474',
                                  'Myxococcia' = '#2E5A51','Polyangia' = '#003029',
                                  'Kiritimatiellae' = '#e3a6ce', 'Campylobacteria' = '#002960',
                                  'Deinococci' = '#ba9864','Fusobacteriia' = '#fb9a99',
                                  'Desulfovibrionia' = '#005c69', 'Desulfobulbia' = '#3D518E',
                                  'Desulfuromonadia' = '#6D7DBE', 'Desulfobacteria' = '#000036', 'Methylomirabilia' = '#ffff00',
                                  'Gemmatimonadetes' = '#c55e5c',
                                  'Anaerolineae' = '#00d198', 'Chloroflexia' = '#009F6A',
                                  'Spirochaetia' = '#5CB1D6', 'Leptospirae' = '#005576', 
                                  'Calditrichia'  = '#8a007a', 'Halobacteria' = '#b79c64', 
                                  'Nitrospiria' = '#41815f',
                                  'Babeliae' = '#5b95e5', 'Saccharimonadia' = '#33af9c', 
                                  'Cloacimonadia' = '#fbed5c',
                                  'Synergistia' = '#ce7800', 'Abditibacteria' = '#87878b', 'Dehalococcoidia' =  '#005c69',                                      
                                  'Thermoplasmata' =     '#5b95e5',
                                  'Deferribacteres' = '#4dbaa9', 'BD2-11terrestrialgroup' = '#ba9864', 'Dadabacteriia' = '#1F78B4',
                                  'NA' = "#ffffff",'others' = '#000000', 'unknown' = "#ffffff",
                                  'unassigned' = '#FBFBFB',
                                  'Marinisomatia'= '#7B968F',
                                  'Planctomycetia' = '#EB9073',
                                  'UBA1144' = '#85DFC9',
                                  'Poseidoniia' = '#6B5C57',
                                  #'Unaligned' = '#000036',
                                  '<NA>' = '#fbed5c') 

# The goal of this code is to observe if the changes in abundances of my bloomers are also observed in the metagenomic data either in mOTUs or in mTAGS----

# -------------------- ################## mTAGs ##################  -------
# upload mTAGs data -----
mtags_table <- read_tsv('data/mtags_2009_2013_ramiro/bbmo.2008.2022.mtags/bbmo.2008.2022.mtags.otu.tsv') |>
  as_tibble() |>
  rename(tax = '#taxpath')

mtags_table_unassigned <-  mtags_table |>
  dplyr::filter(str_detect(tax, 'Unassigned') |
                  str_detect(tax, 'Unaligned')) |>
  dplyr::mutate(domain = case_when(tax == 'Unassigned' ~ 'Unassigned',
                                   tax == 'Unaligned' ~ 'Unaligned'),
                phylum = case_when(tax == 'Unassigned' ~ 'Unassigned',
                                   tax == 'Unaligned' ~ 'Unaligned'),
                class = case_when(tax == 'Unassigned' ~ 'Unassigned',
                                  tax == 'Unaligned' ~ 'Unaligned'),
                order = case_when(tax == 'Unassigned' ~ 'Unassigned',
                                  tax == 'Unaligned' ~ 'Unaligned'), 
                family = case_when(tax == 'Unassigned' ~ 'Unassigned',
                                   tax == 'Unaligned' ~ 'Unaligned'), 
                genus = case_when(tax == 'Unassigned' ~ 'Unassigned',
                                  tax == 'Unaligned' ~ 'Unaligned'), 
                motu_num = case_when(tax == 'Unassigned' ~ 'Unassigned',
                                     tax == 'Unaligned' ~ 'Unaligned')) |>
  dplyr::select(-tax)

mtags_table_unassigned_ed <- mtags_table_unassigned |>
  dplyr::select(-domain, -phylum, -class, -order, -family, -genus)

mtags_tax_unassigned <- mtags_table_unassigned |>
  dplyr::select(domain, phylum, class, order, family, genus, motu_num)

mtags_table <- mtags_table[-c(37240, 37241),]

mtags_table_ed <- mtags_table |> # i add at the end the unassigned rows
  dplyr::select(-tax) |>
  dplyr::mutate(motu_num = paste0('mtags', 1:nrow(mtags_table))) |>
  bind_rows(mtags_table_unassigned_ed)

# mtags_table |>
#   dplyr::select(-tax) |>
#   colSums()

 mtags_tax <- mtags_table |>
  separate(tax, sep = '__', into = c('root1', 'domain',  'phylum', 'class',
                                    'order',  'family',  'genus','otu', 'database'), remove = F) |>
   dplyr::mutate(domain = str_replace(domain, 'Root;','')) |>
   dplyr::mutate(phylum = str_replace(phylum, ';phylum','')) |>
   dplyr::mutate(class = str_replace(class, ';class','')) |>
   dplyr::mutate(order = str_replace(order, ';order','')) |>
   dplyr::mutate(family = str_replace(family, ';family','')) |>
   dplyr::mutate(genus = str_replace(genus, ';genus','')) |>
   dplyr::mutate(otu = str_replace(otu, ';otu','')) |>
   dplyr::select(domain, phylum, class, order, family, genus, otu, tax) |>
   dplyr::mutate(motu_num = paste0('mtags', 1:nrow(mtags_table))) |>
   dplyr::select(-domain, domain = phylum, phylum = class, class = order, order = family, family = genus, genus = otu, motu_num, tax)
 
 mtags_tax <- mtags_tax |>
    bind_rows(mtags_tax_unassigned)
    
mtags_tax |>
   dplyr::filter(domain == 'Eukaryota')  #5829 mtags

 m_mtags <- mtags_table |>
   colnames() |>
   as_tibble_col(column_name = 'sample_id') |>
   dplyr::filter(sample_id != 'tax') |>
   dplyr::mutate(sample_id = str_replace(sample_id, '.bins', '')) |>
   tidyr::separate_wider_position(col = sample_id, widths = c('station' = 2, 'year' = 2, 'month' = 2, 'day' = 2), too_many = 'debug') |>
   dplyr::mutate(date = paste0(day, '-', month, '-', year)) |>
   dplyr::mutate(season = case_when(month %in% c('04', '05') ~ 'spring', ## i would like to add a column with the season ----
                                    month %in% c('01', '02') ~ 'winter',
                                    month %in% c('07', '08') ~ 'summer',
                                    month %in% c('10', '11') ~ 'autumn',
                                    month == '03' & day %in% c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                                                               '13', '14', '15', '16', '17', '18', '19', '20') ~ 'winter',
                                    month == '03' & !day %in% c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                                                                '13', '14', '15', '16', '17', '18', '19', '20') ~ 'spring',
                                    month == '06' & day %in% c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                                                               '13', '14', '15', '16', '17', '18', '19', '20') ~ 'spring',
                                    month == '06' & !day %in% c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                                                                '13', '14', '15', '16', '17', '18', '19', '20') ~ 'summer',
                                    month == '09' & day %in% c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                                                               '13', '14', '15', '16', '17', '18', '19', '20') ~ 'summer',
                                    month == '09' & !day %in% c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                                                                '13', '14', '15', '16', '17', '18', '19', '20') ~ 'autumn',
                                    month == '12' & day %in% c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                                                               '13', '14', '15', '16', '17', '18', '19', '20') ~ 'autumn',
                                    month == '12' & !day %in% c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                                                                '13', '14', '15', '16', '17', '18', '19', '20') ~ 'winter'))
  
  mtags_tax_prok <- mtags_tax |>
   dplyr::filter(domain != 'Eukaryota') |>
   dplyr::filter(order != 'Chloroplast') |>
   dplyr::filter(family != 'Mitochondria')

 mtags_table_ed_pk <- mtags_table_ed |>
   dplyr::filter(motu_num %in%  mtags_tax_prok$motu_num)
 
 mtags_table_ed_pk_rel <- mtags_table_ed_pk |>
   pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'reads') |>
   calculate_rel_abund(group_cols = sample_id)
   
 mtags_table_ed_pk_rel_tax_m <- mtags_table_ed_pk_rel |>
   left_join(mtags_tax) |>
   dplyr::mutate(sample_id = str_replace(sample_id, '.bins', '')) |>
   left_join(m_mtags, by = 'sample_id')
 
# mtags_plot <-  mtags_table_ed_pk_rel_tax_m |>
#    #dplyr::filter(!is.na(phylum))   |>
#    group_by(class, sample_id, date) |>
#    dplyr::reframe(sum_rel = sum(relative_abundance)) |>
#   #left_join(m_mtags, by = c('sample_id') |>
#    dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#    ggplot(aes(date, sum_rel))+
#    geom_area(aes(fill = class, group = class), alpha = 0.9,  position='stack')+
#    scale_fill_manual(values = palette_long)+
#    scale_x_datetime(expand = c(0,0), 
#                     breaks = (by = '1 year'),
#                     date_labels = "%Y")+
#    scale_y_continuous(labels = percent_format(), expand = c(0,0))+
#    labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class')+
#    #facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
#    guides(fill = guide_legend(ncol = 6, size = 10,
#                               override.aes = aes(label = '')),
#           alpha = 'none')+
#    theme_bw()+
#    theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
#          panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
#          legend.position = 'bottom', axis.text.y = element_text(size = 8),
#          axis.title = element_text(size = 8), strip.background = element_blank(), 
#          legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  
#  
# mtags_plot 

# mtags_plot_red <-  mtags_table_ed_pk_rel_tax_m |>
#   #dplyr::filter(!is.na(phylum))   |>
#   group_by(class, sample_id, date) |>
#   dplyr::reframe(sum_rel = sum(relative_abundance)) |>
#   #left_join(m_mtags, by = c('sample_id') |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |>
#   ggplot(aes(date, sum_rel))+
#   geom_area(aes(fill = class, group = class), alpha = 0.9,  position='stack')+
#   #scale_fill_manual(values = palette_long)+
#   scale_x_datetime(expand = c(0,0), 
#                    breaks = (by = '1 year'),
#                    date_labels = "%Y")+
#   scale_y_continuous(labels = percent_format(), expand = c(0,0))+
#   labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class')+
#   #facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
#   guides(fill = guide_legend(ncol = 6, size = 10,
#                              override.aes = aes(label = '')),
#          alpha = 'none')+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
#         legend.position = 'bottom', axis.text.y = element_text(size = 8),
#         axis.title = element_text(size = 8), strip.background = element_blank(), 
#         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside')  
# 
# mtags_plot_red 

## The bloomers (metabarcoding) that I have during those years are ASV11 Alteromonadaceae, ASV15 Clade I, ASV27 Rhodobacteraceae, ASV58 Flavobacteriaceae
## ASV 1 Cyanobacteriaceae, ASV7 Cyanobacteriaceae

## look for my bloomers during those years that match both timeseries -----
# mtags_table_ed_pk_rel_tax_m |>
#   #dplyr::filter(!is.na(phylum))   |>
#   #group_by(class, sample_id, date) |>
#   #dplyr::reframe(sum_rel = sum(relative_abundance)) |>
#   #left_join(m_mtags, by = c('sample_id') |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |> 
#   dplyr::filter(class == 'Gammaproteobacteria') |>
#   dplyr::ungroup() |>
#   arrange(-relative_abundance) |>
#   dplyr::distinct(motu_num)

## LOOKS LIKE MY ASV11 ----
# mtags_table_ed_pk_rel_tax_m |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |> 
#   dplyr::filter(motu_num == 'mtags14723') |>
#   ggplot(aes(date, relative_abundance))+
#   geom_line(aes(group = motu_num))

# ## 
# mtags_max_abund <- mtags_table_ed_pk_rel_tax_m |>
#   #dplyr::filter(!is.na(phylum))   |>
#   #group_by(class, sample_id, date) |>
#   #dplyr::reframe(sum_rel = sum(relative_abundance)) |>
#   #left_join(m_mtags, by = c('sample_id') |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |> 
#   dplyr::filter(family == 'Cyanobacteriales') |>
#   dplyr::ungroup() |>
#   arrange(-relative_abundance) |>
#   dplyr::distinct(motu_num) |>
#   slice_head(n = 10)
# 
# mtags_table_ed_pk_rel_tax_m |>
#   #dplyr::filter(!is.na(phylum))   |>
#   #group_by(class, sample_id, date) |>
#   #dplyr::reframe(sum_rel = sum(relative_abundance)) |>
#   #left_join(m_mtags, by = c('sample_id') |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |> 
#   dplyr::filter(motu_num %in% mtags_max_abund$motu_num) |>
#   ggplot(aes(date, relative_abundance))+
#   geom_line(aes(group = motu_num))+
#   facet_wrap(vars(motu_num))
# 
# mtags_table_ed_pk_rel_tax_m$family |>
#   unique()
# 
# mtags_max_abund <- mtags_table_ed_pk_rel_tax_m |>
#   #dplyr::filter(!is.na(phylum))   |>
#   #group_by(class, sample_id, date) |>
#   #dplyr::reframe(sum_rel = sum(relative_abundance)) |>
#   #left_join(m_mtags, by = c('sample_id') |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |> 
#   dplyr::filter(family == 'Rhodobacterales') |>
#   dplyr::ungroup() |>
#   arrange(-relative_abundance) |>
#   dplyr::distinct(motu_num) |>
#   slice_head(n = 10)
# 
# mtags_table_ed_pk_rel_tax_m |>
#   #dplyr::filter(!is.na(phylum))   |>
#   #group_by(class, sample_id, date) |>
#   #dplyr::reframe(sum_rel = sum(relative_abundance)) |>
#   #left_join(m_mtags, by = c('sample_id') |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
#   dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |> 
#   dplyr::filter(motu_num %in% mtags_max_abund$motu_num) |>
#   ggplot(aes(date, relative_abundance))+
#   geom_line(aes(group = motu_num))+
#   facet_wrap(vars(motu_num))

## Look for potential BLOOMERS in mtags dataset ------
mtags_table_ed_pk_rel_tax_m_red <- mtags_table_ed_pk_rel_tax_m |>
  #dplyr::filter(!is.na(phylum))   |>
  #group_by(class, sample_id, date) |>
  #dplyr::reframe(sum_rel = sum(relative_abundance)) |>
  #left_join(m_mtags, by = c('sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) 

## colored only those taxa that represent > 1% at some point -----
mtags_table_ed_pk_community_class <- mtags_table_ed_pk_rel_tax_m_red |>
  #dplyr::filter(!is.na(class))   |>
  dplyr::mutate(class_ed = case_when(
    class %in% unique((mtags_table_ed_pk_rel_tax_m_red |> 
                         group_by(motu_num) |> 
                         filter(any(relative_abundance > 0.01)) |> 
                         pull(class))) ~ class,
    TRUE ~ 'others'))

mtags_table_ed_pk_rel_tax_m_red <- mtags_table_ed_pk_rel_tax_m_red |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                motu_num_f = as_factor(motu_num))

mtags_table_ed_pk_rel_tax_m_red$class_f <-  factor(mtags_table_ed_pk_rel_tax_m_red$class_f, 
                                        levels=unique(mtags_table_ed_pk_rel_tax_m_red$class_f[order(mtags_table_ed_pk_rel_tax_m_red$phylum_f)]), 
                                        ordered=TRUE)

mtags_table_ed_pk_rel_tax_m_red$order_f <-  factor(mtags_table_ed_pk_rel_tax_m_red$order_f, 
                                        levels=unique(mtags_table_ed_pk_rel_tax_m_red$order_f[order(mtags_table_ed_pk_rel_tax_m_red$phylum_f,
                                                                                         mtags_table_ed_pk_rel_tax_m_red$class_f)]), 
                                        ordered=TRUE)

mtags_table_ed_pk_rel_tax_m_red$family_f <-  factor(mtags_table_ed_pk_rel_tax_m_red$family_f, 
                                         levels=unique(mtags_table_ed_pk_rel_tax_m_red$family_f[order(mtags_table_ed_pk_rel_tax_m_red$phylum_f,
                                                                                           mtags_table_ed_pk_rel_tax_m_red$class_f,
                                                                                           mtags_table_ed_pk_rel_tax_m_red$order_f)]), 
                                         ordered=TRUE)

mtags_table_ed_pk_rel_tax_m_red$motu_num_f <-  factor(mtags_table_ed_pk_rel_tax_m_red$motu_num_f, 
                                          levels=unique(mtags_table_ed_pk_rel_tax_m_red$motu_num_f[order(mtags_table_ed_pk_rel_tax_m_red$phylum_f,
                                                                                             mtags_table_ed_pk_rel_tax_m_red$class_f,
                                                                                             mtags_table_ed_pk_rel_tax_m_red$order_f,
                                                                                             mtags_table_ed_pk_rel_tax_m_red$family_f)]), 
                                          ordered=TRUE)

mtags_table_ed_pk_community_class$class_ed <-  factor(mtags_table_ed_pk_community_class$class_ed, 
                                            levels=unique(mtags_table_ed_pk_community_class$class_ed[order(mtags_table_ed_pk_rel_tax_m_red$class)]), 
                                            ordered=TRUE)

mtags_table_ed_pk_community_class$class_ed <-  factor(mtags_table_ed_pk_community_class$class_ed, 
                                            levels=unique(mtags_table_ed_pk_community_class$class_ed[order(mtags_table_ed_pk_rel_tax_m_red$class)]), 
                                            ordered=TRUE)

mtags_table_ed_pk_community_class$class_ed |>
  unique()

mtags_table_ed_pk_community_class$class_ed <- factor(mtags_table_ed_pk_community_class$class_ed, 
                                                     levels = c(NA, 'others', 'Verrucomicrobiae', 'Gammaproteobacteria', 'Alphaproteobacteria', 'Cyanobacteriia' ,    
                                                               'Dehalococcoidia',     'Rhodothermia' ,       'Bacteroidia',         'Actinobacteria',
                                                               'Thermoplasmata',    'Unassigned',         ' Unaligned'    ))

mtags_table_ed_pk_community_class$class_ed |>
  unique()

mtags_table_ed_pk_community_class_plot <- mtags_table_ed_pk_community_class  |>
  dplyr::mutate(fraction = '0.2') |>
  group_by(class_ed, sample_id, fraction, date) |>
  dplyr::summarize(sum_rel = sum(relative_abundance)) |>
  dplyr::filter(!class_ed %in% c('Unassined', 'Unaligned')) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, sum_rel))+
  geom_area(aes(fill = class_ed, group = fct_rev(class_ed)), alpha = 1,  position='stack')+
  scale_fill_manual(values = palette_class_assigned_mtags)+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.5))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 10),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), strip.placement = 'outside')  

mtags_table_ed_pk_community_class_plot

# 
# ggsave(filename = 'mtags_table_ed_pk_community_class_plot_1perc_ed3.pdf', plot = mtags_table_ed_pk_community_class_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 288, height = 150, units = 'mm')

## with community Evenness ----
mtags_table_ed_pk_community_class_plot <- mtags_table_ed_pk_community_class  |>
  dplyr::mutate(fraction = '0.2') |>
  group_by(class_ed, sample_id, fraction, date) |>
  dplyr::summarize(sum_rel = sum(relative_abundance)) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(!class_ed %in% c('Unassined', 'Unaligned')) |>
  ggplot(aes(date, sum_rel))+
  geom_area(aes(fill = class_ed, group = fct_rev(class_ed)), alpha = 1,  position='stack')+
  scale_fill_manual(values = palette_class_assigned_mtags)+
  geom_line(data = community_eveness_02_red_mtags, aes(date, community_eveness_rar_mtags), linewidth = 2, color = '#DAE4E7')+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~. , name = 'Community Evenness'), limits = c(0,0.5))+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free_y',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 10),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), strip.placement = 'outside')  

mtags_table_ed_pk_community_class_plot

# ggsave(filename = 'mtags_table_ed_pk_community_class_plot_1perc_ed4.pdf', plot = mtags_table_ed_pk_community_class_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 288, height = 150, units = 'mm')

## NMDS
## NMDS general for the community -----
asv_tab_20y_w_rar_rel_w <- mtags_table_ed_pk_rel_tax_m_red |>
  dplyr::select(motu_num, sample_id, relative_abundance) |>
  #dplyr::filter(sample_id != 'BL231114') |> # seq depth higher than the rest for a check
  pivot_wider(id_cols = sample_id, names_from = 'motu_num', values_from = 'relative_abundance')

sample_id_col <-asv_tab_20y_w_rar_rel_w$sample_id

asv_tab_20y_w_rar_rel_w <- asv_tab_20y_w_rar_rel_w |>
  as.data.frame()

row.names(asv_tab_20y_w_rar_rel_w) <- asv_tab_20y_w_rar_rel_w[,1]  

asv_tab_20y_w_rar_rel_w <- asv_tab_20y_w_rar_rel_w[,-1]

data.hel <- asv_tab_20y_w_rar_rel_w |>
  decostand(method="hellinger"); str(data.hel)

data.dist <- vegdist(data.hel, method="bray")
head(data.dist)
data.nmds <- metaMDS(data.dist)                   # càlcul per poder col·locar a l'espai les comparacions entre comunitats
str(data.nmds)                                 # stress num 0.085 (per sota de 20; és acceptable)
data.nmds.points <- data.frame(data.nmds$points)  # convertir dades a data.frame per utilitzar amb qplot
plot(data.nmds.points)
head(data.nmds.points)

data.nmds.points |>
  colnames()

data.nmds.points |>
  dim() ## 64 samples

nmds_mtags <- data.nmds.points |>
  bind_cols(sample_id_col) |>
  dplyr::select(MDS1, MDS2, sample_id = '...3') |>
  tidyr::separate_wider_position(col = sample_id, widths = c('station' = 2, 'year' = 2, 'month' = 2, 'day' = 2), too_many = 'debug') |>
  dplyr::mutate(date = paste0(day, '-', month, '-', year)) |>
  left_join(m_mtags) |>
  dplyr::mutate(date = as.POSIXct(date, format = '%d-%m-%y'))

nmds_mtags$season <- factor(nmds_mtags$season, levels = c('winter', 'spring', 'summer', 'autumn'))

nmds_mtags_plot <- nmds_mtags |>
  dplyr::filter(sample_id != 'BL061107') |>
  ggplot(aes(MDS1, MDS2, color = season))+ # shape = fraction,
  geom_point()+
  #scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_seasons_4)+
  labs( color = 'Season')+
  #geom_text(aes(label = `sample_id`), check_overlap = F, nudge_x = 0.02, nudge_y = 0.04) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 16), panel.grid.minor = element_blank(),
        strip.text = element_text(size = 0),
        legend.position = 'bottom', axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16), strip.background = element_blank(), 
        legend.text = element_text(size = 16), legend.title = element_text(size = 16), 
        strip.placement = 'outside')

nmds_mtags_plot
 
# ggsave(filename = 'nmds_mtags_plot_ed1.pdf', plot = nmds_mtags_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 150, height = 150, units = 'mm')

## calculate rclr ----

# Transform data to rCLR ( FOR NOW I DO NOT PERFORM THIS NORMALIZATION FOR mTAGS )----
# ## transform asv_tab into wider format to go into vegan package
# asv_tab_bbmo_20y_w <- bbmo20y_tab_ed |>
#   as.data.frame()
# 
# rownames(asv_tab_bbmo_20y_w) <-  asv_tab_bbmo_20y_w[,1]
# 
# asv_tab_bbmo_20y_w <- asv_tab_bbmo_20y_w[,-1]
# 
# ## with the deconstant function from vegan using pseudocunt we don't lose samples but, in Coenen 2020 they say that
# ## adding a pseudocount disproportionately affects rare taxa, where the magnitude of differences between samples may 
# ## be similar to the magnitude of the added pseudocount and therefore obscured.
# 
# ## However, when using the rCLR which is similar to the CLR it allows data that contains zeroes. This method does not use
# ## pseudocounts, unlike the standard CLR. Robust clr divides the vales by geometric mean o the observed features; zero values 
# ## are kept as zeroes, and not taken into account. In high dimensional data, the geometric mean of rclr is a good approximation
# ## of the true geometric mean (from deconstand documentation)
# 
# ## The centered log-ratio (CLR) transformation is a CoDa approach that uses the geometric mean of the read counts of all taxa 
# ## within a sample as the reference/denominator for that sample. In this approach all taxon read counts within a sample are 
# ## divided by this geometric mean and the log fold changes in this ratio between samples are compared.
# 
# #### MAYBE IT SHOULD BE TRANSPOSED SO THAT SAMPLES ARE IN ROWS AND ASVS IN COLUMNS ??? (WE NORMALIZE BY THE GEOMETRIC MEAN OF READS IN A SAMPLE)
# ### after checking the results from the different approximations we keep the rCLR transformation.
# asv_tab_bbmo_20y_w_t <- asv_tab_bbmo_20y_w |>
#   t()
# 
# # asv_tab_bbmo_20y_w |>
# #   colnames()
# # 
# # asv_tab_bbmo_20y_w |>
# #   as_tibble() |>
# #   dplyr::filter(any() < 0)
# 
# rclr_df <- decostand(asv_tab_bbmo_20y_w_t, method = 'rclr' )
# #clr_df <- decostand(asv_tab_bbmo_10y_w, method = 'clr' , pseudocount = 1)
# 
# rclr_df_sample_id <- rclr_df |>
#   row.names() |>
#   as_tibble_col(column_name = 'sample_id')
# 
# ##we create two datasets one for FL and one for PA
# asv_tab_20y_rclr <- rclr_df |>
#   as_data_frame() |>
#   bind_cols(rclr_df_sample_id ) |>
#   #rownames_to_column(var = 'motu_num') |>
#   pivot_longer(cols = starts_with('ASV'), names_to = 'motu_num', values_to = 'rclr')
# 
# bbmp_20y_rel_rclr_tax <- mtags_table_ed_pk_rel_tax_m_ed |>
#   left_join(asv_tab_20y_rclr, by = c('motu_num', 'sample_id'))
# 
# bbmp_20y_rel_rclr_tax  |>
#   colnames()
# 
# bbmp_20y_rel_rclr_tax |>
#   dplyr::filter(motu_num == 'ASV1') |>
#   ggplot(aes(rclr, relative_abundance))+
#   geom_point()
# 
# bbmp_20y_rel_rclr_tax |>
#   dplyr::filter(motu_num == 'ASV1') |>
#   ggplot(aes(date, relative_abundance))+
#   geom_point()
# 
# bbmp_20y_rel_rclr_tax |>
#   dplyr::filter(motu_num == 'ASV1') |>
#   ggplot(aes(date, rclr))+
#   geom_point()


## find blooms -----
# Discover anomalies ----
## For each ASVs based on relative abundances (but i also calculate them for rCLR transformation and pseudoabundances to check)-----
#x <- 120*0.75 ##percentage of ASV present in the dataset that we want to subset by (occurrence)

mtags_table_ed_pk_rel_tax_m_red

z_02 <- mtags_table_ed_pk_rel_tax_m_red |>
  group_by(motu_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.025)) |> #estan en format 0-1
  #group_by(motu_num) |>
  dplyr::reframe(#anomalies_ps = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1.96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
    anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)])

## I only filter for those anomalies in relative abundance because pseudoabundance I can only use it for FL not for PA and rclr will be used to explore if the changes observed are true

# Filter the ASV_tab by only those ASVs that have an anomaly at some point of the dataset (IN THIS CASE WE ARE NOT FILTERING BY ANY RELATIVE ABUNDANCE IN THE COMMUITY JUST IF THEY PRESENT AN ANOMALY AT SOME POINT!)----

asv_anom_02 <- find_asv_with_anomalies(anomalies_result = z_02, 
                                       anomaly_in1 = anomalies_ra, 
                                       anomaly_in2 = NULL,  #anomalies_ps,
                                       anomaly_in3 = NULL, #anomalies_clr, 
                                       logic1 = 'TRUE',
                                       logic2 = NULL, 
                                       logic3 = NULL,
                                       asv_col = motu_num)
## 4 bloomers

## Recover z-score for each ASV ----
### I want to highlight anomalies for each ASV to do so I recover z-scores for those ASVs that that have high z-scores
### at some point of the dataset. Easy to observe if those ASVs are having random anomalies or all of them happen at the same time
#x <- 120*0.75
m_mtags <- m_mtags |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_mtags))) 

m_mtags_red <- m_mtags |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) 

m_mtags_red <- m_mtags |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_mtags_red)))
  
z_scores_02 <- mtags_table_ed_pk_rel_tax_m_red |>
  group_by(motu_num) |>
  dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 50% of the samples
  #dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
  group_by(motu_num) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra) |>
  group_by(motu_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_mtags_red))) |>
  left_join(m_mtags_red, by = 'sample_id_num')

## create datasets with the information from 'bloomers' in the mtags timeseries ----
z_scores_02_red <- z_scores_02 |>
  dplyr::select(-sample_id_ok, -sample_id_remainder, -sample_id_width)

mtags_table_ed_pk_rel_tax_m_red_z_scores <- mtags_table_ed_pk_rel_tax_m_red |>
  dplyr::select(-sample_id_ok, -sample_id_remainder, -sample_id_width) |>
  left_join(z_scores_02_red)

mtags_bloo_tax <- asv_anom_02 |>
  as_tibble_col(column_name = 'motu_num') |>
  left_join(mtags_tax)

#write.csv(mtags_bloo_tax , 'data/mtags_bloo_tax.csv')
#write.csv(mtags_table_ed_pk_rel_tax_m_red_z_scores, 'data/mtags_table_ed_pk_rel_tax_m_red_z_scores.csv')

### Explore what happens if we do not use threshold for abundance at least at some point of the dataset, justify my threshold at 10%-----
# Create a function to generate datasets based on the number of ASVs that fulfill the criteria of potential bloomers at different relative abundance thresholds.
source('src/count_number_potential_bloomers_threshold.R')

mtags_table_ed_pk_rel_tax_m_red |>
  colnames()

mtags_table_ed_pk_rel_tax_m_red_ed <- mtags_table_ed_pk_rel_tax_m_red |>
  dplyr::select(everything(), asv_num = motu_num)

z_scores_02_ed <- z_scores_02 |>
  dplyr::select(everything(), asv_num = motus_num)

# I apply a loop for all the thresholds that I'm interested in
# Define a vector of threshold values
threshold_values <- c(0, 0.001, 0.0015, 0.0025, 0.005, 0.0075, 
                      0.01, 0.015, 0.0, 0.025, 0.05, 0.075,
                      0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6)

# Create an empty list to store the results
datasets <- list()

# Iterate over each threshold value and apply the function
for (threshold in threshold_values) {
  dataset_name <- paste0("n_0.2_", threshold * 100)  # Create dataset name
  dataset <- count_num_bloomers(threshold, 0.2, mtags_table_ed_pk_rel_tax_m_red_ed, z_scores_tb = z_scores_02_ed)
  datasets[[dataset_name]] <- dataset  # Store the dataset in the list
}

# Combine all datasets into a single dataframe
result_dataset_02 <- bind_rows(datasets)

### plot the results
blooming_threshold <- result_dataset_02 |>
  ggplot(aes(threshold, num))+
  geom_point(size = 1)+
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = percent_format())+
  scale_y_log10(expand = c(0,0))+
  #scale_y_log10()+
  facet_grid(vars(fraction), labeller = labs_fraction)+
  geom_vline(xintercept = 0.025, linetype = 'dashed')+
  labs(x = 'Relative abundance (%) threshold of the potential blooming event', y = 'Number of potential bloomers detected')+
  geom_line()+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 5),
        axis.ticks = element_blank())

# ggsave(blooming_threshold,  filename = 'blooming_threshold_mTAGs.pdf',
#        path = 'Results/Figures/',
#        width = 88, height = 88, units = 'mm')

## plot the 'bloomers' in the mtags time series -----
mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo <- mtags_table_ed_pk_rel_tax_m_red_z_scores |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(motu_num %in% mtags_bloo_tax$motu_num)

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo  <- mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo  |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                motu_num_f = as_factor(motu_num))

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f <-  factor(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f, 
                                                       levels=unique(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f[order(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$phylum_f)]), 
                                                       ordered=TRUE)

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$order_f <-  factor(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo $order_f, 
                                                       levels=unique(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo $order_f[order(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$phylum_f,
                                                                                                                        mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f)]), 
                                                       ordered=TRUE)

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$family_f <- factor(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$family_f, 
                                                       levels=unique(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$family_f[order(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$phylum_f,
                                                                                                                        mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f,
                                                                                                                        mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$order_f)]), 
                                                       ordered=TRUE)

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$motu_num_f <-  factor(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$motu_num_f, 
                                                         levels=unique(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$motu_num_f[order(mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$phylum_f,
                                                                                                                           mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f,
                                                                                                                           mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$order_f,
                                                                                                                           mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo$family_f)]), 
                                                         ordered=TRUE)

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot <- mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  dplyr::filter(motu_num != 'Unassigned') |> # it does not behave as a bloomer
  dplyr::filter(motu_num != 'Unaligned') |> # it does not behave as a bloomer
  ggplot(aes(date, relative_abundance))+
  #geom_line(aes(group = motu_num))+
  geom_area(aes(fill = order_f))+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  scale_fill_manual(values = palette_order_assigned_all)+
  facet_wrap(vars(interaction(family_f, motu_num)), ncol = 1, scales = 'free_y')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12), strip.background = element_blank(), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12), strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 2, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot
 
# ggsave(filename = 'mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_ed2.pdf', plot = mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 180, height = 250, units = 'mm')

## all together with community evenness ----
mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot <- mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  #dplyr::filter(motu_num != 'ASV1') |> # it does not behave as a bloomer
  dplyr::filter(motu_num != 'Unassigned') |> # it does not behave as a bloomer
  dplyr::filter(motu_num != 'Unaligned') |> # it does not behave as a bloomer
  group_by(date, sample_id) |>
  dplyr::mutate(abund_max = sum(relative_abundance)) |>
  dplyr::group_by(order, date, abund_max) |>
  dplyr::reframe(abund_order = sum(relative_abundance)) |>
  ungroup() |>
  #dplyr::distinct(date, sample_id, order, relative_abundance, abund_max) |>
  #dplyr::left_join(community_eveness_02) |>
  ggplot(aes(date, abund_max))+
  #geom_line(aes( y = abund_max))+
  geom_area(aes(fill = order, y = abund_order), position = 'stack')+
  geom_line(data = community_eveness_02_red, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2004-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot
# # 
# ggsave(filename = 'mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot_ed1.pdf', plot = mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot,
#        path = 'results/figures/',
#        width = 300, height = 180, units = 'mm')

## all together with community evenness ----
mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot <- mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  #dplyr::filter(motu_num != 'ASV1') |> # it does not behave as a bloomer
  dplyr::filter(motu_num != 'Unassigned') |> # it does not behave as a bloomer
  dplyr::filter(motu_num != 'Unaligned') |> # it does not behave as a bloomer
  group_by(date, sample_id) |>
  dplyr::mutate(abund_max = sum(relative_abundance)) |>
  dplyr::group_by(order, date, abund_max) |>
  dplyr::reframe(abund_order = sum(relative_abundance)) |>
  ungroup() |>
  #dplyr::distinct(date, sample_id, order, relative_abundance, abund_max) |>
  #dplyr::left_join(community_eveness_02) |>
  ggplot(aes(date, abund_max))+
  #geom_line(aes( y = abund_max))+
  geom_area(aes(fill = order, y = abund_order), position = 'stack')+
  geom_line(data = community_eveness_02_red, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2008-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot
# # 
# ggsave(filename = 'mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot_ed2.pdf', plot = mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot,
#        path = 'results/figures/',
#        width = 300, height = 180, units = 'mm')

# Calculate diversity parameters ----
## We calculate different diversity parameters to check if the blooming events detected have an effect on the community structure
## Community Evenness----
### We need to apply rarefaction because we have uneven sequencing effort, and this will affect community diversity analysis
### Rarefaction (process that repeats the subsampling many times and tries to overcome the uneven sequencing effort bias)----
 #wider format to go into vegan package 

#### calculate the minimum read size of all samples to rarefy at that number
min_n_seqs <- mtags_table_ed_pk |>
  #rownames_to_column(var = 'motu_num') |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'reads') |>
  group_by(sample_id) |>
  dplyr::summarize(n_seqs = sum(reads)) |>
  dplyr::summarize(min = min(n_seqs)) |>
  pull(min)

# mtags_table_ed_pk_t <- mtags_table_ed_pk |>
#   t() # crec que no ho necessitem en aquest cas

mtags_table_ed_pk_ed <- mtags_table_ed_pk |>
  dplyr::select(-motu_num) |>
  t()

mtags_table_ed_pk_motu_num <- mtags_table_ed_pk |>
  dplyr::select(motu_num) |>
  add_row( .before = 1 ) |>
  dplyr::mutate(motu_num = ifelse(is.na(motu_num), 'sample_id', motu_num))

#this part is commented since is computationally slow so i don't want to repeat this step many times. 
##this function gives us a randomly rarefied community data
# rrarefy(mtags_table_ed_pk_ed, sample = min_n_seqs) |>
#   as_tibble(rownames = 'sample_id') #just rarefying (one simgle subsampling)

## we use rarefaction (which repeats the subsampling step many times)
## perform this a 1000 times to get an empirical diversity values calculating the mean value for each timepoint.
# mtags_table_ed_pk_rar <- rrarefy.perm(mtags_table_ed_pk_ed, # samples in rows!
#                                        sample = min_n_seqs,
#                                        n = 1000,
#                                        round.out = T)

#write.csv2(mtags_table_ed_pk_rar, file = 'data/mtags_table_ed_pk_rar.csv')

mtags_table_ed_pk_rar <- read.csv2('data/mtags_table_ed_pk_rar.csv') |>
  as_tibble()|> 
  rename('sample_id' = 'X') 

mtags_table_ed_pk_rar|>
  colnames() <- mtags_table_ed_pk_motu_num$motu_num

## Rarefied dataset to calculate Community eveness----
source('../../Bloomers/R/community_evenness.R')

community_eveness_02_mtags <- mtags_table_ed_pk_rar |>
  #as.data.frame() |>
  #tibble::rownames_to_column(var = 'sample_id') |>
  dplyr::select(-Unassigned, -Unaligned) |>
  pivot_longer(cols = starts_with(c('mtags')), names_to = 'motu_num', values_to = 'reads_rar') |>
  #dplyr::select(sample_id, reads, motu_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou')) |>
  dplyr::mutate(sample_id = str_replace(sample_id, '.bins', ''))

community_eveness_02_mtags <- community_eveness_02_mtags |>
  left_join(m_mtags, by = 'sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y")))

community_eveness_02_red_mtags <- community_eveness_02_mtags |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d")))

### plot community evenness----
community_eveness_02_red_mtags |>
  #left_join(m_mtags, by = 'sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, community_eveness_rar))+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E')+
  geom_point()+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Community Eveness')+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank())

## Bray Curtis dissimilarity----
### 0 means the two sites have the same composition (that is they share all the species), and 1 means the two sites 
### do not share any species.
source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R') #we need to upload it since we updated the function but I didn't compile the package
# 
# bray_curtis_02 <- dissimilarity_matrix(data = asv_tab_10y_02_rel, 
#                                            sample_id_col = sample_id,
#                                            values_cols_prefix = 'BL')
# 
# bray_curtis_3 <- dissimilarity_matrix(data = asv_tab_10y_3_rel, 
#                                           sample_id_col = sample_id,
#                                           values_cols_prefix = 'BL')

### We need the rarefied table transformed to relative abundances
mtags_table_ed_pk_rel_rar <- mtags_table_ed_pk_rar |>
  # as.data.frame() |>
  # rownames_to_column(var = 'sample_id') |>
  dplyr::select(-Unassigned, -Unaligned) |>
  pivot_longer(cols = starts_with('mtags'), values_to = 'reads', names_to = 'motu_num') |>
  dplyr::mutate(reads = as.numeric(reads)) |>
  calculate_rel_abund(group_cols = sample_id) |>
  dplyr::mutate(sample_id = str_replace(sample_id, '.bins', ''))

mtags_table_ed_pk_rel_rar <- mtags_table_ed_pk_rel_rar %>%
  dplyr::filter(sample_id %in% community_eveness_02_red_mtags$sample_id)

mtags_table_ed_pk_rel_rar_ed <- mtags_table_ed_pk_rel_rar |>
  dplyr::select(sample_id, asv_num = motu_num,  reads, total_reads, relative_abundance)

bray_curtis_02_rar_mtags <- dissimilarity_matrix(data = mtags_table_ed_pk_rel_rar_ed, 
                                           sample_id_col = sample_id,
                                           values_cols_prefix = 'BL')

### plot bray curtis dissimilarity----
#### compare rarified dataset with non-rarified dataset (NO CHANGES, WHY?)
# bray_curtis_02 |>
#   bind_rows(bray_curtis_3) |>
#   dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
#   left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   ggplot(aes(date, bray_curtis_result))+
#   geom_point()+
#   facet_grid(vars(fraction))+
#   scale_x_datetime()+
#   labs(x = 'Time', y = 'Bray Curtis Dissimilarity')+
#   geom_line()+
#   theme_bw()+
#   theme(panel.grid = element_blank(), strip.background = element_blank())

bray_curtis_02_rar_plot <- bray_curtis_02_rar_mtags |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  left_join(m_mtags, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, bray_curtis_result))+
  geom_point()+
  geom_line(aes(date, bray_curtis_result))+
  #facet_grid(vars(fraction))+
  #scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  #geom_smooth(method = 'lm')+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity', color = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

bray_curtis_02_rar_plot

### plot Bray-Curtis dissimilarity and Community Evenness together----
bray_curtis_02_rar_mtags_red <-  bray_curtis_02_rar_mtags |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  left_join(m_mtags, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d")))

community_diversity_plot <- community_eveness_02_red_mtags |> 
  left_join(bray_curtis_02_rar_mtags_red) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  #left_join(m_mtags, by = c('sample_id')) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, value))+
  #geom_point(aes())+
  geom_line(aes(date, value))+
  geom_smooth(method = 'loess', span = 0.06, color = 'black')+
  geom_smooth(method = 'lm', color = 'darkgreen')+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  facet_grid(vars(diversity_index), labeller = labs_diversity)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Time', y = 'Community diversity', color = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom',
        axis.text.x = element_text(size = 6))

community_diversity_plot

# ggsave(filename = 'community_diversity_mtags_plot_ed2.pdf', plot = community_diversity_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 180, height = 180, units = 'mm')

## Explore the threshold that we decide to use for this datset ----
### Explore what happens if we do not use threshold for abundance at least at some point of the dataset, justify my threshold at 10%-----
# Create a function to generate datasets based on the number of ASVs that fulfill the criteria of potential bloomers at different relative abundance thresholds.
source('src/count_number_potential_bloomers_threshold.R')

# I apply a loop for all the thresholds that I'm interested in
# Define a vector of threshold values
threshold_values <- c(0, 0.001, 0.0015, 0.0025, 0.005, 0.0075, 
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
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = percent_format())+
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
        text = element_text(size = 5),
        axis.ticks = element_blank())

# ggsave(blooming_threshold,  filename = 'blooming_threshold_mtags.pdf',
#        path = 'Results/Figures/',
#        width = 88, height = 88, units = 'mm')

## Comparison mTAGs community Evenness with ASVs community Evenness in the FL fraction
community_eveness_02_red_motus

community_eveness_all_m_red <- community_eveness_all_m |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::select(sample_id_ed, community_eveness_rar_asv = community_eveness_rar, date)

## edit the dates since the are some dates that are different fom mOTUs and ASVs (one day difference)
community_eveness_all_m_red <- community_eveness_all_m_red |>
  dplyr::mutate(sample_id_ed = case_when(sample_id_ed == 'BL110705' ~ 'BL110704',
                                         sample_id_ed ==  'BL100414' ~ 'BL100413',
                                         sample_id_ed ==  'BL120511' ~ 'BL120518',
                                         sample_id_ed ==  'BL131215' ~ 'BL131204',
                                         sample_id_ed == 'BL080312' ~ 'BL080311',
                                         T ~ sample_id_ed))

comparison_community_evenness_asvs_motus_tb <- community_eveness_02_red_motus |>
  dplyr::select(-station, -year, -month, -day, -sample_id_ok, -sample_id_remainder, -sample_id_width) |>
  left_join(community_eveness_all_m_red, by = c('sample_id' = 'sample_id_ed'))

comparison_community_evenness_asvs_motus_tb |>
  dplyr::filter(is.na(community_eveness_rar_asv))

# sample_id community_eveness_rar date                season sample_id_ed community_eveness_rar_asv
# <chr>                     <dbl> <dttm>              <chr>  <chr>                            <dbl>
#   1 BL100413              0.8258590 2010-04-13 00:00:00 spring NA                                  NA
# 2 BL110704              0.7983031 2011-07-04 00:00:00 summer NA                                  NA
# 3 BL120518              0.6883892 2012-05-18 00:00:00 spring NA                                  NA
# 4 BL131204              0.8384126 2013-12-04 00:00:00 autumn NA                                  NA

# community_eveness_all_m_red |>
#   dplyr::filter(str_detect(date, '2008'))

correlation_community_evenness_ASvs_mOTUs_plot <- comparison_community_evenness_asvs_motus_tb |>
  ggplot(aes(community_eveness_rar_asv, community_eveness_rar))+
  geom_point(aes(), color = '#344848')+
  labs(x = 'ASVs Community Evenness', y = 'mTAGs Community Evenness')+
  geom_smooth(method = 'lm', color = '#344848')+
  geom_text(aes(label = `sample_id`), check_overlap = F, nudge_x = 0.02, nudge_y = 0.04)+
  theme_bw()+
  theme(strip.background = element_blank(),
        aspect.ratio = 4/4,
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.position = 'bottom')

correlation_community_evenness_ASvs_mOTUs_plot  

# ggsave(correlation_community_evenness_ASvs_mOTUs_plot,  filename = 'correlation_community_evenness_ASvs_mOTUs_plot_labs_outliers.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 188, units = 'mm')

comparison_community_evenness_mTAGs_ASVs_plot <- comparison_community_evenness_asvs_motus_tb |>
  pivot_longer(cols = starts_with('community')) |>
  dplyr::mutate(community_evenness = case_when(!str_detect(name, 'asv') ~ 'mTAGs',
                                               str_detect(name, 'asv') ~ 'ASVs')) |>
  dplyr::mutate(date.x = (as.POSIXct(date.x, format = "%d-%m-%y")))   |>
  ggplot(aes(date.x, value, color = community_evenness, fill = community_evenness))+
  facet_wrap(vars(community_evenness), ncol = 1)+
  scale_color_manual(values = palette_motus_asvs)+
  scale_fill_manual(values = palette_motus_asvs)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(aes(group = community_evenness))+
  labs(x = 'Time', y = 'Community Evenness', fill = 'Community data', color = 'Community data' )+
  geom_smooth(method = 'loess', span = 0.07)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.position = 'bottom')

comparison_community_evenness_mTAGs_ASVs_plot              
               
# ggsave(comparison_community_evenness_mTAGs_ASVs_plot,  filename = 'comparison_community_evenness_mTAGs_ASVs_plot.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 128, units = 'mm')

## correlation ASV11 and mOTUs number ----
motus14723 <- motus_table_ed_pk_rel_tax_m |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  #dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |> 
  dplyr::filter(motu_num == 'motus14723') |>
  dplyr::mutate(sample_id_ed = case_when(sample_id == 'BL110705' ~ 'BL110704',
                                         sample_id ==  'BL100414' ~ 'BL100413',
                                         sample_id ==  'BL120511' ~ 'BL120518',
                                         sample_id ==  'BL131215' ~ 'BL131204',
                                         sample_id == 'BL080312' ~ 'BL080311',
                                         T ~ sample_id)) |>
  ungroup() |>
  dplyr::select(sample_id_ed, date, relative_abundance, variable = motu_num)

asv11 <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new) |>
  ungroup() |>
  #dplyr::distinct(order, family, class)
  left_join(m_bbmo_10y) |>
  dplyr::filter(asv_num == 'asv11') |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::select(sample_id_ed, date, relative_abundance, variable =  asv_num)

motus14723_vs_asv11 <- motus14723 |>
  bind_rows(asv11)

motus_tax |>
  dplyr::filter(motu_num == 'motus14723') %$%
  tax

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num == 'asv11') %$%
  seq |>
  unique()

##before correlation we check normality
shapiro.test(as.numeric(motus14723$relative_abundance)) # =>p-value = 1.884e-15 (NO NORMALITY)
ggqqplot(as.numeric(motus14723$relative_abundance))

shapiro.test(as.numeric(asv11$relative_abundance)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(asv11$relative_abundance))

motus14723_vs_asv11_plot <- motus14723_vs_asv11 |>
  pivot_wider(id_cols = c('sample_id_ed', 'date'), names_from = 'variable', values_from = 'relative_abundance') |>
  ggplot(aes(motus14723, asv11))+
  geom_point()+
 # geom_text(aes(label = `sample_id_ed`), check_overlap = T, nudge_x = -0.02, nudge_y = 0.02)+
  labs(x = 'mTAG 14723', y = 'ASV 11')+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.1,
    #label.y = -0.02,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside', aspect.ratio = 4/4)

motus14723_vs_asv11_plot

ggsave(motus14723_vs_asv11_plot,  filename = 'motus14723_vs_asv11_plot.pdf',
       path = 'Results/Figures/',
       width = 88, height = 88, units = 'mm')

## timeseries ----
timeseries_motus14723_vs_asv11_plot <- motus14723_vs_asv11 |>
  ggplot(aes(date, relative_abundance))+
  geom_line(aes(group = variable, color = variable))+
  #facet_wrap(vars(community_evenness), ncol = 1)+
  scale_color_manual(values = palette_long)+
  scale_fill_manual(values = palette_long)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  #geom_line(aes(group = community_evenness))+
  labs(x = 'Time', y = 'Relative abundance', fill = '', color = '' )+
  geom_smooth(method = 'loess', span = 0.089, aes(group = variable, color = variable, fill = variable))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        axis.ticks = element_blank(),
        legend.position = 'bottom')

timeseries_motus14723_vs_asv11_plot 

# ggsave(timeseries_motus14723_vs_asv11_plot,  filename = 'timeseries_motus14723_vs_asv11_plot_smooth.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 88, units = 'mm')

## Blastn ASV11 with mTAG -----

# >silva_138_complink_cons_otu_93644
# AGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACASGATASCTTGCTA
# TCGCTGACGAGTGGCGGACGGGTGAGTAABACTTASGGATCTGCCTCTGTGTGGGGGATAACTATTGGAAACGATAGCTA
# ATACCGCATAATGTCTMCRGACCAARGMGGGCTTYRGCTCKYGCGCAGAGAGGAACCTAAGCGAGATTAGCTAGTTGGTG
# AGGTAAAGGCTCACCAAGGCGACGATCTCTAGCTGTTCTGARAGGAAGATCAGCCMCMCTGGAACTGARACMCGGTCCAG
# ACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCCMTGCCGCGTGTGTGAAGAAGG
# SCTTCGGGTTGTAAAGCACTTTCAGTTGTGAGGAAAGTTTAGTAGTTAATACCTGCTAGATGTGACGTTAGCAACAGAAG
# AAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCG
# CACGCARRCGGTCTGTTAAGCTAGATGTGAAAGCCCCGCGCTCAACGTGKGAGGGTCATTTAGAACTGGCAGACTAGAGT
# CTTGGAGAGGGGAGTGGAATTCCASGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAACATCAATGGCGAAGGCAACTC
# CCTGGCCAAAGACTGACGCTCATGTGCGAAAGTGTGGGTAGCGAACAGGATTAGATACCCTGGTAGTCCACACCGTAAAC
# GCTGTCTACTAGCTGTTTGTGRATTTAATYCGTGAGTAGCGMAGCTAACGCGATAAGTAGACCGCCTGGGGAGTACGGCC
# GCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGA
# ACCTTACCTACTCTTGACATACTAGAAACTTTTCAGAGATGAATTGGTGCCTTCGGGAATCTAGATACAGGTGCTGCATG
# GCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCCTTAGTTGCCAGCCTTAA
# GTTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGACGACGTCAAGTCATCATGGCCCTTACGAGT
# AGGGMTACACACGTGCTACAATGGMKRGTACAAAGGGATGCRARCCTGCGAAGGTAAGCGGACCCCTTAAAGCYMKTCGT
# AGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTAGGTCAGCATACTACGGTGAATAC
# GTTCCCGGGCCTTGTACACRCCGCCCGTCACACCAYGGGASTGGGATGCAAAAGAAGTAGGTAGCTTAACCTTCGGGATG
# GCGCTTACCACTTTGTGTTTCATGACTGGGGTGAAGTCGTAACAAGG

## Genome of Glaciecola BL_1001_bin.full.234 (RAMIRO) that could belong to the ASV11.

## ASV11 seq 
# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(asv_num == 'asv11') |>
#   distinct(seq)

## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGTTGTGAGGAAAGTTTAGTAGTTAATACCTGCTAGATGTGACGTTAGCAACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTTAAGCTAGATGTGAAAGCCCCGCGCTCAACGTGGGAGGGTCATTTAGAACTGGCAGACTAGAGTCTTGGAGAGGGGAGTGGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAACATCAATGGCGAAGGCAACTCCCTGGCCAAAGACTGACGCTCATGTGCGAAAGTGTGGGTAGCGAACAGG

## Blast-n result 98% of identity. The aligment is in results/reports/glaciecola_asv11




# -------------------- ################## mOTUs ##################  -------

## mOTUs comparison (mOTUs are regions that map not only 16S rRNA) -----
### Fran Latorre data ----
motus_tb <- read.delim2('data/genes_sergio/BBMO_all_sample_profiles_counts.txt', skip = 2,  sep = '\t') |>
  as_tibble()

motus_tb <- motus_tb |>
  rename('tax' = 'X.consensus_taxonomy')

motus_tb |>
  dim()

## organize tax
### upload reference database ---- 
ref_dtbs_motus <- read_tsv('data/genes_sergio/mOTUs_3.1.0_GTDB_tax.tsv', col_names = F) 

ref_dtbs_motus |>
  colnames() <- c('motus_num_ref', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')

ref_dtbs_motus <- ref_dtbs_motus |> 
  dplyr::mutate(domain = str_replace(domain, 'd__','')) |>
  dplyr::mutate(phylum = str_replace(phylum, 'p__','')) |>
  dplyr::mutate(class = str_replace(class, 'c__','')) |>
  dplyr::mutate(order = str_replace(order, 'o__','')) |>
  dplyr::mutate(family = str_replace(family, 'f__','')) |>
  dplyr::mutate(genus = str_replace(genus, 'g__','')) |>
  dplyr::mutate(species = str_replace(species, 's__','')) |>
  distinct(motus_num_ref, domain, phylum, class, order, family, genus, species) |>
  dplyr::filter(!str_detect(motus_num_ref, 'no_mOTU'))
  
motus_tax <- motus_tb |>
  dplyr::select(tax)

motus_tax_unassigned <- motus_tax |>
  separate(tax, into = c('tax_red', 'motus_num_ref'), sep = '\\[', remove = F) |>
  dplyr::filter(str_detect(tax, 'unassigned')) |>
  dplyr::mutate(motus_num_ref = case_when(is.na(motus_num_ref) ~ 'unassigned'))

motus_tax |>
  dim()

motus_tax <- motus_tax[-34341,]

motus_tax_ed <- motus_tax |>
  separate(tax, into = c('tax_red', 'motus_num_ref'), sep = '\\[', remove = F) |>
  dplyr::mutate(motus_num_ref = str_replace(motus_num_ref, ']', '')) |>
  bind_rows(motus_tax_unassigned)
  
motus_tax_ed_complete <- motus_tax_ed |>
  dplyr::mutate(motus_num = paste0('motus_', 1:nrow(motus_tax_ed))) |>
  left_join(ref_dtbs_motus, by = 'motus_num_ref') |>
  #dplyr::filter(tax == 'unassigned') |>
  dplyr::mutate(motus_num = case_when(tax == 'unassigned' ~ 'unassigned',
                                      TRUE ~ motus_num),
                domain = case_when(tax == 'unassigned' ~ 'unassigned',
                                   TRUE ~ domain),
                phylum = case_when(tax == 'unassigned' ~ 'unassigned',
                                   TRUE ~ phylum),
                class = case_when(tax == 'unassigned' ~ 'unassigned',
                                  TRUE ~ class),
                order = case_when(tax == 'unassigned' ~ 'unassigned',
                                  TRUE ~ order), 
                family = case_when(tax == 'unassigned' ~ 'unassigned',
                                   TRUE ~ family), 
                genus = case_when(tax == 'unassigned' ~ 'unassigned',
                                  TRUE ~ genus),
                species = case_when(tax == 'unassigned' ~ 'unassigned',
                                  TRUE ~ species)) 

motus_tax_ed_complete  |>
  dplyr::filter(domain == 'Eukaryota')  #0 mOTUs  

motus_tax_ed_complete_unique <- motus_tax_ed_complete |>
  distinct(tax, tax_red, motus_num_ref, domain, phylum, class, order, motus_num) # in this type of taxonomy there are some duplicated taxonomies that differe at species or genus level which can add some errors. That is why we are workiing with tax (tax consensus)
  
motus_tax_ed_complete_unique |>
  dplyr::filter(tax == 'unassigned')

## calculate relative abundances
motus_tb_l <- motus_tb |>
  dplyr::mutate(motus_num = paste0('motus_', 1:nrow(motus_tax_ed))) |>
  dplyr::mutate(motus_num = case_when(tax == 'unassigned' ~ 'unassigned',
                                      TRUE ~ motus_num)) |>
  dplyr::select(-tax) |>
  pivot_longer(cols = starts_with(c('BL', 'unassigned')), names_to = 'sample_id', values_to = 'counts')

motus_tb_rel <- motus_tb_l |>
  rename('asv_num' = 'motus_num', 'reads' = 'counts') |>
  calculate_rel_abund(group_cols = sample_id) |>
  rename('motus_num' = 'asv_num') |>
  tidyr::separate_wider_position(col = sample_id, widths = c('station' = 2, 'year' = 2, 'month' = 2, 'day' = 2), too_many = 'debug') |>
  dplyr::mutate(date = paste0(day, '-', month, '-', year)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y")))

motus_tb_rel |>
  group_by(sample_id) |>
  reframe(total = sum(relative_abundance)) |>
  dplyr::filter(total != '1')

# motus_tb_rel |>
#   ungroup() |>
#   distinct(tax)
# 
# motus_tb_rel |>
#   colnames()

identify_errors <- motus_tb_rel %$%
  motus_num |>
  unique() == motus_tax_ed_complete_unique  %$%
  motus_num |>
  unique()

identify_errors |>
  summary()

motus_rel_tax <- motus_tb_rel |>
  dplyr::select(-sample_id_width, -sample_id_remainder, -sample_id_ok, -reads, -total_reads) |>
  ungroup() |>
  #distinct(tax) |>
  left_join(motus_tax_ed_complete_unique, by = 'motus_num')

motus_tax_ed_complete_unique |>
  distinct(motus_num) |>
  summary()

motus_tb_rel |>
  ungroup() |>
  distinct(motus_num) |>
  summary()

# motus_tax_ed_complete_unique |>
#   distinct(tax)
# 
# motus_rel_tax |>
#   ungroup() |>
#   distinct(tax)

motus_rel_tax |>
  group_by(sample_id) |>
  reframe(total = sum(relative_abundance)) |>
  dplyr::filter(total != '1')

## colored only those taxa that represent > 1% at some point -----
motus_table_ed_pk_abund_class <- motus_rel_tax |>
  #dplyr::filter(!is.na(class))   |>
  dplyr::mutate(class_ed = case_when(
    class %in% unique((motus_rel_tax |> 
                         group_by(motus_num) |> 
                         filter(any(relative_abundance > 0.01)) |> 
                         pull(class))) ~ class,
    TRUE ~ 'others'))

motus_rel_tax |>
  colnames()

motus_rel_tax <- motus_rel_tax  |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                #family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class))

motus_rel_tax$class_f <-  factor(motus_rel_tax$class_f,
                                                   levels=unique(motus_rel_tax$class_f[order(motus_rel_tax$phylum_f)]),
                                                   ordered=TRUE)

motus_rel_tax$order_f <-  factor(motus_rel_tax$order_f,
                                                   levels=unique(motus_rel_tax$order_f[order(motus_rel_tax$phylum_f,
                                                                                                               motus_rel_tax$class_f)]),
                                                   ordered=TRUE)

# motus_rel_tax$family_f <-  factor(motus_rel_tax$family_f,
#                                                     levels=unique(motus_rel_tax$family_f[order(motus_rel_tax$phylum_f,
#                                                                                                                  motus_rel_tax$class_f,
#                                                                                                                  motus_rel_tax$order_f)]),
#                                                     ordered=TRUE)

# motus_rel_tax$motu_num_f <-  factor(motus_rel_tax$motu_num_f,
#                                                       levels=unique(motus_rel_tax$motu_num_f[order(motus_rel_tax$phylum_f,
#                                                                                                                      motus_rel_tax$class_f,
#                                                                                                                      motus_rel_tax$order_f,
#                                                                                                                      motus_rel_tax$family_f)]),
#                                                       ordered=TRUE)

motus_table_ed_pk_abund_class$class_ed <-  factor( motus_table_ed_pk_abund_class$class_ed,
                                                      levels=unique(motus_rel_tax$class_f[order(motus_rel_tax$class_f)]),
                                                      ordered=TRUE)

motus_table_ed_pk_abund_class$class_ed <-  factor(motus_table_ed_pk_abund_class$class_ed,
                                                      levels=unique(motus_rel_tax$class_f[order(motus_rel_tax$class_f)]),
                                                      ordered=TRUE)

motus_table_ed_pk_abund_class$class_ed |>
  unique()

# motus_table_ed_pk_abund_class$class_ed <- factor(motus_table_ed_pk_abund_class$class_ed,
#                                                      levels = c(NA, 'others',  'Gammaproteobacteria', 'Alphaproteobacteria', 'Cyanobacteriia' ,
#                                                                 'Dehalococcoidia',     'Rhodothermia' ,       'Bacteroidia',         'Actinobacteria',
#                                                                 'Thermoplasmata',  "Marinisomatia"  ,     "Nitrososphaeria"  ,   "Acidimicrobiia"  ,     
#                                                                 "Verrucomicrobiae"  ,  "Kiritimatiellae"   ,    "Poseidoniia" ,  'unassigned'))


motus_table_ed_pk_abund_class$class_ed |>
  unique()

motus_table_ed_pk_abund_motus_plot <- motus_table_ed_pk_abund_class |>
  dplyr::mutate(fraction = '0.2') |>
  group_by(class_ed, sample_id, fraction, date) |>
  dplyr::summarize(sum_rel = sum(relative_abundance)) |>
  dplyr::filter(!class_ed %in% c('unassined')) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, sum_rel))+
  geom_area(aes(fill = class_ed, group = fct_rev(class_ed)), alpha = 1,  position='stack')+
  scale_fill_manual(values = palette_class_assigned_motus)+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 4, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 10),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), strip.placement = 'outside')  

motus_table_ed_pk_abund_motus_plot

# ggsave(filename = 'motus_table_ed_pk_community_class_plot_1perc.pdf', plot = motus_table_ed_pk_abund_motus_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 288, height = 150, units = 'mm')

## general plot -----
motus_rel_tax |>
  ungroup() |>
  distinct(domain)
 
motus_table_ed_pk_domain_plot <- motus_rel_tax  |>
  dplyr::mutate(fraction = '0.2') |>
  group_by(domain, sample_id, fraction, date) |>
  dplyr::summarize(sum_rel = sum(relative_abundance)) |>
  #dplyr::filter(!clas %in% c('unassined')) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, sum_rel))+
  geom_area(aes(fill = domain, group = fct_rev(domain)), alpha = 1,  position='stack')+
  scale_fill_manual(values = palette_long)+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.5)
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 4, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 10),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), strip.placement = 'outside')  

motus_table_ed_pk_domain_plot  

# ggsave(filename = 'motus_table_ed_pk_domain_plot.pdf', plot = motus_table_ed_pk_domain_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 288, height = 150, units = 'mm')

## mOTUs Community Evenness ----
### We need to apply rarefaction because we have uneven sequencing effort, and this will affect community diversity analysis
### Rarefaction (process that repeats the subsampling many times and tries to overcome the uneven sequencing effort bias)----
 #wider format to go into vegan package 

#### calculate the minimum read size of all samples to rarefy at that number
min_n_seqs <- motus_tb  |>
  #rownames_to_column(var = 'motu_num') |>
  pivot_longer(cols = starts_with('BL'), names_to = 'sample_id', values_to = 'reads') |>
  group_by(sample_id) |>
  dplyr::summarize(n_seqs = sum(reads)) |>
  dplyr::summarize(min = min(n_seqs)) |>
  pull(min) #12891

motus_table_ed_pk_ed <- motus_tb |>
  dplyr::select(-tax) |>
  t()

motus_table_tax_col <- motus_tb |>
  dplyr::select(tax) |>
  add_row( .before = 1 ) |>
  dplyr::mutate(tax= ifelse(is.na(tax), 'sample_id', tax))

#this part is commented since is computationally slow so i don't want to repeat this step many times. 
##this function gives us a randomly rarefied community data
# rrarefy(motus_table_ed_pk_ed, sample = min_n_seqs) |>
#   as_tibble(rownames = 'sample_id') #just rarefying (one simgle subsampling)

## we use rarefaction (which repeats the subsampling step many times)
## perform this a 1000 times to get an empirical diversity values calculating the mean value for each timepoint.
# motus_table_ed_pk_rar <- rrarefy.perm(motus_table_ed_pk_ed, # samples in rows!
#                                        sample = min_n_seqs,
#                                        n = 1000,
#                                        round.out = T)

#write.csv2(motus_table_ed_pk_rar, file = 'data/motus_table_ed_pk_rar.csv')

motus_table_ed_pk_rar <- read.csv2('data/motus_table_ed_pk_rar.csv') |>
  as_tibble()|> 
  rename('sample_id' = 'X') 

motus_table_ed_pk_rar |>
  colnames() <- motus_table_tax_col$tax

## Rarefied dataset to calculate Community eveness----
source('../../Bloomers/R/community_evenness.R')

motus_table_ed_pk_rar |>
  dim()

community_eveness_02 <- motus_table_ed_pk_rar |>
  #as.data.frame() |>
  #tibble::rownames_to_column(var = 'sample_id') |>
  dplyr::select(-unassigned) |>
  pivot_longer(cols = !('sample_id'), names_to = 'tax', values_to = 'reads_rar') |>
  #dplyr::select(sample_id, reads, motu_num) |>
  as_tibble() |>
  group_by(sample_id) |>
  dplyr::mutate(reads_rar = as.numeric(reads_rar)) |>
  #ungroup() |>
  dplyr::reframe(community_eveness_rar = community_evenness(abundances = reads_rar, index = 'Pielou')) 

community_eveness_02 <- community_eveness_02 |>
  left_join(m_motus, by = 'sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y")))

community_eveness_02_red_motus <- community_eveness_02 |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d")))

### plot community evenness----
community_eveness_02_red_motus |>
  #left_join(m_motus, by = 'sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, community_eveness_rar))+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E')+
  geom_point()+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Community Eveness')+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank())

## Bray Curtis dissimilarity----
### 0 means the two sites have the same composition (that is they share all the species), and 1 means the two sites 
### do not share any species.
source('../../Bloomers/R/compute_bray_curtis_dissimilariy.R') #we need to upload it since we updated the function but I didn't compile the package
# 
### We need the rarefied table transformed to relative abundances
motus_table_ed_pk_rel_rar <- motus_table_ed_pk_rar |>
  # as.data.frame() |>
  # rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = !('sample_id'), names_to = 'tax', values_to = 'reads') |>
  dplyr::mutate(reads = as.numeric(reads)) |>
  calculate_rel_abund(group_cols = sample_id) 

motus_table_ed_pk_rel_rar <- motus_table_ed_pk_rel_rar %>%
  dplyr::filter(sample_id %in% m_motus_red$sample_id)

motus_table_ed_pk_rel_rar_ed <- motus_table_ed_pk_rel_rar |>
  dplyr::filter(tax != 'unassigned') |>
  dplyr::select(sample_id, asv_num = tax,  reads, total_reads, relative_abundance)

bray_curtis_02_rar_motus <- dissimilarity_matrix(data = motus_table_ed_pk_rel_rar_ed, 
                                           sample_id_col = sample_id,
                                           values_cols_prefix = 'BL')

### plot bray curtis dissimilarity----
bray_curtis_02_rar_plot <- bray_curtis_02_rar_motus |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  left_join(m_motus, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, bray_curtis_result))+
  geom_point()+
  geom_line(aes(date, bray_curtis_result))+
  #facet_grid(vars(fraction))+
  #scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  #geom_smooth(method = 'lm')+
  scale_x_datetime()+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity', color = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom')

bray_curtis_02_rar_plot

### plot Bray-Curtis dissimilarity and Community Evenness together----
bray_curtis_02_rar_motus_red <-  bray_curtis_02_rar_motus |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  left_join(m_motus, by = c('samples' = 'sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d")))

community_diversity_plot <- community_eveness_02_red_motus |> 
  left_join(bray_curtis_02_rar_motus_red) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  #left_join(m_motus, by = c('sample_id')) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, value))+
  #geom_point(aes())+
  geom_line(aes(date, value))+
  geom_smooth(method = 'loess', span = 0.06, color = 'black')+
  geom_smooth(method = 'lm', color = 'darkgreen')+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  facet_grid(vars(diversity_index), labeller = labs_diversity)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Time', y = 'Community diversity', color = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = 'bottom',
        axis.text.x = element_text(size = 6))

community_diversity_plot
# 
# ggsave(filename = 'community_diversity_motus_plot_ed2.pdf', plot = community_diversity_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 180, height = 180, units = 'mm')

## find blooms -----
# Discover anomalies ----
## For each ASVs based on relative abundances (but i also calculate them for rCLR transformation and pseudoabundances to check)-----
#x <- 120*0.75 ##percentage of ASV present in the dataset that we want to subset by (occurrence)

motus_rel_tax 

z_02 <- motus_rel_tax  |>
  group_by(motus_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(relative_abundance >=  0.05)) |> #estan en format 0-1
  #group_by(motu_num) |>
  dplyr::reframe(#anomalies_ps = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1.96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
    anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                 na_rm = TRUE, values = relative_abundance, plotting = FALSE)[c(1,2,3)])

## I only filter for those anomalies in relative abundance because pseudoabundance I can only use it for FL not for PA and rclr will be used to explore if the changes observed are true

# Filter the ASV_tab by only those ASVs that have an anomaly at some point of the dataset (IN THIS CASE WE ARE NOT FILTERING BY ANY RELATIVE ABUNDANCE IN THE COMMUITY JUST IF THEY PRESENT AN ANOMALY AT SOME POINT!)----

asv_anom_02 <- find_asv_with_anomalies(anomalies_result = z_02, 
                                       anomaly_in1 = anomalies_ra, 
                                       anomaly_in2 = NULL,  #anomalies_ps,
                                       anomaly_in3 = NULL, #anomalies_clr, 
                                       logic1 = 'TRUE',
                                       logic2 = NULL, 
                                       logic3 = NULL,
                                       asv_col = motus_num)
## 18 bloomers at 5%.

## Recover z-score for each ASV ----
### I want to highlight anomalies for each ASV to do so I recover z-scores for those ASVs that that have high z-scores
### at some point of the dataset. Easy to observe if those ASVs are having random anomalies or all of them happen at the same time
#x <- 120*0.75
m_motus <- motus_rel_tax |> 
  distinct(date, sample_id, year, month, day) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_motus))) 

m_motus_red <- m_motus |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) 

m_motus_red <- m_motus |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_motus_red)))

z_scores_02 <- motus_rel_tax |>
  dplyr::filter(date %in% m_motus_red$date) |>
  #group_by(tax) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |> ##only anomalies for ASVs that are present in > 50% of the samples
  #dplyr::filter(num_0 < 30) |> ##only anomalies for ASVs that are present in > 25% of the samples
  group_by(motus_num) |>
  dplyr::reframe(z_score_ra = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = relative_abundance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_ra)

duplicated_tax <- z_scores_02 |>
  group_by(motus_num) |>
  reframe(n = n()) |>
  dplyr::filter(n > 72)
  
  motus_rel_tax |>
    dplyr::filter(date %in% m_motus_red$date) |>
    dplyr::filter(motus_num %in% duplicated_tax$motus_num) |>
    dplyr::filter(motus_num == "motus_21") |>
    ggplot(aes(date, relative_abundance))+
    geom_line()+
    geom_point()
  
  # z_scores_02 |>
  #   group_by(motus_num) |>
  # dplyr::mutate(sample_id_num = str_c(1:nrow(m_motus_red))) |>
  #   left_join(m_motus_red, by = 'sample_id_num')

 z_scores_02 <- z_scores_02 |>
    dplyr::filter(!motus_num %in%  duplicated_tax$motus_num) |> # I remove those that are duplicated BUT NEEDS TO BE SOLVED!!!!
    group_by(motus_num) |>
    dplyr::mutate(sample_id_num = str_c(1:nrow(m_motus_red))) |>
    left_join(m_motus_red, by = 'sample_id_num')

## create datasets with the information from 'bloomers' in the mOTUs timeseries ----
z_scores_02_red <- z_scores_02 

motus_table_ed_pk_rel_tax_m_red <- motus_rel_tax |>
  dplyr::filter(date %in%  m_motus_red$date)

motus_table_ed_pk_rel_tax_m_red_z_scores <- motus_table_ed_pk_rel_tax_m_red |>
  #dplyr::select(-sample_id_ok, -sample_id_remainder, -sample_id_width) |>
  left_join(z_scores_02_red)

motus_bloo_tax <- asv_anom_02 |>
  as_tibble_col(column_name = 'motus_num') |>
  left_join(motus_tax_ed_complete_unique) |>
  dplyr::filter(motus_num != 'unassigned') # it's not one taxa

#write.csv(motus_bloo_tax , 'data/motus_bloo_tax.csv')
#write.csv(motus_table_ed_pk_rel_tax_m_red_z_scores, 'data/motus_table_ed_pk_rel_tax_m_red_z_scores.csv')

### Explore what happens if we do not use threshold for abundance at least at some point of the dataset, justify my threshold -----
# Create a function to generate datasets based on the number of ASVs that fulfill the criteria of potential bloomers at different relative abundance thresholds.
source('src/count_number_potential_bloomers_threshold.R')

motus_table_ed_pk_rel_tax_m_red |>
  colnames()

motus_table_ed_pk_rel_tax_m_red_ed <- motus_table_ed_pk_rel_tax_m_red |>
  dplyr::select(everything(), asv_num = motus_num)

z_scores_02_ed <- z_scores_02 |>
  dplyr::select(everything(), asv_num = motus_num)

# I apply a loop for all the thresholds that I'm interested in
# Define a vector of threshold values
threshold_values <- c(0, 0.001, 0.0015, 0.0025, 0.005, 0.0075, 
                      0.01, 0.015, 0.0, 0.025, 0.05, 0.075,
                      0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6)

# Create an empty list to store the results
datasets <- list()

# Iterate over each threshold value and apply the function
for (threshold in threshold_values) {
  dataset_name <- paste0("n_0.2_", threshold * 100)  # Create dataset name
  dataset <- count_num_bloomers(threshold, 0.2, motus_table_ed_pk_rel_tax_m_red_ed, z_scores_tb = z_scores_02_ed)
  datasets[[dataset_name]] <- dataset  # Store the dataset in the list
}

# Combine all datasets into a single dataframe
result_dataset_02 <- bind_rows(datasets)

### plot the results
blooming_threshold <- result_dataset_02 |>
  ggplot(aes(threshold, num))+
  geom_point(size = 1)+
  #scale_x_continuous(expand = c(0,0), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = percent_format())+
  scale_x_log10(expand = c(0,0))+
  #scale_y_log10()+
  facet_grid(vars(fraction), labeller = labs_fraction)+
  geom_vline(xintercept = 0.05, linetype = 'dashed')+
  labs(x = 'Log10 Relative abundance threshold of the potential blooming event', y = 'Number of potential bloomers detected')+
  geom_line()+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 5),
        axis.ticks = element_blank())

blooming_threshold 

# ggsave(blooming_threshold,  filename = 'blooming_threshold_mOTUs_ed1.pdf',
#        path = 'Results/Figures/',
#        width = 88, height = 88, units = 'mm')

## plot the 'bloomers' in the motus time series -----
motus_table_ed_pk_rel_tax_m_red_z_scores_bloo <- motus_table_ed_pk_rel_tax_m_red_z_scores |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  dplyr::filter(motus_num %in% motus_bloo_tax$motus_num)

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo  <- motus_table_ed_pk_rel_tax_m_red_z_scores_bloo  |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                #family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                motus_num_f = as_factor(motus_num))

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f <-  factor(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f, 
                                                                 levels=unique(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f[order(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$phylum_f)]), 
                                                                 ordered=TRUE)

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$order_f <-  factor(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo $order_f, 
                                                                 levels=unique(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo $order_f[order(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$phylum_f,
                                                                                                                                            motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f)]), 
                                                                 ordered=TRUE)

# motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$family_f <- factor(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$family_f, 
#                                                                  levels=unique(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$family_f[order(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$phylum_f,
#                                                                                                                                             motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f,
#                                                                                                                                             motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$order_f)]), 
#                                                                  ordered=TRUE)

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$motus_num_f <-  factor(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$motus_num_f, 
                                                                    levels=unique(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$motus_num_f[order(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$phylum_f,
                                                                                                                                                 motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$class_f,
                                                                                                                                                 motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$order_f)]), 
                                                                    ordered=TRUE)

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot <- motus_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  dplyr::filter(motus_num != 'unasigned') |> # it does not behave as a bloomer
  #dplyr::filter(motus_num != 'Unaligned') |> # it does not behave as a bloomer
  ggplot(aes(date, relative_abundance))+
  #geom_line(aes(group = motus_num))+
  geom_area(aes(fill = order_f))+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  #scale_fill_manual(values = palette_order_assigned_all)+
  facet_wrap(vars(interaction(order_f, motus_num)), ncol = 1, scales = 'free_y')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12), strip.background = element_blank(), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12), strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 2, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot

# ggsave(filename = 'motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_ed2.pdf', plot = motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot,
#        path = 'results/figures/relationship_genes_blooms/',
#        width = 180, height = 250, units = 'mm')

## all together with community evenness ----
motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$order_f |>
  unique()

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot <- motus_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  #dplyr::filter(motu_num != 'ASV1') |> # it does not behave as a bloomer
  dplyr::filter(motus_num != 'unassigned') |> # it does not behave as a bloomer
 # dplyr::filter(motu_num != 'Unaligned') |> # it does not behave as a bloomer
  group_by(date, sample_id) |>
  dplyr::mutate(abund_max = sum(relative_abundance)) |>
  dplyr::group_by(order, date, abund_max) |>
  dplyr::reframe(abund_order = sum(relative_abundance)) |>
  ungroup() |>
  #dplyr::distinct(date, sample_id, order, relative_abundance, abund_max) |>
  #dplyr::left_join(community_eveness_02) |>
  ggplot(aes(date, abund_max))+
  #geom_line(aes( y = abund_max))+
  geom_area(aes(fill = order, y = abund_order), position = 'stack')+
  geom_line(data = community_eveness_02_red_motus, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2004-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot
# # 
# ggsave(filename = 'motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot.pdf', plot = motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot,
#        path = 'results/figures/',
#        width = 300, height = 180, units = 'mm')

## all together with community evenness ----
motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot <- motus_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  #dplyr::filter(motu_num != 'ASV1') |> # it does not behave as a bloomer
  dplyr::filter(motus_num != 'unassigned') |> # it does not behave as a bloomer
  #dplyr::filter(motu_num != 'Unaligned') |> # it does not behave as a bloomer
  group_by(date, sample_id) |>
  dplyr::mutate(abund_max = sum(relative_abundance)) |>
  dplyr::group_by(order, date, abund_max) |>
  dplyr::reframe(abund_order = sum(relative_abundance)) |>
  ungroup() |>
  #dplyr::distinct(date, sample_id, order, relative_abundance, abund_max) |>
  #dplyr::left_join(community_eveness_02) |>
  ggplot(aes(date, abund_max))+
  #geom_line(aes( y = abund_max))+
  geom_area(aes(fill = order, y = abund_order), position = 'stack')+
  geom_line(data = community_eveness_02_red, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2008-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot
# # 
# ggsave(filename = 'motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot_ed2.pdf', plot = motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot,
#        path = 'results/figures/',
#        width = 300, height = 180, units = 'mm')


# -------------------- ################## COMPOSITION PLOTs ##################  ----------

## arrange the community of the three different approximations class level----

### mTAGs ----
mtags_table_ed_pk_community_class_plot <- mtags_table_ed_pk_community_class  |>
  dplyr::mutate(fraction = '0.2') |>
  group_by(class_ed, sample_id, fraction, date) |>
  dplyr::summarize(sum_rel = sum(relative_abundance)) |>
  dplyr::filter(!class_ed %in% c('Unassigned', 'Unaligned')) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, sum_rel))+
  geom_area(aes(fill = class_ed, group = fct_rev(class_ed)), alpha = 1,  position='stack')+
  scale_fill_manual(values = palette_class_assigned_mtags_motus)+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class', title = 'mTAGs')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 10),
        legend.position = 'none', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), strip.placement = 'outside')  

mtags_table_ed_pk_community_class_plot

###  ASVs -----
asv_tab_10y_02_rel_tax_m_red <- asv_tab_10y_02_rel  |>
  left_join(tax_bbmo_10y_new) |>
  left_join(m_02) |>
  #dplyr::filter(!is.na(phylum))   |>
  #group_by(class, sample_id, date) |>
  #dplyr::reframe(sum_rel = sum(relative_abundance)) |>
  #left_join(m_mtags, by = c('sample_id') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))  |>
  dplyr::filter(date > (as.POSIXct('2008-01-01', format = "%Y-%m-%d")))

###  colored only those taxa that represent > 1% at some point 
asv_tab_10y_02_rel_tax_m_red_community_class <- asv_tab_10y_02_rel_tax_m_red |>
  #dplyr::filter(!is.na(class))   |>
  dplyr::mutate(class_ed = case_when(
    class %in% unique((asv_tab_10y_02_rel_tax_m_red  |> 
                         group_by(asv_num) |> 
                         filter(any(relative_abundance > 0.01)) |> 
                         pull(class))) ~ class,
    TRUE ~ 'others'))

asv_tab_10y_02_rel_tax_m_red_community_class$class_ed <- factor(asv_tab_10y_02_rel_tax_m_red_community_class$class_ed, 
                                                      levels=unique(asv_tab_10y_02_rel_tax_m_red_community_class$class_ed[order(asv_tab_10y_02_rel_tax_m_red_community_class$class)]), 
                                                      ordered=TRUE)

asv_tab_10y_02_rel_tax_m_red_community_class$class_ed <-  factor(asv_tab_10y_02_rel_tax_m_red_community_class$class_ed, 
                                                      levels=unique(asv_tab_10y_02_rel_tax_m_red_community_class$class_ed[order(asv_tab_10y_02_rel_tax_m_red_community_class$class)]), 
                                                      ordered=TRUE)

asv_tab_10y_02_rel_tax_m_red_community_class$class_ed |>
  unique()

asv_tab_10y_02_rel_tax_m_red_community_class$class_ed <- factor(asv_tab_10y_02_rel_tax_m_red_community_class$class_ed, 
                                                     levels = c(NA, 'others', 'Verrucomicrobiae', 'Gammaproteobacteria', 'Alphaproteobacteria', 'Cyanobacteriia' ,    
                                                                'Dehalococcoidia',    "Planctomycetes" ,  'Rhodothermia' ,       'Bacteroidia',  'Actinobacteria', "Desulfuromonadia",
                                                                'Thermoplasmata',  "Campylobacteria" ,   "Clostridia", "Bacilli"  ,"Acidimicrobiia"  ))

asv_tab_10y_02_rel_tax_m_red_community_class$class_ed |>
  unique()

asv_tab_10y_02_rel_tax_m_red_community_class_plot <- asv_tab_10y_02_rel_tax_m_red_community_class  |>
  dplyr::mutate(fraction = '0.2') |>
  group_by(class_ed, sample_id, fraction, date) |>
  dplyr::summarize(sum_rel = sum(relative_abundance)) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, sum_rel))+
  geom_area(aes(fill = class_ed, group = fct_rev(class_ed)), alpha = 1,  position='stack')+
  scale_fill_manual(values = palette_class_assigned_mtags_motus)+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class', title = 'ASVs')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 10),
        legend.position = 'none', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), strip.placement = 'outside')  

asv_tab_10y_02_rel_tax_m_red_community_class_plot

### mOTUs ----
### colored only those taxa that represent > 1% at some point 
motus_table_ed_pk_abund_class <- motus_rel_tax |>
  #dplyr::filter(!is.na(class))   |>
  dplyr::mutate(class_ed = case_when(
    class %in% unique((motus_rel_tax |> 
                         group_by(motus_num) |> 
                         filter(any(relative_abundance > 0.01)) |> 
                         pull(class))) ~ class,
    TRUE ~ 'others'))

motus_rel_tax |>
  colnames()

motus_rel_tax <- motus_rel_tax  |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                #family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class))

motus_rel_tax$class_f <-  factor(motus_rel_tax$class_f,
                                 levels=unique(motus_rel_tax$class_f[order(motus_rel_tax$phylum_f)]),
                                 ordered=TRUE)

motus_rel_tax$order_f <-  factor(motus_rel_tax$order_f,
                                 levels=unique(motus_rel_tax$order_f[order(motus_rel_tax$phylum_f,
                                                                           motus_rel_tax$class_f)]),
                                 ordered=TRUE)

# motus_rel_tax$family_f <-  factor(motus_rel_tax$family_f,
#                                                     levels=unique(motus_rel_tax$family_f[order(motus_rel_tax$phylum_f,
#                                                                                                                  motus_rel_tax$class_f,
#                                                                                                                  motus_rel_tax$order_f)]),
#                                                     ordered=TRUE)

# motus_rel_tax$motu_num_f <-  factor(motus_rel_tax$motu_num_f,
#                                                       levels=unique(motus_rel_tax$motu_num_f[order(motus_rel_tax$phylum_f,
#                                                                                                                      motus_rel_tax$class_f,
#                                                                                                                      motus_rel_tax$order_f,
#                                                                                                                      motus_rel_tax$family_f)]),
#                                                       ordered=TRUE)

motus_table_ed_pk_abund_class$class_ed <-  factor( motus_table_ed_pk_abund_class$class_ed,
                                                   levels=unique(motus_rel_tax$class_f[order(motus_rel_tax$class_f)]),
                                                   ordered=TRUE)

motus_table_ed_pk_abund_class$class_ed <-  factor(motus_table_ed_pk_abund_class$class_ed,
                                                  levels=unique(motus_rel_tax$class_f[order(motus_rel_tax$class_f)]),
                                                  ordered=TRUE)

motus_table_ed_pk_abund_class$class_ed |>
  unique()

# motus_table_ed_pk_abund_class$class_ed <- factor(motus_table_ed_pk_abund_class$class_ed,
#                                                      levels = c(NA, 'others',  'Gammaproteobacteria', 'Alphaproteobacteria', 'Cyanobacteriia' ,
#                                                                 'Dehalococcoidia',     'Rhodothermia' ,       'Bacteroidia',         'Actinobacteria',
#                                                                 'Thermoplasmata',  "Marinisomatia"  ,     "Nitrososphaeria"  ,   "Acidimicrobiia"  ,     
#                                                                 "Verrucomicrobiae"  ,  "Kiritimatiellae"   ,    "Poseidoniia" ,  'unassigned'))

motus_table_ed_pk_abund_class$class_ed |>
  unique()

motus_table_ed_pk_abund_motus_plot <- motus_table_ed_pk_abund_class |>
  dplyr::filter(date < (as.POSIXct('2014-01-01', format = "%Y-%m-%d"))) |>
  dplyr::mutate(fraction = '0.2') |>
  group_by(class_ed, sample_id, fraction, date) |>
  dplyr::summarize(sum_rel = sum(relative_abundance)) |>
  dplyr::filter(!class_ed %in% c('unassigned')) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y"))) |>
  ggplot(aes(date, sum_rel))+
  geom_area(aes(fill = class_ed, group = fct_rev(class_ed)), alpha = 1,  position='stack')+
  scale_fill_manual(values = palette_class_assigned_motus)+
  scale_x_datetime(expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Class', title = 'mOTUs')+
  facet_wrap(vars(fraction), dir = 'v', scales = 'free',  labeller = labs_fraction)+
  guides(fill = guide_legend(ncol = 4, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 10),
        legend.position = 'none', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), strip.placement = 'outside')  

motus_table_ed_pk_abund_motus_plot

### create a plot to extract the legend ----
motus_table_ed_pk_abund_class_red <- motus_table_ed_pk_abund_class |>
  dplyr::select(class_ed, relative_abundance, date)
  
asv_tab_10y_02_rel_tax_m_red_community_class_red <- asv_tab_10y_02_rel_tax_m_red_community_class |>
  ungroup() |> 
  dplyr::select(class_ed, relative_abundance, date)
  
mtags_table_ed_pk_community_class_red <-  mtags_table_ed_pk_community_class |>
  ungroup() |> 
  dplyr::select(class_ed, relative_abundance, date)

dataset_all <- motus_table_ed_pk_abund_class_red |>
  dplyr::mutate(class_ed = as.character(class_ed)) |>
  bind_rows(asv_tab_10y_02_rel_tax_m_red_community_class_red |>
              dplyr::mutate(class_ed = as.character(class_ed))) |>
  bind_rows(mtags_table_ed_pk_community_class_red |>
              dplyr::mutate(class_ed = as.character(class_ed)))

dataset_all |>
  distinct(class_ed) |>
  as_vector()

dataset_all$class_ed <- factor(dataset_all$class_ed, 
                               levels = c('Gammaproteobacteria', 'Alphaproteobacteria', 'Cyanobacteriia' , "Bacteroidia", 
                                          'Dehalococcoidia',  'Rhodothermia' ,  'Actinobacteria',
                                          "Kiritimatiellae", "Desulfuromonadia", "Campylobacteria" , "Clostridia",
                                          'Thermoplasmata',  "Marinisomatia"  , "Nitrososphaeria"  ,   "Acidimicrobiia"  ,     
                                          "Verrucomicrobiae"  ,     "Poseidoniia" , "Planctomycetia", "Planctomycetes",
                                          "UBA1144", "Bacilli" , 'others',
                                          'unassigned',  "Unassigned"))

all_plot_for_legend <- dataset_all |>
  ggplot(aes(date, relative_abundance))+
           geom_area(aes(fill = class_ed))+
  scale_fill_manual(values = palette_class_assigned_mtags_motus)+
  labs(fill = 'Class')+
  theme(legend.position = 'bottom', legend.text = element_text(size = 8))+
  guides(guide_legend(fill = guides(ncol = 4)))
  
legend_class_asv_motus_mtags <- get_legend(all_plot_for_legend)

### all together ----
motus_mtags_asvs_community_class_1perc <- grid.arrange(asv_tab_10y_02_rel_tax_m_red_community_class_plot,
             mtags_table_ed_pk_community_class_plot,
             motus_table_ed_pk_abund_motus_plot,
             legend_class_asv_motus_mtags,
             ncol = 1,
             widths = c(1), heights = c(1, 1, 1, 0.5)
)

# ggsave(filename = 'motus_mtags_asvs_community_class_1perc_plot.pdf', plot = motus_mtags_asvs_community_class_1perc,
#        path = 'results/figures/',
#        width = 180, height = 350, units = 'mm')

## Arrange the potential blooms detected in the three different approximations ----
mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot <- mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  #dplyr::filter(motu_num != 'ASV1') |> # it does not behave as a bloomer
  dplyr::filter(motu_num != 'Unassigned') |> # it does not behave as a bloomer
  dplyr::filter(motu_num != 'Unaligned') |> # it does not behave as a bloomer
  ungroup() |>
  dplyr::mutate(order = case_when(order == 'Alteromonadales' ~ 'Enterobacterales',
                                   order != 'Alteromonadales'  ~ order)) |>
  group_by(date, sample_id) |>
  dplyr::mutate(abund_max = sum(relative_abundance)) |>
  dplyr::group_by(order, date, abund_max) |>
  dplyr::reframe(abund_order = sum(relative_abundance)) |>
  ungroup() |>
  #dplyr::distinct(date, sample_id, order, relative_abundance, abund_max) |>
  #dplyr::left_join(community_eveness_02) |>
  ggplot(aes(date, abund_max))+
  #geom_line(aes( y = abund_max))+
  geom_area(aes(fill = order, y = abund_order), position = 'stack')+
  geom_line(data = community_eveness_02_red, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2008-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'none', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot

asv_table_02_bloo_plot_evenness_plot <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(!asv_num  %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  group_by(date, sample_id) |>
  dplyr::mutate(abund_max = sum(abundance_value)) |>
  dplyr::group_by(order, date, abund_max) |>
  dplyr::reframe(abund_order = sum(abundance_value)) |>
  ungroup() |>
  #dplyr::distinct(date, sample_id, order, relative_abundance, abund_max) |>
  #dplyr::left_join(community_eveness_02) |>
  ggplot(aes(date, abund_max))+
  #geom_line(aes( y = abund_max))+
  geom_area(aes(fill = order, y = abund_order), position = 'stack')+
  geom_line(data = community_eveness_02_red, aes(date, community_eveness_rar))+
  scale_x_datetime(limits = c(as.POSIXct('2008-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*1 , name = 'Community Evenness'))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'none', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot <- motus_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  #dplyr::filter(motu_num != 'ASV1') |> # it does not behave as a bloomer
  dplyr::filter(motus_num != 'unassigned') |> # it does not behave as a bloomer
  # dplyr::filter(motu_num != 'Unaligned') |> # it does not behave as a bloomer
  group_by(date, sample_id) |>
  dplyr::mutate(abund_max = sum(relative_abundance)) |>
  dplyr::group_by(order, date, abund_max) |>
  dplyr::reframe(abund_order = sum(relative_abundance)) |>
  ungroup() |>
  #dplyr::distinct(date, sample_id, order, relative_abundance, abund_max) |>
  #dplyr::left_join(community_eveness_02) |>
  ggplot(aes(date, abund_max))+
  #geom_line(aes( y = abund_max))+
  geom_area(aes(fill = order, y = abund_order), position = 'stack')+
  geom_line(data = community_eveness_02_red_motus, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2008-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'none', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot

## general legend for all the plots 
motus_table_ed_pk_abund_class_red <- motus_table_ed_pk_rel_tax_m_red_z_scores_bloo|>
  dplyr::select(order_f, relative_abundance, date) |>
  dplyr::filter(order_f != 'unassigned')

asv_tab_all_bloo_z_tax_red <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  ungroup() |> 
  dplyr::select(order_f = order, relative_abundance = abundance_value, date)

mtags_table_ed_pk_community_class_red <-  mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo |>
  ungroup() |> 
  dplyr::select(order_f, relative_abundance, date) |>
  dplyr::filter(!order_f %in% c( "Unassigned", 'Unaligned'))

dataset_all <- motus_table_ed_pk_abund_class_red |>
  dplyr::mutate(order_f = as.character(order_f)) |>
  bind_rows(asv_tab_all_bloo_z_tax_red |>
              dplyr::mutate(order_f = as.character(order_f))) |>
  bind_rows(mtags_table_ed_pk_community_class_red |>
              dplyr::mutate(order_f = as.character(order_f)))
# 
# levels(motus_table_ed_pk_rel_tax_m_red_z_scores_bloo$order_f)
# levels(mtags_table_ed_pk_community_class_red$order_f)
# unique(asv_tab_all_bloo_z_tax_red$order)

dataset_all |>
  distinct(order_f) |>
  as_vector()

dataset_all$order_f <- factor(dataset_all$order_f, 
                               levels = c("SAR11 clade",  "Sphingomonadales", "Rhizobiales"  ,
                                          'Rhodobacterales', "Rhodospirillales"  , 'Pelagibacterales', 'Puniceispirillales',
                                           "SAR86" ,"Cellvibrionales", "Pseudomonadales", "Alteromonadales", "Enterobacterales",
                                          'Nitrososphaerales', 'PCC-6307', 
                                          'Poseidoniales', "Thiotrichales"  , 'PS1', "Synechococcales" ,  "Flavobacteriales" ,  "Bacteriovoracales",
                                          "Verrucomicrobiales"
                                          
                                           ))

all_plot_for_legend <- dataset_all |>
  ggplot(aes(date, relative_abundance))+
  geom_area(aes(fill = order_f))+
  scale_fill_manual(values = palette_order_assigned_all)+
  labs(fill = 'Order')+
  theme(legend.position = 'bottom', legend.text = element_text(size = 8))+
  guides(guide_legend(fill = guides(ncol = 4)))

legend_order <- get_legend(all_plot_for_legend)

## all together
bloomers_asvs_mtags_motus <- grid.arrange(asv_table_02_bloo_plot_evenness_plot,
             mtags_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot,
             motus_table_ed_pk_rel_tax_m_red_z_scores_bloo_plot_evenness_plot,
             legend_order,
            ncol = 1,
             widths = c(1), heights = c(1, 1, 1,0.5)
            )
# ggsave(filename = 'bloomers_asvs_mtags_motus_plot.pdf', plot = bloomers_asvs_mtags_motus,
#        path = 'results/figures/',
#        width = 180, height = 350, units = 'mm')

## Bray Curtis with the three different approximations ----
bray_curtis_02_motus <- bray_curtis_02_rar_motus_red |>
  dplyr::select(sample_id_ed = samples, bray_curtis_result_motus = bray_curtis_result, date)

bray_curtis_02_mtags <- bray_curtis_02_rar_mtags_red |>
  dplyr::select(sample_id_ed = samples, bray_curtis_result_mtags = bray_curtis_result, date) 

bray_curtis_02_asvs <- bray_curtis_rar_all_m |>  
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(sample_id_ed = str_replace(samples, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::select(sample_id_ed, bray_curtis_result_asv = bray_curtis_result, date) |>
  dplyr::mutate(sample_id_ed = case_when(sample_id_ed == 'BL110705' ~ 'BL110704',
                                         sample_id_ed ==  'BL100414' ~ 'BL100413',
                                         sample_id_ed ==  'BL120511' ~ 'BL120518',
                                         sample_id_ed ==  'BL131215' ~ 'BL131204',
                                         sample_id_ed == 'BL080312' ~ 'BL080311',
                                         T ~ sample_id_ed))

comparison_bray_curtis_asvs_motus_tb <- bray_curtis_02_mtags |>
  left_join(bray_curtis_02_asvs, by = c('sample_id_ed')) |>
  left_join(bray_curtis_02_motus, by = c('sample_id_ed'))

comparison_bray_curtis_asvs_motus_tb_l <- comparison_bray_curtis_asvs_motus_tb |>
  pivot_longer(cols = starts_with('bray')) |>
  dplyr::mutate(bray_curtis = case_when(str_detect(name, 'mtags') ~ 'mTAGs',
                                               str_detect(name, 'asv') ~ 'ASVs',
                                               str_detect(name, 'motus') ~ 'mOTUs')) |>
  dplyr::mutate(date.x = (as.POSIXct(date.x, format = "%d-%m-%y"))) 

comparison_bray_curtis_asvs_motus_tb_l$bray_curtis <- factor(comparison_bray_curtis_asvs_motus_tb_l$bray_curtis,
                                                                           levels = c('mTAGs', 'ASVs', 'mOTUs'))

comparison_bray_curtis_mTAGs_ASVs_plot <-  comparison_bray_curtis_asvs_motus_tb_l |>
  ggplot(aes(date.x, value, color = bray_curtis, fill = bray_curtis))+
  facet_wrap(vars(bray_curtis), ncol = 1)+
  scale_color_manual(values = palette_motus_asvs_mtags)+
  scale_fill_manual(values = palette_motus_asvs_mtags)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(aes(group = bray_curtis), linewidth = 1)+
  labs(x = 'Time', y = 'Bray Curtis Dissimilarity', fill = 'Community data', color = 'Community data' )+
  #geom_smooth(method = 'loess', span = 0.07)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.position = 'bottom')

comparison_bray_curtis_mTAGs_ASVs_plot              
# 
# ggsave(comparison_bray_curtis_mTAGs_ASVs_plot,  filename = 'comparison_bray_curtis_mTAGs_ASVs_mOTUs_plot.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 128, units = 'mm')

### I plot mTAGs, mOTUS and ASVs commnity Evenness together -----
## datasets 
community_eveness_02_red_motus
community_eveness_02_mtags
community_evenness_02_asvs <- community_eveness_all_m |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::select(sample_id_ed, community_eveness_rar_asv = community_eveness_rar, date)

community_eveness_02_red_mtags <- community_eveness_02_mtags |>
  dplyr::select(-station, -year, -month, -day, -sample_id_ok, -sample_id_remainder, -sample_id_width, community_eveness_rar_mtags = community_eveness_rar)

## edit the dates since the are some dates that are different fom mOTUs and ASVs (one day difference)
community_evenness_02_asvs <- community_evenness_02_asvs |>
  dplyr::mutate(sample_id_ed = case_when(sample_id_ed == 'BL110705' ~ 'BL110704',
                                         sample_id_ed ==  'BL100414' ~ 'BL100413',
                                         sample_id_ed ==  'BL120511' ~ 'BL120518',
                                         sample_id_ed ==  'BL131215' ~ 'BL131204',
                                         sample_id_ed == 'BL080312' ~ 'BL080311',
                                         T ~ sample_id_ed))

comparison_community_evenness_asvs_motus_tb <- community_eveness_02_red_motus |>
  dplyr::select(-station, -year, -month, -day, -sample_id_ok, -sample_id_remainder, -sample_id_width, community_eveness_rar_motus = community_eveness_rar) |>
  left_join(community_eveness_all_m_red, by = c('sample_id' = 'sample_id_ed')) |>
  left_join(community_eveness_02_red_mtags)

comparison_community_evenness_asvs_motus_tb |>
  dplyr::filter(is.na(community_eveness_rar_asv))

# sample_id community_eveness_rar date                season sample_id_ed community_eveness_rar_asv
# <chr>                     <dbl> <dttm>              <chr>  <chr>                            <dbl>
#   1 BL100413              0.8258590 2010-04-13 00:00:00 spring NA                                  NA
# 2 BL110704              0.7983031 2011-07-04 00:00:00 summer NA                                  NA
# 3 BL120518              0.6883892 2012-05-18 00:00:00 spring NA                                  NA
# 4 BL131204              0.8384126 2013-12-04 00:00:00 autumn NA                                  NA

# community_eveness_all_m_red |>
#   dplyr::filter(str_detect(date, '2008'))

correlation_community_evenness_ASvs_mOTUs_plot <- comparison_community_evenness_asvs_motus_tb |>
  ggplot(aes(community_eveness_rar_asv, community_eveness_rar_mtags))+
  geom_point(aes(), color = '#344848')+
  labs(x = 'ASVs Community Evenness', y = 'mTAGs Community Evenness')+
  geom_smooth(method = 'lm', color = '#344848')+
  geom_text(aes(label = `sample_id`), check_overlap = F, nudge_x = 0.02, nudge_y = 0.04)+
  theme_bw()+
  theme(strip.background = element_blank(),
        aspect.ratio = 4/4,
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.position = 'bottom')

correlation_community_evenness_ASvs_mOTUs_plot  

# ggsave(correlation_community_evenness_ASvs_mOTUs_plot,  filename = 'correlation_community_evenness_ASvs_mOTUs_plot_labs_outliers.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 188, units = 'mm')

comparison_community_evenness_asvs_motus_tb_l <- comparison_community_evenness_asvs_motus_tb |>
  pivot_longer(cols = starts_with('community')) |>
  dplyr::mutate(community_evenness = case_when(str_detect(name, 'mtags') ~ 'mTAGs',
                                               str_detect(name, 'asv') ~ 'ASVs',
                                               str_detect(name, 'motus') ~ 'mOTUs')) |>
  dplyr::mutate(date.x = (as.POSIXct(date.x, format = "%d-%m-%y"))) 

comparison_community_evenness_asvs_motus_tb_l$community_evenness <- factor(comparison_community_evenness_asvs_motus_tb_l$community_evenness,
                                                                           levels = c('mTAGs', 'ASVs', 'mOTUs'))

comparison_community_evenness_mTAGs_ASVs_plot <-  comparison_community_evenness_asvs_motus_tb_l |>
  ggplot(aes(date.x, value, color = community_evenness, fill = community_evenness))+
  facet_wrap(vars(community_evenness), ncol = 1)+
  scale_color_manual(values = palette_motus_asvs_mtags)+
  scale_fill_manual(values = palette_motus_asvs_mtags)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_line(aes(group = community_evenness), linewidth = 1)+
  labs(x = 'Time', y = 'Community Evenness', fill = 'Community data', color = 'Community data' )+
  #geom_smooth(method = 'loess', span = 0.07)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.position = 'bottom')

comparison_community_evenness_mTAGs_ASVs_plot              
# 
# ggsave(comparison_community_evenness_mTAGs_ASVs_plot,  filename = 'comparison_community_evenness_mTAGs_ASVs_mOTUs_plot_ed2.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 128, units = 'mm')


## Cyanobacteria differences between the different approximations -----
### Comparison abundances mTAGs, ASVs and flow cytometry of Cyanobacteria ----

cianos_mtags <- mtags_table_ed_pk_rel_tax_m |>
  dplyr::filter(phylum == 'Cyanobacteria') |>
  ungroup() |>
  dplyr::mutate(sample_id_ed = case_when(sample_id == 'BL110705' ~ 'BL110704',
                                         sample_id ==  'BL100414' ~ 'BL100413',
                                         sample_id ==  'BL120511' ~ 'BL120518',
                                         sample_id ==  'BL131215' ~ 'BL131204',
                                         sample_id == 'BL080312' ~ 'BL080311',
                                         T ~ sample_id)) |>
  #dplyr::distinct(order, family, class)
  dplyr::select(date, sample_id_ed, relative_abundance) |>
  group_by(date, sample_id_ed) |>
  dplyr::reframe(cianos_total = sum(relative_abundance)) |>
  dplyr::mutate(type = 'mtags') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y")))

cianos_motus <- motus_rel_tax |>
  dplyr::filter(phylum == 'Cyanobacteria') |>
  ungroup() |>
  dplyr::mutate(sample_id_ed = case_when(sample_id == 'BL110705' ~ 'BL110704',
                                         sample_id ==  'BL100414' ~ 'BL100413',
                                         sample_id ==  'BL120511' ~ 'BL120518',
                                         sample_id ==  'BL131215' ~ 'BL131204',
                                         sample_id == 'BL080312' ~ 'BL080311',
                                         T ~ sample_id)) |>
  #dplyr::distinct(order, family, class)
  dplyr::select(date, sample_id_ed, relative_abundance) |>
  group_by(date, sample_id_ed) |>
  dplyr::reframe(cianos_total = sum(relative_abundance)) |>
  dplyr::mutate(type = 'motus') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%d-%m-%y")))

asv_tab_10y_02_rel |> 
  left_join(tax_bbmo_10y_new) |>
  ungroup() |>
  dplyr::distinct(phylum) |>
  as_vector()

cianos_asvs <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new) |>
  ungroup() |>
  dplyr::filter(phylum == 'Cyanobacteria') |>
  #dplyr::distinct(order, family, class)
  left_join(m_bbmo_10y) |>
  dplyr::select(date, relative_abundance, sample_id) |>
  group_by(date, sample_id) |>
  dplyr::reframe(cianos_total = sum(relative_abundance)) |>
  dplyr::mutate(type = 'asvs') |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(-sample_id)

m_bbmo_10y |>
  colnames()

cianos_fcm <- m_bbmo_10y |>
  dplyr::filter(fraction != '3') |>
  dplyr::select(date, prochlorococcus_FC, synechococcus,  sample_id) |>
  pivot_longer(cols = !c('date', 'sample_id')) |>
  group_by(date, sample_id) |>
  dplyr::reframe(cianos_total = sum(value)) |>
  dplyr::select(date, cianos_total, sample_id) |>
  dplyr::mutate(type = 'fcm') |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(-sample_id)

total_bac_fcm <-  m_bbmo_10y |>
  dplyr::filter(fraction != '3') |>
  dplyr::select(date, bacteria_joint, sample_id) |>
  pivot_longer(cols = !c('date', 'sample_id')) |>
  group_by(date, sample_id) |>
  dplyr::reframe(bac_total = sum(value)) |>
  dplyr::select(date, bac_total, sample_id) |>
  dplyr::mutate(type = 'fcm') |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(-sample_id)

cianos_fcm <- cianos_fcm |>
  left_join(total_bac_fcm) |>
  dplyr::mutate(cianos_total = cianos_total/bac_total) |>
  dplyr::select(-bac_total)

## change citometry data by relative citometry abundance (Pep counts) ----
citometry_data_pep <- read_xlsx('data/env_data/Dades_cianos _per _Ona.xlsx', skip = 1)

citometry_data_pep |>
  colnames() <- c('bacteria', 'synecho', 'prochloro', 'per_year', 'decimal_date', 'ratio_ciano_all')

citometry_data_pep |>
  dim()

citometry_data_pep <- citometry_data_pep |>
  dplyr::mutate(decimal_date_round = round(decimal_date, digits = 2)) |>
  right_join(m_02, by = c('decimal_date_round' = 'decimal_date')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

citometry_data_pep |>
  ggplot(aes(date, ratio_ciano_all))+
  geom_line()+
  geom_smooth(method = 'loess', span = 0.088, color = 'black')+
  scale_x_datetime(limits = c(as.POSIXct('2004-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0),
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  # scale_x_continuous(limits = c(2004, 2018))+
  labs(y = 'Ratio Cyanobacteria / Prokaryotes')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

citometry_data_pep |>
  dplyr::select(date, ratio_ciano_all) |>
  left_join(cianos_fcm) |>
  ggplot(aes(ratio_ciano_all, cianos_total))+
  geom_point()

cianos_fcm <- citometry_data_pep |>
  dplyr::mutate(type = 'fcm') |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::select(date, cianos_total = ratio_ciano_all, type, sample_id_ed )

## create a dataset with all the cianobacteria information -----
cianos_data <- cianos_asvs |>
  #bind_rows(cianos_motus) |>
  bind_rows(cianos_fcm) |>
  #bind_rows(cianos_mtags) |>
  dplyr::select(-date) |>
  distinct(cianos_total, type, sample_id_ed)

##before correlation we check normality
shapiro.test(as.numeric(cianos_mtags$cianos_total)) # => p-value = 3.391e-11 (NO NORMALITY)
ggqqplot(as.numeric(cianos_mtags$cianos_total))

shapiro.test(as.numeric(cianos_asvs$cianos_total)) # =>p-value = 1.028e-06 (NO NORMALITY)
ggqqplot(as.numeric(cianos_asvs$cianos_total))

shapiro.test(as.numeric(cianos_motus$cianos_total)) # => p-value = 1e-12 (NO NORMALITY)
ggqqplot(as.numeric(cianos_motus$cianos_total))

shapiro.test(as.numeric(cianos_fcm$cianos_total)) # => p-value = 1.373e-08 (NO NORMALITY)
ggqqplot(as.numeric(cianos_fcm$cianos_total))

# Calculate correlations for each type
cianos_data_w <- cianos_data |>
  ungroup() |> 
  distinct(sample_id_ed, type,.keep_all = TRUE) |>
  ungroup() |>
  pivot_wider(id_cols = c( 'sample_id_ed'), values_from = 'cianos_total', names_from = 'type', values_fill = NA) |>
  as_tibble() |>
  dplyr::mutate(asvs = as.numeric(asvs), 
               # motus = as.numeric(motus),
                fcm = as.numeric(fcm) #, mtags = as.numeric(mtags)
               )

cianos_asvs_mtags_plot <- cianos_data_w |>
  ggplot(aes(asvs, mtags))+
  geom_point()+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.15,
    #label.y = -0.02,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  labs(x = 'ASVs Cyanobacteria', y = 'mTAGs Cyanobacteria')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside', aspect.ratio = 4/4)

cianos_asvs_mtags_plot

cianos_asvs_fcm_plot <- cianos_data_w |>
  ggplot(aes(asvs, fcm))+
  geom_point()+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.25,
    #label.y = -0.02,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  labs(x = 'ASVs Cyanobacteria', y = 'Flow cytometry Cyanobacteria')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside', aspect.ratio = 4/4,
        panel.border = element_blank())

cianos_asvs_fcm_plot <- cianos_data_w |>
  ggplot(aes(asvs, fcm))+
  geom_abline(slope = 1, intercept = 0, color = '#0F3331', linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.25,
    #label.y = -0.02,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.25, label.y = 0.1, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  labs(x = 'Ratio ASVs Cyanobacteria / Total Bacteria ', y = 'Ratio Flow cytometry Cyanobacteria / Total Bacteria')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 10), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

cianos_asvs_fcm_plot

# ggsave(cianos_asvs_fcm_plot,  filename = 'cianos_asvs_fcm_plot_ed4.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 188, units = 'mm')

cianos_asvs_motus_plot <- cianos_data_w |>
  ggplot(aes(asvs, motus))+
  geom_point()+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.25,
    #label.y = -0.02,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  labs(x = 'ASVs Cyanobacteria', y = 'mOTUs Cyanobacteria')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside', aspect.ratio = 4/4)

cianos_mtags_fcm_plot <- cianos_data_w |>
  ggplot(aes(mtags, fcm))+
  geom_point()+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.0,
    label.y = 0.1,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  labs(x = 'mTAGs Cyanobacteria relative', y = 'Flow cytometry Cyanobacteria')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside', aspect.ratio = 4/4)

cianos_motus_fcm_plot <- cianos_data_w |>
  ggplot(aes(motus, fcm))+
  geom_point()+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.0,
    label.y = 0.1,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  labs(x = 'mOTUs Cyanobacteria relative', y = 'Flow cytometry Cyanobacteria relative')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside', aspect.ratio = 4/4)

cianos_motus_mtags_plot <- cianos_data_w |>
  ggplot(aes(mtags, motus))+
  geom_point()+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.0,
    label.y = 0.1,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  labs(x = 'mTAGs Cyanobacteria', y = 'mOTUs Cyanobacteria')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside', aspect.ratio = 4/4)

correlation_cianobacteria <- grid.arrange(cianos_asvs_mtags_plot,
                                          cianos_asvs_fcm_plot,
                                          cianos_asvs_motus_plot,
                                          cianos_mtags_fcm_plot,
                                          cianos_motus_fcm_plot,
                                          cianos_motus_mtags_plot,
                                          ncol = 3)

correlation_cianobacteria

# ggsave(correlation_cianobacteria,  filename = 'correlation_cianobacteria_complete_plot.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 138, units = 'mm')

cianos_data$type <- factor(cianos_data$type, levels = c('fcm', 'asvs', 'mtags', 'motus'), labels = c('Flow cytometry', 'ASVs', 'mTAGs','mOTUs'))

bbmo_cianos_plot <- cianos_data |>
  ggplot(aes(date, cianos_total))+
  facet_wrap(vars(type), scales = 'free_y', ncol =1)+
  geom_line()+
  # geom_line(data = community_eveness_02_red, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2004-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  # scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Time', y = 'Cyanobacteria')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside') +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')
# 
# ggsave(bbmo_cianos_plot,  filename = 'bbmo_cianos_plot_ed2.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 188, units = 'mm')


## Only the dynamics over the timeseries ---
cianos_data %$%
  unique(type)

cianos_asv_fcm_time_plot <- cianos_data |>
  dplyr::filter(type %in% c('ASVs', 'Flow cytometry')) |>
  tidyr::separate_wider_position(col = sample_id_ed, widths = c('station' = 2, 'year' = 2, 'month' = 2, 'day' = 2), too_many = 'debug') |>
  dplyr::mutate(date = paste0(year, '-', month, '-', day)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%y-%m-%d"))) |>
  ggplot(aes(date, cianos_total))+
  #facet_wrap(vars(type), scales = 'free_y', ncol =1)+
  geom_line(aes(group = type, linetype = type))+
  # geom_line(data = community_eveness_02_red, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2004-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  # scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Date', y = 'Cyanobacteria / Total Bacteria', linetype = '')+
  scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 6), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), 
        strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside',
        panel.border = element_blank()) +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

cianos_asv_fcm_time_plot

# ggsave(cianos_asv_fcm_time_plot,  filename = 'cianos_asv_fcm_time_plot.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 110, units = 'mm')

### test with cianos in PA too ----
cianos_asvs <- asv_tab_10y_3_rel |>
  bind_rows(asv_tab_10y_02_rel) |>
  left_join(tax_bbmo_10y_new) |>
  ungroup() |>
  dplyr::filter(phylum == 'Cyanobacteria') |>
  #dplyr::distinct(order, family, class)
  left_join(m_bbmo_10y) |>
  dplyr::select(date, relative_abundance, sample_id, fraction) |>
  group_by(date, sample_id, fraction) |>
  dplyr::reframe(cianos_total = sum(relative_abundance)) |>
  dplyr::mutate(type = case_when(fraction == '3' ~ 'asvs_3',
                                 fraction == '0.2' ~ 'asvs_02')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_0.2_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_0.2_4408', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id, '_3_4020', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_3_4408', '')) |>
  dplyr::mutate(sample_id_ed = str_replace(sample_id_ed, '_3_4021', '')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::select(-sample_id)

cianos_asv_fcm_time_plot <- cianos_asvs |>
  bind_rows(cianos_fcm) |>
  tidyr::separate_wider_position(col = sample_id_ed, widths = c('station' = 2, 'year' = 2, 'month' = 2, 'day' = 2), too_many = 'debug') |>
  dplyr::mutate(date = paste0(year, '-', month, '-', day)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%y-%m-%d"))) |>
  ggplot(aes(date, cianos_total))+
  #facet_wrap(vars(type), scales = 'free_y', ncol =1)+
  geom_line(aes(group = type, linetype = type))+
  # geom_line(data = community_eveness_02_red, aes(date, community_eveness_rar/2))+
  scale_x_datetime(limits = c(as.POSIXct('2004-01-01'), as.POSIXct('2014-01-01')),
                   expand = c(0,0), 
                   breaks = (by = '1 year'),
                   date_labels = "%Y")+
  # scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~.*2 , name = 'Community Evenness'))+
  labs(x = 'Date', y = 'Cyanobacteria / Total Bacteria', linetype = '')+
  #scale_fill_manual(values = palette_order_assigned_all)+
  #facet_wrap(vars(interaction(family_f, motu_num)), ncol = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 6), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), 
        strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 8), 
        strip.placement = 'outside',
        panel.border = element_blank()) +
  guides(fill = guide_legend(ncol = 3, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')

cianos_asv_fcm_time_plot
# 
# ggsave(cianos_asv_fcm_time_plot,  filename = 'cianos_asv_fcm_time_frac_plot.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 110, units = 'mm')

###
cianos_asvs_fcm_plot <- cianos_asvs |>
  dplyr::filter(type == 'asvs_3') |>
  dplyr::select(-fraction) |>  
  bind_rows(cianos_fcm) |>
  ungroup() |> 
  distinct(sample_id_ed, type,.keep_all = TRUE) |>
  ungroup() |>
  pivot_wider(id_cols = c( 'sample_id_ed'), values_from = 'cianos_total', names_from = 'type', values_fill = NA) |>
  as_tibble() |>
  dplyr::mutate(asvs = as.numeric(asvs_3), 
                # motus = as.numeric(motus),
                fcm = as.numeric(fcm) #, mtags = as.numeric(mtags)
  ) |>
  ggplot(aes(asvs_3, fcm))+
  geom_abline(slope = 1, intercept = 0, color = '#0F3331', linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  #geom_smooth(method = 'loess', color = 'grey')+
  geom_smooth(method = 'lm', color = 'black')+
  stat_cor(aes( #color = 'black',
    label =   paste(..p.label..)), label.x = 0.25,
    #label.y = -0.02,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman'#,
    #position = position_jitter(0.0)
  )+
  stat_cor(aes( label = paste0(..r.label..)),label.x = 0.25, label.y = 0.1, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'spearman')+
  labs(x = 'Ratio ASVs Cyanobacteria / Total Bacteria ', y = 'Ratio Flow cytometry Cyanobacteria / Total Bacteria')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank(),
        #panel.grid.major.y = element_blank(), strip.text = element_text(size = 12),
        legend.position = 'bottom', axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10), strip.background = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 10), 
        panel.border = element_blank(),
        strip.placement = 'outside', aspect.ratio = 12/12)

cianos_asvs_fcm_plot

# ggsave(cianos_asvs_fcm_plot,  filename = 'cianos_asvs_fcm_PA_plot_ed1.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 188, units = 'mm')

### ----  Comparison CARD - FISH counts Clara Ruiz with ASV counts #### ----
card_fish_perc <- read_xlsx('data/env_data/CARDFISH_BLANES 2008-2010 for ONA_ed1.xlsx') |>
  as_tibble()

asv_tab_10y_02_rel |>
  colnames()

tax_bbmo_10y_new |>
  distinct(order, family ) |>
  dplyr::filter(str_detect(order, 'SAR11 clade'))

## SAR 11 clade -----
data <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(str_detect(order, "SAR11 clade") &
                  family %in% c('Clade I', 'Clade II') #', NA
  ) |>
  dplyr::group_by(order, sample_id) |>
  dplyr::reframe(SAR11_abund = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(card_fish_perc, by = c('date'))

## test normality before correlation
shapiro.test(as.numeric(data$SAR11_abund)) # => NORMAL
ggqqplot(as.numeric(data$SAR11_abund))

shapiro.test(as.numeric(data$SAR11_perc)) # => NORMAL
ggqqplot(as.numeric(data$SAR11_perc)) 

SAR11_plot <- data |>
  ggplot(aes(SAR11_abund, SAR11_perc/100))+
  geom_point()+
  stat_cor( aes( 
                                                    label =   paste(..p.label..)), label.x = 0.05,  
           label.y = 0.4,
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'pearson', color = 'black'
           #position = position_jitter(0.0)
  )+
  stat_cor( aes( 
    label =   paste(..r.label..)), label.x = 0.05,  
    label.y = 0.5,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'pearson', color = 'black'
    #position = position_jitter(0.0)
  )+
  geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
  labs(x = 'SAR11 Clade\nAmplicon metabarcoding relative', y = 'SAR11 Clade CARD-FISH relative')+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 6),
        aspect.ratio = 4/4)

SAR11_plot

## GAMMMA -----
data <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(class == 'Gammaproteobacteria') |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(gam_abund = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(card_fish_perc, by = c('date'))

## test normality before correlation
shapiro.test(as.numeric(data$gam_abund)) # => NO NORMAL
ggqqplot(as.numeric(data$gam_abund))
 
gamma_plot <- data |>
  ggplot(aes(gam_abund, GAM_perc/100))+
  geom_point()+
  stat_cor( aes( 
    label =   paste(..p.label..)), label.x = 0.05,  
    label.y = 0.3,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  stat_cor( aes( 
    label =   paste(..r.label..)), label.x = 0.05,  
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
  labs(x = 'Gammaproteobacteria\nAmplicon metabarcoding relative', y = 'Gammaproteobacteria\nCARD-FISH relative')+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 6),
        aspect.ratio = 4/4)

gamma_plot 

## ROS probe ----
data <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family %in% c('Rhodobacteraceae')) |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(ros_abund = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(card_fish_perc, by = c('date'))

## test normality before correlation
shapiro.test(as.numeric(data$ros_abund)) # =>  NORMAL p-value = 0.042
ggqqplot(as.numeric(data$ros_abund))

shapiro.test(as.numeric(data$ROS_perc)) # =>  NO NORMAL p-value = 0.007405
ggqqplot(as.numeric(data$ROS_perc))

rhodo_plot <- data |>
  ggplot(aes(ros_abund, ROS_perc/100))+
  geom_point()+
  stat_cor( aes( 
    label =   paste(..p.label..)), label.x = 0.05,  
    label.y = 0.15,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  stat_cor( aes( 
    label =   paste(..r.label..)), label.x = 0.05,  
    label.y = 0.1,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
  labs(x = 'Rhodobacteraceae\nAmplicon metabarcoding relative', y = 'ROS CARD-FISH relative')+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 6),
        aspect.ratio = 4/4)

rhodo_plot

##   NOR5_perc -----
data <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family == 'Halieaceae') |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(nor5_abund = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(card_fish_perc, by = c('date'))

## test normality before correlation
shapiro.test(as.numeric(data$nor5_abund)) # =>  No NORMAL p-value = 0.007302
ggqqplot(as.numeric(data$nor5_abund))

nor5_plot <-   data  |>
  ggplot(aes(nor5_abund, NOR5_perc/100))+
  geom_point()+
  stat_cor( aes( 
    label =   paste(..p.label..)), label.x = 0.005,  
    label.y = 0.04,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  stat_cor( aes( 
    label =   paste(..r.label..)), label.x = 0.005,  
    label.y = 0.035,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
  labs(x = 'Halieaceae\nAmplicon metabarcoding relative', y = 'NOR5 CARD-FISH relative')+
  theme_bw()+
  theme(panel.border = element_blank(), axis.text = element_text(size = 6),
        aspect.ratio = 4/4)

nor5_plot

# CF_perc  ----
data <-  asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family %in% c('Flavobacteriaceae', 'Cytophagaceae', 'NS11-12 marine group')) |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(cf_abund = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(card_fish_perc, by = c('date'))

## test normality before correlation
shapiro.test(as.numeric(data$cf_abund)) # =>  NO NORMAL p-value = 0.006002
ggqqplot(as.numeric(data$cf_abund))

cf_plot <- data |>
  ggplot(aes(cf_abund, CF_perc/100))+
  geom_point()+
  stat_cor( aes( 
    label =   paste(..p.label..)), label.x = 0.05,  
    label.y = 0.25,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  stat_cor( aes( 
    label =   paste(..r.label..)), label.x = 0.05,  
    label.y = 0.22,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
  labs(x = 'Flavobacteriaceae, Cytophagaceae,\nNS11-12 marine group\nAmplicon metabarcoding relative', y = 'CF319 CARD-FISH relative')+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 5),
        aspect.ratio = 4/4)
  
cf_plot
  
#  CYA_perc

## composition of the results 
corr_card_fish_amplicons <- plot_grid(
  SAR11_plot, gamma_plot,
  rhodo_plot,
  nor5_plot, cf_plot,
  rel_heights = c(1,1,1),
  ncol = 2,
  labels = c('A', 'B', 'C', 'D', 'E'), label_fontface = 'plain',
  label_size = 9,
  hjust = 0.1
)

corr_card_fish_amplicons
# 
# ggsave( plot = corr_card_fish_amplicons, filename = 'corr_card_fish_amplicons.pdf',
#         path = 'results/figures/',
#         width = 180, height = 200, units = 'mm')

## Pseudoabundances of CARD-FISH counts and pseudoabundances of Amplicon Metabarcoding Data ----
citometry_data_pep <- citometry_data_pep |>
  dplyr::mutate(decimal_date_ed = round(decimal_date, digits = 2))

pseudo_card_fish <- card_fish_perc |>
  left_join(m_02, by = c( 'date')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::select('GAM_perc', 'CF_perc', 'SAR11_perc',  'ROS_perc' , 'NOR5_perc', date, bacteria ) |>
  pivot_longer(cols = c('GAM_perc', 'CF_perc', 'SAR11_perc',  'ROS_perc' , 'NOR5_perc'), names_to = 'card_fish_probe', values_to = 'perc') |>
  dplyr::mutate(pseudo_abund_card_fish = as.numeric(perc)*as.numeric(bacteria)/100) |>
  dplyr::select(date, pseudo_abund_card_fish, card_fish_probe)

## pseudoabundances for the amplicon metabarcoding data ---
pseudoabund_amplicon_sar11 <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(str_detect(order, "SAR11 clade") &
                  family %in% c('Clade I', 'Clade II') #', NA
  ) |>
  dplyr::group_by(order, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'SAR11_perc')

pseudo_card_fish |>
  dplyr::filter(card_fish_probe == 'SAR11_perc') |>
  left_join(pseudoabund_amplicon_sar11) |>
  ggplot(aes(pseudo_abund_amplicon, pseudo_abund_card_fish))+
  geom_point()+
  stat_cor( aes( 
    label =   paste(..p.label..)), label.x = 0.05,  
    label.y = 0.4,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'pearson', color = 'black'
    #position = position_jitter(0.0)
  )+
  stat_cor( aes( 
    label =   paste(..r.label..)), label.x = 0.05,  
    label.y = 0.5,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'pearson', color = 'black'
    #position = position_jitter(0.0)
  )+
  geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
  labs(x = 'SAR11 Clade\nAmplicon metabarcoding relative', y = 'SAR11 Clade CARD-FISH relative')+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 6),
        aspect.ratio = 4/4)

pseudoabund_amplicon_altero <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(str_detect(family, "Alteromonadaceae") #&
                 # family %in% c('Clade I', 'Clade II') #', NA
  ) |>
  dplyr::group_by(order, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'GAM_perc')

# other probes
pseudoabund_amplicon_gam <- pseudoabund_amplicon_gamm <- asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(class == 'Gammaproteobacteria') |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'GAM_perc')

pseudoabund_amplicon_ros <-  asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family %in% c('Rhodobacteraceae')) |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'ROS_perc')

pseudoabund_amplicon_nor5 <-  asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family == 'Halieaceae') |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'NOR5_perc')

pseudoabund_amplicon_cf <-  asv_tab_10y_02_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family %in% c('Flavobacteriaceae', 'Cytophagaceae', 'NS11-12 marine group')) |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'CF_perc')

## bind datasets from amplicon pseudocounts ---
pseudoabund_amplicon_tb <- pseudoabund_amplicon_sar11 |>
  bind_rows(pseudoabund_amplicon_gam) |>
  #bind_rpws(pseudoabund_amplicon_altero)
  bind_rows(pseudoabund_amplicon_ros) |>
  bind_rows(pseudoabund_amplicon_nor5) |>
  bind_rows(pseudoabund_amplicon_cf)

pseudoabund_amplicon_card_fish_tb <- pseudoabund_amplicon_tb |>
  left_join(pseudo_card_fish ) 

## without computing correlations to observe 
pseudoabund_amplicon_card_fish_corr_plot <- pseudoabund_amplicon_card_fish_tb |>
  ggplot(aes(pseudo_abund_amplicon, pseudo_abund_card_fish))+
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  geom_point()+
  scale_x_log10(labels = scientific_format())+
  scale_y_log10(labels = scientific_format())+
  stat_cor( aes( 
    label =   paste(..p.label..)), #label.x = 0.05,  
    #label.y = 0.4,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  facet_wrap(vars(card_fish_probe), scales = 'free', ncol = 2)+
  # stat_cor( aes( 
  #   label =   paste(..r.label..)), #label.x = 0.05,  
  #   #label.y = 0.5,
  #   p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'pearson', color = 'black'
  #   #position = position_jitter(0.0)
  # )+
  geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
  labs(x = 'log10(Amplicon metabarcoding pseudoabund)', y = 'log10(CARD-FISH pseudoabund)')+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 6),
        aspect.ratio = 4/4, strip.background = element_rect(fill = 'transparent'))

pseudoabund_amplicon_card_fish_corr_plot

ggsave( plot = pseudoabund_amplicon_card_fish_corr_plot, filename = 'pseudoabund_amplicon_card_fish_corr_plot.pdf',
        path = 'results/figures/',
        width = 180, height = 200, units = 'mm')

## Gammas only compared with Alteromonas ----- 
# pseudoabund_amplicon_tb <- pseudoabund_amplicon_sar11 |>
#   #bind_rows(pseudoabund_amplicon_gam) |>
#   bind_rows(pseudoabund_amplicon_altero) |> # compare Gamma probe with only alteromonas (worst)
#   bind_rows(pseudoabund_amplicon_ros) |>
#   bind_rows(pseudoabund_amplicon_nor5) |>
#   bind_rows(pseudoabund_amplicon_cf)
# 
# pseudoabund_amplicon_card_fish_tb <- pseudoabund_amplicon_tb |>
#   left_join(pseudo_card_fish ) 
# 
# ## without computing correlations to observe 
# pseudoabund_amplicon_card_fish_corr_plot <- pseudoabund_amplicon_card_fish_tb |>
#   ggplot(aes(pseudo_abund_amplicon, pseudo_abund_card_fish))+
#   geom_point()+
#   scale_x_log10()+
#   scale_y_log10()+
#   stat_cor( aes( 
#     label =   paste(..p.label..)), #label.x = 0.05,  
#     #label.y = 0.4,
#     p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
#     #position = position_jitter(0.0)
#   )+
#   facet_wrap(vars(card_fish_probe), scales = 'free', ncol = 2)+
#   # stat_cor( aes( 
#   #   label =   paste(..r.label..)), #label.x = 0.05,  
#   #   #label.y = 0.5,
#   #   p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'pearson', color = 'black'
#   #   #position = position_jitter(0.0)
#   # )+
#   geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
#   labs(x = 'log10(Amplicon metabarcoding pseudoabund)', y = 'log10(CARD-FISH pseudoabund)')+
#   theme_bw()+
#   theme(panel.border = element_blank(),
#         axis.text = element_text(size = 6),
#         aspect.ratio = 4/4, strip.background = element_rect(fill = 'transparent'))
# 
# pseudoabund_amplicon_card_fish_corr_plot

## Correct by number of 16S copy number and then compare again CARD-FISH and AMPLICON abundances -----
## SCRIPT D'EXEMPLE MARKEL:

# Load tables:
# fpath <- "C:/Users/Usuario/Documents/marta/tax.RDS"
# tax <- readRDS(fpath)

#fpath <- "C:/Users/Usuario/Documents/marta/tax.RDS"
tax <- tax_bbmo_10y_new #mirar d'on surt aquesta taula

fpath <- "../../REMEI/REMEI_data/intermediate_files/16S_rna_copy_num_normalization/cnadjusted_remei_asvs.hier.tsv"
cn <- read.table(file = fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cn |>
  colnames()

cn$remei_asvs.fasta

# We need to add $cnumber to tax, using the data that gets as close as possible to the taxonomy of each ASV in tax.

# We need to somehow pair the data from both table considering the taxonomic resolution available for each ASV.
# If phylum, phylum; if genus, genus; etc. So if we have domainA_phylumB, we need to match the cnumber 
# corresponding to domainA_phylumB, not only to domainA.

# Create a modified lineage variable where only the actual taxonomy is left (remove "root", "domain", etc.)
cn$id.lin <- gsub(pattern = "Root;rootrank;|domain;|phylum;|class;|order;|family;|genus;", replacement = "", x = cn$lineage)
# check:
# cn$id.lin[4]

# Remove the last ";" and replace the others with "_" so that the id is "NAME_NAME_NAME" etc.
# Well no there are some taxa with "_" already in their name so better use mmm "." maybe?
cn$id.lin <- gsub(";", ".", gsub(";$", "", cn$id.lin))
# check:
# cn$id.lin[4]

# Similarly, merge the taxonomic assignations in tax the same way:
tax %>%
  colnames()

tax$id.tax <- apply(tax, MARGIN = 1, FUN = function(x){paste(x[c("domain", "phylum", "class", "order", "family", "genus")], collapse = ".")})
# check:
# tax$id.tax[1]

# Remove the ".NA"s introduced by unidentified taxonomy levels:
# check:
# tax$id.tax[4]
# gsub(".NA", "", tax$id.tax[4])
tax$id.tax <- gsub(".NA", "", tax$id.tax)

# Now they have shared IDs.
#
# Add cnumber to taxa:
tax$cnumber <- cn$cnumber[match(tax$id.tax, cn$id.lin)]

#no funciona el match faig aquesta funció
tax <- tax %>%
  left_join(cn, by = c('id.tax' = 'id.lin')) %>%
  dplyr::mutate(cnumber = remei_asvs.fasta)

# How many ASVs have cnumber?
table(!is.na(tax$cnumber))
# From the top 200:
table(!is.na(tax$cnumber[1:200]))

# And in terms of unique taxa? (to avoid repetitions as ASVs may share the same taxonomy)
table(!is.na(tax$cnumber[!duplicated(tax$id.tax)]))

# So many are missing... An option would be that if there is no match, seek values of previous taxonomic
# levels. E.g., if a genus has no cnumber, see if its family has.

# But first, write down which ASVs have exact matches:
tax$cnumber.exact <- !is.na(tax$cnumber)

cn <- 
  cn %>%
  mutate(cnumber = remei_asvs.fasta)

# Find close results. Now this is more tricky.
# We can maybe do a loop and go entry by entry. If no exact match, crop the id one level, which equals going one up.
# That could work.
# can check with i = 3 (exact match) 6 (genus to family), 46 (genus to order)
for(i in 1:nrow(tax)){
  
  # Is there a cnumber entry for ASV i?
  if(tax$id.tax[i] %in% cn$id.lin){
    
    # If there is, get the cnumber:
    tax$cnumber[i] <- cn$cnumber[cn$id.lin==tax$id.tax[i]]
    
    # Write down match type. In this case, it is an exact match:
    tax$cnumber.match[i] <- "exact"
    
  } else {
    
    # If not, find the closest value.
    # We can do that by cropping the ID one level back and checking again.
    # But never going over "order".
    # To know if we exceed order or not, we can check the id length.
    # Ids look like:
    # domain.phylum.class.order.family.genus
    # So when split by ".", it has to have at least four elements (to be at least order):
    # domain.phylum.class.order
    
    # How long is the ID we are checking?    
    l <- length(strsplit(tax$id.tax[i], split = "\\.")[[1]])
    
    # We can only seek for alternatives if its last level is family or genus, i.e., l>4
    if(l==5){
      
      # If it is family (l=5), we can only try going one level up until reaching order:
      alt <- paste(strsplit(tax$id.tax[i], split = "\\.")[[1]][1:4], collapse = ".")
      
      # Does its order have an entry?
      if(alt %in% cn$id.lin){ 
        
        # if yes:
        tax$cnumber[i] <- cn$cnumber[cn$id.lin==alt]
        tax$cnumber.match[i] <- "from.order"
        
      } else {
        
        # if not, set default value to 1.
        tax$cnumber[i] <- 1
        tax$cnumber.match[i] <- "default.to.1"
        
      }
      
    } else if(l==6){
      
      # If it is genus (l=6), we can try going up twice.
      # First, try with family:
      alt <- paste(strsplit(tax$id.tax[i], split = "\\.")[[1]][1:5], collapse = ".")
      
      # Does its family have an entry?
      if(alt %in% cn$id.lin){ 
        
        # if yes, 
        tax$cnumber[i] <- cn$cnumber[cn$id.lin==alt]
        tax$cnumber.match[i] <- "from.family"
        
      } else {
        
        # if not, try its order:
        alt2 <- paste(strsplit(tax$id.tax[i], split = "\\.")[[1]][1:4], collapse = ".")
        
        if(alt2 %in% cn$id.lin){ 
          
          # if yes:
          tax$cnumber[i] <- cn$cnumber[cn$id.lin==alt2]
          tax$cnumber.match[i] <- "from.order"
          
        } else {
          
          # if not, set default value to 1.
          tax$cnumber[i] <- 1
          tax$cnumber.match[i] <- "default.to.1"
          
        }
      }
    } else {
      
      # If it already was order (or class/phylum) (l<=4), do not seek alternatives.
      # Set instead to 1 as default.
      tax$cnumber[i] <- 1
      tax$cnumber.match[i] <- "default.to.1"
      
    }
  }
}


## add this to find information for genus that have a different taxonomic rute -------
tax$cnumber_match_genus<-cn[match(tax$genus, cn$name),7]
cn_curated<-NULL
for (i in 1:dim(tax)[1]){
  if(is.na(tax[i,12]))
    cn_curated[i]<-tax[i, 9]
  else  cn_curated[i]<-tax[i,12]
}
cn_curated
tax$cn_curated<-cn_curated
cnumber_match_curated<-NULL
for (i in 1:dim(tax)[1]){
  if(is.na(tax[i,12]))
    cnumber_match_curated[i]<-tax[i, 11]
  else  cnumber_match_curated[i]<- "exact"
}
cnumber_match_curated
tax$cnumber_match_curated <- cnumber_match_curated

## in case cnumber match curated !na use this value, otherwise use cnumber
tax %>%
  colnames()

tax <- tax %>%
  mutate(cn_final = case_when(!is.na(cnumber_match_genus) ~ cnumber_match_genus,
                              is.na(cnumber_match_genus) ~ cnumber))

## count nº of ASVs that have 1 as copy number (default)
tax %>%
  filter(cn_final == 1) %$%
  asv_num %>%
  unique() #3687

tax %>%
  filter(cn_final != 1) %$%
  asv_num %>%
  unique() #4162

## now normalize reads by the number of 16S rRNA copy numbers 
tax |>
  colnames()

tax_filt <- tax %>%
  dplyr::select(asv_num, cn_final)

asv_tab_cn <- tax %>%
  dplyr::select(asv_num, cn_final) %>%
  right_join(asv_tab_10y_02_rel, by = c('asv_num' ))

asv_tab_cn %>%
  dim() == asv_tab_10y_02_rel |>
  dim()

asv_tab_cn %>%
  glimpse()

asv_tab_cn_nor <- asv_tab_cn %>%
  dplyr::select(-relative_abundance) |>
  # Convert all columns except 'asv_num' and 'sample_id' to numeric
  dplyr::mutate(across(-c(asv_num, sample_id), as.numeric)) %>%
  # Normalize numeric columns, excluding 'asv_num' and 'cn_final'
  mutate(across(where(is.numeric) & -c(asv_num, cn_final), ~ . / cn_final)) %>%
  # Retain 'sample_id' explicitly
  dplyr::select(sample_id, everything(), -asv_num) %>%
  # Combine with 'tax_filt' after ensuring both data frames have compatible dimensions
  cbind(tax_filt)

##  ----- relative abundance for reads normalized by 16S rRNA copy number -----
asv_tab_cn_nor %>%
  dim()

asv_tab_cn_nor %$%
  unique(asv_num)

asv_tab_cn_nor_rel <- asv_tab_cn %>%
  dplyr::select(-relative_abundance) |>
  dplyr::mutate(reads = reads/cn_final) |>
  calculate_rel_abund(group_cols = sample_id)

asv_tab_cn_nor_rel %>%
  dplyr::group_by(sample_id) |>
  reframe(n = sum(relative_abundance)) ## !! it's not 1
  # as_tibble() %>%
  # dplyr::select(-V309) %>%
 
# other probes
pseudoabund_amplicon_gam <- asv_tab_cn_nor_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(class == 'Gammaproteobacteria') |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'GAM_perc')

pseudoabund_amplicon_altero <- asv_tab_cn_nor_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(str_detect(genus, "Alteromonas") #&
                # family %in% c('Clade I', 'Clade II') #', NA
  ) |>
  dplyr::group_by(order, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'GAM_perc')

pseudoabund_amplicon_ros <-  asv_tab_cn_nor_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family %in% c('Rhodobacteraceae')) |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'ROS_perc')

pseudoabund_amplicon_nor5 <-  asv_tab_cn_nor_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family == 'Halieaceae') |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'NOR5_perc')

pseudoabund_amplicon_cf <-  asv_tab_cn_nor_rel |>
  left_join(tax_bbmo_10y_new, by = c('asv_num')) |>
  dplyr::filter(family %in% c('Flavobacteriaceae', 'Cytophagaceae', 'NS11-12 marine group')) |>
  dplyr::group_by(class, sample_id) |>
  dplyr::reframe(relative_abundance = sum(relative_abundance)) |>
  left_join(m_02, by = c('sample_id')) |>
  dplyr::filter(year %in% c('2008', '2009', '2010')) |>
  left_join(citometry_data_pep, by = c('decimal_date'= 'decimal_date_ed')) |>
  dplyr::mutate(pseudo_abund_amplicon = as.numeric(relative_abundance)*as.numeric(bacteria)) |>
  dplyr::select(date, pseudo_abund_amplicon) |>
  dplyr::mutate(card_fish_probe = 'CF_perc')

## bind datasets from amplicon pseudocounts ---
pseudoabund_amplicon_tb <- pseudoabund_amplicon_sar11 |>
  bind_rows(pseudoabund_amplicon_gam) |>
  #bind_rows(pseudoabund_amplicon_altero) |> #it didn't improve filtering by only Alteromonas
  bind_rows(pseudoabund_amplicon_ros) |>
  bind_rows(pseudoabund_amplicon_nor5) |>
  bind_rows(pseudoabund_amplicon_cf)

pseudoabund_amplicon_card_fish_tb <- pseudoabund_amplicon_tb |>
  left_join(pseudo_card_fish ) 

## without computing correlations to observe 
pseudoabund_amplicon_card_fish_corr_plot <- pseudoabund_amplicon_card_fish_tb |>
  ggplot(aes(pseudo_abund_amplicon, pseudo_abund_card_fish))+
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  geom_point()+
  scale_x_log10(labels = scientific_format())+
  scale_y_log10(labels = scientific_format())+
  stat_cor( aes( 
    label =   paste(..p.label..)), #label.x = 0.05,  
    #label.y = 0.4,
    p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman', color = 'black'
    #position = position_jitter(0.0)
  )+
  facet_wrap(vars(card_fish_probe), scales = 'free', ncol = 2)+
  # stat_cor( aes( 
  #   label =   paste(..r.label..)), #label.x = 0.05,  
  #   #label.y = 0.5,
  #   p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'pearson', color = 'black'
  #   #position = position_jitter(0.0)
  # )+
  geom_smooth(method = 'lm', color = 'black', fill = 'grey')+
  labs(x = 'log10(Amplicon metabarcoding pseudoabund)', y = 'log10(CARD-FISH pseudoabund)')+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 6),
        aspect.ratio = 4/4, strip.background = element_rect(fill = 'transparent'))

pseudoabund_amplicon_card_fish_corr_plot

# ggsave( plot = pseudoabund_amplicon_card_fish_corr_plot, filename = 'pseudoabund_amplicon_card_fish_corr_plot_16S_nor.pdf',
#         path = 'results/figures/',
#         width = 180, height = 200, units = 'mm')

bloo_all_types_summary_tax |>
  dplyr::mutate(recurrency = case_when(
    asv_num == 'asv15' & fraction == '0.2' ~ 'recurrent',
    TRUE ~ recurrency
  )) |>
  group_by(fraction, recurrency) |>
  dplyr::reframe(n = n())


get_citations(("WaveletComp"), out.dir ='../../Bibliografia/')
get_citations(("grateful"), out.dir ='../../Bibliografia/')

## Which blooms did we identify in the FL fraction during 2009-2014 ----

## Do we find any gene that has an anomaly in their abundance? -----
