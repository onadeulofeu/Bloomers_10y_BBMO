# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# upload packages 
library(tidyverse)
library(phyloseq)
library(speedyseq)

##update taxonomy with database 138.1 silva
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DECIPHER")

# setwd("~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/")
bbmo_10y <-readRDS("~/Documentos/Doctorat/BBMO/BBMO_bloomers/data/blphy10years.rds") 
head(bbmo_10y)

# bbmo_10y@otu_table |>
#   View()

#new taxonomy created with the database SILVA 138
new_tax <-  readRDS('data/03_tax_assignation/raw/devotes_all_tax_assignation.rds') |>
  as_tibble(rownames = 'sequence')

tax_bbmo_10y_old <- bbmo_10y@tax_table |>
  #mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(bbmo_10y@otu_table))) |> #ja té uns nº d'OTUs
  as_tibble()

colnames(tax_bbmo_10y_old) <- c("asv_num", "kingdom", "phylum", "class", "order", "family", "genus",
                                "species", "curated", "otu_corr","seq")

tax_bbmo_10y_new <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax, by = c('seq' = 'sequence'))
# 
# ##UPLOAD BLOOMERS DATA-----
# asv_tab_all_bloo_z_tax_old <- read.csv2('data/asv_tab_all_bloo_z_tax.csv')
# asv_tab_all_bloo_z_tax <- read.csv2('data/asv_tab_all_bloo_z_tax_new.csv')
# 
# #I get ASVs which don't have domain, which is strange, as they get classified with blast.
# asv_tab_all_bloo_z_tax_old |>
#   dplyr::filter(is.na(phylum))
# 
# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(is.na(phylum))
# 
# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(asv_num == 'asv80') |>
#   distinct(seq) == asv_tab_all_bloo_z_tax_old |>
#   dplyr::filter(asv_num == 'asv80') |>
#   distinct(seq)
# 
# asv_tab_all_bloo_z_tax_old |>
#   dplyr::filter(asv_num == 'asv80') |>
# distinct(seq,  phylum, class, order, family)
# 
# #Podria ser que el decihper sigui molt més restrictiu del que és el assign tax i per això em queden més ASVs sense classificar. 
# # De totes maneres no tenir domain em sembla bastant... Podria provar de fer una assignació amb assign_tax de dada2 i veure quins 
# # resultats obtenim.
# 
# ##Test that the match performed between the two databases is correct. 
# y <- asv_tab_all_bloo_z_tax |> 
#   dplyr::select(seq, asv_num) |>
#   distinct(seq, asv_num) |>
#   arrange(asv_num)
# 
# y |>
#   dim()
# 
# x <- asv_tab_all_bloo_z_tax_old |> 
#   dplyr::select(seq, asv_num) |>
#   distinct(seq, asv_num) |>
#   arrange(asv_num)
# 
# y |>
#   dim()
# 
# x == y 
# 
# ## it seems that the new taxonomy is correct, and that I mantained the asv_num for the new dataset. 
# ## check that Glaciecola (asv11) is present in the dataset.
# asv_tab_all_bloo_z_tax |>
#   dplyr::filter( str_detect(seq, 
#                             pattern = 'TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGTTGTGAGGAAAGTTTAGTAGTTAATACCTGCTAGATGTGACGTTAGCAACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTTAAGCTAGATGTGAAAGCCCCGCGCTCAACGTGGGAGGGTCATTTAGAACTGGCAGACTAGAGTCTTGGAGAGGGGAGTGGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAACATCAATGGCGAAGGCAACTCCCTGGCCAAAGACTGACGCTCATGTGCGAAAGTGTGGGTAGCGAACAGG',
#   ))
# 
# asv_tab_all_bloo_z_tax$asv_num
# 
# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(asv_num == 'asv11') |>
#   distinct(seq)

#  We USED ASSIGN TAX INSTEAD OF DECIPHER -----------------
# I used the assign tax function from dada2 to assign the taxonomy. For this purpose I don't need a .RData
# in this case I use the database SILVA 138.1 the last one.
new_tax_assign <-  readRDS('data/03_tax_assignation/devotes_all_assign_tax_assignation_v2.rds') |>
  as_tibble(rownames = 'sequence')

new_tax_assign |>
  colnames()

## I compare it with the old one
new_tax <-  readRDS('data/03_tax_assignation/raw/devotes_all_tax_assignation.rds') |>
  as_tibble(rownames = 'sequence')

new_tax |>
  colnames()

differences_new_tax_dada2_decipher <- new_tax_assign |>
  left_join(new_tax, by = 'sequence') |>
  dplyr::mutate(domain_dif = (Kingdom == domain),
                phylum_dif = (Phylum == phylum),
                class_dif = (Class == class),
                order_dif = (Order == order),
                family_dif = (Family == family),
                genus_dif = (Genus == genus)) |>
  summary()

## Conclusions: after phylum there's FALSE assignation between both strategies and appear NAs. Let's explore the unclassifieds at the different taxonomic
## levels.
tax_bbmo_10y_new_assign <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax_assign, by = c('seq' = 'sequence'))

tax_bbmo_10y_new <- tax_bbmo_10y_old |>
  dplyr::select(asv_num, seq) |>
  left_join(new_tax, by = c('seq' = 'sequence'))

tax_bbmo_10y_new |>
  dplyr::filter(is.na(phylum))

tax_bbmo_10y_new_assign |>
  dplyr::filter(is.na(Phylum))

tax_bbmo_10y_new |>
  dplyr::filter(is.na(class))

tax_bbmo_10y_new_assign |>
  dplyr::filter(is.na(Class))

tax_bbmo_10y_new_assign |>
  dplyr::filter(asv_num == 'asv24713') %$%
  seq

##calculate the percentage o NAs between one table and the other, I keep both datasets just in case.

# #prepare an asv_tab that I can use to update the taxonomy
# ## add seqs rather than asv_num
# asv_tab_bbmo_10y_l <- bbmo_10y@otu_table |>
#   as_tibble()
# 
# bbmo_10y_seqtab <- asv_tab_bbmo_10y_l |>
#   left_join(tax_bbmo_10y, by = c('.otu' = 'asv_num')) |>
#   dplyr::select(seq, .sample, .abundance) |>
#   pivot_wider(id_cols = .sample, names_from =seq, values_from = .abundance) |>
#   as.data.frame()
# 
# tax_bbmo_10y %$%
#   seq |>
#   duplicated() |>
#   summary()
# 
# ##gestionar el tema dels duplicats
# asv_tab_bbmo_10y_l |> 
#   dplyr::group_by(.sample, .otu) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   dplyr::filter(n > 1) 
# 
# asv_tab_bbmo_10y_l |> 
#   left_join(tax_bbmo_10y, by = c('.otu' = 'asv_num')) |>
#   dplyr::select(seq, .sample, .abundance) |>
#   dplyr::group_by(.sample, seq) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   dplyr::filter(n > 1) 
# 
# asv_tab_bbmo_10y_l |>
#   left_join(tax_bbmo_10y, by = c('.otu' = 'asv_num')) |>
#   dplyr::select(seq, .sample, .abundance) |>
#   pivot_wider(id_cols = .sample, names_from = seq, values_from = .abundance, values_fn = list) |>
#   as.data.frame()
# 
# row.names(bbmo_10y_seqtab) <- bbmo_10y_seqtab$.sample
# 
# bbmo_10y_seqtab <- bbmo_10y_seqtab |>
#   #rownames_to_column() |>
#   dplyr::select(-.sample)
# 
# bbmo_10y_seqtab |>
# class()
# bbmo_10y_seqtab_final <- bbmo_10y_seqtab_final |>
#   as.matrix()
# 
# # asv_tab_bbmo_10y_l %>%
# #   left_join(tax_bbmo_10y, by = c('.otu' = 'asv_num')) |>
# #   dplyr::group_by(seq, .sample) %>%
# #   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
# #   dplyr::filter(n > 1L) 
# 
# 
# saveRDS(bbmo_10y_seqtab, 'bbmo_10y_seqtab_final.rds')
# 
# ##create data with the same format as the one we used in REMEI see if it works
# src(remei_1_2_pool_seqtab_final)
# src(bbmo_10y_seqtab_final)
# class(remei_1_2_pool_seqtab_final)
# class(bbmo_10y_seqtab_final)
# 
# bbmo_10y_seqtab_final
