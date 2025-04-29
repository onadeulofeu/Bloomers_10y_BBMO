# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# upload packages ------
library(tidyverse)
library(ggtreeExtra)
library(ggtree)
library(ggplot2)
library(tidytree)
library(ggpubr)
library(ggridges)

# packages version ----
# R version 4.2.3 (2023-03-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS 14.5
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# dendextend_1.17.1 ggridges_0.5.4    ggpubr_0.6.0      tidytree_0.4.6    ggtree_3.11.0     ggtreeExtra_1.8.1 ape_5.7-1        
# lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2       readr_2.1.4       tidyr_1.3.1      
# tibble_3.2.1      ggplot2_3.4.4     tidyverse_2.0.0   phyloseq_1.42.0  
# 
# loaded via a namespace (and not attached):
# nlme_3.1-163           bitops_1.0-7           fs_1.6.3               GenomeInfoDb_1.34.9    backports_1.4.1        tools_4.2.3           
# utf8_1.2.4             R6_2.5.1               vegan_2.6-4            lazyeval_0.2.2         BiocGenerics_0.44.0    mgcv_1.9-0            
# colorspace_2.1-0       permute_0.9-7          rhdf5filters_1.10.1    ade4_1.7-22            withr_3.0.0            tidyselect_1.2.0      
# gridExtra_2.3          compiler_4.2.3         cli_3.6.2              Biobase_2.58.0         FD_1.0-12.3            scales_1.3.0          
# digest_0.6.34          yulab.utils_0.1.4      XVector_0.38.0         pkgconfig_2.0.3        fastmap_1.1.1          rlang_1.1.3           
# rstudioapi_0.15.0      gridGraphics_0.5-1     generics_0.1.3         jsonlite_1.8.8         car_3.1-2              RCurl_1.98-1.12       
# magrittr_2.0.3         ggplotify_0.1.2        GenomeInfoDbData_1.2.9 biomformat_1.26.0      patchwork_1.2.0        Matrix_1.6-1.1        
# Rcpp_1.0.12            munsell_0.5.0          S4Vectors_0.36.2       Rhdf5lib_1.20.0        fansi_1.0.6            viridis_0.6.4         
# abind_1.4-5            ggnewscale_0.4.10      lifecycle_1.0.4        stringi_1.8.3          carData_3.0-5          MASS_7.3-60           
# zlibbioc_1.44.0        rhdf5_2.42.1           plyr_1.8.9             grid_4.2.3             parallel_4.2.3         crayon_1.5.2          
# lattice_0.22-5         Biostrings_2.66.0      cowplot_1.1.1          splines_4.2.3          multtest_2.54.0        hms_1.1.3             
# pillar_1.9.0           igraph_1.5.1           ggsignif_0.6.4         reshape2_1.4.4         codetools_0.2-19       stats4_4.2.3          
# magic_1.6-1            glue_1.7.0             ggfun_0.1.4            data.table_1.14.8      vctrs_0.6.5            treeio_1.27.0.002     
# tzdb_0.4.0             foreach_1.5.2          gtable_0.3.4           cachem_1.0.8           broom_1.0.5            rstatix_0.7.2         
# viridisLite_0.4.2      survival_3.5-7         geometry_0.4.7         iterators_1.0.14       aplot_0.2.2            memoise_2.0.1         
# IRanges_2.32.0         cluster_2.1.4          timechange_0.2.0  

## Closely related bloomers distance ----
tax_bbmo_10y_new |>
  dplyr::filter(asv_num == 'asv7') %$%
  seq

tax_bbmo_10y_new |>
  dplyr::filter(asv_num == 'asv1') %$%
  seq

## 3 different bp

tax_bbmo_10y_new |>
  dplyr::filter(asv_num == 'asv4') %$%
  seq

tax_bbmo_10y_new |>
  dplyr::filter(asv_num == 'asv31') %$%
  seq

## 1 different bp

tax_bbmo_10y_new |>
  dplyr::filter(asv_num == 'asv17') %$%
  seq

tax_bbmo_10y_new |>
  dplyr::filter(asv_num == 'asv77') %$%
  seq

## 1 different bp

# Previously I need to create a fasta file with the sequences of my potential bloomers----
# Aligment
##https://www.genome.jp/tools-bin/clustalw

# ## add fasta seqs, DNA format, and Select Weight Matrix: CLUSTALW (for DNA)
# 
#  asv_num_seq_bloo <- tax |>
#   left_join(tax_bbmo_10y_new, by = 'asv_num') |>
#   dplyr::select(seq, asv_num) |>
#   dplyr::mutate(seq = as.character(seq),
#                 asv_num = as.character(asv_num))
# 
# # Assuming your tibble is named 'sequences_tibble' with columns 'sequence_name' and 'sequence'
# 
# # Open a connection to write the output to a file
# output_file <- "sequences.fasta"
# file_conn <- file(output_file, "w")
# 
# # Iterate over each row in the tibble and write to the file in FASTA format
# for (i in 1:nrow(asv_num_seq_bloo)) {
#   # Extract sequence name and sequence
#   seq_name <- as.character(asv_num_seq_bloo$asv_num[i])  # Extract character string from the tibble
#   seq <- as.character(asv_num_seq_bloo$seq[i])  # Ensure seq is a character vector
#   
#   # Write to file in FASTA format
#   cat(">", seq_name, "\n", seq, "\n", file = file_conn, sep = "")
# }
# 
# # Close the file connection
# close(file_conn)
# 
# # Print a message indicating the file has been written
# cat("FASTA file '", output_file, "' has been created.\n")
# 
# ## Run the code in marbits or in local, look at README file in the folder----
# 
# # module load
# # 
# # module load raxml-ng/0.9.0
# # 
# # raxml-ng --check --msa /alignment_remei_100_long_trimmed.fasta --model GTR+I+G -seed 123 --prefix remei100 ## for checking sequences 
# # 
# # raxml-ng --msa alignment_remei_100_long_trimmed.fasta --model GTR+I+G —-outgroup NR_074309.1_Synechococcus_elongatus_PCC_6301 -seed 123 --threads 2 --prefix remei100 ## for doing the main tree
# # 
# # raxml-ng --rfdist --tree remei100.raxml.mlTrees --prefix remei100 ## para comprobar topologias
# # 
# # raxml-ng --bootstrap --msa alignment_remei_100_long_trimmed.fasta --model GTR+I+G --outgroup NR_074309.1_Synechococcus_elongatus_PCC_6301 -seed 123 --threads 2 --bs-trees 600 --prefix remei100 ## for making bootstraps
# # 
# # raxml-ng --bsconverge --bs-trees remei100.raxml.bootstraps --prefix remei100 -seed 123 --threads 1 --bs-cutoff 0.01 ## para ver si las bootstraps convergen
# # 
# # raxml-ng --support --tree remei100.raxml.bestTree --bs-trees remei100.raxml.bootstraps --prefix remei100 --threads 2 ## para juntar las bootstraps con el tree

# Visualization of our tree ----
tree <- ggtree::read.tree('data/raxml/bloo_bbmo.raxml.support')

tree |>
  ggplot() + 
  geom_tree() + 
  theme_tree2() +
  #geom_treescale()+
  geom_tiplab(align = T)
  #geom_tippoint()

# tree |>
#   ggplot() + 
#   geom_tree() + 
#   theme_tree2() +
#   #geom_treescale()+
#   geom_tiplab(align = T)
# 
# tax <- asv_tab_all_bloo_z_tax |>
#   dplyr::select(asv_num, phylum, class, order, family, genus) |>
#   distinct()
# 
# metadata <- bloo_02 |>
#   dplyr::mutate(fraction = '0.2') |>
#   bind_rows(bloo_3) |>
#   dplyr::mutate(fraction = case_when(is.na(fraction)~ '3',
#                                      !is.na(fraction) ~ fraction)) |>
#   group_by(value) |>
#   dplyr::mutate(detect = n()) |>
#   dplyr::mutate(detect_both = case_when(detect == '2' ~ 'both',
#                                         detect == '1' ~ fraction)) |>
#   group_by( value) |>
#   dplyr::select(detect_both, value) |>
#   distinct() |>
#   rename(asv_num = value) |>
#   left_join(tax)
# 
# merged_data <- merge(as_tibble_col(tree$tip.label, column_name = 'asv_num'), metadata, by = "asv_num") |>
#   dplyr::mutate(family_num = paste0(family,' ', asv_num))
# 
#  tree_plot <- ggtree(tree, layout='dendrogram', branch.length='none') %<+% ## this is used to assign new data to the plot
#   merged_data +
#   geom_tree(aes(color=class))+
#   scale_color_manual(values = palette_class_assigned)+
#   geom_tiplab( aes(label=family), size=2, align=TRUE) +
#     labs(color = 'Class')+
#   #geom_treescale()+
#   #geom_tiplab(align = T)+
#   theme_tree2()
#   #geom_treescale(x = 10, y = 360)
#  
#  ## prepare the data if they are blooming on the FL or the PA fraction
#  detect_both <- merged_data |>
#    dplyr::select(detect_both, asv_num) |>
#    dplyr::mutate(blooming_fraction = 1,
#                  detect_both = str_replace(detect_both, '3', 'PA')) |>
#    dplyr::mutate(detect_both = str_replace(detect_both, '0.2', 'FL')) |>
#    pivot_wider(id_cols = asv_num, values_fill = 0, names_from = detect_both, values_from = blooming_fraction)
#  
#  detect_both <- detect_both |>
#    dplyr::mutate_all(as.character) |>
#    dplyr::mutate(PA = case_when(both == '1' ~ '1',
#                                 both == '0' ~ PA),
#                  FL =  case_when(both == '1' ~ '1',
#                                  both == '0' ~ FL)) |>
#    dplyr::select(-both)   |>
#    column_to_rownames(var = 'asv_num') 
#  
#  palette_fraction_bw <- c(  '0' = 'white', '1' = 'black' )
#  
 tips_to_remove <- c("asv2", "asv3", "asv5", "asv8") ## discarded as potential bloomers for our dataset
 tree <- drop.tip(tree, tips_to_remove)

 # Plot the modified tree
 plot(tree)

 # Plot the dendrogram with colored branches based on 'class'
 tree_plot <- ggtree(tree, branch.length='none') %<+% ## this is used to assign new data to the plot
  merged_data +
  geom_tree() + #aes(color=class)
  scale_color_manual(values = palette_order_assigned_bloo) +
   geom_tiplab(aes(label = family_num), size=3.5, align=TRUE) + 
   geom_tippoint(aes(color=order), alpha = 0.8)+
  labs(color = 'Order') +
  theme_tree2()+
   theme(legend.position = "none",  # Remove legend
         panel.grid.major = element_blank(),  # Remove grid lines
         panel.grid.minor = element_blank(),  # Remove grid lines
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.line.x = element_blank(),
         axis.ticks.x = element_blank()) 

# Create the heatmap colored by 'detect_both' variable
heatmap <- gheatmap(tree_plot, detect_both, offset=7, width = .1, color=NA, font.size = 3) + 
  scale_fill_manual(values = palette_fraction_bw, name="Potential\nblooming\nin fraction", labels=c("0" = "No", "1" = "Yes"))+
  theme(legend.position = 'bottom', text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines")
        )

print(heatmap)

#### Based on wavelets analysis I create a label for each taxa and add the information to the phylogenetic tree ----
bloo_type_biased_all ##here I have their maximum coefficient. Maybe the lowest ones are not significant and other criteria need to be followed to decide what do we trust

bloo_type_biased_all <- read.csv('data/bloo_type_biased_all_checked.csv') # not sure if this is the table that i would like to use
wavelets_result_02 <- read.csv('data/wavelets_analysis/wavelets_result_ed_tibble_tax_02_biased_red.csv')
wavelets_result_3 <- read.csv('data/wavelets_analysis/wavelets_result_ed_tibble_tax_3_biased_red.csv')

wavelets_result_ed_tibble_tax_02_biased_red_coeff <- wavelets_result_02 %>%
  group_by(asv_num, wavelets_transformation) %>%
  dplyr::filter(!is.na(wavelets_result_ed)) |>
  dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) |>
  dplyr::mutate(fraction = '0.2') |>
  dplyr::mutate(wavelets_fraction = paste0(wavelets_transformation,'_', fraction)) |>
  dplyr::select(-fraction, -wavelets_transformation)

wavelets_result_ed_tibble_tax_3_biased_red_coeff <- wavelets_result_3 %>%
  group_by(asv_num, wavelets_transformation) %>%
  dplyr::filter(!is.na(wavelets_result_ed)) |>
  dplyr::summarize(coefficients = sqrt(sum(wavelets_result_ed^2))) |>
  dplyr::mutate(fraction = '3') |>
  dplyr::mutate(wavelets_fraction = paste0(wavelets_transformation,'_', fraction)) |>
  dplyr::select(-fraction, -wavelets_transformation)

wavelets_result_ed_tibble_biased_red_coeff_all <-  wavelets_result_ed_tibble_tax_3_biased_red_coeff |>
  bind_rows(wavelets_result_ed_tibble_tax_02_biased_red_coeff) |>
  pivot_wider(id_cols = asv_num, values_from = coefficients, names_from = wavelets_fraction, values_fill = 0)

## I prepare a table similar to the one I have for the fraction at which they bloom or not bloom.
type_bloo <- bloo_type_biased_all |>
  dplyr::mutate(type_fraction = paste0(fraction, '.', bloomer_type)) |>
  dplyr::select(asv_num, coefficients, type_fraction) |>
  pivot_wider(id_cols = asv_num, values_from = coefficients, names_from = type_fraction, values_fill = 0) |>
  as_tibble() |>
  column_to_rownames(var = 'asv_num') 

palette_gradient_bw <- c(  'white', 'black' ) 

type_bloo_heatmap <- gheatmap(tree_plot, type_bloo, offset=12, width = .1, color=NA, font.size = 3) + 
  scale_fill_gradientn(colours = palette_gradient_bw,
                       name="Bloomer's type")+
  theme(legend.position = 'bottom', text = element_text(size = 6), # labels=c("0" = "No", "1" = "Yes")
        legend.key.size = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 1, angle =180))

combined_plot <- tree_plot +
  heatmap+
type_bloo_heatmap

wavelets_result_ed_tibble_biased_red_coeff_all <- wavelets_result_ed_tibble_biased_red_coeff_all |>
  as_tibble() |>
  column_to_rownames(var = 'asv_num') 

# Create the tree plot
tree_plot <- ggtree(tree, branch.length='none') %<+% 
  merged_data +
  geom_tree() + #aes(color=class)
  scale_color_manual(values = palette_order_assigned_bloo) +
  geom_tiplab(aes(label = family_num), size=2.5, align=TRUE) + 
  geom_tippoint(aes(color=order), alpha = 1)+
  labs(color = 'Order') +
  theme_tree2()+
  theme(legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),  # Remove grid lines
        panel.grid.minor = element_blank(),  # Remove grid lines
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 5))

# Create the first heatmap
heatmap1 <- gheatmap(tree_plot, wavelets_result_ed_tibble_biased_red_coeff_all, offset=11, width = 2, color=NA, font.size = 2) + 
  scale_fill_gradientn(colours = palette_gradient_bw,
                       name="Bloomer's\nwavelets\ncoefficients")+
  theme(legend.position = 'bottom', text = element_text(size = 5), # labels=c("0" = "No", "1" = "Yes")
        legend.key.size = unit(0.5, "lines"))

# Print the combined plot
print(heatmap1)

# ggsave(heatmap1, filename = 'phylogenetic_tree_wavelets_coeff.pdf',
#        path = 'Results/Figures/',
#        width = 230, height = 200, units = 'mm')

### Remove d2 and d4 signals (no information)----
wavelets_result_ed_tibble_biased_red_coeff_all_red <- wavelets_result_ed_tibble_biased_red_coeff_all |>
  dplyr::select(-d2_3, -d4_3, -d2_0.2, -d4_0.2)

# Create the tree plot
tree_plot <- ggtree(tree, branch.length='none') %<+% 
  merged_data +
  geom_tree() + #aes(color=class)
  scale_color_manual(values = palette_order_assigned_bloo) +
  geom_tiplab(aes(label = family_num), size=2.5, align=T) + 
  geom_tippoint(aes(color=order), alpha = 1)+
  labs(color = 'Order') +
  theme_tree2()+
  theme(legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),  # Remove grid lines
        panel.grid.minor = element_blank(),  # Remove grid lines
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 5))

# Create the first heatmap
wavelets_result_ed_tibble_biased_red_coeff_all_red <- wavelets_result_ed_tibble_biased_red_coeff_all_red |>
column_to_rownames(var = 'asv_num') 

heatmap1 <- gheatmap(tree_plot, wavelets_result_ed_tibble_biased_red_coeff_all_red, offset=8, width = 1, color=NA, font.size = 2) + 
  scale_x_discrete(labels = labs_wavelets_fract) + 
  scale_fill_gradientn(colours = palette_gradient_bw,
                       name="Bloomer's\nwavelets\ncoefficients\nmagnitude")+

  #geom_text(data ="Particle attached (3-20 um)")+
  #annotate(geom = "Particle attached (3-20 um)", x = 0.5,   hjust = 0.5) +  # Adjust x position as needed
  theme(legend.position = 'bottom', text = element_text(size = 5), # labels=c("0" = "No", "1" = "Yes")
        legend.key.size = unit(0.5, "lines"))

# Print the combined plot
print(heatmap1)

# ggsave(heatmap1, filename = 'phylogenetic_tree_wavelets_coeff_red_nosar11.pdf',
#        path = 'Results/Figures/',
#        width = 230, height = 200, units = 'mm')

#### i can not change the labs names because they are in the column names not as a variable
labs_wavelets_fract <- as_labeller(c('d1_3' = 'Fine-scale PA',
                               #'d2_3' = 'Half-yearly',
                               'd3_3' = 'Seasonal',
                               #'d4_3' = 'Year-to-year',
                               's4_3' = 'Inter-annual',
                               'd1_0.2' = 'Fine-scale',
                               #'d2_0.2' = 'Half-yearly',
                               'd3_0.2' = 'Seasonal',
                               #'d4_0.2' = 'Year-to-year',
                               's4_0.2' = 'Inter-annual'))

heatmap1 |>
    scale_x_discrete(labels = c('d1_3' = 'Fine-scale',
                               'd2_3' = 'Half-yearly',
                               'd3_3' = 'Seasonal',
                               'd4_3' = 'Year-to-year',
                               's4_3' = 'Inter-annual',
                               'd1_0.2' = 'Fine-scale',
                               'd2_0.2' = 'Half-yearly',
                               'd3_0.2' = 'Seasonal',
                               'd4_0.2' = 'Year-to-year',
                               's4_0.2' = 'Inter-annual'))

### Phyogenetic tree with variances for each coefficient------
wavelets_variance <- wavelets_result_ed_tibble_tax_3_biased_red |> 
  dplyr::mutate(fraction = '3') |>
  bind_rows( wavelets_result_ed_tibble_tax_02_biased_red) |>
  dplyr::filter(!is.na(wavelets_result_ed)) |>
  dplyr::group_by(asv_num, wavelets_transformation, class, fraction) |>
  dplyr::summarize(variance = var(wavelets_result_ed)) |>
  ungroup() |>
  dplyr::mutate(transfromation = paste0(wavelets_transformation,'_',fraction)) |>
  dplyr::select(-fraction, wavelets_transformation) |>
  pivot_wider(id_cols = asv_num, names_from = transfromation, values_from = variance, values_fill = 0  ) |>
  as_tibble() |>
  column_to_rownames(var = 'asv_num') |>
  dplyr::select("d1_3", "d2_3", "d3_3", "d4_3" ,  "s4_3", "d1_0.2" ,"d2_0.2", "d3_0.2", "d4_0.2", "s4_0.2")

wavelets_result_ed_tibble_biased_red_coeff_all |>
  colnames()

# Create the tree plot
tree_plot <- ggtree(tree, branch.length='none') %<+% 
  merged_data +
  geom_tree() + #aes(color=class)
  scale_color_manual(values = palette_order_assigned_bloo) +
  geom_tiplab(aes(label = family_num), size=2.5, align=T) + 
  geom_tippoint(aes(color=order), alpha = 1)+
  labs(color = 'Order') +
  theme_tree2()+
  theme(legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),  # Remove grid lines
        panel.grid.minor = element_blank(),  # Remove grid lines
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 5))

# Create the first heatmap
heatmap1 <- gheatmap(tree_plot, wavelets_variance, offset=8, width = 1, color=NA, font.size = 2) + 
  scale_x_discrete(labels = labs_wavelets_fract) + 
  scale_fill_gradientn(colours = palette_gradient_bw,
                       name="Bloomer's\nwavelets\ncoefficients\nvariance")+
  
  #geom_text(data ="Particle attached (3-20 um)")+
  #annotate(geom = "Particle attached (3-20 um)", x = 0.5,   hjust = 0.5) +  # Adjust x position as needed
  theme(legend.position = 'bottom', text = element_text(size = 5), # labels=c("0" = "No", "1" = "Yes")
        legend.key.size = unit(0.5, "lines"))

# Print the combined plot
print(heatmap1)

# ggsave(heatmap1, filename = 'phylogenetic_tree_wavelets_coeff_variance.pdf',
#        path = 'Results/Figures/',
#        width = 230, height = 200, units = 'mm')


## I analyse the relationship between different types of bloomers and phylogeny-----
### For this I use the phylogenetic tree and hierarchical clutsering analysis from the wavelets
library(dendextend)

# Define the taxa to be removed (replace "taxon1", "taxon2", etc. with actual taxon names)
taxa_to_remove <- bloo_02 |>
  anti_join(bloo_3) |>
  as_vector()

# Prune the phylogenetic tree to remove the specified taxa
pruned_phylo_tree <- drop.tip(tree, taxa_to_remove)

# Plot the pruned phylogenetic tree
plot(pruned_phylo_tree)

# Convert phylogenetic tree to dendrogram
# Extract branch lengths from the phylogenetic tree
branch_lengths <- cophenetic(pruned_phylo_tree)

# Create a dendrogram object using the branch lengths
dend <- as.dendrogram(hclust(as.dist(branch_lengths)))

## edit the labels that we have in the phylogenetic tree so that they match the ones we have in the hierarchical clustering 
old_label_factor <- factor(c(  "asv163" , "asv69"  , "asv80"  , "asv237" , "asv563" , "asv182" , "asv58"  , "asv126" , "asv223" , "asv116" , "asv105" , "asv249" , "asv62"  , "asv219" , "asv471" , "asv752" , "asv72"  , "asv27"  , "asv385",
                               "asv84"  , "asv43"  , "asv192" , "asv17"  , "asv77"  , "asv153" , "asv311" , "asv555" , "asv118" , "asv8"   , "asv38"  , "asv225" , "asv3"   , "asv15"  , "asv2"   , "asv264" , "asv200" , "asv5"   , "asv276",
                               "asv49"  , "asv42"  , "asv113" , "asv511" , "asv302" , "asv317" , "asv114" , "asv11"  , "asv559" , "asv282" , "asv100" , "asv22"  , "asv194" , "asv178" , "asv179" , "asv25"  , "asv23"  , "asv85"  , "asv7" , 
                               "asv1"   , "asv4"   , "asv31"  , "asv28"),
                           levels = c(
                             "asv163", "asv69", "asv80", "asv237", "asv563", "asv182", "asv58", "asv126", "asv223", "asv116",
                             "asv105", "asv249", "asv62", "asv219", "asv471", "asv752", "asv72", "asv27", "asv385", "asv84",
                             "asv43", "asv192", "asv17", "asv77", "asv153", "asv311", "asv555", "asv118", "asv8", "asv38",
                             "asv225", "asv3", "asv15", "asv2", "asv264", "asv200", "asv5", "asv276", "asv49", "asv42",
                             "asv113", "asv511", "asv302", "asv317", "asv114", "asv11", "asv559", "asv282", "asv100", "asv22",
                             "asv194", "asv178", "asv179", "asv25", "asv23", "asv85", "asv7", "asv1", "asv4", "asv31", "asv28"
                           )) 

old_label_factor <- droplevels(old_label_factor[!old_label_factor %in% taxa_to_remove])

tax_phylo_tree <- tax |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  arrange(factor(asv_num, levels = levels(old_label_factor))) |>
  dplyr::filter(!asv_num %in% taxa_to_remove)

# Create a tibble with old and new labels
new_label = tax_phylo_tree$asv_f

dend_ed_pa <- dend |>
  set("labels", new_label)

# Print the updated dendrogram
print(dend_ed_pa)

dend_ed_pa |>
  labels()

##change names to differenciate between fl and pa fraction
hc_pa_1 <- hc1
hc_pa_2 <- hc2
hc_pa_3 <- hc3
hc_pa_4 <- hc4
hc_pa_5 <- hc5

# Plot tanglegram----
dend_list <- dendlist(dend_ed_pa, hc1)
tanglegram(dend_ed_pa, hc1, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2, 
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

# Plot tanglegram
dend_list <- dendlist(dend_ed_pa, hc2)
tanglegram(dend_ed_pa, hc2, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2, 
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

# Plot tanglegram
dend_list <- dendlist(dend_ed_pa, hc3)
tanglegram(dend_ed_pa, hc3, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2, 
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

# Plot tanglegram
dend_list <- dendlist(dend_ed_pa, hc4)
tanglegram(dend_ed_pa, hc4, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2, 
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

# Plot tanglegram
dend_list <- dendlist(dend_ed_pa, hc5)
tanglegram(dend_ed_pa, hc5, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2, 
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

### The same for the FL fraction-----
# Define the taxa to be removed (replace "taxon1", "taxon2", etc. with actual taxon names)
taxa_to_remove <- bloo_3 |>
  anti_join(bloo_02) |>
  as_vector()

# Prune the phylogenetic tree to remove the specified taxa
pruned_phylo_tree <- drop.tip(tree, taxa_to_remove)

# Plot the pruned phylogenetic tree
plot(pruned_phylo_tree)

# Convert phylogenetic tree to dendrogram
# Extract branch lengths from the phylogenetic tree
branch_lengths <- cophenetic(pruned_phylo_tree)

# Create a dendrogram object using the branch lengths
dend <- as.dendrogram(hclust(as.dist(branch_lengths)))

## edit the labels that we have in the phylogenetic tree so that they match the ones we have in the hierarchical clustering 
old_label_factor <- factor(c(  "asv163" , "asv69"  , "asv80"  , "asv237" , "asv563" , "asv182" , "asv58"  , "asv126" , "asv223" , "asv116" , "asv105" , "asv249" , "asv62"  , "asv219" , "asv471" , "asv752" , "asv72"  , "asv27"  , "asv385",
                               "asv84"  , "asv43"  , "asv192" , "asv17"  , "asv77"  , "asv153" , "asv311" , "asv555" , "asv118" , "asv8"   , "asv38"  , "asv225" , "asv3"   , "asv15"  , "asv2"   , "asv264" , "asv200" , "asv5"   , "asv276",
                               "asv49"  , "asv42"  , "asv113" , "asv511" , "asv302" , "asv317" , "asv114" , "asv11"  , "asv559" , "asv282" , "asv100" , "asv22"  , "asv194" , "asv178" , "asv179" , "asv25"  , "asv23"  , "asv85"  , "asv7" , 
                               "asv1"   , "asv4"   , "asv31"  , "asv28"),
                           levels = c(
                             "asv163", "asv69", "asv80", "asv237", "asv563", "asv182", "asv58", "asv126", "asv223", "asv116",
                             "asv105", "asv249", "asv62", "asv219", "asv471", "asv752", "asv72", "asv27", "asv385", "asv84",
                             "asv43", "asv192", "asv17", "asv77", "asv153", "asv311", "asv555", "asv118", "asv8", "asv38",
                             "asv225", "asv3", "asv15", "asv2", "asv264", "asv200", "asv5", "asv276", "asv49", "asv42",
                             "asv113", "asv511", "asv302", "asv317", "asv114", "asv11", "asv559", "asv282", "asv100", "asv22",
                             "asv194", "asv178", "asv179", "asv25", "asv23", "asv85", "asv7", "asv1", "asv4", "asv31", "asv28"
                           )) 

old_label_factor <- droplevels(old_label_factor[!old_label_factor %in% taxa_to_remove])

tax_phylo_tree <- tax |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  arrange(factor(asv_num, levels = levels(old_label_factor))) |>
  dplyr::filter(!asv_num %in% taxa_to_remove)

# Create a tibble with old and new labels
label_changes <- tibble(
  old_label = old_label_factor,  # Replace with actual old labels
  new_label = tax_phylo_tree$asv_f   # Replace with corresponding new labels
)

dend_ed_fl <- dend |>
  set("labels", label_changes$new_label)

# Print the updated dendrogram
hc1 |>
  plot()

print(dend_ed)

dend_ed |>
  labels()

hc1 |>
  labels()

##change names to differenciate between fl and pa fraction
hc_fl_1 <- hc1
hc_fl_2 <- hc2
hc_fl_3 <- hc3
hc_fl_4 <- hc4
hc_fl_5 <- hc5

# Plot tanglegram
dend_list <- dendlist(dend_ed_fl, hc_fl_1)
tanglegram(dend_ed_fl, hc_fl_1, 
           common_subtrees_color_lines = TRUE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2,
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

# Plot tanglegram
dend_list <- dendlist(dend_ed_fl, hc_fl_2)
tanglegram(dend_ed_fl, hc_fl_2, 
           common_subtrees_color_lines = TRUE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2,
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

# Plot tanglegram
dend_list <- dendlist(dend_ed_fl, hc_fl_3)

tanglegram(dend_ed_fl, hc_fl_3, 
           common_subtrees_color_lines =  TRUE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2,
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

# Plot tanglegram
dend_list <- dendlist(dend_ed_fl, hc_fl_4)
tanglegram(dend_ed_fl, hc_fl_4, 
           common_subtrees_color_lines = TRUE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2,
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

# Plot tanglegram
dend_list <- dendlist(dend_ed_fl, hc_fl_5)
tanglegram(dend_ed_fl, hc_fl_5, 
           common_subtrees_color_lines = TRUE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner = 7,
           lwd = 2, 
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

##Phylogenetic tree with the category of each bloomer-----
bloo_all_types_summary <- read.csv( 'results/tables/bloo_all_types_summary_tb_tax.csv')

## I prepare the table with the category of each blooming event
bloo_all_types_summary_ed_02 <- bloo_all_types_summary |>
  dplyr::select(-X.1, -X) |>
  dplyr::filter(fraction == '0.2') 

#### i can not change the labs names because they are in the column names not as a variable
labs_wavelets_fract <- as_labeller(c('d1_3' = 'Fine-scale PA',
                                     #'d2_3' = 'Half-yearly',
                                     'd3_3' = 'Seasonal',
                                     #'d4_3' = 'Year-to-year',
                                     's4_3' = 'Inter-annual',
                                     'd1_0.2' = 'Fine-scale',
                                     #'d2_0.2' = 'Half-yearly',
                                     'd3_0.2' = 'Seasonal',
                                     #'d4_0.2' = 'Year-to-year',
                                     's4_0.2' = 'Inter-annual'))

heatmap1 |>
  scale_x_discrete(labels = c('d1_3' = 'Fine-scale',
                              'd2_3' = 'Half-yearly',
                              'd3_3' = 'Seasonal',
                              'd4_3' = 'Year-to-year',
                              's4_3' = 'Inter-annual',
                              'd1_0.2' = 'Fine-scale',
                              'd2_0.2' = 'Half-yearly',
                              'd3_0.2' = 'Seasonal',
                              'd4_0.2' = 'Year-to-year',
                              's4_0.2' = 'Inter-annual'))

# Define the taxa to be removed (replace "taxon1", "taxon2", etc. with actual taxon names)
taxa_to_remove <- bloo_3 |>
  anti_join(bloo_02) |>
  as_vector()

# Prune the phylogenetic tree to remove the specified taxa
pruned_phylo_tree <- drop.tip(tree, taxa_to_remove)

# Plot the pruned phylogenetic tree
plot(pruned_phylo_tree)

# Convert phylogenetic tree to dendrogram
# Extract branch lengths from the phylogenetic tree
branch_lengths <- cophenetic(pruned_phylo_tree)

# Create a dendrogram object using the branch lengths
dend <- as.dendrogram(hclust(as.dist(branch_lengths)))

## edit the labels that we have in the phylogenetic tree so that they match the ones we have in the hierarchical clustering 
old_label_factor <- factor(c(  "asv163" , "asv69"  , "asv80"  , "asv237" , "asv563" , "asv182" , "asv58"  , "asv126" , "asv223" , "asv116" , "asv105" , "asv249" , "asv62"  , "asv219" , "asv471" , "asv752" , "asv72"  , "asv27"  , "asv385",
                               "asv84"  , "asv43"  , "asv192" , "asv17"  , "asv77"  , "asv153" , "asv311" , "asv555" , "asv118" , "asv8"   , "asv38"  , "asv225" , "asv3"   , "asv15"  , "asv2"   , "asv264" , "asv200" , "asv5"   , "asv276",
                               "asv49"  , "asv42"  , "asv113" , "asv511" , "asv302" , "asv317" , "asv114" , "asv11"  , "asv559" , "asv282" , "asv100" , "asv22"  , "asv194" , "asv178" , "asv179" , "asv25"  , "asv23"  , "asv85"  , "asv7" , 
                               "asv1"   , "asv4"   , "asv31"  , "asv28"),
                           levels = c(
                             "asv163", "asv69", "asv80", "asv237", "asv563", "asv182", "asv58", "asv126", "asv223", "asv116",
                             "asv105", "asv249", "asv62", "asv219", "asv471", "asv752", "asv72", "asv27", "asv385", "asv84",
                             "asv43", "asv192", "asv17", "asv77", "asv153", "asv311", "asv555", "asv118", "asv8", "asv38",
                             "asv225", "asv3", "asv15", "asv2", "asv264", "asv200", "asv5", "asv276", "asv49", "asv42",
                             "asv113", "asv511", "asv302", "asv317", "asv114", "asv11", "asv559", "asv282", "asv100", "asv22",
                             "asv194", "asv178", "asv179", "asv25", "asv23", "asv85", "asv7", "asv1", "asv4", "asv31", "asv28"
                           )) 

old_label_factor <- droplevels(old_label_factor[!old_label_factor %in% taxa_to_remove])

tax_phylo_tree <- tax |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  arrange(factor(asv_num, levels = levels(old_label_factor))) |>
  dplyr::filter(!asv_num %in% taxa_to_remove)

# Create a tibble with old and new labels
label_changes <- tibble(
  old_label = old_label_factor,  # Replace with actual old labels
  new_label = tax_phylo_tree$asv_f   # Replace with corresponding new labels
)

dend_ed_fl <- dend |>
  set("labels", label_changes$new_label)

# Print the updated dendrogram
hc1 |>
  plot()

print(dend_ed)

dend_ed |>
  labels()

hc1 |>
  labels()

heatmap(bloo_all_types_summary_ed_02)

bloo_all_types_summary_ed_02 |>
  dplyr::select(asv_num, cluster_fr, fraction) 

# ggsave(heatmap1, filename = 'phylogenetic_tree_wavelets_coeff_variance.pdf',
#        path = 'Results/Figures/',
#        width = 230, height = 200, units = 'mm')

############ TREE FOR ALL MY TAXA ######## -----------------
## I need to create a fasta file with the sequences of all the taxa ----

# Aligment
##https://www.genome.jp/tools-bin/clustalw

## add fasta seqs, DNA format, and Select Weight Matrix: CLUSTALW (for DNA)
# 
# asv_num_seq_bbmo <- bbmo_10y@tax_table |>
#   as_tibble() |>
#   dplyr::select(asv_num = .otu, seq ) |>
#   left_join(tax_bbmo_10y_new, by = c('asv_num', 'seq')) |>
#   #dplyr::select(seq, asv_num) |>
#   dplyr::mutate(seq = as.character(seq),
#                 asv_num = as.character(asv_num))
# 
# # Assuming your tibble is named 'sequences_tibble' with columns 'sequence_name' and 'sequence'
# 
# # Open a connection to write the output to a file
# output_file <- "sequences_bbmo.fasta"
# file_conn <- file(output_file, "w")
# 
# # Iterate over each row in the tibble and write to the file in FASTA format
# for (i in 1:nrow(asv_num_seq_bbmo)) {
#   # Extract sequence name and sequence
#   seq_name <- as.character(asv_num_seq_bbmo$asv_num[i])  # Extract character string from the tibble
#   seq <- as.character(asv_num_seq_bbmo$seq[i])  # Ensure seq is a character vector
#   
#   # Write to file in FASTA format
#   cat(">", seq_name, "\n", seq, "\n", file = file_conn, sep = "")
# }
# 
# # Close the file connection
# close(file_conn)
# 
# # Print a message indicating the file has been written
# cat("FASTA file '", output_file, "' has been created.\n")

## Is being a bloomer a phylogenetic related trait? ----

## we compute a phylogenetic tree with all taxa present in the BBMO during those 10-Y

### upload data 
tree_complete <- ggtree::read.tree('data/raxml/complete_tree_cesga/bbmo.raxml.support')

## plot the tree
tree_complete |>
  ggplot() + 
  geom_tree(layout = 'circular')+
  theme_tree2() +
  # #geom_treescale()+
  geom_tiplab(align = T)
#geom_tippoint()

ggtree(tree_complete, layout="circular")

## I plot a tree with the following information 
### highlight bloomers in PA or FL differently
### which type of bloomer are they? recurrent or not, in the community broad, narrow or intermediate
### recurrency pattern: seasonal - chaotic

#prepare the information
tax <- tax_bbmo_10y_new |>
  dplyr::select(asv_num, phylum, class, order, family, genus) |>
  distinct()

metadata <- bloo_02 |>
  dplyr::mutate(fraction = '0.2') |>
  bind_rows(bloo_3) |>
  dplyr::mutate(fraction = case_when(is.na(fraction)~ '3',
                                     !is.na(fraction) ~ fraction)) |>
  group_by(value) |>
  dplyr::mutate(detect = n()) |>
  dplyr::mutate(detect_both = case_when(detect == '2' ~ 'both',
                                        detect == '1' ~ fraction)) |>
  group_by( value) |>
  dplyr::select(detect_both, value) |>
  distinct() |>
  rename(asv_num = value) |>
  left_join(tax)

metadata <- tax_bbmo_10y_new

merged_data <- merge(as_tibble_col(tree_complete$tip.label, column_name = 'asv_num'), metadata, by = "asv_num") |>
  dplyr::mutate(family_num = paste0(family,' ', asv_num)) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                !asv_num %in% bloo_taxonomy$asv_num ~ 'no_bloomer'))

 ## prepare the data if they are blooming on the FL or the PA fraction
 detect_both <- merged_data |>
   dplyr::select(detect_both, asv_num) |>
   dplyr::mutate(blooming_fraction = 1,
                 detect_both = str_replace(detect_both, '3', 'PA')) |>
   dplyr::mutate(detect_both = str_replace(detect_both, '0.2', 'FL')) |>
   pivot_wider(id_cols = asv_num, values_fill = 0, names_from = detect_both, values_from = blooming_fraction)

 detect_both <- detect_both |>
   dplyr::mutate_all(as.character) |>
   dplyr::mutate(PA = case_when(both == '1' ~ '1',
                                both == '0' ~ PA),
                 FL =  case_when(both == '1' ~ '1',
                                 both == '0' ~ FL)) |>
   dplyr::select(-both)   |>
   column_to_rownames(var = 'asv_num')

 palette_fraction_bw <- c(  '0' = 'white', '1' = 'black' )

 # extract the clade label information. Because some nodes of tree are
 # annotated to genera, which can be displayed with high light using ggtree.
 nodeids <- nodeid(tree_complete, tree_complete$node.label[nchar(tree_complete$node.label)>4])
 nodedf <- data.frame(node=nodeids)
 nodelab <- gsub("[\\.0-9]", "", tree_complete$node.label[nchar(tree_complete$node.label)>4])
 # The layers of clade and hightlight
 poslist <- c(1.6, 1.4, 1.6, 0.8, 0.1, 0.25, 1.6, 1.6, 1.2, 0.4,
              1.2, 1.8, 0.3, 0.8, 0.4, 0.3, 0.4, 0.4, 0.4, 0.6,
              0.3, 0.4, 0.3)
 labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)
 
 palette_order_assigned_bloo <-  c("SAR11 clade" =      "#ca6094",
                                   "Rhodospirillales"  = '#FFA180', 
                                   "Sphingomonadales"  = '#8C000A', 
                                   "Puniceispirillales" = '#4cb50f',
                                   'Rhizobiales' = '#B31722',    
                                   "Rhodobacterales" = '#2d373b',
                                   "Verrucomicrobiales"= '#005c69',
                                   "Opitutales"   =   '#74B9C8',
                                   "Phycisphaerales"  = '#e3a6ce', 
                                   "Flavobacteriales"   =  '#0051BF', 
                                   'Chitinophagales' = '#92ABFF', 
                                   "Synechococcales"  = '#009F6A', 
                                   "Bacteriovoracales" = '#8C789D',
                                   "Pseudomonadales"  = '#FF8E00', 
                                   "Enterobacterales" = "#f1c510",
                                   'Thiotrichales' =  "#000000")
 
 palette_bloomer <- c('bloomer' =  '#e3a6ce',
                      'no_bloomer' = '#ffffff')
 
 ggtree(tree_complete, layout="fan", size=0.05, open.angle=1) %<+% 
   merged_data +
   geom_tree(aes(color = order)) + #aes(color=class)
   scale_color_manual(values = palette_order_assigned_bloo) +
   #geom_tiplab(aes(), size=2.5, align=T) + 
   #geom_tippoint(aes(color=bloomer), alpha = 1, size = 1 )+
   labs(color = 'Order') +
   geom_fruit(
     geom = geom_tile,
     aes(fill = bloomer),
     color = NA,
     offset = 0.01,
     width = 1
   ) +
   scale_fill_manual(values = palette_bloomer) +
   #theme_tree2()+
   theme(
     legend.position = "right",
     panel.grid.major = element_blank(),  # Remove major grid lines
     panel.grid.minor = element_blank(),  # Remove minor grid lines
     axis.text.x = element_blank(),
     axis.text.y = element_blank(),
     axis.line.x = element_blank(),
     axis.line.y = element_blank(),
     axis.ticks = element_blank(),
     plot.background = element_blank(),
     panel.background = element_blank(),
     legend.text = element_text(size = 5)
   ) +
   scale_alpha_continuous(range = c(0, 0.7), guide = 'none')
 
 # Assume merged_data has a column named 'order' which contains the order information
 unique(bloo_taxonomy$order_f) 
 
 desired_order <- "Enterobacterales"  # Replace with the desired order
 desired_order <- "Thiotrichales"  # Replace with the desired order
 
 # Subset the data to include only the desired order
 subset_data <- merged_data[merged_data$order == desired_order, ] |>
   dplyr::filter(!is.na(order))
 
 # Subset the tree based on the tips that belong to the desired order
 subset_tree <- keep.tip(tree_complete, subset_data$asv_num)
 
 # Plot the subsetted tree
 ggtree(subset_tree, layout="fan", size=0.05, open.angle=1) %<+%
   subset_data +
   geom_tree() +  # Color branches by 'bloomer'
   scale_color_manual(values = palette_bloomer) +
   geom_fruit(
     geom = geom_tile,
     aes(fill = bloomer),
     color = NA,
     offset = 0.03,
     width = 0.05
   ) +
   scale_fill_manual(values = palette_bloomer) +
   labs(color = 'Order', fill = 'Bloomer Status') +
   theme(
     legend.position = "right",
     panel.grid.major = element_blank(),  # Remove major grid lines
     panel.grid.minor = element_blank(),  # Remove minor grid lines
     axis.text.x = element_blank(),
     axis.text.y = element_blank(),
     axis.line.x = element_blank(),
     axis.line.y = element_blank(),
     axis.ticks = element_blank(),
     plot.background = element_blank(),
     panel.background = element_blank(),
     legend.text = element_text(size = 5)
   ) +
   scale_alpha_continuous(range = c(0, 0.7), guide = 'none')  # Make points transparent where family_num is NA
 
 ## i try to loop
 unique_orders <- unique(bloo_taxonomy$order_f)  
 
 palette_bloomer <- c('bloomer' =  '#e3a6ce',
                      'no_bloomer' = '#ffffff')
 
 # Loop over each unique order
 for (desired_order in unique_orders) {
   # Subset the data to include only the desired order
   subset_data <- merged_data |>
     filter(order == desired_order & !is.na(order))
   
   # Subset the tree based on the tips that belong to the desired order
   subset_tree <- keep.tip(tree_complete, subset_data$asv_num)
   
   
   # Plot the subsetted tree
   plot <-  ggtree(subset_tree, layout="fan", size=0.05, open.angle=1) %<+%
     subset_data +
     geom_tree() +  # Color branches by 'bloomer'
     scale_color_manual(values = palette_bloomer) +
     geom_fruit(
       geom = geom_tile,
       aes(fill = bloomer),
       color = NA,
       offset = 0.03,
       width = 0.05
     ) +
     scale_fill_manual(values = palette_bloomer) +
     labs(color = 'Order', fill = 'Bloomer Status') +
     theme(
       legend.position = "right",
       panel.grid.major = element_blank(),  # Remove major grid lines
       panel.grid.minor = element_blank(),  # Remove minor grid lines
       axis.text.x = element_blank(),
       axis.text.y = element_blank(),
       axis.line.x = element_blank(),
       axis.line.y = element_blank(),
       axis.ticks = element_blank(),
       plot.background = element_blank(),
       panel.background = element_blank(),
       legend.text = element_text(size = 5)
     ) +
     scale_alpha_continuous(range = c(0, 0.7), guide = 'none')
   
   # Save or print the plot
   ggsave(paste0(desired_order, "_tree_plot.png"), plot, width = 10, height = 8)
 }

## There are too much taxa in the tree I will filter to plot only those orders that have the potential to bloom in our dataset-----
 library(ggtreeExtra)
 library(ggtree)
 library(treeio)
 library(tidytree)
 library(ggstar)
 library(ggnewscale)
 
 ### highlight bloomers in PA or FL differently
 ### which type of bloomer are they? recurrent or not, in the community broad, narrow or intermediate
 ### recurrency pattern: seasonal - chaotic
 
 palette_bloomer <- c('bloomer' =  '#e3a6ce',
                      'bloomer_0.2' =  '#93397a',
                      'bloomer_3' =  "#220000",
                      'no_bloomer' = '#ffffff',
                      'no_bloomer_NA' =  '#ffffff')
 
 merged_data <- merged_data |>
   left_join(bloo_types_summary, by = 'asv_num') |>
   dplyr::mutate(bloomer_frac = paste0(bloomer, '_', fraction))
 
 # Assume merged_data has a column named 'order' which contains the order information
 unique(bloo_taxonomy$order_f) 

 # Subset the data to include only the desired order
 subset_data <- merged_data |>
  # dplyr::filter(order %in%  unique(bloo_taxonomy$order_f))
 dplyr::filter(family %in% unique(bloo_taxonomy$family_f))

 # Subset the tree based on the tips that belong to the desired order
 subset_tree <- keep.tip(tree_complete, subset_data$asv_num)
 
 # Plot the subseted tree
 tree_plot_complete <- ggtree(subset_tree, layout="fan", size=0.05, open.angle=1) %<+%
   subset_data +
   geom_tree() +  # Color branches by 'bloomer'
   #scale_color_manual(values = palette_bloomer) +
   #geom_nodepoint(aes(color = family))+
   #scale_color_manual(values = palette_family_assigned_bloo)+
   geom_fruit(
     geom = geom_tile,
     aes(fill = bloomer),
     #color = NA,
     offset = 0.01,
     width = 0.05
   ) +
   geom_fruit(
     geom = geom_tile,
     aes(fill = bloomer_frac),
     #color = NA,
     offset = 0.01,
     width = 0.05
   ) +
   scale_fill_manual(values = palette_bloomer) +
   labs(color = 'Order', fill = 'Bloomer Status') +
   theme(
     legend.position = "right",
     panel.grid.major = element_blank(),  # Remove major grid lines
     panel.grid.minor = element_blank(),  # Remove minor grid lines
     axis.text.x = element_blank(),
     axis.text.y = element_blank(),
     axis.line.x = element_blank(),
     axis.line.y = element_blank(),
     axis.ticks = element_blank(),
     plot.background = element_blank(),
     panel.background = element_blank(),
     legend.text = element_text(size = 5)
   ) +
   scale_alpha_continuous(range = c(0, 0.7), guide = 'none')  +# Make points transparent where family_num is NA
   new_scale_fill()+
   geom_fruit(
     geom = geom_tile,
     aes(fill = occurrence_category),
     #color = NA,
     offset = 0.01,
     width = 0.05
   )+
   scale_fill_manual(
     values = c("#0000FF", "#FFA500", "#FF0000", "#800000", "#006400", "#800080"),
     na.value = "white",  # Set NA values to white
     guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
   )+
   geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
   theme(legend.position=c(0.93, 0.5),
         legend.background=element_rect(fill=NA),
         legend.title=element_text(size=6.5),
         legend.text=element_text(size=4.5),
         legend.spacing.y = unit(0.02, "cm"),
   )
 
tree_plot_complete 

ggsave(paste0(tree_plot_complete, "_tree_plot.png"), tree_plot_complete, width = 188, height = 188, units = 'mm')

 merged_data |>
   colnames()
 
### Hypothesis 2: more close phylogenetic taxa in the community leads to less blooms (for each taxa) -----

#### calculate phylogenetic distance form a tree----
# Example of calculating distances between taxa in a phylogenetic tree
# Replace 'tree' with your actual phylogenetic tree object

# Calculate pairwise distances between all taxa in the tree
pairwise_dist <- cophenetic(tree_complete)

# Display the pairwise distances
print(pairwise_dist)

row_names_tb <- pairwise_dist |>
  rownames() |>
  as_tibble_col(column_name = 'asv_num_1')

phylogenetic_distances_tb <- pairwise_dist |>
  as_tibble() |>
  bind_cols(row_names_tb) |>
  #rownames_to_column(var = 'rows_asv_num') |>
  pivot_longer(cols = -c('asv_num_1'), values_to = 'phylogenetic_distance', names_to = 'asv_num_2')

phylogenetic_distances_tb |>
  dplyr::reframe(mean = mean(phylogenetic_distance))

phylogenetic_distances_tb |>
  dplyr::filter(asv_num_1 == 'asv43') |>
  dplyr::filter(phylogenetic_distnace < 0.1) |>
  left_join(bloo_taxonomy, by = c('asv_num_2' = 'asv_num_f'))

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% c('asv43', 'asv192')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  ggplot(aes(date, abundance_value))+
  #geom_point()+
  geom_line(aes(group = asv_num, color = asv_num_f))+
  facet_wrap(vars(fraction), scales = 'free_y')+
  theme_bw()

1-0.988 ## less than 5 nulceotides of difference between those ASVs

asv_num_close <- phylogenetic_distances_tb |>
  dplyr::filter(asv_num_1 == 'asv1') |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  left_join(bloo_taxonomy, by = c('asv_num_2' = 'asv_num')) |>
  dplyr::select(asv_num_2)

closely_phylogenetically_related_example1 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% asv_num_close$asv_num_2) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  ggplot(aes(date, abundance_value))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(group = interaction(asv_num, family), color =interaction(asv_num, family)))+
  scale_color_manual(values = c("#2d373b", "#4cb76a"))+
  scale_y_continuous(labels = percent_format())+
  labs(y = 'Relative abundance (%)', color = 'ASV', x = 'Time')+
  facet_wrap(vars(fraction), scales = 'free_y', labeller = labs_fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

closely_phylogenetically_related_example1 

# ggsave(closely_phylogenetically_related_example1, filename = 'closely_phylogenetically_related_example1.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 80, units = 'mm')

closely_phylogenetically_related_example1 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% asv_num_close$asv_num_2) |>
  dplyr::filter(abundance_type == 'rclr') |>
  ggplot(aes(date, abundance_value))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(group = interaction(asv_num, family), color =interaction(asv_num, family)))+
  scale_color_manual(values = c("#2d373b", "#4cb76a"))+
  labs(y = 'rCLR', color = 'ASV', x = 'Time')+
  facet_wrap(vars(fraction), scales = 'free_y', labeller = labs_fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

closely_phylogenetically_related_example1 

# ggsave(closely_phylogenetically_related_example1, filename = 'closely_phylogenetically_related_example1_rclr.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 80, units = 'mm')

## look for those taxa that are closely phylogenetically related----
### new palette
palete_close_asvs <- c("asv1" =     '#2D373B' ,
                       'asv7' = '#009F6A',
                       'asv4' = "#686868",
                       'asv31' =  "#ae60ca",
                       'asv17'= '#9C0000',
                       'asv77' = '#74B9C8',
                       'asv200' = "#ca6094",
                       'asv264' = '#2D373B',
                       'asv42' = '#FFA200',
                       'asv49' = '#92ABFF')

labs_family_fraction <- as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                                                         '3' = 'Particle attached (3-20 um)',
                                      'asv4.asv31'='Cyanobiaceae',
                                      'asv7.asv1' = 'Cyanobiaceae',
                                      'asv264.asv200' = 'Clade I',
                                      'asv49.asv42' = 'SAR86 clade',
                                      'asv17.asv77' = 'Sphingomonadaceae'))
                          
close_bloomers <- phylogenetic_distances_tb |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  dplyr::filter(!asv_num_1 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(!asv_num_2 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  # dplyr::filter(asv_num_1 %in% bloo_all_types_summary_tb_tax$asv_num &
  #                 asv_num_2 %in% bloo_all_types_summary_tb_tax$asv_num) |>
  dplyr::filter(asv_num_1 != asv_num_2) |>
  dplyr::select(asv_num_2)
  
close_groups <- phylogenetic_distances_tb |>
  dplyr::filter(!asv_num_1 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(!asv_num_2 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  dplyr::filter(asv_num_1 %in% bloo_all_types_summary_tb_tax$asv_num &
                  asv_num_2 %in% bloo_all_types_summary_tb_tax$asv_num) |>
  dplyr::filter(asv_num_1 != asv_num_2) |>
  dplyr::mutate(close_group = paste0(asv_num_1, '.', asv_num_2)) |>
  dplyr::filter(asv_num_1 %in% c('asv17', 'asv264', 'asv49', 'asv7', 'asv4')) |>
  pivot_longer(cols = starts_with('asv_num'))

closely_phylogenetically_related_example2 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% close_bloomers$asv_num_2) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::left_join(close_groups, by = c('asv_num' = 'value')) |>
  group_by(close_group, date) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, abundance_value))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(group = asv_num, color = asv_num))+
  #geom_area(aes(fill = asv_num, color = asv_num), position = 'stack', alpha = 0.2)+
  scale_fill_manual( values = palete_close_asvs)+
  scale_color_manual( values = palete_close_asvs)+
  scale_y_continuous(labels = percent_format())+
  labs(y = 'Relative abundance (%)', fill = 'ASV', x = 'Time', color = 'ASV')+
  facet_grid(close_group~fraction, scales = 'free_y', labeller = labs_family_fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

closely_phylogenetically_related_example2

# ggsave(closely_phylogenetically_related_example2, filename = 'closely_phylogenetically_related_example2.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 220, units = 'mm')

closely_phylogenetically_related_example2_data <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% close_bloomers$asv_num_2) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::left_join(close_groups, by = c('asv_num' = 'value'))
  
closely_phylogenetically_related_example2_data$asv_num <- factor(closely_phylogenetically_related_example2_data$asv_num,
                                                                 levels = c('asv17', 'asv77',
                                                                            'asv264', 'asv200',
                                                                            'asv4', 'asv31',
                                                                            'asv42', 'asv49',
                                                                            'asv1', 'asv7'))
  
closely_phylogenetically_related_example2 <-  
  closely_phylogenetically_related_example2_data |>
  ggplot(aes(date, abundance_value))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(color = asv_num), alpha = 0.5, linewidth = 0.2)+
  geom_smooth(method = 'loess', aes(color = asv_num), span = 0.086)+
  #geom_area(aes(fill = asv_num, color = asv_num), position = 'stack', alpha = 0.2)+
  scale_fill_manual( values = palete_close_asvs)+
  scale_color_manual( values = palete_close_asvs)+
  #scale_y_continuous(labels = percent_format())+
  labs(y = 'rCLR', fill = 'ASV', x = 'Time', color = 'ASV')+
  facet_grid(close_group~fraction, scales = 'free_y', labeller = labs_family_fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(4, 'mm'))+
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.05)),
    color = guide_legend(override.aes = list(alpha = 0.05))
  )

closely_phylogenetically_related_example2

# ggsave(closely_phylogenetically_related_example2, filename = 'closely_phylogenetically_related_example2_rclr_ed.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 220, units = 'mm')

#### I filter for the distances based on the phylogenetic tree from the complete dataset ------
close_bloomers <- phylogenetic_distances_tb_com |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  dplyr::filter(!asv_num_1 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(!asv_num_2 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(asv_num_1 %in% bloo_all_types_summary_tb_tax$asv_num &
                  asv_num_2 %in% bloo_all_types_summary_tb_tax$asv_num) |>
  dplyr::filter(asv_num_1 != asv_num_2) |>
  dplyr::select(asv_num_2)

close_groups <- phylogenetic_distances_tb_com |>
  dplyr::filter(!asv_num_1 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(!asv_num_2 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  dplyr::filter(asv_num_1 %in% bloo_all_types_summary_tb_tax$asv_num &
                  asv_num_2 %in% bloo_all_types_summary_tb_tax$asv_num) |>
  dplyr::filter(asv_num_1 != asv_num_2) |>
  dplyr::mutate(close_group = paste0(asv_num_1, '.', asv_num_2)) |>
  dplyr::filter(asv_num_1 %in% c('asv17', 'asv264', 'asv49', 'asv7', 'asv4')) |>
  pivot_longer(cols = starts_with('asv_num'))

closely_phylogenetically_related_example2 <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% close_bloomers$asv_num_2) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::left_join(close_groups, by = c('asv_num' = 'value')) |>
  group_by(close_group, date) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, abundance_value))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(group = asv_num, color = asv_num))+
  #geom_area(aes(fill = asv_num, color = asv_num), position = 'stack', alpha = 0.2)+
  scale_fill_manual( values = palete_close_asvs)+
  scale_color_manual( values = palete_close_asvs)+
  scale_y_continuous(labels = percent_format())+
  labs(y = 'Relative abundance (%)', fill = 'ASV', x = 'Time', color = 'ASV')+
  facet_grid(close_group~fraction, scales = 'free_y', labeller = labs_family_fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

closely_phylogenetically_related_example2

# ggsave(closely_phylogenetically_related_example2, filename = 'closely_phylogenetically_related_example2.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 220, units = 'mm')

closely_phylogenetically_related_example2_data <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num_f %in% close_bloomers$asv_num_2) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::left_join(close_groups, by = c('asv_num' = 'value'))

closely_phylogenetically_related_example2_data$asv_num <- factor(closely_phylogenetically_related_example2_data$asv_num,
                                                                 levels = c('asv17', 'asv77',
                                                                            #'asv264', 'asv200',
                                                                            'asv4', 'asv31'
                                                                            #'asv42', 'asv49',
                                                                            #'asv1', 'asv7'
                                                                            ))

closely_phylogenetically_related_example2 <-  
  closely_phylogenetically_related_example2_data |>
  ggplot(aes(date, abundance_value))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(color = asv_num), alpha = 0.5, linewidth = 0.2)+
  geom_smooth(method = 'loess', aes(color = asv_num), span = 0.086)+
  #geom_area(aes(fill = asv_num, color = asv_num), position = 'stack', alpha = 0.2)+
  scale_fill_manual( values = palete_close_asvs)+
  scale_color_manual( values = palete_close_asvs)+
  #scale_y_continuous(labels = percent_format())+
  labs(y = 'rCLR', fill = 'ASV', x = 'Time', color = 'ASV')+
  facet_grid(close_group~fraction, scales = 'free_y', labeller = labs_family_fraction)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(4, 'mm'))+
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.05)),
    color = guide_legend(override.aes = list(alpha = 0.05))
  )

closely_phylogenetically_related_example2

# ggsave(closely_phylogenetically_related_example2, filename = 'closely_phylogenetically_related_example3_rclr_ed.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 220, units = 'mm')

# i do it for all taxa (not only those that are potential bloomers): -------
tree_complete

# Calculate pairwise distances between all taxa in the tree
pairwise_dist_com <- cophenetic(tree_complete)

row_names_tb <- pairwise_dist_com |>
  rownames() |>
  as_tibble_col(column_name = 'asv_num_1')

phylogenetic_distances_tb_com <- pairwise_dist_com |>
  as_tibble() |>
  bind_cols(row_names_tb) |>
  #rownames_to_column(var = 'rows_asv_num') |>
  pivot_longer(cols = -c('asv_num_1'), values_to = 'phylogenetic_distance', names_to = 'asv_num_2')

1-0.988 ## less than 5 nulceotides of difference between those ASVs

# asv_num_close <- phylogenetic_distances_tb_com |>
#   dplyr::filter(phylogenetic_distance < 0.012) |>
#   left_join(bloo_taxonomy, by = c('asv_num_2' = 'asv_num_f')) |>
#   dplyr::select(asv_num_2)

close_bloomers_com <- phylogenetic_distances_tb_com |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  # dplyr::filter(!asv_num_1 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  # dplyr::filter(!asv_num_2 %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(asv_num_1 != asv_num_2)

close_bloomers_com_f <- close_bloomers_com |>
  dplyr::filter(asv_num_1 %in% bloo_02$value |
                  asv_num_1 %in% bloo_3$value |
                  asv_num_2 %in% bloo_02$value |
                  asv_num_2 %in% bloo_3$value)

closely_phylogenetically_related_complete_data <- asv_tab_10y_l_rel |>
  dplyr::filter(asv_num %in% close_bloomers_com_f$asv_num_1) |>
  left_join(m_bbmo_10y) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num')
  
close_bloomers_com_f |>
  dplyr::filter(asv_num_1 == 'asv1')
  
closely_phylogenetically_related_complete_data |>  
    dplyr::filter(asv_num %in% c('asv1', 'asv3073', 'asv343', 'asv26831', 'asv5224', 'asv5327', 'asv4780', 'asv5283', 'asv7590')) |>
  ggplot(aes(date, relative_abundance))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(group = interaction(asv_num, family), color =interaction(asv_num, family)))+
  #scale_color_manual(values = c("#2d373b", "#4cb76a"))+
  scale_y_continuous(labels = percent_format())+
  labs(y = 'Relative abundance (%)', color = 'ASV', x = 'Time')+
  facet_wrap(vars(fraction), scales = 'free_y')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

## I plot rclr values 
asv_tab_10y_02_rclr <- rclr_df |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'rclr') |>
  dplyr::filter(str_detect(sample_id, '_0.2_')) 

asv_tab_10y_3_rclr <- rclr_df |>
  rownames_to_column(var = 'sample_id') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'rclr') |>
  dplyr::filter(str_detect(sample_id, '_3_'))

asv_tab_10y_rclr <- asv_tab_10y_02_rclr |>
  bind_rows(asv_tab_10y_3_rclr)

asv_tab_10y_rclr_f <- asv_tab_10y_rclr |>
  dplyr::filter(asv_num %in% close_bloomers_com_f$asv_num_1)

closely_phylogenetically_related_complete_data <- asv_tab_10y_rclr_f  |>
  left_join(m_bbmo_10y) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num')

close_bloomers_com_f |>
  dplyr::filter(asv_num_1 == 'asv1')

closely_phylogenetically_related_complete_data |>  
  dplyr::filter(asv_num %in% c('asv1', 'asv3073', 'asv343', 'asv26831', 'asv5224', 'asv5327', 'asv4780', 'asv5283', 'asv7590')) |>
  ggplot(aes(date, rclr))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(group = interaction(asv_num, family), color =interaction(asv_num, family)))+
  #scale_color_manual(values = c("#2d373b", "#4cb76a"))+
  #scale_y_continuous(labels = percent_format())+
  labs(y = 'rCLR', color = 'ASV', x = 'Time')+
  facet_wrap(vars(fraction), scales = 'free_y')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

filter_tb <- close_bloomers_com_f |>
  dplyr::filter(asv_num_1 == 'asv17')

closely_phylogenetically_related_complete_data |>  
  dplyr::filter(asv_num %in% filter_tb$asv_num_1 |
                  asv_num %in% filter_tb$asv_num_2) |>
  ggplot(aes(date, rclr))+
  #geom_point(aes(shape = asv_num_f))+
  geom_line(aes(group = interaction(asv_num, family), color =interaction(asv_num, family)))+
  #scale_color_manual(values = c("#2d373b", "#4cb76a"))+
  #scale_y_continuous(labels = percent_format())+
  labs(y = 'rCLR', color = 'ASV', x = 'Time')+
  facet_wrap(vars(fraction), scales = 'free_y')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

## competition over time (asv77)
closely_phylogenetically_related_complete_data |>
  colnames()

closely_phylogenetically_related_complete_data$date

closely_phylogenetically_related_complete_data |>  
  separate(date, into = c('year', 'month', 'day'), sep = '-', remove = F) |>
  dplyr::filter(asv_num %in% filter_tb$asv_num_1 |
                  asv_num %in% filter_tb$asv_num_2) |>
  pivot_wider(id_cols = c('date', 'fraction', 'year'), values_from = 'rclr', names_from = 'asv_num') |>
  pivot_longer(cols = c('asv2300', 'asv955', 'asv557', 'asv3052', 'asv2929', 'asv1185', 'asv876', 'asv77')) |>
  dplyr::group_by(date, fraction, asv17, year) |>
  dplyr::reframe(competition = sum(value)) |>
  ggplot(aes(asv17, competition))+
  geom_point(aes(color = year))+
  geom_smooth(method = 'loess')+
  scale_color_manual(values = palette_years)+
  facet_wrap(vars(fraction))

closely_phylogenetically_related_complete_data |>  
  separate(date, into = c('year', 'month', 'day'), sep = '-', remove = F) |>
  dplyr::filter(asv_num %in% filter_tb$asv_num_1 |
                  asv_num %in% filter_tb$asv_num_2) |>
  pivot_wider(id_cols = c('date', 'fraction', 'year'), values_from = 'rclr', names_from = 'asv_num') |>
  pivot_longer(cols = c('asv2300', 'asv955', 'asv557', 'asv3052', 'asv2929', 'asv1185', 'asv876', 'asv17')) |>
  dplyr::group_by(date, fraction, asv77, year) |>
  dplyr::reframe(competition = sum(value)) |>
  ggplot(aes(asv77, competition))+
  geom_point(aes(color = year))+
  geom_smooth(method = 'loess')+
  scale_color_manual(values = palette_years)+
  facet_wrap(vars(fraction))

## asv7
filter_tb <- close_bloomers_com_f |>
  dplyr::filter(asv_num_1 == 'asv7')

closely_phylogenetically_related_complete_data |>  
  separate(date, into = c('year', 'month', 'day'), sep = '-', remove = F) |>
  dplyr::filter(asv_num %in% filter_tb$asv_num_1 |
                  asv_num %in% filter_tb$asv_num_2) |>
  pivot_wider(id_cols = c('date', 'fraction', 'year'), values_from = 'rclr', names_from = 'asv_num') |>
  pivot_longer(cols = c('asv31034', 'asv109', 'asv974')) |>
  dplyr::group_by(date, fraction, asv7, year) |>
  dplyr::reframe(competition = sum(value)) |>
  ggplot(aes(competition, asv7))+
  geom_point(aes(color = year))+
  geom_smooth(method = 'loess')+
  scale_color_manual(values = palette_years)+
  facet_wrap(fraction)

closely_phylogenetically_related_complete_data |>  
  separate(date, into = c('year', 'month', 'day'), sep = '-', remove = F) |>
  dplyr::filter(asv_num %in% filter_tb$asv_num_1 |
                  asv_num %in% filter_tb$asv_num_2) |>
  pivot_wider(id_cols = c('date', 'fraction', 'year'), values_from = 'rclr', names_from = 'asv_num') |>
  pivot_longer(cols = c('asv31034', 'asv109', 'asv974')) |>
  dplyr::group_by(date, fraction, asv7, year) |>
  ggplot(aes(value, asv7))+
  geom_point(aes(color = year))+
  geom_smooth(method = 'loess', aes(group = year, color = year))+
  scale_color_manual(values = palette_years)+
  facet_wrap(fraction~name)

## asv4
filter_tb <- close_bloomers_com_f |>
  dplyr::filter(asv_num_1 == 'asv4')

closely_phylogenetically_related_complete_data |>  
  dplyr::filter(asv_num %in% filter_tb$asv_num_1 |
                  asv_num %in% filter_tb$asv_num_2) |>
  pivot_wider(id_cols = c('date', 'fraction'), values_from = 'rclr', names_from = 'asv_num') |>
  pivot_longer(cols = filter_tb$asv_num_2) |>
  dplyr::group_by(date, fraction, unique(filter_tb$asv_num_1)) |>
  dplyr::reframe(competition = sum(value)) |>
  ggplot(aes(competition, unique(filter_tb$asv_num_1)))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap(vars(fraction))

## asv1
bloo_02_filt <- bloo_02 |>
  dplyr::filter(value %in% c('asv2','asv3', 'asv5', 'asv8'))

filter_tb <- close_bloomers_com_f |>
  dplyr::mutate(competition_group = case_when(
    asv_num_1 %in% bloo_02_filt$value ~ asv_num_1,
    asv_num_1 %in% bloo_3$value ~ asv_num_1,
    TRUE ~ NA_character_ # This line ensures that all cases are handled
  )) |>
  dplyr::filter(!is.na(competition_group))

closely_phylogenetically_related_complete_data_f <- closely_phylogenetically_related_complete_data |>
  dplyr::select(asv_num, bloomer_rclr = rclr , date, fraction) 

closely_phylogenetically_related_complete_data |>  
  separate(date, into = c('year', 'month', 'day'), sep = '-', remove = F) |>
  dplyr::filter(asv_num %in% filter_tb$asv_num_2) |>
  left_join(filter_tb, by = c('asv_num' = 'asv_num_2'), relationship = "many-to-many") |>
  dplyr::select(rclr_competitor = rclr, date, fraction, asv_num, competition_group, year) |>
  # dplyr::group_by(date, fraction, competition_group, year) |>
  # dplyr::reframe(competition = sum(rclr)) |>
  left_join(closely_phylogenetically_related_complete_data_f, by = c('date', 'fraction', 'competition_group' = 'asv_num')) |>
  dplyr::filter(rclr_competitor != 0 |
                  bloomer_rclr != 0) |> # when there is NO competition we do not expect it to be affecting bloomers abundance
  group_by(competition_group) |>
  dplyr::filter(n() > 10) |>
  ggplot(aes(rclr_competitor, bloomer_rclr))+
  geom_point(aes(color = year))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_smooth(aes(group = competition_group),method = 'loess')+
  facet_wrap(fraction~competition_group, scales = 'free')+
  scale_color_manual(values = palette_years)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/11,
        panel.grid.minor  = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        panel.border = element_blank())

closely_phylogenetically_related_complete_data %$%
  unique(asv_num)

## plot the relationship between bloomers and their closely phylogenetically related taxa-----
paired_bona <- c("#002300","#919991",
                 "#00639a","#b8df83","#fb8a90","#adcecd",
                 "#31a02d","#d13e2f","#fdbf1f",
                 "#fe7f00","#cab4d6","#560000",
                 "#e3e317","#6a1d9a","#28b153","#ffd047","#a6cee3",
                 "#009e73","#fa344e","#9a57b1","#1f78b4","#5bb50a","#fb9a99",
                 "#002900","#666666","#e69f00","#009e73","#cc6666","#6f41d1","#237300","#f15e70",
                 "#252626","#3bb6cc","#e41a1c","#626591","#8d4c6a","#377eb8","#6f473a","#419681","#377eb8","#6f473a","#419681","#002300","#00639a","#919991",
                 "#b8df83","#fb8a90","#adcecd","#e31a1c","#fdbf1f",
                 "#fe7f00","#cab4d6","#c86f75","#6a1d9a","#e3e317","#b15928","#e69f00","#009e73",
                 "#cc6666","#a6cee3","#d14141","#1f78b4","#5bb50a","#fb9a99","#cab2d6","#6a3d9a",
                 "#ffff99","#b15928","#8c8888","#002900","#666666","#e69f00","#009e73","#cc6666","#6f41d1","#237300","#f15e70","#252626","#3bb6cc","#e41a1c",
                 "#626591","#8d4c6a","#377eb8","#6f473a","#419681","#377eb8","#6f473a","#419681","#002300","#00639a","#919991","#b8df83","#fb8a90","#adcecd",
                 "#31a02d","#e31a1c","#fdbf1f","#fe7f00","#cab4d6","#c86f75","#6a1d9a","#e3e317","#b15928","#e69f00","#009e73","#cc6666","#a6cee3","#d14141",
                 "#1f78b4","#5bb50a","#fb9a99","#cab2d6","#6a3d9a","#ffff99","#b15928","#8c8888","#002900","#666666","#e69f00","#009e73","#cc6666",
                 "#6f41d1","#237300","#f15e70","#252626","#3bb6cc","#e41a1c","#626591","#8d4c6a","#377eb8","#6f473a","#419681","#377eb8","#6f473a","#419681")

closely_phylogenetically_related_complete_data <- closely_phylogenetically_related_complete_data |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num))

closely_phylogenetically_related_complete_data$class_f <-  factor(closely_phylogenetically_related_complete_data$class_f, 
                                          levels=unique(closely_phylogenetically_related_complete_data$class_f[order(closely_phylogenetically_related_complete_data$phylum_f)]), 
                                          ordered=TRUE)

closely_phylogenetically_related_complete_data$order_f <-  factor(closely_phylogenetically_related_complete_data$order_f, 
                                          levels=unique(closely_phylogenetically_related_complete_data$order_f[order(closely_phylogenetically_related_complete_data$phylum_f,
                                                                                             closely_phylogenetically_related_complete_data$class_f)]), 
                                          ordered=TRUE)

closely_phylogenetically_related_complete_data$family_f <-  factor(closely_phylogenetically_related_complete_data$family_f, 
                                           levels=unique(closely_phylogenetically_related_complete_data$family_f[order(closely_phylogenetically_related_complete_data$phylum_f,
                                                                                               closely_phylogenetically_related_complete_data$class_f,
                                                                                               closely_phylogenetically_related_complete_data$order_f)]), 
                                           ordered=TRUE)

closely_phylogenetically_related_complete_data$asv_num_f <-  factor(closely_phylogenetically_related_complete_data$asv_num_f, 
                                            levels=unique(closely_phylogenetically_related_complete_data$asv_num_f[order(closely_phylogenetically_related_complete_data$phylum_f,
                                                                                                 closely_phylogenetically_related_complete_data$class_f,
                                                                                                 closely_phylogenetically_related_complete_data$order_f,
                                                                                                 closely_phylogenetically_related_complete_data$family_f)]), 
                                            ordered=TRUE)

closely_phylogenetically_related_complete_data |>  
  dplyr::filter(asv_num %in% filter_tb$asv_num_2) |>
  left_join(filter_tb, by = c('asv_num' = 'asv_num_2'), relationship = "many-to-many") |>
  dplyr::select(rclr, date, fraction, asv_num_f, competition_group, asv_num) |>
  # dplyr::group_by(date, fraction, competition_group) |>
  # dplyr::reframe(competition = sum(rclr)) |>
  left_join(closely_phylogenetically_related_complete_data_f, by = c('date', 'fraction', 'competition_group' = 'asv_num')) |>
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(rclr >=  0.1)) |>
  group_by(competition_group) |>
  #dplyr::filter(any(bloomer_rclr >=  0.1)) |>
  dplyr::left_join(occurrence_bloo_bbmo, by = c('competition_group' = 'asv_num', 'fraction')) |>
  dplyr::left_join(bloo_taxonomy, by = c('competition_group' = 'asv_num_f')) |>
  dplyr::filter(occurrence_category != 'narrow') |>
  dplyr::filter(rclr != 0 &
                  bloomer_rclr != 0) |> # when there is NO competition we do not expect it to be affecting bloomers abundance
  group_by(asv_num, competition_group, fraction) |>
  dplyr::filter(n() > 8) |>
  ungroup() |>
  # group_by(competition_group) |>
  # dplyr::filter(n() > 10) |>
  ggplot(aes(rclr, bloomer_rclr))+
  geom_point(aes(color = asv_num_f, shape = fraction))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  #geom_smooth(aes(rclr, bloomer_rclr, group = asv_num_f), method = 'loess', se = F)+
  #scale_color_viridis(discrete = T)+
  scale_color_manual(values = paired_bona)+
  scale_shape_discrete(labels = labs_fraction)+
  geom_smooth(aes(group = interaction(fraction, asv_num_f), color = asv_num_f), method = 'lm')+
  #stat_poly_line(aes(group = interaction(fraction, asv_num_f), color = asv_num_f), fm.values = TRUE)+
  stat_poly_eq(aes(group = interaction(fraction, asv_num_f), color = asv_num_f,
    label =  paste(after_stat(p.value.label))), 
    #position = 'stack',
    #formula = x ~ y,
    method = stats::lm,
    p.digits = 2, 
    coef.keep.zeros = T, #npcx = 1, 
   #npcy = 1,
    na.rm = FALSE)+
  facet_wrap(vars(interaction(competition_group, family_f)), scales = 'fixed', ncol = 3)+
  labs(color = 'ASV', x = 'rCLR closely phylogenetically related taxa with a potential bloomer',
       y = 'Potential bloomer rCLR', shape = 'Fraction')+
  theme_bw()+
  theme(strip.background = element_blank(),#element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/6,
        panel.grid.minor  = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),   
        panel.border = element_blank(),
        legend.key.size = unit(4, 'mm'))+
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.05)),
    color = guide_legend(override.aes = list(alpha = 0.05)),
    shape = guide_legend(ncol = 1)
  )

### I divide them for those that present significant correlations with their closely photogenically related ASVs (either negative or positive) -------

# Define a custom function to compute correlation and p-value
compute_correlation <- function(data) {
  cor_test <- cor.test(data$rclr, data$bloomer_rclr)
  tibble(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    significant = cor_test$p.value < 0.05
  )
}

# Data processing and correlation calculation pipeline
correlation_results <- closely_phylogenetically_related_complete_data %>%
  filter(asv_num %in% filter_tb$asv_num_2) %>%
  left_join(filter_tb, by = c('asv_num' = 'asv_num_2'), relationship = "many-to-many") %>%
  select(rclr, date, fraction, asv_num_f, competition_group, asv_num) %>%
  left_join(closely_phylogenetically_related_complete_data_f, by = c('date', 'fraction', 'competition_group' = 'asv_num')) %>%
  group_by(asv_num) %>%
  filter(any(rclr >= 0.1)) %>%
  group_by(competition_group) %>%
  left_join(occurrence_bloo_bbmo, by = c('competition_group' = 'asv_num', 'fraction')) %>%
  left_join(bloo_taxonomy, by = c('competition_group' = 'asv_num_f')) %>%
  filter(occurrence_category != 'narrow') %>%
  filter(rclr != 0 & bloomer_rclr != 0) %>%
  group_by(asv_num, competition_group, fraction) %>%
  filter(n() > 8) %>%
  ungroup() %>%
  group_by(asv_num, competition_group, fraction) %>%
  summarise(correlation_data = list(compute_correlation(pick(everything()))), .groups = 'drop') %>%
  unnest(correlation_data)

# View the results
print(correlation_results)

#filter those bloomers that have at least a significant correlation result ----
correlation_results_f <- correlation_results |>
  group_by(competition_group) |>
  dplyr::filter(any(significant) == 'TRUE') |>
  dplyr::mutate(relationship = case_when(significant == 'FALSE' ~ 'Neutral',
                                         significant == 'TRUE' & correlation > 0 ~ 'Cooperation',
                                         significant == 'TRUE' & correlation < 0 ~ 'Competition'))

correlation_results_plot <- correlation_results_f |>
  dplyr::left_join(bloo_taxonomy, by = c('competition_group' = 'asv_num_f')) |>
  group_by(competition_group,relationship, fraction) |>
  #dplyr::reframe(relationships_num = n()) |>
  ggplot(aes(interaction(family_f, competition_group)))+
  geom_bar(aes(fill = relationship))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c(Competition = "#b82a2e", Cooperation =  "#7bd596", Neutral = "#323232"))+
  facet_wrap(vars(fraction), labeller = labs_fraction)+
  theme_bw()+
  labs(fill = 'Relationship', y = 'Number',
       x = 'ASV', shape = 'Fraction')+
  theme(strip.background = element_blank(),#element_rect(fill = 'transparent'),
        legend.position = 'bottom',
       aspect.ratio = 6/7,
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor  = element_blank(), text = element_text(size = 6),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),   
        panel.border = element_blank(),
        legend.key.size = unit(4, 'mm'))+
  guides(
    fill = guide_legend(override.aes = list(alpha = 1), ncol = 3),
    color = guide_legend(override.aes = list(alpha = 1), ncol = 3),
    shape = guide_legend(ncol = 1)
  )

correlation_results_plot

5/21 ## 25% of the evaluated relationships were affected by competition
6/21  ## 14% of the evaluated relationships were cooperating
9/21 ## 67.86

# ggsave(correlation_results_plot, filename = 'rcorrelation_results_plot.pdf',
#        path = 'Results/Figures/',
#        width = 120, height = 100, units = 'mm')
  
## for those taxa that are not included in the plot (they do not have any significant correlation with their closely phylogenetically related taxa-----
correlation_results_f <- correlation_results |>
  group_by(competition_group) |>
  #dplyr::filter(!any(significant) == 'TRUE') |>
  dplyr::mutate(relationship = case_when(significant == 'FALSE' ~ 'Neutral',
                                         significant == 'TRUE' & correlation > 0 ~ 'Cooperation',
                                         significant == 'TRUE' & correlation < 0 ~ 'Competition'))

correlation_results_plot <- correlation_results_f |>
  dplyr::left_join(bloo_taxonomy, by = c('competition_group' = 'asv_num_f')) |>
  group_by(competition_group,relationship, fraction) |>
  #dplyr::reframe(relationships_num = n()) |>
  ggplot(aes(interaction(family_f, competition_group)))+
  geom_bar(aes(fill = relationship))+
  scale_fill_manual(values = c(Competition = "#b82a2e", Cooperation =  "#7bd596", Neutral = "#323232"))+
  facet_wrap(vars(fraction), labeller = labs_fraction)+
  theme_bw()+
  labs(fill = 'Relationship', y = 'Number',
       x = 'ASV', shape = 'Fraction')+
  theme(strip.background = element_blank(),#element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 7/8,
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor  = element_blank(), text = element_text(size = 5),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),   
        panel.border = element_blank(),
        legend.key.size = unit(4, 'mm'))+
  guides(
    fill = guide_legend(override.aes = list(alpha = 1), ncol = 3),
    color = guide_legend(override.aes = list(alpha = 1), ncol = 3),
    shape = guide_legend(ncol = 1)
  )

correlation_results_plot

# ggsave(correlation_results_plot, filename = 'rcorrelation_results_plot_all.pdf',
#        path = 'Results/Figures/',
#        width = 120, height = 100, units = 'mm')

correlation_results_f |>
  group_by(relationship) |>
  dplyr::reframe(n = n()) |>
  dplyr::mutate(n_perc = n/sum(n))
  
16.1+22.6+61.3

# filter the plot by significant correlations
relationship_with_close_tax <- closely_phylogenetically_related_complete_data |>
  dplyr::filter(asv_num %in% filter_tb$asv_num_2) |>
  left_join(filter_tb, by = c('asv_num' = 'asv_num_2'), relationship = "many-to-many") |>
  dplyr::select(rclr, date, fraction, asv_num_f, competition_group, asv_num) |>
  # dplyr::group_by(date, fraction, competition_group) |>
  # dplyr::reframe(competition = sum(rclr)) |>
  left_join(closely_phylogenetically_related_complete_data_f, by = c('date', 'fraction', 'competition_group' = 'asv_num')) |>
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  dplyr::filter(any(rclr >=  0.1)) |>
  group_by(competition_group) |>
  #dplyr::filter(any(bloomer_rclr >=  0.1)) |>
  dplyr::left_join(occurrence_bloo_bbmo, by = c('competition_group' = 'asv_num', 'fraction')) |>
  dplyr::left_join(bloo_taxonomy, by = c('competition_group' = 'asv_num_f')) |>
  dplyr::filter(occurrence_category != 'narrow') |>
  dplyr::filter(rclr != 0 &
                  bloomer_rclr != 0) |> # when there is NO competition we do not expect it to be affecting bloomers abundance
  group_by(asv_num, competition_group, fraction) |>
  dplyr::filter(n() > 8) |>
  ungroup() |>
  dplyr::filter(competition_group %in% correlation_results_f$competition_group ) |>
  ggplot(aes(rclr, bloomer_rclr))+
  geom_point(aes(color = asv_num_f, shape = fraction))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  #geom_smooth(aes(rclr, bloomer_rclr, group = asv_num_f), method = 'loess', se = F)+
  #scale_color_viridis(discrete = T)+
  scale_color_manual(values = paired_bona)+
  scale_shape_discrete(labels = labs_fraction)+
  geom_smooth(aes(group = interaction(fraction, asv_num_f), color = asv_num_f), method = 'lm')+
  #stat_poly_line(aes(group = interaction(fraction, asv_num_f), color = asv_num_f), fm.values = TRUE)+
  stat_poly_eq(aes(group = interaction(fraction, asv_num_f), color = asv_num_f,
                   label =  paste(after_stat(p.value.label))), 
               #position = 'stack',
               #formula = x ~ y,
               method = stats::lm,
               p.digits = 2, 
               coef.keep.zeros = T, #npcx = 1, 
               #npcy = 1,
               na.rm = FALSE)+
  facet_wrap(vars(interaction(competition_group, family_f)), scales = 'fixed', ncol = 3)+
  labs(color = 'ASV', x = 'rCLR closely phylogenetically related taxa with a potential bloomer',
       y = 'Potential bloomer rCLR', shape = 'Fraction')+
  theme_bw()+
  theme(strip.background = element_blank(),#element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 8/8,
        panel.grid.minor  = element_blank(), text = element_text(size = 6),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),   
        panel.border = element_blank(),
        legend.key.size = unit(4, 'mm'))+
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.05), ncol = 3),
    color = guide_legend(override.aes = list(alpha = 0.05), ncol = 3),
    shape = guide_legend(ncol = 1)
  )

relationship_with_close_tax

# ggsave(relationship_with_close_tax, filename = 'relationship_with_close_tax.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 180, units = 'mm')

## unfiltered----
closely_phylogenetically_related_complete_data |>  
  dplyr::filter(asv_num %in% filter_tb$asv_num_2) |>
  left_join(filter_tb, by = c('asv_num' = 'asv_num_2'), relationship = "many-to-many") |>
  dplyr::select(rclr, date, fraction, asv_num, competition_group) |>
  # dplyr::group_by(date, fraction, competition_group) |>
  # dplyr::reframe(competition = sum(rclr)) |>
  left_join(closely_phylogenetically_related_complete_data_f, by = c('date', 'fraction', 'competition_group' = 'asv_num')) |>
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  # dplyr::filter(any(rclr >=  0.1)) |>
  # group_by(competition_group) |>
  # dplyr::filter(any(bloomer_rclr >=  0.1)) |>
  dplyr::left_join(occurrence_bloo_bbmo, by = c('competition_group' = 'asv_num', 'fraction')) |>
  dplyr::left_join(bloo_taxonomy, by = c('competition_group' = 'asv_num_f')) |>
  dplyr::filter(occurrence_category == 'narrow') |>
  # dplyr::filter(rclr != 0 &
  #                 bloomer_rclr != 0) |> # when there is NO competition we do not expect it to be affecting bloomers abundance
  # group_by(competition_group) |>
  # dplyr::filter(n() > 10) |>
  ggplot(aes(rclr, bloomer_rclr))+
  geom_point(aes(color = asv_num))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  #geom_smooth(aes(rclr, bloomer_rclr, group = asv_num), method = 'loess', se = F)+
  #scale_color_viridis(discrete = T)+
  scale_color_manual(values = paired_bona)+
  geom_smooth(aes(group = asv_num, color = asv_num), method = 'lm')+
  facet_wrap(fraction~interaction(competition_group, family_f), scales = 'free')+
  labs(color = 'ASV', x = 'rCLR closely phylogenetically related taxa with a potential bloomer',
       y = 'Potential bloomer rCLR')+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        aspect.ratio = 6/6,
        panel.grid.minor  = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),   
        panel.border = element_blank(),
        legend.key.size = unit(4, 'mm'))+
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.05)),
    color = guide_legend(override.aes = list(alpha = 0.05))
  )

### for the supplementary material unfiltered ----
closely_phylogenetically_related_complete_data |>  
  separate(date, into = c('year', 'month', 'day'), sep = '-', remove = F) |>
  dplyr::filter(asv_num %in% filter_tb$asv_num_2) |>
  left_join(filter_tb, by = c('asv_num' = 'asv_num_2'), relationship = "many-to-many") |>
  dplyr::select(rclr, date, fraction, asv_num, competition_group, year) |>
  # dplyr::group_by(date, fraction, competition_group) |>
  # dplyr::reframe(competition = sum(rclr)) |>
  left_join(closely_phylogenetically_related_complete_data_f, by = c('date', 'fraction', 'competition_group' = 'asv_num')) |>
  group_by(asv_num) |>
  #dplyr::mutate(num_0 = sum(relative_abundance == 0)) |>
  #dplyr::filter(num_0 <= x) |>
  #as_tibble() |>
  # dplyr::filter(any(rclr >=  0.1)) |>
  # group_by(competition_group) |>
  # dplyr::filter(any(bloomer_rclr >=  0.1)) |>
  dplyr::left_join(occurrence_bloo_bbmo, by = c('competition_group' = 'asv_num', 'fraction')) |>
  dplyr::left_join(bloo_taxonomy, by = c('competition_group' = 'asv_num_f')) |>
 # dplyr::filter(occurrence_category != 'narrow') |>
  dplyr::filter(rclr != 0 &
                  bloomer_rclr != 0) |> # when there is NO competition we do not expect it to be affecting bloomers abundance
  group_by(asv_num, competition_group, fraction) |>
  dplyr::filter(n() > 10) |>
  ungroup() |>
  ggplot(aes(rclr, bloomer_rclr))+
  geom_point(aes(color = year, shape = fraction))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = palette_years)+
  scale_shape_discrete(labels = labs_fraction)+
  geom_smooth(aes(group = interaction(fraction, asv_num)), method = 'loess')+
  #stat_poly_line(aes(group = interaction(fraction, asv_num), color = asv_num), fm.values = TRUE)+
  # stat_poly_eq(aes(group = interaction(fraction, asv_num), color = asv_num,
  #                  label =  paste(after_stat(p.value.label))),
  #              #position = 'stack',
  #              #formula = x ~ y,
  #              method = stats::lm,
  #              p.digits = 2,
  #              coef.keep.zeros = T, #npcx = 1,
  #              #npcy = 1,
  #              na.rm = FALSE)+
  # stat_cor(aes(group = interaction(fraction, asv_num), color = 'black',
  #              label =   paste(..p.label..)), #label.x = 0.5, 
  #          #label.y = 5,
  #          p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman',
  #          position = position_jitter(0.025))+
  facet_wrap(vars(interaction(competition_group, family_f)), scales = 'fixed', ncol = 3)+
  labs(color = 'ASV', x = 'rCLR closely phylogenetically related taxa with a potential bloomer',
       y = 'Potential bloomer rCLR', shape = 'Fraction')+
  theme_bw()+
  theme(strip.background = element_blank(),#element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 6/6,
        panel.grid.minor  = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),   
        panel.border = element_blank(),
        legend.key.size = unit(4, 'mm'))+
  guides(
    fill = 'none',
    color = 'none'
  )

  ## loop for all bloooming ASVs -----
# Assuming 'filter_tb' and 'closely_phylogenetically_related_complete_data' are your data frames

library(dplyr)
library(ggplot2)
library(purrr)

# Function to create plot for a given asv_num
create_plot <- function(asv) {
  # Filter the tibble for the given asv_num
  filter_tb <- close_bloomers_com_f |>
    dplyr::filter(asv_num_1 == asv)
  
  # Filter the main data for the asv_num in filter_tb
  plot_data <- closely_phylogenetically_related_complete_data |>  
    dplyr::filter(asv_num %in% filter_tb$asv_num_1 |
                    asv_num %in% filter_tb$asv_num_2)
  
  # Create the plot
  p <- ggplot(plot_data, aes(date, rclr)) +
    geom_line(aes(group = interaction(asv_num, family), color = interaction(asv_num, family))) +
    labs(y = 'rCLR', color = 'ASV', x = 'Time') +
    facet_wrap(vars(fraction), scales = 'free_y', labeller = labs_fraction) +
    theme_bw() +
    theme(strip.background = element_rect(fill = 'transparent'),
          legend.position = 'bottom',
          aspect.ratio = 6/11,
          panel.grid = element_blank(), text = element_text(size = 8),
          strip.text = element_text(margin = margin(2, 2, 2, 2)),
          plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
          legend.key.size = unit(3, 'mm'))
  
  # Return the plot
  return(p)
}

# Get a vector of unique asv_nums from bloo_taxonomy
unique_asv_nums <- unique(shared_asvs$value)

shared_asvs <- bloo_02 |>
  left_join(bloo_3) |>
  filter(!is.na(recurrency))

# Create a list of plots
plots <- map(unique_asv_nums, create_plot)

# Save each plot with a unique filename
walk2(plots, unique_asv_nums, ~ggsave(filename = paste0("plot_", .y, ".png"), plot = .x, width = 10, height = 8))

# Phylogenetic distance and blooming events vs non-blooming events ##### ----
library(FD)

## functcomp returns the functional composition of a set of communities, as measured by the community-level weighted means of trait values (CWM; e.g. Lavorel et al. 2008)

phylogenetic_distances_tb_com

bloo_02 |>
  head()

comp <- phylogenetic_distances_tb_com |>
  dplyr::filter(asv_num_1 ==  'asv38') |> #con asv.id y phylogenetic distance
  dplyr::select(-asv_num_1)

dim(comp)

abund <- asv_tab_10y_02_rel |> # con asv.id, abund, y muestra
  dplyr::select(-total_reads) |>
  pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') |>
  ungroup() |>
  dplyr::select(-sample_id)

dim(abund)

# Ensure species labels match and are ordered
comp <- comp %>%
  arrange(species)

abund <- abund %>%
  arrange(species)

# Convert tibbles to matrices
comp_matrix <- comp |>
  dplyr::select(-asv_num_2) |>
  as.matrix()

rownames(comp_matrix) <- comp$asv_num_2

abund_matrix <- as.matrix(abund)

all(comp$asv_num_2  == colnames(abund))

weighted_traits <- FD::functcomp(x = as.matrix(comp_matrix), 
                                a = as.matrix(abund_matrix))

m_02 |>
  colnames()

weighted_traits |>
  dplyr::select(asv_38_phylogenetic_distance = phylogenetic_distance ) |>
  as_tibble() |>
  dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits))) |>
  left_join(m_02, by = 'sample_id_num') |>
  ggplot(aes(date, asv_38_phylogenetic_distance))+
  geom_line()

asv_tab_10y_02_rel |>
  dplyr::filter(asv_num == 'asv38' &
                  str_detect(sample_id, '_0.2_')) |>
  left_join(m_02) |>
  ggplot(aes(date, relative_abundance))+
  geom_line()

asv38 <- weighted_traits |>
  dplyr::select(asv_38_phylogenetic_distance = phylogenetic_distance ) |>
  as_tibble() |>
  dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits))) |>
  left_join(m_02, by = 'sample_id_num') 

asv_tab_10y_02_rel |>
  dplyr::filter(asv_num == 'asv38' &
                  str_detect(sample_id, '_0.2_')) |>
  left_join(m_02)  |>
  left_join(asv38) |>
  ggplot(aes(relative_abundance, asv_38_phylogenetic_distance))+
  geom_point()

## another try 
comp <- phylogenetic_distances_tb_com |>
  dplyr::filter(asv_num_1 ==  'asv23') |> #con asv.id y phylogenetic distance
  dplyr::select(-asv_num_1)

dim(comp)

abund <- asv_tab_10y_02_rel |> # con asv.id, abund, y muestra
  dplyr::select(-total_reads) |>
  pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') |>
  ungroup() |>
  dplyr::select(-sample_id)

dim(abund)

# Ensure all(comp$asv_num_2 == colnames(abund))
if (all(comp$asv_num_2 %in% colnames(abund))) {
  # Reorder columns in abund
  abund <- abund[, comp$asv_num_2]
} else {
  stop("Not all elements in comp$asv_num_2 match colnames in abund.")
}

# Convert tibbles to matrices
comp_matrix <- comp |>
  dplyr::select(-asv_num_2) |>
  as.matrix()

rownames(comp_matrix) <- comp$asv_num_2

abund_matrix <- as.matrix(abund)

all(comp$asv_num_2  == colnames(abund))

weighted_traits <- FD::functcomp(x = as.matrix(comp_matrix), 
                                 a = as.matrix(abund_matrix))

m_02 |>
  colnames()

weighted_traits |>
  dplyr::select(asv_38_phylogenetic_distance = phylogenetic_distance ) |>
  as_tibble() |>
  dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits))) |>
  left_join(m_02, by = 'sample_id_num') |>
  ggplot(aes(date, asv_38_phylogenetic_distance))+
  geom_line()

asv_tab_10y_02_rel |>
  dplyr::filter(asv_num == 'asv23' &
                  str_detect(sample_id, '_0.2_')) |>
  left_join(m_02) |>
  ggplot(aes(date, relative_abundance))+
  geom_line()

asv23 <- weighted_traits |>
  dplyr::select(asv_38_phylogenetic_distance = phylogenetic_distance ) |>
  as_tibble() |>
  dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits))) |>
  left_join(m_02, by = 'sample_id_num') 

asv_tab_10y_02_rel |>
  dplyr::filter(asv_num == 'asv23' &
                  str_detect(sample_id, '_0.2_')) |>
  left_join(m_02)  |>
  left_join(asv23) |>
  ggplot(aes(relative_abundance, asv_38_phylogenetic_distance))+
  geom_point()

## I loop over all the potential bloomers that i have in the FL fraction  --------
library(dplyr)
library(tidyr)
library(ggplot2)
library(FD) ## weighted traits

# Example data
# Assuming `bloo_02` is a tibble with `asv_num` column
bloo_taxonomy <- asv_tab_all_bloo_z_tax %>%
  dplyr::select(phylum, class, order, family, asv_num) |>
  unique()

bloo_02 <- bloo_02 |>
  dplyr::filter(!value %in% c('asv2', 'asv3', 'asv5', 'asv8'))

bloo_taxonomy <- bloo_taxonomy |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8'))

# # Loop over each `asv_num` in `bloo_02`
# for (asv in bloo_02$value) {
#   # Filter comp based on the current asv_num
#   comp <- phylogenetic_distances_tb_com %>%
#     dplyr::filter(asv_num_1 == asv) %>%
#     dplyr::select(-asv_num_1)
#   
#   # Prepare abund tibble
#   abund <- asv_tab_10y_02_rel %>%
#     dplyr::select(-total_reads) %>%
#     pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') %>%
#     ungroup() %>%
#     dplyr::select(-sample_id)
#   
#   # Ensure all species labels match and are ordered
#   comp <- comp %>%
#     arrange(asv_num_2)
#   
#   if (all(comp$asv_num_2 %in% colnames(abund))) {
#     # Reorder columns in abund
#     abund <- abund[, comp$asv_num_2]
#   } else {
#     message(paste("Not all elements in comp$asv_num_2 match colnames in abund for", asv))
#     next
#   }
#   
#   # Convert tibbles to matrices
#   comp_matrix <- comp %>%
#     dplyr::select(-asv_num_2) %>%
#     as.matrix()
#   
#   rownames(comp_matrix) <- comp$asv_num_2
#   
#   abund_matrix <- as.matrix(abund)
#   
#   # Ensure the column names of abund_matrix match the row names of comp_matrix
#   if (all(comp$asv_num_2 == colnames(abund_matrix))) {
#     # Call the functcomp function
#     weighted_traits <- FD::functcomp(x = comp_matrix, a = abund_matrix)
#     
#     # Rename the column
#     colname <- paste0(asv, "_phylogenetic_distance")
#     weighted_traits <- weighted_traits %>%
#       as_tibble() %>%
#       rename(!!colname := phylogenetic_distance) %>%
#       dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits)))
#     
#     # Example plot with m_02
#     p1 <- weighted_traits %>%
#       left_join(m_02, by = 'sample_id_num') %>%
#       ggplot(aes(date, !!sym(colname))) +
#       geom_line() +
#       ggtitle(paste("Phylogenetic Distance for", asv)) +
#       theme_minimal()
#     ggsave(filename = paste0(asv, "_phylogenetic_distance_plot.png"), plot = p1)
#     
#     # Example plot with relative abundance
#     p2 <- asv_tab_10y_02_rel %>%
#       filter(asv_num == asv & str_detect(sample_id, '_0.2_')) %>%
#       left_join(m_02, by = 'sample_id') %>%
#       ggplot(aes(date, relative_abundance)) +
#       geom_line() +
#       ggtitle(paste("Relative Abundance for", asv)) +
#       theme_minimal()
#     ggsave(filename = paste0(asv, "_relative_abundance_plot.png"), plot = p2)
#     
#     # Save the combined data for further analysis
#     combined_data <- asv_tab_10y_02_rel %>%
#       filter(asv_num == asv & str_detect(sample_id, '_0.2_')) %>%
#       left_join(m_02, by = 'sample_id') %>%
#       left_join(weighted_traits, by = 'sample_id_num')
#     
#     p3 <- combined_data %>%
#       ggplot(aes(relative_abundance, !!sym(colname))) +
#       geom_point() +
#       ggtitle(paste("Abundance vs Phylogenetic Distance for", asv)) +
#       theme_minimal()
#     ggsave(filename = paste0(asv, "_abundance_vs_distance_plot.png"), plot = p3)
#   } else {
#     message(paste("Species labels in 'comp' and 'abund' are not identical and ordered for", asv))
#   }
# }

##### I first calculate the weighted phylogenetic distances for each blooming ASV and all the samples ----

## I would like to add some information to the Order plot
#### where are the blooming events
asv_tab_all_bloo_z_tax <- read.delim2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv', sep = ';') 
## add blooming events or not 
bloom_event <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(bloom_event = case_when(abundance_type == 'relative_abundance' &
                                          abundance_value >= 0.1 &
                                          z_score_ra > 1.96 ~ 'bloom',
                                        TRUE ~ 'no-bloom')) |>
  filter(bloom_event != 0) |>
  dplyr::select(date, asv_num, bloom_event, fraction) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ungroup() |>
  distinct()

## to distinguish between bloom events and no-blooming events

## FL 
# Create an empty list to store combined data for each ASV

m_02_red <- m_02 |>
  dplyr::select(sample_id, date, fraction, sample_id_num) |>
  distinct()

# combined_data_list <- list()
# 
# # Loop over each `asv_num` in the datset that is defined as a potential blooming event either in one fraction or the other 
# for (asv in bloo_taxonomy$asv_num_f) {
#   # dplyr::filter comp based on the current asv_num
#   comp <- phylogenetic_distances_tb_com |>
#     dplyr::filter(asv_num_1 == asv) |>
#     dplyr::select(-asv_num_1)
#   
#   # Prepare abund tibble
#   abund <- asv_tab_10y_02_rel |>
#     dplyr::select(-total_reads) |>
#     pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') %>%
#     ungroup() %>%
#     dplyr::select(-sample_id)
#   
#   # Ensure all species labels match and are ordered
#   comp <- comp |>
#     arrange(asv_num_2)
#   
#   if (all(comp$asv_num_2 %in% colnames(abund))) {
#     # Reorder columns in abund
#     abund <- abund[, comp$asv_num_2]
#   } else {
#     message(paste("Not all elements in comp$asv_num_2 match colnames in abund for", asv))
#     next
#   }
#   
#   # Convert tibbles to matrices
#   comp_matrix <- comp |>
#     dplyr::select(-asv_num_2) |>
#     as.matrix()
#   
#   rownames(comp_matrix) <- comp$asv_num_2
#   
#   abund_matrix <- as.matrix(abund)
#   
#   # Ensure the column names of abund_matrix match the row names of comp_matrix
#   if (all(comp$asv_num_2 == colnames(abund_matrix))) {
#     # Call the functcomp function
#     weighted_traits <- FD::functcomp(x = comp_matrix, a = abund_matrix)
#     
#     # Rename the column
#     colname <- paste0(asv, "_phylogenetic_distance")
#     weighted_traits <- weighted_traits |>
#       as_tibble() |>
#       dplyr::mutate(asv_num_phylo = asv) |>
#       dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits)))
#     
#     # Save the combined data for further analysis
#     combined_data <- asv_tab_10y_02_rel |>
#       dplyr::filter(asv_num == asv & str_detect(sample_id, '_0.2_')) |>
#       left_join(m_02_red, by = 'sample_id') |>
#       left_join(weighted_traits, by = 'sample_id_num') |>
#       dplyr::mutate(asv_num = asv)
#     
#     # Append the combined data to the list
#     combined_data_list[[asv]] <- combined_data
#   } else {
#     message(paste("Species labels in 'comp' and 'abund' are not identical and ordered for", asv))
#   }
# }
# 
# # Combine the list elements into a single tibble
# combined_data_all_02 <- bind_rows(combined_data_list) |>
#   dplyr::left_join(bloom_event) |>
#   distinct() |>
#   left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'))
# 
# # Generate plots outside the loop
# combined_data_all_02 |>
#   colnames()
# 
# ## plot the relationship between relative abundance and weighted phylogenetic distance with the community -----
# #plot
# combined_data_all_02 |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family))+
#   geom_point(aes(color = family_f, shape = bloom_event, alpha = if_else(bloom_event == 'bloom', 1, 0.5)))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 16),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# combined_data_all_02 |>
#   dplyr::filter(bloom_event == 'bloom') |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family))+
#   geom_point(aes(color = family_f, shape = bloom_event))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 8),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# combined_data_all_02 |>
#   dplyr::filter(bloom_event == 'no-bloom') |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family))+
#   geom_point(aes(color = family_f, shape = bloom_event))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 8),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# ## The same for PA-------
# # Create an empty list to store combined data for each ASV
# 
# combined_data_list <- list()
# 
# m_3_red <- m_3 |>
#   dplyr::select(date, fraction, sample_id, sample_id_num)
# 
# # Loop over each `asv_num` in `bloo_3`
# for (asv in bloo_3$value) {
#   # dplyr::filter comp based on the current asv_num
#   comp <- phylogenetic_distances_tb_com %>%
#     dplyr::filter(asv_num_1 == asv) %>%
#     dplyr::select(-asv_num_1)
#   
#   # Prepare abund tibble
#   abund <- asv_tab_10y_3_rel |>
#     dplyr::select(-total_reads) |>
#     pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') %>%
#     ungroup() %>%
#     dplyr::select(-sample_id)
#   
#   # Ensure all species labels match and are ordered
#   comp <- comp %>%
#     arrange(asv_num_2)
#   
#   if (all(comp$asv_num_2 %in% colnames(abund))) {
#     # Reorder columns in abund
#     abund <- abund[, comp$asv_num_2]
#   } else {
#     message(paste("Not all elements in comp$asv_num_2 match colnames in abund for", asv))
#     next
#   }
#   
#   # Convert tibbles to matrices
#   comp_matrix <- comp %>%
#     dplyr::select(-asv_num_2) %>%
#     as.matrix()
#   
#   rownames(comp_matrix) <- comp$asv_num_2
#   
#   abund_matrix <- as.matrix(abund)
#   
#   # Ensure the column names of abund_matrix match the row names of comp_matrix
#   if (all(comp$asv_num_2 == colnames(abund_matrix))) {
#     # Call the functcomp function
#     weighted_traits <- FD::functcomp(x = comp_matrix, a = abund_matrix)
#     
#     # Rename the column
#     colname <- paste0(asv, "_phylogenetic_distance")
#     weighted_traits <- weighted_traits |>
#       as_tibble() |>
#       dplyr::mutate(asv_num_phylo = asv) |>
#       dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits)))
#     
#     # Save the combined data for further analysis
#     combined_data <- asv_tab_10y_3_rel |>
#       dplyr::filter(asv_num == asv & str_detect(sample_id, '_3_')) |>
#       left_join(m_3_red, by = 'sample_id') |>
#       left_join(weighted_traits, by = 'sample_id_num') |>
#       dplyr::mutate(asv_num = asv)
#     
#     # Append the combined data to the list
#     combined_data_list[[asv]] <- combined_data
#   } else {
#     message(paste("Species labels in 'comp' and 'abund' are not identical and ordered for", asv))
#   }
# }
# 
# # Combine the list elements into a single tibble
# combined_data_all_3 <- bind_rows(combined_data_list) |>
#   dplyr::left_join(bloom_event) |>
#   distinct() |>
#   left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'))
# 
# # Generate plots outside the loop
# combined_data_all_3 |>
#   colnames()
# 
# ## plot the relationship between relative abundance and weighted phylogenetic distance with the community -----
# #plot
# combined_data_all_3 |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family))+
#   geom_point(aes(color = family_f, shape = bloom_event, alpha = if_else(bloom_event == 'bloom', 1, 0.5)))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 16),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# combined_data_all_3 |>
#   dplyr::filter(bloom_event == 'bloom') |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family))+
#   geom_point(aes(color = family_f, shape = bloom_event))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 8),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# combined_data_all_3 |>
#   dplyr::filter(bloom_event == 'no-bloom') |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family))+
#   geom_point(aes(color = family_f, shape = bloom_event))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 8),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# ## FL and PA at the same time----
# combined_data_all_3 |>
#   bind_rows(combined_data_all_02) |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family_f))+
#   geom_point(aes(color = family_f, shape = bloom_event, alpha = if_else(bloom_event == 'bloom', 1, 0.5)))+
#   geom_smooth(method = 'loess', aes(group = asv_num, color = family_f))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 14),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# ## density plot 
occurrence_bloo_bbmo <- read.delim2('data/occurrence_bloo_bbmo.csv', sep = ',')

occurrence_bloo_bbmo_red <-  occurrence_bloo_bbmo |>
  dplyr::select(occurrence_perc, fraction, date, asv_num) |>
  dplyr::mutate(occurrence_category = ifelse(occurrence_perc > 2/3, 'broad',
                                             ifelse(occurrence_perc < 1/3, 'narrow',
                                                    'intermediate'))) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))
# 
# combined_data_all_3 |>
#   bind_rows(combined_data_all_02) |>
#   dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
#   left_join(occurrence_bloo_bbmo_red) |>
#   dplyr::filter(occurrence_perc > 2/3) |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family_f))+
#   geom_point(aes(color = family_f, shape = bloom_event, alpha = if_else(bloom_event == 'bloom', 1, 0.5)))+
#   geom_smooth(method = 'loess', aes(group = asv_num, color = family_f))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num), scales = 'free')+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 14),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# combined_data_all_3 |>
#   dplyr::group_by(sample_id) |>
#   reframe(n = n())
# 
# combined_data_all_02 |>
#   dplyr::group_by(sample_id) |>
#   reframe(n = n())
# 
# phylogenetic_weighted_tb <- combined_data_all_3 |>
#   bind_rows(combined_data_all_02) 
# 
# phylogenetic_weighted_tb <- phylogenetic_weighted_tb |>
#   dplyr::mutate(phylum_f = as_factor(phylum_f),
#                 family_f = as_factor(family_f),
#                 order_f = as_factor(order_f),
#                 class_f = as_factor(class_f),
#                 asv_num_f = as_factor(asv_num))
# 
# phylogenetic_weighted_tb$class_f <-  factor(phylogenetic_weighted_tb$class_f, 
#                                                                   levels=unique(phylogenetic_weighted_tb$class_f[order(phylogenetic_weighted_tb$phylum_f)]), 
#                                                                   ordered=TRUE)
# 
# phylogenetic_weighted_tb$order_f <-  factor(phylogenetic_weighted_tb$order_f, 
#                                                                   levels=unique(phylogenetic_weighted_tb$order_f[order(phylogenetic_weighted_tb$phylum_f,
#                                                                                                                                              phylogenetic_weighted_tb$class_f)]), 
#                                                                   ordered=TRUE)
# 
# phylogenetic_weighted_tb$family_f <-  factor(phylogenetic_weighted_tb$family_f, 
#                                                                    levels=unique(phylogenetic_weighted_tb$family_f[order(phylogenetic_weighted_tb$phylum_f,
#                                                                                                                                                phylogenetic_weighted_tb$class_f,
#                                                                                                                                                phylogenetic_weighted_tb$order_f)]), 
#                                                                    ordered=TRUE)
# 
# phylogenetic_weighted_tb$asv_num_f <-  factor(phylogenetic_weighted_tb$asv_num_f, 
#                                                                     levels=unique(phylogenetic_weighted_tb$asv_num_f[order(phylogenetic_weighted_tb$phylum_f,
#                                                                                                                                                  phylogenetic_weighted_tb$class_f,
#                                                                                                                                                  phylogenetic_weighted_tb$order_f,
#                                                                                                                                                  phylogenetic_weighted_tb$family_f)]), 
#                                                                     ordered=TRUE)
# 
# phylogenetic_weighted_tb$bloom_event <- factor(phylogenetic_weighted_tb$bloom_event)
# 
# more_blooms <- phylogenetic_weighted_tb |>
#   dplyr::group_by(bloom_event, asv_num, fraction) |>
#   dplyr::reframe(n = n()) |>
#   dplyr::filter(bloom_event == 'bloom' &
#       n > 1)
# 
# boxplot_phylogenetic_distance <- phylogenetic_weighted_tb  |>
#   dplyr::filter(asv_num %in% c(occurrence_bloo_bbmo_red |>
#                                  dplyr::filter(occurrence_perc > 2/3) %$%
#                                  asv_num)) |>
#   dplyr::filter(asv_num %in% more_blooms$asv_num) |>
#   dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
#   ggplot(aes(bloom_event, phylogenetic_distance, group = as.factor(asv_num_f), color = family_f))+
#   geom_point(aes(color = family_f, shape = fraction), position = position_jitter(width = 0.1))+ #alpha = if_else(bloom_event == 'bloom', 1, 0.5)
#   scale_x_discrete(labels = labs_bloom_events)+
#   geom_boxplot(aes(group = bloom_event, alpha = 0.8))+
#   scale_shape_discrete(labels = labs_fraction)+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   labs(x='Blooming event', y = 'Weighted phylogenetic distance', color = 'Family', shape = 'Fraction')+
#   facet_wrap(vars(interaction(asv_num_f, family_f)))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 14),
#         strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))+
#   guides(color = guide_legend(ncol = 1),
#          alpha = 'none',
#          shape = guide_legend(ncol = 1))
# 
# boxplot_phylogenetic_distance
# 
# ggsave(boxplot_phylogenetic_distance, file = 'boxplot_phylogenetic_distance_v2.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 188, units = 'mm')
# 
# labs_bloom_events <- as_labeller(c('bloom' = 'Yes',
#                                    'no-bloom' = 'No'))
# 
# phylogenetic_weighted_tb |>
#   ungroup() |>
#   distinct(asv_num)
# 
# weigth_phylogenetic_distance_blooms_plot <- phylogenetic_weighted_tb  |>
#   #dplyr::filter(occurrence_perc > 2/3) |>
#   dplyr::filter(asv_num %in% more_blooms$asv_num) |>
#   dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
#   ggplot(aes( x= as.numeric(phylogenetic_distance), y =interaction(family_f, asv_num_f) ,fill = as.factor(bloom_event)), 
#          group = as.factor(bloom_event))+
#   geom_density_ridges(alpha = 0.8)+
#   #scale_color_manual(values = palette_family_assigned_bloo)+
#   labs(y ='Density', x = 'Weighted phylogenetic distance', fill = 'Bloom event')+
#   #facet_wrap(vars(asv_num))+
#   scale_fill_manual(values = c('bloom' = '#00808F' , 'no-bloom' = '#454545'), labels = labs_bloom_events)+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 14),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'),
#           panel.border = element_blank())
# 
# weigth_phylogenetic_distance_blooms_plot

## weighted phylogenetic distances removing the ASV of interest from the analysis ----
## FL 
m_02_red <- m_02 |>
  dplyr::select(sample_id, date, fraction, sample_id_num) |>
  distinct()

# Create an empty list to store combined data for each ASV

combined_data_list <- list()

# Loop over each `asv_num` in `bloo_02`
for (asv in bloo_taxonomy$asv_num) {
  # dplyr::filter comp based on the current asv_num
  comp <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv) |>
    dplyr::filter(asv_num_2 != asv) |> # remove the ASV of interest from the analysis
    dplyr::select(-asv_num_1)
  
  # I need to transform 0 by NA so that the functcomp function does not consider them
  asv_tab_10y_02_rel_ed <- asv_tab_10y_02_rel |>
    dplyr::mutate(relative_abundance = case_when(
      relative_abundance == 0 ~ NA_real_,  # Replace 0 with NA (numeric NA)
      TRUE ~ relative_abundance            # Keep the original value otherwise
    ))
  
  # Prepare abund tibble
  abund <- asv_tab_10y_02_rel_ed |>
    dplyr::select(-total_reads) |>
    pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') %>%
    dplyr::select(-asv) |> # remove the ASV of interest from the analysis
    ungroup() %>%
    dplyr::select(-sample_id)
  
  # Ensure all species labels match and are ordered
  comp <- comp |>
    arrange(asv_num_2)
  
  if (all(comp$asv_num_2 %in% colnames(abund))) {
    # Reorder columns in abund
    abund <- abund[, comp$asv_num_2]
  } else {
    message(paste("Not all elements in comp$asv_num_2 match colnames in abund for", asv))
    next
  }
  
  # Convert tibbles to matrices
  comp_matrix <- comp |>
    dplyr::select(-asv_num_2) |>
    as.matrix()
  
  rownames(comp_matrix) <- comp$asv_num_2
  
  abund_matrix <- as.matrix(abund)
  
  # Ensure the column names of abund_matrix match the row names of comp_matrix
  if (all(comp$asv_num_2 == colnames(abund_matrix))) {
    # Call the functcomp function
    weighted_traits <- FD::functcomp(x = comp_matrix, a = abund_matrix)
    
    # Rename the column
    colname <- paste0(asv, "_phylogenetic_distance")
    weighted_traits <- weighted_traits |>
      as_tibble() |>
      dplyr::mutate(asv_num_phylo = asv) |>
      dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits)))
    
    # Save the combined data for further analysis
    combined_data <- asv_tab_10y_02_rel |>
      dplyr::filter(asv_num == asv & str_detect(sample_id, '_0.2_')) |>
      left_join(m_02_red, by = 'sample_id') |>
      left_join(weighted_traits, by = 'sample_id_num') |>
      dplyr::mutate(asv_num = asv)
    
    # Append the combined data to the list
    combined_data_list[[asv]] <- combined_data
  } else {
    message(paste("Species labels in 'comp' and 'abund' are not identical and ordered for", asv))
  }
}

# Combine the list elements into a single tibble
combined_data_all_02 <- bind_rows(combined_data_list) |>
  dplyr::left_join(bloom_event) |>
  #dplyr::mutate(fraction = '0.2') |>
  distinct() |>
  dplyr::select(asv_num, sample_id, relative_abundance, fraction, date, phylogenetic_distance, bloom_event) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num'))

# Generate plots outside the loop
combined_data_all_02 |>
  colnames()

## plot the relationship between relative abundance and weighted phylogenetic distance with the community -----
#plot
combined_data_all_02 |>
  ggplot(aes(phylogenetic_distance, relative_abundance, group = asv_num, color = family))+
  geom_point(aes(color = family, shape = bloom_event, alpha = if_else(bloom_event == 'bloom', 1, 0.5)))+
  #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
  scale_color_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(asv_num))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 16),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

combined_data_all_02 |>
  dplyr::filter(bloom_event == 'bloom') |>
  ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family_f))+
  geom_point(aes(color = family_f, shape = bloom_event))+
  #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
  scale_color_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(asv_num))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

combined_data_all_02 |>
  dplyr::filter(bloom_event == 'no-bloom') |>
  ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family_f))+
  geom_point(aes(color = family_f, shape = bloom_event))+
  #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
  scale_color_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(asv_num))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 8),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))

## The same for PA-------
# Create an empty list to store combined data for each ASV

combined_data_list <- list()

m_3_red <- m_3 |>
  dplyr::select(date, fraction, sample_id, sample_id_num)

# Loop over each `asv_num` in `bloo_3`
for (asv in bloo_taxonomy$asv_num) {
  # dplyr::filter comp based on the current asv_num
  comp <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv) |>
    dplyr::filter(asv_num_2 != asv) |> # remove the ASV of interest from the analysis
    dplyr::select(-asv_num_1)
  
  # I need to transform 0 by NA so that the functcomp function does not consider them
  asv_tab_10y_3_rel_ed <- asv_tab_10y_3_rel |>
    dplyr::mutate(relative_abundance = case_when(
      relative_abundance == 0 ~ NA_real_,  # Replace 0 with NA (numeric NA)
      TRUE ~ relative_abundance            # Keep the original value otherwise
    ))
  
  # Prepare abund tibble
  abund <- asv_tab_10y_3_rel_ed |>
    dplyr::select(-total_reads) |>
    pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') %>%
    dplyr::select(-asv) |> # remove the ASV of interest from the analysis
    ungroup() |>
    dplyr::select(-sample_id)
  
  # Ensure all species labels match and are ordered
  comp <- comp %>%
    arrange(asv_num_2)
  
  if (all(comp$asv_num_2 %in% colnames(abund))) {
    # Reorder columns in abund
    abund <- abund[, comp$asv_num_2]
  } else {
    message(paste("Not all elements in comp$asv_num_2 match colnames in abund for", asv))
    next
  }
  
  # Convert tibbles to matrices
  comp_matrix <- comp %>%
    dplyr::select(-asv_num_2) %>%
    as.matrix()
  
  rownames(comp_matrix) <- comp$asv_num_2
  
  abund_matrix <- as.matrix(abund)
  
  # Ensure the column names of abund_matrix match the row names of comp_matrix
  if (all(comp$asv_num_2 == colnames(abund_matrix))) {
    # Call the functcomp function
    weighted_traits <- FD::functcomp(x = comp_matrix, a = abund_matrix)
    
    # Rename the column
    colname <- paste0(asv, "_phylogenetic_distance")
    weighted_traits <- weighted_traits |>
      as_tibble() |>
      dplyr::mutate(asv_num_phylo = asv) |>
      dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits)))
    
    # Save the combined data for further analysis
    combined_data <- asv_tab_10y_3_rel |>
      dplyr::filter(asv_num == asv & str_detect(sample_id, '_3_')) |>
      left_join(m_3_red, by = 'sample_id') |>
      left_join(weighted_traits, by = 'sample_id_num') |>
      dplyr::mutate(asv_num = asv)
    
    # Append the combined data to the list
    combined_data_list[[asv]] <- combined_data
  } else {
    message(paste("Species labels in 'comp' and 'abund' are not identical and ordered for", asv))
  }
}

# Combine the list elements into a single tibble
combined_data_all_3 <- bind_rows(combined_data_list) |>
  dplyr::left_join(bloom_event)  |>
  #dplyr::mutate(fraction = '3') |>
  distinct() |>
  dplyr::select(asv_num, sample_id, relative_abundance, fraction, date, phylogenetic_distance, bloom_event) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num'))
  
combined_data_all_3 |>
  colnames()

combined_data_all_02 |>
  colnames()

phylogenetic_weighted_tb <- combined_data_all_3 |>
  bind_rows(combined_data_all_02) |>
  dplyr::mutate(bloom_event = case_when(is.na(bloom_event) ~ 'no-bloom',
                                        TRUE ~ bloom_event))

phylogenetic_weighted_tb <- phylogenetic_weighted_tb |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num))

phylogenetic_weighted_tb$class_f <-  factor(phylogenetic_weighted_tb$class_f, 
                                            levels=unique(phylogenetic_weighted_tb$class_f[order(phylogenetic_weighted_tb$phylum_f)]), 
                                            ordered=TRUE)

phylogenetic_weighted_tb$order_f <-  factor(phylogenetic_weighted_tb$order_f, 
                                            levels=unique(phylogenetic_weighted_tb$order_f[order(phylogenetic_weighted_tb$phylum_f,
                                                                                                 phylogenetic_weighted_tb$class_f)]), 
                                            ordered=TRUE)

phylogenetic_weighted_tb$family_f <-  factor(phylogenetic_weighted_tb$family_f, 
                                             levels=unique(phylogenetic_weighted_tb$family_f[order(phylogenetic_weighted_tb$phylum_f,
                                                                                                   phylogenetic_weighted_tb$class_f,
                                                                                                   phylogenetic_weighted_tb$order_f)]), 
                                             ordered=TRUE)

phylogenetic_weighted_tb$asv_num_f <-  factor(phylogenetic_weighted_tb$asv_num_f, 
                                              levels=unique(phylogenetic_weighted_tb$asv_num_f[order(phylogenetic_weighted_tb$phylum_f,
                                                                                                     phylogenetic_weighted_tb$class_f,
                                                                                                     phylogenetic_weighted_tb$order_f,
                                                                                                     phylogenetic_weighted_tb$family_f)]), 
                                              ordered=TRUE)

phylogenetic_weighted_tb$bloom_event <- factor(phylogenetic_weighted_tb$bloom_event)

more_blooms <- phylogenetic_weighted_tb |>
  dplyr::group_by(bloom_event, asv_num, fraction) |>
  dplyr::reframe(n = n()) |>
  dplyr::filter(bloom_event == 'bloom' &
                  n > 5) |>
  ungroup()

boxplot_phylogenetic_distance <- phylogenetic_weighted_tb  |>
  dplyr::filter(asv_num %in% c(occurrence_bloo_bbmo_red |>
                                 dplyr::filter(occurrence_perc > 2/3) %$%
                                 asv_num)) |>
  dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes(bloom_event, phylogenetic_distance, group = as.factor(asv_num_f), color = family_f))+
  geom_point(aes(color = family_f, shape = fraction), position = position_jitter(width = 0.1))+ #alpha = if_else(bloom_event == 'bloom', 1, 0.5)
  scale_x_discrete(labels = labs_bloom_events)+
  geom_boxplot(aes(group = bloom_event, alpha = 0.8))+
  scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  labs(x='Blooming event', y = 'Weighted phylogenetic distance', color = 'Family', shape = 'Fraction')+
  facet_wrap(vars(interaction(asv_num_f, family_f)))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend(ncol = 1),
         alpha = 'none',
         shape = guide_legend(ncol = 1))

boxplot_phylogenetic_distance

# ggsave(boxplot_phylogenetic_distance, file = 'boxplot_phylogenetic_distance_remove_asv_interest_remove0_v3.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 188, units = 'mm')

phylogenetic_weighted_tb |>
  ungroup() |>
  distinct(asv_num)

weigth_phylogenetic_distance_blooms_plot <- phylogenetic_weighted_tb  |>
  #dplyr::filter(occurrence_perc > 2/3) |>
  dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes( x= as.numeric(phylogenetic_distance), y =interaction(family_f, asv_num_f) ,fill = as.factor(bloom_event)), 
         group = as.factor(bloom_event))+
  geom_density_ridges(alpha = 0.8)+
  #scale_color_manual(values = palette_family_assigned_bloo)+
  labs(y ='Density', x = 'Weighted phylogenetic distance', fill = 'Bloom event')+
  #facet_wrap(vars(asv_num))+
  scale_fill_manual(values = c('bloom' = '#00808F' , 'no-bloom' = '#454545'), labels = labs_bloom_events)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        panel.border = element_blank())

weigth_phylogenetic_distance_blooms_plot

### I plot the distance between the community in the different fractions PA and FL -----
phylogenetic_weighted_tb_pa_fl_plot <- phylogenetic_weighted_tb  |>
  dplyr::filter(asv_num %in% c(occurrence_bloo_bbmo_red |>
                                 dplyr::filter(occurrence_perc > 2/3) %$%
                                 asv_num)) |>
  #dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  #dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes(interaction(fraction, asv_num_f), phylogenetic_distance, group = as.factor(asv_num_f), color = family_f))+
  geom_point(aes(color = family_f, shape = fraction), position = position_jitter(width = 0.1))+ #alpha = if_else(bloom_event == 'bloom', 1, 0.5)
  #scale_x_discrete(labels = labs_bloom_events)+
  geom_violin(aes(group = interaction(fraction, asv_num_f), alpha = 0.8))+
  scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  labs(x='Fraction', y = 'Weighted phylogenetic distance', color = 'Family', shape = 'Fraction')+
  #facet_wrap(vars(fraction))+
  coord_flip()+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend(ncol = 1),
         alpha = 'none',
         shape = guide_legend(ncol = 1))

phylogenetic_weighted_tb_pa_fl_plot


phylogenetic_weighted_tb_pa_fl_plot <- phylogenetic_weighted_tb  |>
  dplyr::filter(asv_num %in% c(occurrence_bloo_bbmo_red |>
                                 dplyr::filter(occurrence_perc > 2/3) %$%
                                 asv_num)) |>
  #dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  #dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes(interaction(fraction, asv_num_f), phylogenetic_distance, group = as.factor(asv_num_f), color = family_f))+
  geom_point(aes(color = family_f, shape = fraction, alpha = if_else(bloom_event == 'bloom', 1, 0.5)), position = position_jitter(width = 0.1))+ #
  #scale_x_discrete(labels = labs_bloom_events)+
  geom_violin(aes(group = interaction(fraction, asv_num_f), alpha = 0.8))+
  scale_shape_discrete(labels = labs_fraction)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  labs(x='Fraction', y = 'Weighted phylogenetic distance', color = 'Family', shape = 'Fraction')+
  #facet_wrap(vars(fraction))+
  coord_flip()+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend(ncol = 1),
         alpha = 'none',
         shape = guide_legend(ncol = 1))

phylogenetic_weighted_tb_pa_fl_plot

phylogenetic_weighted_tb_pa_fl_plot <- phylogenetic_weighted_tb  |>
  dplyr::filter(asv_num %in% c(occurrence_bloo_bbmo_red |>
                                 dplyr::filter(occurrence_perc > 2/3) %$%
                                 asv_num)) |>
  #dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  #dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes(x = phylogenetic_distance, y =  asv_num))+
  #geom_point(aes(color = fraction, shape = fraction, alpha = if_else(bloom_event == 'bloom', 1, 0.5)), position = position_jitter(width = 0.25))+ #
  #scale_x_discrete(labels = labs_bloom_events)+
  #geom_violin(aes(group = interaction(fraction, asv_num_f), alpha = 0.8))+
  geom_density_ridges(aes( y =  asv_num, fill = fraction, alpha = 0.5))+
  scale_shape_discrete(labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  labs(y ='ASV number', x = 'Weighted phylogenetic distance', shape = 'Fraction', fill = 'Fraction', color = 'Fraction')+
  #facet_wrap(vars(fraction))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend(ncol = 1),
         alpha = 'none',
         shape = guide_legend(ncol = 1))

phylogenetic_weighted_tb_pa_fl_plot

ggsave(phylogenetic_weighted_tb_pa_fl_plot, file = 'phylogenetic_weighted_tb_pa_fl_plot_abundant.pdf',
       path = 'Results/Figures/',
       width = 188, height = 188, units = 'mm')

phylogenetic_weighted_tb_pa_fl_plot <- phylogenetic_weighted_tb  |>
  dplyr::filter(asv_num %in% c(occurrence_bloo_bbmo_red |>
                                 dplyr::filter(occurrence_perc < 2/3) %$%
                                 asv_num)) |>
  #dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  #dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes(x = phylogenetic_distance, y =  asv_num))+
  #geom_point(aes(color = fraction, shape = fraction, alpha = if_else(bloom_event == 'bloom', 1, 0.5)), position = position_jitter(width = 0.25))+ #
  #scale_x_discrete(labels = labs_bloom_events)+
  #geom_violin(aes(group = interaction(fraction, asv_num_f), alpha = 0.8))+
  geom_density_ridges(aes( y =  asv_num, fill = fraction, alpha = 0.5))+
  scale_shape_discrete(labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  labs(y ='ASV number', x = 'Weighted phylogenetic distance', shape = 'Fraction', fill = 'Fraction', color = 'Fraction')+
  #facet_wrap(vars(fraction))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend(ncol = 1),
         alpha = 'none',
         shape = guide_legend(ncol = 1))

phylogenetic_weighted_tb_pa_fl_plot

ggsave(phylogenetic_weighted_tb_pa_fl_plot, file = 'phylogenetic_weighted_tb_pa_fl_plot_narrow_inter.pdf',
       path = 'Results/Figures/',
       width = 188, height = 188, units = 'mm')

## clear differences between PA and FL 

phylogenetic_weighted_tb  |>
  # dplyr::filter(asv_num %in% c(occurrence_bloo_bbmo_red |>
  #                                dplyr::filter(occurrence_perc < 2/3) %$%
  #                                asv_num)) |>
  dplyr::filter(asv_num %in% c('asv264', 'asv225', 'asv200', 'asv15')) |>
  #dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  #dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes(x = phylogenetic_distance, y =  asv_num))+
  #geom_point(aes(color = fraction, shape = fraction, alpha = if_else(bloom_event == 'bloom', 1, 0.5)), position = position_jitter(width = 0.25))+ #
  #scale_x_discrete(labels = labs_bloom_events)+
  #geom_violin(aes(group = interaction(fraction, asv_num_f), alpha = 0.8))+
  geom_density_ridges(aes(y = asv_num, fill = fraction, alpha = 0.5))+
  scale_shape_discrete(labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  labs(y ='ASV number', x = 'Weighted phylogenetic distance', shape = 'Fraction', fill = 'Fraction', color = 'Fraction')+
  facet_wrap(vars(interaction(order_f, family_f)))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'right',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend(ncol = 1),
         alpha = 'none',
         shape = guide_legend(ncol = 1))

## only previous to blooming event and bloom event --------

phylogenetic_weighted_tb |>
  colnames()

m_red <- m_02_red |>
  bind_rows(m_3_red)

previous_bloom_event <- phylogenetic_weighted_tb |>
  left_join(m_red) |>
  dplyr::group_by(asv_num_f, fraction, bloom_event) |>
  arrange(-as.numeric(sample_id_num)) |>
  dplyr::filter(bloom_event == 'bloom') |>
  dplyr::select(sample_id_num) |>
  dplyr::mutate(sample_id_num_previous = (as.numeric(sample_id_num)-1)) |>
  ungroup() |>
  dplyr::select(-sample_id_num, -bloom_event) 

previous_bloom_event_phylodistance <- phylogenetic_weighted_tb |>
  left_join(m_red) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)) |>
  ungroup() |> 
  dplyr::select(asv_num_f, fraction, sample_id_num, phylum_f, class_f, order_f, family_f, relative_abundance, phylogenetic_distance) |>
  right_join(previous_bloom_event, by = c('asv_num_f', 'fraction', 'sample_id_num' = 'sample_id_num_previous')) |>
  dplyr::mutate(bloom_event = 'pre-bloom')

phylogenetic_weighted_tb_bloom <- phylogenetic_weighted_tb |>
  left_join(m_red) |>
  ungroup() |> 
  dplyr::filter(bloom_event == 'bloom') |>
  dplyr::select(asv_num_f, fraction, sample_id_num, phylum_f, class_f, order_f, family_f, relative_abundance, phylogenetic_distance, bloom_event) |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)) |>
  bind_rows(previous_bloom_event_phylodistance)

# Assuming phylogenetic_weighted_tb is your tibble
weigth_phylogenetic_distance_blooms_plot <- phylogenetic_weighted_tb_bloom  |>
  #dplyr::filter(occurrence_perc > 2/3) |>
  #dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  #dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes( x= as.numeric(phylogenetic_distance), y =interaction(family_f, asv_num_f), fill = as.factor(bloom_event)), 
         group = as.factor(bloom_event))+
  geom_density_ridges(alpha = 0.8)+
  #scale_color_manual(values = palette_family_assigned_bloo)+
  labs(y ='Density', x = 'Weighted phylogenetic distance', fill = 'Bloom event')+
  #facet_wrap(vars(asv_num))+
  scale_fill_manual(values = c('bloom' = '#00808F' , 'pre-bloom' = '#454545'))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        panel.border = element_blank())

weigth_phylogenetic_distance_blooms_plot

phylogenetic_weighted_tb_bloom$bloom_event <- factor(phylogenetic_weighted_tb_bloom$bloom_event, levels = c('pre-bloom', 'bloom'))

phylogenetic_weighted_tb_bloom |>
dplyr::filter(asv_num_f %in% more_blooms$asv_num) |>
  ungroup() |>
  ggplot(aes(relative_abundance, phylogenetic_distance))+
  geom_point(aes(size = relative_abundance, shape = fraction,
                 color = bloom_event))+
  #geom_boxplot(alpha = 0.4)+
  geom_smooth(method = 'loess', color = 'black')+
  facet_wrap(vars(interaction(asv_num_f, family_f)))+
  scale_color_manual(values = c('bloom' = '#00808F' , 'pre-bloom' = '#454545'))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        panel.border = element_blank())

phylogenetic_distance_bloom_events_plot <- phylogenetic_weighted_tb |>
  #dplyr::filter(order_f == 'Flavobacteriales') |>
  #dplyr::filter(asv_num_f == 'asv237') |>
  dplyr::filter(!class_f %in% c('Phycisphaerae', 'Bdellovibrionia')) |>
  #dplyr::filter(asv_num_f %in% more_blooms$asv_num) |>
  #left_join(occurrence_bloo_bbmo, by = c('fraction', 'asv_num_f' = 'asv_num')) |>
  #dplyr::filter(order_f == 'Rhodobacterales') |>
  ggplot(aes(date, phylogenetic_distance))+
  #geom_line(aes(), alpha = 0.5)+
  facet_grid(class_f~fraction)+
  #facet_wrap(vars(asv_num_f))+
  geom_smooth(aes(group = asv_num_f, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.2), shape = fraction, size = relative_abundance))+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(x = 'Time', y = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol = 2),
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 13),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_blank(),
       # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
       legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_distance_bloom_events_plot

# ggsave(phylogenetic_distance_bloom_events_plot, file = 'phylogenetic_distance_bloom_events_plot_v2.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 250, units = 'mm')

phylogenetic_distance_bloom_events_plot <- phylogenetic_weighted_tb |>
  dplyr::filter(class_f %in% c('Phycisphaerae', 'Bdellovibrionia')) |>
  #dplyr::filter(asv_num_f %in% more_blooms$asv_num) |>
  #left_join(occurrence_bloo_bbmo_ed, by = c('fraction', 'asv_num_f' = 'asv_num')) |>
  #dplyr::filter(order_f == 'Rhodobacterales') |>
  ggplot(aes(date, phylogenetic_distance))+
  #geom_line(aes(), alpha = 0.5)+
  facet_grid(vars(order_f))+
  geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.2), shape = fraction, size = relative_abundance))+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(x = 'Time', y = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol = 2),
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 13),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_blank(),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_distance_bloom_events_plot

phylogenetic_distance_bloom_events_plot_ciano <- phylogenetic_weighted_tb |>
  dplyr::filter(class_f %in% c('Cyanobacteriia')) |>
  #dplyr::filter(asv_num_f %in% more_blooms$asv_num) |>
  #left_join(occurrence_bloo_bbmo_ed, by = c('fraction', 'asv_num_f' = 'asv_num')) |>
  #dplyr::filter(order_f == 'Rhodobacterales') |>
  ggplot(aes(date, phylogenetic_distance))+
  facet_grid(vars(asv_num_f))+
  geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.4), shape = fraction, size = relative_abundance))+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(x = 'Time', y = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol = 2),
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 13),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_blank(),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_distance_bloom_events_plot_ciano

phylogenetic_weighted_tb |>
  colnames()

## do bloom events for each ASV happen at the same phylogenetic distance with the community? ----
phylogenetic_weighted_tb |>
  dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  ggplot(aes(bloom_event, phylogenetic_distance))+
  geom_point(aes(shape = 'fraction'))+
  facet_wrap(vars(asv_num_f))


phylogenetic_weighted_tb |>
  dplyr::mutate(high_abund = case_when(relative_abundance > 0.1 ~ 'high-abund',
                TRUE ~ 'low-abund')) |>
  dplyr::  mutate(phylo_distance_cat = case_when(
    phylogenetic_distance < 1.5 ~ 'low',
    phylogenetic_distance >= 1.5 & phylogenetic_distance < 2.75 ~ 'mid',
    phylogenetic_distance >= 2.75 ~ 'super-high'
  )) |>
  #dplyr::filter(!is.na(year)) |>
  #dplyr::filter(order == 'Rhodobacterales') |>
  ggplot(aes(phylo_distance_cat, relative_abundance))+
  #geom_line(aes(group = asv_num_f), alpha = 0.5)+
  geom_boxplot()+
  geom_point(aes(alpha = if_else(high_abund == 'high-abund', 1, 0.2), shape = fraction))+
  facet_wrap(vars(fraction))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 14),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        panel.border = element_blank())

## new analysis
library(betapart)

### environmental distance vs phylogenetic distance
euclidean_distance_tb

phylogenetic_weighted_tb |>
  colnames()

asv_tab_all_bloo_z_tax |>
  colnames()

decimal_date_tb_3 <- asv_tab_10y_3_pseudo |>
  dplyr::select(decimal_date, sample_id) |>
  distinct()

decimal_date_tb <- asv_tab_10y_02_pseudo |>
  dplyr::select(decimal_date, sample_id) |>
  distinct() |>
  bind_rows(decimal_date_tb_3)

phylogenetic_distance_vs_env_distance_plot <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(euclidean_distance_tb) |>
  ggplot(aes(phylogenetic_distance, euclidean_distance))+
  #geom_point(aes(size = relative_abundance))+
  facet_wrap(vars(interaction(family_f,asv_num_f)))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.2), shape = fraction, #size = relative_abundance/2, 
                 color = order_f))+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Euclidean distance previous sampling point', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol = 2),
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 13),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #strip.text.x = element_blank(),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_distance_vs_env_distance_plot
  
# ggsave(phylogenetic_distance_vs_env_distance_plot, file = 'phylogenetic_distance_vs_env_distance_plot.pdf',
#        path = 'Results/Figures/',
#        width = 250, height = 250, units = 'mm')

## witn NMDS env variables ----
nmds_10y_red <- nmds_10y_env |>
  dplyr::select(MDS1, MDS2, decimal_date)

phylogenetic_distance_vs_env_distance_plot <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(nmds_10y_red) |>
  ggplot(aes(phylogenetic_distance, MDS1))+
  #geom_point(aes(size = relative_abundance))+
  facet_wrap(vars(interaction(family_f,asv_num_f)))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.2), shape = fraction, #size = relative_abundance/2, 
                 color = order_f))+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Euclidean distance previous sampling point', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol = 2),
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 13),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #strip.text.x = element_blank(),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_distance_vs_env_distance_plot

phylogenetic_distance_vs_env_distance_plot <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(nmds_10y_red) |>
  ggplot(aes(phylogenetic_distance, MDS2))+
  #geom_point(aes(size = relative_abundance))+
  facet_wrap(vars(interaction(family_f,asv_num_f)))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.2), shape = fraction, #size = relative_abundance/2, 
                 color = order_f))+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Euclidean distance previous sampling point', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol = 2),
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 13),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        #strip.text.x = element_blank(),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_distance_vs_env_distance_plot

### nmds with all data ----
nmds_10y_red <- nmds_10y_env |>
  dplyr::select(MDS1, MDS2, decimal_date)

phylogenetic_distance_vs_env_distance_plot <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(nmds_10y_red) |>
  ggplot(aes(phylogenetic_distance, MDS1))+
  #geom_point(aes(size = relative_abundance))+
  facet_wrap(vars(interaction(family_f,asv_num_f)), ncol = 4)+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.2), shape = fraction, #size = relative_abundance/2, 
                 color = order_f))+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Environmental NMDS axis', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol = 2),
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 13),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_distance_vs_env_distance_plot

# ggsave(phylogenetic_distance_vs_env_distance_plot, file = 'phylogenetic_distance_vs_env_distance_plot_nmds.pdf',
#        path = 'Results/Figures/',
#        width = 250, height = 450, units = 'mm')

## another idea of plot 
phylogenetic_distance_vs_env_distance_plot <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(nmds_10y_red) |>
  ggplot(aes(phylogenetic_distance, MDS1))+
  #geom_point(aes(size = relative_abundance))+
  #facet_grid(vars(interaction(family_f,asv_num_f), cols = 4 ))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.2), shape = fraction, #size = relative_abundance/2, 
                 color = order_f))+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Environmental NMDS axis', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol = 2),
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 13),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_distance_vs_env_distance_plot

box_plot_phylogenetic_distance <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  #left_join(nmds_10y_red) |>
  ggplot(aes(phylogenetic_distance, order_f))+
  #geom_point(aes(size = relative_abundance))+
  #facet_grid(vars(interaction(family_f,asv_num_f), cols = 4 ))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.15), shape = fraction, #size = relative_abundance/2, 
                 color = order_f),
             position = position_jitter(width = 0.2))+
  geom_violin(aes(group = order_f), alpha = 0.1)+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Order', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = 'none',
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(1, 'mm'),
        legend.text = element_text(size = 9),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

box_plot_phylogenetic_distance

# ggsave(box_plot_phylogenetic_distance, file = 'box_plot_phylogenetic_distance.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 200, units = 'mm')

## edited to compare with the CCM results ----
phylogenetic_weighted_tb |>
  ungroup() |>
  distinct(asv_num)
### I need to recover the information for all those taxa that are not BLOOMERS!!! 

box_plot_phylogenetic_distance <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  #left_join(nmds_10y_red) |>
  ggplot(aes(phylogenetic_distance, order_f))+
  #geom_point(aes(size = relative_abundance))+
  facet_grid(vars(fraction))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.15), shape = fraction, #size = relative_abundance/2, 
                 color = order_f),
             position = position_jitter(width = 0.2))+
  geom_violin(aes(group = order_f), alpha = 0.1)+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Order', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = 'none',
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(1, 'mm'),
        legend.text = element_text(size = 9),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

box_plot_phylogenetic_distance

ggsave(box_plot_phylogenetic_distance, file = 'box_plot_phylogenetic_distance_v2.pdf',
       path = 'Results/Figures/',
       width = 188, height = 200, units = 'mm')

box_plot_phylogenetic_distance <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(nmds_10y_red) |>
  ggplot(aes(phylogenetic_distance, family_f))+
  #geom_point(aes(size = relative_abundance))+
  #facet_grid(vars(interaction(family_f,asv_num_f), cols = 4 ))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.15), shape = fraction, #size = relative_abundance/2, 
                 color = order_f),
             position = position_jitter(width = 0.2))+
  geom_violin(aes(group = family_f), alpha = 0.1)+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Order', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = 'none',
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(1, 'mm'),
        legend.text = element_text(size = 9),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

box_plot_phylogenetic_distance

ggsave(box_plot_phylogenetic_distance, file = 'box_plot_phylogenetic_distance_family.pdf',
       path = 'Results/Figures/',
       width = 188, height = 230, units = 'mm')
  
box_plot_phylogenetic_distance <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(nmds_10y_red) |>
  ggplot(aes(phylogenetic_distance, interaction(family_f, asv_num_f)))+
  #geom_point(aes(size = relative_abundance))+
  #facet_grid(vars(interaction(family_f,asv_num_f), cols = 4 ))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.15), shape = fraction, #size = relative_abundance/2, 
                 color = order_f),
             position = position_jitter(width = 0.2))+
  geom_violin(aes(group = asv_num_f), alpha = 0.1)+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Order', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = 'none',
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(1, 'mm'),
        legend.text = element_text(size = 9),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

box_plot_phylogenetic_distance

# ggsave(box_plot_phylogenetic_distance, file = 'box_plot_phylogenetic_distance_asv_num_f.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 260, units = 'mm')

## env variables vs phylogenetic distance ----
phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(env_data_interpolated_values_all) |>
  ggplot(aes(temperature_no_nas, phylogenetic_distance))+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.15), shape = fraction, #size = relative_abundance/2, 
                 color = order_f),
             position = position_jitter(width = 0.2))+
  #geom_violin(aes(group = asv_num_f), alpha = 0.1)+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  facet_wrap(vars(asv_num), scales = 'free')+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Order', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = 'none',
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(1, 'mm'),
        legend.text = element_text(size = 9),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

## only non-recurrent blooms relationship with the community -----
bloo_all_types_summary_tb2 <- bloo_all_types_summary_tb |>
  dplyr::select(-X) |>
  dplyr::mutate(fraction = as.factor(fraction))

phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(bloo_all_types_summary_tb2) |>
  dplyr::filter(recurrency == 'no') |>
  ggplot(aes(phylogenetic_distance, interaction(family_f, asv_num_f)))+
  #geom_point(aes(size = relative_abundance))+
  #facet_grid(vars(interaction(family_f,asv_num_f), cols = 4 ))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(#alpha = if_else(bloom_event == 'bloom', 1, 0.15), 
                 shape = fraction, size = relative_abundance/2, 
                 color = order_f),
             position = position_jitter(width = 0.2))+
  geom_violin(aes(group = asv_num_f), alpha = 0.1)+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Order', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = 'none',
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(1, 'mm'),
        legend.text = element_text(size = 9),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())
  
library(ggridges)

density_ridges_plot_phylogenetic_distance <- phylogenetic_weighted_tb |>
  left_join(decimal_date_tb) |>
  left_join(bloo_all_types_summary_tb2) |>
  dplyr::filter(recurrency == 'no') |>
  ggplot(aes(y = interaction(family_f, asv_num_f), x = phylogenetic_distance, fill = order_f, group = asv_num_f, color = order_f))+
  #geom_point(aes(size = relative_abundance))+
  #facet_grid(vars(interaction(family_f,asv_num_f), cols = 4 ))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.15),
    shape = fraction, #size = relative_abundance/2,
    color = order_f),
    position = position_jitter(width = 0.15))+
  #geom_density_ridges(aes(group = asv_num_f, alpha = 0.3))+
  geom_density_ridges(alpha = 0.3, panel_scaling = TRUE, #scale = 1,
                      #jittered_points = TRUE,
                      #point_shape = fraction, 
                      #point_size = 0.2, point_alpha = if_else(bloom_event == 'bloom', 1, 0.15),
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  #geom_violin(aes(group = asv_num_f), alpha = 0.1)+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Density', 
       x = 'ASV weighted phylogenetic distance', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = 'none',
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(1, 'mm'),
        legend.text = element_text(size = 9),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

density_ridges_plot_phylogenetic_distance

# ggsave(density_ridges_plot_phylogenetic_distance, file = 'density_ridges_plot_phylogenetic_distance.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 280, units = 'mm')

## phylogenetic distances in rCLR and for ALL taxa either potential bloomers or not ---------
## weighted phylogenetic distances removing the ASV of interest from the analysis for all the community ----
## FL ---
m_02_red <- m_02 |>
  dplyr::select(sample_id, date, fraction, sample_id_num) |>
  distinct()

asv_tab_10y_02_rel

asv_tab_10y_02_rclr

## I do it for high occurrent taxa in the dataset occurrence > 2/3 either in the FL or the PA fraction ----
occurrence_bbmo <- asv_tab_10y_02_rel |> 
  group_by(asv_num) |>
  dplyr::mutate(n_occurrence = case_when(relative_abundance > 0 ~ 1,
                                         relative_abundance == 0 ~ 0)) |>
  dplyr::mutate(n_samples = n_distinct(sample_id)) |>
  group_by(asv_num) |>
  dplyr::mutate(occurrence_perc = sum(n_occurrence)/n_samples)

occurrence_bbmo_abund_02 <- occurrence_bbmo |>
  dplyr::filter(occurrence_perc > 1/3) |>
  distinct(asv_num)

occurrence_bbmo <- asv_tab_10y_3_rel |> 
  group_by(asv_num) |>
  dplyr::mutate(n_occurrence = case_when(relative_abundance > 0 ~ 1,
                                         relative_abundance == 0 ~ 0)) |>
  dplyr::mutate(n_samples = n_distinct(sample_id)) |>
  group_by(asv_num) |>
  dplyr::mutate(occurrence_perc = sum(n_occurrence)/n_samples)

occurrence_bbmo_abund_3 <- occurrence_bbmo |>
  dplyr::filter(occurrence_perc > 1/3) |>
  distinct(asv_num)

occurrence_bbmo_abund <- occurrence_bbmo_abund_02 |>
  bind_rows(occurrence_bbmo_abund_3) |>
  distinct(asv_num) 

bbmo_10y_asv_num_tax <- asv_tab_10y_02_rel |>
  ungroup() |>
  distinct(asv_num) |>
  left_join(tax_bbmo_10y_new) 

bbmo_10y_asv_num_tax_abund <- bbmo_10y_asv_num_tax |>
  dplyr::filter(asv_num %in% occurrence_bbmo_abund$asv_num)

# Create an empty list to store combined data for each ASV (this is with rCLR)

combined_data_list_all <- list()

# Loop over each `asv_num` in `bloo_02`
for (asv in bbmo_10y_asv_num_tax_abund$asv_num) {
  # dplyr::filter comp based on the current asv_num
  comp <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv) |>
    dplyr::filter(asv_num_2 != asv) |> # remove the ASV of interest from the analysis
    dplyr::select(-asv_num_1)
  
  # I need to transform 0 by NA so that the functcomp function does not consider them
  asv_tab_10y_02_rel_ed <- asv_tab_10y_02_rel |>
    dplyr::mutate(relative_abundance = case_when(
      relative_abundance == 0 ~ NA_real_,  # Replace 0 with NA (numeric NA)
      TRUE ~ relative_abundance            # Keep the original value otherwise
    ))
  
  # Prepare abund tibble
  abund <- asv_tab_10y_02_rel_ed |>
    dplyr::select(-total_reads) |>
    pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') %>%
    dplyr::select(-asv) |> # remove the ASV of interest from the analysis
    ungroup() %>%
    dplyr::select(-sample_id)
  
  # Ensure all species labels match and are ordered
  comp <- comp |>
    arrange(asv_num_2)
  
  if (all(comp$asv_num_2 %in% colnames(abund))) {
    # Reorder columns in abund
    abund <- abund[, comp$asv_num_2]
  } else {
    message(paste("Not all elements in comp$asv_num_2 match colnames in abund for", asv))
    next
  }
  
  # Convert tibbles to matrices
  comp_matrix <- comp |>
    dplyr::select(-asv_num_2) |>
    as.matrix()
  
  rownames(comp_matrix) <- comp$asv_num_2
  
  abund_matrix <- as.matrix(abund)
  
  # Ensure the column names of abund_matrix match the row names of comp_matrix
  if (all(comp$asv_num_2 == colnames(abund_matrix))) {
    # Call the functcomp function
    weighted_traits <- FD::functcomp(x = comp_matrix, a = abund_matrix)
    
    # Rename the column
    colname <- paste0(asv, "_phylogenetic_distance")
    weighted_traits <- weighted_traits |>
      as_tibble() |>
      dplyr::mutate(asv_num_phylo = asv) |>
      dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits)))
    
    # Save the combined data for further analysis
    combined_data <- asv_tab_10y_02_rel |>
      dplyr::filter(asv_num == asv & str_detect(sample_id, '_0.2_')) |>
      left_join(m_02_red, by = 'sample_id') |>
      left_join(weighted_traits, by = 'sample_id_num') |>
      dplyr::mutate(asv_num = asv)
    
    # Append the combined data to the list
    combined_data_list_all[[asv]] <- combined_data
  } else {
    message(paste("Species labels in 'comp' and 'abund' are not identical and ordered for", asv))
  }
}

# Combine the list elements into a single tibble
combined_data_all_02_rel <- bind_rows(combined_data_list_all) |>
  dplyr::left_join(bloom_event) |>
  #dplyr::mutate(fraction = '0.2') |>
  distinct() |>
  dplyr::select(asv_num, sample_id, relative_abundance, fraction, date, phylogenetic_distance, bloom_event) |>
  left_join(bbmo_10y_asv_num_tax, by = c('asv_num' = 'asv_num'))

combined_data_all_02_rel |>
  ungroup() |>
  distinct(asv_num)

# Generate plots outside the loop----
combined_data_all_02_rel |>
  colnames()

## plot the relationship between relative abundance and weighted phylogenetic distance with the community -----
#plot
# combined_data_all_02_rel |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = asv_num, color = family))+
#   geom_point(aes(color = family, shape = bloom_event, alpha = if_else(bloom_event == 'bloom', 1, 0.5)))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   #facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 16),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# combined_data_all_02_rel |>
#   dplyr::filter(bloom_event == 'bloom') |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family_f))+
#   geom_point(aes(color = family_f, shape = bloom_event))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 8),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# combined_data_all_02_rel |>
#   dplyr::filter(bloom_event == 'no-bloom') |>
#   ggplot(aes(phylogenetic_distance, relative_abundance, group = as.factor(asv_num), color = family_f))+
#   geom_point(aes(color = family_f, shape = bloom_event))+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 8),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))
# 
# combined_data_all_02_rel |>
#   dplyr::mutate(bloomer = case_when(asv_num %in% bloo_02$value ~ 'bloomer',
#                                     TRUE ~ 'no-bloomer')) |>
#   dplyr::filter(family %in% bloo_taxonomy$family_f) |>
#   dplyr::mutate(bloom_event = case_when(
#     bloom_event == NA_character_ ~ 'no-bloom',  # Replace  NA with 0 (numeric NA)
#     TRUE ~ bloom_event        # Keep the original value otherwise
#   )) |>
#   ggplot(aes(phylogenetic_distance, family, group = as.factor(asv_num), color = family))+
#   geom_point(aes(color = family, shape = bloom_event))+
#   geom_density_ridges(aes(fill = family, group = asv_num), alpha = 0.2)+
#   #geom_smooth(method = 'loess', aes(group = asv_num, color = family))+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(bloomer))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = 'transparent'),
#         legend.position = 'bottom',
#         panel.grid = element_blank(), text = element_text(size = 8),
#         strip.text = element_text(margin = margin(2, 2, 2, 2)),
#         plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
#         legend.key.size = unit(3, 'mm'))

## The same for PA-------
# Create an empty list to store combined data for each ASV

### I filter by those taxa that have more tan 2/3 of occurrence in at least FL or PA
bbmo_10y_asv_num_tax <- asv_tab_10y_3_rel |>
  ungroup() |>
  distinct(asv_num) |>
  left_join(tax_bbmo_10y_new) 

bbmo_10y_asv_num_tax_abund <- bbmo_10y_asv_num_tax |>
  dplyr::filter(asv_num %in% occurrence_bbmo_abund$asv_num)

combined_data_list_all <- list()

m_3_red <- m_3 |>
  dplyr::select(date, fraction, sample_id, sample_id_num)

# asv_num_3 <- asv_tab_10y_3_rel |>
#   ungroup() |>  
#   dplyr::filter(relative_abundance > 0) |>
#   distinct(asv_num)

# Loop over each `asv_num` in `bloo_3`
for (asv in bbmo_10y_asv_num_tax_abund$asv_num) {
  # dplyr::filter comp based on the current asv_num
  comp <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv) |>
    dplyr::filter(asv_num_2 != asv) |> # remove the ASV of interest from the analysis
    dplyr::select(-asv_num_1)
  
  # I need to transform 0 by NA so that the functcomp function does not consider them
  asv_tab_10y_3_rel_ed <- asv_tab_10y_3_rel |>
    dplyr::mutate(relative_abundance = case_when(
      relative_abundance == 0 ~ NA_real_,  # Replace 0 with NA (numeric NA)
      TRUE ~ relative_abundance            # Keep the original value otherwise
    ))
  
  # Prepare abund tibble
  abund <- asv_tab_10y_3_rel_ed |>
    dplyr::select(-total_reads) |>
    pivot_wider(id_cols = 'sample_id', values_from = 'relative_abundance', names_from = 'asv_num') %>%
    dplyr::select(-asv) |> # remove the ASV of interest from the analysis
    ungroup() |>
    dplyr::select(-sample_id)
  
  # Ensure all species labels match and are ordered
  comp <- comp %>%
    arrange(asv_num_2)
  
  if (all(comp$asv_num_2 %in% colnames(abund))) {
    # Reorder columns in abund
    abund <- abund[, comp$asv_num_2]
  } else {
    message(paste("Not all elements in comp$asv_num_2 match colnames in abund for", asv))
    next
  }
  
  # Convert tibbles to matrices
  comp_matrix <- comp %>%
    dplyr::select(-asv_num_2) %>%
    as.matrix()
  
  rownames(comp_matrix) <- comp$asv_num_2
  
  abund_matrix <- as.matrix(abund)
  
  # Ensure the column names of abund_matrix match the row names of comp_matrix
  if (all(comp$asv_num_2 == colnames(abund_matrix))) {
    # Call the functcomp function
    weighted_traits <- FD::functcomp(x = comp_matrix, a = abund_matrix)
    
    # Rename the column
    colname <- paste0(asv, "_phylogenetic_distance")
    weighted_traits <- weighted_traits |>
      as_tibble() |>
      dplyr::mutate(asv_num_phylo = asv) |>
      dplyr::mutate(sample_id_num = as.character(1:nrow(weighted_traits)))
    
    # Save the combined data for further analysis
    combined_data <- asv_tab_10y_3_rel |>
      dplyr::filter(asv_num == asv & str_detect(sample_id, '_3_')) |>
      left_join(m_3_red, by = 'sample_id') |>
      left_join(weighted_traits, by = 'sample_id_num') |>
      dplyr::mutate(asv_num = asv)
    
    # Append the combined data to the list
    combined_data_list_all[[asv]] <- combined_data
  } else {
    message(paste("Species labels in 'comp' and 'abund' are not identical and ordered for", asv))
  }
}

# Combine the list elements into a single tibble
combined_data_all_3 <- bind_rows(combined_data_list_all) |>
  dplyr::left_join(bloom_event)  |>
  #dplyr::mutate(fraction = '3') |>
  distinct() |>
  dplyr::select(asv_num, sample_id, relative_abundance, fraction, date, phylogenetic_distance, bloom_event) |>
  left_join(bbmo_10y_asv_num_tax, by = c('asv_num' = 'asv_num'))

phylogenetic_weighted_tb_all <- combined_data_all_3 |>
  bind_rows(combined_data_all_02_rel) |>
  dplyr::mutate(bloom_event = case_when(is.na(bloom_event) ~ 'no-bloom',
                                        TRUE ~ bloom_event))

# phylogenetic_weighted_tb_all$asv_num |>
#   unique()

phylogenetic_weighted_tb_all <- phylogenetic_weighted_tb_all |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num)) |>
  dplyr::mutate(asv_num_f_fa = paste0(asv_num_f,'.', family_f)) |>
  dplyr::mutate(asv_num_f_fa = as.factor(asv_num_f_fa))

phylogenetic_weighted_tb_all$class_f <-  factor(phylogenetic_weighted_tb_all$class_f, 
                                            levels=unique(phylogenetic_weighted_tb_all$class_f[order(phylogenetic_weighted_tb_all$phylum_f)]), 
                                            ordered=TRUE)

phylogenetic_weighted_tb_all$order_f <-  factor(phylogenetic_weighted_tb_all$order_f, 
                                            levels=unique(phylogenetic_weighted_tb_all$order_f[order(phylogenetic_weighted_tb_all$phylum_f,
                                                                                                 phylogenetic_weighted_tb_all$class_f)]), 
                                            ordered=TRUE)

phylogenetic_weighted_tb_all$family_f <-  factor(phylogenetic_weighted_tb_all$family_f, 
                                             levels=unique(phylogenetic_weighted_tb_all$family_f[order(phylogenetic_weighted_tb_all$phylum_f,
                                                                                                   phylogenetic_weighted_tb_all$class_f,
                                                                                                   phylogenetic_weighted_tb_all$order_f)]), 
                                             ordered=TRUE)

phylogenetic_weighted_tb_all$asv_num_f <-  factor(phylogenetic_weighted_tb_all$asv_num_f, 
                                              levels=unique(phylogenetic_weighted_tb_all$asv_num_f[order(phylogenetic_weighted_tb_all$phylum_f,
                                                                                                     phylogenetic_weighted_tb_all$class_f,
                                                                                                     phylogenetic_weighted_tb_all$order_f,
                                                                                                     phylogenetic_weighted_tb_all$family_f)]), 
                                              ordered=TRUE)

phylogenetic_weighted_tb_all$asv_num_f_fa <-  factor(phylogenetic_weighted_tb_all$asv_num_f_fa, 
                                                  levels=unique(phylogenetic_weighted_tb_all$asv_num_f_fa[order(phylogenetic_weighted_tb_all$phylum_f,
                                                                                                             phylogenetic_weighted_tb_all$class_f,
                                                                                                             phylogenetic_weighted_tb_all$order_f,
                                                                                                             phylogenetic_weighted_tb_all$family_f)]), 
                                                  ordered=TRUE)

phylogenetic_weighted_tb_all$bloom_event <- factor(phylogenetic_weighted_tb_all$bloom_event)

## plots weigthed differences between fractions and blooms ----
phylogenetic_weighted_tb_all <- 
  phylogenetic_weighted_tb_all |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'Bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num ~ 'No-bloomer'))

phylogenetic_weighted_tb_pa_fl_plot <- phylogenetic_weighted_tb_all  |>
  # dplyr::filter(asv_num %in% c(occurrence_bloo_bbmo_red |>
  #                                dplyr::filter(occurrence_perc > 2/3) %$%
  #                                asv_num)) |>
  #dplyr::filter(asv_num %in% more_blooms$asv_num) |>
  #dplyr::filter(asv_num != 'asv15') |> # it covaries with asv7
  ggplot(aes(x = phylogenetic_distance, y =  asv_num_f_fa))+
  #geom_point(aes(color = fraction, shape = fraction, alpha = if_else(bloom_event == 'bloom', 1, 0.5)), position = position_jitter(width = 0.25))+ #
  #scale_x_discrete(labels = labs_bloom_events)+
  #geom_violin(aes(group = interaction(fraction, asv_num_f), alpha = 0.8))+
  geom_density_ridges(aes( y =  asv_num_f_fa, fill = fraction, alpha = 0.5), scale = 3)+
  scale_shape_discrete(labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  labs(y ='ASV number', x = 'Weighted phylogenetic distance', shape = 'Fraction', fill = 'Fraction', color = 'Fraction')+
  facet_wrap(vars(bloomer))+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2), size = 10),
        #plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        axis.text.x = element_text(size = 8),
        legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend(ncol = 1),
         alpha = 'none',
         shape = guide_legend(ncol = 2))

phylogenetic_weighted_tb_pa_fl_plot

# ggsave(phylogenetic_weighted_tb_pa_fl_plot, file = 'density_phylogenetic_weighted_tb_pa_fl_plot_1-3_occ.pdf',
#        path = 'Results/Figures/',
#        width = 188, height = 180, units = 'mm')

## we calculate weighted unifrac distance to observe the dynamics of it ----
# We introduce here a new method for computing differences between microbial communities based on 
# phylogenetic information. This method, UniFrac, measures the phylogenetic distance between sets 
# of taxa in a phylogenetic tree as the fraction of the branch length of the tree that leads 
# to descendants from either one environment or the other, but not both. UniFrac can be used 
# to determine whether communities are significantly different, to compare many communities 
# simultaneously using clustering and ordination techniques, and to measure the relative 
# contributions of different factors, such as chemistry and geography, to similarities between samples. (Lozupone et al., 2005) 

weighted_unifrac <- UniFrac(bbmo_10y, weighted = T)

# Assuming you have a phyloseq object named `physeq`
# Calculate UniFrac distance matrix (unweighted or weighted)
unifrac_dist <- phyloseq::distance(bbmo_10y, 
                                   method = "wunifrac") # or "wunifrac" for weighted

# Convert the distance matrix to a tidy tibble
unifrac_tibble <- as.matrix(unifrac_dist) |>
  as.data.frame() |>
  rownames_to_column(var = "sample_id_1") %>%
  pivot_longer(-sample_id_1, names_to = "sample_id_2", 
               values_to = "wunifrac_distance")

# Convert to tibble
unifrac_tibble <- unifrac_tibble |>
  as_tibble()

# plot weighted unifrac distances -----
unifrac_tibble_m <- unifrac_tibble |>
  left_join(m_red, by = c('sample_id_1' = 'sample_id')) |>
  dplyr::rename(sample_id_num_1 = sample_id_num) |>
  #dplyr::select(-date) |>
  left_join(m_red, by = c('sample_id_2' = 'sample_id')) |>
  dplyr::rename(sample_id_num_2 = sample_id_num) |>
  dplyr::mutate(sample_distance = (as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))) |>
  dplyr::filter(sample_distance == 1 &
                  fraction.x == fraction.y) |>
  dplyr::mutate(date = (as.POSIXct(date.x, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  dplyr::select(-date.y, -date.x)
  
weigthed_unifrac_plot <- unifrac_tibble |>
  left_join(m_red, by = c('sample_id_1' = 'sample_id')) |>
  dplyr::rename(sample_id_num_1 = sample_id_num) |>
  dplyr::select(-date) |>
  left_join(m_red, by = c('sample_id_2' = 'sample_id')) |>
  dplyr::rename(sample_id_num_2 = sample_id_num) |>
  dplyr::mutate(sample_distance = (as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))) |>
  dplyr::filter(sample_distance == 1 &
                  fraction.x == fraction.y) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  ggplot(aes(date, wunifrac_distance))+
  #facet_wrap(vars(fraction.x))+
  geom_line(aes(date, group = fraction.x, color = fraction.x, linetype = fraction.x))+
  scale_linetype_discrete(labels = labs_fraction)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = 'Weighted UNIFRAC distance', color = 'Fraction', linetype = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'none',
        panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
        # axis.text.x = element_text(size = 7), 
        panel.grid.major.y = element_blank())

weigthed_unifrac_plot

### compostion plot with env euclidean distance, evenness and phylogenetic distance + CCM
euclidean_distance_plot <- euclidean_distance_tb |>  
  ggplot(aes(date, euclidean_distance))+
  geom_line()+
  # facet_wrap(vars(season), labeller = labs_fraction_s, scales = 'free_x')+
  #scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(y = 'Environmental Euclidean distance', x = 'Date')+
  theme(legend.position = "right", panel.grid.minor = element_blank(),
       # axis.text.x = element_text(size = 7), 
        panel.grid.major.y = element_blank(), 
       strip.text = element_text(size = 7),
        #axis.text.y = element_text(size = 8),
        #axis.title = element_text(size = 8), 
       strip.background = element_blank(), 
       #  legend.text = element_text(size = 7), 
       # legend.title = element_text(size = 8), 
       axis.title = element_text(size = 10),
       strip.placement = 'outside')+
  guides(shape = 'none',
         color = guide_legend(ncol =1, size = 9,
                              override.aes = aes(label = ''))) 

euclidean_distance_plot

bray_curtis_m_tb <- community_eveness_all |> 
  left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(diversity_index == 'bray_curtis_result')

beta_diversity_plot <- community_eveness_all |> 
  left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(diversity_index == 'bray_curtis_result') |>
  ggplot(aes(date, value))+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  #geom_point(aes(shape = fraction, color = fraction))+
  geom_line(aes(date, value, group = fraction, color = fraction, linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  facet_grid(vars(diversity_index), labeller = labs_diversity)+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = 'Bray Curtis Dissimilarity', color = 'Fraction', linetype = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 10),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(), 
        strip.background = element_blank(), legend.position = 'none',
        strip.text = element_blank())

beta_diversity_plot

#legend_fraction <- get_legend(beta_diversity_plot)

beta_diversity_phylo_env_plot <- grid.arrange(beta_diversity_plot, 
             weigthed_unifrac_plot,
             legend_fraction,
             euclidean_distance_plot,
             heights=c(1,1, 0.25, 1),
             ncol = 1)

# ggsave(beta_diversity_phylo_env_plot,
#        filename = 'beta_diversity_phylo_env_plot.pdf',
#               path = 'Results/Figures/',
#               width = 180, height = 160, units = 'mm'
#        )

## Relationship env. distance phylogenetic distance and community distance ----
### is there a relationship between the environmental distance between samples and the 
### phylogenetic distance or the bray curtis distance between them?

## I have three different dataset with euclidean distances euclidean_distance_tb (all variables), euclidean_distance_tb_bio (bio parammeters)
## euclidean_distance_tb_phch (with temperature, day_length and nutrients)

m_bbmo_10y$season <- factor(m_bbmo_10y$season, levels = c('winter', 'spring', 'summer', 'autumn'))

euclidean_distance_tb$season

euclidean_distance_tb <- euclidean_distance_tb |>
  dplyr::select(#-fraction, 
               # -sample_id,
                -name_complete) |>
  ungroup()

## by season ----
euclidean_distance_tb_red <- euclidean_distance_tb |>
  dplyr::select(date, euclidean_distance)

community_eveness_all |> 
  left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(diversity_index == 'bray_curtis_result') |>
  #dplyr::filter(fraction == '3') |>
  left_join(euclidean_distance_tb_red) |>
 # dplyr::filter(is.na(euclidean_distance))
  ggplot(aes( euclidean_distance, value, color = season))+
  geom_point()+
  facet_wrap(vars(season))+
  geom_smooth(aes(group = season, color = season),method = 'loess')+
  labs(x = 'Bray Curtis Dissimilarity', y = 'Euclidean Environmental Distance', color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(aspect.ratio = 4/4,
        panel.grid.minor = element_blank(),
        axis.title  = element_text(size = 10),
        # axis.text.x = element_text(size = 7), 
        #panel.grid.major.y = element_blank(), 
        strip.background = element_blank(), legend.position = 'none',
        strip.text = element_blank())

unifrac_tibble |>
  left_join(m_red, by = c('sample_id_1' = 'sample_id')) |>
  dplyr::rename(sample_id_num_1 = sample_id_num) |>
  dplyr::select(-date) |>
  left_join(m_red, by = c('sample_id_2' = 'sample_id')) |>
  dplyr::rename(sample_id_num_2 = sample_id_num) |>
  dplyr::mutate(sample_distance = (as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))) |>
  dplyr::filter(sample_distance == 1 &
                  fraction.x == fraction.y) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  left_join(euclidean_distance_tb) |>
  ungroup() |>
  left_join(m_bbmo_10y, by = c('date', 'season')) |>
  dplyr::filter(!is.na(season)) |>
  ggplot(aes(wunifrac_distance, euclidean_distance, color = season))+
  facet_wrap(vars(season))+
  geom_point()+
  geom_smooth(aes(group = season, color = season),method = 'loess')+
  scale_color_manual(values = palette_seasons_4)+
  labs(x = 'Weigthed UniFrac distance', y = 'Euclidean Environmental Distance', color = 'Season')+
  theme_bw()+
  theme(aspect.ratio = 4/4,
        panel.grid.minor = element_blank(),
        axis.title  = element_text(size = 10),
        # axis.text.x = element_text(size = 7), 
        #panel.grid.major = element_blank(), 
        strip.background = element_blank(), legend.position = 'none',
        strip.text = element_blank())

m_bbmo_10y |>
  colnames()

  unifrac_tibble_ed <- unifrac_tibble |>
  left_join(m_red, by = c('sample_id_1' = 'sample_id')) |>
  dplyr::rename(sample_id_num_1 = sample_id_num) |>
  dplyr::select(-date) |>
  left_join(m_red, by = c('sample_id_2' = 'sample_id')) |>
  dplyr::rename(sample_id_num_2 = sample_id_num) |>
  dplyr::mutate(sample_distance = (as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))) |>
  dplyr::filter(sample_distance == 1 &
                  fraction.x == fraction.y) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
    dplyr::select(fraction = fraction.x, wunifrac_distance, date)

  community_eveness_all |> 
    left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
    dplyr::select(-row_index_2) |>
    pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
    left_join(m_bbmo_10y, by = c('sample_id')) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::filter(diversity_index == 'bray_curtis_result') |> 
    left_join(unifrac_tibble_ed) |>
    ggplot(aes(wunifrac_distance, value))+
    geom_point(aes(color = season))+
    labs(x = 'Weighted Unifrac Distance', y = 'Bray Curtis Dissimilarity')+
    scale_color_manual(values = palette_seasons_4)+
    geom_smooth(aes(group = season, color = season),method = 'loess')+
    theme_bw()+
    theme(aspect.ratio = 4/4,
          panel.grid.minor = element_blank(),
          axis.title  = element_text(size = 10),
          # axis.text.x = element_text(size = 7), 
          #panel.grid.major.y = element_blank(), 
          strip.background = element_blank(), legend.position = 'none',
          strip.text = element_blank())
  
## by year ----
  m_bbmo_10y$year <- factor(m_bbmo_10y$year, levels = c('2004', '2005', '2006', '2007',
                                                        '2008', '2009', '2010', '2011',
                                                        '2012', '2013'))
  
  euclidean_distance_tb$year <- factor(euclidean_distance_tb$year, levels = c('2004', '2005', '2006', '2007',
                                                                       '2008', '2009', '2010', '2011',
                                                                       '2012', '2013'))
  community_eveness_all |> 
    left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
    dplyr::select(-row_index_2) |>
    pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
    left_join(m_bbmo_10y, by = c('sample_id')) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::filter(diversity_index == 'bray_curtis_result') |>
    left_join(euclidean_distance_tb_red) |>
    ggplot(aes(value, euclidean_distance, color = year, shape = season))+
    geom_point()+
    facet_grid(fraction~year)+
    labs(x = 'Bray Curtis Dissimilarity', y = 'Euclidean Environmental Distance', color = 'year')+
    scale_color_manual(values = palette_years)+
    theme_bw()+
    theme(aspect.ratio = 4/4,
          panel.grid.minor = element_blank(),
          axis.title  = element_text(size = 10),
          # axis.text.x = element_text(size = 7), 
          #panel.grid.major.y = element_blank(), 
          strip.background = element_blank(), legend.position = 'none'#,
         # strip.text = element_blank()
          )
  
  unifrac_tibble |>
    left_join(m_red, by = c('sample_id_1' = 'sample_id')) |>
    dplyr::rename(sample_id_num_1 = sample_id_num) |>
    dplyr::select(-date) |>
    left_join(m_red, by = c('sample_id_2' = 'sample_id')) |>
    dplyr::rename(sample_id_num_2 = sample_id_num) |>
    dplyr::mutate(sample_distance = (as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))) |>
    dplyr::filter(sample_distance == 1 &
                    fraction.x == fraction.y) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::filter(!str_detect(date, '2003-')) |>
    left_join(euclidean_distance_tb) |>
    ungroup() |>
    left_join(m_bbmo_10y) |>
    ggplot(aes(wunifrac_distance, euclidean_distance, color = year))+
    geom_point()+
    #facet_wrap(fraction.x~year)+
    scale_color_manual(values = palette_years)+
    labs(x = 'Weigthed UniFrac distance', y = 'Euclidean Environmental Distance', color = 'year')+
    theme_bw()+
    theme(aspect.ratio = 4/4,
          panel.grid.minor = element_blank(),
          axis.title  = element_text(size = 10),
          # axis.text.x = element_text(size = 7), 
          #panel.grid.major = element_blank(), 
          strip.background = element_blank(), legend.position = 'none',
         # strip.text = element_blank()
          )
  
  m_bbmo_10y |>
    colnames()
  
  unifrac_tibble_ed <- unifrac_tibble |>
    left_join(m_red, by = c('sample_id_1' = 'sample_id')) |>
    dplyr::rename(sample_id_num_1 = sample_id_num) |>
    dplyr::select(-date) |>
    left_join(m_red, by = c('sample_id_2' = 'sample_id')) |>
    dplyr::rename(sample_id_num_2 = sample_id_num) |>
    dplyr::mutate(sample_distance = (as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))) |>
    dplyr::filter(sample_distance == 1 &
                    fraction.x == fraction.y) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::filter(!str_detect(date, '2003-')) |>
    dplyr::select(fraction = fraction.x, wunifrac_distance, date)
  
  community_eveness_all |> 
    left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
    dplyr::select(-row_index_2) |>
    pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
    left_join(m_bbmo_10y, by = c('sample_id')) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::filter(diversity_index == 'bray_curtis_result') |> 
    left_join(unifrac_tibble_ed) |>
    ggplot(aes(wunifrac_distance, value))+
    geom_point(aes(color = year, shape = season))+
    labs(x = 'Weighted Unifrac Distance', y = 'Bray Curtis Dissimilarity')+
    scale_color_manual(values = palette_years)+
    #geom_smooth(method = 'loess')+
    facet_grid(fraction~year)+
    theme_bw()+
    theme(aspect.ratio = 4/4,
          panel.grid.minor = element_blank(),
          axis.title  = element_text(size = 10),
          # axis.text.x = element_text(size = 7), 
          #panel.grid.major.y = element_blank(), 
          strip.background = element_blank(), 
          legend.position = 'none',
          #strip.text = element_blank()
          )
  
## At the level of community I would like to detect anomalies with the different approximations ---------
  ### Are they coupled?
  
### Weigthed phylogenetic anomalies -----
wunifrac_anomalies <- unifrac_tibble_ed |>
  group_by(fraction) |>
  dplyr::reframe(#anomalies_ab = get_anomalies(time_lag = 3, negative = FALSE, na_rm = TRUE, cutoff = 1,96, values = pseudoabundance, plotting = FALSE)[c(1,2,3)],
    anomalies_ra = get_anomalies(time_lag = 2, negative = FALSE, 
                                 cutoff = 1.96, na_rm = TRUE, 
                                 values = wunifrac_distance, 
                                 plotting = TRUE)[c(1,2,3)])
  
z_scores_wunifrac_anomalies <- unifrac_tibble_ed |>
  group_by(fraction) |>
  dplyr::reframe(z_score_wu = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = wunifrac_distance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols =z_score_wu)

z_scores_wunifrac_anomalies_02 <-   
  z_scores_wunifrac_anomalies |>
  dplyr::filter(fraction == '0.2') |>
  as_tibble() |>
  unnest(cols = z_score_wu) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')
  
z_scores_wunifrac_anomalies_3 <-   
  z_scores_wunifrac_anomalies |>
  dplyr::filter(fraction == '3') |>
  as_tibble() |>
  unnest(cols = z_score_wu) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

### Euclidean distance anomalies ----
z_environment <- euclidean_distance_tb |>
  #ungroup() %>%
  #group_by(sample_id) %>%
  dplyr::reframe(anomalies_euclidean = get_anomalies(time_lag = 2, values = euclidean_distance, plotting = TRUE)[c(1,2,3)])

m_02_plus_1 <- m_02 |>
  dplyr::mutate(sample_id_num = as.numeric(sample_id_num)+1)

z_scores_environment <- euclidean_distance_tb |>
  dplyr::reframe(z_score_euclidean = get_anomalies(time_lag = 2, values = euclidean_distance, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_euclidean) |>
  dplyr::filter(!is.na(z_score_euclidean)) |>
  dplyr::mutate(sample_id_num = str_c(4:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

## At the level of community, we use the Bray Curtis dissimilarity ----
z_diversity <- bray_curtis_02_rar |>
  dplyr::right_join(community_eveness_02, by = join_by("samples" == "sample_id")) |> 
  #ungroup() %>%
  #group_by(sample_id) %>%
  dplyr::reframe(anomalies_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = TRUE)[c(1,2,3)],# ),
                 anomalies_eveness = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = TRUE)[c(1,2,3)])

## Recover z_scores of diversity----
z_scores_div_02_bray <- bray_curtis_02_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_div_3_bray <- bray_curtis_3_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_bray <- z_scores_div_02_bray |>
  bind_rows(z_scores_div_3_bray)
  
## At the level of community, we use the Evenness result and Bray Curtis dissimilarity ----
z_diversity <- bray_curtis_02_rar |>
  dplyr::right_join(community_eveness_02, by = join_by("samples" == "sample_id")) |> 
  #ungroup() %>%
  #group_by(sample_id) %>%
  dplyr::reframe(anomalies_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = TRUE)[c(1,2,3)],# ),
                 anomalies_eveness = get_anomalies(time_lag = 2, values = community_eveness_rar, plotting = TRUE)[c(1,2,3)])

z_diversity %>%
  str()

## Recover z_scores of diversity----
z_scores_div_02_bray <- bray_curtis_02_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

z_scores_div_3_bray <- bray_curtis_3_rar |>
  dplyr::filter(bray_curtis_result != is.na(bray_curtis_result)) |>
  dplyr::reframe(z_score_bray = get_anomalies(time_lag = 2, values = bray_curtis_result, plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_bray) |>
  dplyr::mutate(sample_id_num = str_c(2:nrow(m_3))) |>
  left_join(m_3, by = 'sample_id_num')

z_scores_bray <- z_scores_div_02_bray |>
  bind_rows(z_scores_div_3_bray)
  
### in these plots we colour the anomaly points ----
z_scores_wunifrac_red <- z_scores_wunifrac_anomalies_02  |>
  dplyr::select(fraction = fraction.x, z_score_wu, sample_id)

z_scores_wunifrac_all <- z_scores_wunifrac_anomalies_3 |>
  dplyr::select(fraction = fraction.x, z_score_wu, sample_id) |>
  bind_rows(z_scores_wunifrac_red)

weigthed_unifrac_plot <- unifrac_tibble |>
  left_join(m_red, by = c('sample_id_1' = 'sample_id')) |>
  dplyr::rename(sample_id_num_1 = sample_id_num) |>
  dplyr::select(-date) |>
  left_join(m_red, by = c('sample_id_2' = 'sample_id')) |>
  dplyr::rename(sample_id_num_2 = sample_id_num) |>
  dplyr::mutate(sample_distance = (as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))) |>
  dplyr::filter(sample_distance == 1 &
                  fraction.x == fraction.y) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(!str_detect(date, '2003-')) |>
  left_join(z_scores_wunifrac_all, by = c('fraction.x' = 'fraction', 'sample_id_1' = 'sample_id')) |>
  dplyr::mutate(anomaly_color = if_else(abs(z_score_wu ) >= 1.96,  '#9F0011', '#080808', missing = '#080808')) |>
  ggplot(aes(date, wunifrac_distance))+
  #facet_wrap(vars(fraction.x))+
  geom_line(aes(date, group = fraction.x, color = fraction.x, linetype = fraction.x))+
  scale_linetype_discrete(labels = labs_fraction)+
  geom_point(aes(alpha = if_else(z_score_wu >= 1.96, 1, 0, missing = 0)))+
  #scale_color_identity()+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = 'Weighted UNIFRAC distance', color = 'Fraction', linetype = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank())

weigthed_unifrac_plot

### composition plot with env euclidean distance, evenness and phylogenetic distance
z_scores_environment_red <- z_scores_environment |>
  dplyr::select(date, z_score_euclidean) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

euclidean_distance_tb |>
  colnames()

euclidean_distance_tb <- euclidean_distance_tb |>
  #dplyr::select(-z_score_euclidean.x, -z_score_euclidean.y) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  left_join(z_scores_environment_red, by = c('date'))

euclidean_distance_plot <- euclidean_distance_tb |>  
  ggplot(aes(date, euclidean_distance))+
  geom_line()+
  # facet_wrap(vars(season), labeller = labs_fraction_s, scales = 'free_x')+
  #scale_y_continuous(expand = c(0,0))+
  geom_point(aes(alpha = if_else(as.numeric(z_score_euclidean) >= 1.96, 1, 0, missing = 0)))+
  #geom_smooth( method = 'loess', span = 0.019)+
  theme_bw()+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(y = 'Environmental Euclidean distance', x = 'Date')+
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        # axis.text.x = element_text(size = 7), 
        panel.grid.major.y = element_blank(), 
        strip.text = element_text(size = 7),
        #axis.text.y = element_text(size = 8),
        #axis.title = element_text(size = 8), 
        strip.background = element_blank(), 
        #  legend.text = element_text(size = 7), 
        # legend.title = element_text(size = 8), 
        axis.title = element_text(size = 10),
        strip.placement = 'outside')+
  guides(shape = 'none',
         color = guide_legend(ncol =1, size = 9,
                              override.aes = aes(label = ''))) 

euclidean_distance_plot

## beta diversity plot with z_scores
z_scores_bray_red <- z_scores_bray |>
  dplyr::select(date, fraction, z_score_bray) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) 

beta_diversity_plot <- 
  community_eveness_all |> 
  left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(diversity_index == 'bray_curtis_result') |>
  left_join(z_scores_bray_red, by = c('fraction', 'date')) |> 
  ggplot(aes(date, value))+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  #geom_point(aes(shape = fraction, color = fraction))+
  geom_line(aes(date, value, group = fraction, color = fraction, linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  facet_grid(vars(diversity_index), labeller = labs_diversity)+
  geom_point(aes(alpha = if_else(z_score_bray >= 1.96, 1, 0, missing = 0)))+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  #geom_smooth(aes(date, value, group = fraction, color = fraction, linetype = fraction), method = 'loess', span = 0.019)+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = 'Bray Curtis Dissimilarity', color = 'Fraction', linetype = 'Fraction')+
  guides(shape = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 10),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(), 
    strip.background = element_blank(), legend.position = 'none',
    strip.text = element_blank())

beta_diversity_plot

beta_diversity_plot_legend <- 
  community_eveness_all |> 
  left_join(bray_curtis_rar_all, by = c('sample_id' = 'samples')) |>
  dplyr::select(-row_index_2) |>
  pivot_longer(cols = c('community_eveness_rar', 'bray_curtis_result'), names_to = 'diversity_index') |>
  left_join(m_bbmo_10y, by = c('sample_id')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(diversity_index == 'bray_curtis_result') |>
  left_join(z_scores_bray_red, by = c('fraction', 'date')) |> 
  ggplot(aes(date, value))+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  #                                                   ymin = -Inf, ymax = Inf), fill = '#94969E', alpha = 0.6)+
  #geom_point(aes(shape = fraction, color = fraction))+
  geom_line(aes(date, value, group = fraction, color = fraction, linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  #facet_wrap(diversity_index~., labeller = labs_diversity)+
  facet_grid(vars(diversity_index), labeller = labs_diversity)+
  geom_point(aes(alpha = if_else(z_score_bray >= 1.96, 1, 0, missing = 0)))+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  #geom_smooth(aes(date, value, group = fraction, color = fraction, linetype = fraction), method = 'loess', span = 0.019)+
  geom_vline(xintercept = as.numeric(as.Date("2005-01-01")), color = '#000000')+
  labs(x = 'Date', y = 'Bray Curtis Dissimilarity', color = 'Fraction', linetype = 'Fraction')+
  guides(shape = 'none', alpha = 'none')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 10),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    strip.text = element_blank())

legend_fraction <- get_legend(beta_diversity_plot_legend)

### Composition plot with anomalies ----

beta_diversity_phylo_env_anomalies_plot <- grid.arrange(beta_diversity_plot, 
                                              weigthed_unifrac_plot,
                                              legend_fraction,
                                              euclidean_distance_plot,
                                              heights=c(1,1, 0.25, 1),
                                              ncol = 1)

# ggsave(beta_diversity_phylo_env_anomalies_plot,
#        filename = 'beta_diversity_phylo_env_plot_v3.pdf',
#               path = 'Results/Figures/',
#               width = 180, height = 160, units = 'mm'
#        )

## Calculate weigthed anomalies in the phylogenetic distance with the community -----
### Are these anomalies the coupled with environmental anomalies in some cases?
### Or are they related with communtiy anamalies

phylogenetic_weighted_tb |>
  colnames()

phylogenetic_weighted_tb_02 <- phylogenetic_weighted_tb |>
  dplyr::select(asv_num, sample_id, date, fraction, phylogenetic_distance) |>
  dplyr::filter(fraction == '0.2') |>
  ungroup() |>
  distinct(asv_num, sample_id, date, fraction, phylogenetic_distance)

z_02_wphy <- phylogenetic_weighted_tb_02 |>
  dplyr::reframe(anomalies_wphy = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1.96, 
                                                values = phylogenetic_distance, plotting = FALSE)[c(1,2,3)])

z_scores_wphy_02 <- phylogenetic_weighted_tb_02 |>
  group_by(asv_num) |>
  dplyr::reframe(z_score_wphy = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                            na_rm = TRUE, values = phylogenetic_distance, 
                                            plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_wphy) |>
  group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num')

phylogenetic_weighted_tb_3 <- phylogenetic_weighted_tb |>
  dplyr::select(asv_num, sample_id, date, fraction, phylogenetic_distance) |>
  dplyr::filter(fraction == '3') |>
  ungroup() |>
  distinct(asv_num, sample_id, date, fraction, phylogenetic_distance)

z_02_wphy <- phylogenetic_weighted_tb_3 |>
  dplyr::reframe(anomalies_wphy = get_anomalies(time_lag = 2, negative = FALSE, na_rm = TRUE, cutoff = 1.96, 
                                                values = phylogenetic_distance, plotting = FALSE)[c(1,2,3)])

z_scores_wphy_3 <- phylogenetic_weighted_tb_02 |>
  group_by(asv_num) |>
  dplyr::reframe(z_score_wphy = get_anomalies(time_lag = 2, negative = FALSE, cutoff = 1.96, 
                                              na_rm = TRUE, values = phylogenetic_distance, 
                                              plotting = FALSE)[c(2)]) |>
  as_tibble() |>
  unnest(cols = z_score_wphy) |>
  group_by(asv_num) |>
  dplyr::mutate(sample_id_num = str_c(1:nrow(m_02))) |>
  left_join(m_02, by = 'sample_id_num') |>
  dplyr::mutate(fraction = '3')

## join datasets from both fractions and plot them 
z_scores_wphy_02_red <- z_scores_wphy_02 |>
  dplyr::select(asv_num, z_score_wphy, fraction, date)

z_scores_wphy_all <- z_scores_wphy_3 |>
  dplyr::select(asv_num, z_score_wphy, fraction, date) |>
  bind_rows(z_scores_wphy_02_red)
  
phylogenetic_weighted_tb_z_scores <- phylogenetic_weighted_tb |>
  left_join(z_scores_wphy_all, by = c('asv_num', 'fraction', 'date')) 
  
## Do bloom events match in anomalies of phylogenetic distance with the community?

box_plot_phylogenetic_distance <- phylogenetic_weighted_tb_z_scores |>
  #left_join(decimal_date_tb) |>
  #left_join(nmds_10y_red) |>
  #dplyr::filter(bloom_event == 'bloom') |>
  ggplot(aes(phylogenetic_distance, interaction(asv_num, family_f)))+
  #geom_point(aes(size = relative_abundance))+
  #facet_grid(vars(interaction(family_f,asv_num_f), cols = 4 ))+
  #geom_smooth(aes(group = fraction, fill = order_f, color = order_f),  method = 'loess', span = 0.09)+
  # geom_point(aes(alpha = if_else(bloom_event == 'bloom', 1, 0.15), shape = fraction, #size = relative_abundance/2, 
  #                color = order_f),
  #            position = position_jitter(width = 0.2))+
  #geom_point()+
  geom_point(aes(alpha = if_else(z_score_wphy >= 1.96, 1, 0.25, missing = 0), shape = fraction, #size = relative_abundance/2, 
                 color = bloom_event),
             position = position_jitter(width = 0.2))+
  #facet_wrap(vars(fraction), labeller = labs_fraction)+
  facet_wrap(fraction~bloom_event)+
#scale_alpha_discrete()+
  #geom_violin(aes(group = asv_num_f), alpha = 0.1)+
  #geom_boxplot(aes(group = order_f, color = order_f), alpha = 0.5)+
  #scale_fill_manual(values = palette_order_assigned_bloo)+
  #scale_color_manual(values = palette_order_assigned_bloo)+
  labs(y= 'Order', x = 'ASV weighted phylogenetic distance with the community', 
       size = 'Relative abundance', 
       alpha = 'Bloom event', 
       fill = 'Order',
       color = 'Order',
       shape = 'Fraction')+
  scale_shape_discrete(labels = labs_fraction)+
  theme_bw()+
  guides(alpha = 'none', 
         size = 'none',
         shape = guide_legend(ncol = 1),
         color = 'none',
         fill = 'none')+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.x = element_text(size = 10),
        # plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(1, 'mm'),
        legend.text = element_text(size = 9),
        legend.key = element_rect(color = 'transparent'),
        panel.border = element_blank())

phylogenetic_weighted_tb_z_scores |>
  ungroup() |>
  distinct(class_f)

phylogenetic_weighted_tb_z_scores |>
  dplyr::filter(class_f == 'Cyanobacteriia') |>
  ggplot(aes(date, phylogenetic_distance))+
  geom_line(aes(group = asv_num))+
  geom_point(aes(alpha = if_else(z_score_wphy >= 1.96, 1, 0.25, missing = 0), shape = fraction, #size = relative_abundance/2, 
                 color = bloom_event),
             position = position_jitter(width = 0.2))+
  facet_wrap(vars(fraction))+
  theme(legend.position = 'none')

box_plot_phylogenetic_distance

## Correlations between Euclidean distance (environmental and Bray Curtis) ----
euclidean_distance_tb_bio <- euclidean_distance_tb_bio |> 
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_bio')

euclidean_distance_tb_phch <- euclidean_distance_tb_phch |> 
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_phch')
  
euclidean_distance_tb_all <- euclidean_distance_tb_g |>
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_g') |>
  bind_rows(euclidean_distance_tb_bio) |>
  bind_rows(euclidean_distance_tb_phch) |>
  pivot_wider(id_cols = ! starts_with('euclidean'), names_from = 'euclidean_distance_type', values_from = 'euclidean_distance') 

bray_curtis_rar_all_m <- bray_curtis_rar_all |>
  dplyr::select(-row_index_2) |>
  # dplyr::mutate(fraction = case_when(str_detect(samples, '0.2_') ~ '0.2',
  #                                    str_detect(samples, '_3_') ~ '3')) |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) 

euc_bray <- euclidean_distance_tb_all |>
  full_join(bray_curtis_rar_all_m, by = 'decimal_date') 

bloo_all_types_summary_ed <- bloo_all_types_summary |>
  dplyr::mutate(fraction = as.character(fraction)) |>
  distinct(asv_num, frequency, fraction)

bloom_event_type_dec_date <- bloom_event |>
  dplyr::mutate(bloom_event = case_when(asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'no-bloom',
                                        !asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ bloom_event)) |>
  left_join(m_bbmo_10y) |>
  dplyr::select(asv_num, bloom_event, decimal_date, fraction) |>
  dplyr::filter(bloom_event == 'bloom') |>
  left_join(bloo_all_types_summary_ed)

euc_bray_bloom_event <- euc_bray |>
  left_join(bloom_event_dec_date)

euc_bray_bloom_event  |>
  ggplot(aes(bray_curtis_result, euclidean_distance_g, shape = fraction, color = fraction))+
  geom_point(aes())+
  facet_wrap(bloom_event~frequency)+
  geom_smooth(method = 'lm')

euc_bray_bloom_event |>
  ggplot(aes(bray_curtis_result, euclidean_distance_phch, shape = fraction, color = fraction))+
  geom_point(aes())+
  facet_wrap(bloom_event~frequency)+
  geom_smooth(method = 'lm')

euc_bray_bloom_event  |>
  ggplot(aes(bray_curtis_result, euclidean_distance_bio, shape = fraction, color = fraction))+
  geom_point(aes())+
  facet_wrap(bloom_event~frequency)+
  geom_smooth(method = 'loess')

## type of bloom ----
euclidean_distance_tb_bio <- euclidean_distance_tb_bio |> 
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_bio')

euclidean_distance_tb_phch <- euclidean_distance_tb_phch |> 
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_phch')

euclidean_distance_tb_all <- euclidean_distance_tb_g |>
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_g') |>
  bind_rows(euclidean_distance_tb_bio) |>
  bind_rows(euclidean_distance_tb_phch) |>
  pivot_wider(id_cols = ! starts_with('euclidean'), names_from = 'euclidean_distance_type', values_from = 'euclidean_distance') |>
  dplyr::mutate(across(starts_with('euclidean'), scale))

bray_curtis_rar_all_m <- bray_curtis_rar_all |>
  dplyr::select(-row_index_2) |>
  # dplyr::mutate(fraction = case_when(str_detect(samples, '0.2_') ~ '0.2',
  #                                    str_detect(samples, '_3_') ~ '3')) |>
  left_join(m_bbmo_10y, by = c('samples' = 'sample_id')) 

euc_bray <- euclidean_distance_tb_all |>
  full_join(bray_curtis_rar_all_m, by = 'decimal_date') 

bloo_all_types_summary_ed <- bloo_all_types_summary |>
  dplyr::mutate(fraction = as.character(fraction)) |>
  distinct(asv_num, frequency, fraction)

bloom_event_type_dec_date <- bloom_event |>
  dplyr::mutate(bloom_event = case_when(asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'no-bloom',
                                        !asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ bloom_event)) |>
  left_join(m_bbmo_10y) |>
  dplyr::select(asv_num, bloom_event, decimal_date, fraction) |>
  dplyr::filter(bloom_event == 'bloom') |>
  left_join(bloo_all_types_summary_ed) |>
  distinct(decimal_date, bloom_event, fraction, frequency)

euc_bray_bloom_event <- euc_bray |>
  dplyr::select(starts_with('euclidean'), bray_curtis_result, decimal_date, fraction) |>
  left_join(bloom_event_type_dec_date, by = c('decimal_date', 'fraction')) |>
  dplyr::mutate(bloom_event = case_when(is.na(bloom_event) ~ 'no-bloom',
                                        !is.na(bloom_event) ~ bloom_event), 
                frequency = case_when(is.na(frequency) ~ 'no-bloom',
                                      !is.na(frequency) ~ frequency))

labs_frequency <- as_labeller(c('no-bloom' = 'No Bloom Event', 'seasonal' = 'Seasonal' , 'stochastic' = 'Chaotic' ))

euclidean_g_bray_bloom_plot <- euc_bray_bloom_event  |>
  dplyr::filter(!is.na(fraction)) |>
 # dplyr::filter(euclidean_distance_g > 0) |>
  ggplot(aes(bray_curtis_result, euclidean_distance_g))+
  geom_point(aes(shape = fraction, color = fraction))+
  facet_wrap(vars(frequency), labeller = labs_frequency)+
  geom_smooth(method = 'loess', aes(group = fraction, color = fraction, fill = fraction))+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_fill_manual(values= palette_fraction, labels = labs_fraction)+
  scale_shape_discrete( labels = labs_fraction)+
  #scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Bray Curtis Dissimilarity', y = 'Euclidean Environmental Distance Scaled', 
       color = 'Fraction', 
       linetype = 'Fraction', shape = 'Fraction', fill = 'Fraction')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank())

euclidean_g_bray_bloom_plot

ggsave(euclidean_g_bray_bloom_plot,
       filename = 'euclidean_g_bray_bloom_sm_plot.pdf',
              path = 'Results/Figures/',
              width = 180, height = 160, units = 'mm'
       )

euc_bray_bloom_event  |>
  dplyr::filter(!is.na(fraction)) |>
  ggplot(aes(bray_curtis_result, euclidean_distance_g))+
  geom_point(aes(shape = fraction, color = fraction))+
  facet_wrap(vars(frequency))+
  geom_smooth(method = 'lm', aes(group = fraction, color = fraction))+
  #cale_color_manual(values= palette_fraction, labels = labs_fraction)+
  #scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Bray Curtis Dissimilarity', y = 'Euclidean Environmental Distance', color = 'Fraction', linetype = 'Fraction')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank())

euclidean_phch_bray_bloom_plot <- euc_bray_bloom_event |>
  dplyr::filter(!is.na(fraction)) |>
  # dplyr::filter(euclidean_distance_g > 0) |>
  ggplot(aes(bray_curtis_result, euclidean_distance_phch))+
  geom_point(aes(shape = fraction, color = fraction))+
  facet_wrap(vars(frequency), labeller = labs_frequency)+
  geom_smooth(method = 'lm', aes(group = fraction, color = fraction, fill = fraction))+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_fill_manual(values= palette_fraction, labels = labs_fraction)+
  scale_shape_discrete( labels = labs_fraction)+
  #scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Bray Curtis Dissimilarity', y = 'Euclidean Environmental Distance Scaled', 
       color = 'Fraction', 
       linetype = 'Fraction', shape = 'Fraction', fill = 'Fraction')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank())

euclidean_bio_bray_bloom_plot

ggsave(euclidean_bio_bray_bloom_plot,
       filename = 'euclidean_phch_bray_bloom_lm_plot.pdf',
       path = 'Results/Figures/',
       width = 180, height = 160, units = 'mm'
)

euc_bray_bloom_event  |>
  ggplot(aes(bray_curtis_result, euclidean_distance_g, shape = fraction, color = fraction))+
  geom_point(aes())+
  facet_wrap(bloom_event~frequency)+
  geom_smooth(method = 'lm')

euclidean_bio_bray_bloom_plot <- euc_bray_bloom_event |>
  dplyr::filter(!is.na(fraction)) |>
  # dplyr::filter(euclidean_distance_g > 0) |>
  ggplot(aes(bray_curtis_result, euclidean_distance_bio))+
  geom_point(aes(shape = fraction, color = fraction))+
  facet_wrap(vars(frequency), labeller = labs_frequency)+
  geom_smooth(method = 'loess', aes(group = fraction, color = fraction, fill = fraction))+
  scale_color_manual(values= palette_fraction, labels = labs_fraction)+
  scale_fill_manual(values= palette_fraction, labels = labs_fraction)+
  scale_shape_discrete( labels = labs_fraction)+
  #scale_x_datetime(date_breaks = 'year', date_labels = '%Y')+
  labs(x = 'Bray Curtis Dissimilarity', y = 'Euclidean Environmental Distance Scaled', 
       color = 'Fraction', 
       linetype = 'Fraction', shape = 'Fraction', fill = 'Fraction')+
  #scale_y_continuous(limits = c(0.15, 1.0))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank())

## I add a lag, expecting that when there is an abrupt change in env conditions then we see an abrupt change in the community ----
euclidean_distance_tb_bio <- euclidean_distance_tb_bio |> 
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_bio')

euclidean_distance_tb_phch <- euclidean_distance_tb_phch |> 
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_phch')

euclidean_distance_tb_all <- euclidean_distance_tb_g |>
  dplyr::mutate(euclidean_distance_type = 'euclidean_distance_g') |>
  bind_rows(euclidean_distance_tb_bio) |>
  bind_rows(euclidean_distance_tb_phch) |>
  pivot_wider(id_cols = ! starts_with('euclidean'), names_from = 'euclidean_distance_type', values_from = 'euclidean_distance') 

m_ed <- m_02 |>
  bind_rows(m_3) |>
  dplyr::select(fraction, decimal_date, sample_id_num)

bray_curtis_rar_all_m <- bray_curtis_rar_all |>
  dplyr::mutate(row_index_lag = (row_index_2 + 1)) |>
  dplyr::mutate(row_index_lag = as.character(row_index_lag)) |>
  dplyr::mutate(fraction = case_when(str_detect(samples, '0.2_') ~ '0.2',
                                     str_detect(samples, '_3_') ~ '3')) |>
  dplyr::select(-row_index_2, -samples) |>
  left_join(m_ed, by = c('row_index_lag' = 'sample_id_num', 'fraction')) 

euc_bray <- euclidean_distance_tb_all |>
  full_join(bray_curtis_rar_all_m, by = 'decimal_date') 

bloo_all_types_summary_ed <- bloo_all_types_summary |>
  dplyr::mutate(fraction = as.character(fraction)) |>
  distinct(asv_num, frequency, fraction)

bloom_event_type_dec_date <- bloom_event |>
  dplyr::mutate(bloom_event = case_when(asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'no-bloom',
                                        !asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ bloom_event)) |>
  left_join(m_bbmo_10y) |>
  dplyr::select(asv_num, bloom_event, decimal_date, fraction) |>
  dplyr::filter(bloom_event == 'bloom') |>
  left_join(bloo_all_types_summary_ed) |>
  distinct(decimal_date, bloom_event, fraction, frequency)

euc_bray_bloom_event <- euc_bray |>
  dplyr::select(starts_with('euclidean'), bray_curtis_result, decimal_date, fraction) |>
  left_join(bloom_event_type_dec_date, by = c('decimal_date', 'fraction')) |>
  dplyr::mutate(bloom_event = case_when(is.na(bloom_event) ~ 'no-bloom',
                                        !is.na(bloom_event) ~ bloom_event), 
                frequency = case_when(is.na(frequency) ~ 'no-bloom',
                                      !is.na(frequency) ~ frequency))

euc_bray_bloom_event  |>
  ggplot(aes(bray_curtis_result, euclidean_distance_g, shape = fraction, color = fraction))+
  geom_point(aes())+
  facet_wrap(bloom_event~frequency)+
  geom_smooth(method = 'loess')

euc_bray_bloom_event |>
  ggplot(aes(bray_curtis_result, euclidean_distance_phch, shape = fraction, color = fraction))+
  geom_point(aes())+
  facet_wrap(bloom_event~frequency)+
  geom_smooth(method = 'loess')

euc_bray_bloom_event  |>
  ggplot(aes(bray_curtis_result, euclidean_distance_bio, shape = fraction, color = fraction))+
  geom_point(aes())+
  facet_wrap(bloom_event~frequency)+
  geom_smooth(method = 'loess')





# Another idea. I would like to observe if the relationship between taxa and their closely related taxa -----
## evolves in some way between them 

distances_asv7 <- phylogenetic_distances_tb_com |>
  dplyr::filter(asv_num_1 == 'asv7') |>
  dplyr::filter(asv_num_2 != 'asv7')

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num == 'asv7') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(fraction == '0.2')

asv_tab_10y_02_rel |>
  full_join(distances_asv7, by = c('asv_num' = 'asv_num_2')) |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  dplyr::mutate(weighted_distance = relative_abundance*phylogenetic_distance) |>
  ggplot(aes(sample_id, weighted_distance))+
  facet_wrap(vars(asv_num))+
  geom_line(aes(group = asv_num))

asv_tab_10y_02_rel |>
  dplyr::filter(asv_num %in% c('asv109', 'asv7')) |>
  ggplot(aes(sample_id, relative_abundance))+
  #facet_wrap(vars(asv_num))+
  geom_line(aes(group = asv_num, color = asv_num))+
  geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09)+
  scale_color_manual(values = palette_long)+
  scale_fill_manual(values = palette_long)+
  theme_bw()

##
distances_asv22 <- phylogenetic_distances_tb_com |>
  dplyr::filter(asv_num_1 == 'asv22') |>
  dplyr::filter(asv_num_2 != 'asv22')

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num == 'asv22') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::filter(fraction == '0.2')

asv_tab_10y_02_rel |>
  full_join(distances_asv22, by = c('asv_num' = 'asv_num_2')) |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  dplyr::mutate(weighted_distance = relative_abundance*phylogenetic_distance) |>
  ggplot(aes(sample_id, weighted_distance))+
  facet_wrap(vars(asv_num))+
  geom_line(aes(group = asv_num))

asv_tab_10y_02_rel |>
  dplyr::filter(asv_num %in% c('asv47', 'asv22', 'asv73', 'asv757')) |>
  ggplot(aes(sample_id, relative_abundance))+
  #facet_wrap(vars(asv_num))+
  geom_line(aes(group = asv_num, color = asv_num))+
  geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09)+
  scale_color_manual(values = palette_long)+
  scale_fill_manual(values = palette_long)+
  theme_bw()

##
close_distances_asv80 <- phylogenetic_distances_tb_com |>
  dplyr::filter(asv_num_1 == 'asv80') |>
  dplyr::filter(asv_num_2 != 'asv80') |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)
# 
# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(asv_num == 'asv80') |>
#   dplyr::filter(abundance_type == 'rclr') |>
#   dplyr::filter(fraction == '0.2')
# 
# asv_tab_10y_02_rel |>
#   full_join(distances_asv80, by = c('asv_num' = 'asv_num_2')) |>
#   dplyr::filter(phylogenetic_distance < 0.012) |>
#   dplyr::mutate(weighted_distance = relative_abundance*phylogenetic_distance) |>
#   ggplot(aes(sample_id, weighted_distance))+
#   facet_wrap(vars(asv_num))+
#   geom_line(aes(group = asv_num))

asv_tab_10y_02_rel |>
  dplyr::filter(asv_num %in% close_distances_asv80) |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, relative_abundance))+
  #facet_wrap(vars(asv_num))+
  labs(x = 'Time', y = 'Relative abundance', color = 'ASV', fill = 'ASV')+
  geom_line(aes(group = asv_num, color = asv_num))+
  geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09)+
  scale_color_manual(values = palette_long)+
  scale_fill_manual(values = palette_long)+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    strip.background = element_blank(), legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    # axis.text.x = element_text(size = 7), 
    panel.grid.major.y = element_blank())

## I loop to observe how what happens with all my bloomers ----
# Define the vector of ASV numbers
bloo_02$value

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_02$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_02_rel |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_02) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, relative_abundance)) +
    labs(x = 'Time', y = 'Relative abundance', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num)) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_long) +
    scale_fill_manual(values = palette_long) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

# Print the plots
for (asv_num in bloo_02$value) {
  print(results[[asv_num]])
}

## the same with rCLR data----
## I loop to observe how what happens with all my bloomers ----
# Define the vector of ASV numbers
bloo_02$value

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_02$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_02_rclr |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_02) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, rclr)) +
    labs(x = 'Time', y = 'rCLR', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num)) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_long) +
    scale_fill_manual(values = palette_long) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_rclr.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

# Print the plots
for (asv_num in bloo_02$value) {
  print(results[[asv_num]])
}

## the same for PA bloomers ----
## I loop to observe how what happens with all my bloomers ----
# Define the vector of ASV numbers
bloo_3$value

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_3$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_3_rel |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_3) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, relative_abundance)) +
    labs(x = 'Time', y = 'Relative abundance', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num)) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_long) +
    scale_fill_manual(values = palette_long) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_pa.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

# Print the plots
for (asv_num in bloo_3$value) {
  print(results[[asv_num]])
}

## the same with rCLR data----
## I loop to observe how what happens with all my bloomers ----
# Define the vector of ASV numbers
bloo_3$value

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_3$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_3_rclr |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_3) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, rclr)) +
    labs(x = 'Time', y = 'rCLR', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num)) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_long) +
    scale_fill_manual(values = palette_long) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_pa_rclr.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

# Print the plots
for (asv_num in bloo_3$value) {
  print(results[[asv_num]])
}


# I would like to highlight the BLOOMERS in the plots and observe what happens -----
## I loop to observe how what happens with all my bloomers ----
# Define the vector of ASV numbers
bloo_02_filt <- bloo_02 |>
  dplyr::filter(!value %in% c('asv5', 'asv8', 'asv2', 'asv3'))

bloo_02_filt$value

asv_tab_joined |>
  colnames()

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_02_filt$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_02_rel |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_02) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                                      !asv_num %in% bloo_taxonomy$asv_num ~ 'no-bloomer'))
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, relative_abundance)) +
    labs(x = 'Time', y = 'Relative abundance', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_genes) +
    scale_fill_manual(values = palette_genes) +
    theme_bw() +
    guides(alpha = 'none')+
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_bloo.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

# Print the plots
for (asv_num in bloo_02$value) {
  print(results[[asv_num]])
}

## the same with rCLR data----
## I loop to observe how what happens with all my bloomers ----
# Define the vector of ASV numbers
bloo_02$value

# Create an empty list to store the results
results <- list()

bloo_taxonomy <- bloo_taxonomy |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8'))

# Loop over the ASV numbers
for (asv_num in bloo_02_filt$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_02_rclr |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_02) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                                      !asv_num %in% bloo_taxonomy$asv_num ~ 'no-bloomer'))
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, rclr)) +
    labs(x = 'Time', y = 'rCLR', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_genes) +
    scale_fill_manual(values = palette_genes) +
    theme_bw() +
    guides(alpha = 'none')+
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_rclr_ed.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

# Print the plots
for (asv_num in bloo_02$value) {
  print(results[[asv_num]])
}

## the same for PA bloomers ----
## I loop to observe how what happens with all my bloomers ----
# Define the vector of ASV numbers
bloo_3$value

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_3$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_3_rel |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_3) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                                      !asv_num %in% bloo_taxonomy$asv_num ~ 'no-bloomer'))
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, relative_abundance)) +
    labs(x = 'Time', y = 'Relative abundance', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_genes) +
    scale_fill_manual(values = palette_genes) +
    theme_bw() +
    guides(alpha = 'none')+
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_pa_ed.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

# Print the plots
for (asv_num in bloo_3$value) {
  print(results[[asv_num]])
}

## the same with rCLR data----
## I loop to observe how what happens with all my bloomers ----
# Define the vector of ASV numbers
bloo_3$value

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_3$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_3_rclr |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_3) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                                      !asv_num %in% bloo_taxonomy$asv_num ~ 'no-bloomer'))
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, rclr)) +
    labs(x = 'Time', y = 'rCLR', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_genes) +
    scale_fill_manual(values = palette_genes) +
    theme_bw() +
    guides(alpha = 'none')+
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_pa_rclr_ed.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

# Print the plots
for (asv_num in bloo_3$value) {
  print(results[[asv_num]])
}

# I would like to select those closely related taxa that have two bloomers (for the main text) ----

## FL ----
# Define the vector of ASV numbers
bloo_02_filt$value

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_02_filt$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_02_rclr |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_02) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                                      !asv_num %in% bloo_taxonomy$asv_num ~ 'no-bloomer'))
  
  # Filter out only bloomers and check for at least two distinct ASV with 'bloomer' status
  asv_tab_joined_bloomers <- asv_tab_joined |> 
    dplyr::filter(bloomer == "bloomer") |> 
    dplyr::distinct(asv_num, bloomer)
  
  # Check if there are two or more distinct bloomers
  if (nrow(asv_tab_joined_bloomers) < 2) {
    next  # Skip this iteration if fewer than two bloomers are present
  }
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, rclr)) +
    labs(x = 'Time', y = 'rCLR', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_genes) +
    scale_fill_manual(values = palette_genes) +
    theme_bw() +
    guides(alpha = 'none')+
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_fl_rclr_ed1.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}

## PA ----
# Define the vector of ASV numbers
bloo_3$value

# Create an empty list to store the results
results <- list()

# Loop over the ASV numbers
for (asv_num in bloo_3$value) {
  # Filter the data for the current ASV number
  close_distances <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2 != asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value)
  
  # Filter the data for the current ASV number
  asv_tab_filtered <- asv_tab_10y_3_rclr |>
    dplyr::filter(asv_num %in% close_distances)
  
  # Join the data with the metadata
  asv_tab_joined <- asv_tab_filtered |>
    left_join(m_3) |>
    dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
    dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                                      !asv_num %in% bloo_taxonomy$asv_num ~ 'no-bloomer'))
  
  # Filter out only bloomers and check for at least two distinct ASV with 'bloomer' status
  asv_tab_joined_bloomers <- asv_tab_joined |> 
    dplyr::filter(bloomer == "bloomer") |> 
    dplyr::distinct(asv_num, bloomer)
  
  # Check if there are two or more distinct bloomers
  if (nrow(asv_tab_joined_bloomers) < 2) {
    next  # Skip this iteration if fewer than two bloomers are present
  }
  
  # Create the plot
  plot <- ggplot(asv_tab_joined, aes(date, rclr)) +
    labs(x = 'Time', y = 'rCLR', color = 'ASV', fill = 'ASV') +
    geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
    geom_smooth(aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
    scale_color_manual(values = palette_genes) +
    scale_fill_manual(values = palette_genes) +
    theme_bw() +
    guides(alpha = 'none')+
    theme(
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  # Add the plot to the results list
  results[[asv_num]] <- plot
  
  # Save the plot to a file
  ggsave(paste0(asv_num, "closely_related_taxa_pa_rclr_ed1.pdf"), plot = plot, 
         path = 'Results/Figures/',
         width = 180, height = 160, units = 'mm')
}


## I design a nice plot for the groups that share two bloomers -----
## ASV 4 and ASV31 -----
palette_asv4_asv31 <- c("#de93a6",
                        "#6673a9",
                        "#67aa6c",
                        "#2a0000", 
                        "#9c3354")

close_distances <- phylogenetic_distances_tb_com |>
  dplyr::filter(asv_num_1 == 'asv4') |>
  dplyr::filter(asv_num_2 != 'asv4') |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

## detailed bloomers taxonomy ----
# detailed_bloo_tax <- tax_bbmo_10y_new |>
#   dplyr::filter(asv_num %in% bloo_taxonomy$asv_num_f)
# 
# detailed_bloo_tax |>
#   distinct(family)

# write.csv(detailed_bloo_tax,
#           'results/tables/detailed_bloo_tax.csv', row.names = F)

# Join the data with the metadata 
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

plot_asv31_asv4 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Time', y = 'rCLR',color = 'Synechococcus\nCC9902\nASV', fill = 'Synechococcus\nCC9902\nASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.5) +
  scale_color_manual(values = palette_asv4_asv31) +
  scale_fill_manual(values =palette_asv4_asv31) +
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1)+
  theme_bw() +
  guides(alpha = 'none',
         color = guide_legend(nrow = 3),
         fill = guide_legend(nrow = 3))+
  theme(
    strip.background = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 4),
    panel.grid.major.y = element_blank()
  )

legend_asv31_asv4 <- get_legend(plot_asv31_asv4)

# Create the plot
plot_asv31_asv4 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Time', y = 'rCLR', color = 'ASV', fill = 'ASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.5) +
  scale_color_manual(values = palette_asv4_asv31) +
  scale_fill_manual(values =palette_asv4_asv31) +
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1)+
  theme_bw() +
  guides(alpha = 'none')+
  theme(
    strip.background = element_blank(),
    legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    panel.grid.major.y = element_blank()
  )

plot_asv31_asv4

## ASV 17 and ASV77 ------
palette_asv17_asv77 <- c('asv17' = "#6a62b3",
                         "asv955" = "#5c6873",
                         "asv557" = "#1c1a33",
                         "asv77"  = "#f2c549",
                         "asv2929" = "#6bb362",
                         "asv1185" = "#1a2e79",
                         "asv2300" = "#b36a62",
                         "asv876" = "#a8001d",
                         "asv3052" = "#93d8e5")

close_distances <- phylogenetic_distances_tb_com |>
  dplyr::filter(asv_num_1 == 'asv17') |>
  dplyr::filter(asv_num_2 != 'as17') |>
  dplyr::filter(phylogenetic_distance < 0.012) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

# Join the data with the metadata
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

#create the legend 
plot_asv17_asv77 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Time', y = 'rCLR', color = 'Sphingomonadaceae\nASV', fill = 'Sphingomonadaceae\nASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.5) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
  scale_color_manual(values =  palette_asv17_asv77) +
  scale_fill_manual(values = palette_asv17_asv77) +
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1)+
  theme_bw() +
  guides(alpha = 'none',
         color = guide_legend(nrow = 3),
         fill = guide_legend(nrow = 3))+
  theme(
    strip.background = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    legend.title  = element_text(size = 4),
    legend.text  = element_text(size = 4),
    panel.grid.major.y = element_blank()
  )

legend_asv17_asv77 <- get_legend(plot_asv17_asv77)

# Create the plot
plot_asv17_asv77 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Time', y = 'rCLR', color = 'ASV', fill = 'ASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.5) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
  scale_color_manual(values =  palette_asv17_asv77) +
  scale_fill_manual(values = palette_asv17_asv77) +
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1)+
  theme_bw() +
  guides(alpha = 'none')+
  theme(
    strip.background = element_blank(),
    legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    panel.grid.major.y = element_blank()
  )

plot_asv17_asv77 

closely_related_bloomers_plot <- plot_grid(
  plot_asv31_asv4,
  plot_asv17_asv77,
  #legend_asv31_asv4,
  #legend_asv17_asv77,
  tree_plot_asv17_asv77,
  ncol = 2,                # One column layout for the main grid
  rel_heights = c(2, 2, 1, 1),
  labels = c('A', 'B'), label_fontface = 'plain'
)

# Print the final plot
print(closely_related_bloomers_plot)

# ggsave( plot = closely_related_bloomers_plot,
#         filename = 'closely_related_bloomers_plot.pdf',
#         path = 'results/figures/',
#         width = 180, height = 180, units = 'mm')

## ------- Different phloygenetic distance based on hamming distances  ## ----------

### Remember that the taxonomy in the phyloseq object is the old one.


# Prune the phyloseq object by keeping only the ASVs present in `seqs$asv_num`
tax <- as(tax_table(bbmo_10y@tax_table), 'matrix') %>% 
  data.frame() %>% 
  rownames_to_column(var = 'asv') %>% 
  mutate_if(.predicate = is.factor,.funs = as.character)

comparison_hamm <- function(chr){
  chr %>%
    cross2(., ., .filter = `==`) %>% 
    map(setNames, c("seq1", "seq2")) %>% 
    bind_rows() %>% 
    mutate( hammdist = stringdist(seq1, seq2, method = "hamming"))
}

## update taxonomy 
tax |>
  colnames()

tax <- tax |>
  dplyr::select( seq) |>
  left_join(tax_bbmo_10y_new, by = c( 'seq'))

## family hamm ----
tax |>
  dim() # 7849

# tax |>
#   dplyr::mutate(length = nchar(seq)) |>
#   ggplot(aes(length, y = Family))+
#   geom_density_ridges()+
#   theme_ridges(center_axis_labels = T)
  
family.hamm <- tax %>% 
  dplyr::filter(family %in% unique(bloo_taxonomy$family_f)) |>
  group_by(family) %>% 
  summarise(n = n()) %>% 
  filter( n > 1) %>% 
  pull(family)

hammdist.all <- tax %>%
  dplyr::filter(family %in% family.hamm) %>%
  split(.$family) %>%
  map('seq') %>%
  map(~comparison_hamm(.)) 

hammdist.all |>
  names() ## 13
 
## calculate p-distance and filter for those that have a closely related bloomer ----

tax <- asv_tab_all_bloo_z_tax |>
  dplyr::distinct(asv_num, order_f, family_f, class_f, seq)

closely_related_bloomers  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  mutate(length_1 = nchar(seq1),
         length_2 = nchar(seq2)) %>%
  mutate(p_distance = hammdist/length_1) %>%
  left_join(tax, by = c('seq1' = 'seq')) %>%
  left_join(tax, by = c('seq2' = 'seq')) %>%
  dplyr::filter(p_distance < 0.012) %>%
  distinct(asv_num.x, p_distance, asv_num.y) |>
  dplyr::select(asv_num = asv_num.x, p_distance, asv_num.y) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(asv_num %in% c('asv1', 'asv7', 'asv4', 'asv31', 'asv42', 'asv49', 'asv200', 'asv15', 'asv264',
                                'asv17', 'asv77')) # I remove those that were anlyzed before 

closely_related_bloomers_vs_bloo  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  mutate(length_1 = nchar(seq1),
         length_2 = nchar(seq2)) %>%
  mutate(p_distance = hammdist/length_1) %>%
  left_join(tax, by = c('seq1' = 'seq')) %>%
  left_join(tax, by = c('seq2' = 'seq')) %>%
  dplyr::filter(p_distance < 0.012) %>%
  distinct(asv_num.x, p_distance, asv_num.y) |>
  dplyr::select(asv_num = asv_num.x, p_distance, asv_num.y) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(asv_num %in% bloo_taxonomy$asv_num_f) |> # I remove those that were anlyzed before 
  dplyr::filter(asv_num.y %in% bloo_taxonomy$asv_num_f)
  
## 
closely_related_families <- closely_related_bloomers |>
  distinct(family) # 28

closely_related_bloomers_and_others  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  dplyr::mutate(length_1 = nchar(seq1),
         length_2 = nchar(seq2)) %>%
  dplyr::mutate(p_distance = hammdist/length_1) %>%
  left_join(tax_bbmo_10y_new, by = c('seq1' = 'seq')) %>%
  left_join(tax_bbmo_10y_new, by = c('seq2' = 'seq')) %>%
  filter(p_distance < 0.012) %>%
  distinct(asv_num.x, asv_num.y, p_distance) |>
  dplyr::select(asv_num = asv_num.x, asv_num.y, p_distance) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(!asv_num %in% c('asv1', 'asv7', 'asv4', 'asv31', 'asv42', 'asv49', 'asv200', 'asv15', 'asv264',
                                'asv17', 'asv77')) # I remove those that were anlyzed before 

## now i create the groups and plot it like I did before 
closely_related_bloomers_and_others |>
  dplyr::filter(asv_num %in% bloo_taxonomy$asv_num_f) |>
  distinct(asv_num)  # n = 11 ASVs 
# # A tibble: 11 × 1
# asv_num
# 1 asv77  
# 2 asv17  
# 3 asv4   
# 4 asv31  
# 5 asv7   
# 6 asv1   
# 7 asv42  
# 8 asv49  
# 9 asv200 
# 10 asv15  
# 11 asv264 

# are they all seasonal??? No
closely_related_bloomers_and_others |>
  dplyr::filter(asv_num %in% bloo_taxonomy$asv_num_f) |>
  distinct(asv_num) |>
  left_join(bloo_all_types_summary_tb)

### PLOT the groups of closely related taxa based on p-distances not the tree ------
###  Those clusters that have two bloomers 
##### Synechococcus has two groups one ASV4-ASV31 and ASV1-ASV7. I will plot them separately but construct one only phylogenetic tree -----
closely_related_bloomers_and_others  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  dplyr::mutate(length_1 = nchar(seq1),
                length_2 = nchar(seq2)) %>%
  dplyr::mutate(p_distance = hammdist/length_1) %>%
  left_join(tax_bbmo_10y_new, by = c('seq1' = 'seq')) %>%
  left_join(tax_bbmo_10y_new, by = c('seq2' = 'seq')) %>%
  filter(p_distance < 0.012) %>%
  distinct(asv_num.x, asv_num.y, p_distance) |>
  dplyr::select(asv_num = asv_num.x, asv_num.y, p_distance) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(asv_num %in% c('asv1', 'asv7', 'asv4', 'asv31')) 

closely_related_bloomers_and_others |>
  distinct(asv_num.y) |>
  as_vector()

# create a beautiful palette for Synechococcus closely related ---- 
palette_cianobac <- c('asv5934' = "#de93a6",
                      'asv1337' =  "#6673a9",
                      'asv34' =  "#67aa6c",
                      "asv31"=   "#2a0000", 
                      "asv4"  = "#9c3354", 
                      # second cluster of synechococcus
                      'asv1' = '#E6CF40',
                      "asv7" = '#605991',
                      "asv974" = '#57806B',
                      'asv92' = '#E55440',
                      "asv109" = '#332A32',
                      "asv3533" =  '#805B57',
                      "asv343" = '#FFB6DA',
                      "asv14"  ='#2D4A2F',
                      "asv3073" = '#8A6276',
                      "asv5327" = '#2B3147',
                      "asv5283" = '#8DA0EB',
                      "asv4780"   = '#86AAB0',
                      "asv8058" ='#A6A6A6',
                      "asv90"    = '#FFFBBC',
                      "asv10"  =  '#CE0000' ,  
                      "asv698" =    '#C0B6F5')

## closely related bloomers asv4 and asv31 ----
closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv4') |>
  arrange(p_distance) 

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv31') |>
  arrange(p_distance) 

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv7') |>
  arrange(p_distance)

close_distances <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv4') |>
  slice_min(n = 3, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

# Join the data with the metadata
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

asv_tab_filtered |>
  distinct(asv_num)

asv_tab_filtered$asv_num <- factor(asv_tab_filtered$asv_num, levels = c('asv31', 'asv4', 'asv34', 'asv5834'))

plot_asv31_asv4 <- ggplot(asv_tab_filtered, aes(date, rclr)) + 
  labs(x = 'Time', y = 'rCLR', color = 'Synechococcus\nCC9902\nASV', fill = 'Synechococcus\nCC9902\nASV') + 
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8, 0.6))) + 
  # Smooth for 'bloomer'
  geom_smooth(data = asv_tab_filtered |> filter(bloomer == 'bloomer'), 
              aes(group = asv_num, fill = asv_num), 
              method = 'loess', color = NA, 
              alpha = 0.3, span = 0.09) + 
  # Smooth for non-bloomer
  geom_smooth(data = asv_tab_filtered |> filter(bloomer != 'bloomer'),
              aes(group = asv_num, fill = asv_num),
              method = 'loess', color = NA,
              alpha = 0.5, span = 0.09) +
  scale_color_manual(values = c("asv31" = "#2a0000",
                                "asv4" = "#9c3354",
                                'asv5934' = '#C2AFB3',
                                'asv34' = '#C2AFB3'),
                     labels = c("asv31",
                                "asv4",
                                'asv34, asv5934')) +
  scale_fill_manual(values = c("asv31" = "#2a0000",
                               "asv4" = "#9c3354",
                               'asv5934' = '#C2AFB3',
                               'asv34' = '#C2AFB3'),
                    labels = c("asv31",
                               "asv4",
                               'asv34, asv5934'))+
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1) + 
  theme_bw() + 
  guides(alpha = 'none', 
         color = guide_legend(ncol = 5), 
         fill = guide_legend(ncol = 5)) + 
  theme(strip.background = element_blank(), 
        legend.position = 'bottom', 
        panel.grid.minor = element_blank(), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), 
        panel.grid.major.y = element_blank())

legend_asv31_asv4 <- get_legend(plot_asv31_asv4)

# Create the plot
# plot_asv31_asv4 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
#   labs(x = 'Date', y = 'rCLR', color = 'ASV', fill = 'ASV') +
#   geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
#   geom_smooth(data = asv_tab_filtered |>
#                 dplyr::filter(bloomer != 'bloomer'),
#               aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.3) +
#   geom_smooth(data = asv_tab_filtered |>
#                 dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num),  linewidth = 0.5, span = 0.09) +
#   
#   scale_color_manual(values = palette_cianobac) +
#   scale_fill_manual(values =palette_cianobac) +
#   facet_grid(vars(fraction), labeller = labs_fraction)+
#   theme_bw() +
#   guides(alpha = 'none')+
#   theme(
#     strip.background = element_blank(),
#     legend.position = 'none',
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(size = 5),
#     strip.text = element_text(size = 6),
#     legend.title  = element_text(size = 6),
#     legend.text  = element_text(size = 4),
#     panel.grid.major.y = element_blank()
#   )
# 
# plot_asv31_asv4

# Create the plot highlighting bloomers ----
plot_asv31_asv4 <- ggplot(asv_tab_filtered, aes(date, rclr)) + 
  labs(x = 'Date', y = 'rCLR', color = 'Synechococcus\nCC9902\nASV', fill = 'Synechococcus\nCC9902\nASV') + 
  
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8, 0.5))) + 
  
  # Smooth for 'bloomer'
  geom_smooth(data = asv_tab_filtered |> filter(bloomer == 'bloomer'), 
              aes(group = asv_num, fill = asv_num), 
              method = 'loess', color = NA, 
              alpha = 0.3, span = 0.09) + 
  # Smooth for non-bloomer
  geom_smooth(data = asv_tab_filtered |> filter(bloomer != 'bloomer'),
              aes(group = asv_num, fill = asv_num),
              method = 'loess', color = NA,
              alpha = 0.2, span = 0.09) +
  scale_color_manual(values = c("asv31" = "#2a0000",
                                "asv4" = "#9c3354",
                                'asv5934' = '#C2AFB3',
                                'asv34' = '#C2AFB3'),
                     labels = c("asv31",
                                "asv4",
                                'asv34, asv5934')) +
  scale_fill_manual(values = c("asv31" = "#2a0000",
                               "asv4" = "#9c3354",
                               'asv5934' = '#C2AFB3',
                               'asv34' = '#C2AFB3'),
                    labels = c("asv31",
                               "asv4",
                               'asv34, asv5934'))+
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1) + 
  theme_bw() + 
  guides(alpha = 'none', 
         color = guide_legend(ncol = 5), 
         fill = guide_legend(ncol = 5)) + 
  theme(strip.background = element_blank(), 
        legend.position = 'none', 
        panel.grid.minor = element_blank(), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), 
        panel.grid.major.y = element_blank())

plot_asv31_asv4

## next group asv1 and asv7 ----
closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv1') |>
  arrange(p_distance)

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv7') |>
  arrange(p_distance) # they are only 4  # 0.007389163

close_distances <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv1') |>
  slice_min(n = 13, order_by = p_distance) |> # to have both bloomers in the cluster
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

# Join the data with the metadata
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

asv_tab_filtered %$%
  asv_num |>
  unique()

asv_tab_filtered <- asv_tab_filtered |>
  dplyr::mutate(new_labs = case_when(asv_num == 'asv1' ~ 'asv1',
                                     asv_num == 'asv7' ~ 'asv7', 
                                     asv_num %in% c("asv343", 'asv14', 'asv3533', 'asv3073', 'asv5327', 'asv5283', 
                                                    'asv4780', 'asv92', 'asv8058', 'asv90', 'asv109', "asv34") ~ 
                                        "asv343, asv14, asv3533, asv3073, asv5327, asv5283, asv4780, asv92, asv8058, asv90, asv109, asv34"))

asv_tab_filtered$new_labs <- factor(asv_tab_filtered$new_labs, c("asv1", "asv7", 
                                                                 "asv343, asv14, asv3533, asv3073, asv5327, asv5283, asv4780, asv92, asv8058, asv90, asv109, asv34"))

# Since data are not smooth then I can't add bloom events points because they are in another position than the one in the plot
# bloom_events_asv1_asv7 <- asv_tab_all_bloo_z_tax |>
#   dplyr::filter(fraction == '0.2') |>
#   dplyr::filter(abundance_type == 'rclr') |>
#   dplyr::filter((
#     asv_num == 'asv1' &
#     date %in% (
#       asv_tab_all_bloo_z_tax |>
#         dplyr::filter(asv_num %in% c('asv1')) |>
#         dplyr::filter(
#           abundance_value > 0.1 &
#             abundance_type == 'relative_abundance' &
#             z_score_ra > 1.96
#         ) |>
#         dplyr::select(date) |>
#         dplyr::pull(date))) |
#       (asv_num == 'asv7' &
#       date %in% (
#         asv_tab_all_bloo_z_tax |>
#           dplyr::filter(asv_num %in% c( 'asv7')) |>
#           dplyr::filter(
#             abundance_value > 0.1 &
#               abundance_type == 'relative_abundance' &
#               z_score_ra > 1.96
#           ) |>
#           dplyr::select(date) |>
#           dplyr::pull(date))) ) |>
#         dplyr::filter(asv_num %in% c('asv1', 'asv7')) |>
#         dplyr::select(date, asv_num, abundance_value) |>
#   dplyr::mutate(fraction = '0.2')
# 
# bloom_events_asv1_asv7_3 <- asv_tab_all_bloo_z_tax |>
#   dplyr::filter(fraction == '3') |>
#   dplyr::filter(abundance_type == 'rclr') |>
#   dplyr::filter((
#     asv_num == 'asv1' &
#       date %in% (
#         asv_tab_all_bloo_z_tax |>
#           dplyr::filter(asv_num %in% c('asv1')) |>
#           dplyr::filter(
#             abundance_value > 0.1 &
#               abundance_type == 'relative_abundance' &
#               z_score_ra > 1.96
#           ) |>
#           dplyr::select(date) |>
#           dplyr::pull(date))) |
#       (asv_num == 'asv7' &
#          date %in% (
#            asv_tab_all_bloo_z_tax |>
#              dplyr::filter(asv_num %in% c( 'asv7')) |>
#              dplyr::filter(
#                abundance_value > 0.1 &
#                  abundance_type == 'relative_abundance' &
#                  z_score_ra > 1.96
#              ) |>
#              dplyr::select(date) |>
#              dplyr::pull(date))) ) |>
#   dplyr::filter(asv_num %in% c('asv1', 'asv7')) |>
#   dplyr::select(date, asv_num, abundance_value) |>
#   dplyr::mutate(fraction = '3')
# 
# bloom_events_asv1_asv7 <- bloom_events_asv1_asv7 |>
#   bind_rows(bloom_events_asv1_asv7_3)

asv_tab_all_bloo_z_tax |>
  colnames()

plot_asv1_asv7 <-  ggplot(asv_tab_filtered, aes(date, rclr, fill = new_labs)) + 
  geom_line(aes(group = new_labs, color = new_labs, alpha = ifelse(bloomer == 'bloomer', 0.9, 0.3))) +  
  #geom_point(data = bloom_events_asv1_asv7, aes(date, abundance_value), fill = NA, color = 'black', shape = 5)+
  geom_smooth(data = asv_tab_filtered |> 
                dplyr::filter(bloomer != 'bloomer'), 
              aes(group = asv_num, fill = new_labs), 
              method = 'loess', alpha = 0.1, span = 0.09, color = NA) +   
  geom_smooth(data = asv_tab_filtered |> 
                dplyr::filter(bloomer == 'bloomer'), 
              aes(group = new_labs, fill = new_labs), 
              method = 'loess', alpha = 0.5, span = 0.09, color = NA) +    
  scale_color_manual(values = c(
    "asv1"  = "#E6CF40",
    "asv7" = "#605991", 
    "asv343, asv14, asv3533, asv3073, asv5327, asv5283, asv4780, asv92, asv8058, asv90, asv109, asv34" = '#C2AFB3')) +    
  scale_fill_manual(values = c(
    "asv1"  = "#E6CF40", 
    "asv7" = "#605991", 
    "asv1"  = "#E6CF40",
    "asv7" = "#605991", 
    "asv343, asv14, asv3533, asv3073, asv5327, asv5283, asv4780, asv92, asv8058, asv90, asv109, asv34" = '#C2AFB3')) + 
  
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1) +   
  labs(x = 'Date', y = 'rCLR', color = 'Synechococcus\nCC9902\nASV', fill = 'Synechococcus\nCC9902\nASV') +   
  theme_bw() +   
  guides(alpha = 'none', 
         color = guide_legend(ncol = 5), 
         fill = guide_legend(ncol = 5)) + 
  theme(strip.background = element_blank(), 
        legend.position = 'bottom', 
        panel.grid.minor = element_blank(), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), 
        panel.grid.major.y = element_blank())

plot_asv1_asv7

legend_plot_asv1_asv7 <- get_legend(plot_asv1_asv7 )

# Create the plot
plot_asv1_asv7 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  geom_point(data = bloom_events_asv1_asv7, aes(date, abundance_value), fill = NA, color = 'black', shape = 5)+
  labs(x = 'Date', y = 'rCLR', color = 'ASV', fill = 'ASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.9, 0.3))) + 
  # Smooth for non-bloomer
  geom_smooth(data = asv_tab_filtered |> filter(bloomer != 'bloomer'),
              aes(group = asv_num, fill = asv_num),
              method = 'loess', color = NA,
              alpha = 0.1, span = 0.09) +
  # Smooth for 'bloomer'
  geom_smooth(data = asv_tab_filtered |> filter(bloomer == 'bloomer'), 
              aes(group = asv_num, fill = asv_num), 
              method = 'loess', color = NA, 
              alpha = 0.5, span = 0.09) + 
  scale_color_manual(values = c(
    "asv1"  = "#E6CF40",
    "asv7" = "#605991", 
    "asv343" = '#C2AFB3', 
    "asv14"  = '#C2AFB3', 
    'asv3533' = '#C2AFB3', 
    'asv3073' = '#C2AFB3', 
    'asv5327'  = '#C2AFB3', 
    'asv5283' = '#C2AFB3', 
    'asv4780' = '#C2AFB3', 
    'asv92' = '#C2AFB3', 
    'asv8058' = '#C2AFB3', 
    'asv90' = '#C2AFB3', 
    'asv109'  = '#C2AFB3', 
    'asv34'  = '#C2AFB3'),
    labels = c("asv1", "asv7", "asv343, asv14, asv3533, asv3073, asv5327, asv5283, asv4780, asv92, asv8058, asv90, asv109, asv34")) +    
  scale_fill_manual(values = c(
    "asv1"  = "#E6CF40", 
    "asv7" = "#605991", 
    "asv343" = '#C2AFB3', 
    "asv14"  = '#C2AFB3', 
    'asv3533' = '#C2AFB3', 
    'asv3073' = '#C2AFB3', 
    'asv5327'  = '#C2AFB3', 
    'asv5283' = '#C2AFB3', 
    'asv4780' = '#C2AFB3', 
    'asv92' = '#C2AFB3', 
    'asv8058' = '#C2AFB3', 
    'asv90' = '#C2AFB3', 
    'asv109'  = '#C2AFB3', 
    'asv34'  = '#C2AFB3'), 
                    labels = c("asv1", "asv7",
                               "asv343, asv14, asv3533, asv3073, asv5327, asv5283, asv4780, asv92, asv8058, asv90, asv109, asv34")
  ) + 
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1) +
  theme_bw() +
  guides(alpha = 'none', 
         color = guide_legend(ncol = 5), 
         fill = guide_legend(ncol = 5)) + 
  theme(strip.background = element_blank(), 
        legend.position = 'none', 
        panel.grid.minor = element_blank(), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), 
        panel.grid.major.y = element_blank())
plot_asv1_asv7

## Cyanobacteria tree ----
palette_cianobac_v2 <- c('Synechococcus CC9902 asv5934' = "#de93a6",
                      'Synechococcus CC9902 asv1337' =  "#6673a9",
                      'Synechococcus CC9902 asv34' =  "#67aa6c",
                      "Synechococcus CC9902 asv31"=   "#2a0000", 
                      "Synechococcus CC9902 asv4"  = "#9c3354", 
                      # second cluster of synechococcus
                      'Synechococcus CC9902 asv1' = '#E6CF40',
                      "Synechococcus CC9902 asv7" = '#605991',
                      "Synechococcus CC9902 asv974" = '#57806B',
                      'Synechococcus CC9902 asv92' = '#E55440',
                      "Synechococcus CC9902 asv109" = '#332A32',
                      "Synechococcus CC9902 asv3533" =  '#805B57',
                      "Synechococcus CC9902 asv343" = '#FFB6DA',
                      "Synechococcus CC9902 asv14"  ='#2D4A2F',
                      "Synechococcus CC9902 asv3073" = '#8A6276',
                      "Synechococcus CC9902 asv5327" = '#2B3147',
                      "Synechococcus CC9902 asv5283" = '#8DA0EB',
                      "Synechococcus CC9902 asv4780"   = '#86AAB0',
                      "Synechococcus CC9902 asv8058" ='#A6A6A6',
                      "Synechococcus CC9902 asv90"    = '#FFFBBC',
                      "Synechococcus CC9902 asv10"  =  '#CE0000' ,  
                      "Synechococcus CC9902 asv698" =    '#C0B6F5')

## I want to plot all Synechococcus closely related to my bloomers
close_distances <- closely_related_bloomers_and_others |>
  distinct(asv_num.y)

tips_to_keep <- close_distances$asv_num.y %in% tree_complete$tip.label

tree <- keep.tip(tree_complete, close_distances$asv_num.y)

# Plot the filtered tree
tree |>
  str()

merged_data <- merge(as_tibble_col(tree$tip.label, column_name = 'asv_num'), tax_bbmo_10y_new, by = "asv_num") |>
  as_tibble() |>
  dplyr::mutate(genus_num = paste0(genus,' ', asv_num),
                asv_num_ed = asv_num)

tree_plot_cianobac <- ggtree(tree, branch.length = T) %<+%   
  merged_data +   
  geom_tiplab(hjust = -0.2, size = 3.5, align = T) +   
  geom_tippoint(aes(color = genus_num), alpha = 0.8, size = 2) +  # Use color for points
  scale_color_manual(values = palette_cianobac_v2) +   # Corrected to scale color instead of fill
  labs(
    title = paste0(unique(merged_data$family), ' ', unique(merged_data$genus)),
    color = 'ASV nº',
    x = 'Phylogenetic Distance'
  ) + 
  guides(color = guide_legend(ncol = 3))+
  theme_tree2() +   
  theme(
    legend.position = "none",  # Show legend
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),  # Remove grid lines
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    #axis.line.x = element_blank(),
    #axis.ticks.x = element_blank(),
    legend.text = element_text(size = 5),
    title = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5)
  )+ 
  xlim(0, 0.25)  # Adjust this value to add space

tree_plot_cianobac

## ASV 17 and ASV77 ------
### beautiful palette for Spingomonadaceae family phylogenetically closely related ----
palette_asv17_asv77 <- c('asv17' = "#6a62b3",
                         "asv955" = "#5c6873",
                         "asv557" = "#1c1a33",
                         "asv77"  = "#f2c549",
                         "asv2929" = "#6bb362",
                         "asv1185" = "#1a2e79",
                         "asv2300" = "#b36a62",
                         "asv876" = "#a8001d",
                         "asv3052" = "#93d8e5",
                         'asv2017' = "#005b00",
                         'asv2115' =   "#93e6bc",
                         "asv6292" = "#e6a093",
                         "asv628"  = "#a69c94",
                         "asv594" =    "#f29f49")

### prepare the data to plot ----
closely_related_bloomers_and_others  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  dplyr::mutate(length_1 = nchar(seq1),
                length_2 = nchar(seq2)) %>%
  dplyr::mutate(p_distance = hammdist/length_1) %>%
  left_join(tax_bbmo_10y_new, by = c('seq1' = 'seq')) %>%
  left_join(tax_bbmo_10y_new, by = c('seq2' = 'seq')) %>%
  filter(p_distance < 0.012) %>%
  distinct(asv_num.x, asv_num.y, p_distance) |>
  dplyr::select(asv_num = asv_num.x, asv_num.y, p_distance) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(asv_num %in% c('asv17', 'asv77')) 

closely_related_bloomers_and_others |>
  arrange(p_distance) |>
  distinct(asv_num.y) |>
  as_vector() ## they are 11 but i will only plot 8 until distance 0.004950495

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv17') |>
  arrange(p_distance)

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv77') |>
  arrange(p_distance) ## they are 13 but i will only plot 8 until distance  0.007425743

x <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv17') |>
  arrange(p_distance)

y <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv77') |>
  arrange(p_distance) ## they are 13 but i will only plot 8 until distance  0.007425743

sphingo_filter_tree <- x |>
  bind_rows(y) |>
  distinct(asv_num.y)

## in the tree I will plot all of them
close_distances <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv17') |>
  slice_min(n = 8, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances_asv77 <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv77') |>
  slice_min(n = 5, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances <- c(close_distances, close_distances_asv77)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

# Join the data with the metadata
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

asv_tab_filtered$asv_num |>
  unique()

palette_asv17_asv77

asv_tab_filtered <- asv_tab_filtered |> 
  dplyr::mutate(new_labs = case_when(
    asv_num == 'asv17' ~ 'asv17',
    asv_num == 'asv77' ~ 'asv77',
    asv_num %in% c('asv2300', 'asv955', 'asv557', 'asv3052', 'asv2929', 'asv1185', 'asv876', 'asv2115', 'asv2017') ~ 
      "asv2300, asv955, asv557, asv3052, asv2929, asv1185, asv876, asv2115, asv2017"))

asv_tab_filtered$new_labs <- factor(asv_tab_filtered$new_labs, levels = c('asv17', 'asv77', 
                                                                          "asv2300, asv955, asv557, asv3052, asv2929, asv1185, asv876, asv2115, asv2017"))
## get the legend ----
plot_asv17_asv77 <- asv_tab_filtered |> 
  ggplot(aes(date, rclr)) + 
  # Smoothed curves for non-bloomer ASVs
  geom_smooth(data = asv_tab_filtered |> dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, fill = new_labs, color = new_labs),  # Map both color and fill
              span = 0.09, alpha = 0.3, linewidth = 0.5) + 
  # Smoothed curves for bloomer ASVs
  geom_smooth(data = asv_tab_filtered |> dplyr::filter(bloomer == 'bloomer'), 
              aes(group = asv_num, fill = new_labs, color = NA),  # Map both color and fill
              span = 0.09) + 
  # Line plot with transparency based on 'bloomer' status
  geom_line(aes(group = asv_num, color = new_labs, alpha = ifelse(bloomer == 'bloomer', 0.8, 0.3))) + 
  # Scale for fill and color, combined under 'new_labs'
  scale_color_manual(values = c(
    'asv17' = "#f29f49", 
    'asv77' = "#005b00", 
    "asv2300, asv955, asv557, asv3052, asv2929, asv1185, asv876, asv2115, asv2017" = '#C2AFB3')) + 
  
  scale_fill_manual(values = c(
    'asv17' = "#f29f49", 
    'asv77' = "#005b00", 
    "asv2300, asv955, asv557, asv3052, asv2929, asv1185, asv876, asv2115, asv2017" = '#C2AFB3')) + 
  
  # Labels
  labs(x = 'Date', y = 'rCLR', fill = 'Sphingomonadaceae\nASV', color = 'Sphingomonadaceae\nASV') + 
  
  # Facet and theme settings
  facet_grid(vars(fraction), labeller = labs_fraction) + 
  theme_bw() + 
  
  # Combined legend for color and fill, arranged in 3 columns
  guides(alpha = 'none', 
         fill = guide_legend(ncol = 3), 
         color = guide_legend(ncol = 3)) + 
  
  # Legend and plot theme
  theme(
    strip.background = element_blank(), 
    legend.position = 'bottom', 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 5), 
    strip.text = element_text(size = 6), 
    legend.title = element_text(size = 6), 
    legend.text = element_text(size = 5), 
    panel.grid.major.y = element_blank()
  )

plot_asv17_asv77

legend_asv17_asv77 <- get_legend(plot_asv17_asv77)

# Create the plot ---
plot_asv17_asv77 <- asv_tab_filtered |> 
  ggplot(aes(date, rclr)) + 
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.9, 0.3))) + 
  # Smooth for non-bloomer
  geom_smooth(data = asv_tab_filtered |> filter(bloomer != 'bloomer'),
              aes(group = asv_num, fill = asv_num),
              method = 'loess', color = NA,
              alpha = 0.1, span = 0.09) +
  # Smooth for 'bloomer'
  geom_smooth(data = asv_tab_filtered |> filter(bloomer == 'bloomer'), 
              aes(group = asv_num, fill = asv_num), 
              method = 'loess', color = NA, 
              alpha = 0.5, span = 0.09) + 
  # Scale for fill, used for both fill and color mapping in legend
  scale_color_manual(values = c(
    'asv17' = "#f29f49", 
    'asv77' = "#005b00", 
    "asv2300, asv955, asv557, asv3052, asv2929, asv1185, asv876, asv2115, asv2017" = '#C2AFB3')) + 
  scale_fill_manual(values = c(
    'asv17' = "#f29f49", 
    'asv77' = "#005b00", 
    "asv2300, asv955, asv557, asv3052, asv2929, asv1185, asv876, asv2115, asv2017" = '#C2AFB3')) + 
  labs(x = 'Date', y = 'rCLR', fill = 'Sphingomonadaceae\nASV', color = 'Sphingomonadaceae\nASV') + 
  facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1) +
  theme_bw() + 
  guides(alpha = 'none', 
         color = guide_legend(ncol = 5), 
         fill = guide_legend(ncol = 5)) + 
  theme(strip.background = element_blank(), 
        legend.position = 'none', 
        panel.grid.minor = element_blank(), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), 
        panel.grid.major.y = element_blank())

plot_asv17_asv77

## plot the phylogenetic tree between those closely related taxa ----
palette_asv17_asv77_v2 <- c('Erythrobacter asv17' = "#6a62b3",
                            "Erythrobacter asv955" = "#5c6873",
                            "Erythrobacter asv557" = "#1c1a33",
                            "Erythrobacter asv77"  = "#f2c549",
                            "Erythrobacter asv2929" = "#6bb362",
                            "Erythrobacter asv1185" = "#1a2e79",
                            "Erythrobacter asv2300" = "#b36a62",
                            "Erythrobacter asv876" = "#a8001d",
                            "Erythrobacter asv3052" = "#93d8e5",
                            'Erythrobacter asv2017' = "#005b00",
                            'Erythrobacter asv2115' =   "#93e6bc",
                            "Erythrobacter asv6292" = "#e6a093",
                            "Erythrobacter asv628"  = "#a69c94",
                            "Erythrobacter asv594" =    "#f29f49")

tips_to_keep <- sphingo_filter_tree$asv_num.y %in% tree_complete$tip.label

tree <- keep.tip(tree_complete, sphingo_filter_tree$asv_num.y)

# Plot the modified tree
tree |>
  str()

merged_data <- merge(as_tibble_col(tree$tip.label, column_name = 'asv_num'), tax_bbmo_10y_new, by = "asv_num") |>
  as_tibble() |>
  dplyr::mutate(genus_num = paste0(genus,' ', asv_num),
                asv_num_ed = asv_num)

tree_plot_asv17_asv77 <- ggtree(tree, branch.length = T) %<+%   
  merged_data +   
  geom_tiplab(hjust = -0.2, size = 3.5, align = TRUE) +   
  geom_tippoint(aes(color = genus_num), alpha = 0.8, size = 4) +  # Use color for points
  scale_color_manual(values = palette_asv17_asv77_v2) +   # Corrected to scale color instead of fill
  labs(
    title = paste0(unique(merged_data$family), ' ', unique(merged_data$genus)),
    color = 'ASV nº',
    x = 'Phylogenetic Distance'
  ) + 
  guides(color = guide_legend(ncol = 4))+
  theme_tree2() +   
  theme(
    legend.position = "none",  # Show legend
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),  # Remove grid lines
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    #axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    title = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )+
  xlim(0, 0.25)
  
tree_plot_asv17_asv77

## COMPOSITION PLOT ------
closely_related_bloomers_plot_v2  <- plot_grid(
  plot_asv1_asv7,
  plot_asv17_asv77,
  plot_asv31_asv4,
  tree_plot_asv17_asv77,
  tree_plot_cianobac,
  ncol = 2,                # One column layout for the main grid
  rel_heights = c(2, 2, 2),
  labels = c('A', 'B', 'C', 'D', 'E'), label_fontface = 'plain',
  align = "hv",           # Align plots horizontally and vertically
  axis = "tblr",          # Display axis labels on all sides
  rel_widths = c(1, 1)    # Set relative widths of columns
)

# Print the final plot
print(closely_related_bloomers_plot_v2)
# 
# ggsave( plot = closely_related_bloomers_plot_v2,
#         filename = 'closely_related_bloomers_plot_v2.pdf',
#         path = 'results/figures/',
#         width = 180, height = 230, units = 'mm')

## Composition plot without the trees -----
closely_related_bloomers_plot_v3  <- plot_grid(
  plot_asv1_asv7,
  legend_plot_asv1_asv7,
  plot_asv31_asv4,
  legend_asv31_asv4,
  plot_asv17_asv77,
  legend_asv17_asv77,
  ncol = 1,                # One column layout for the main grid
  rel_heights = c(2, 0.25, 2, 0.25, 2, 0.25),
  labels = c('A','', 'B', '', 'C'), label_fontface = 'plain',
  label_size = 10,
  align = "hv",           # Align plots horizontally and vertically
  axis = "tblr",          # Display axis labels on all sides
  rel_widths = c(1, 1)    # Set relative widths of columns
)

# Print the final plot
print(closely_related_bloomers_plot_v3)

# ggsave( plot = closely_related_bloomers_plot_v3,
#         filename = 'closely_related_bloomers_plot_v5.pdf',
#         path = 'results/figures/',
#         width = 180, height = 230, units = 'mm')

## Trees go the supplementary material ----
closely_related_bloomers_trees  <- plot_grid(
  tree_plot_cianobac,
  tree_plot_asv17_asv77,
  ncol = 2,                # One column layout for the main grid
  rel_heights = c(1),
  labels = c('A', 'B'), label_fontface = 'plain',label_size = 8, label_x = -0.01,
  align = "hv",           # Align plots horizontally and vertically
  axis = "tblr",          # Display axis labels on all sides
  rel_widths = c(1, 1)    # Set relative widths of columns
)

# Print the final plot
print(closely_related_bloomers_trees)

# ggsave( plot = closely_related_bloomers_trees,
#         filename = 'closely_related_bloomers_trees_v2.pdf',
#         path = 'results/figures/',
#         width = 180, height = 150, units = 'mm')

# Supplementary examples of other closely related bloomers -----
#### Group ASV15, ASV200 and ASV264 #### ----
### beautiful palette for Spingomonadaceae family phylogenetically closely related ----
# palette_asv15_asv2XX <- c(  "asv5"= "#4965f2",    "asv3368" = "#384eba",
#                             "asv86"  = "#3245a6",  
#                             "asv4160" = "#1f2c69",  "asv406" =  "#192354",
#                             "asv9"  = "#9d646c",  
#                             "asv3" = "#649d6f",    "asv1877" = "#5c6073",
#                             "asv1340"  = "#f2c949", "asv1103" = "#7a4e54",
#                             "asv200"  = "#e0909a",  "asv169"  = "#664146",
#                             "asv12" = "#90e0a0",  
#                             "asv1683" = "#416649", "asv2820" = "#83cc91", 
#                             "asv2"  = "#cfaa3e",   
#                             "asv272"  ="#9787f5","asv12265"="#b3626e",
#                             "asv3455"  ="#6e62b3","asv2159" ="#64607f",
#                             "asv529" ="#9db362",  "asv593"="#b3628e",   "asv4053" ="#777f60", 
#                             "asv1770" ="#2a3312", "asv4840" ="#b36d62", "asv15" ="#b3a09e",
#                             "asv264" = "#3f3866",  "asv472" = '#332C26')

palette_asv15_asv2XX <- c(  "asv5"= "#C2AFB3",    "asv3368" = "#C2AFB3", # I will only highlight those blooming taxa
                            "asv86"  = "#C2AFB3",  
                            "asv4160" = "#C2AFB3",  "asv406" =  "#C2AFB3",
                            "asv9"  = "#C2AFB3",  
                            "asv3" = "#C2AFB3",    "asv1877" = "#C2AFB3",
                            "asv1340"  = "#C2AFB3", "asv1103" = "#C2AFB3",
                            "asv200"  = "#e0909a",  "asv169"  = "#C2AFB3",
                            "asv12" = "#C2AFB3",  
                            "asv1683" = "#C2AFB3", "asv2820" = "#C2AFB3", 
                            "asv2"  = "#C2AFB3",   
                            "asv272"  ="#C2AFB3","asv12265"="#C2AFB3",
                            "asv3455"  ="#C2AFB3","asv2159" ="#C2AFB3",
                            "asv529" ="#C2AFB3",  "asv593"="#C2AFB3",   "asv4053" ="#C2AFB3", 
                            "asv1770" ="#C2AFB3", "asv4840" ="#C2AFB3", "asv15" ="#2a3312",
                            "asv264" = "#3f3866",  "asv472" = '#C2AFB3')

### prepare the data to plot ----
closely_related_bloomers_and_others  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  dplyr::mutate(length_1 = nchar(seq1),
                length_2 = nchar(seq2)) %>%
  dplyr::mutate(p_distance = hammdist/length_1) %>%
  left_join(tax_bbmo_10y_new, by = c('seq1' = 'seq')) %>%
  left_join(tax_bbmo_10y_new, by = c('seq2' = 'seq')) %>%
  filter(p_distance < 0.012) %>%
  distinct(asv_num.x, asv_num.y, p_distance) |>
  dplyr::select(asv_num = asv_num.x, asv_num.y, p_distance) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(asv_num %in% c('asv15', 'asv200', 'asv264')) 

closely_related_bloomers_and_others |>
  arrange(p_distance) |>
  distinct(asv_num.y) |>
  as_vector() ## they are 32 but i will only plot 8 until distance 0.004950495

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv15') |>
  arrange(p_distance) ## it's too far from the others it is at 0.009. Therefore I will only plot it in the tree.

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv200') |>
  arrange(p_distance) ## they are 20 but i will only plot 11 until distance  0.007425743

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv264') |>
  arrange(p_distance) ## they are 26 but i will only plot x until distance  0.007425743

x <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv264') |>
  arrange(p_distance)

y <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv200') |>
  arrange(p_distance) ## they are 13 but i will only plot 8 until distance  0.007425743

sar11_filter_tree <- x |>
  bind_rows(y) |>
  distinct(asv_num.y)

## in the tree I will plot all of them
close_distances <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv200') |>
  slice_min(n = 11, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances_asv264 <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv264') |>
  slice_min(n = 17, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances <- c(close_distances, close_distances_asv264)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

# Join the data with the metadata
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

#create the legend
plot_asv15_asv2XX <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Date', y = 'rCLR', color = 'Sphingomonadaceae\nASV', fill = 'Sphingomonadaceae\nASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.5) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
  scale_color_manual(values =  palette_asv15_asv2XX) +
  scale_fill_manual(values = palette_asv15_asv2XX) +
  facet_grid(vars(fraction), labeller = labs_fraction)+
  theme_bw() +
  guides(alpha = 'none',
         color = guide_legend(ncol = 6),
         fill = guide_legend(ncol = 6))+
  theme(
    strip.background = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 5),
    strip.text = element_text(size = 6),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 5),
    panel.grid.major.y = element_blank()
  )

legend_asv15_asv2XX <- get_legend(plot_asv15_asv2XX)

# Create the plot
plot_asv15_asv2XX <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Date', y = 'rCLR', color = 'ASV', fill = 'ASV') +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = NA, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0) +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0. , 0.003))) +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.9 , 0.0001))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = NA, fill = asv_num),linewidth = 0, 
              span = 0.09, alpha = 0.5) +
  scale_color_manual(values =  palette_asv15_asv2XX) +
  scale_fill_manual(values = palette_asv15_asv2XX) +
  facet_grid(vars(fraction), labeller = labs_fraction)+
  theme_bw() +
  guides(alpha = 'none')+
  theme(
    strip.background = element_blank(),
    legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size =5),
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 4),
    panel.grid.major.y = element_blank()
  )

plot_asv15_asv2XX 

## plot the phylogenetic tree between those closely related taxa ----
# palette_asv15_asv2XX_v2 <- c(  "Clade Ia asv5"= "#4965f2",    "asv3368" = "#384eba",
#                                "Clade Ia asv86"  = "#3245a6",  
#                                "Clade Ia asv4160" = "#1f2c69",  "Clade Ia asv406" =  "#192354",
#                                "Clade Ia asv9"  = "#9d646c",  
#                                "Clade Ia asv3" = "#649d6f",    "Clade Ia asv1877" = "#5c6073",
#                                "Clade Ia asv1340"  = "#f2c949", "Clade Ia asv1103" = "#7a4e54",
#                                "Clade Ia asv200"  = "#e0909a",  "Clade Ia asv169"  = "#664146",
#                                "Clade Ia asv12" = "#90e0a0",  
#                                "Clade Ia asv1683" = "#416649", "Clade Ia asv2820" = "#83cc91", 
#                                "Clade Ia asv2"  = "#cfaa3e",   
#                                "Clade Ia asv272"  ="#9787f5","Clade Ia asv12265"="#b3626e",
#                                "NA asv3455"  ="#6e62b3","Clade Ia asv2159" ="#64607f",
#                                "Clade Ia asv529" ="#9db362",  "Clade Ia asv593"="#b3628e",   "Clade Ia asv4053" ="#777f60", 
#                                "Clade Ia asv1770" ="#2a3312", "Clade Ia asv4840" ="#b36d62", "Clade Ia asv15" ="#b3a09e",
#                                "Clade Ia asv264" = "#3f3866",  "Clade Ia asv472" = '#332C26')

palette_asv15_asv2XX_v2 <- c(  "Clade Ia asv5"= "#C2AFB3",    "Clade Ia asv3368" = "#C2AFB3", # I will only highlight those blooming taxa
    "Clade Ia asv86"  = "#C2AFB3",  
    "Clade Ia asv4160" = "#C2AFB3",  "Clade Ia asv406" =  "#C2AFB3",
    "Clade Ia asv9"  = "#C2AFB3",  
    "Clade Ia asv3" = "#C2AFB3",    "Clade Ia asv1877" = "#C2AFB3",
    "Clade Ia asv1340"  = "#C2AFB3", "Clade Ia asv1103" = "#C2AFB3",
    "Clade Ia asv200"  = "#e0909a",  "Clade Ia asv169"  = "#C2AFB3",
    "Clade Ia asv12" = "#C2AFB3",  
    "Clade Ia asv1683" = "#C2AFB3", "Clade Ia asv2820" = "#C2AFB3", 
    "Clade Ia asv2"  = "#C2AFB3",   
    "Clade Ia asv272"  ="#C2AFB3","Clade Ia asv12265"="#C2AFB3",
    "Clade Ia asv3455"  ="#C2AFB3","Clade Ia asv2159" ="#C2AFB3",
    "Clade Ia asv529" ="#C2AFB3",  "Clade Ia asv593"="#C2AFB3",   "Clade Ia asv4053" ="#C2AFB3", 
    "Clade Ia asv1770" ="#C2AFB3", "Clade Ia asv4840" ="#C2AFB3", "Clade Ia asv15" ="#2a3312",
    "Clade Ia asv264" = "#3f3866",  "Clade Ia asv472" = '#C2AFB3')

tips_to_keep <- sar11_filter_tree$asv_num.y %in% tree_complete$tip.label

tree <- keep.tip(tree_complete, sar11_filter_tree$asv_num.y)

# Plot the modified tree
tree |>
  str()

merged_data <- merge(as_tibble_col(tree$tip.label, column_name = 'asv_num'), tax_bbmo_10y_new, by = "asv_num") |>
  as_tibble() |>
  dplyr::mutate(genus_num = paste0(genus,' ', asv_num),
                asv_num_ed = asv_num)

tree_plot_asv15_asv2XX <- ggtree(tree, branch.length = T) %<+%   
  merged_data +   
  geom_tiplab(hjust = -0.25, size = 2.5, align = TRUE) +   
  geom_tippoint(aes(color = genus_num), alpha = 0.8, size = 3) +  # Use color for points
  scale_color_manual(values = palette_asv15_asv2XX_v2) +   # Corrected to scale color instead of fill
  labs(
    title = paste0(unique(merged_data$order),' ',
                   unique(merged_data$family), ' ', unique(merged_data$genus)),
    color = 'ASV nº'
  ) + 
  guides(color = guide_legend(ncol = 4))+
  theme_tree2() +   
  theme(
    legend.position = "none",  # Show legend
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),  # Remove grid lines
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    #axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    title = element_text(size = 6),
    plot.margin = margin(2.5, 5, 5, 2.5)
  )+
  xlim(0, 0.25)

tree_plot_asv15_asv2XX

## create a composition with the time series and ----

closely_related_bloomers_sar11cladea_plot  <- plot_grid(
  plot_asv15_asv2XX,
  tree_plot_asv15_asv2XX,
  ncol = 2,                # One column layout for the main grid
  rel_heights = c(1),
  labels = c('A', 'B'), label_fontface = 'plain',
  label_size = 10,
  align = "hv",           # Align plots horizontally and vertically
  axis = "tblr",          # Display axis labels on all sides
  rel_widths = c(2, 1)    # Set relative widths of columns
)

# ggsave( plot = closely_related_bloomers_sar11cladea_plot,
#         filename = 'closely_related_bloomers_sar11cladea_plot_v2.pdf',
#         path = 'results/figures/',
#         width = 180, height = 150, units = 'mm')

#### Group ASV42 and ASV49 #### -----
### beautiful palette for Spingomonadaceae family phylogenetically closely related ----
# palette_asv42_asv49 <- c(  "asv1271"= "#4965f2",    "asv87"   = "#384eba",
#                            "asv49"   = "#3245a6",  
#                            "asv1819" = "#1f2c69",   "asv1942"   =  "#192354",
#                            "asv5485"  = "#9d646c",  
#                            "asv477" = "#649d6f",    "asv679"  = "#5c6073",
#                            "asv1051" = "#f2c949",  "asv1021"   = "#7a4e54",
#                            "asv329"  = "#e0909a", "asv531"   = "#664146",
#                            "asv10095" = "#90e0a0",  
#                            "asv1294"  = "#416649", "asv808"    = "#83cc91", 
#                             "asv258"  = '#A7A4C4',   
#                            "asv1004"  ="#9787f5", "asv520"="#b3626e",
#                            "asv414"   ="#6e62b3","asv4588"  ="#64607f",
#                            "asv132" ="#9db362",  "asv124" ="#b3628e",   "asv2516" ="#777f60", 
#                            "asv670" ="#2a3312", "asv6017"  ="#b36d62", "asv1061" ="#b3a09e",
#                            "asv281" = "#3f3866",   "asv1822" = '#332C26',
#                            "asv500" = '#B2211B',   
#                              "asv42" = "#cfaa3e",
#                              "asv2077" =  '#7A1612')

palette_asv42_asv49 <- c(  "asv1271"= "#C2AFB3",    "asv87"   = "#C2AFB3",
                           "asv49"   = "#3245a6",  "asv42" = "#cfaa3e",
                           "asv1819" = "#C2AFB3",   "asv1942"   =  "#C2AFB3",
                           "asv5485"  = "#C2AFB3",  
                           "asv477" = "#C2AFB3",    "asv679"  = "#C2AFB3",
                           "asv1051" = "#C2AFB3",  "asv1021"   = "#C2AFB3",
                           "asv329"  = "#C2AFB3", "asv531"   = "#C2AFB3",
                           "asv10095" = "#C2AFB3",  
                           "asv1294"  = "#C2AFB3", "asv808"    = "#C2AFB3", 
                           "asv258"  = '#C2AFB3',   
                           "asv1004"  ="#C2AFB3", "asv520"= "#C2AFB3",
                           "asv414"   ="#C2AFB3","asv4588"  ="#C2AFB3",
                           "asv132" = "#C2AFB3",  "asv124" ="#C2AFB3",   "asv2516" ="#C2AFB3", 
                           "asv670" = "#C2AFB3", "asv6017"  ="#C2AFB3", "asv1061" ="#b3a09e",
                           "asv281" = "#C2AFB3",   "asv1822" = '#C2AFB3',
                           "asv500" = '#C2AFB3',   "asv2077" =  '#C2AFB3')

### prepare the data to plot ----
closely_related_bloomers_and_others  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  dplyr::mutate(length_1 = nchar(seq1),
                length_2 = nchar(seq2)) %>%
  dplyr::mutate(p_distance = hammdist/length_1) %>%
  left_join(tax_bbmo_10y_new, by = c('seq1' = 'seq')) %>%
  left_join(tax_bbmo_10y_new, by = c('seq2' = 'seq')) %>%
  filter(p_distance < 0.012) %>%
  distinct(asv_num.x, asv_num.y, p_distance) |>
  dplyr::select(asv_num = asv_num.x, asv_num.y, p_distance) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(asv_num %in% c('asv42', 'asv49')) 

closely_related_bloomers_and_others |>
  arrange(p_distance) |>
  distinct(asv_num.y) |>
  as_vector() ## they are 29 but i will only plot until distance 0.002331002

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv42') |>
  arrange(p_distance) ## it's too far from the others it is at 0.009. Therefore I will only plot it in the tree.

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv49') |>
  arrange(p_distance) ## they are 20 but i will only plot 11 until distance  0.002331002

x <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv42') |>
  arrange(p_distance)

y <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv49') |>
  arrange(p_distance) ## they are 13 but i will only plot 8 until distance  0.002331002

sar86_filter_tree <- x |>
  bind_rows(y) |>
  distinct(asv_num.y)

## in the tree I will plot all of them
close_distances <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv42') |>
  slice_min(n = 3, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances_asv49 <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv49') |>
  slice_min(n = 3, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances <- c(close_distances, close_distances_asv49)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

# Join the data with the metadata
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

#create the legend
plot_asv42_asv49 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Date', y = 'rCLR', color = 'Sphingomonadaceae\nASV', fill = 'Sphingomonadaceae\nASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.5) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
  scale_color_manual(values =  palette_asv42_asv49) +
  scale_fill_manual(values = palette_asv42_asv49) +
  facet_grid(vars(fraction), labeller = labs_fraction)+
  theme_bw() +
  guides(alpha = 'none',
         color = guide_legend(ncol = 6),
         fill = guide_legend(ncol = 6))+
  theme(
    strip.background = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 5),
    strip.text = element_text(size = 6),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 5),
    panel.grid.major.y = element_blank()
  )

legend_asv42_asv49 <- get_legend(plot_asv42_asv49)

# Create the plot
plot_asv42_asv49 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Date', y = 'rCLR', color = 'ASV', fill = 'ASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = NA, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.0) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = NA, fill = asv_num),linewidth = 0.0, span = 0.09) +
  scale_color_manual(values =  palette_asv42_asv49) +
  scale_fill_manual(values = palette_asv42_asv49) +
  facet_grid(vars(fraction), labeller = labs_fraction)+
  theme_bw() +
  guides(alpha = 'none')+
  theme(
    strip.background = element_blank(),
    legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size =5),
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 4),
    panel.grid.major.y = element_blank()
  )

plot_asv42_asv49 

## plot the phylogenetic tree between those closely related taxa ----
# palette_asv42_asv49_v2 <- c(  "SAR86 clade asv1271"= "#4965f2",    "SAR86 clade asv87"   = "#384eba",
#                               "SAR86 clade asv49"   = "#3245a6",  
#                               "SAR86 clade asv1819" = "#1f2c69",   "SAR86 clade asv1942"   =  "#192354",
#                               "SAR86 clade asv5485"  = "#9d646c",  
#                               "SAR86 clade asv477" = "#649d6f",    "SAR86 clade asv679"  = "#5c6073",
#                               "SAR86 clade asv1051" = "#f2c949",  "SAR86 clade asv1021"   = "#7a4e54",
#                               "SAR86 clade asv329"  = "#e0909a", "SAR86 clade asv531"   = "#664146",
#                               "SAR86 clade asv10095" = "#90e0a0",  
#                               "SAR86 clade asv1294"  = "#416649", "SAR86 clade asv808"    = "#83cc91", 
#                               "SAR86 clade asv258"  = '#A7A4C4',   
#                               "SAR86 clade asv1004"  ="#9787f5", "SAR86 clade asv520"="#b3626e",
#                               "SAR86 clade asv414"   ="#6e62b3","SAR86 clade asv4588"  ="#64607f",
#                               "SAR86 clade asv132" ="#9db362",  "SAR86 clade asv124" ="#b3628e",   "SAR86 clade asv2516" ="#777f60", 
#                               "SAR86 clade asv670" ="#2a3312", "SAR86 clade asv6017"  ="#b36d62", "SAR86 clade asv1061" ="#b3a09e",
#                               "SAR86 clade asv281" = "#3f3866",   "SAR86 clade asv1822" = '#332C26',
#                               "SAR86 clade asv500" = '#B2211B',   
#                               "SAR86 clade asv42" = "#cfaa3e",
#                               "SAR86 clade asv2077" =  '#7A1612')

palette_asv42_asv49_v2 <- c(  "SAR86 clade asv1271"= "#C2AFB3",    "SAR86 clade asv87"   = "#C2AFB3",
                           "SAR86 clade asv49"   = "#3245a6",  "SAR86 clade asv42" = "#cfaa3e",
                           "SAR86 clade asv1819" = "#C2AFB3",   "SAR86 clade asv1942"   =  "#C2AFB3",
                           "SAR86 clade asv5485"  = "#C2AFB3",  
                           "SAR86 clade asv477" = "#C2AFB3",    "SAR86 clade asv679"  = "#C2AFB3",
                           "SAR86 clade asv1051" = "#C2AFB3",  "SAR86 clade asv1021"   = "#C2AFB3",
                           "SAR86 clade asv329"  = "#C2AFB3", "SAR86 clade asv531"   = "#C2AFB3",
                           "SAR86 clade asv10095" = "#C2AFB3",  
                           "SAR86 clade asv1294"  = "#C2AFB3", "SAR86 clade asv808"    = "#C2AFB3", 
                           "SAR86 clade asv258"  = '#C2AFB3',   
                           "SAR86 clade asv1004"  ="#C2AFB3", "SAR86 clade asv520"= "#C2AFB3",
                           "SAR86 clade asv414"   ="#C2AFB3","SAR86 clade asv4588"  ="#C2AFB3",
                           "SAR86 clade asv132" = "#C2AFB3",  "SAR86 clade asv124" ="#C2AFB3",   "SAR86 clade asv2516" ="#C2AFB3", 
                           "SAR86 clade asv670" = "#C2AFB3", "SAR86 clade asv6017"  ="#C2AFB3", "SAR86 clade asv1061" ="#b3a09e",
                           "SAR86 clade asv281" = "#C2AFB3",   "SAR86 clade asv1822" = '#C2AFB3',
                           "SAR86 clade asv500" = '#C2AFB3',   "SAR86 clade asv2077" =  '#C2AFB3')

tips_to_keep <- sar86_filter_tree$asv_num.y %in% tree_complete$tip.label

tree <- keep.tip(tree_complete, sar86_filter_tree$asv_num.y)

merged_data$genus_num

# Plot the modified tree
tree |>
  str()

merged_data <- merge(as_tibble_col(tree$tip.label, column_name = 'asv_num'), tax_bbmo_10y_new, by = "asv_num") |>
  as_tibble() |>
  dplyr::mutate(genus_num = paste0(family,' ', asv_num),
                asv_num_ed = asv_num)

tree_plot_asv42_asv49 <- ggtree(tree, branch.length = T) %<+%   
  merged_data +   
  geom_tiplab(hjust = -0.25, size = 2.5, align = TRUE) +   
  geom_tippoint(aes(color = genus_num), alpha = 0.8, size = 3) +  # Use color for points
  scale_color_manual(values = palette_asv42_asv49_v2) +   # Corrected to scale color instead of fill
  labs(
    title = paste0(unique(merged_data$order),' ',
                   unique(merged_data$family)),
    color = 'ASV nº'
  ) + 
  guides(color = guide_legend(ncol = 4))+
  theme_tree2() +   
  theme(
    legend.position = "none",  # Show legend
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),  # Remove grid lines
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    #axis.line.x = element_blank(),
    #axis.ticks.x = element_blank(),
    title = element_text(size = 6),
    plot.margin = margin(2.5, 5, 5, 2.5)
  )+
  xlim(0, 0.25)

tree_plot_asv42_asv49

## create a composition with the time series and ----
closely_related_bloomers_sar86clade_plot  <- plot_grid(
  plot_asv42_asv49,
  tree_plot_asv42_asv49,
  ncol = 2,                # One column layout for the main grid
  rel_heights = c(1),
  labels = c('C', 'D'), label_fontface = 'plain',
  label_size = 10,
  align = "hv",           # Align plots horizontally and vertically
  axis = "tblr",          # Display axis labels on all sides
  rel_widths = c(2, 1)    # Set relative widths of columns
)

# ggsave( plot = closely_related_bloomers_sar86clade_plot,
#         filename = 'closely_related_bloomers_sar86clade_plot_V2.pdf',
#         path = 'results/figures/',
#         width = 180, height = 150, units = 'mm')


## I will add another example ASV11 -----
### beautiful palette for Glaciecolas ----
palette_asv11 <- c('asv11' = "#f2c949",  'asv1368' = '#C2AFB3',  'asv2977' = '#C2AFB3',  'asv18888' = '#C2AFB3',
                   'asv764'= '#C2AFB3',  'asv10327' = '#C2AFB3', 'asv7413' = '#C2AFB3')

### prepare the data to plot ----
closely_related_bloomers_and_others  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  dplyr::mutate(length_1 = nchar(seq1),
                length_2 = nchar(seq2)) %>%
  dplyr::mutate(p_distance = hammdist/length_1) %>%
  left_join(tax_bbmo_10y_new, by = c('seq1' = 'seq')) %>%
  left_join(tax_bbmo_10y_new, by = c('seq2' = 'seq')) %>%
  filter(p_distance < 0.012) %>%
  distinct(asv_num.x, asv_num.y, p_distance) |>
  dplyr::select(asv_num = asv_num.x, asv_num.y, p_distance) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(asv_num %in% c('asv11'))

# asv11 asv1368 asv2977 asv18888 asv764 asv10327 asv7413
closely_related_bloomers_and_others |>
  arrange(p_distance) |>
  distinct(asv_num.y) |>
  as_vector() ## they are 32 but i will only plot 8 until distance 0.004950495

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv11') |>
  arrange(p_distance) ## it's too far from the others it is at 0.009. Therefore I will only plot it in the tree.

x <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv11') |>
  arrange(p_distance) |>
  dplyr::select(asv_num = asv_num.y, asv_num.y = asv_num, p_distance, seq, domain, phylum, class, order, family , genus)

y <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv11') |>
  arrange(p_distance) ## they are 13 but i will only plot 8 until distance  0.007425743

asv11_filter_tree <- x |>
  bind_rows(y) |>
  distinct(asv_num.y)

## in the tree I will plot all of them
close_distances <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv11') |>
  slice_min(n = 11, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num'))

close_distances |>
  distinct(family, genus)

close_distances <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv11') |>
  slice_min(n = 11, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances <- c(close_distances)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

# Join the data with the metadata
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num_f ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num_f ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

#create the legend
plot_asv11 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Date', y = 'rCLR', color = 'Sphingomonadaceae\nASV', fill = 'Sphingomonadaceae\nASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.5) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
  scale_color_manual(values =  palette_asv15_asv2XX) +
  scale_fill_manual(values = palette_asv15_asv2XX) +
  facet_grid(vars(fraction), labeller = labs_fraction)+
  theme_bw() +
  guides(alpha = 'none',
         color = guide_legend(ncol = 6),
         fill = guide_legend(ncol = 6))+
  theme(
    strip.background = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 5),
    strip.text = element_text(size = 6),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 5),
    panel.grid.major.y = element_blank()
  )

legend_asv11 <- get_legend(plot_asv11)

# Create the plot
plot_asv11 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Date', y = 'rCLR', color = 'ASV', fill = 'ASV') +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = NA, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0) +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0. , 0.003))) +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.9 , 0.0001))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = NA, fill = asv_num),linewidth = 0, 
              span = 0.09, alpha = 0.5) +
  scale_color_manual(values =  palette_asv11) +
  scale_fill_manual(values = palette_asv11) +
  facet_grid(vars(fraction), labeller = labs_fraction)+
  theme_bw() +
  guides(alpha = 'none')+
  theme(
    strip.background = element_blank(),
    legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size =5),
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 4),
    panel.grid.major.y = element_blank()
  )

plot_asv11

## plot the phylogenetic tree between those closely related taxa ----
asv_tab_filtered |>
  unique(family, genus)

asv_tab_filtered |>
  colnames()

palette_asv11_V2 <- c('Glaciecola asv11' = "#f2c949",  "Glaciecola asv1368" = '#C2AFB3',  
                      "Glaciecola asv2977" = '#C2AFB3',  "Glaciecola asv1888" = '#C2AFB3',
                                                "Glaciecola asv764" = '#C2AFB3',  "Glaciecola asv10327" = '#C2AFB3', 
                      "Glaciecola asv7413" = '#C2AFB3')

tips_to_keep <- asv11_filter_tree$asv_num.y %in% tree_complete$tip.label

tree <- keep.tip(tree_complete, asv11_filter_tree$asv_num.y)

# Plot the modified tree
tree |>
  str()

merged_data <- merge(as_tibble_col(tree$tip.label, column_name = 'asv_num'), tax_bbmo_10y_new, by = "asv_num") |>
  as_tibble() |>
  dplyr::mutate(genus_num = paste0(genus,' ', asv_num),
                asv_num_ed = asv_num)

tree_plot_asv11 <- ggtree(tree, branch.length = T) %<+%   
  merged_data +   
  geom_tiplab(hjust = -0.25, size = 2.5, align = TRUE) +   
  geom_tippoint(aes(color = genus_num), alpha = 0.8, size = 3) +  # Use color for points
  scale_color_manual(values = palette_asv11_V2) +   # Corrected to scale color instead of fill
  labs(
    title = paste0(unique(merged_data$order),' ',
                   unique(merged_data$family), ' ', unique(merged_data$genus)),
    color = 'ASV nº'
  ) + 
  guides(color = guide_legend(ncol = 4))+
  theme_tree2() +   
  theme(
    legend.position = "none",  # Show legend
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),  # Remove grid lines
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    #axis.line.x = element_blank(),
    #axis.ticks.x = element_blank(),
    title = element_text(size = 6),
    plot.margin = margin(2.5, 5, 5, 2.5)
  )+
  xlim(0, 0.25)

tree_plot_asv11

## create a composition with the time series and ----

closely_related_bloomers_glaciecola11_plot  <- plot_grid(
  plot_asv11,
  tree_plot_asv11,
  ncol = 2,                # One column layout for the main grid
  rel_heights = c(1),
  labels = c('E', 'F'), label_fontface = 'plain',
  label_size = 10,
  align = "hv",           # Align plots horizontally and vertically
  axis = "tblr",          # Display axis labels on all sides
  rel_widths = c(2, 1)    # Set relative widths of columns
)

# ggsave( plot = closely_related_bloomers_glaciecola11_plot,
#         filename = 'closely_related_bloomers_glaciecola11_plot_v2.pdf',
#         path = 'results/figures/',
#         width = 180, height = 150, units = 'mm')

closely_related_bloomers_supplementary  <- plot_grid(
  closely_related_bloomers_sar11cladea_plot,
  closely_related_bloomers_sar86clade_plot ,
  closely_related_bloomers_glaciecola11_plot,
  ncol = 1,                # One column layout for the main grid
  rel_heights = c(1),
  #labels = c('C', 'D'), label_fontface = 'plain',
  label_size = 10,
  align = "hv",           # Align plots horizontally and vertically
  axis = "tblr",          # Display axis labels on all sides
  rel_widths = c(2, 1)    # Set relative widths of columns
)

# ggsave( plot = closely_related_bloomers_supplementary,
#         filename = 'closely_related_bloomers_supplementary_v1.pdf',
#         path = 'results/figures/',
#         width = 180, height = 250, units = 'mm')

#### Group ASV42 and ASV49 #### -----
### beautiful palette for Spingomonadaceae family phylogenetically closely related ----
# palette_asv42_asv49 <- c(  "asv1271"= "#4965f2",    "asv87"   = "#384eba",
#                            "asv49"   = "#3245a6",  
#                            "asv1819" = "#1f2c69",   "asv1942"   =  "#192354",
#                            "asv5485"  = "#9d646c",  
#                            "asv477" = "#649d6f",    "asv679"  = "#5c6073",
#                            "asv1051" = "#f2c949",  "asv1021"   = "#7a4e54",
#                            "asv329"  = "#e0909a", "asv531"   = "#664146",
#                            "asv10095" = "#90e0a0",  
#                            "asv1294"  = "#416649", "asv808"    = "#83cc91", 
#                             "asv258"  = '#A7A4C4',   
#                            "asv1004"  ="#9787f5", "asv520"="#b3626e",
#                            "asv414"   ="#6e62b3","asv4588"  ="#64607f",
#                            "asv132" ="#9db362",  "asv124" ="#b3628e",   "asv2516" ="#777f60", 
#                            "asv670" ="#2a3312", "asv6017"  ="#b36d62", "asv1061" ="#b3a09e",
#                            "asv281" = "#3f3866",   "asv1822" = '#332C26',
#                            "asv500" = '#B2211B',   
#                              "asv42" = "#cfaa3e",
#                              "asv2077" =  '#7A1612')

palette_asv42_asv49 <- c(  "asv1271"= "#C2AFB3",    "asv87"   = "#C2AFB3",
                           "asv49"   = "#3245a6",  "asv42" = "#cfaa3e",
                           "asv1819" = "#C2AFB3",   "asv1942"   =  "#C2AFB3",
                           "asv5485"  = "#C2AFB3",  
                           "asv477" = "#C2AFB3",    "asv679"  = "#C2AFB3",
                           "asv1051" = "#C2AFB3",  "asv1021"   = "#C2AFB3",
                           "asv329"  = "#C2AFB3", "asv531"   = "#C2AFB3",
                           "asv10095" = "#C2AFB3",  
                           "asv1294"  = "#C2AFB3", "asv808"    = "#C2AFB3", 
                           "asv258"  = '#C2AFB3',   
                           "asv1004"  ="#C2AFB3", "asv520"= "#C2AFB3",
                           "asv414"   ="#C2AFB3","asv4588"  ="#C2AFB3",
                           "asv132" = "#C2AFB3",  "asv124" ="#C2AFB3",   "asv2516" ="#C2AFB3", 
                           "asv670" = "#C2AFB3", "asv6017"  ="#C2AFB3", "asv1061" ="#b3a09e",
                           "asv281" = "#C2AFB3",   "asv1822" = '#C2AFB3',
                           "asv500" = '#C2AFB3',   "asv2077" =  '#C2AFB3')

### prepare the data to plot ----
closely_related_bloomers_and_others  <- hammdist.all %>%
  purrr::reduce(bind_rows) %>%
  dplyr::mutate(length_1 = nchar(seq1),
                length_2 = nchar(seq2)) %>%
  dplyr::mutate(p_distance = hammdist/length_1) %>%
  left_join(tax_bbmo_10y_new, by = c('seq1' = 'seq')) %>%
  left_join(tax_bbmo_10y_new, by = c('seq2' = 'seq')) %>%
  filter(p_distance < 0.012) %>%
  distinct(asv_num.x, asv_num.y, p_distance) |>
  dplyr::select(asv_num = asv_num.x, asv_num.y, p_distance) |>
  left_join(tax_bbmo_10y_new) |>
  dplyr::filter(asv_num %in% c('asv42', 'asv49')) 

closely_related_bloomers_and_others |>
  arrange(p_distance) |>
  distinct(asv_num.y) |>
  as_vector() ## they are 29 but i will only plot until distance 0.002331002

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv42') |>
  arrange(p_distance) ## it's too far from the others it is at 0.009. Therefore I will only plot it in the tree.

closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv49') |>
  arrange(p_distance) ## they are 20 but i will only plot 11 until distance  0.002331002

x <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv42') |>
  arrange(p_distance)

y <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv49') |>
  arrange(p_distance) ## they are 13 but i will only plot 8 until distance  0.002331002

sar86_filter_tree <- x |>
  bind_rows(y) |>
  distinct(asv_num.y)

## in the tree I will plot all of them
close_distances <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv42') |>
  slice_min(n = 3, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances_asv49 <- closely_related_bloomers_and_others |>
  dplyr::filter(asv_num == 'asv49') |>
  slice_min(n = 3, order_by = p_distance) |>
  pivot_longer(cols = starts_with('asv_num')) %$%
  unique(value)

close_distances <- c(close_distances, close_distances_asv49)

# Filter the data for the current ASV number
asv_tab_filtered_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_distances)

asv_tab_filtered_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_distances)

# Join the data with the metadata
asv_tab_joined_02 <- asv_tab_filtered_02 |>
  left_join(m_02) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num ~ 'no-bloomer'))

asv_tab_joined_3 <- asv_tab_filtered_3 |>
  left_join(m_3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(bloomer = case_when(asv_num %in% bloo_taxonomy$asv_num ~ 'bloomer',
                                    !asv_num %in% bloo_taxonomy$asv_num ~ 'no-bloomer'))

asv_tab_filtered <- asv_tab_joined_3 |>
  bind_rows(asv_tab_joined_02)

#create the legend
plot_asv42_asv49 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Date', y = 'rCLR', color = 'Sphingomonadaceae\nASV', fill = 'Sphingomonadaceae\nASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.5) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = asv_num, fill = asv_num), span = 0.09) +
  scale_color_manual(values =  palette_asv42_asv49) +
  scale_fill_manual(values = palette_asv42_asv49) +
  facet_grid(vars(fraction), labeller = labs_fraction)+
  theme_bw() +
  guides(alpha = 'none',
         color = guide_legend(ncol = 6),
         fill = guide_legend(ncol = 6))+
  theme(
    strip.background = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 5),
    strip.text = element_text(size = 6),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 5),
    panel.grid.major.y = element_blank()
  )

legend_asv42_asv49 <- get_legend(plot_asv42_asv49)

# Create the plot
plot_asv42_asv49 <- ggplot(asv_tab_filtered, aes(date, rclr)) +
  labs(x = 'Date', y = 'rCLR', color = 'ASV', fill = 'ASV') +
  geom_line(aes(group = asv_num, color = asv_num, alpha = ifelse(bloomer == 'bloomer', 0.8 , 0.3))) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer != 'bloomer'),
              aes(group = asv_num, color = NA, fill = asv_num), span = 0.09, alpha = 0.3, linewidth = 0.0) +
  geom_smooth(data = asv_tab_filtered |>
                dplyr::filter(bloomer == 'bloomer'), aes(group = asv_num, color = NA, fill = asv_num),linewidth = 0.0, span = 0.09) +
  scale_color_manual(values =  palette_asv42_asv49) +
  scale_fill_manual(values = palette_asv42_asv49) +
  facet_grid(vars(fraction), labeller = labs_fraction)+
  theme_bw() +
  guides(alpha = 'none')+
  theme(
    strip.background = element_blank(),
    legend.position = 'none',
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size =5),
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    legend.title  = element_text(size = 6),
    legend.text  = element_text(size = 4),
    panel.grid.major.y = element_blank()
  )

plot_asv42_asv49 

## plot the phylogenetic tree between those closely related taxa ----
# palette_asv42_asv49_v2 <- c(  "SAR86 clade asv1271"= "#4965f2",    "SAR86 clade asv87"   = "#384eba",
#                               "SAR86 clade asv49"   = "#3245a6",  
#                               "SAR86 clade asv1819" = "#1f2c69",   "SAR86 clade asv1942"   =  "#192354",
#                               "SAR86 clade asv5485"  = "#9d646c",  
#                               "SAR86 clade asv477" = "#649d6f",    "SAR86 clade asv679"  = "#5c6073",
#                               "SAR86 clade asv1051" = "#f2c949",  "SAR86 clade asv1021"   = "#7a4e54",
#                               "SAR86 clade asv329"  = "#e0909a", "SAR86 clade asv531"   = "#664146",
#                               "SAR86 clade asv10095" = "#90e0a0",  
#                               "SAR86 clade asv1294"  = "#416649", "SAR86 clade asv808"    = "#83cc91", 
#                               "SAR86 clade asv258"  = '#A7A4C4',   
#                               "SAR86 clade asv1004"  ="#9787f5", "SAR86 clade asv520"="#b3626e",
#                               "SAR86 clade asv414"   ="#6e62b3","SAR86 clade asv4588"  ="#64607f",
#                               "SAR86 clade asv132" ="#9db362",  "SAR86 clade asv124" ="#b3628e",   "SAR86 clade asv2516" ="#777f60", 
#                               "SAR86 clade asv670" ="#2a3312", "SAR86 clade asv6017"  ="#b36d62", "SAR86 clade asv1061" ="#b3a09e",
#                               "SAR86 clade asv281" = "#3f3866",   "SAR86 clade asv1822" = '#332C26',
#                               "SAR86 clade asv500" = '#B2211B',   
#                               "SAR86 clade asv42" = "#cfaa3e",
#                               "SAR86 clade asv2077" =  '#7A1612')

palette_asv42_asv49_v2 <- c(  "SAR86 clade asv1271"= "#C2AFB3",    "SAR86 clade asv87"   = "#C2AFB3",
                              "SAR86 clade asv49"   = "#3245a6",  "SAR86 clade asv42" = "#cfaa3e",
                              "SAR86 clade asv1819" = "#C2AFB3",   "SAR86 clade asv1942"   =  "#C2AFB3",
                              "SAR86 clade asv5485"  = "#C2AFB3",  
                              "SAR86 clade asv477" = "#C2AFB3",    "SAR86 clade asv679"  = "#C2AFB3",
                              "SAR86 clade asv1051" = "#C2AFB3",  "SAR86 clade asv1021"   = "#C2AFB3",
                              "SAR86 clade asv329"  = "#C2AFB3", "SAR86 clade asv531"   = "#C2AFB3",
                              "SAR86 clade asv10095" = "#C2AFB3",  
                              "SAR86 clade asv1294"  = "#C2AFB3", "SAR86 clade asv808"    = "#C2AFB3", 
                              "SAR86 clade asv258"  = '#C2AFB3',   
                              "SAR86 clade asv1004"  ="#C2AFB3", "SAR86 clade asv520"= "#C2AFB3",
                              "SAR86 clade asv414"   ="#C2AFB3","SAR86 clade asv4588"  ="#C2AFB3",
                              "SAR86 clade asv132" = "#C2AFB3",  "SAR86 clade asv124" ="#C2AFB3",   "SAR86 clade asv2516" ="#C2AFB3", 
                              "SAR86 clade asv670" = "#C2AFB3", "SAR86 clade asv6017"  ="#C2AFB3", "SAR86 clade asv1061" ="#b3a09e",
                              "SAR86 clade asv281" = "#C2AFB3",   "SAR86 clade asv1822" = '#C2AFB3',
                              "SAR86 clade asv500" = '#C2AFB3',   "SAR86 clade asv2077" =  '#C2AFB3')

tips_to_keep <- sar86_filter_tree$asv_num.y %in% tree_complete$tip.label

tree <- keep.tip(tree_complete, sar86_filter_tree$asv_num.y)

merged_data$genus_num

# Plot the modified tree
tree |>
  str()

merged_data <- merge(as_tibble_col(tree$tip.label, column_name = 'asv_num'), tax_bbmo_10y_new, by = "asv_num") |>
  as_tibble() |>
  dplyr::mutate(genus_num = paste0(family,' ', asv_num),
                asv_num_ed = asv_num)

tree_plot_asv42_asv49 <- ggtree(tree, branch.length = 'none') %<+%   
  merged_data +   
  geom_tiplab(hjust = -0.25, size = 2.5, align = TRUE) +   
  geom_tippoint(aes(color = genus_num), alpha = 0.8, size = 3) +  # Use color for points
  scale_color_manual(values = palette_asv42_asv49_v2) +   # Corrected to scale color instead of fill
  labs(
    title = paste0(unique(merged_data$order),' ',
                   unique(merged_data$family)),
    color = 'ASV nº'
  ) + 
  guides(color = guide_legend(ncol = 4))+
  theme_tree2() +   
  theme(
    legend.position = "none",  # Show legend
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),  # Remove grid lines
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    title = element_text(size = 6),
    plot.margin = margin(2.5, 5, 5, 2.5)
  )+
  xlim(0, 20)

tree_plot_asv42_asv49

## create a composition with the time series and ----

closely_related_bloomers_sar86clade_plot  <- plot_grid(
  plot_asv42_asv49,
  tree_plot_asv42_asv49,
  ncol = 2,                # One column layout for the main grid
  rel_heights = c(1),
  labels = c('A', 'B'), label_fontface = 'plain',
  label_size = 10,
  align = "hv",           # Align plots horizontally and vertically
  axis = "tblr",          # Display axis labels on all sides
  rel_widths = c(2, 1)    # Set relative widths of columns
)

# ggsave( plot = closely_related_bloomers_sar86clade_plot,
#         filename = 'closely_related_bloomers_sar86clade_plot_V2.pdf',
#         path = 'results/figures/',
#         width = 180, height = 150, units = 'mm')


## what happens with those ASVs that do not have any closely related taxa? Are they driven by different factors? ----

bloo_taxonomy |>
  dplyr::filter(!asv_num %in% unique(closely_related_bloomers_and_others$asv_num)) |>
  left_join(bloo_all_types_summary_tb_v2) ## 45

## transform the previous plots to heatmaps -----
### i need to create a new variable to facet_wrap which is colse taxa

### dataset with all closely related taxa group ----
bloo_02_filt <- bloo_02 |>
  dplyr::filter(value != c('asv2', 'asv3', 'asv8', 'asv5'))

asv_num <- bloo_02_filt$value[[1]]

  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- temp_tb
  
  asv_num <- bloo_02_filt$value[[2]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 

  asv_num <- bloo_02_filt$value[[3]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[4]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_02_filt$value[[5]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[6]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 

  asv_num <- bloo_02_filt$value[[7]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 

  asv_num <- bloo_02_filt$value[[8]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[9]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  asv_num <- bloo_02_filt$value[[10]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[11]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[12]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[13]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[14]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[15]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_02_filt$value[[16]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  # the same for those ASVs that appear only in PA
  bloo_3_exclusive <- bloo_3 |>
    anti_join(bloo_02_filt)
  
  asv_num <- bloo_3_exclusive$value[[1]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb) 
  
  asv_num <- bloo_3_exclusive$value[[2]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[3]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[4]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[5]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[6]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[7]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[8]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[9]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[10]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[11]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[12]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[13]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[14]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[15]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[16]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[17]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[18]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[19]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[20]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[21]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[22]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[23]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[24]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[25]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[26]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[27]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[28]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[29]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[30]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
  asv_num <- bloo_3_exclusive$value[[31]]
  
  temp_tb <- phylogenetic_distances_tb_com |>
    dplyr::filter(asv_num_1 == asv_num) |>
    dplyr::filter(asv_num_2!= asv_num) |>
    dplyr::filter(phylogenetic_distance < 0.012) |>
    pivot_longer(cols = starts_with('asv_num')) %$%
    unique(value) |>
    as_tibble_col(column_name = 'asv_num_var') |>
    dplyr::mutate(close_taxa_group = paste0('close_taxa_', asv_num))
  
  close_taxa_tb <- close_taxa_tb |>
    bind_rows(temp_tb)
  
## heat map with all closely related taxa groups ----
close_taxa_asv_num <- close_taxa_tb |>
    distinct(asv_num_var)

asv_tab_close_tax_tb_02 <- asv_tab_10y_02_rclr |>
  dplyr::filter(asv_num %in% close_taxa_asv_num$asv_num_var) |>
  full_join(close_taxa_tb, by = c('asv_num' = 'asv_num_var')) |>
  left_join(m_02) |>
  dplyr::select(date, sample_id_num, asv_num, close_taxa_group, rclr, fraction, order_f)

asv_tab_close_tax_tb_3 <- asv_tab_10y_3_rclr |>
  dplyr::filter(asv_num %in% close_taxa_asv_num$asv_num_var) |>
  full_join(close_taxa_tb, by = c('asv_num' = 'asv_num_var')) |>
  left_join(m_3) |>
  dplyr::select(date, sample_id_num, asv_num, close_taxa_group, rclr, fraction, order_f)

asv_tab_close_tax_tb_all <- asv_tab_close_tax_tb_02 |>
  bind_rows(asv_tab_close_tax_tb_3)

#write.csv(asv_tab_close_tax_tb_all , 'data/asv_tab_close_tax_tb_all.csv')

asv_tab_close_tax_tb$rclr |>
  range()

closely_related_groups_plot <- asv_tab_close_tax_tb_all |>
  ggplot(aes(sample_id_num, interaction(asv_num, fraction))) +
  geom_tile(aes(fill = rclr))+
  facet_wrap(vars(close_taxa_group), scales = 'free', ncol = 2)+
  #scale_fill_gradientn(colours = palette_gradient_bw)+
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-4.00, 6.6), na.value = '#D7D6D3')+
  theme_bw() +
  theme(
    panel.border = element_blank(),
    strip.background = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    text = element_text(size = 4)
  )

closely_related_groups_plot 
# 
# ggsave('closely_related_groups_plot.pdf', plot = closely_related_groups_plot, 
#        path = 'Results/Figures/',
#        width = 180, height = 500, units = 'mm')

## i do it with them individually ---
# Get unique values of close_taxa_group
groups <- unique(asv_tab_close_tax_tb_all$close_taxa_group)

# Loop over each group and create a plot
for (group in groups) {
  # Filter the data for the current group
  group_data <- asv_tab_close_tax_tb_all |>
    dplyr::filter(close_taxa_group == group)
  
  # Create the plot
  plot <- group_data |> 
    ggplot(aes(as.numeric(sample_id_num), asv_num)) + 
    geom_tile(aes(fill = rclr)) + 
    scale_x_continuous()+
    scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-4.00, 6.6), na.value = '#D7D6D3') + 
    facet_wrap(vars(fraction), labeller = labs_fraction, ncol = 1)+
    labs(x = 'Time', y = 'Closely Related ASVs', fill = 'rCLR', title = unique(group_data$close_taxa_group) )+
      theme_bw() + 
    theme(
      axis.text.x = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      axis.title  = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      text = element_text(size = 14),
      title = element_text(size = 10)
    )
  
  # Save the plot to a file
  ggsave(filename = paste0("plot_", group, ".pdf"),  path = 'Results/Figures/',
         plot = plot, width = 180, height = 120, units = 'mm')
}


## heatmap complete

### remove duplicated groups
 asv_tab_close_tax_tb_all |>
   dplyr::filter(close_taxa_group %in% c('close_taxa_asv17' , 'close_taxa_asv77')) |>
  distinct(asv_num, close_taxa_group) |>
   arrange(close_taxa_group) ## not completely the same
 
 asv_tab_close_tax_tb_all |>
   dplyr::filter(close_taxa_group %in% c('close_taxa_asv4' , 'close_taxa_asv31')) |>
   distinct(asv_num, close_taxa_group) |>
   arrange(close_taxa_group) ## completely the same
 
 asv_tab_close_tax_tb_all |>
   colnames()

closely_related_groups_plot <- asv_tab_close_tax_tb_all |>
  left_join(bloo_taxonomy, by = 'asv_num') |>
  dplyr::filter(close_taxa_group !=  'close_taxa_asv31') |>
  ggplot(aes(sample_id_num, asv_num)) +
  geom_tile(aes(fill = rclr))+
  scale_y_discrete()+
  labs(x = 'Time', y = 'Closely Related ASVs', fill = 'rCLR')+
  facet_grid(close_taxa_group~fraction, scales = 'free_y', labeller = labs_fraction)+
  #scale_fill_gradientn(colours = palette_gradient_bw)+
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-4.00, 6.6), na.value = '#D7D6D3')+
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    strip.background = element_rect(aes(fill = order_f)),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    axis.title  = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    text = element_text(size = 12),
    strip.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.key.size = unit(0.5, "lines")
  )+
  scale_fill_manual(name = "Order", values = palette_order_assigned_bloo)

closely_related_groups_plot 
# 
# ggsave('closely_related_groups_plot_ed1.pdf', plot = closely_related_groups_plot,
#        path = 'Results/Figures/',
#        width = 180, height = 380, units = 'mm')
