
library(tidyverse)

## I need to create a fasta file with the sequences of my potential bloomers----

# Aligment
##https://www.genome.jp/tools-bin/clustalw

## add fasta seqs, DNA format, and Select Weight Matrix: CLUSTALW (for DNA)

 asv_num_seq_bloo <- tax |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::select(seq, asv_num) |>
  dplyr::mutate(seq = as.character(seq),
                asv_num = as.character(asv_num))

# Assuming your tibble is named 'sequences_tibble' with columns 'sequence_name' and 'sequence'

# Open a connection to write the output to a file
output_file <- "sequences.fasta"
file_conn <- file(output_file, "w")

# Iterate over each row in the tibble and write to the file in FASTA format
for (i in 1:nrow(asv_num_seq_bloo)) {
  # Extract sequence name and sequence
  seq_name <- as.character(asv_num_seq_bloo$asv_num[i])  # Extract character string from the tibble
  seq <- as.character(asv_num_seq_bloo$seq[i])  # Ensure seq is a character vector
  
  # Write to file in FASTA format
  cat(">", seq_name, "\n", seq, "\n", file = file_conn, sep = "")
}

# Close the file connection
close(file_conn)

# Print a message indicating the file has been written
cat("FASTA file '", output_file, "' has been created.\n")

## Run the code in marbits or in local look at README file in the folder----

# module load
# 
# module load raxml-ng/0.9.0
# 
# raxml-ng --check --msa /alignment_remei_100_long_trimmed.fasta --model GTR+I+G -seed 123 --prefix remei100 ## for checking sequences 
# 
# raxml-ng --msa alignment_remei_100_long_trimmed.fasta --model GTR+I+G —-outgroup NR_074309.1_Synechococcus_elongatus_PCC_6301 -seed 123 --threads 2 --prefix remei100 ## for doing the main tree
# 
# raxml-ng --rfdist --tree remei100.raxml.mlTrees --prefix remei100 ## para comprobar topologias
# 
# raxml-ng --bootstrap --msa alignment_remei_100_long_trimmed.fasta --model GTR+I+G --outgroup NR_074309.1_Synechococcus_elongatus_PCC_6301 -seed 123 --threads 2 --bs-trees 600 --prefix remei100 ## for making bootstraps
# 
# raxml-ng --bsconverge --bs-trees remei100.raxml.bootstraps --prefix remei100 -seed 123 --threads 1 --bs-cutoff 0.01 ## para ver si las bootstraps convergen
# 
# raxml-ng --support --tree remei100.raxml.bestTree --bs-trees remei100.raxml.bootstraps --prefix remei100 --threads 2 ## para juntar las bootstraps con el tree


## Visualization of our tree----
#devtools::install_github("GuangchuangYu/ggtree")
library(ggtree)
library(ggplot)
library(tidyverse)

tree <- ggtree::read.tree('data/raxml/bloo_bbmo.raxml.support')

tree |>
  ggplot() + 
  geom_tree() + 
  theme_tree2() +
  #geom_treescale()+
  geom_tiplab(align = T)
  #geom_tippoint()

# I try to add a column with information of the types of blooms that they form (from the wavelets analysis). 


# And another column with the information of in which fraction do they bloom.
tree |>
  ggplot() + 
  geom_tree() + 
  theme_tree2() +
  #geom_treescale()+
  geom_tiplab(align = T)

tax <- asv_tab_all_bloo_z_tax |>
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

merged_data <- merge(as_tibble_col(tree$tip.label, column_name = 'asv_num'), metadata, by = "asv_num") |>
  dplyr::mutate(family_num = paste0(family,' ', asv_num))

 tree_plot <- ggtree(tree, layout='dendrogram', branch.length='none') %<+% ## this is used to assign new data to the plot
  merged_data +
  geom_tree(aes(color=class))+
  scale_color_manual(values = palette_class_assigned)+
  geom_tiplab( aes(label=family), size=2, align=TRUE) +
    labs(color = 'Class')+
  #geom_treescale()+
  #geom_tiplab(align = T)+
  theme_tree2()
  #geom_treescale(x = 10, y = 360)
 
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

#### Based on wavelets analysis I create a label for each taxa and add the information to the phylogenetic tree----
bloo_type_biased_all ##here I have their maximum coefficient. Maaybe the lowest ones are not significant and other criteria need to be followed to decide what do we trust

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

ggsave(heatmap1, filename = 'phylogenetic_tree_wavelets_coeff.pdf',
       path = 'Results/Figures/',
       width = 230, height = 200, units = 'mm')

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

ggsave(heatmap1, filename = 'phylogenetic_tree_wavelets_coeff_red.pdf',
       path = 'Results/Figures/',
       width = 230, height = 200, units = 'mm')

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
                       name="Bloomer's\nwavelets\ncoefficients\nmagnitude")+
  
  #geom_text(data ="Particle attached (3-20 um)")+
  #annotate(geom = "Particle attached (3-20 um)", x = 0.5,   hjust = 0.5) +  # Adjust x position as needed
  theme(legend.position = 'bottom', text = element_text(size = 5), # labels=c("0" = "No", "1" = "Yes")
        legend.key.size = unit(0.5, "lines"))

# Print the combined plot
print(heatmap1)

ggsave(heatmap1, filename = 'phylogenetic_tree_wavelets_coeff_variance.pdf',
       path = 'Results/Figures/',
       width = 230, height = 200, units = 'mm')

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




