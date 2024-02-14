
## I need to create a fasta file with the sequences of my potential bloomers----

# Aligment
##https://www.genome.jp/tools-bin/clustalw

## add fasta seqs, DNA format, and Select Weight Matrix: CLUSTALW (for DNA)

 asv_num_seq_bloo <- tax |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::select(seq, asv_num)


# Assuming your tibble is named 'sequences_tibble' with columns 'sequence_name' and 'sequence'

# Open a connection to write the output to a file
output_file <- "sequences.fasta"
file_conn <- file(output_file, "w")

# Iterate over each row in the tibble and write to the file in FASTA format
for (i in 1:nrow( asv_num_seq_bloo)) {
  # Extract sequence name and sequence
  seq_name <-  asv_num_seq_bloo[i, "asv_num"]
  seq <-  asv_num_seq_bloo[i, "seq"]
  
  # Write to file in FASTA format
  cat(">", seq_name, "\n", seq, "\n", file = file_conn, sep = "")
}

# Close the file connection
close(file_conn)

# Print a message indicating the file has been written
cat("FASTA file '", output_file, "' has been created.\n")

## Run thi code in marbits or in local look at README file in the folder----

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
 
 palette_fraction <- c(  '0' = 'white', '1' = 'black' )

 # Plot the dendrogram with colored branches based on 'class'
 tree_plot <- ggtree(tree, branch.length='none') %<+% ## this is used to assign new data to the plot
  merged_data +
  geom_tree() + #aes(color=class)
  #scale_color_manual(values = palette_class_assigned) +
   geom_tiplab(aes(label = family_num), size=2, align=TRUE) + 
  labs(color = 'Class') +
  theme_tree2()+
   theme(legend.position = "none",  # Remove legend
         panel.grid.major = element_blank(),  # Remove grid lines
         panel.grid.minor = element_blank(),  # Remove grid lines
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.line.x = element_blank(),
         axis.ticks.x = element_blank()) 


# Create the heatmap colored by 'detect_both' variable
heatmap <- gheatmap(tree_plot, detect_both, offset=7, width = .1, color=NA, font.size = 2) + 
  scale_fill_manual(values = palette_fraction, name="Potential\nblooming\nin fraction", labels=c("0" = "No", "1" = "Yes"))+
  theme(legend.position = 'bottom', text = element_text(size = 4),
        legend.key.size = unit(0.5, "lines"),
        )

print(heatmap)
