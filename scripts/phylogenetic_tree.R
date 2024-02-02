
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

## Después corremos este trozo de código en marbits o en local look at README file in the folder----

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


## Visualization of our tree
#devtools::install_github("GuangchuangYu/ggtree")
library(ggtree)
library(ggtree)

tree <- ggtree::read.tree('raxml/bloo_bbmo.raxml.support')

tree |>
  ggplot() + 
  geom_tree() + 
  theme_tree2() +
  #geom_treescale()+
  geom_tiplab(align = T)

  #geom_tippoint()
  


