library(seqinr)
library(ape)
library(bio3d)

## upload data

asv_tab_all_bloo_z_tax <- read.csv('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') ## bloomers table

asv_tab_10y_02_pseudo_rclr <- read.csv('data/asv_tab_10y_02_pseudo_rclr.csv')  ## community table
asv_tab_10y_3_pseudo_rclr <- read.csv('data/asv_tab_10y_3_pseudo_rclr.csv')  ## community table

# Read FASTA sequences
sequences <- read.fasta(file = "data/asv_seqs_bbmo10y.fasta")

# Perform multiple sequence alignment
alignment <- seqaln(sequences, method = "ebi", email = "odeulofeu@icm.csic.es")

# Calculate distance matrix
dist_matrix <- dist.alignment(alignment)

# Build phylogenetic tree
tree <- nj(dist_matrix)

# Visualize the tree (optional)
plot(tree)

# Identify closely related sequences (you may need to adjust the threshold)
closely_related <- ape::cophenetic(tree) < threshold

## When the competition increases (more taxa phylogenetically related are present in the dataset) is it difficult for them to bloom? ----
### This hypothesis will be confirmed in case that the community has higher influence on the stochastic blooms than the environment

## I test it for ASV11 and see

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num == 'asv11') |>
  dplyr::select(family, genus)

asv_tab_all_bloo_z_tax |>
  dplyr::filter(family == 'Alteromonadaceae') |>
  distinct(asv_num)
