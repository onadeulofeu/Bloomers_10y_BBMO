library(seqinr)
library(ape)
library(bio3d)


library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)

palette_fraction <- c('0.2' = '#00808F', '3' = '#454545')

## upload data
asv_tab_all_bloo_z_tax <- read.csv('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') ## bloomers table

asv_tab_10y_02_pseudo_rclr <- read.csv('data/asv_tab_10y_02_pseudo_rclr.csv')  ## community table
asv_tab_10y_3_pseudo_rclr <- read.csv('data/asv_tab_10y_3_pseudo_rclr.csv')  ## community table


# FL and PA taxa in the whole community and in the bloomer's community-----
## I will use a Venn Diagram to illustrate the taxa that are exclusive of one fraction or in both 
### for bloomer's commmunity----

df <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::select(asv_num, fraction) |>
  distinct(asv_num, fraction) |>
  dplyr::mutate(value = 1) |>
  pivot_wider(id_cols = asv_num, names_from = fraction, values_from = value, values_fill = 0)
  
# Assuming your data is stored in a data frame called 'df'
set_0.2 <- df$asv_num[df$`0.2` == 1]
set_3 <- df$asv_num[df$`3` == 1]
set_shared <- intersect(set_0.2, set_3)

# Create a list of sets
sets <- list(`0.2` = set_0.2, `3` = set_3)

# Create the Venn diagram
venn_diagram <- venn.diagram(
  x = list(`0.2` = set_0.2, `3` = set_3),
  category.names = c("Free living\n(0.2-3 um)", "Particle\nattached\n(3-20 um)"),
  fill = c(alpha("#00808F",0.3), alpha('#454545',0.3)),
  col = c(alpha("#00808F",0.3), alpha('#454545',0.3)),
  cex = 2,
  cat.cex = 2,
  fontface = 1,
  cat.dist = 0.05,
  cat.pos = 5,  # Adjust category name position as needed
  filename = NULL)

# Save the diagram to a file
png("venn_diagram.png", width = 800, height = 600)  # Adjust width and height as needed
grid.draw(venn_diagram)
dev.off()

### whole_community----
asv_tab_10y_02_rel_02 <- asv_tab_10y_02_rel |>
  dplyr::mutate(fraction = '0.2')

asv_tab_10y_3_rel_3 <- asv_tab_10y_3_rel |>
  dplyr::mutate(fraction = '3')

df <- asv_tab_10y_02_rel_02 |>
  bind_rows(asv_tab_10y_3_rel_3) |>
  dplyr::filter(relative_abundance > 0) |>
  ungroup() |>
  dplyr::select(-sample_id) |>
  dplyr::select(asv_num, fraction) |>
  distinct(asv_num, fraction) |>
  dplyr::mutate(value = 1) |>
  pivot_wider(id_cols = asv_num, names_from = fraction, values_from = value, values_fill = 0)

# Assuming your data is stored in a data frame called 'df'
set_0.2 <- df$asv_num[df$`0.2` == 1]
set_3 <- df$asv_num[df$`3` == 1]
set_shared <- intersect(set_0.2, set_3)

# Create a list of sets
sets <- list(`0.2` = set_0.2, `3` = set_3)

# Create the Venn diagram
venn_diagram <- venn.diagram(
  x = list(`0.2` = set_0.2, `3` = set_3),
  category.names = c("Free living\n(0.2-3 um)", "Particle\nattached\n(3-20 um)"),
  fill = c(alpha("#00808F",0.3), alpha('#454545',0.3)),
  col = c(alpha("#00808F",0.3), alpha('#454545',0.3)),
  cex = 2,
  cat.cex = 2,
  fontface = 1,
  cat.dist = 0.05,
  cat.pos = 5,  # Adjust category name position as needed
  filename = NULL
)

# Save the diagram to a file
png("venn_diagram_whole_community.png", width = 800, height = 600)  # Adjust width and height as needed
grid.draw(venn_diagram)
dev.off()


# Read FASTA sequences------
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

### look Marta's script differential_abundance_corncob_forONA.Rmd

library(phyloseq)
library(readxl)
library(phangorn)
library(Biostrings)
library(DECIPHER)

#read all data 
#read all data 
asv_tax<-read.csv("data/MALASPINA_16S_VP_amplicons_dada2_from_megarun_4R_NoChloroplastMitochondria_min2000.txt",sep="\t",header=TRUE,row.names = 1,stringsAsFactors = TRUE)
env<-read.table("data/env_vertical_profiles_forR_all.txt",sep="\t", header=TRUE, check.names = FALSE)

colnames(asv_tax)
tax<-asv_tax[,184:190]
asv<-asv_tax[,colnames(asv_tax)%in%env$sample_id]

colnames(asv)


## PHYLOGENETIC DISTANCE MARTA SCRIPT------
#load always in *this* order:
library(plyr)
library(tidyr)
library(dplyr)
library(tidyverse)

seqs<-asv_tax$ASV
seqs<-as.data.frame(seqs)
rownames(asv) <- str_c( 'asv' , 1:nrow(asv))
rownames(seqs)<-rownames(asv)
fasta<-seqs
rm(seqs)
fasta<-rownames_to_column(fasta)
colnames(fasta)[1]<-"asv.id"
asv$ASV<-NULL

samdata <-env[match(colnames(asv), env$sample_id),]
asv.table<-asv

# lets do the tree

seqs <- DNAStringSet(fasta$seqs)
names(seqs) <- fasta$asv.id
alignment <- AlignSeqs(seqs, anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)



## negative edges length changed to 0!
## this takes a loooong time...
fitGTR <- update(fit, k=4, inv=0.2)

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

# now we have to integrate all the data.frames/matrices into a single object 

saveRDS(fitGTR,"output/fitGTR.rds")

OTU <- otu_table(as.matrix(asv.table), taxa_are_rows = TRUE)

#next is the taxonomy 
# i will create a new column called seq 

all(rownames(tax)==fasta$asv.id)
tax$seq<-fasta$seqs

TAX <- tax_table(as.matrix(tax))

#Samdata same procedure 
rownames(samdata) <- samdata$sample_id

DAT <- sample_data(samdata)


#We join all the data together 
Profilesphy <- phyloseq(OTU, TAX,DAT, phy_tree(fitGTR$tree))

saveRDS(object = Profilesphy, file = "data/Profiles_phyloseq_clean.rds")