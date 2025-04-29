# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## packages ---- 
library(tidyverse)
library(fields) ##add the legend
library(magrittr)
library(ggdendro)
library(gridExtra) ## combine the dendogram and the heatmap
library(factoextra) ## visualize hierarchical clusters

## palettes----
palette_clustering <- c("#FFE355",
                        "#EF8D00",
                        "#AE659B",
                        "#2b347a",
                        "#c23939",
                        "#4cb76a",
                        "#804c90",
                        "#006cb0",
                        "#518535",
                        "#2d373b",
                        "#6a6964",
                        "#57a9a8",
                        "#c9acb8",
                        "#8289ba",
                        "#bec735",
                        "#c5ebff",
                        "#ca6094"
                        )

palete_seasonal_bloo <- c('#2466BF', '#5B9B57', '#FFD700', '#C73F4E', '#AD7CE6')

palette_clustering_assigned <- c('cl_9_3' ="#FFD700",
                                 "#EF8D00",
                                 "#AE659B",
                                 "#2b347a",
                                 "#c23939",
                                 'cl_8_3' = "#4cb76a",
                                 "#804c90",
                                 'cl_5_3'= '#2466BF',
                                 'cl_3_0.2' = '#2466BF', 
                                 'cl_5_3'= "#518535",
                                 "#2d373b",
                                 "#6a6964",
                                 "#57a9a8",
                                 "#c9acb8",
                                 "#8289ba",
                                 "#bec735",
                                 "#c5ebff",
                                 "#ca6094",
                                 'asv27_3'=  '#C73F4E' , 
                                 'asv28_3'= '#AD7CE6',
                                 'asv27_0.2'=  '#C73F4E' , 
                                
                                 'asv62_0.2' = '#5B9B57' , 
                                 'asv1_0.2'= '#2F2F2C')

palette_occurrence <- c(narrow = "#AE659B",
                        intermidiate = "#3e3e3e",
                        broad = "#57a9a8")

## List of my clusters and which ASVs belong in there ----
## (cl_1 = 'habor restoration 1' ('asv311', 'asv302', 'asv511')
#                                   cl_2 = 'habor restoration 2', c('asv194', 'asv105', 'asv559')
#                                   cl_3 = 'habor restoration 3', c('asv22', 'asv85', 'asv163', 'asv219', 'asv80', 'asv192', 'asv49')
#                                   cl_4 = 'harbor restoration 4', c('asv276', 'asv264', 'asv223', 'asv471', 'asv752')
#                                   cl_5 = 'seasonal 1', c('asv7', 'asv15')
#                                   cl_6 = 'harbor restoration 5', c('asv317', 'asv200', 'asv113')
#                                   cl_7 = 'recurrent random', c('asv116', 'asv182', 'asv84')
#                                   cl_8 = 'seasonal 2', c('asv100', 'asv25', 'asv72', 'asv42')
#                                   cl_9 = 'seasonal 3', c('asv23', 'asv1', 'asv31', 'asv4')
#                                   cl_10 = 'harbor restoration 6', c('asv69', 'asv225')
#                                   cl_11 = 'harbor restoration 7', c('asv153', 'asv77')
#                                   unclear = 'ungruped' asv179, asv11, asv385, asv27, asv17, asv43, asv118, asv126, asv28
# cl_1_0.2 = 'recurrent random', c('asv58', 'asv178')
# cl_2_0.2 = 'ephemeral random', c('asv555', 'asv114', 'asv249', 'asv237', 'asv563', 'asv282')
# cl_3_0.2 = 'seasonal 1', c('asv15', 'asv7')
# cl_4_0.2 = 'SAR11 cluster', c('asv2', 'asv3', 'asv5', 'asv8') 
# unclear_0.2 = 'ungrouped' 'asv11' 'asv27' 'asv38' 'asv1' 'asv17' 'asv62'

palete_seasonal_bloo <- c('cl_5_3'= '#2466BF', 'cl_8_3' =  '#5B9B57','cl_9_3' = '#FFD700', 'asv27_3'=  '#C73F4E' , 'asv28_3'= '#AD7CE6',
                          'asv27_0.2'=  '#C73F4E' , 'cl_3_0.2' = '#2466BF', 'asv62_0.2' = '#5B9B57' , 'asv1_0.2'= '#2F2F2C')

labs_clusters_pa_fl <-   as_labeller(c(cl_1_3 = 'habor restoration 1', ## clusters meaning
                                       cl_2_3 = 'habor restoration 2',
                                       cl_3_3 = 'habor restoration 3',
                                       cl_4_3 = 'harbor restoration 4',
                                       cl_5_3 = 'seasonal 1',
                                       cl_6_3 = 'harbor restoration 5', 
                                       cl_7_3 = 'recurrent random',
                                       cl_8_3 = 'seasonal 2',
                                       cl_9_3 = 'seasonal 3',
                                       cl_10_3 = 'harbor restoration 6',
                                       cl_11_3 = 'harbor restoration 7',
                                       unclear_3 = 'ungruped', 
                                       cl_1_0.2 = 'recurrent random', 
                                       cl_2_0.2 = 'ephemeral random',
                                       cl_3_0.2 = 'seasonal 1',
                                       cl_4_0.2 = 'SAR11 cluster',
                                       unclear_0.2 = 'ungrouped'))

## labels----
year_labels <- c("2004", "2005", "2006", '2007', '2008', '2009', '2010',
                 '2011', '2012', '2013') ## create the labels with correspond with the sample num


labs_occurrence <- as_labeller(c(narrow = "Narrow (<1/3)",
                                 intermidiate = "Intermediate (1/3-2/3)",
                                 broad = "Broad (>2/3)"))


## upload some datasets that I will need for the plots----
##upload occurrence data -----
occurrence_bloo_bbmo <- read.delim2('data/occurrence_bloo_bbmo.csv', sep = ',')

occurrence_bloo_bbmo |>
  head()

occurence_perc_3 <- occurrence_bloo_bbmo |>
  dplyr::filter(fraction == '3') |>
  dplyr::select(asv_num, fraction, occurrence_perc) |>
  distinct(asv_num, fraction, occurrence_perc)

occurence_perc_02 <- occurrence_bloo_bbmo |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::select(asv_num, fraction, occurrence_perc) |>
  distinct(asv_num, fraction, occurrence_perc) 

## add occurrence categories 
occurrence_perc_3 <- occurence_perc_3 |>
dplyr::mutate(occurrence_category = ifelse(occurrence_perc > 2/3, 'broad',
                                           ifelse(occurrence_perc < 1/3, 'narrow',
                                                  'intermediate')))

occurrence_perc_02 <- occurence_perc_02 |>
  dplyr::mutate(occurrence_category = ifelse(occurrence_perc > 2/3, 'broad',
                                             ifelse(occurrence_perc < 1/3, 'narrow',
                                                    'intermediate')))

occurrence_bloo_bbmo <- occurence_perc_3 |>
  bind_rows(occurrence_perc_02)

bloo_all_types_summary_tb <- read.csv('results/tables/summary_types_of_blooms.csv')

bloo_all_types_summary_tb_tax <- bloo_all_types_summary_tb |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::select(-seq) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8'))

#write.csv(bloo_all_types_summary_tb_tax, 'results/tables/bloo_all_types_summary_tb_tax.csv')

summary_types_blooms_order <- bloo_all_types_summary_tb_tax |>
  group_by(fraction, occurrence_category, recurrency, order) |>
  reframe(n = n())

bloo_all_types_summary_tb_tax |>
  group_by(fraction, recurrency) |>
  reframe(n = n())

#write.csv(summary_types_blooms_order, 'results/tables/summary_types_blooms_order.csv')

bloo_all_types_summary_tb_tax |>
  distinct(asv_num)

bloo_all_types_summary_tb_tax |>
  group_by(asv_num, fraction) |>
  reframe(n = n()) |>
  group_by(fraction) |>
  reframe(sun = sum(n))

bloo_all_types_summary_tb_tax |>
  group_by(asv_num) |>
  filter(n() > 1) |>
  distinct(asv_num)

### (I'm not using this approximation)------
## Fuzzy C-Means
### FCM is a soft clustering alogrithm proposed by Bezdek. Unlike K-means algorithm in which each data object is the
### the member of only one cluster, a data object is the member of all clusters with varyinh degrees of fuzzy membership
### between 0 and 1 in FCM. Hence, the data objects closer to the centers of clusters have higher degrees of membership 
### than objects scattered in the borders of clusters

# library(ppclust)

## a hierarchical clustering method is one which works by partitioning the data into groups with
## increasing similar features.

## prior to clustering we need to decide the adequate number of clusters we need for our dataset:
## for this purpose we not to implement a clustering algorithm using a range of possible numbers 
## of cluster, and then comparison of these indices will indicate which number has a high degree
## of fit without over-fitting.
## Although there are many internal indexes that have originally been proposed
## for working with hard membership degrees preoduced by the K-means and its variants,
## most of these indexes cannot be used for fuzzy clustering results.

# asv_tab_all_bloo_z_tax |>
#   colnames()
# 
# # I want to select only the bloomers that bloom en each particular fraction
# bloo_3 # vector with ASV number of potential bloomers in PA fraction
# bloo_02 # vector with ASV number of potential bloomers in FL fraction
# 
# bloo_3_w <- asv_tab_all_bloo_z_tax |>
#   dplyr::filter(fraction == '3' &
#                   asv_num %in% bloo_3) |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   dplyr::select(asv_num, sample_id, abundance_value) |>
#   pivot_wider(id_cols = sample_id, values_from = abundance_value, names_from = asv_num) |>
#   as.data.frame(row.names = sample_id)
# 
# bloo_02_w <- asv_tab_all_bloo_z_tax |>
#   dplyr::filter(fraction == '0.2' &
#                   asv_num %in% bloo_02) |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   dplyr::select(asv_num, sample_id, abundance_value) |>
#   pivot_wider(id_cols = sample_id, values_from = abundance_value, names_from = asv_num) |>
#   as.data.frame(row.names = sample_id)
# 
# # test |>
# #   class()
# # 
# # test |>
# #   glimpse()
# 
# par(mar = c(1, 1, 1, 1)) # we use this to fix the error with margins too large. 
# 
# cluster_3_15c <- fcm(bloo_3_w[,-1], centers = 15, m = 4)
# cluster_3_10c <- fcm(bloo_3_w[,-1], centers = 10)
# cluster_3_5c <- fcm(bloo_3_w[,-1], centers = 5)
# cluster_3_4c <- fcm(bloo_3_w[,-1], centers = 4)
# cluster_3_3c <- fcm(bloo_3_w[,-1], centers = 3)
# 
# summary(cluster_3_15c)
# summary(cluster_3_10c)
# summary(cluster_3_5c)
# summary(cluster_3_4c)
# summary(cluster_3_3c)
# 
# plotcluster(cluster_3_15c, cp=5)
# plotcluster(cluster_3_10c, cp=1)
# plotcluster(cluster_3_5c, cp=1)
# plotcluster(cluster_3_4c, cp=1)
# plotcluster(cluster_3_3c, cp=1)
# 
# res.fcm2_15c <- ppclust2(cluster_3_15c, "kmeans")
# res.fcm2_10c <- ppclust2(cluster_3_10c, "kmeans")
# res.fcm2_5c <- ppclust2(cluster_3_5c, "kmeans")
# res.fcm2_4c <- ppclust2(cluster_3_4c, "kmeans")
# res.fcm2_3c <- ppclust2(cluster_3_3c, "kmeans")
# 
# factoextra::fviz_cluster(res.fcm2_15c, data = bloo_3_w[,-1], 
#                          ellipse.type = "convex",
#                          palette = "jco",
#                          repel = TRUE)
# 
# factoextra::fviz_cluster(res.fcm2_10c, data = bloo_3_w[,-1], 
#                          ellipse.type = "convex",
#                          palette = "jco",
#                          repel = TRUE)
# 
# factoextra::fviz_cluster(res.fcm2_5c, data = bloo_3_w[,-1], 
#                          ellipse.type = "convex",
#                          palette = "jco",
#                          repel = TRUE)
# 
# factoextra::fviz_cluster(res.fcm2_4c, data = bloo_3_w[,-1], 
#                          ellipse.type = "convex",
#                          palette = "jco",
#                          repel = TRUE)
# 
# factoextra::fviz_cluster(res.fcm2_3c, data = bloo_3_w[,-1], 
#                          ellipse.type = "convex",
#                          palette = "jco",
#                          repel = TRUE)
# 
# library(cluster)
# res.fcm3_15 <- ppclust2(cluster_3_15c, "fanny")
# res.fcm3_10 <- ppclust2(cluster_3_10c, "fanny")
# res.fcm3_5 <- ppclust2(cluster_3_5c, "fanny")
# res.fcm3_4 <- ppclust2(cluster_3_4c, "fanny")
# res.fcm3_3 <- ppclust2(cluster_3_3c, "fanny")
# 
# # cluster::clusplot(scale(test[,-1]), res.fcm3$cluster,  
# #                   main = "Cluster plot of Bloomers in the PA fraction",
# #                   color=TRUE, labels = 2, lines = 2, cex=1)
# 
# cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_15$cluster,  
#                   main = "Cluster plot of Bloomers in the PA fraction",
#                   color=TRUE, labels = 2, lines = 2, cex=1)
# 
# cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_10$cluster,  
#                   main = "Cluster plot of Bloomers in the PA fraction",
#                   color=TRUE, labels = 2, lines = 2, cex=1)
# 
# cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_5$cluster,  
#                   main = "Cluster plot of Bloomers in the PA fraction",
#                   color=TRUE, labels = 2, lines = 2, cex=1)
# 
# cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_4$cluster,  
#                   main = "Cluster plot of Bloomers in the PA fraction",
#                   color=TRUE, labels = 2, lines = 2, cex=1)
# 
# cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_3$cluster,  
#                   main = "Cluster plot of Bloomers in the PA fraction",
#                   color=TRUE, labels = 2, lines = 2, cex=1)
# 
# ## Validation of the clustering results find the best clustering option
# library(fclust)
# 
# ## 15 clusters
# res.fcm4 <- ppclust2(cluster_3_15c, "fclust")
# idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
# idxpe <- PE(res.fcm4$U)
# idxpc <- PC(res.fcm4$U)
# idxmpc <- MPC(res.fcm4$U)
# 
# cat("Partition Entropy: ", idxpe)
# cat("Partition Coefficient: ", idxpc)
# cat("Modified Partition Coefficient: ", idxmpc)
# ### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
# ### the index takes the maximum value.
# cat("Fuzzy Silhouette Index: ", idxsf)
# 
# ## 10 clusters
# res.fcm4 <- ppclust2(cluster_3_10c, "fclust")
# idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
# idxpe <- PE(res.fcm4$U)
# idxpc <- PC(res.fcm4$U)
# idxmpc <- MPC(res.fcm4$U)
# 
# cat("Partition Entropy: ", idxpe)
# cat("Partition Coefficient: ", idxpc)
# cat("Modified Partition Coefficient: ", idxmpc)
# ### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
# ### the index takes the maximum value.
# cat("Fuzzy Silhouette Index: ", idxsf)
# 
# ## 5 clusters
# res.fcm4 <- ppclust2(cluster_3_5c, "fclust")
# idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
# idxpe <- PE(res.fcm4$U)
# idxpc <- PC(res.fcm4$U)
# idxmpc <- MPC(res.fcm4$U)
# 
# cat("Partition Entropy: ", idxpe)
# cat("Partition Coefficient: ", idxpc)
# cat("Modified Partition Coefficient: ", idxmpc)
# ### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
# ### the index takes the maximum value.
# cat("Fuzzy Silhouette Index: ", idxsf)
# 
# ## 4 clusters
# res.fcm4 <- ppclust2(cluster_3_4c, "fclust")
# idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
# idxpe <- PE(res.fcm4$U)
# idxpc <- PC(res.fcm4$U)
# idxmpc <- MPC(res.fcm4$U)
# 
# cat("Partition Entropy: ", idxpe)
# cat("Partition Coefficient: ", idxpc)
# cat("Modified Partition Coefficient: ", idxmpc)
# ### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
# ### the index takes the maximum value.
# cat("Fuzzy Silhouette Index: ", idxsf)
# 
# ## 3 clusters
# res.fcm4 <- ppclust2(cluster_3_3c, "fclust")
# idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
# idxpe <- PE(res.fcm4$U)
# idxpc <- PC(res.fcm4$U)
# idxmpc <- MPC(res.fcm4$U)
# 
# cat("Partition Entropy: ", idxpe)
# cat("Partition Coefficient: ", idxpc)
# cat("Modified Partition Coefficient: ", idxmpc)
# ### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
# ### the index takes the maximum value.
# cat("Fuzzy Silhouette Index: ", idxsf)
# 
# ##group by clustering ASVs and then plot using streamchart or area-----
# ### 3 clusters
# cluster_membership <- cluster_3_4c$
# cluster_3_3c$u
# cluster_membership <- cluster_3_3c$v |>
#   as.data.frame() |>
#   rownames_to_column(var = 'cluster_id') |>
#   t() 
# 
# cluster_id <- cluster_membership[1,]
# 
# cluster_membership |>
#   colnames() <- cluster_id 
# 
# cluster_membership <- cluster_membership[-1,]
# 
# cluster_membership <- cluster_membership |>
#   as.data.frame() |>
#   rownames_to_column(var = 'asv_num') |>
#   as_tibble() |>
#   dplyr::mutate(across(-asv_num, as.numeric)) |>
#   dplyr::mutate(cluster_belonging = case_when("Cluster 1" > "Cluster 2" & "Cluster 1" > "Cluster 3" ~ 'cluster_1', 
#                                               "Cluster 2" > "Cluster 1" &  "Cluster 2" > "Cluster 3" ~ 'cluster_2', 
#                                               "Cluster 3" > "Cluster 1" & "Cluster 3" > "Cluster 2" ~ 'cluster_3'))
# ## problem when they are equal then it decides for cluster 3...
# 
# 
#   # pivot_longer(cols = cluster_id, values_to = 'cluster_member_ship_degree', names_to = 'cluster_id') |>
#   # #group_by(asv_num) |>
#   # dplyr::filter(case_when(cluster_member_ship_degree) )
# ### 4 clusters
# cluster_3_4c$u
# cluster_membership <- cluster_3_4c$v |>
#   as.data.frame() |>
#   rownames_to_column(var = 'cluster_id') |>
#   t() 
# 
# cluster_id <- cluster_membership[1,]
# 
# cluster_membership |>
#   colnames() <- cluster_id 
# 
# cluster_membership <- cluster_membership[-1,]
# 
# cluster_membership <- cluster_membership |>
#   as.data.frame() |>
#   rownames_to_column(var = 'asv_num') |>
#   as_tibble() |>
#   dplyr::mutate(across(-asv_num, as.numeric)) |>
#   dplyr::mutate(cluster_belonging = case_when("Cluster 1" > "Cluster 2" & "Cluster 1" > "Cluster 3" ~ 'cluster_1', 
#                                               "Cluster 2" > "Cluster 1" &  "Cluster 2" > "Cluster 3" ~ 'cluster_2', 
#                                               "Cluster 3" > "Cluster 1" & "Cluster 3" > "Cluster 2" ~ 'cluster_3',
#                                               "Cluster 4" > "Clutser 1" & "Cluster 4" > "Clutser 2" & "Cluster 4" > "Clutser 3" ~ 'cluster_4'))



### Another way of clustering could be using the Wavelet results----
### Using the Euclidean distances then the Ward Hierarchical clustering using the results obtained from the wavelets and then cross correlate them----
## Cross-correlation functions provide a measure of association between signals. When two time series data sets are cross-correlated, a measure
## of temporal similarity is achieved. 

#upload wavelets results 

# wavelets_result_tibble_tax_3_biased <- read.csv2('data/wavelets_result_tibble_tax_3_biased.csv')
# wavelets_result_tibble_tax_02_biased <- read.csv2('data/wavelets_result_tibble_tax_02_biased.csv')

# I use the dataset that has the results with wavelets analysis without removing all samples that could be affected by the margin effects but some were
## trying to reduce the bias but not completely

wavelets_result_tibble_tax_3_biased <- read.csv('data/wavelets_analysis/wavelets_result_ed_tibble_tax_3_biased_red.csv')
wavelets_result_tibble_tax_02_biased <- read.csv('data/wavelets_analysis/wavelets_result_ed_tibble_tax_02_biased_red.csv')

asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') ##using dada2 classifier assign tax with silva 138.1

tax_bbmo_10y_new <- asv_tab_all_bloo_z_tax |>
  dplyr::select(asv_num, seq, domain, phylum, class, order, family, genus) |>
  distinct()

bloo_02 <- read.csv('data/detect_bloo/bloo_02.csv') |>
  as_tibble()

bloo_3 <- read.csv('data/detect_bloo/bloo_3.csv') |>
  as_tibble()

##reorder taxonomy as factors ----
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

## transform date into a date format for all the plots
asv_tab_all_bloo_z_tax <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

# packages
# library(cluster)    # clustering algorithms
# library(factoextra) # clustering visualization
# library(dendextend) # for comparing two dendrograms

## Prior to cross correlate the results of the Ward Hierachical clustering I explore the wavelets results at different decompositions to understand why they
## clustered together

### PA ------
time_series_1 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd1' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_2 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd2' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_3 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd3' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_4 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_5 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 's4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

# time_series_1 |>
#   pivot_longer(cols = c(-decimal_date)) |>
#   ggplot(aes(decimal_date, name))+
#   geom_tile(aes(fill = value))+
#   #scale_fill_gradientn(colours = scale_fill_viridis_b)+
#   theme_minimal()

# Dissimilarity matrix
distances_1 <- time_series_1 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_2 <- time_series_2 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_3 <- time_series_3 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_4 <- time_series_4 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_5 <- time_series_5 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

# Hierarchical clustering
hc1 <- hclust(distances_1, method = "ward.D" )  |>
  as.dendrogram()

hc2 <- hclust(distances_2, method = "ward.D" )  |>
  as.dendrogram()

hc3 <- hclust(distances_3, method = "ward.D" )  |>
  as.dendrogram()

hc4 <- hclust(distances_4, method = "ward.D" )  |>
  as.dendrogram()

hc5 <- hclust(distances_5, method = "ward.D" )  |>
  as.dendrogram()

plot(hc1)
plot(hc2)
plot(hc3)
plot(hc4)
plot(hc5)

## I plot the dendogram and the wavelets results in a heatmap and observe the clustering analysis----

year_labels <- c("2004", "2005", "2006", '2007', '2008', '2009', '2010',
                 '2011', '2012', '2013') ## create the labels with correspond with the sample num
### d1------
heatmap_data_l <- time_series_1 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc1)

dendro <- ggdendrogram(data = hc1, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(6, 0.1, 6, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_d1_PA_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()

### d2------
heatmap_data_l <- time_series_2 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc2)

dendro <- ggdendrogram(data = hc2, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(6, 0.1, 6, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_d2_PA_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()

### d3------
heatmap_data_l <- time_series_3 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc3)

dendro <- ggdendrogram(data = hc3, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(6, 0.1, 6, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_d3_PA_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()

### d4------
heatmap_data_l <- time_series_4 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc4)

dendro <- ggdendrogram(data = hc4, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(6, 0.1, 6, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_d4_PA_line.pdf", width = 8, height = 4)  # Adjust width and height as neededƒ

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()

### s4------
heatmap_data_l <- time_series_5 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc5)

dendro <- ggdendrogram(data = hc5, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(6, 0.1, 6, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_s4_PA_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()


## Visualize the clusters and try to identify some groups to label them----
# Cut tree into groups
sub_grp_1 <- cutree(hc1, k = 13)
sub_grp_2 <- cutree(hc2, k = 7)
sub_grp_3 <- cutree(hc3, k = 9)
sub_grp_4 <- cutree(hc4, k = 6)
sub_grp_5 <- cutree(hc4, k = 11)

# Number of members in each cluster
table(sub_grp_1)
table(sub_grp_2)
table(sub_grp_3)
table(sub_grp_4)
table(sub_grp_5)

df_1 <- time_series_1 |>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_2 <- time_series_2 |>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_3 <- time_series_3|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_4 <- time_series_4|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_5 <- time_series_5|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

fviz_cluster(list(data = df_1, cluster = sub_grp_1),
             geom = "point", 
             aes(color = as.factor(sub_grp_1)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_1)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_2 , cluster = sub_grp_2),
             geom = "point", 
             aes(color = as.factor(sub_grp_2)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_2)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_3 , cluster = sub_grp_3),
             geom = "point", 
             aes(color = as.factor(sub_grp_3)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_3)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_4 , cluster = sub_grp_4),
             geom = "point", 
             aes(color = as.factor(sub_grp_4)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_4)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_5 , cluster = sub_grp_5),
             geom = "point", 
             aes(color = as.factor(sub_grp_5)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_5)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

## keep the information, which cluster do each taxa belong to?-----
d1_PA_clusters <- tibble(names = names(sub_grp_1), values = sub_grp_1[names(sub_grp_1)]) |>
  rename(asv_f = names, cluster_1 = values)

d2_PA_clusters <- tibble(names = names(sub_grp_2), values = sub_grp_2[names(sub_grp_2)]) |>
  rename(asv_f = names, cluster_2 = values)

d3_PA_clusters <- tibble(names = names(sub_grp_3), values = sub_grp_3[names(sub_grp_3)]) |>
  rename(asv_f = names, cluster_3 = values)

d4_PA_clusters <- tibble(names = names(sub_grp_4), values = sub_grp_4[names(sub_grp_4)]) |>
  rename(asv_f = names, cluster_4 = values)

s5_PA_clusters <- tibble(names = names(sub_grp_5), values = sub_grp_5[names(sub_grp_5)]) |>
  rename(asv_f = names, cluster_5 = values)

clutsers_results <- d1_PA_clusters |>
  left_join(d2_PA_clusters) |>
  left_join(d3_PA_clusters) |>
  left_join(d4_PA_clusters) |>
  left_join(s5_PA_clusters)

clutsers_results_red <- clutsers_results |>
  dplyr::select(-cluster_2, -cluster_4)

cluster_1_ed <- clutsers_results_red %>%
  select(asv_f, cluster_1)

cluster_3_ed <- clutsers_results_red %>%
  select(asv_f, cluster_3)

cluster_5_ed <- clutsers_results_red %>%
  select(asv_f, cluster_5)

# Create a matrix with unique values of asv_f
asv_f_values <- sort(unique(cluster_1_ed$asv_f))

# Initialize an empty matrix with dimensions based on the number of unique asv_f values
mat_pa_1 <- matrix(0, nrow = length(asv_f_values), ncol = length(asv_f_values), dimnames = list(asv_f_values, asv_f_values))

# Fill in the matrix based on whether the asv_f have the same value for cluster_1
for(i in 1:nrow(mat_pa_1)) {
  for(j in 1:ncol(mat_pa_1)) {
    mat_pa_1[i, j] <- ifelse(cluster_1_ed[cluster_1_ed$asv_f == rownames(mat_pa_1)[i], "cluster_1"] == cluster_1_ed[cluster_1_ed$asv_f == colnames(mat_pa_1)[j], "cluster_1"], 1, 0)
  }
}

dim(mat_pa_1)

# Create a matrix with unique values of asv_f
asv_f_values <- sort(unique(cluster_3_ed$asv_f))

# Initialize an empty matrix with dimensions based on the number of unique asv_f values
mat_pa_3 <- matrix(0, nrow = length(asv_f_values), ncol = length(asv_f_values), dimnames = list(asv_f_values, asv_f_values))

# Fill in the matrix based on whether the asv_f have the same value for cluster_3
for(i in 1:nrow(mat_pa_3)) {
  for(j in 1:ncol(mat_pa_3)) {
    mat_pa_3[i, j] <- ifelse(cluster_3_ed[cluster_3_ed$asv_f == rownames(mat_pa_3)[i], "cluster_3"] == cluster_3_ed[cluster_3_ed$asv_f == colnames(mat_pa_3)[j], "cluster_3"], 1, 0)
  }
}

dim(mat_pa_3)

# Create a matrix with unique values of asv_f
asv_f_values <- sort(unique(cluster_5_ed$asv_f))

# Initialize an empty matrix with dimensions based on the number of unique asv_f values
mat_pa_5 <- matrix(0, nrow = length(asv_f_values), ncol = length(asv_f_values), dimnames = list(asv_f_values, asv_f_values))

# Fill in the matrix based on whether the asv_f have the same value for cluster_5
for(i in 1:nrow(mat_pa_5)) {
  for(j in 1:ncol(mat_pa_5)) {
    mat_pa_5[i, j] <- ifelse(cluster_5_ed[cluster_5_ed$asv_f == rownames(mat_pa_5)[i], "cluster_5"] == cluster_5_ed[cluster_5_ed$asv_f == colnames(mat_pa_5)[j], "cluster_5"], 1, 0)
  }
}

dim(mat_pa_5)

sum_matrix <- mat_pa_1 + mat_pa_3 + mat_pa_5

sum_matrix |>
  as.data.frame() |>
  rownames_to_column(var = 'asv_f') |>
  as_tibble() |>
  pivot_longer(cols = -asv_f, names_to = 'asv_f_2', values_to = 'consistency') |>
  dplyr::filter(asv_f == asv_f_2) |>
  dplyr::filter(consistency != 3) ## check that the same asv compared to itself is always 3 to be sure that the code is working correcly

results_clusters_consistency_pa <- sum_matrix |>
  as.data.frame() |>
  rownames_to_column(var = 'asv_f') |>
  as_tibble() |>
  pivot_longer(cols = -asv_f, names_to = 'asv_f_2', values_to = 'consistency') |>
  dplyr::filter(asv_f != asv_f_2 &
                consistency != 0 & 
                  consistency != 1 &
                  consistency != 2) 

results_clusters_consistency_pa |>
  dplyr::summarise(n_unique_chars_asv_f = n_distinct(asv_f),
            n_unique_chars_asv_f_2 = n_distinct(asv_f)) ## from 47 bloomers in PA I have 27 that always group together at distance 10 (revise the number) and i got 11 different clusters

### plot the grouped ASVs in the initial dataset find what brings them together------
asv_tab_all_bloo_z_tax |>
  colnames()

## harbor responsive group1
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv302', 'asv511', 'asv311')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.25))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## harbor responsive group 2?
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv194', 'asv559')) |> # 'asv105', at the threshold distance 10 now is not included
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

# ## harbor responsive group 3? at the distance 10 they are not clustered together anymore.
# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(asv_num %in% c('asv22', 'asv85', 'asv163', 'asv219', 'asv80', 'asv192', 'asv49')) |>
#   dplyr::filter(abundance_type == 'relative_abundance' &
#                   fraction == '3') |>
#   group_by(date, fraction) |>
#   dplyr::mutate(max_abund = sum(abundance_value)) |>
#   ungroup() |>
#   ggplot(aes(date, max_abund))+
#   scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
#                    #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
#                    #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
#   )+
#   # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
#   # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
#   scale_y_continuous(labels = percent_format(), expand = c(0,0))+
#   geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   #facet_wrap(vars(asv_num_f))+
#   labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Family')+
#   guides(fill = guide_legend(ncol = 6, size = 10,
#                              override.aes = aes(label = '')),
#          alpha = 'none')+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(), strip.text = element_text(size = 7),
#         legend.position = 'bottom', axis.text.y = element_text(size = 8),
#         axis.title = element_text(size = 8), strip.background = element_blank(), 
#         legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
#         plot.margin = margin(2,5,0,5)) 

## harbor responsive group 4?
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv276', 'asv223')) |> # 'asv264',, 'asv471', 'asv752'
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## recurrent
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv7', 'asv15')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## ephemeral
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv317', 'asv200', 'asv113')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## recurrent
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv116', 'asv182', 'asv84')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## recurrent
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv100', 'asv72', 'asv42')) |> #'asv25', 
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## recurrent
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c( 'asv1',  'asv4')) |> #'asv31', 'asv23',
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## ephemeral
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv69', 'asv225')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## harbor affected but more persistent blooms
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv163', 'asv219')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous( expand = c(0,0))+ #labels = percent_format(),
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## harbor affected but more persistent blooms
asv_tab_all_bloo_z_tax |>
  #dplyr::filter(asv_num %in% c('asv179')) |>
  dplyr::filter(asv_num %in% c('asv264', 'asv471', 'asv752')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## study all ASVs that do not have always the same group ----
asvs_3_with_cluster <- results_clusters_consistency_pa |>
  separate(asv_f_2, into = c('family', 'asv_num'), sep = '\\.', remove = F) |>
  distinct(asv_num)

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(!asv_num %in% asvs_3_with_cluster$asv_num) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous( expand = c(0,0))+ #, limits = c(0,0.25) #labels = percent_format(),
  geom_area(aes(date, abundance_value, fill = order_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(asv_num_f), scales = 'free')+
  labs(x = 'Time', y = 'rCLR', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

### CLUSTERS IN THE PARTICLE ATTACHED FRACTION (EUCLIDEAN DISTANCE 10) -----
#### I have 11 different clusters
# c('asv302', 'asv511', 'asv311')
# c('asv194', 'asv559')
# c('asv276', 'asv223')
# c('asv7', 'asv15')
# c('asv317', 'asv200', 'asv113')
# c('asv116', 'asv182', 'asv84')
# c('asv100', 'asv72', 'asv42')
# c( 'asv1',  'asv4')
# c('asv69', 'asv225')
# c('asv163', 'asv219')
# c('asv264', 'asv471', '182')

###table with the label corresponding to their signals-------
### for each transformation used each bloomer has a different label only if it gives a clear strong signal 
#### this table will contain the asv num + for each transformation which label did it got

###The most important signals for our dataset have been the d1, d3, and s5 they are the ones that we will use for 
### labeling the bloomers.

## With this analysis we will pick a representative of each group to perform the random tree analysis-----
## asvs with less clear clutser asv179, asv11, asv385, asv27, asv17, asv43, asv118, asv126, asv28

## when we cut at distance 10 then we have more ASVs with a less clear cluster: 
## asv179 asv385 asv27  asv153 asv17  asv77  asv43  asv192 asv118 asv23  asv85  asv25  asv80  asv126 asv105 asv28  asv31  asv22  asv11  asv49
# asv_tab_all_bloo_z_tax |>
#   dplyr::filter(asv_num %in% bloo_3$value) |>
#   dplyr::filter(!asv_num %in% asvs_3_with_cluster$asv_num) |>
#   distinct(asv_num_f) |>
#   as.vector()

wavelets_result_ed_tibble_tax_3_biased_red_coeff |>
  dplyr::filter(asv_num %in% c('asv179', 'asv11', 'asv385'))

wavelets_result_ed_tibble_tax_3_biased_red_coeff |>
  dplyr::filter(asv_num %in% c('asv27', 'asv17', 'asv43'))

wavelets_result_ed_tibble_tax_3_biased_red_coeff |>
  dplyr::filter(asv_num %in% c('asv118', 'asv126', 'asv28'))

bloo_3$value

bloo_3_types_summary <- bloo_3 |>
  rename('asv_num' = 'value') |>
  dplyr::mutate(
    recurrency = case_when(asv_num %in% c('asv311', 'asv302', 'asv511') ~ 'no',
                           asv_num %in% c('asv194', 'asv105', 'asv559') ~ 'no',
                           asv_num %in% c('asv22', 'asv85', 'asv163', 'asv219', 'asv80', 'asv192', 'asv49')  ~ 'no',
                           asv_num %in% c('asv276', 'asv264', 'asv223', 'asv471', 'asv752')  ~ 'no',
                           asv_num %in% c('asv7', 'asv15') ~ 'yes',
                           asv_num %in% c('asv317', 'asv200', 'asv113') ~ 'no',
                           asv_num %in% c('asv116', 'asv182', 'asv84') ~ 'no',
                           asv_num %in% c('asv100', 'asv25', 'asv72', 'asv42') ~ 'yes',
                           asv_num %in% c('asv23', 'asv1', 'asv31', 'asv4')  ~ 'yes',
                           asv_num %in% c('asv69', 'asv225') ~ 'no',
                           asv_num %in% c('asv153', 'asv77') ~ 'no',
                           asv_num == 'asv179'~ 'no',
                           asv_num == 'asv11' ~ 'no',
                           asv_num == 'asv385' ~ 'no',
                           asv_num == 'asv27' ~ 'yes',
                           asv_num == 'asv17' ~ 'no',
                           asv_num == 'asv43' ~ 'no',
                           asv_num == 'asv118'~ 'no',
                           asv_num == 'asv126'~ 'no',
                           asv_num == 'asv28'~ 'yes'
                           
    ),
    frequency = case_when(asv_num %in% c('asv311', 'asv302', 'asv511') ~ 'stochastic',
                          asv_num %in% c('asv194', 'asv105', 'asv559') ~ 'stochastic',
                          asv_num %in% c('asv276', 'asv264', 'asv223', 'asv471', 'asv752') ~ 'stochastic',
                          asv_num %in% c('asv22', 'asv85', 'asv163', 'asv219', 'asv80', 'asv192', 'asv49')  ~ 'stochastic',
                          asv_num %in% c('asv7', 'asv15') ~ 'seasonal',
                          asv_num %in% c('asv317', 'asv200', 'asv113') ~ 'stochastic',
                          asv_num %in% c('asv116', 'asv182', 'asv84') ~ 'stochastic',
                          asv_num %in% c('asv100', 'asv25', 'asv72', 'asv42') ~ 'seasonal',
                          asv_num %in% c('asv23', 'asv1', 'asv31', 'asv4')  ~ 'seasonal',
                          asv_num %in% c('asv69', 'asv225') ~ 'stochastic',
                          asv_num %in% c('asv153', 'asv77') ~ 'stochastic',
                          asv_num == 'asv179' ~ 'stochastic',
                          asv_num == 'asv11'~ 'stochastic',
                          asv_num == 'asv385'~ 'stochastic',
                          asv_num == 'asv27'~ 'seasonal',
                          asv_num == 'asv17'~ 'stochastic',
                          asv_num == 'asv43'~ 'stochastic',
                          asv_num == 'asv118'~ 'stochastic',
                          asv_num == 'asv126'~ 'stochastic',
                          asv_num == 'asv28'~ 'seasonal'
    ))

bloo_3_types_summary <- bloo_3_types_summary |>
  dplyr::left_join(occurence_perc_3)|>
  dplyr::mutate(fraction = '3')

bloo_3_clustering_results <- bloo_3 |>
  rename('asv_num' = 'value') |>
  dplyr::mutate(
    # I update the clustering groups for those at euclidean distance 10.
    clustering_group = case_when(asv_num %in% c('asv302', 'asv511', 'asv311') ~ 'cl_1',
                                 #asv_num %in% c('asv311', 'asv302', 'asv511') ~ 'cl_1',
                                 # asv_num %in% c('asv194', 'asv105', 'asv559') ~ 'cl_2',
                                 # asv_num %in% c('asv22', 'asv85', 'asv163', 'asv219', 'asv80', 'asv192', 'asv49')  ~ 'cl_3',
                                 # asv_num %in% c('asv276', 'asv264', 'asv223', 'asv471', 'asv752') ~ 'cl_4',
                                 # asv_num %in% c('asv7', 'asv15')  ~ 'cl_5',
                                 # asv_num %in% c('asv317', 'asv200', 'asv113') ~ 'cl_6',
                                 # asv_num %in% c('asv116', 'asv182', 'asv84') ~ 'cl_7',
                                 # asv_num %in% c('asv100', 'asv25', 'asv72', 'asv42') ~ 'cl_8',
                                 # asv_num %in% c('asv23', 'asv1', 'asv31', 'asv4') ~ 'cl_9',
                                 # asv_num %in% c('asv69', 'asv225') ~ 'cl_10',
                                 # asv_num %in% c('asv153', 'asv77') ~ 'cl_11',
                                 asv_num %in% c('asv194', 'asv559') ~ 'cl_2',
                                 asv_num %in% c('asv276', 'asv223') ~ 'cl_3',
                                 asv_num %in% c('asv7', 'asv15') ~ 'cl_4',
                                 asv_num %in% c('asv317', 'asv200', 'asv113') ~ 'cl_5',
                                 asv_num %in% c('asv116', 'asv182', 'asv84') ~ 'cl_6',
                                 asv_num %in% c('asv100', 'asv72', 'asv42') ~ 'cl_7',
                                 asv_num %in% c( 'asv1',  'asv4')  ~ 'cl_8',
                                 asv_num %in% c('asv69', 'asv225') ~ 'cl_9',
                                 asv_num %in% c('asv163', 'asv219') ~ 'cl_10',
                                 asv_num %in% c('asv264', 'asv471', 'asv752') ~ 'cl_11',
                                 asv_num %in% c('asv179', 'asv385', 'asv27',  'asv153', 'asv17', 'asv77',  'asv43',  'asv192', 'asv118', 
                                                'asv23' , 'asv85',  'asv25',  'asv80',  'asv126', 'asv105', 'asv28', 'asv31' , 'asv22',  'asv11',  'asv49') ~ 'unclear'))|>
  dplyr::mutate(fraction = '3')

## clusters meaning ----
labs_clusters_pa <- as_labeller(c(cl_1 = 'habor restoration 1', 
                                  cl_2 = 'habor restoration 2',
                                  cl_3 = 'habor restoration 3',
                                  cl_4 = 'harbor restoration 4',
                                  cl_5 = 'seasonal 1',
                                  cl_6 = 'harbor restoration 5', 
                                  cl_7 = 'recurrent random',
                                 cl_8 = 'seasonal 2',
                               cl_9 = 'seasonal 3',
                                cl_10 = 'harbor restoration 6',
                               cl_11 = 'harbor restoration 7',
                               unclear = 'ungruped'))

## add the occurrence and the relative abundance of this taxa----
occurrence_bloo_bbmo |>
  colnames()

### PLOT THE RESULTS OF THE CLUSTERING ANALYSIS ---------
bloo_3_types_summary |>
  colnames()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::left_join(bloo_3_clustering_results) |>
  dplyr::filter(clustering_group != 'unclear') |>
  dplyr::group_by(clustering_group) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(clustering_group), scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::left_join(bloo_3_clustering_results) |>
  dplyr::filter(clustering_group != 'unclear') |>
  #dplyr::group_by(clustering_group, asv_num) |>
  #dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, abundance_value))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 1,  position = 'stack', outline.type = 'upper')+
  #geom_line(aes(date, abundance_value, group = asv_num), color = 'black', linewidth = 0.1) +
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(clustering_group), scales = 'free')+
  labs(x = 'Time', y = 'rCLR', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  dplyr::left_join(bloo_3_types_summary) |>
  group_by(date, fraction, recurrency, frequency) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(recurrency~frequency, scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  dplyr::left_join(bloo_3_types_summary) |>
  dplyr::left_join(bloo_3_clustering_results) |>
  group_by(date, fraction, frequency) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(frequency), scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))

## explore the examples by the different types ----
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(!asv_num %in% asvs_3_with_cluster$asv_num) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::left_join(bloo_3_clustering_results) |>
  dplyr::left_join(bloo_3_types_summary) |>
  dplyr::filter(frequency == 'stochastic') |>
  group_by(date, fraction,  occurrence_category) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(vars(occurrence_category), scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::left_join(bloo_3_types_summary) |>
  dplyr::left_join(bloo_3_clustering_results) |>
  dplyr::filter(recurrency == 'yes') |>
  group_by(date, fraction, recurrency) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #facet_wrap(vars(asv_num_f), scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::left_join(bloo_3_types_summary) |>
  dplyr::left_join(bloo_3_clustering_results) |>
  dplyr::filter(recurrency == 'yes') |>
  group_by(date, fraction, recurrency) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_grid(vars(clustering_group))+
  geom_smooth(method = 'loess', se = F, span = 0.3)+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))

#### SEASONAL BLOOMERS -----
seasonal_clusters_labs <- as_labeller(c( cl_5 = 'first bloom type 1',
                                                             cl_8 = 'second bloom',
                                                             cl_9 = 'third bloom',
                                                            asv27 = 'first bloom type 2',
                                                             asv28 = 'between second and third bloom'))

bloo_3_types_summary |>
  left_join(bloo_3_clustering_results) |>
  dplyr::filter(recurrency == 'yes',
                clustering_group == 'unclear') 

subset <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::left_join(bloo_3_types_summary) |>
  left_join(bloo_3_clustering_results) |>
  dplyr::filter(recurrency == 'yes') |>
  dplyr::mutate(clustering_group = case_when(clustering_group == 'unclear' &
                                               asv_num == 'asv27' ~ 'asv27',
                                             clustering_group == 'unclear' &
                                               asv_num == 'asv28' ~ 'asv28',
                                             T ~ clustering_group)) |>
  group_by(day_of_year, clustering_group) |> 
  dplyr::mutate(max_abund = sum(abundance_value))

subset |>
  dplyr::filter(!clustering_group %in% c('asv28', 'asv27')) |>
  ggplot(aes(day_of_year, max_abund, group = clustering_group, fill = clustering_group))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  scale_x_continuous(expand = c(0,0), limits = c(0, 365))+
  facet_grid(vars(asv_num))+
  #geom_area(aes(day_of_year, y = abundance_value, group = asv_num_f, fill = clustering_group), alpha = 0.3,  position= 'stack')+
  scale_fill_manual(values = palete_seasonal_bloo)+
  #geom_area(aes(day_of_year,max_abund), alpha = 0.5,  position='identity')+
  scale_color_manual(values = palete_seasonal_bloo)+
  geom_point(aes(color = clustering_group, y = abundance_value), alpha = 0.2)+
  geom_smooth(aes(color = clustering_group, y = abundance_value), span = 0.8)+
  #facet_grid(year~clustering_group)+
  #geom_smooth(aes(group = clustering_group), method = 'loess', se = F, span = 0.6)+
  labs(x = 'Day of the year', y = 'Relative abundance (%)', fill = 'Clustering group')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         color = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))

### Stochastic blooms and persitent-----
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::mutate(relative_abundance_category = case_when(mean(abundance_value) > 0.02 ~ 'abundant',
                                                        mean(abundance_value) < 0.02 ~ 'mid',
                                                        mean(abundance_value) < 0.01 ~ 'rare')) |>
  dplyr::left_join(bloo_3_types_summary) |>
  dplyr::filter(recurrency == 'yes') |>
  dplyr::mutate(clustering_group = case_when(clustering_group == 'unclear' &
                                               asv_num == 'asv27' ~ 'asv27',
                                             clustering_group == 'unclear' &
                                               asv_num == 'asv28' ~ 'asv28',
                                             T ~ clustering_group)) |>
  group_by(day_of_year, clustering_group) |> 
  dplyr::mutate(max_abund = sum(abundance_value))

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::mutate(relative_abundance_category = case_when(mean(abundance_value) > 0.02 ~ 'abundant',
                                                        mean(abundance_value) < 0.02 ~ 'mid',
                                                        mean(abundance_value) < 0.01 ~ 'rare')) |>
  dplyr::left_join(bloo_3_types_summary) |>
  dplyr::filter(recurrency == 'no' &
                  frequency == 'stochastic' &
                  type_of_bloom == 'ephemeral') |>
  dplyr::filter(!clustering_group %in% c('asv28', 'asv27')) |>
  ggplot(aes(day_of_year, max_abund, group = clustering_group, fill = clustering_group))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  scale_x_continuous(expand = c(0,0), limits = c(0, 365))+
  #facet_grid(vars(asv_num))+
  geom_area(aes(day_of_year, y = abundance_value, group = asv_num_f, fill = clustering_group), alpha = 0.3,  position= 'identity')+
  #scale_fill_manual(values = palete_seasonal_bloo)+
  #geom_area(aes(day_of_year,max_abund), alpha = 0.5,  position='identity')+
  #scale_color_manual(values = palete_seasonal_bloo)+
  #geom_point(aes(color = clustering_group, y = abundance_value), alpha = 0.2)+
  #geom_smooth(aes(color = clustering_group, y = abundance_value), span = 0.8)+
  facet_grid(vars(year))+
  #geom_smooth(aes(group = clustering_group), method = 'loess', se = F, span = 0.6)+
  labs(x = 'Day of the year', y = 'Relative abundance (%)', fill = 'Clustering group')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         color = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  group_by(asv_num) |>
  dplyr::mutate(relative_abundance_category = case_when(mean(abundance_value) > 0.02 ~ 'abundant',
                                                        mean(abundance_value) < 0.02 ~ 'mid',
                                                        mean(abundance_value) < 0.01 ~ 'rare')) |>
  dplyr::left_join(bloo_3_types_summary) |>
  dplyr::filter(recurrency == 'no' &
                  frequency == 'stochastic' &
                  type_of_bloom == 'persistent') |>
  dplyr::filter(!clustering_group %in% c('asv28', 'asv27')) |>
  group_by(clustering_group, day_of_year, year) |>
  dplyr::summarize(abundance_cluster = sum(abundance_value)) |>
  ggplot(aes(day_of_year, max_abund, group = clustering_group, fill = clustering_group))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  scale_x_continuous(expand = c(0,0), limits = c(0, 365))+
  #facet_grid(vars(asv_num))+
  geom_area(aes(day_of_year, y = abundance_cluster, group = clustering_group, fill = clustering_group), alpha = 0.3,  position= 'identity')+
  #scale_fill_manual(values = palete_seasonal_bloo)+
  #geom_area(aes(day_of_year,max_abund), alpha = 0.5,  position='identity')+
  #scale_color_manual(values = palete_seasonal_bloo)+
  #geom_point(aes(color = clustering_group, y = abundance_value), alpha = 0.2)+
  #geom_smooth(aes(color = clustering_group, y = abundance_value), span = 0.8)+
  facet_grid(vars(year))+
  #geom_smooth(aes(group = clustering_group), method = 'loess', se = F, span = 0.6)+
  labs(x = 'Day of the year', y = 'Relative abundance (%)', fill = 'Clustering group')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         color = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))

### The same for the FL fraction----
time_series_1 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |> # I remove the SAR11 cluster from the analyisis as they as discarded as true bloomers
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd1' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_2 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |> # I remove the SAR11 cluster from the analyisis as they as discarded as true bloomers
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd2' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_3 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |> # I remove the SAR11 cluster from the analyisis as they as discarded as true bloomers
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd3' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_4 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |> # I remove the SAR11 cluster from the analyisis as they as discarded as true bloomers
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_5 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |> # I remove the SAR11 cluster from the analyisis as they as discarded as true bloomers
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = case_when(!is.na(family) ~ paste0(family,'.',asv_num),
                                  is.na(family) & !is.na(order) ~ paste0(order, '.', asv_num),
                                  is.na(family) & is.na(order) ~ paste0(class, '.', asv_num),
                                  is.na(family) & is.na(order) & is.na(class) ~ paste0(phylum, '.', asv_num))) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 's4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

# time_series_1 |>
#   pivot_longer(cols = c(-decimal_date)) |>
#   ggplot(aes(decimal_date, name))+
#   geom_tile(aes(fill = value))+
#   #scale_fill_gradientn(colours = scale_fill_viridis_b)+
#   theme_minimal()

# Dissimilarity matrix
distances_1 <- time_series_1 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_2 <- time_series_2 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_3 <- time_series_3 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_4 <- time_series_4 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_5 <- time_series_5 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

# Hierarchical clustering
### when analysing clustering results and not ploting do not run as.dendogram
hc1 <- hclust(distances_1, method = "ward.D" )  |>
  as.dendrogram()

hc2 <- hclust(distances_2, method = "ward.D" )  |>
  as.dendrogram()

hc3 <- hclust(distances_3, method = "ward.D" )  |>
  as.dendrogram()

hc4 <- hclust(distances_4, method = "ward.D" )  |>
  as.dendrogram()

hc5 <- hclust(distances_5, method = "ward.D" )  |>
  as.dendrogram()

plot(hc1)
plot(hc2)
plot(hc3)
plot(hc4)
plot(hc5)

## I plot the dendogram and the wavelets results in a heatmap and observe the clustering analysis----
### d1------
heatmap_data_l <- time_series_1 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc1)

dendro <- ggdendrogram(data = hc1, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(7, 0.1, 7, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_d1_FL_no_sar11_cluster_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()

### d2------
heatmap_data_l <- time_series_2 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc2)

dendro <- ggdendrogram(data = hc2, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(7, 0.1, 7, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_d2_FL_no_sar11_cluster_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()

### d3------
heatmap_data_l <- time_series_3 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc3)

dendro <- ggdendrogram(data = hc3, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(7, 0.1, 7, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_d3_FL_no_sar11_cluster_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()

### d4------
heatmap_data_l <- time_series_4 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

asv_order <- order.dendrogram(hc4)

dendro <- ggdendrogram(data = hc4, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(7, 0.1, 7, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-5.21, 5.21), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_d4_FL_no_sar11_cluster_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()

### s4------
heatmap_data_l <- time_series_5 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'asv_num')

heatmap_data_l |>
  dplyr::filter(!is.na(value)) %$%
  value |>
  range()

asv_order <- order.dendrogram(hc5)

dendro <- ggdendrogram(data = hc5, rotate = T)+
  scale_y_reverse(expand = c(0,0))+
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  labs(y = 'Distance')+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5),
        plot.margin = unit(c(7, 0.1, 7, 5), "mm"))

heatmap_data_l$asv_num <- factor(x = heatmap_data_l$asv_num,
                                 levels = heatmap_data_l$asv_num[asv_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = asv_num)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-8.87, 8.87), na.value = '#D7D6D3')+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'),
        plot.margin = unit(c(5, 1, 5, 0), "mm"))

# Create a new page
pdf("results/figures/hierarchical_clustering_scales/clust_s4_FL_no_sar11_cluster_line.pdf", width = 8, height = 4)  # Adjust width and height as needed

# Arrange and print the plots
composition_of_plots <- grid.arrange(dendro, 
                                     heatmap_plot,
                                     ncol = 2,
                                     widths = c(1, 3))

# Close the PDF device
dev.off()



## Visualize the clusters and try to identify some groups to label them----
# Cut tree into 4 groups
sub_grp_1 <- cutree(hc1, k = 6)
sub_grp_2 <- cutree(hc2, k = 6)
sub_grp_3 <- cutree(hc3, k = 5)
sub_grp_4 <- cutree(hc4, k = 6)
sub_grp_5 <- cutree(hc4, k = 4)

# Number of members in each cluster
table(sub_grp_1)
table(sub_grp_2)
table(sub_grp_3)
table(sub_grp_4)
table(sub_grp_5)

df_1 <- time_series_1 |>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_2 <- time_series_2 |>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_3 <- time_series_3|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_4 <- time_series_4|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_5 <- time_series_5|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

fviz_cluster(list(data = df_1, cluster = sub_grp_1),
             geom = "point", 
             aes(color = as.factor(sub_grp_1)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_1)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_2 , cluster = sub_grp_2),
             geom = "point", 
             aes(color = as.factor(sub_grp_2)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_2)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_3 , cluster = sub_grp_3),
             geom = "point", 
             aes(color = as.factor(sub_grp_3)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_3)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_4 , cluster = sub_grp_4),
             geom = "point", 
             aes(color = as.factor(sub_grp_4)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_4)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_5 , cluster = sub_grp_5),
             geom = "point", 
             aes(color = as.factor(sub_grp_5)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_5)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

## Visualize the clusters and try to identify some groups to label them----

# Cut tree into 4 groups
sub_grp_1 <- cutree(hc1, k = 6)
sub_grp_2 <- cutree(hc2, k = 5)
sub_grp_3 <- cutree(hc3, k = 5)
sub_grp_4 <- cutree(hc4, k = 9)
sub_grp_5 <- cutree(hc4, k = 4)

# Number of members in each cluster
table(sub_grp_1)
table(sub_grp_2)
table(sub_grp_3)
table(sub_grp_4)
table(sub_grp_5)

df_1 <- time_series_1 |>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_2 <- time_series_2 |>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_3 <- time_series_3|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_4 <- time_series_4|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

df_5 <- time_series_5|>
  dplyr::select(-decimal_date) |>
  na.omit() %>%
  t()

fviz_cluster(list(data = df_1, cluster = sub_grp_1),
             geom = "point", 
             aes(color = as.factor(sub_grp_1)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_1)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_2 , cluster = sub_grp_2),
             geom = "point", 
             aes(color = as.factor(sub_grp_1)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_1)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_3 , cluster = sub_grp_3),
             geom = "point", 
             aes(color = as.factor(sub_grp_1)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_1)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_4 , cluster = sub_grp_4),
             geom = "point", 
             aes(color = as.factor(sub_grp_1)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_1)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_5 , cluster = sub_grp_5),
             geom = "point", 
             aes(color = as.factor(sub_grp_1)),  # Specify color aesthetic
             ggtheme = theme_light())+
  scale_color_manual(values = palette_clustering)+
  geom_text(aes(label = rownames(df_1)), check_overlap = TRUE, size = 2, vjust = 1.5) +
  scale_fill_manual(values = palette_clustering)+
  labs(title = '')+
  theme(text = element_text(size = 6))

## keep the information, which cluster do each taxa belong to?
d1_fl_clusters <- tibble(names = names(sub_grp_1), values = sub_grp_1[names(sub_grp_1)]) |>
  rename(asv_f = names, cluster_1 = values)

d2_fl_clusters <- tibble(names = names(sub_grp_2), values = sub_grp_2[names(sub_grp_2)]) |>
  rename(asv_f = names, cluster_2 = values)

d3_fl_clusters <- tibble(names = names(sub_grp_3), values = sub_grp_3[names(sub_grp_3)]) |>
  rename(asv_f = names, cluster_3 = values)

d4_fl_clusters <- tibble(names = names(sub_grp_4), values = sub_grp_4[names(sub_grp_4)]) |>
  rename(asv_f = names, cluster_4 = values)

s5_fl_clusters <- tibble(names = names(sub_grp_5), values = sub_grp_5[names(sub_grp_5)]) |>
  rename(asv_f = names, cluster_5 = values)

clutsers_results <- d1_fl_clusters |>
  left_join(d2_fl_clusters) |>
  left_join(d3_fl_clusters) |>
  left_join(d4_fl_clusters) |>
  left_join(s5_fl_clusters)

clutsers_results_red <- clutsers_results |>
  dplyr::select(-cluster_2, -cluster_4)

cluster_1_ed <- clutsers_results_red %>%
  select(asv_f, cluster_1)

cluster_3_ed <- clutsers_results_red %>%
  select(asv_f, cluster_3)

cluster_5_ed <- clutsers_results_red %>%
  select(asv_f, cluster_5)

# Create a matrix with unique values of asv_f-----
asv_f_values <- sort(unique(cluster_1_ed$asv_f))

# Initialize an empty matrix with dimensions based on the number of unique asv_f values
mat_fl_1 <- matrix(0, nrow = length(asv_f_values), ncol = length(asv_f_values), dimnames = list(asv_f_values, asv_f_values))

# Fill in the matrix based on whether the asv_f have the same value for cluster_1
for(i in 1:nrow(mat_fl_1)) {
  for(j in 1:ncol(mat_fl_1)) {
    mat_fl_1[i, j] <- ifelse(cluster_1_ed[cluster_1_ed$asv_f == rownames(mat_fl_1)[i], "cluster_1"] == cluster_1_ed[cluster_1_ed$asv_f == colnames(mat_fl_1)[j], "cluster_1"], 1, 0)
  }
}

dim(mat_fl_1)

# Create a matrix with unique values of asv_f
asv_f_values <- sort(unique(cluster_3_ed$asv_f))

# Initialize an empty matrix with dimensions based on the number of unique asv_f values
mat_fl_3 <- matrix(0, nrow = length(asv_f_values), ncol = length(asv_f_values), dimnames = list(asv_f_values, asv_f_values))

# Fill in the matrix based on whether the asv_f have the same value for cluster_3
for(i in 1:nrow(mat_fl_3)) {
  for(j in 1:ncol(mat_fl_3)) {
    mat_fl_3[i, j] <- ifelse(cluster_3_ed[cluster_3_ed$asv_f == rownames(mat_fl_3)[i], "cluster_3"] == cluster_3_ed[cluster_3_ed$asv_f == colnames(mat_fl_3)[j], "cluster_3"], 1, 0)
  }
}

dim(mat_fl_3)

# Create a matrix with unique values of asv_f
asv_f_values <- sort(unique(cluster_5_ed$asv_f))

# Initialize an empty matrix with dimensions based on the number of unique asv_f values
mat_fl_5 <- matrix(0, nrow = length(asv_f_values), ncol = length(asv_f_values), dimnames = list(asv_f_values, asv_f_values))

# Fill in the matrix based on whether the asv_f have the same value for cluster_5
for(i in 1:nrow(mat_fl_5)) {
  for(j in 1:ncol(mat_fl_5)) {
    mat_fl_5[i, j] <- ifelse(cluster_5_ed[cluster_5_ed$asv_f == rownames(mat_fl_5)[i], "cluster_5"] == cluster_5_ed[cluster_5_ed$asv_f == colnames(mat_fl_5)[j], "cluster_5"], 1, 0)
  }
}

dim(mat_fl_5)

sum_matrix <- mat_fl_1 + mat_fl_3 + mat_fl_5

sum_matrix |>
  as.data.frame() |>
  rownames_to_column(var = 'asv_f') |>
  as_tibble() |>
  pivot_longer(cols = -asv_f, names_to = 'asv_f_2', values_to = 'consistency') |>
  dplyr::filter(asv_f == asv_f_2) |>
  dplyr::filter(consistency != 3) ## check that the same asv compared to itself is always 3 to be sure that the code is working correcly

results_clusters_consistency_fl <- sum_matrix |>
  as.data.frame() |>
  rownames_to_column(var = 'asv_f') |>
  as_tibble() |>
  pivot_longer(cols = -asv_f, names_to = 'asv_f_2', values_to = 'consistency') |>
  dplyr::filter(asv_f != asv_f_2 &
                  consistency != 0 & 
                  consistency != 1 &
                  consistency != 2) 

results_clusters_consistency_fl |>
  summarise(n_unique_chars_asv_f = n_distinct(asv_f),
            n_unique_chars_asv_f_2 = n_distinct(asv_f)) ## from 20 bloomers I have 10 that always group together and form 3 clusters

### plot the grouped ASVs in the initial dataset find what brings them together-----
asv_tab_all_bloo_z_tax |>
  colnames()

## recurrent random bloom----
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv58', 'asv178')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '0.2') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
                                                   # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous( expand = c(0,0))+ #labels = percent_format(),
  geom_area(aes(date, abundance_value, fill = order_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  #facet_wrap(vars(asv_num_f))+
  labs(x = 'Time', y = 'rCLR', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## ephemeral BLOOMS -----
results_clusters_consistency_fl |>
  dplyr::filter(asv_f == 'AEGEAN-169 marine group.asv555')

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv555', 'asv114', 'asv249', 'asv237', 'asv563', 'asv282')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '0.2') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #geom_point(aes(y = abundance_value))+
  labs(x = 'Time', y = 'rCLR', fill = 'Order')+
  #facet_wrap(vars(asv_num_f))+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## Recurrent bloom -----
results_clusters_consistency_fl |>
  dplyr::filter(asv_f == 'Clade I.asv15')

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv15', 'asv7')) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '0.2') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous( expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = order_f, group = asv_num), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  labs(x = 'Time', y = 'rCLR', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## SAR11 GROUP (CLADE-I and II) (DISCARDED NO REAL BLOOMERS) -----
results_clusters_consistency_fl |>
  dplyr::filter(asv_f == 'Clade I.asv2')

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  facet_wrap(vars(asv_num_f))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.25))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## study all ASVs that do not have always the same group ----
asvs_02_with_cluster <- results_clusters_consistency_fl |>
  separate(asv_f_2, into = c('family', 'asv_num'), sep = '\\.', remove = F) 

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(!asv_num %in% asvs_02_with_cluster$asv_num) |>
  distinct(asv_num) |>
  as_vector()

##  "asv38"    "asv8"    "asv5"    "asv3"    "asv2"   "asv27"   "asv17"   "asv62"    "asv1"   "asv11"

bloo_02 <- bloo_02 |>
  as_tibble()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(!asv_num %in% asvs_02_with_cluster$asv_num) |>
  dplyr::filter(abundance_type == 'rclr' &
                  fraction == '0.2') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = order_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  facet_wrap(vars(asv_num_f), scales = 'free')+
  labs(x = 'Time', y = 'rCLR', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## winter blooms?----
## they are not clearly clustered
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(asv_num %in% c('asv17', 'asv62')) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = order_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_order_assigned_bloo)+
  #facet_wrap(vars(asv_num_f), scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## winter vs summer blooms
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% c('asv15', 'asv7', 'asv17', 'asv62')) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  dplyr::mutate(seasonality = case_when(asv_num %in% c('asv15', 'asv7', 'asv62') ~ 'spring_bloom',
                                        TRUE ~ 'winter_bloom')) |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(decimal_date, max_abund))+
  # scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
  #                  #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
  #                  #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  # )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.25))+
  #geom_area(aes(date, abundance_value, fill = order_f, group = asv_num), alpha = 0.8,  position='stack')+
  geom_area(aes(day_of_year, abundance_value, fill = family_f), alpha = 0.5, position = 'stack')+
  #geom_smooth(method = 'loess', se = F, aes(day_of_year, abundance_value, group = asv_num, color = family_f))+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  facet_grid(vars(seasonality))+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

wavelets_result_ed_tibble_tax_02_biased_red |>
  dplyr::filter(asv_num %in% c('asv15', 'asv7', 'asv17', 'asv62')) |>
  dplyr::mutate(seasonality = case_when(asv_num %in% c('asv15', 'asv7', 'asv62') ~ 'spring_bloom',
                                        TRUE ~ 'winter_bloom')) |>
  dplyr::filter(!wavelets_transformation %in% c('d1','d2', 'd4')) |>
  ggplot(aes(decimal_date, wavelets_result_ed))+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  #scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.25))+
  #geom_area(aes(date, abundance_value, fill = order_f, group = asv_num), alpha = 0.8,  position='stack')+
  geom_line(aes(group = asv_num, color = seasonality))+
  geom_point(aes(shape = asv_num, color = seasonality))+
  facet_grid(vars(wavelets_transformation))+
  #scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  


###table with the label corresponding to their signals for FL and PA -------

### for each transformation used each bloomer has a different label only if it gives a clear strong signal 
#### this table will contain the asv num + for each transformation which label did it got
###The most important signals for our dataset have been the d1, d3, and s5 they are the ones that we will use for 
### labeling the bloomers.

bloo_all_types_summary_tb <- bloo_all_types_summary |>
  group_by(recurrency, frequency, type_of_bloom, clustering_group, fraction, occurrence_category) %>%
  mutate(group_asv = paste(asv_num, collapse = ", ")) %>%
  ungroup() |>
  dplyr::select(-asv_num, -occurrence_perc) |>
  distinct(recurrency, frequency, type_of_bloom, clustering_group, fraction, occurrence_category, group_asv)

#write.csv(bloo_all_types_summary_tb, 'results/tables/bloo_all_types_summary_tb.csv')

## With this analysis we will pick a representative of each group to perform the random tree analysis-----
## asvs with less clear clutser asv179, asv11, asv385, asv27, asv17, asv43, asv118, asv126, asv28
wavelets_result_ed_tibble_tax_02_biased_red_coeff |>
  dplyr::filter(asv_num %in% c('asv11', 'asv27', 'asv38'))

wavelets_result_ed_tibble_tax_02_biased_red_coeff |>
  dplyr::filter(asv_num %in% c('asv17', 'asv62', 'asv1'))

bloo_3$value

# unclear group
#asv11, asv38, asv27, asv17, asv62, asv1

bloo_02_types_summary <- bloo_02 |>
  as_tibble() |>
  rename('asv_num' = 'value') |>
  dplyr::mutate(
    #ocurrence
    #mean_relative_abundance
    #sd_relative_abundance
    recurrency = case_when(asv_num %in%  c('asv58', 'asv178') ~ 'yes',
                           asv_num %in% c('asv555', 'asv114', 'asv249', 'asv237', 'asv563', 'asv282') ~ 'no',
                           asv_num %in% c('asv15', 'asv7')  ~ 'yes',
                           asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')  ~ 'no',
                           asv_num == 'asv11'~ 'no',
                           asv_num == 'asv27' ~ 'yes',
                           asv_num == 'asv38' ~ 'yes',
                           asv_num == 'asv1' ~ 'yes',
                           asv_num == 'asv17' ~ 'no',
                           asv_num == 'asv62' ~ 'yes'
    ),
    frequency = case_when(#asv_num %in% c('asv15', 'asv7', 'asv16', 'asv116', 'asv182', 'asv84', 'asv100', 'asv25', 'asv72', 'av42') ~ 'seasonal',
      asv_num %in% c('asv58', 'asv178') ~ 'stochastic',
      asv_num %in% c('asv555', 'asv114', 'asv249', 'asv237', 'asv563', 'asv282') ~ 'stochastic',
      asv_num %in% c('asv15', 'asv7') ~ 'seasonal',
      asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')  ~ 'stochastic',
      asv_num == 'asv11' ~ 'stochastic',
      asv_num == 'asv27' ~ 'seasonal',
      asv_num == 'asv38' ~ 'stochastic', ## the coeff are almost the same for d1 and d3 6.81 vs 6.75 respectively, so it has also a seasonal signal
      asv_num == 'asv1' ~ 'seasonal',
      asv_num == 'asv17' ~ 'stochastic',
      asv_num == 'asv62' ~ 'seasonal',
    ) ## our type of data is not adequate to define the length of these blooming events
    # type_of_bloom = case_when(asv_num %in% c('asv58', 'asv178') ~ 'persistent',
    #                           asv_num %in% c('asv555', 'asv114', 'asv249', 'asv237', 'asv563', 'asv282') ~ 'ephemeral',
    #                           asv_num %in% c('asv15', 'asv7')  ~ 'persistent',
    #                           asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'persistent',
    #                           asv_num == 'asv11' ~ 'ephemeral',
    #                           asv_num == 'asv27' ~ 'persistent',
    #                           asv_num == 'asv38' ~ 'persistent',
    #                           asv_num == 'asv1' ~ 'persistent',
    #                           asv_num == 'asv17' ~ 'persistent',
    #                           asv_num == 'asv62' ~ 'persistent',
    )
    
bloo_02_clustering_results <- bloo_02 |>
  as_tibble() |>
  rename('asv_num' = 'value') |>
  dplyr::mutate(
    clustering_group = case_when(asv_num %in% c('asv58', 'asv178') ~ 'cl_1',
                                 asv_num %in% c('asv555', 'asv114', 'asv249', 'asv237', 'asv563', 'asv282') ~ 'cl_2',
                                 asv_num %in% c('asv15', 'asv7')  ~ 'cl_3',
                                 asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'cl_4', ## we need to remove this cluster
                                 asv_num == 'asv11' ~ 'unclear',
                                 asv_num == 'asv27' ~ 'unclear',
                                 asv_num == 'asv38' ~ 'unclear',
                                 asv_num == 'asv1' ~ 'unclear',
                                 asv_num == 'asv17' ~ 'unclear',
                                 asv_num == 'asv62' ~ 'unclear'))  |>
  dplyr::mutate(fraction ='0.2')

bloo_02_types_summary <- bloo_02_types_summary |>
  dplyr::left_join(occurence_perc_02) |> ##categories 1/3 narrow, 1-2/3 intermidiate, 2/3 broad
  dplyr::mutate(fraction = '0.2')

## clusters meaning ----
labs_clusters_fl <- as_labeller(c(cl_1 = 'recurrent random', 
                                  cl_2 = 'ephemeral random',
                                  cl_3 = 'seasonal 1',
                                  cl_4 = 'SAR11 cluster',
                                  unclear = 'ungrouped'
))

### GENERAL TABLES FOR PA AND FL TYPES OF BLOOMS AND CLUSTERING RESULTS-------

## Clustering results at distance 10 in d1, d3, s4
bloo_clustering_results <- bloo_3_clustering_results |>
  bind_rows(bloo_02_clustering_results)

#write.csv(bloo_clustering_results, 'results/tables/clustering_results_d10_d12s4.csv')  

## types of blooms FL and PA together
bloo_types_summary <- bloo_3_types_summary |>
  bind_rows(bloo_02_types_summary)

#write.csv(bloo_types_summary, 'results/tables/summary_types_of_blooms.csv')  

### PLOT THE RESULTS OF THE CLUSTERING ANALYSIS ---------
bloo_clustering_results <- read.csv('results/tables/clustering_results_d10_d12s4.csv')|>
  dplyr::mutate(fraction = as.character(fraction))

bloo_all_types_summary_tax <- bloo_all_types_summary_tax |>
  dplyr::mutate(fraction = as.character(fraction))

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  group_by(asv_num) |>
  left_join(bloo_all_types_summary_tax) |>
  left_join(bloo_clustering_results) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  group_by(date, fraction, clustering_group) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(recurrency~occurrence_category, scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## recurrent blooms
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  group_by(asv_num) |>
  left_join(bloo_02_types_summary) |>
  left_join(bloo_02_clustering_results) |>
  dplyr::filter(recurrency == 'yes') |>
  group_by(date, fraction, recurrency, occurrence_category) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, abundance_value, fill = family_f), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  facet_wrap(recurrency~occurrence_category, scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## seasonal blooms
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  group_by(asv_num) |>
  left_join(bloo_02_types_summary) |>
  left_join(bloo_02_clustering_results) |>
  dplyr::filter(recurrency == 'yes' &
                  frequency == 'seasonal') |>
  group_by(day_of_year, fraction, year) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  group_by(day_of_year, fraction, year, clustering_group, max_abund) |>
  dplyr::summarize(cluster_abund = sum(abundance_value)) |>
  ggplot(aes(day_of_year, max_abund))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(y = cluster_abund, fill = clustering_group), alpha = 0.7,  position='stack')+
  scale_fill_manual(values = palette_clustering)+
  facet_grid(vars(year), scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  group_by(asv_num) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  group_by(asv_num) |>
  left_join(bloo_02_types_summary) |>
  left_join(bloo_02_clustering_results) |>
  dplyr::filter(recurrency == 'yes' &
                  frequency == 'seasonal') |>
  #group_by(day_of_year, fraction, year) |>
  #dplyr::mutate(max_abund = sum(abundance_value)) |>
  #group_by(day_of_year, fraction, year, clustering_group, max_abund) |>
  #dplyr::summarize(cluster_abund = sum(abundance_value)) |>
  ggplot(aes(day_of_year, abundance_value, color = case_when(clustering_group != 'unclear' ~ clustering_group,
                                                             clustering_group == 'unclear' ~ asv_num)))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  #geom_area(aes(y = cluster_abund, fill = clustering_group), alpha = 0.7,  position='stack')+
  geom_point(aes())+
  scale_color_manual(values = palette_clustering)+
  geom_smooth(method = 'loess', aes(color = case_when(clustering_group != 'unclear' ~ clustering_group,
                                                      clustering_group == 'unclear' ~ asv_num)), span = 0.6)+
  #facet_grid(vars(year), scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order', color = 'Seasonal group')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

##seasonal from PA
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '3') |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  group_by(asv_num) |>
  left_join(bloo_3_types_summary) |>
  left_join(bloo_3_clustering_results) |>
  dplyr::filter(recurrency == 'yes' &
                  frequency == 'seasonal') |>
  ggplot(aes(day_of_year, abundance_value, color = case_when(clustering_group != 'unclear' ~ clustering_group,
                                                             clustering_group == 'unclear' ~ asv_num)))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  #geom_area(aes(y = cluster_abund, fill = clustering_group), alpha = 0.7,  position='stack')+
  geom_point(aes())+
  scale_color_manual(values = palette_clustering)+
  geom_smooth(method = 'loess', aes(color = case_when(clustering_group != 'unclear' ~ clustering_group,
                                                      clustering_group == 'unclear' ~ asv_num)), span = 0.6)+
  #facet_grid(vars(year), scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order', color = 'Seasonal group')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))  

## Non recurrent blooms
asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_02$value) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'relative_abundance' &
                  fraction == '0.2') |>
  left_join(bloo_02_types_summary) |>
  left_join(bloo_02_clustering_results) |>
  group_by(asv_num) |>
  dplyr::filter(recurrency == 'no') |>
  group_by(date, fraction) |>
  dplyr::mutate(max_abund = sum(abundance_value)) |>
  group_by(date, fraction, clustering_group, max_abund, asv_num) |>
  dplyr::summarize(cluster_abund = sum(abundance_value)) |>
  ggplot(aes(date, max_abund))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, cluster_abund, fill = case_when(clustering_group != 'unclear' ~ clustering_group,
                                                      clustering_group == 'unclear' ~ asv_num)), alpha = 0.8,  position='stack')+
  scale_fill_manual(values = palette_clustering)+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Order')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5)) 

## Plot all the clusters ----
## general plot by clusters PA and FL
bloo_clustering_results <- bloo_clustering_results |>
  dplyr::mutate(cluster_fr = paste0(clustering_group, '_', fraction))

# labs_clusters_pa_fl <-   as_labeller(c(cl_1_3 = 'habor restoration 1', ## clusters meaning
#                                        cl_2_3 = 'habor restoration 2',
#                                        cl_3_3 = 'habor restoration 3',
#                                        cl_4_3 = 'harbor restoration 4',
#                                        cl_5_3 = 'seasonal 1',
#                                        cl_6_3 = 'harbor restoration 5', 
#                                        cl_7_3 = 'recurrent random',
#                                        cl_8_3 = 'seasonal 2',
#                                        cl_9_3 = 'seasonal 3',
#                                        cl_10_3 = 'harbor restoration 6',
#                                        cl_11_3 = 'harbor restoration 7',
#                                        unclear_3 = 'ungruped', 
#                                        cl_1_0.2 = 'recurrent random', 
#                                        cl_2_0.2 = 'ephemeral random',
#                                        cl_3_0.2 = 'seasonal 1',
#                                        cl_4_0.2 = 'SAR11 cluster',
#                                        unclear_0.2 = 'ungrouped'))

asv_tab_all_bloo_z_tax |>
  colnames()

asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::left_join(bloo_clustering_results, by = c('asv_num' = 'asv_num','fraction')) |>
  dplyr::filter(clustering_group != 'unclear') |>
  group_by(date, clustering_group, fraction, cluster_fr) |>
  dplyr::reframe(abundance_cluster = sum(abundance_value)) |>
  ggplot(aes(date, abundance_cluster))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  scale_x_datetime(expand = c(0,0))+
  geom_area(aes(date, y = abundance_cluster, group = cluster_fr, fill = cluster_fr), alpha = 1,  position= 'stack')+
  scale_fill_manual(values = palette_clustering)+ #, labels = labs_clusters_pa_fl
  facet_wrap(vars(fraction), labeller = labs_fraction)+
  labs(x = 'Day of the year', y = 'Relative abundance (%)', fill = 'Clustering group')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         color = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5))

labs_clusters_pa_fl <-   as_labeller(c(cl_1_3 = 'Habor Restoration 1', ## clusters meaning
                                       cl_2_3 = 'Habor Restoration 2',
                                       cl_3_3 = 'Habor restoration 3',
                                       cl_4_3 = 'Harbor restoration 4',
                                       cl_5_3 = 'Seasonal 1',
                                       cl_6_3 = 'Harbor restoration 5', 
                                       cl_7_3 = 'Recurrent Stochastic',
                                       cl_8_3 = 'Seasonal 2',
                                       cl_9_3 = 'Seasonal 3',
                                       cl_10_3 = 'Harbor restoration 6',
                                       cl_11_3 = 'Harbor restoration 7',
                                       unclear_3 = 'Ungruped', 
                                       cl_1_0.2 = 'Recurrent Stochastic', 
                                       cl_2_0.2 = 'Ephemeral Stochastic',
                                       cl_3_0.2 = 'Seasonal 1',
                                       cl_4_0.2 = 'SAR11 cluster',
                                       unclear_0.2 = 'Ungrouped',
                                       "stochastic_3" = 'Stochastic unclear' , 
                                       "stochastic_0.2" = 'Stochastic unclear',
                                       "seasonal_3" = 'Sesonal unclear' , 
                                       "seasonal_0.2" = 'Seasonal unclear',
                                       '0.2' = 'Free living (0.2-3 um)',
                                       '3' = 'Particle attached (3-20 um)'
                                       
))

bloo_clustering_results <- bloo_clustering_results |>
  dplyr::mutate(fraction = as.character(fraction))

bloo_clustering_results |>
  colnames()

asv_tab_all_bloo_z_tax_cl <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::left_join(bloo_types_summary, by = c('asv_num' = 'asv_num','fraction')) |>
  left_join(bloo_clustering_results) |>
  group_by(date, clustering_group, fraction,cluster_fr) |>
  dplyr::mutate(abundance_cluster = sum(abundance_value)) |>
  dplyr::mutate(facet_variable = case_when(
    !str_detect(cluster_fr, 'unclear') ~ cluster_fr,
    str_detect(cluster_fr, 'unclear') ~ paste0(frequency, '_', fraction) 
  )) 

asv_tab_all_bloo_z_tax_cl$facet_variable |>
  unique()

asv_tab_all_bloo_z_tax_cl$facet_variable <- factor(asv_tab_all_bloo_z_tax_cl$facet_variable, levels = c(     
                                                                                                 "cl_5_3" , "cl_8_3", "cl_9_3", "cl_3_0.2"  ,   ##seasonal
                                                                                                 "cl_4_0.2"  , "seasonal_3" , "seasonal_0.2",
                                                                                                    "cl_7_3"  ,   "cl_1_0.2" ,  "cl_2_0.2" , #stochastic
                                                                                                 "cl_1_3",    "cl_2_3", "cl_3_3" , 
                                                                                                 "cl_4_3" , "cl_6_3",
                                                                                                 "cl_10_3" ,   "cl_11_3" , #harbor restoration clusters
                                                                                                 "stochastic_3" , "stochastic_0.2" 
                                                                                                ))
  
bloo_bbmo_clusters <- asv_tab_all_bloo_z_tax_cl |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, abundance_cluster))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = family_f),   position= 'stack')+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(facet_variable~fraction, 
             labeller = labs_clusters_pa_fl, ncol = 3)+ #, labels = labs_fraction
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         color = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
         legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'bloo_bbmo_clusters.pdf', plot = bloo_bbmo_clusters,
#        path = 'results/figures/',
#        width = 188, height = 220, units = 'mm')

## clusters without SAR11
bloo_bbmo_clusters_no_sar_cluster <- asv_tab_all_bloo_z_tax_cl |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  #dplyr::filter(order != 'SAR11 clade') |>
  dplyr::filter(cluster_fr != 'cl_4_0.2') |> #no sar11 clade but the others yes
  ggplot(aes(date, abundance_cluster))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = family_f),   position= 'stack')+
  #scale_fill_manual(values = palette_clustering, labels = labs_clusters_pa_fl)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(facet_variable~fraction, 
             labeller = labs_clusters_pa_fl, ncol = 3)+ #, labels = labs_fraction
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         color = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'bloo_bbmo_clusters_no_sar_cluster.pdf', plot = bloo_bbmo_clusters_no_sar_cluster,
#        path = 'results/figures/',
#        width = 188, height = 200, units = 'mm')

## plot rclr not relative abundance ----
bloo_clustering_results %$%
  cluster_fr |>
  unique()

labs_clusters_pa_fl_ed <-   as_labeller(c(  cl_4_3 = 'ASV7, ASV15', 
                                            cl_3_0.2 = 'ASV7, ASV15',
                                            cl_8_3 = 'ASV1, ASV4' ,
                                            cl_9_3 =  'ASV69, ASV225', 
                                            cl_2_0.2 = 'ASV555, ASV114, ASV249, ASV237, ASV563, ASV282',
                                            cl_2_3 = 'ASV194, ASV559', 
                                            cl_1_0.2 = 'ASV58, ASV178', 
                                            cl_11_3 = 'ASV264, ASV471, ASV752', 
                                            cl_5_3 =  'ASV317, ASV200, ASV113',
                                            cl_7_3 = 'ASV100, ASV72, ASV42',  
                                            cl_6_3 = 'ASV116, ASV182, ASV84', 
                                            cl_1_3 = 'ASV302, ASV511, ASV311', 
                                            cl_10_3 = 'ASV163, ASV219',
                                            cl_3_3 = 'ASV276, ASV223',  
                                       '0.2' = 'Free living (0.2-3 um)',
                                       '3' = 'Particle attached (3-20 um)'
                                       
))

bloo_types_summary <- bloo_all_types_summary_tax |>
  dplyr::mutate(fraction = as.character(fraction))

asv_tab_all_bloo_z_tax_cl <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::left_join(bloo_types_summary, by = c('asv_num' = 'asv_num','fraction')) |>
  left_join(bloo_clustering_results) |>
  dplyr::filter(clustering_group != 'unclear') |>
  group_by(date, clustering_group, fraction, cluster_fr) |>
  dplyr::mutate(abundance_cluster = sum(abundance_value))

asv_tab_all_bloo_z_tax_cl$cluster_fr |>
  unique()
asv_tab_all_bloo_z_tax_cl |>
  group_by(cluster_fr) |>
  dplyr::reframe(max = max(abundance_value)) |>
  arrange(max)

asv_tab_all_bloo_z_tax_cl$cluster_fr <- factor(asv_tab_all_bloo_z_tax_cl$cluster_fr, levels = c(     
  "cl_8_3" ,   "cl_9_3" , 
   "cl_4_3",  "cl_3_0.2",
  "cl_7_3",  "cl_2_3"  , 
  "cl_5_3"  , "cl_1_0.2", 
  "cl_2_0.2",  "cl_6_3" ,
  "cl_10_3" ,  "cl_3_3" ,
   "cl_1_3" ,"cl_11_3"
))

# bloo_bbmo_clusters <- asv_tab_all_bloo_z_tax_cl |>  
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   ggplot(aes(date, abundance_cluster))+
#   scale_y_continuous( expand = c(0,0))+ #, limits = c(0,0.25)
#   geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = family_f),   position= 'stack')+
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
#   facet_wrap(facet_variable~fraction, 
#              labeller = labs_clusters_pa_fl, ncol = 3)+ #, labels = labs_fraction
#   labs(x = 'Date', y = 'rCLR', fill = 'Family')+
#   guides(fill = guide_legend(ncol = 7, size = 8,
#                              override.aes = aes(label = '')),
#          color = 'none',
#          alpha = 'none')+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
#         legend.position = 'bottom', axis.text.y = element_text(size = 4),
#         axis.title = element_text(size = 7), strip.background = element_blank(), 
#         legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
#         plot.margin = margin(2,5,0,5),
#         legend.key.size = unit(3, 'mm'))

asv_tab_all_bloo_z_tax_cl |>
  colnames()

bloo_bbmo_clusters_no_sar_cluster <- asv_tab_all_bloo_z_tax_cl |> 
  #dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  #dplyr::filter(order != 'SAR11 clade') |>
  dplyr::filter(cluster_fr != 'cl_4_0.2') |> #no sar11 clade but the others yes c('asv2', 'asv3', 'asv5', 'asv8'))
  dplyr::filter(clustering_group != 'unclear') |>

  ggplot(aes(date, abundance_cluster))+
  scale_y_continuous(expand = c(0,0))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = family_f),   position= 'stack')+
  #scale_fill_manual(values = palette_clustering, labels = labs_clusters_pa_fl)+
  geom_hline(yintercept = 0, alpha = 0.2)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(cluster_fr~fraction, 
             labeller = labs_clusters_pa_fl_ed, ncol = 2)+ #, labels = labs_fraction
  labs(x = 'Date', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 1, size = 8,
                             override.aes = aes(label = '')),
         color = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 6), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 6, margin = margin(1, 2, 2, 2)),
        legend.position = 'right', axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 6), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'))

bloo_bbmo_clusters_no_sar_cluster

# ggsave(filename = 'bloo_bbmo_clusters_no_sar_cluster_rclr_no_unclear.pdf', plot = bloo_bbmo_clusters_no_sar_cluster,
#        path = 'results/figures/',
#        width = 188, height = 200, units = 'mm')

## Is SAR11 clade a real bloom or are they responding to compositional constrictions? -----
asv_tab_all_bloo_z_tax_cl |>
  colnames()

asv15_seq <- common_bloomers_tax |>
  dplyr::filter(value == 'asv15') |>
  dplyr::select(seq) |>
  as.vector()

asv2_seq <- common_bloomers_tax |>
  dplyr::filter(value == 'asv2') |>
  dplyr::select(seq) |>
  as.vector()

bloo_bbmo_sar11_cluster <- asv_tab_all_bloo_z_tax_cl |> 
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(blooms = case_when(cluster_fr == 'cl_4_0.2' ~ 'SAR11_cluster',
                                   cluster_fr != 'cl_4_0.2' ~ 'Other blooms')) |>
  dplyr::group_by(blooms, date) |>
  dplyr::reframe(blooms_abund = sum(abundance_value)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, blooms_abund))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.75))+ #, limits = c(0,0.25)
  scale_color_manual(values = c('SAR11_cluster' = '#C73F4E', 'Other blooms' = 'black'))+
  scale_linetype_discrete()+
  geom_line(aes(group = blooms, color = blooms, linetype = blooms), linewidth = 0.25)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = '', linetype = '')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

bloo_bbmo_sar11_cluster 
# ggsave(filename = 'bloo_bbmo_sar11_cluster.pdf', plot = bloo_bbmo_sar11_cluster ,
#        path = 'results/figures/',
#        width = 88, height = 60, units = 'mm')

m_02_HNA <- m_02 |>
  dplyr::select(date, LNA, HNA) |>
  dplyr::mutate(perc_HNA = HNA/(HNA+LNA))|>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

bloo_bbmo_sar11_cluster <- asv_tab_all_bloo_z_tax |> 
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(blooms = case_when(asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'SAR11_cluster',
                                   !asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'Other blooms')) |>
  dplyr::group_by(blooms, date) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::reframe(blooms_abund = sum(abundance_value)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, blooms_abund))+
  #geom_area(data = m_02_HNA, aes(date, perc_HNA))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.75))+ #, limits = c(0,0.25)
  scale_color_manual(values = c('SAR11_cluster' = '#C73F4E', 'Other blooms' = 'black'))+
  scale_linetype_discrete()+
  geom_line(data = asv_tab_all_bloo_z_tax |> 
              dplyr::filter(fraction == '0.2') |>
              dplyr::mutate(blooms = case_when(asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'SAR11_cluster',
                                               !asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'Other blooms')) |>
              dplyr::group_by(blooms, date) |>
              dplyr::filter(abundance_type == 'relative_abundance') |>
              dplyr::reframe(blooms_abund = sum(abundance_value)) |>
              dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
              dplyr::filter(blooms == 'SAR11_cluster'),
              aes(group = blooms, color = blooms), linewidth = 0.25, alpha = 0.2)+
  geom_line(aes(group = blooms, color = blooms, linetype = blooms), linewidth = 0.25)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = '', linetype = '')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

bloo_bbmo_sar11_cluster 
# ggsave(filename = 'bloo_bbmo_sar11_cluster.pdf', plot = bloo_bbmo_sar11_cluster ,
#        path = 'results/figures/',
#        width = 88, height = 60, units = 'mm')

sar11_abundace_cluster <- asv_tab_all_bloo_z_tax |> 
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(blooms = case_when(asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'SAR11_cluster',
                                   !asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'Other blooms')) |>
  dplyr::group_by(blooms, date) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::reframe(blooms_abund = sum(abundance_value)) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(blooms == 'SAR11_cluster')

m_02_HNA |>
  dplyr::mutate(date = as.character.Date(date)) |>
  left_join(sar11_abundace_cluster) |>
  ggplot(aes(perc_HNA, blooms_abund))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme(aspect.ratio = 4/4)

m_02_HNA |>
  dplyr::mutate(date = as.character.Date(date)) |>
  left_join(sar11_abundace_cluster) |>
  ggplot(aes(LNA, blooms_abund))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme(aspect.ratio = 4/4)

no_sar11_abundace_cluster <- asv_tab_all_bloo_z_tax |> 
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(blooms = case_when(asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'SAR11_cluster',
                                   !asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8') ~ 'Other blooms')) |>
  dplyr::group_by(blooms, date) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::reframe(blooms_abund = sum(abundance_value)) |>
  #dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(blooms != 'SAR11_cluster')

m_02_HNA |>
  dplyr::mutate(date = as.character.Date(date)) |>
  left_join(no_sar11_abundace_cluster) |>
  ggplot(aes(perc_HNA, blooms_abund))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme(aspect.ratio = 4/4)
library(ggpmisc)

hna_blooms <- m_02_HNA |>
  dplyr::mutate(date = as.character.Date(date)) |>
  left_join(no_sar11_abundace_cluster)

hna_blooms |>
  ggplot(aes(HNA, blooms_abund))+
  geom_point()+
  stat_poly_eq(aes(
                   label =  paste(after_stat(p.value.label)), size = 3), 
               #position = 'stack',
               #formula = x ~ y,
               method = stats::lm,
               p.digits = 2, 
               coef.keep.zeros = T, #npcx = 1, 
               #npcy = 1,
               na.rm = FALSE)+
  stat_poly_eq(aes(
                   label =  paste(after_stat(rr.label)), size = 3), 
               #position = 'stack',
               #formula = x ~ y,
               method = stats::lm,
               rr.digits = 2, 
               coef.keep.zeros = T, #npcx = 1, 
               npcy = 50,
               na.rm = FALSE)+
  scale_x_continuous(labels = scientific_format())+
  scale_y_continuous(labels = percent_format())+
  geom_smooth(method = 'lm')+
  theme(aspect.ratio = 4/4)



##not only SAR11 clade but all SAR11
bloo_bbmo_sar11_order <- asv_tab_all_bloo_z_tax_cl |> 
  dplyr::filter(fraction == '0.2') |>
  # dplyr::mutate(blooms = case_when(cluster_fr == 'cl_4_0.2' ~ 'SAR11_cluster',
  #                                  cluster_fr != 'cl_4_0.2' ~ 'Other blooms')) |>
  dplyr::mutate(blooms = case_when(order == 'SAR11 clade' ~ 'SAR11_cluster',
                                   order != 'SAR11 clade' ~ 'Other blooms')) |>
  dplyr::group_by(blooms, date) |>
  dplyr::reframe(blooms_abund = sum(abundance_value)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, blooms_abund))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.75))+ #, limits = c(0,0.25)
  scale_color_manual(values = c('SAR11_cluster' = '#C73F4E', 'Other blooms' = 'black'))+
  scale_linetype_discrete()+
  geom_line(aes(group = blooms, color = blooms, linetype = blooms), linewidth = 0.25)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = '', linetype = '')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

asv_tab_all_bloo_z_tax_cl |> 
  dplyr::filter(fraction == '0.2') |>
  # dplyr::mutate(blooms = case_when(cluster_fr == 'cl_4_0.2' ~ 'SAR11_cluster',
  #                                  cluster_fr != 'cl_4_0.2' ~ 'Other blooms')) |>
  dplyr::mutate(blooms = case_when(order == 'SAR11 clade' ~ 'SAR11 clade',
                                   order != 'SAR11 clade' ~ 'Other blooms')) |>
  dplyr::group_by(blooms, date) |>
  dplyr::reframe(blooms_abund = sum(abundance_value)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, blooms_abund))+
  geom_area(aes(fill = blooms), position = 'stack')+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.75))+ #, limits = c(0,0.25)
  scale_fill_manual(values = c('SAR11 clade' = '#C73F4E', 'Other blooms' = 'black'))+
  scale_linetype_discrete()+
  #geom_line(aes(group = blooms, color = blooms, linetype = blooms), linewidth = 0.25)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = '', linetype = '')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

asv_tab_all_bloo_z_tax_cl |> 
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(blooms = case_when(cluster_fr == 'cl_4_0.2' ~ 'SAR11_cluster',
                                   cluster_fr != 'cl_4_0.2' ~ 'Other blooms')) |>
  # dplyr::mutate(blooms = case_when(order == 'SAR11 clade' ~ 'SAR11_cluster',
  #                                  order != 'SAR11 clade' ~ 'Other blooms')) |>
  dplyr::group_by(blooms, date) |>
  dplyr::reframe(blooms_abund = sum(abundance_value)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, blooms_abund))+
  geom_area(aes(fill = blooms), position = 'stack')+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.75))+ #, limits = c(0,0.25)
  scale_fill_manual(values = c('SAR11_cluster' = '#C73F4E', 'Other blooms' = 'black'))+
  scale_linetype_discrete()+
  #geom_line(aes(group = blooms, color = blooms, linetype = blooms), linewidth = 0.25)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = '', linetype = '')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

asv_tab_all_bloo_z_tax_cl |> 
  dplyr::filter(fraction == '0.2') |>
  dplyr::mutate(blooms = case_when(cluster_fr == 'cl_4_0.2' ~ 'SAR11_cluster',
                                   cluster_fr != 'cl_4_0.2' ~ 'Other blooms')) |>
  # dplyr::mutate(blooms = case_when(order == 'SAR11 clade' ~ 'SAR11_cluster',
  #                                  order != 'SAR11 clade' ~ 'Other blooms')) |>
  dplyr::group_by(blooms, date) |>
  dplyr::reframe(blooms_abund = sum(abundance_value)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, blooms_abund))+
  geom_area(aes(fill = blooms), position = 'stack')+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.75))+ #, limits = c(0,0.25)
  scale_fill_manual(values = c('SAR11_cluster' = '#C73F4E', 'Other blooms' = 'black'))+
  scale_linetype_discrete()+
  #geom_line(aes(group = blooms, color = blooms, linetype = blooms), linewidth = 0.25)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = '', linetype = '')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

bloo_all_types_summary_tax <- bloo_all_types_summary |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'))

### We have different SAR11 the group that is clustered in the SAR11 clade and the ones that are not clustered in there explore and decide if I need to remove all of them or just
### those that were in the SAR11 cluster

## SAR11 clade ASVs study them in detail 
# 1   asv38
# 2    asv8 (cluster SAR11)
# 3  asv225
# 4  asv264
# 5    asv5 (cluster SAR11)
# 6  asv200
# 7    asv3 (cluster SAR11)
# 8    asv2(cluster SAR11)
# 9   asv15

asv_tab_all_bloo_z_tax_cl |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(order == 'SAR11 clade') |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
  #geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = family_f),   position= 'stack')+
  geom_line(aes(y = abundance_value, group = asv_num_f, color = family_f))+
  geom_line(aes(y = bacteria_joint/10^8))+
  #scale_fill_manual(values = palette_clustering, labels = labs_clusters_pa_fl)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(facet_variable~fraction~asv_num, 
             ncol = 3)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl
  labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         color = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

### Seasonal  -----
palete_seasonal_bloo <- c('cl_5_3'= '#2466BF', 'cl_8_3' =  '#5B9B57','cl_9_3' = '#FFD700', 'asv27_3'=  '#C73F4E' , 'asv28_3'= '#AD7CE6',
                          'asv27_0.2'=  '#C73F4E' , 'cl_3_0.2' = '#2466BF', 'asv62_0.2' = '#5B9B57' , 'asv1_0.2'= '#2F2F2C')

seasonal_clusters_labs <- as_labeller(c( cl_5_3 = '1st bloom (ASV15, ASV7)',
                                         cl_8_3 = '2nbloom (ASV72, ASV25, ASV100, ASV42)',
                                         cl_9_3 = '3rd bloom (ASV23, ASV1, ASV4, ASV31)',
                                         asv27_3 = 'ASV27',
                                         asv28_3 = 'ASV28', #between second and third bloom
                                         cl_3_0.2 = '1st bloom (ASV15, ASV7)',
                                         asv62_0.2 = 'ASV62',
                                           asv1_0.2 = 'ASV1',
                                         asv27_0.2   = 'ASV27',
                                         '2004' = '2004',
                                         '2005' = '2005',
                                         '2006' = '2006',
                                         '2007' = '2007', 
                                         '2008' = '2008',
                                         '2009' = '2009',
                                         '2010' = '2010',
                                         '2011' = '2011',
                                         '2012' = '2012',
                                         '2013' = '2013',
                                         '0.2' = 'Free living (0.2-3 um)',
                                         '3' = 'Particle attached (3-20 um)'))

#### add which ASVs belong to each cluster
##geom line for seasonal blooms
asv_tab_all_bloo_z_tax |>
  str()

bloo_all_types_summary_tb <- bloo_all_types_summary_tb_tax_v2 |>
  dplyr::mutate(fraction = as.factor(fraction))

data <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(fraction = as.factor(fraction)) |>
  dplyr::left_join(bloo_all_types_summary_tb, by = c('asv_num_f' = 'asv_num','fraction', 'order', 'class', 'phylum', 'family')) |> 
  dplyr::filter(type_of_bloomer == 'Seasonal' )
# |>
#   dplyr::mutate(facet_variable = case_when(
#     !str_detect(cluster_fr, 'unclear') ~ cluster_fr,
#     str_detect(cluster_fr, 'unclear') ~ paste0(asv_num,'_',fraction))) 

# data |> 
#   arrange(cluster_fr) |>
#   distinct(cluster_fr, asv_num)

# asv_tab_all_bloo_z_tax  |>
#   dplyr::filter(asv_num %in% bloo_3$value &
#                   fraction == '3' |
#                   asv_num %in% bloo_02$value &
#                   fraction == '0.2') |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   dplyr::left_join(bloo_all_types_summary, by = c('asv_num_f' = 'asv_num','fraction')) |> 
#   dplyr::filter(frequency == 'seasonal') |>
#   dplyr::mutate(facet_variable = case_when(
#     !str_detect(cluster_fr, 'unclear') ~ cluster_fr,
#     str_detect(cluster_fr, 'unclear') ~ paste0(asv_num,'_',fraction))) |>
#   distinct(facet_variable)
# 
# data$facet_variable <- data$facet_variable |>
#   factor(levels = c( 'cl_5_3', 'cl_8_3',  'cl_9_3',  
#                      'cl_3_0.2',  'asv27_3', 'asv27_0.2', 
#                      'asv1_0.2', 'asv28_3', 'asv62_0.2'))
  
seasonal_bloo <- data |>
  ggplot(aes(day_of_year, abundance_value))+
  geom_point(aes(color = family_f), alpha = 0.8)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  geom_smooth(method = 'loess', color = 'black', span = 0.7)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, 
  facet_wrap(recurrency~fraction)+ #, labeller = seasonal_clusters_labs 
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance (%)' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(0, 0, 2, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'seasonal_bloo.pdf', plot = seasonal_bloo ,
#        path = 'results/figures/',
#        width = 188, height = 150, units = 'mm')

seasonal_bloo_y <- data |>
  ggplot(aes(day_of_year, abundance_value))+
  geom_point(aes(color = family_f))+
  scale_color_manual(values = palette_family_assigned_bloo)+
  geom_smooth(method = 'loess', color = 'black', span = 0.7)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, 
  facet_grid(asv_num~fraction~year)+ #, labeller = seasonal_clusters_labs
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance (%)' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'seasonal_bloo_y.pdf', plot = seasonal_bloo_y ,
#        path = 'results/figures/',
#        width = 288, height = 250, units = 'mm')

seasonal_bloo_frac <- data |>
  ggplot(aes(day_of_year, abundance_value, group = facet_variable))+
  geom_point(aes(color = facet_variable), alpha = 0.8)+
  #scale_color_manual(values = palette_family_assigned_bloo)+
  scale_color_manual(values = palete_seasonal_bloo, labels = seasonal_clusters_labs)+
  scale_fill_manual(values = palete_seasonal_bloo)+
  geom_smooth(method = 'loess', aes(fill = facet_variable, color = facet_variable), span = 0.7)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, 
  facet_wrap(vars(fraction))+ #, labeller = seasonal_clusters_labs 
  guides(fill = 'none')+
  theme_bw()+
  labs(color = 'Clustering group', x = 'Day of the year', y = 'Relative abundance (%)' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(0, 0, 2, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'seasonal_bloo_frac.pdf', plot = seasonal_bloo_frac,
#        path = 'results/figures/',
#        width = 150, height = 100, units = 'mm')

## I would like to illustrate that seasonal bloomers do not bloom every year -----
data |>
  distinct(asv_num)

data |>
  distinct(asv_num)

data$asv_num_f

data <- data |>
  dplyr::mutate(phylum_f = as_factor(phylum),
                family_f = as_factor(family),
                order_f = as_factor(order),
                class_f = as_factor(class),
                asv_num_f = as_factor(asv_num))

data$class_f <-  factor(data$class_f, 
                                          levels=unique(data$class_f[order(data$phylum_f)]), 
                                          ordered=TRUE)

data$order_f <-  factor(data$order_f, 
                                          levels=unique(data$order_f[order(data$phylum_f,
                                                                                             data$class_f)]), 
                                          ordered=TRUE)

data$family_f <-  factor(data$family_f, 
                                           levels=unique(data$family_f[order(data$phylum_f,
                                                                                               data$class_f,
                                                                                               data$order_f)]), 
                                           ordered=TRUE)

data$asv_num_f <-  factor(data$asv_num_f, 
                                            levels=unique(data$asv_num_f[order(data$phylum_f,
                                                                                                 data$class_f,
                                                                                                 data$order_f,
                                                                                                 data$family_f)]), 
                                            ordered=TRUE)

data <- data |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) 

data$facet_labels <- factor(data$facet_labels, 
                            levels = unique(data$facet_labels[order(
                                                                 data$family_f,
                                                                 data$asv_num_f)]), 
                            ordered=TRUE)

seasonal_bloo_y <- data |>
  ggplot(aes(day_of_year, abundance_value, color = family_f))+
  geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_point(aes(color = family_f), size = 0.2)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = fraction, linetype = fraction), linewidth = 1)+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype(labels = labs_fraction)+
  labs(linetype = 'Fraction', color = 'Family')+
  facet_grid(facet_labels ~ year, scales = 'free_y')+
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance' )+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text.x = element_text(size = 7, margin = margin(5, 5, 5, 5)),
        strip.text.y = element_text(size = 5, margin = margin(1, 1, 1, 1)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 10), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))+
  guides(linetype = guide_legend(ncol = 1))

seasonal_bloo_y

ggsave(filename = 'results/figures/seasonal_bloo_y_fraction_ed2.pdf', plot = seasonal_bloo_y,
      # path = 'results/figures/',
       width = 180, height = 200, units = 'mm')

## I would like to plot the chl-a concentrations 
chla_years <- data |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) |>
  dplyr::select(chla_total, chla_3um, date, day_of_year, year) |>
  distinct(chla_total, chla_3um, date, day_of_year, year) |>
  pivot_longer(cols = starts_with('chl'), values_to = 'chla_concentrations', names_to = 'chla_fraction') |>
  ggplot(aes(day_of_year, chla_concentrations))+
  geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  #geom_point(aes(color = chla_fraction), size = 0.2)+
  scale_color_manual(values = palette_clustering)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = chla_fraction, linetype = chla_fraction))+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  #scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype()+
  labs(linetype = 'Chl-a', color = 'Chl-a fraction')+
  facet_wrap(vars(year), nrow = 1)+
  theme_bw()+
  labs(color = 'Chl-a', x = 'Day of the year', y = 'Concentration' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text.x = element_text(size = 7, margin = margin(5, 5, 5, 5)),
        strip.text.y = element_text(size = 4, margin = margin(1, 1, 1, 1)),
        aspect.ratio = 5/7,
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 4), legend.title = element_text(size = 6), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))

ggsave(filename = 'results/figures/chla_years.pdf', plot = chla_years,
      # path = 'results/figures/',
       width = 180, height = 40, units = 'mm')


data <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::mutate(fraction = as.factor(fraction)) |>
  dplyr::left_join(bloo_all_types_summary_tb, by = c('asv_num_f' = 'asv_num','fraction')) |> 
  dplyr::filter(frequency == 'seasonal' )

seasonal_bloo_y <- data |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) |>
  ggplot(aes(day_of_year, abundance_value, color = family_f))+
  #geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_point(aes(color = family_f), size = 0.2)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = fraction, linetype = fraction))+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  #scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype(labels = labs_fraction)+
  labs(linetype = 'Fraction', color = 'Family')+
  facet_grid(facet_labels ~ year, scales = 'free_y')+
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(5, 5, 5, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 6), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))

seasonal_bloo_y

# ggsave(filename = '../results/figures/seasonal_bloo_y_fraction_rclr.pdf', plot = seasonal_bloo_y,
#        # path = 'results/figures/',
#        width = 180, height = 200, units = 'mm')


## I do the same plot for my non seasonal bloomers -----

data <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(fraction = as.factor(fraction)) |>
  dplyr::left_join(bloo_all_types_summary_tb, by = c('asv_num_f' = 'asv_num','fraction')) |> 
  dplyr::filter(frequency != 'seasonal' ) |>
  dplyr::filter(!asv_num %in% c('asv2', 'asv3', 'asv5', 'asv8'))

## I divide them because I have too much taxa to be properly observed in one plot----
data_shared <- data |>
  group_by(fraction, asv_num) |>
  distinct(fraction, asv_num) |> 
  group_by(asv_num) |>
  dplyr::filter(n() == 2)

data_broad_shared <- data |>
  dplyr::filter(occurrence_category == 'broad') 

data_intermediate_non_shared <- data |>
  dplyr::filter(occurrence_category == 'intermediate') |>
  dplyr::filter(!asv_num %in% data_shared$asv_num)

data_nar_non_shared <- data |>
  dplyr::filter(occurrence_category == 'narrow') |>
  dplyr::filter(!asv_num %in% data_shared$asv_num)

# non_seasonal_bloo_y_shared <- data |>
#   dplyr::filter(asv_num %in% data_shared$asv_num) |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) |>
#   ggplot(aes(day_of_year, abundance_value, color = family_f))+
#   geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
#   geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
#   geom_point(aes(color = family_f), size = 0.2)+
#   scale_color_manual(values = palette_family_assigned_bloo)+
#   # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
#   geom_line(aes(group = fraction, linetype = fraction))+
#   scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
#   scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
#   scale_linetype(labels = labs_fraction)+
#   labs(linetype = 'Fraction', color = 'Family')+
#   facet_grid(facet_labels ~ year, scales = 'free_y')+
#   theme_bw()+
#   labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance (%)' )+
#   theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(5, 5, 5, 5)),
#         legend.position = 'bottom', axis.text.y = element_text(size = 4),
#         axis.title = element_text(size = 6), strip.background = element_blank(), #element_rect(fill = 'transparent')
#         legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
#         plot.margin = margin(10,10,10,10),
#         legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))
# 
# non_seasonal_bloo_y_shared
# ggsave(filename = 'results/figures/non_seasonal_bloo_y_shared.pdf', plot = non_seasonal_bloo_y_shared,
#        # path = 'results/figures/',
#        width = 180, height = 100, units = 'mm')

non_seasonal_bloo_y_inter <- data |>
  dplyr::filter(asv_num %in% data_intermediate_non_shared$asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) |>
  ggplot(aes(day_of_year, abundance_value, color = family_f))+
  geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_point(aes(color = family_f), size = 0.2)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = fraction, linetype = fraction))+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype(labels = labs_fraction)+
  labs(linetype = 'Fraction', color = 'Family')+
  facet_grid(facet_labels ~ year, scales = 'free_y')+
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(5, 5, 5, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 6), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))

non_seasonal_bloo_y_inter

# ggsave(filename = 'results/figures/non_seasonal_bloo_y_inter.pdf', plot = non_seasonal_bloo_y_inter,
#        # path = 'results/figures/',
#        width = 180, height = 120, units = 'mm')

non_seasonal_bloo_y_borad <- data |>
  dplyr::filter(asv_num %in% data_broad_shared$asv_num) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) |>
  ggplot(aes(day_of_year, abundance_value, color = family_f))+
  geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_point(aes(color = family_f), size = 0.2)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = fraction, linetype = fraction))+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  scale_y_continuous(expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype(labels = labs_fraction)+
  labs(linetype = 'Fraction', color = 'Family')+
  facet_grid(facet_labels ~ year, scales = 'free_y')+
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(5, 5, 5, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 6), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))

non_seasonal_bloo_y_borad

# ggsave(filename = 'results/figures/non_seasonal_bloo_y_borad_shared.pdf', plot = non_seasonal_bloo_y_borad,
#        # path = 'results/figures/',
#        width = 180, height = 80, units = 'mm')

### i divide the narrow group because i have too many to observe them well
unique_asv <- data_nar_non_shared$asv_num |>
  unique()

# Determine the size of each chunk
chunk_size <- ceiling(length(unique_asv) / 3)

# Split unique_asv into three chunks
chunks <- split(unique_asv, cut(seq_along(unique_asv), breaks = 3, labels = FALSE))

# Create three vectors
vector1 <- chunks[[1]]
vector2 <- chunks[[2]]
vector3 <- chunks[[3]]

non_seasonal_bloo_y_nar <- data |>
  dplyr::filter(asv_num %in% vector1) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) |>
  ggplot(aes(day_of_year, abundance_value, color = family_f))+
  geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_point(aes(color = family_f), size = 0.2)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = fraction, linetype = fraction))+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype(labels = labs_fraction)+
  labs(linetype = 'Fraction', color = 'Family')+
  facet_grid(facet_labels ~ year, scales = 'free_y')+
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(5, 5, 5, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 6), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))

non_seasonal_bloo_y_nar

# ggsave(filename = 'results/figures/non_seasonal_bloo_y_nar1.pdf', plot = non_seasonal_bloo_y_nar,
#        # path = 'results/figures/',
#        width = 180, height = 220, units = 'mm')

non_seasonal_bloo_y_nar <- data |>
  dplyr::filter(asv_num %in% vector2) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) |>
  ggplot(aes(day_of_year, abundance_value, color = family_f))+
  geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_point(aes(color = family_f), size = 0.2)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = fraction, linetype = fraction))+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype(labels = labs_fraction)+
  labs(linetype = 'Fraction', color = 'Family')+
  facet_grid(facet_labels ~ year, scales = 'free_y')+
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(5, 5, 5, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 6), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))

non_seasonal_bloo_y_nar

# ggsave(filename = 'results/figures/non_seasonal_bloo_y_nar2.pdf', plot = non_seasonal_bloo_y_nar,
#        # path = 'results/figures/',
#        width = 180, height = 220, units = 'mm')

non_seasonal_bloo_y_nar <- data |>
  dplyr::filter(asv_num %in% vector3) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(facet_labels = paste(asv_num_f, "\n", family_f)) |>
  ggplot(aes(day_of_year, abundance_value, color = family_f))+
  geom_hline(yintercept = 0.1, color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_vline(xintercept = c(80,173,267,354), color = '#BFBFC3', linewidth = 0.2, alpha = 0.7)+
  geom_point(aes(color = family_f), size = 0.2)+
  scale_color_manual(values = palette_family_assigned_bloo)+
  # geom_smooth(aes(group = fraction, linetype = fraction), method = 'loess', color = 'black', span = 0.7)+
  geom_line(aes(group = fraction, linetype = fraction))+
  scale_x_continuous(expand = c(0,0), breaks = c(80,173,267,354))+
  scale_y_continuous( expand = c(0,0), labels = function(x) round(x, 1), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+ #, labels = percent_format(),
  scale_linetype(labels = labs_fraction)+
  labs(linetype = 'Fraction', color = 'Family')+
  facet_grid(facet_labels ~ year, scales = 'free_y')+
  theme_bw()+
  labs(color = 'Family', x = 'Day of the year', y = 'Relative abundance' )+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(5, 5, 5, 5)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 6), strip.background = element_blank(), #element_rect(fill = 'transparent')
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(3, 'mm'),panel.spacing.x = unit(0.2, "lines"))

non_seasonal_bloo_y_nar

# ggsave(filename = 'results/figures/non_seasonal_bloo_y_nar3.pdf', plot = non_seasonal_bloo_y_nar,
#        # path = 'results/figures/',
#        width = 180, height = 220, units = 'mm')

## taxonomy of my seasonal blooms------
data |>
  colnames()

data |>
  distinct(clustering_group)

labs_clustering_tax <- as_labeller(c(cl_3 = '1st bloom (ASV15, ASV7)',
                                  unclear = 'Unclear',
                                  cl_5 = '1st bloom (ASV15, ASV7)',
                                  cl_8 = '2nbloom (ASV72, ASV25, ASV100, ASV42)',
                                  cl_9 = '3rd bloom (ASV23, ASV1, ASV4, ASV31)' ))

seasonal_tax <- data |>
  group_by(asv_num, clustering_group, family_f, fraction) |>
  distinct(asv_num, facet_variable) |>
  group_by(clustering_group, family_f,   fraction) |>
  dplyr::reframe(n = n()) |>
  ungroup() |>
  ggplot(aes(clustering_group, n))+
  geom_col(aes(fill = family_f))+
  scale_x_discrete(labels = labs_clustering_tax )+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(fill= 'Family', x = 'Custering group', y = 'Number of taxa' )+
  theme_bw()+
  facet_grid(vars(fraction), labeller = labs_fraction, scales = 'free_x')+
  coord_flip()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), strip.text = element_text(size = 4, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4), strip.background = element_blank(), 
        legend.text = element_text(size = 3), legend.title = element_text(size = 4), strip.placement = 'outside',
        plot.margin = margin(2,5,0,0),
        legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'seasonal_tax.pdf', plot = seasonal_tax,
#        path = 'results/figures/',
#        width = 88, height = 88, units = 'mm')

## how many different types of blooming events do i have?----
distinct_bloo_events <- bloo_all_types_summary |>
  group_by(recurrency, occurrence_category, frequency, type_of_bloom) |>
  distinct(recurrency, occurrence_category, frequency, type_of_bloom) |>
  arrange(occurrence_category)

#install.packages("ggpattern")
library(ggpattern)

#write.csv(bloo_all_types_summary, 'results/tables/bloo_all_types_summary.csv')

bloo_all_types_summary |>
  dplyr::mutate(present = 1) |>
  ungroup() |>
  dplyr::select(-clustering_group, -recurrency, -type_of_bloom, -frequency, -fraction) |>
  pivot_wider( values_from = present, names_from = asv_num)

bloo_all_types_summary |>
  dplyr::mutate(present = 1) |>
  ungroup() |>
  dplyr::select(-asv_num) |>
  distinct(cluster_fr, recurrency, type_of_bloom, frequency, occurrence_category) |>
  pivot_wider( values_from = present, names_from = asv_num)

bloo_all_types_summary <- bloo_all_types_summary |>
  dplyr::mutate(type_of_bloom = str_replace(type_of_bloom, 'ephimeral', 'ephemeral'))

bloo_all_types_summary_tab <- bloo_all_types_summary |>
  group_by(recurrency, occurrence_category, frequency, type_of_bloom, fraction) |>
  dplyr::reframe(n = n()) |>
  arrange(occurrence_category) |>
  pivot_wider(id_cols = c(recurrency, occurrence_category, frequency, type_of_bloom), names_from = fraction, values_from = n, values_fill = 0)

#write.csv(bloo_all_types_summary_tab, 'results/tables/bloo_all_types_summary.csv')
# 
# bloo_all_types_summary |>
#   ggplot(aes(type_of_bloom, occurrence_category, pattern = recurrency))+
#   scale_pattern_manual(values = c(
#     "Texture1" = "stripe",
#     "Texture2" = "crosshatch"
#   ))+
#   geom_tile(aes(fill = n))+
#   facet_grid(frequency~fraction)

####stream graph-----
#remotes::install_github("hrbrmstr/streamgraph")

# library(streamgraph)
# data <- asv_tab_all_bloo_z_tax |>
#   dplyr::filter(asv_num %in% bloo_02$value) |>
#   dplyr::filter(abundance_type == 'relative_abundance' &
#                   fraction == '0.2') |>
#   group_by(asv_num) |>
#   dplyr::mutate(relative_abundance_category = case_when(mean(abundance_value) > 0.01 ~ 'abundant',
#                                                         mean(abundance_value) < 0.01 ~ 'mid',
#                                                         mean(abundance_value) < 0.001 ~ 'rare')) |>
#   dplyr::left_join(bloo_02_types_summary) |>
#   dplyr::filter(recurrency == 'yes' &
#                   frequency == 'seasonal') |>
#   group_by(date, fraction, clustering_group) |>
#   dplyr::summarize(cluster_abund = sum(abundance_value)) 
# 
#   streamgraph(data = data, key = clustering_group, value = cluster_abund, date = date, offset="zero") %>%
#   sg_fill_brewer("BuPu")

#####PLOT FOR EACH EXAMPLE OF TYPE OF BLOOM-------

##before the examples for the presentation of the blooms----
## I reduced my types of blooms to 6

# asv_tab_all_bloo_z_tax_examples <- asv_tab_all_bloo_z_tax |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   dplyr::filter(asv_num %in% bloo_3$value &
#                   fraction == '3' |
#                   asv_num %in% bloo_02$value &
#                   fraction == '0.2') |>
#   dplyr::filter(abundance_type == 'relative_abundance') |>
#   dplyr::left_join(bloo_types_summary, by = c('asv_num_f' = 'asv_num','fraction')) |>
#   dplyr::filter(order != 'SAR11 clade') |>
#   #dplyr::filter(asv_num_f %in% c('asv17', 'asv23', 'asv11', 'asv84', 'asv62', 'asv58', 'asv555', 'asv559')) |>
#   filter((asv_num_f %in% c('asv17', 'asv23', 'asv11', 'asv84', 'asv62', 'asv58', 'asv555', 'asv559') & fraction == '3') |
#            (asv_num_f %in% c('asv62', 'asv58', 'asv555') & fraction == '0.2') |
#            (asv_num_f == 'asv11' & fraction == '0.2')) |>
#   dplyr::mutate(type_of_bloom = str_replace(type_of_bloom, 'ephimeral', 'ephemeral')) |>
#   dplyr::mutate(facet_variable_2 = paste0(occurrence_category, '.', recurrency, '.', frequency, '.', type_of_bloom)) |>
#   dplyr::mutate(facet_var_3 = paste0(occurrence_category, '-', frequency),
#                 facet_var_4 = paste0(type_of_bloom, '-', recurrency)) |>
#   dplyr::mutate(facet_var_5 = paste0( recurrency, '-', frequency),
#                 facet_var_6 = paste0(type_of_bloom)) |>
#   dplyr::mutate(bloom = case_when (z_score_ra >= 1.96 &
#                                      abundance_value >= 0.1~ 'Bloom',
#                                    TRUE ~ 'No-bloom')) |>
#   dplyr::mutate(bloom = as.factor(bloom))
#   
#   asv_summary_bloom <- asv_tab_all_bloo_z_tax_examples |>
#   distinct(asv_num, occurrence_category, type_of_bloom, recurrency, frequency) |>
#   left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f'))
# 
# asv_tab_all_bloo_z_tax_examples$asv_num |>
#   unique() 
# 
# asv_tab_all_bloo_z_tax_examples |>
#   ungroup() |>
#   select(facet_variable_2, asv_num_f) |>
#   distinct()
# 
# asv_tab_all_bloo_z_tax_examples$facet_variable_2  <- factor(asv_tab_all_bloo_z_tax_examples$facet_variable_2, levels = c(     
#   "broad.no.stochastic.persistent" ,
#   "broad.yes.seasonal.persistent" ,
#   "intermidiate.no.stochastic.ephemeral" ,
#   "intermidiate.yes.seasonal.persistent",
#   "intermidiate.yes.stochastic.persistent",
#   "intermidiate.no.stochastic.persistent" ,
#     "narrow.no.stochastic.ephemeral"     ,    
#     "narrow.no.stochastic.persistent"  
# ))
# 
# labels_type_bloom  <- as_labeller(c(     
#   broad.no.stochastic.persistent = 'Broad Stochastic Persistent',
#   broad.yes.seasonal.persistent = 'Broad Recurrent Seasonal Persistent',
#   intermidiate.no.stochastic.ephemeral = 'Intermediate Seasonal Ephemeral',
#   intermidiate.yes.seasonal.persistent = 'Intermediate Recurrent Seasonal Persistent',
#   intermidiate.yes.stochastic.persistent = 'Intermediate Recurrent Stochastic Persistent',
#   intermidiate.no.stochastic.persistent = 'Intermediate Stochastic Persistent',
#   narrow.no.stochastic.ephemeral = 'Narrow Stochastic ephemeral'   ,    
#   narrow.no.stochastic.persistent  = 'Narrow Stochastic persistent'
# ))
# 
# asv_tab_all_bloo_z_tax_examples <- asv_tab_all_bloo_z_tax_examples |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

# bloo_bbmo_clusters_asv <- asv_tab_all_bloo_z_tax_examples |> 
#   # dplyr::filter(facet_variable_2 %in% c(#"broad.no.stochastic.persistent"# ,
#   #   #"narrow.no.stochastic.persistent"
#   #   #
#   #   'intermidiate.yes.stochastic.persistent'
#   #                                       )) |>
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   ggplot(aes(date, abundance_value))+
#   geom_area(aes(date, y = abundance_value, group = fraction, fill = family_f),   position= 'stack')+
#   geom_line(aes(group = fraction), linewidth = 0.25)+
#   scale_linetype_discrete()+
#   # geom_point(data = asv_tab_all_bloo_z_tax_examples |>
#   #              dplyr::filter(z_score_ra >= 1.96 &
#   #                              abundance_value >= 0.1),  
#   #            aes(date, abundance_value, color =  '#9F0011', alpha = 1))+ 
#   scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   facet_wrap(vars(asv_num_f), 
#              ncol = 3, labeller = labels_type_bloom)+ #, labels = labs_fraction
#   scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
#   labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Family')+
#   guides(fill = guide_legend(ncol = 7, size = 8,
#                              override.aes = aes(label = '')),
#          color = 'none',
#          alpha = 'none')+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(), strip.text = element_text(size = 8, margin = margin(0, 0, 2, 0)),
#         legend.position = 'bottom', axis.text.y = element_text(size = 4),
#         axis.title = element_text(size = 7), strip.background = element_blank(), 
#         legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
#         plot.margin = margin(2,5,0,5),
#         legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'bloo_bbmo_clusters.pdf', plot = bloo_bbmo_clusters,
#        path = 'results/figures/',
#        width = 188, height = 220, units = 'mm')

## SAR11 clade ASVs study them in detail 
# asv_num
# 1   asv38
# 2    asv8
# 3  asv225
# 4  asv264
# 5    asv5
# 6  asv200
# 7    asv3
# 8    asv2
# 9   asv15
asv_tab_all_bloo_z_tax |>
  colnames()

1/3

types_of_bloomers

wavelet_results_category

bloo_types_summary |>
  dplyr::filter(recurrency == 'yes' &
                  fraction == '3' &
                  occurrence_category == 'narrow') ## I am looking for a narrow recurrent bloomer to use as example

bloo_types_summary |>
  dplyr::filter((asv_num %in% c('asv17', 'asv23', 'asv11', 'asv62', 'asv58', 'asv555', 'asv72') & fraction == '3') |
           (asv_num %in% c('asv62', 'asv58', 'asv555') & fraction == '0.2') |
           (asv_num == 'asv11' & fraction == '0.2'))

asv_tab_all_bloo_z_tax_examples <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::mutate(fraction = as.double(fraction)) |>
  dplyr::select(-order,-family, -phylum, -class) |>
  dplyr::left_join(bloo_all_types_summary_tb_tax, by = c( 'asv_num','fraction')) |>
  dplyr::filter(order != 'SAR11 clade') |>
  #dplyr::filter(asv_num_f %in% c('asv17', 'asv23', 'asv11', 'asv84', 'asv62', 'asv58', 'asv555', 'asv559')) |>
  dplyr::filter((asv_num %in% c('asv17', 'asv23', 'asv11', 'asv62', 'asv58', 'asv555', 'asv72') & fraction == 3) |
                  (asv_num %in% c('asv62', 'asv58', 'asv555') & fraction == 0.2) |
                  (asv_num == 'asv11' & fraction == 0.2)) |>
  #dplyr::left_join(occurrence_bloo_bbmo) |>
  dplyr::mutate(bloom = case_when (z_score_ra >= 1.96 &
                                     abundance_value >= 0.1~ 'Bloom',
                                   TRUE ~ 'No-bloom')) |>
  dplyr::mutate(bloom = as.factor(bloom)) |>
  dplyr::mutate(asv_num_f = asv_num)

  asv_summary_bloom <- asv_tab_all_bloo_z_tax_examples |>
  distinct(asv_num, occurrence_category, recurrency, frequency) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num'))
  


# examples_of_blooms <- asv_tab_all_bloo_z_tax_examples |>  
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   #dplyr::filter(asv_num_f %in% c('asv17', 'asv23', 'asv11', 'asv84', 'asv62', 'asv58', 'asv555', 'asv558')) |>
#   ggplot(aes(date, abundance_value))+
#   scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
#   #scale_y_log10( expand = c(0,0))+ #, limits = c(0,0.25)
#   geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'stack')+
#   geom_point(data = asv_tab_all_bloo_z_tax_examples |>
#                dplyr::filter(z_score_ra >= 1.96 &
#                                abundance_value >= 0.1),
#              aes(date, abundance_value, color =  '#9F0011', alpha = 1))+
#   scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
#   # facet_wrap(vars(facet_variable_2), 
#   #            ncol = 2, labeller = labels_type_bloom,
#   #            scales = 'free')+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
#   facet_grid(facet_var_3~facet_var_4, 
#               #labeller = labels_type_bloom,
#              scales = 'free',
#               drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   geom_hline(yintercept = 0.1, linetype = 'dashed')+
#   labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Family')+
#   guides(fill = guide_legend(ncol = 7, size = 8,
#                              override.aes = aes(label = '')),
#          colour = 'none',
#          alpha = 'none')+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(), strip.text = element_text(size = 6, margin = margin(0, 0, 2, 0)),
#         legend.position = 'bottom', axis.text.y = element_text(size = 5),
#         axis.title = element_text(size = 7), #strip.background = element_blank(), 
#         legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
#         plot.margin = margin(2,5,0,5),
#         legend.key.size = unit(3, 'mm'))
# 
# examples_of_blooms
# examples_of_blooms_no_occurrence <- asv_tab_all_bloo_z_tax_examples |>  
#   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
#   #dplyr::filter(asv_num_f %in% c('asv17', 'asv23', 'asv11', 'asv84', 'asv62', 'asv58', 'asv555', 'asv558')) |>
#   ggplot(aes(date, abundance_value))+
#   scale_y_continuous(labels = percent_format(), expand = c(0,0))+ #, limits = c(0,0.25)
#   geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
#   geom_point(data = asv_tab_all_bloo_z_tax_examples |>
#                dplyr::filter(z_score_ra >= 1.96 &
#                                abundance_value >= 0.1),
#              aes(date, abundance_value, color =  '#9F0011', alpha = 1))+
#   scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
#   # facet_wrap(vars(facet_variable_2), 
#   #            ncol = 2, labeller = labels_type_bloom,
#   #            scales = 'free')+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
#   facet_grid(facet_var_5~facet_var_6, 
#              #labeller = labels_type_bloom,
#              scales = 'free',
#              drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
#   scale_fill_manual(values = palette_family_assigned_bloo)+
#   geom_hline(yintercept = 0.1, linetype = 'dashed')+
#   labs(x = 'Date', y = 'Relative abundance (%)', fill = 'Family')+
#   guides(fill = guide_legend(ncol = 7, size = 8,
#                              override.aes = aes(label = '')),
#          colour = 'none',
#          alpha = 'none')+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(), strip.text = element_text(size = 12, margin = margin(4, 4, 4, 4)),
#         legend.position = 'bottom', axis.text.y = element_text(size = 5),
#         axis.title = element_text(size = 7), #strip.background = element_blank(), 
#         legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
#         plot.margin = margin(2,5,0,5),
#         legend.key.size = unit(3, 'mm'))
# 
# examples_of_blooms_no_occurrence

### CONCEPTUAL FIGURE OF THE TYPES OF BLOOMERS THAT WE IDENTIFIED IN THE BBMO ------
### i will prepare each graph individually so that I can organize them in a tidy way easy to follow 

palette_occurrence <- c(narrow = "#AE659B",
                        intermediate = "#3e3e3e",
                        broad = "#57a9a8")

labs_occurrence <- as_labeller(c(narrow = "Narrow\n(<1/3)",
                                 intermediate = "Intermediate\n(1/3 < x < 2/3)",
                                 broad = "Broad\n(>2/3)"))

asv_tab_all_bloo_z_tax_examples |>
  dplyr::filter(asv_num == 'asv62') |>
  distinct(asv_num, fraction)

##plot1 
plot1 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv23')) |>
  ggplot(aes(date, abundance_value))+
  geom_hline(yintercept = 0.1, linetype = 'dashed', linewidth = 0.3)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.3))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),  position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv23'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+ #, title = 'Recurrents'
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 0),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 7), strip.background = element_rect(fill = "#57a9a8",color = "#57a9a8" ), #strip.text.x = element_blank(),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6, title = element_text(size = 8))

##plot2
plot2 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv62')) |>
  ggplot(aes(date, abundance_value))+
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.2))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv62'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 0),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 7), strip.background = element_rect(fill = "#3e3e3e",color = "#3e3e3e"), #strip.text.x = element_blank(),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6, title = element_text(size = 8))
# plot 3
plot3 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv72')) |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0),  limits = c(0,0.3))+ #
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv72'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 0),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 7), strip.background = element_rect(fill = "#AE659B", color = "#AE659B"), #strip.text.x = element_blank(),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6, title = element_text(size = 8))

# plot 4
plot4 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv17')) |>
  ggplot(aes(date, abundance_value))+
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.35))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv17'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+ #, title = 'Non-Recurrents'
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#57a9a8", color = "#57a9a8" ),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6)

# plot 5
plot5 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv84')) |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.2))+ #, limits = c(0,0.25)
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv84'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#3e3e3e",color = "#3e3e3e" ),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6)

# plot 6
plot6 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv559')) |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.15))+ #, limits = c(0,0.25)
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv559'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#AE659B", color = "#AE659B"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6)

# plot 7
plot7 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv11')) |>
dplyr::filter(fraction == '0.2') |>
  ggplot(aes(date, abundance_value))+
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.5))+ #, limits = c(0,0.25)
  geom_area(aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv11',
                             fraction == '0.2'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+ #, title = 'Ephemeral'
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#3e3e3e", color = "#3e3e3e"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6)

# plot 8
plot8 <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num_f %in% c('asv555')) |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.2))+ #, limits = c(0,0.25)
  geom_hline(yintercept = 0.1, linetype = 'dashed', size = 0.3)+
  geom_area( aes(date, y = abundance_value, group = asv_num_f, fill = 'grey'),   position= 'identity')+
  geom_point(data = asv_tab_all_bloo_z_tax_examples |>
               dplyr::filter(bloom == 'Bloom' &
                               asv_num == 'asv555'), aes(color = bloom), size = 0.5) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  facet_wrap(vars(occurrence_category), 
             #labeller = labels_type_bloom,
             scales = 'free',
             drop = T)+ #, labels = labs_fraction  labeller = labs_clusters_pa_fl,
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         colour = 'none',
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(10, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#AE659B", color = "#AE659B"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'), aspect.ratio = 4/6)

# legend
legend_plot <- asv_tab_all_bloo_z_tax_examples |>  
  dplyr::mutate(bloom = case_when (z_score_ra >= 1.96 &
                  abundance_value >= 0.1~ 'Bloom',
                  TRUE ~ 'No-bloom')) |>
  dplyr::mutate(bloom = as.factor(bloom)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, abundance_value)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, 0.2)) +
  geom_area(aes(group = asv_num_f, fill = occurrence_category), position = 'identity') +
  geom_point(aes(color = bloom), size = 1) +  # Specify shape aesthetic for points
  scale_color_manual(values = c( 'Bloom' = '#9F0011', 'No-bloom' = 'white')) +
  scale_fill_manual(values = palette_occurrence, labels = labs_occurrence) +
  labs(x = 'Time (Y)', y = 'Relative abundance (%)', fill = 'Occurrence category', color = 'Bloom')+
  guides(fill = guide_legend(ncol = 3, size = 6),
         color = guide_legend(ncol = 2, size = 6),  # Hide points in the color legend
         alpha = 'none') +
  theme_bw()+
  theme(axis.text.x = element_text(size = 0), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 0, margin = margin(4, 4, 4, 4)),
        legend.position = 'bottom', axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "#AE659B", color = "transparent"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'))

#legend_only  <- cowplot::get_legend(legend_plot)

#library(cowplot)
#library(gridExtra)
# examples_bloomers_types_titles <- grid.arrange(plot1, plot2, plot3,
#              plot4, plot5, plot6,
#              legend_only, plot7, plot8,
#              ncol = 3)

# # Add text to the combined plot
# text1 <- textGrob("Non-recurrent", x = 0.05, y = 1, gp = gpar(fontsize = 8))
# text2 <- textGrob("Non-recurrent", x = 0.4, y = 1, gp = gpar(fontsize = 8))
# text3 <- textGrob("Non-recurrent", x = 0.7, y = 1, gp = gpar(fontsize = 8))
# 
# # Add text to the combined plot
# examples_bloomers_types_with_text <- examples_bloomers_types +
#   annotate("text", x = c(0.05, 0.38, 0.71), y = 1, label = c("Non-recurrent", "Non-recurrent", "Non-recurrent"), 
#            vjust = 1, size = 3.5)

#Save the combined plot
# ggsave(filename = 'examples_bloomers_types.pdf', plot = examples_bloomers_types_titles,
#        path = 'results/figures/',
#        width = 188, height = 160, units = 'mm')

#### I decide to keep only 6 blooming examples ------
legend_only  <- cowplot::get_legend(legend_plot)

#library(cowplot)
#library(gridExtra)

# Define the layout matrix
layout_mat <- rbind(c(1, 2, 3),
                    c(4, 5, 6),
                    c(7, 7, 7))  # This specifies that legend_only occupies all columns in the last row

# Define heights for each row
heights <- c(1, 1, 0.25)  # Adjust the height of the last row to be shorter

# Arrange the plots using the layout matrix and heights
examples_bloomers_types_titles <- grid.arrange(plot1, plot2, plot3,
                                               plot4, plot7, plot8,
                                               legend_only,
                                               layout_matrix = layout_mat,
                                               heights = heights)
# Add text to each row
grid.text("Recurrent", x = unit(0.02, "npc"), y = unit(1, "npc") - unit(0.5 * heights[1], "cm"), just = "left", gp = gpar(fontsize = 8))
grid.text("Non-recurrent", x = unit(0.02, "npc"), y = unit(1, "npc") - unit(heights[1] + 3 * heights[2], "cm"), just = "left", gp = gpar(fontsize = 8))
# Add grid line between rows 1 and 2 (still tpe be improved)
grid.lines(x = unit(c(0, 1), "npc"), y = unit(1, "npc") - unit(heights[1], "npc"), gp = gpar(col = "black", lwd = 1))

#Save the combined plot (trobar la manera de guardar-ho amb el text!!)
ggsave(filename = 'examples_bloomers_types_red.pdf', plot = examples_bloomers_types_titles,
       path = 'results/figures/',
       width = 188, height = 120, units = 'mm')

# ggsave(filename = 'examples_bloomers_types_red.svg', 
#        plot = examples_bloomers_types_titles,
#        path = 'results/figures/poster_svg_format/',
#        width = 188, height = 120, units = 'mm')

### EXAMPLES OF BLOOMING EVENTS ---- 
labs_blooming_events <-  as_labeller(c('0.2' = 'Free living (0.2-3 um)',
                                       '3' = 'Particle attached (3-20 um)',
                                       seasonal = 'Seasonal',
                                       stochastic = 'Stochastic'))

bloom_events_bbmo_description <- asv_tab_all_bloo_z_tax |>
  plyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::left_join(bloo_all_types_summary, by = c('asv_num_f' = 'asv_num','fraction')) |>
  dplyr::filter(!asv_num_f %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(type_of_bloom = str_replace(type_of_bloom, 'ephimeral', 'ephemeral')) |>
  dplyr::mutate(bloom = case_when (z_score_ra >= 1.96 &
                                     abundance_value >= 0.1~ 'Bloom',
                                   TRUE ~ 'No-bloom')) |>
  dplyr::mutate(bloom = as.factor(bloom)) |>
  dplyr::filter(bloom == 'Bloom') |>
  group_by(date, recurrency, frequency, fraction) |>
  dplyr::reframe(bloom_abund = sum(abundance_value)) |>
  ggplot(aes(date, bloom_abund))+
  facet_grid(fraction~frequency, labeller = labs_blooming_events)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  geom_smooth(method = 'loess', span = 0.7, color = 'black')+
  labs(x = 'Time (Y)', y = 'Relative abundance (%)')+ #, fill = 'Occurrence category', color = 'Bloom'
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), #strip.text = element_text(size = 0, margin = margin(4, 4, 4, 4)),
        legend.position = 'right', axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "white", color = "transparent"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'))

# ggsave(filename = 'bloom_events_bbmo_description.pdf', plot = bloom_events_bbmo_description,
#        path = 'results/figures/',
#        width = 188, height = 130, units = 'mm')

bloom_events_bbmo_description <- asv_tab_all_bloo_z_tax |>
  plyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% bloo_3$value &
                  fraction == '3' |
                  asv_num %in% bloo_02$value &
                  fraction == '0.2') |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::left_join(bloo_all_types_summary, by = c('asv_num_f' = 'asv_num','fraction')) |>
  dplyr::filter(!asv_num_f %in% c('asv2', 'asv3', 'asv5', 'asv8')) |>
  dplyr::mutate(type_of_bloom = str_replace(type_of_bloom, 'ephimeral', 'ephemeral')) |>
  dplyr::mutate(bloom = case_when (z_score_ra >= 1.96 &
                                     abundance_value >= 0.1~ 'Bloom',
                                   TRUE ~ 'No-bloom')) |>
  dplyr::mutate(bloom = as.factor(bloom)) |>
  dplyr::filter(bloom == 'Bloom') |>
  group_by(date, recurrency, frequency, fraction) |>
  dplyr::reframe(bloom_abund = sum(abundance_value)) |>
  ggplot(aes(date, bloom_abund))+
  facet_grid(fraction~frequency, labeller = labs_blooming_events)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  geom_smooth(method = 'loess', span = 0.7, color = 'black')+
  labs(x = 'Time (Y)', y = 'rCLR')+ #, fill = 'Occurrence category', color = 'Bloom'
  scale_y_continuous(expand = c(0, 0)) +
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), #strip.text = element_text(size = 0, margin = margin(4, 4, 4, 4)),
        legend.position = 'right', axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 7), strip.background = element_rect(fill = "white", color = "transparent"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 7),# strip.placement = 'outside',
        plot.margin = margin(5,5,5,5),
        legend.key.size = unit(3, 'mm'))

### ASV11 is persistent in FL and ephemeral in PA go a little bit further------
#### It is always ephemeral
asv_tab_all_bloo_z_tax_cl |> 
  dplyr::filter(asv_num_f == 'asv11') |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, abundance_value))+
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,0.6))+ #, limits = c(0,0.25)
  scale_color_manual(values = palette_fraction)+
  scale_linetype_discrete()+
  geom_line(aes(group = fraction, color = fraction, linetype = fraction), linewidth = 0.25)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = '', linetype = '')+
  guides(fill = guide_legend(ncol = 7, size = 8,
                             override.aes = aes(label = '')),
         
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 5, margin = margin(0, 0, 2, 0)),
        legend.position = 'bottom', axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 7), strip.background = element_blank(), 
        legend.text = element_text(size = 5), legend.title = element_text(size = 7), strip.placement = 'outside',
        plot.margin = margin(2,5,0,5),
        legend.key.size = unit(3, 'mm'))

## Interesting ASVs blooming together-----
### asv7 and asv15 cluster----
asv7_asv15_cluster <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv7', 'asv15')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  # group_by(date, fraction) |>
  # dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, abundance_value))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_color_manual(values = palette_family_assigned_bloo)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #geom_line(aes(color = family_f, group = asv_num_f))+
  facet_wrap(vars(fraction), labeller = labs_fraction, scales = 'free')+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(5,15,5,15)) 

asv7_asv15_cluster
# ggsave(filename = 'asv7_asv15_cluster_ed.pdf', plot = asv7_asv15_cluster,
#        path = 'results/figures/',
#        width = 188, height = 80, units = 'mm')

asv7_asv15_cluster <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv7', 'asv15')) |>
  dplyr::filter(abundance_type == 'rclr') |>
  # group_by(date, fraction) |>
  # dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, abundance_value))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous( expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 1,  position='stack')+
  scale_color_manual(values = palette_family_assigned_bloo)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #geom_line(aes(color = family_f, group = asv_num_f))+
  facet_wrap(vars(fraction), labeller = labs_fraction)+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(5,15,5,15)) 

asv7_asv15_cluster
# ggsave(filename = 'asv7_asv15_cluster_rclr_v1.pdf', plot = asv7_asv15_cluster,
#        path = 'results/figures/',
#        width = 188, height = 80, units = 'mm')

ggsave(filename = 'asv7_asv15_cluster_ed_rclr.svg', plot = asv7_asv15_cluster,
       path = 'results/figures/poster_svg_format/',
       width = 188, height = 80, units = 'mm')


##cluster seasonal 2 in the PA ----
c('asv100', 'asv25', 'asv72', 'asv42')

seasonal2_cluster_pa <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv100', 'asv25', 'asv72', 'asv42')) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  # group_by(date, fraction) |>
  # dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, abundance_value))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous(labels = percent_format(), expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 0.8,  position='stack')+
  scale_color_manual(values = palette_family_assigned_bloo)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #geom_line(aes(color = family_f, group = asv_num_f))+
  facet_wrap(vars(fraction), labeller = labs_fraction)+
  labs(x = 'Time', y = 'Relative abundance (%)', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(5,15,5,15)) 

seasonal2_cluster_pa
# 
# ggsave(filename = 'seasonal2_cluster_pa_ed.pdf', plot = seasonal2_cluster_pa,
#        path = 'results/figures/',
#        width = 188, height = 80, units = 'mm')

seasonal2_cluster_pa <- asv_tab_all_bloo_z_tax |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::filter(asv_num %in% c('asv100', 'asv25', 'asv72', 'asv42')) |>
  dplyr::filter(abundance_type == 'rclr') |>
  # group_by(date, fraction) |>
  # dplyr::mutate(max_abund = sum(abundance_value)) |>
  ungroup() |>
  ggplot(aes(date, abundance_value))+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0)
                   #limits = c(min(asv_tab_all_bloo_z_tax$date), max(asv_tab_all_bloo_z_tax$date),
                   #limits = c(as.POSIXct(2004-01-26, origin = '2004-01-26'), as.POSIXct(2014-01-01, origin = '2014-01-01'))
  )+
  # geom_rect(data = harbour_restoration, mapping=aes(xmin = date_min, xmax = date_max, x=NULL, y=NULL,
  # ymin = -Inf, ymax = Inf), fill = '#C7C7C7', alpha = 0.5)+
  scale_y_continuous( expand = c(0,0))+
  geom_area(aes(date, abundance_value, fill = family_f, group = asv_num_f), alpha = 1,  position='stack')+
  scale_color_manual(values = palette_family_assigned_bloo)+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  #geom_line(aes(color = family_f, group = asv_num_f))+
  facet_wrap(vars(fraction), labeller = labs_fraction)+
  labs(x = 'Time', y = 'rCLR', fill = 'Family')+
  guides(fill = guide_legend(ncol = 6, size = 10,
                             override.aes = aes(label = '')),
         alpha = 'none')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 7),
        legend.position = 'bottom', axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8), strip.background = element_blank(), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 8), strip.placement = 'outside',
        plot.margin = margin(5,15,5,15)) 

seasonal2_cluster_pa
# 
# ggsave(filename = 'seasonal2_cluster_pa_ed_rclr.pdf', plot = seasonal2_cluster_pa,
#        path = 'results/figures/',
#        width = 188, height = 80, units = 'mm')

#### SUMMARY OF TYPES OF BLOOMS AND THE CLUSTER THEY BELONG TO FOR PA AND FL----
  
  ## relative abundance category is related to the rarefied dataset < 1% criteria from Alonso Saez Alonso-Sáez L, Díaz-Pérez L, Morán XAG. The hidden seasonality of the rare biosphere in coastal marine bacterioplankton. Environ Microbiol. 2015;17:3766–80
  asv_tab_10y_3_rel_rar |>
    colnames()
  
  # asv_tab_10y_3_rel_rar |>
  #   dplyr::filter(relative_abundance < 0.01) |>
  #   group_by(asv_num) |>
  #   dplyr::reframe(n = n()) |>
  #   dplyr::filter(n == 117) |>
  #   dplyr::filter(asv_num %in% bloo_3$value)
library(ggridges)
  
  asv_tab_all_bloo_z_tax |>
    dplyr::filter(abundance_type == 'relative_abundance') |>
    group_by(asv_num, fraction) |>
    dplyr::reframe(mean_abund = mean(abundance_value)) |>
    ggplot(aes(x = mean_abund, fill = fraction))+
    scale_y_discrete(labels = labs_fraction)+
    geom_density_ridges2(aes(y = fraction), alpha = 0.8, scale = 1)+
    scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
    labs(fill = 'Fraction')+
    theme_ridges(font_size = 6, center_axis_labels = F, grid = F)+
    theme(axis.text.y = element_blank())
  
 asv_tab_all_bloo_z_tax_summary_02 <- asv_tab_all_bloo_z_tax |>
    dplyr::filter(asv_num %in% bloo_02$value) |>
    dplyr::filter(abundance_type == 'relative_abundance' &
                    fraction == '0.2') |>
    group_by(asv_num) |>
    dplyr::mutate(relative_abundance_category = case_when(mean(abundance_value) > 0.01 ~ 'abundant',
                                                          mean(abundance_value) <= 0.01 ~ 'rare')) |>
   ungroup() |>
    dplyr::left_join(bloo_02_types_summary)
 
 asv_tab_all_bloo_z_tax_summary_3 <- asv_tab_all_bloo_z_tax |>
    dplyr::filter(asv_num %in% bloo_3$value) |>
    dplyr::filter(abundance_type == 'relative_abundance' &
                    fraction == '3') |>
    group_by(asv_num) |>
    dplyr::mutate(relative_abundance_category = case_when(mean(abundance_value) > 0.01 ~ 'abundant',
                                                          mean(abundance_value) <= 0.01 ~ 'rare')) |>
   ungroup() |>
    dplyr::left_join(bloo_3_types_summary)
 
 asv_tab_all_bloo_z_tax_summary_3 <-  asv_tab_all_bloo_z_tax_summary_3 |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))

 asv_tab_all_bloo_z_tax_summary_02  <- asv_tab_all_bloo_z_tax_summary_02  |>
   dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d")))
   
asv_tab_all_bloo_z_tax_summary_all <- asv_tab_all_bloo_z_tax_summary_02 |>
   bind_rows(asv_tab_all_bloo_z_tax_summary_3)
  
 # write_csv(asv_tab_all_bloo_z_tax_summary_all, 'data/asv_tab_all_bloo_z_tax_summary.csv')
 
 #### ANSWERING DIFFERENT QUESTIONS DEPENDING ON THE TYPE OF BLOOM THAT THEY ARE-----
bloo_02_types_summary <- bloo_02_types_summary |>
   dplyr::mutate(fraction = '0.2')
 
bloo_3_types_summary <- bloo_3_types_summary |>
   dplyr::mutate(fraction = '3')
 
bloo_all_types_summary <- bloo_02_types_summary |>
  bind_rows(bloo_3_types_summary)
 
bloo_all_types_summary |>
  group_by(recurrency, fraction) |>
  dplyr::reframe(n = n())

bloo_all_types_summary |>
  group_by(recurrency, frequency, fraction) |>
  dplyr::reframe(n = n())

bloo_all_types_summary |>
  group_by(recurrency, type_of_bloom, fraction) |>
  dplyr::reframe(n = n())

bloo_all_types_summary |>
  group_by(clustering_group, fraction) |>
  dplyr::reframe(n = n())

bloo_all_types_summary |>
  group_by(occurrence_category, fraction) |>
  dplyr::reframe(n = n())

## check if bloomers have the same behavior in both size fractions-----

### look at Marta's script differential_abundance_corncob_forONA.Rmd to do it more quantiatively!
shared_blooms <- bloo_all_types_summary |>
  dplyr::filter(duplicated(asv_num)) 

bloo_all_types_summary |>
  dplyr::filter(asv_num %in% shared_blooms$asv_num) |>
  arrange(asv_num)

##for those that are always considered bloomers
shared_bloo_fractions <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::filter(asv_num %in% shared_blooms$asv_num) |>
  dplyr::mutate(asv_fam = paste0(family_f, ' ', asv_num)) |>
  ggplot(aes(date, abundance_value))+
  geom_line(aes(color = fraction, linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  facet_wrap(vars(asv_fam), ncol = 3)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  scale_y_continuous(label = percent_format(), expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = 'Fraction', linetype = 'Fraction')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 6),
        legend.position = 'bottom', axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6), strip.background = element_blank(),
        legend.text = element_text(size = 5), strip.placement = 'outside',
        legend.title = element_text(size = 5),
        aspect.ratio = 3/7)

# ggsave(filename = 'shared_bloo_fractions.pdf', plot = shared_bloo_fractions,
#        path = 'results/figures/relationship_bloo_02_3/',
#        width = 188, height = 100, units = 'mm')

##for those that are only bloomers in one size fraction
m_bbmo_10y_1 <- m_bbmo_10y |>
  dplyr::select(date, sample_id)

exclusive_02 <- bloo_02 |>
  anti_join(bloo_3)

exclusive_3 <- bloo_3 |>
  anti_join(bloo_02)

exclusive_02_in_3 <- asv_tab_10y_3_rel |>
  dplyr::mutate(fraction = '3') |>
  dplyr::filter(asv_num %in% exclusive_02$value)

exclusive_02_in_3_plot <- asv_tab_10y_02_rel |>
  dplyr::mutate(fraction = '0.2') |>
  dplyr::filter(asv_num %in% exclusive_02$value) |>
  bind_rows(exclusive_02_in_3) |>
  left_join(m_bbmo_10y_1) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  dplyr::mutate(asv_fam = paste0(family_f, ' ', asv_num)) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  ggplot(aes(date, relative_abundance))+
  geom_line(aes(color = fraction, linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  facet_wrap(vars(asv_fam), ncol = 3)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  scale_y_continuous(label = percent_format(), expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = 'Fraction', linetype = 'Fraction')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 6),
        legend.position = 'bottom', axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6), strip.background = element_blank(),
        legend.text = element_text(size = 5), strip.placement = 'outside',
        legend.title = element_text(size = 5),
        aspect.ratio = 3/7)

# ggsave(filename = 'exclusive_02_in_3_plot.pdf', plot = exclusive_02_in_3_plot,
#        path = 'results/figures/relationship_bloo_02_3/',
#        width = 188, height = 150, units = 'mm')

exclusive_3_in_02 <- asv_tab_10y_02_rel |>
  dplyr::mutate(fraction = '0.2') |>
  dplyr::filter(asv_num %in% exclusive_3$value)

asv_num_gr1 <- exclusive_3_in_02  |>
  ungroup() |>
  distinct(asv_num) |>
  slice_max(n = 15, order_by = asv_num)

asv_num_gr2 <- exclusive_3_in_02  |>
  dplyr::filter(!asv_num %in% asv_num_gr1$asv_num) |>
  ungroup() |>
  distinct(asv_num) |>
  slice_max(n = 15, order_by = asv_num)

asv_num_gr3 <- exclusive_3_in_02  |>
  dplyr::filter(!asv_num %in% asv_num_gr1$asv_num &
                  !asv_num %in% asv_num_gr2$asv_num ) |>
  ungroup() |>
  distinct(asv_num) 

exclusive_3_in_02_plot_1 <- asv_tab_10y_3_rel |>
  dplyr::mutate(fraction = '3') |>
  dplyr::filter(asv_num %in% exclusive_3$value) |>
  bind_rows(exclusive_3_in_02) |>
  dplyr::filter(asv_num %in% asv_num_gr1$asv_num) |>
  left_join(m_bbmo_10y_1) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(asv_fam = paste0(family_f, ' ', asv_num)) |>
  ggplot(aes(date, relative_abundance))+
  geom_line(aes(color = fraction, linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  facet_wrap(vars(asv_fam), ncol = 3)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  scale_y_continuous(label = percent_format(), expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = 'Fraction', linetype = 'Fraction')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 6),
        legend.position = 'bottom', axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6), strip.background = element_blank(),
        legend.text = element_text(size = 5), strip.placement = 'outside',
        legend.title = element_text(size = 5),
        aspect.ratio = 3/7)

exclusive_3_in_02_plot_2 <- asv_tab_10y_3_rel |>
  dplyr::mutate(fraction = '3') |>
  dplyr::filter(asv_num %in% exclusive_3$value) |>
  bind_rows(exclusive_3_in_02) |>
  dplyr::filter(asv_num %in% asv_num_gr2$asv_num) |>
  left_join(m_bbmo_10y_1) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(asv_fam = paste0(family_f, ' ', asv_num)) |>
  ggplot(aes(date, relative_abundance))+
  geom_line(aes(color = fraction, linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  facet_wrap(vars(asv_fam), ncol = 3)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  scale_y_continuous(label = percent_format(), expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = 'Fraction', linetype = 'Fraction')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 6),
        legend.position = 'bottom', axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6), strip.background = element_blank(),
        legend.text = element_text(size = 5), strip.placement = 'outside',
        legend.title = element_text(size = 5),
        aspect.ratio = 3/7)

exclusive_3_in_02_plot_3 <- asv_tab_10y_3_rel |>
  dplyr::mutate(fraction = '3') |>
  dplyr::filter(asv_num %in% exclusive_3$value) |>
  bind_rows(exclusive_3_in_02) |>
  dplyr::filter(asv_num %in% asv_num_gr3$asv_num) |>
  left_join(m_bbmo_10y_1) |>
  left_join(bloo_taxonomy, by = c('asv_num' = 'asv_num_f')) |>
  dplyr::mutate(date = (as.POSIXct(date, format = "%Y-%m-%d"))) |>
  dplyr::mutate(asv_fam = paste0(family_f, ' ', asv_num)) |>
  ggplot(aes(date, relative_abundance))+
  geom_line(aes(color = fraction, linetype = fraction))+
  scale_linetype_discrete(labels = labs_fraction)+
  facet_wrap(vars(asv_fam), ncol = 3)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0))+
  scale_y_continuous(label = percent_format(), expand = c(0,0))+
  labs(x = 'Date', y = 'Relative abundance (%)', color = 'Fraction', linetype = 'Fraction')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 4), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 6),
        legend.position = 'bottom', axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6), strip.background = element_blank(),
        legend.text = element_text(size = 5), strip.placement = 'outside',
        legend.title = element_text(size = 5),
        aspect.ratio = 3/7)
# 
# ggsave(filename = 'exclusive_3_in_02_plot1.pdf', plot = exclusive_3_in_02_plot_1,
#        path = 'results/figures/relationship_bloo_02_3/',
#        width = 188, height = 150, units = 'mm')
# 
# ggsave(filename = 'exclusive_3_in_02_plot2.pdf', plot = exclusive_3_in_02_plot_2,
#        path = 'results/figures/relationship_bloo_02_3/',
#        width = 188, height = 150, units = 'mm')
# 
# ggsave(filename = 'exclusive_3_in_02_plot3.pdf', plot = exclusive_3_in_02_plot_3,
#        path = 'results/figures/relationship_bloo_02_3/',
#        width = 188, height = 150, units = 'mm')
##try to plot a heatmap with labels 

## strong signal on s4 understand why -----

### Compare two dendograms----

# Custom these kendo, and place them in a list
# dl <- dendlist(
#   hc1 %>% 
#     set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
#     set("branches_lty", 1) %>%
#     set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
#   hc3 %>% 
#     set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
#     set("branches_lty", 1) %>%
#     set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3)
# )
# 
# dl <- dendlist(
#   hc1 %>% 
#     set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
#     set("branches_lty", 1) %>%
#     set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
#   hc5 %>% 
#     set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
#     set("branches_lty", 1) %>%
#     set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3)
# )
# 
# dl <- dendlist(
#   hc3 %>% 
#     set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
#     set("branches_lty", 1) %>%
#     set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
#   hc5 %>% 
#     set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
#     set("branches_lty", 1) %>%
#     set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3)
# )
# 
# # Plot them together
# tanglegram(dl, 
#            common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
#            margin_inner=7,
#            lwd=2
# )

## dendogram + heatmap
library(pheatmap) ## for heatmap generation
##tibble with information
anot <- time_series_3 |>
  colnames() |>
  as_tibble_col(column_name = 'asv_num') |>
  dplyr::filter(asv_num != c('decimal_date'))

anot <- time_series_3 |>
  dplyr::select(decimal_date)

distances |>
  dim()

pheatmap(distances_1,scale="row", #annotation_col = anot$decimal_date,
         color=colorRampPalette(c("navy", "white", "red"))(25))

pheatmap(distances_2,scale="row", #annotation_col = anot$decimal_date,
         color=colorRampPalette(c("navy", "white", "red"))(25))

pheatmap(distances_3,scale="row", #annotation_col = dfh,
         color=colorRampPalette(c("navy", "white", "red"))(25))

pheatmap(distances_5,scale="row", #annotation_col = dfh,
         color=colorRampPalette(c("navy", "white", "red"))(25))

## Cross-correlation for the PA fraction----
## data preparation rows are observations and columns are variables----
time_series_1 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd1' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_2 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd2' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_3 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd3' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_4 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_5 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 's4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

# Dissimilarity matrix
distances_1 <- time_series_1 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_2 <- time_series_2 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_3 <- time_series_3 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_5 <- time_series_5 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

# Hierarchical clustering using ward.d
distances_1 |>
  str()

hc1 <- hclust(distances_1, method = "ward.D" )  |>
  as.dendrogram()

# hc2 <- hclust(distances_2, method = "ward.D" ) |>
#   as.dendrogram()
# 
# hc3 <- hclust(distances_3, method = "ward.D" ) |>
#   as.dendrogram()
# 
# hc5 <- hclust(distances_5, method = "ward.D" ) |>
#   as.dendrogram()
# 
# # Plot the obtained dendrogram
# plot(hc1, cex = 0.6, hang = -1)
# plot(hc2, cex = 0.6, hang = -1)
# 
# # Compute with agnes
# hc2 <- agnes(distances, method = 'ward.D')
# 
# # Agglomerative coefficient
# hc2$ac
# 
# # methods to assess
# m <- c( "average", "single", "complete", "ward")
# names(m) <- c( "average", "single", "complete", "ward")
# 
# # function to compute coefficient
# ac <- function(x) {
#   agnes(distances, method = x)$ac
# }
# 
# map_dbl(m, ac)
# 
# hc3 <- agnes(distances, method = "ward")
# pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes")

### Cross-correlation prepare the dataset------
distances_1 <- distances_1 |>
  as.matrix()

distances_2 <- distances_2 |>
  as.matrix()

distances_3 <- distances_3 |>
  as.matrix()

distances_5 <- distances_5 |>
  as.matrix()

# List of distances matrices
distance_matrices <- list(distances_1, distances_2, distances_3, distances_5)

# Initialize a list to store cross-correlation matrices and hierarchical clustering objects
cross_correlation_matrices <- list()
hc_list <- list()

# Perform cross-correlations and hierarchical clustering
for (i in 1:(length(distance_matrices) - 1)) {
  for (j in (i + 1):length(distance_matrices)) {
    name <- paste("cross_correlation_", i, ".", j, sep = "")
    cross_correlation_matrices[[name]] <- cor(distance_matrices[[i]], distance_matrices[[j]])
    hc_list[[name]] <- hclust(as.dist(1 - cross_correlation_matrices[[name]]), method = 'ward.D')
  }
}

# Plot heatmaps with hierarchical clustering
for (i in 1:(length(distance_matrices) - 1)) {
  for (j in (i + 1):length(distance_matrices)) {
    name <- paste("cross_correlation_", i, ".", j, sep = "")
    heatmap(cross_correlation_matrices[[name]], Rowv = as.dendrogram(hc_list[[name]]), Colv = as.dendrogram(hc_list[[name]]),
            col = colorRampPalette(c("navy", "white", "red"))(25),
            scale = "none",
            main = paste("Cross-correlation (3)", name),
            xlab = '', ylab = '')
  }
}

# Add color legend
image.plot(x = c(-1, 1), y = c(-1, 1), z = matrix(seq(-1, 1, length.out = 100), nrow = 1),
           col = colorRampPalette(c("navy", "white", "red"))(100),
           zlim = c(-1, 1),
           legend.only = TRUE, legend.shrink = 0.6, legend.width = 0.5,
           legend.mar = 5, legend.lab = "Correlation values")

### LOOP TO SAVE THEM
# Define output directory
output_dir <- "results/figures/heatmaps_cross_corr_3_biased_red/"

# Create the directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Plot heatmaps and save them
for (i in 1:(length(distance_matrices) - 1)) {
  for (j in (i + 1):length(distance_matrices)) {
    name <- paste("cross_correlation_", i, ".", j, sep = "")
    heatmap_file <- paste(output_dir, "heatmap_", name, ".png", sep = "")
    
    heatmap(cross_correlation_matrices[[name]], Rowv = as.dendrogram(hc_list[[name]]), Colv = as.dendrogram(hc_list[[name]]),
            col = colorRampPalette(c("navy", "white", "red"))(25),
            scale = "none",
            main = paste("Cross-correlation (3)", name),
            xlab = '', ylab = '')
    
    # Save the heatmap as a PNG file
    dev.print(png, file = heatmap_file, width = 800, height = 600)
    dev.off()  # Close the PNG device
  }
}

# cross_correlation_1.2 <- cor(distances_1, distances_2)
# 
# cross_correlation_1.3 <- cor(distances_1, distances_3)
# 
# cross_correlation_1.5 <- cor(distances_1, distances_5)
# 
# cross_correlation_2.5 <- cor(distances_2, distances_5)
# 
# cross_correlation_3.5 <- cor(distances_3, distances_5)
# 
# cross_correlation_2.3 <- cor(distances_2, distances_3)
# 
# distances |>
#   str()
# 
# dim(cross_correlation)
# 
# hc <- hclust(as.dist(1 - cross_correlation_1.2), method = 'ward.D')
# hc <- hclust(as.dist(1 - cross_correlation), method = 'ward.D')
# 
# # Plot the heatmap
# #library(ggplot2)
# # heatmap(cross_correlation, Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc),
# #         col = colorRampPalette(c("navy", "white", "red"))(45))
# #   
# # legend("bottomright", legend = seq(min(cross_correlation), max(cross_correlation), length.out = 5), fill = colorRampPalette(c("navy", "white", "red"))(25), 
# #        title = "Cross-correlation\nvalues")
# 
# # heatmap(cross_correlation, 
# #         Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc),
# #         xlab = "", ylab = "", 
# #         main = "",
# #         scale = "column",
# #         #margins = c(-1,100,40,20),
# #         grid_color = "white",
# #         grid_width = 0.00001,
# #         titleX = FALSE,
# #         hide_colorbar = TRUE,
# #         #branches_lwd = 0.1,
# #         label_names = c("Country", "Feature:", "Value"),
# #         #fontsize_row = 4, fontsize_col = 4,
# #        # labCol = colnames(),
# #         #labRow = rownames(mat),
# #         heatmap_layers = theme(axis.line=element_blank()))
# 
# heatmap(cross_correlation, Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc),
#         col = colorRampPalette(c("navy", "white", "red"))(25),
#         scale = "none", # prevent scaling of the data
#         #main = "Cross-correlation Heatmap",
#         xlab = 'Fine-scale', ylab = 'Half-yearly')
# 
# # Add color legend using image.plot()
# 
# image.plot(x = c(-1, 1), y = c(-1, 1), z = matrix(seq(-1, 1, length.out = 100), nrow = 1),
#            col = colorRampPalette(c("navy", "white", "red"))(100),
#            zlim = c(-1, 1),
#            legend.only = TRUE, legend.shrink = 0.8, legend.width = 1,
#            legend.mar = 10, legend.lab = "Correlation values")
# 
# heatmap(cross_correlation, Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc),
#         col = colorRampPalette(c("navy", "white", "red"))(25),
#         scale = "none", # prevent scaling of the data
#         main = "Cross-correlation Heatmap",
#         xlab = 'Fine-scale', ylab = 'Seasonal')

### cross correlation for the FL fraction ------
time_series_1 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd1' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_2 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd2' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_3 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd3' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_4 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 'd4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

time_series_5 <- wavelets_result_tibble_tax_02_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  left_join(tax_bbmo_10y_new, by = 'asv_num') |>
  dplyr::mutate(asv_f = paste0(family,'.',asv_num)) |>
  dplyr::select(-c(family, asv_num, seq, class, order, family, domain, phylum, genus)) |>
  dplyr::filter(wavelets_transformation == 's4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_f, values_from = wavelets_result_ed)

# Dissimilarity matrix
distances_1 <- time_series_1 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_2 <- time_series_2 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_3 <- time_series_3 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_4 <- time_series_4 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_5 <- time_series_5 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

# # Convert cross-correlation matrix to data frame
# cross_correlation_df <- melt(as.matrix(cross_correlation))
# 
# # Plot heatmap using ggplot with color legend
# ggplot(cross_correlation_df, aes(Var1, Var2, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "navy", mid = "white", high = "red", 
#                        midpoint = 0, limits = c(-1, 1), name = "Correlation values") +
#   labs(x = "Variable 1", y = "Variable 2", title = "Cross-correlation Heatmap")

### Cross-correlation------
distances_1 <- distances_1 |>
  as.matrix()

distances_2 <- distances_2 |>
  as.matrix()

distances_3 <- distances_3 |>
  as.matrix()

distances_4 <- distances_4 |>
  as.matrix()

distances_5 <- distances_5 |>
  as.matrix()

# cross_correlation_1.2 <- cor(distances_1, distances_2)
# 
# cross_correlation_1.3 <- cor(distances_1, distances_3)
# 
# cross_correlation_1.5 <- cor(distances_1, distances_5)
# 
# cross_correlation_2.5 <- cor(distances_2, distances_5)
# 
# cross_correlation_3.5 <- cor(distances_3, distances_5)
# 
# cross_correlation_2.3 <- cor(distances_2, distances_3)
# 
# distances |>
#   str()
# 
# dim(cross_correlation)
# 
# hc <- hclust(as.dist(1 - cross_correlation_1.2), method = 'ward.D')
# hc <- hclust(as.dist(1 - cross_correlation_1.3), method = 'ward.D')
# hc <- hclust(as.dist(1 - cross_correlation), method = 'ward.D')
# 
# # Plot the heatmap
# #library(ggplot2)
# # heatmap(cross_correlation, Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc),
# #         col = colorRampPalette(c("navy", "white", "red"))(45))
# #   
# # legend("bottomright", legend = seq(min(cross_correlation), max(cross_correlation), length.out = 5), fill = colorRampPalette(c("navy", "white", "red"))(25), 
# #        title = "Cross-correlation\nvalues")
# 
# # heatmap(cross_correlation, 
# #         Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc),
# #         xlab = "", ylab = "", 
# #         main = "",
# #         scale = "column",
# #         #margins = c(-1,100,40,20),
# #         grid_color = "white",
# #         grid_width = 0.00001,
# #         titleX = FALSE,
# #         hide_colorbar = TRUE,
# #         #branches_lwd = 0.1,
# #         label_names = c("Country", "Feature:", "Value"),
# #         #fontsize_row = 4, fontsize_col = 4,
# #        # labCol = colnames(),
# #         #labRow = rownames(mat),
# #         heatmap_layers = theme(axis.line=element_blank()))
# 
# heatmap(cross_correlation_1.3, Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc),
#         col = colorRampPalette(c("navy", "white", "red"))(25),
#         scale = "none", # prevent scaling of the data
#         #main = "Cross-correlation Heatmap",
#         xlab = 'Fine-scale', ylab = 'Seasonal')

# # Add color legend using image.plot()
# 
# image.plot(x = c(-1, 1), y = c(-1, 1), z = matrix(seq(-1, 1, length.out = 100), nrow = 1),
#            col = colorRampPalette(c("navy", "white", "red"))(100),
#            zlim = c(-1, 1),
#            legend.only = TRUE, legend.shrink = 0.8, legend.width = 0.5,
#            legend.mar = 10, legend.lab = "Correlation values")


## CREAE A LOOP TO PERFORM THE CROSS-CORRELATIONS -----
#### now it sais incompatible dimensions, this could be because when I remove the samples most affected by margins effects, depending on the decomposition 
#### we remove more or less samples.

# List of distances matrices
distance_matrices <- list(distances_1, distances_2, distances_3, distances_5)

# Initialize a list to store cross-correlation matrices and hierarchical clustering objects
cross_correlation_matrices <- list()
hc_list <- list()

# Perform cross-correlations and hierarchical clustering
for (i in 1:(length(distance_matrices) - 1)) {
  for (j in (i + 1):length(distance_matrices)) {
    name <- paste("cross_correlation_", i, ".", j, sep = "")
    cross_correlation_matrices[[name]] <- cor(distance_matrices[[i]], distance_matrices[[j]])
    hc_list[[name]] <- hclust(as.dist(1 - cross_correlation_matrices[[name]]), method = 'ward.D')
  }
}

# Plot heatmaps with hierarchical clustering
for (i in 1:(length(distance_matrices) - 1)) {
  for (j in (i + 1):length(distance_matrices)) {
    name <- paste("cross_correlation_", i, ".", j, sep = "")
    heatmap(cross_correlation_matrices[[name]], Rowv = as.dendrogram(hc_list[[name]]), Colv = as.dendrogram(hc_list[[name]]),
            col = colorRampPalette(c("navy", "white", "red"))(25),
            scale = "none",
            main = paste("Cross-correlation (0.2)", name),
            xlab = '', ylab = '')
  }
}

# Add color legend
image.plot(x = c(-1, 1), y = c(-1, 1), z = matrix(seq(-1, 1, length.out = 100), nrow = 1),
           col = colorRampPalette(c("navy", "white", "red"))(100),
           zlim = c(-1, 1),
           legend.only = TRUE, legend.shrink = 0.6, legend.width = 0.5,
           legend.mar = 5, legend.lab = "Correlation values")

### LOOP TO SAVE THEM
# Define output directory
output_dir <- "results/figures/heatmaps_cross_corr_02_biased_red/"

# Create the directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Plot heatmaps and save them
for (i in 1:(length(distance_matrices) - 1)) {
  for (j in (i + 1):length(distance_matrices)) {
    name <- paste("cross_correlation_", i, ".", j, sep = "")
    heatmap_file <- paste(output_dir, "heatmap_", name, ".png", sep = "")
    
    heatmap(cross_correlation_matrices[[name]], Rowv = as.dendrogram(hc_list[[name]]), Colv = as.dendrogram(hc_list[[name]]),
            col = colorRampPalette(c("navy", "white", "red"))(25),
            scale = "none",
            main = paste("Cross-correlation (0.2)", name),
            xlab = '', ylab = '')
    
    # Save the heatmap as a PNG file
    dev.print(png, file = heatmap_file, width = 800, height = 600)
    dev.off()  # Close the PNG device
  }
}

# Convert cross-correlation matrix to data frame-----
library(reshape2)
library(dendextend)
cross_correlation_df <- melt(as.matrix(cross_correlation))

# Plot heatmap using ggplot with dendrogram
heatmap_plot <- ggplot(cross_correlation_df, aes(Var1, Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradient(low = "blue", high = "red", name = "Cross-correlation values") +
  labs(x = "Variable 1", y = "Variable 2", title = "Cross-correlation Heatmap")

# Add dendrogram to the plot
heatmap_with_dendrogram <- heatmap_plot + 
  geom_hline(yintercept = hc$hclust$height[hc$hclust$order], color = "gray") +
  geom_vline(xintercept = hc$hclust$height[hc$hclust$order], color = "gray")

# Print the heatmap with dendrogram
print(heatmap_with_dendrogram)

# In the dendrogram displayed above, each leaf corresponds to one observation. As we move up the tree, observations that are similar to each other are combined into branches, which are themselves fused at a higher height.
# The height of the fusion, provided on the vertical axis, indicates the (dis)similarity between two observations. The higher the height of the fusion, the less similar the observations are. Note that, conclusions about the proximity of two observations can be drawn only based on the height where branches containing those two observations first are fused. We cannot use the proximity of two observations along the horizontal axis as a criteria of their similarity.
# The height of the cut to the dendrogram controls the number of clusters obtained.

# Ward's method
hc5 <- hclust(distances, method = "ward.D2" )

# Cut tree into 4 groups
sub_grp <- cutree(hc5, k = 4)

# Number of members in each cluster
table(sub_grp)

USArrests %>%
  dplyr::mutate(cluster = sub_grp) %>%
  head

plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 4, border = 2:5)

fviz_cluster(list(data = distances, cluster = sub_grp))

# Compute euclidean distance and ward hierarchical clustering-----
dist_matrix <- dist(rbind(ts1, ts2)) 

time_series_1 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  dplyr::filter(wavelets_transformation == 'd1' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_num, values_from = wavelets_result)

time_series_2 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed, asv_num, wavelets_transformation) |>
  dplyr::filter(wavelets_transformation == 'd3' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_num, values_from = wavelets_result)

# Get unique taxa
unique_taxa <- unique(wavelets_result_tibble_tax_3_biased$asv_num)

# Initialize lists to store cross-correlation results
cross_corr_results_ts1 <- list()
cross_corr_results_ts2 <- list()

# Iterate over each taxon
for (taxon in unique_taxa) {
  # Subset data for the current taxon
  subset_data <- time_series_1 |>
    dplyr::filter(asv_num == taxon)  # Remove the 'taxa' column# Remove the 'taxa' column
  
  # Compute cross-correlation for time series 1 and current taxon
  cross_corr_ts1 <- ccf(time_series_1, subset_data)
  cross_corr_results_ts1[[taxon]] <- cross_corr_ts1
  
  # Compute cross-correlation for time series 2 and current taxon
  cross_corr_ts2 <- ccf(time_series_2, subset_data)
  cross_corr_results_ts2[[taxon]] <- cross_corr_ts2
}



##### CLUSTERING ENVIRONMENTAL VARIABLES AND CROSSCOREATE THEM WITH ASVs (USING WAVELETS RESULTS)------
wavelets_result_env_tibble_red <- wavelets_result_env_tibble_red |>
  dplyr::mutate(environmental_variable = str_replace(environmental_variable,'_no_nas',''))

### Prepare the environmental data----
time_series_env_1 <- wavelets_result_env_tibble_red |>
  dplyr::select(decimal_date, wavelets_result_ed, environmental_variable, wavelets_transformation) |>
  dplyr::filter(wavelets_transformation == 'd1' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = environmental_variable, values_from = wavelets_result_ed)

time_series_env_2 <-wavelets_result_env_tibble_red |>
  dplyr::select(decimal_date, wavelets_result_ed, environmental_variable, wavelets_transformation) |>
  dplyr::filter(wavelets_transformation == 'd2' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = environmental_variable, values_from = wavelets_result_ed)

time_series_env_3 <- wavelets_result_env_tibble_red |>
  dplyr::select(decimal_date, wavelets_result_ed, environmental_variable, wavelets_transformation) |>
  dplyr::filter(wavelets_transformation == 'd3' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = environmental_variable, values_from = wavelets_result_ed)

time_series_env_4 <-wavelets_result_env_tibble_red |>
  dplyr::select(decimal_date, wavelets_result_ed, environmental_variable, wavelets_transformation) |>
  dplyr::filter(wavelets_transformation == 'd4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = environmental_variable, values_from = wavelets_result_ed)

time_series_env_5 <- wavelets_result_env_tibble_red |>
  dplyr::select(decimal_date, wavelets_result_ed, environmental_variable, wavelets_transformation) |>
  dplyr::filter(wavelets_transformation == 's4' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = environmental_variable, values_from = wavelets_result_ed)

# time_series_env_1 |>
#   pivot_longer(cols = c(-decimal_date)) |>
#   ggplot(aes(decimal_date, name))+
#   geom_tile(aes(fill = value))+
#   #scale_fill_gradientn(colours = scale_fill_viridis_b)+
#   theme_minimal()

# Dissimilarity matrix
distances_env_1 <- time_series_env_1 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_env_2 <- time_series_env_2 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_env_3 <- time_series_env_3 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_env_4 <- time_series_env_4 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

distances_env_5 <- time_series_env_5 |>
  dplyr::select(-decimal_date) |>
  t() |>
  stats::dist( method = "euclidean")

# Hierarchical clustering
hc1_env <- hclust(distances_env_1, method = "ward.D" )  |>
  as.dendrogram()

hc2_env <- hclust(distances_env_2, method = "ward.D" )  |>
  as.dendrogram()

hc3_env <- hclust(distances_env_3, method = "ward.D" )  |>
  as.dendrogram()

hc4_env <- hclust(distances_env_4, method = "ward.D" )  |>
  as.dendrogram()

hc5_env <- hclust(distances_env_5, method = "ward.D" )  |>
  as.dendrogram()

plot(hc1_env)
plot(hc2_env)
plot(hc3_env)
plot(hc4_env)
plot(hc5_env)

## I plot the dendogram and the wavelets results in a heatmap and observe the clustering analysis----

year_labels <- c("2004", "2005", "2006", '2007', '2008', '2009', '2010',
                 '2011', '2012', '2013') ## create the labels with correspond with the sample num
### d1------
heatmap_data_l <- time_series_env_1 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'env_var')

env_order <- order.dendrogram(hc1_env)

dendro <- ggdendrogram(data = hc1_env, rotate = T)+
  scale_y_reverse()+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5))

heatmap_data_l$env_var <- factor(x = heatmap_data_l$env_var,
                                 levels = heatmap_data_l$env_var[env_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = env_var)) +
  geom_tile(aes(fill = as.numeric(value))) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0)+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  scale_y_discrete(labels = labs_env)+
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text.y = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'))

### d2------
heatmap_data_l <- time_series_env_1 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'env_var')

env_order <- order.dendrogram(hc2_env)

dendro <- ggdendrogram(data = hc2_env, rotate = T)+
  scale_y_reverse()+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5))

heatmap_data_l$env_var <- factor(x = heatmap_data_l$env_var,
                                 levels = heatmap_data_l$env_var[env_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = env_var)) +
  geom_tile(aes(fill = as.numeric(value))) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0)+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  scale_y_discrete(labels = labs_env)+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text.y = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'))

### d3------
heatmap_data_l <- time_series_env_3 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'env_var')

env_order <- order.dendrogram(hc3_env)

dendro <- ggdendrogram(data = hc3_env, rotate = T)+
  scale_y_reverse()+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5))

heatmap_data_l$env_var <- factor(x = heatmap_data_l$env_var,
                                 levels = heatmap_data_l$env_var[env_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = env_var)) +
  geom_tile(aes(fill = as.numeric(value))) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0)+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  scale_y_discrete(labels = labs_env)+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text.y = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'))

### d4------
heatmap_data_l <- time_series_env_4 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'env_var')

env_order <- order.dendrogram(hc4_env)

dendro <- ggdendrogram(data = hc4_env, rotate = T)+
  scale_y_reverse()+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5))

heatmap_data_l$env_var <- factor(x = heatmap_data_l$env_var,
                                 levels = heatmap_data_l$env_var[env_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = env_var)) +
  geom_tile(aes(fill = as.numeric(value))) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0)+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  scale_y_discrete(labels = labs_env)+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text.y = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'))

### s4------
heatmap_data_l <- time_series_env_5 |>
  dplyr::mutate(sample_num = row_number()) |>
  pivot_longer(cols = -c(decimal_date, sample_num), names_to = 'env_var')

env_order <- order.dendrogram(hc5_env)

dendro <- ggdendrogram(data = hc5_env, rotate = T)+
  scale_y_reverse()+
  theme(axis.text.y = element_text(size = 0), text = element_text(size = 5))

heatmap_data_l$env_var <- factor(x = heatmap_data_l$env_var,
                                 levels = heatmap_data_l$env_var[env_order], 
                                 ordered = TRUE)

heatmap_plot <- heatmap_data_l |>
  ggplot( aes(x = as.numeric(sample_num), y = env_var)) +
  geom_tile(aes(fill = as.numeric(value))) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0)+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_continuous(expand = c(0, 0), breaks = c(1,   13,   25,   37,   49,   61,   73,   85,   97,  109), 
                     labels =  year_labels) +
  labs(x = 'Time (Years)', y = '', fill = 'Wavelet\nresult')+
  scale_y_discrete(labels = labs_env)+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'right',
        axis.text.y = element_text(size = 5), text = element_text(size = 5),
        panel.border = element_blank(),  legend.key.size = unit(4, 'mm'))


# PERFORM THE CROSS-CORRELATIONS BETWEEN ENVIRONMENTAL VARIABLES AND ASVs AT DIFFERENT TRANSFORMATIONS-----

## FL----
time_series_fl_1 <- time_series_1 

time_series_fl_2 <- time_series_2

time_series_fl_3 <- time_series_3 

time_series_fl_4 <- time_series_4 

time_series_fl_5 <- time_series_5 

### D1 comparison-----
time_series_1$decimal_date == time_series_env_1$decimal_date # we sure that the tibbles are ordered in the same way. 

# time_series_1_mat <- time_series_1 |>
#   dplyr::select(-decimal_date) |>
#   dplyr::filter(if_all(everything(), ~ !is.na(.))) |>
#   as.matrix()
# 
# time_series_env_1 <- time_series_env_1 |>
#   mutate_if(is.character, as.numeric) 
# 
# time_series_env_1_mat <- time_series_env_1 |>
#   dplyr::select(-decimal_date) |>
#   mutate_if(is.character, as.numeric) |>
#   dplyr::filter(if_all(everything(), ~ !is.na(.))) |>
#   as.matrix()

time_series_env_1_ed <-  time_series_env_1 |>
  mutate_if(is.character, as.numeric) |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_fl_1_ed <- time_series_fl_1 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_1_ed) == dim(time_series_fl_1_ed) ##the first needs to be true

cor_1 <- cor(time_series_fl_1_ed, time_series_env_1_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_fl_1_ed, time_series_env_1_ed)

pdf(file = "results/figures/cross_corr_wavelets_env/cross_corr_env_fl_d1.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_fl_d1 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_fl_1), 
                                Colv = as.dendrogram(hc1_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels

dev.off()

# Add color legend
# image.plot(x = c(-1, 1), y = c(-1, 1), z = matrix(seq(-1, 1, length.out = 10), nrow = 1),
#            col = colorRampPalette(c("navy", "white", "red"))(100),
#            zlim = c(-1, 1),
#            legend.only = TRUE, legend.shrink = 0.6, legend.width = 0.5,
#            legend.mar = 5, legend.lab = "Correlation values")


#### D2 comparison----
time_series_fl_2$decimal_date == time_series_env_2$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_2_ed <-  time_series_env_2 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_fl_2_ed <- time_series_fl_2 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_2_ed) == dim(time_series_fl_2_ed) ##the first needs to be true

cor_2 <- cor(time_series_fl_2_ed, time_series_env_2_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_fl_2_ed, time_series_env_2_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_fl_2), Colv = as.dendrogram(hc2_env),
        col = colorRampPalette(c("navy", "white", "red"))(25),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.2)", name),
        xlab = '', ylab = '')

pdf(file = "results/figures/cross_corr_wavelets_env/cross_corr_env_fl_d2.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_fl_d2 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_fl_2), 
                                Colv = as.dendrogram(hc2_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()

#### D3 comparison----
time_series_fl_3$decimal_date == time_series_env_3$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_3_ed <-  time_series_env_3 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_fl_3_ed <- time_series_fl_3 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_3_ed) == dim(time_series_fl_3_ed) ##the first needs to be true

cor_3 <- cor(time_series_fl_3_ed, time_series_env_3_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_fl_3_ed, time_series_env_3_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_fl_3), Colv = as.dendrogram(hc3_env),
        col = colorRampPalette(c("navy", "white", "red"))(35),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.3)", name),
        xlab = '', ylab = '')

pdf(file = "../results/figures/cross_corr_wavelets_env/cross_corr_env_fl_d3.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_fl_d3 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_fl_3), 
                                Colv = as.dendrogram(hc3_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()

#### D4 comparison----
time_series_fl_4$decimal_date == time_series_env_4$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_4_ed <-  time_series_env_4 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_fl_4_ed <- time_series_fl_4 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_4_ed) == dim(time_series_fl_4_ed) ##the first needs to be true

cor_4 <- cor(time_series_fl_4_ed, time_series_env_4_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_fl_4_ed, time_series_env_4_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_fl_4), Colv = as.dendrogram(hc4_env),
        col = colorRampPalette(c("navy", "white", "red"))(45),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.4)", name),
        xlab = '', ylab = '')+
  theme(text = element_text(size = 6))+
  theme_bw()

pdf(file = "../results/figures/cross_corr_wavelets_env/cross_corr_env_fl_d4.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_fl_d4 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_fl_4), 
                                Colv = as.dendrogram(hc4_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()

#### D5 comparison----
time_series_fl_5$decimal_date == time_series_env_5$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_5_ed <-  time_series_env_5 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_fl_5_ed <- time_series_fl_5 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_5_ed) == dim(time_series_fl_5_ed) ##the first needs to be true

cor_5 <- cor(time_series_fl_5_ed, time_series_env_5_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_fl_5_ed, time_series_env_5_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_fl_5), Colv = as.dendrogram(hc5_env),
        col = colorRampPalette(c("navy", "white", "red"))(55),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.5)", name),
        xlab = '', ylab = '')+
  theme(text = element_text(size = 6))+
  theme_bw()

pdf(file = "../results/figures/cross_corr_wavelets_env/cross_corr_env_fl_s5.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_fl_s5 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_fl_5), 
                                Colv = as.dendrogram(hc5_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()


## PA----
time_series_pa_1 <- time_series_1 

time_series_pa_2 <- time_series_2

time_series_pa_3 <- time_series_3 

time_series_pa_4 <- time_series_4 

time_series_pa_5 <- time_series_5 

#### D1 comparison----
time_series_pa_1$decimal_date == time_series_env_1$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_1_ed <-  time_series_env_1 |>
  arrange(decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_pa_1_ed <- time_series_pa_1 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_1_ed) == dim(time_series_pa_1_ed) ##the first needs to be true

cor_1 <- cor(time_series_pa_1_ed, time_series_env_1_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_pa_1_ed, time_series_env_1_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_pa_1), Colv = as.dendrogram(hc1_env),
        col = colorRampPalette(c("navy", "white", "red"))(25),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.2)", name),
        xlab = '', ylab = '')

# Add color legend
image.plot(x = c(-1, 1), y = c(-1, 1), z = matrix(seq(-1, 1, length.out = 10), nrow = 1),
           col = colorRampPalette(c("navy", "white", "red"))(100),
           zlim = c(-1, 1),
           legend.only = TRUE, legend.shrink = 0.6, legend.width = 0.5,
           legend.mar = 5, legend.lab = "Correlation values")

pdf(file = "results/figures/cross_corr_wavelets_env/cross_corr_env_pa_d1.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_pa_d1 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_pa_1), 
                                Colv = as.dendrogram(hc1_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()

#### D2 comparison----
time_series_pa_2$decimal_date == time_series_env_2$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_2_ed <-  time_series_env_2 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_pa_2_ed <- time_series_pa_2 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_2_ed) == dim(time_series_pa_2_ed) ##the first needs to be true

cor_2 <- cor(time_series_pa_2_ed, time_series_env_2_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_pa_2_ed, time_series_env_2_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_pa_2), Colv = as.dendrogram(hc2_env),
        col = colorRampPalette(c("navy", "white", "red"))(25),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.2)", name),
        xlab = '', ylab = '')

pdf(file = "../results/figures/cross_corr_wavelets_env/cross_corr_env_pa_d2.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_pa_d2 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_pa_2), 
                                Colv = as.dendrogram(hc2_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()

#### D3 comparison----
time_series_pa_3$decimal_date == time_series_env_3$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_3_ed <-  time_series_env_3 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_pa_3_ed <- time_series_pa_3 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_3_ed) == dim(time_series_pa_3_ed) ##the first needs to be true

cor_3 <- cor(time_series_pa_3_ed, time_series_env_3_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_pa_3_ed, time_series_env_3_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_pa_3), Colv = as.dendrogram(hc3_env),
        col = colorRampPalette(c("navy", "white", "red"))(35),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.3)", name),
        xlab = '', ylab = '')
pdf(file = "../results/figures/cross_corr_wavelets_env/cross_corr_env_pa_d3.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_pa_d3 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_pa_3), 
                                Colv = as.dendrogram(hc3_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()


#### D4 comparison----
time_series_pa_4$decimal_date == time_series_env_4$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_4_ed <-  time_series_env_4 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_pa_4_ed <- time_series_pa_4 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_4_ed) == dim(time_series_pa_4_ed) ##the first needs to be true

cor_4 <- cor(time_series_pa_4_ed, time_series_env_4_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_pa_4_ed, time_series_env_4_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_pa_4), Colv = as.dendrogram(hc4_env),
        col = colorRampPalette(c("navy", "white", "red"))(45),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.4)", name),
        xlab = '', ylab = '')+
  theme(text = element_text(size = 6))+
  theme_bw()

pdf(file = "results/figures/cross_corr_wavelets_env/cross_corr_env_pa_d4.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_pa_d4 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_pa_4), 
                                Colv = as.dendrogram(hc4_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()

#### S5 comparison----
time_series_pa_5$decimal_date == time_series_env_5$decimal_date # we sure that the tibbles are ordered in the same way. 

time_series_env_5_ed <-  time_series_env_5 |>
  dplyr::select(-decimal_date) |>
  mutate_if(is.character, as.numeric) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

time_series_pa_5_ed <- time_series_pa_5 |>
  dplyr::select(-decimal_date) |>
  dplyr::filter(if_all(everything(), ~ !is.na(.)))

dim(time_series_env_5_ed) == dim(time_series_pa_5_ed) ##the first needs to be true

cor_5 <- cor(time_series_pa_5_ed, time_series_env_5_ed) |>
  melt()

cross_correlation_matrices <- cor(time_series_pa_5_ed, time_series_env_5_ed)

heatmap(cross_correlation_matrices, Rowv = as.dendrogram(hc_pa_5), Colv = as.dendrogram(hc5_env),
        col = colorRampPalette(c("navy", "white", "red"))(55),
        scale = "none",
        #labCol = labs_env, 
        #main = paste("Cross-correlation (0.5)", name),
        xlab = '', ylab = '')+
  theme(text = element_text(size = 6))+
  theme_bw()

pdf(file = "../results/figures/cross_corr_wavelets_env/cross_corr_env_pa_s5.pdf", width = 14, height = 14)  # Adjust width and height as needed

cross_corr_env_pa_s5 <- heatmap(cross_correlation_matrices, 
                                Rowv = as.dendrogram(hc_pa_5), 
                                Colv = as.dendrogram(hc5_env),
                                col = colorRampPalette(c("navy", "white", "red"))(25),
                                scale = "none",
                                #labCol = labs_env, 
                                #main = paste("Cross-correlation (0.2)", name),
                                xlab = '', 
                                ylab = '',
                                margins = c(6, 12),
                                cexRow = 0.8,  # Adjust the size of row labels
                                cexCol = 0.8)  # Adjust the size of column labels
dev.off()






# ----------- ######## -------- CLUSTER ASVs BASASED ON THE CONTINUOUS WAVELET TRANSFROMATION ------- ########--------
### The cross-wavelet spectrum and coherency spectrum of two time series can be analyzed with function coherency analyze.

## my data needs to be in wide format 
wavelet_02_df_date_w <- wavelet_02_df_date |>
  pivot_wider(id_cols = 'date', values_from = 'abundance_value', names_from = 'asv_num')

wavelet_3_df_date_w <- wavelet_3_df_date |>
  pivot_wider(id_cols = 'date', values_from = 'abundance_value', names_from = 'asv_num')

library(WaveletComp)
library(utils)

bloo_02_filt <- bloo_02 |>
  dplyr::filter(!value %in% c('asv2', 'asv3', 'asv5', 'asv8'))

pairs <- combn(bloo_02_filt$value, 2)

# Initialize a list to store the coherency results
coherency_results <- list()

# Loop through each pair of time series
for (i in 1:ncol(pairs)) {
  ts1 <- wavelet_02_df_date_w[, pairs[1, i]]
  ts2 <- wavelet_02_df_date_w[, pairs[2, i]]
  
  # Analyze coherency between the two time series
  coherency_result <- analyze.coherency(my.data = cbind(ts1, ts2),
                                        loess.span = 0, 
                                        dt = 1, # number of observations per time unit
                                        dj = 1/12, 
                                        #type = "xwavelet",
                                        lowerPeriod = 2, 
                                        upperPeriod = 32, 
                                        make.pval = TRUE, method = "white.noise", params = NULL,
                                        n.sim = 100, 
                                        date.format = '%Y-%M-%d', date.tz = NULL, 
                                        verbose = TRUE)
  
  # Store the result
  coherency_results <- c(coherency_results, list(coherency_result))
}

library(biwavelet)

ts1 <- wavelet_02_df_date_w[, pairs[1]]
ts2 <- wavelet_02_df_date_w[, pairs[2]]

cross_wavelet_results <- crossWavelet(ts1, ts2, 
                                      scales = 1:32, 
                                      boundary = "reflection", 
                                      dj = 1/12)

## Now I extract the results from the coherency_results list ----
date_col <- wavelet_02_df_date_w$date |>
  as_tibble_col(column_name = 'date') |>
  dplyr::mutate(sample_id_num = 1:nrow(wavelet_02_df_date_w))

coherency_result[[1]]
coherency_result$Coherency

pairs_ed <- pairs |>
  as_tibble() |>
  t() |>
  as_tibble() |>
  dplyr::mutate(coherence_asvs = paste0(V1,'.', V2))

coherency_result_tibble <- coherency_result$Coherency |> 
  as_tibble()

colnames(coherency_result_tibble) <- pairs_ed$coherence_asvs

coherency_result_tibble <- coherency_result_tibble |>
  #dplyr::mutate(values = 'coherency')|>
  pivot_longer(cols= starts_with('asv'), values_to = 'coherency') 
  
coherency_result_pval_tibble <- coherency_result$Coherence.pval |>
  as_tibble() 

colnames(coherency_result_pval_tibble) <- pairs_ed$coherence_asvs

coherency_result_pval_tibble <- coherency_result_pval_tibble |>
  #dplyr::mutate(values = 'pvalue') |>
  pivot_longer(cols= starts_with('asv'), values_to = 'pvalue')

coherency_results_pvalue_tb_02 <- coherency_result_tibble |>
  bind_cols(coherency_result_pval_tibble) 

period_tb <- coherency_results$Period

significant_02_coherency <- coherency_results_pvalue_tb_02 |>
  dplyr::filter(pvalue < 0.05) |>
  dplyr::filter(!str_detect(coherency, 'NaN')) |>
  separate(coherency, into = c('coherency_real', 'radiant'),  sep = "[1-9](?=[+-][0-1.])", remove = FALSE)

significant_02_coherency |> 
  dplyr::mutate(radiant = str_replace(radiant, 'i', '')) |>
  dplyr::filter(as.numeric(radiant) < 0.2 &
                  as.numeric(radiant) > 0)

# # Loop through each coherency result
# for (i in 1:ncol(pairs)) {
#   # Extract the period and coherence values
#   period <- coherency_result$Period[[i]]
#   coherence <- coherency_result$Coherency[[i]]
#   
#   # Plot the coherence as a function of period
#   plot(period, coherence, type = "l", xlab = "Period (days)", ylab = "Coherence", main = paste("Coherency Plot - Pair", i), 
#        xlim = c(0, max(period)), ylim = c(0, 1))
# }
# 
# # Extract the coherency result for a specific pair of time series
# coherency_result <- coherency_results[[1]]
# 
# # Extract the period and coherence values
# period <- coherency_result$Period
# coherence <- coherency_result$Coherency
# significance <- coherency_result$Coherency.pval
# 
# # Plot the wavelet coherence spectra and their significance
# plot.coherency(period, coherence, significance, xlab = "Period (days)", ylab = "Coherence")

# # Loop through each coherency result
# for (i in length(1:ncol(pairs))) {
#   # Extract the period, coherence, and significance values
#   period <- coherency_results[[i]]$Period
#   coherency <- coherency_results[[i]]$Coherency
#   coherency <- coherency |>
#     as_tibble() |>
#     separate(coherency, into = c('coherency_real', 'radiant'),  sep = "[1-9](?=[+-][0-1.])", remove = FALSE)
#   
#   significance <- coherency_results[[i]]$Coherence.pval
#   
#   # Plot the coherence as a function of period
#   plot <- ggplot(data.frame(period, coherence, significance), aes(x = period, y = coherence)) + 
#     geom_line() + 
#     geom_hline(yintercept = significance, linetype = "dashed") + 
#     labs(x = "Period (days)", y = "Coherence", title = paste("Coherency Plot - Pair", i))
#   
#   plot
# }

## i would like to plot the results 
wt.image(coherency_results[[1]], 
         my.series = 1,
         color.key = "interval", 
         main = paste0("Wavelet Coherence - ", pairs_ed$coherence_asvs[[1]]), 
         periodlab = "Period (Months)",
         legend.params = list(width = 4, shrink = 0.9, mar = 5.1,
                              n.ticks = 6,
                              label.digits = 1, label.format = "f",
                              lab = "Coherence levels", lab.line = 2.5),
         #color.palette = 'rainbow(n.levels, start = 0, end = 0.7)',
         date.format = '%Y-%M-%d', 
         label.time.axis = 'Time',
         date.tz = NULL,
         verbose = TRUE,
         show.date = TRUE)

wt.image(coherency_results[[1]], 
         my.series = 2,
         color.key = "interval", 
         main = paste0("Wavelet Coherence - ", pairs_ed$coherence_asvs[[1]]), 
         periodlab = "Period (Months)",
         legend.params = list(width = 4, shrink = 0.9, mar = 5.1,
                              n.ticks = 6,
                              label.digits = 1, label.format = "f",
                              lab = "Coherence levels", lab.line = 2.5),
         #color.palette = 'rainbow(n.levels, start = 0, end = 0.7)',
         date.format = '%Y-%M-%d', 
         label.time.axis = 'Time',
         date.tz = NULL,
         verbose = TRUE,
         show.date = TRUE)

wt.avg(coherency_results[[1]],
       my.series = 1,
       siglvl = 0.05, sigcol = "red", 
       periodlab = "period (days)")

# Loop through each coherency result
for (i in 1:length(coherency_results)) {
  
  fraction <-  '0.2'
  # Extract the period, coherence, and significance values

  # Plot the coherence as a function of period using wt.image
  wt.image(coherency_results[[i]], color.key = "interval", 
           main = paste0("Wavelet Coherence - Pair", i), 
           legend.params = list(lab = "coherence levels"),
           periodlab = "period (days)" )
  
  # Plot the wavelet average of coherence using wt.avg
  wt.avg(coherency_results[[i]], siglvl = 0.05, sigcol = "red", 
         periodlab = "period (days)")
  
  # Create the directory if it doesn't exist
  dir.create(paste0("results/figures/wavelets_coherence_", fraction, "_ed/"), recursive = TRUE)
  
  # Save the wavelet coherence plot
  pdf(file = paste0("results/figures/wavelets_coherence_", fraction, "_ed/", i, "_wavelet_coherence.pdf"), 
      width = 40, height = 34
      # , units = "in"
  )
  wt.image(coherency_results[[i]], color.key = "interval", 
           main = paste0("Wavelet Coherence - Pair", i), 
           legend.params = list(lab = "coherence levels"),
           periodlab = "period (days)",
           show.date = T)
  dev.off()
  
  # Save the wavelet average of coherence plot
  pdf(file = paste0("results/figures/wavelets_coherence_", fraction, "_ed/", i, "_wavelet_average_coherence.pdf"), 
      width = 40, height = 34#, units = "in"
  )
  wt.avg(coherency_results[[i]], siglvl = 0.05, sigcol = "red", 
         periodlab = "period (days)")
  dev.off()
}

## pa fraction ----
 # Generate all possible pairs of time series columns 
pairs <- combn(names(wavelet_3_df_date_w)[2:ncol(wavelet_3_df_date_w)], 2)

# Initialize a list to store the coherency results
coherency_results <- list()

# Loop through each pair of time series
for (i in 1:nrow(pairs)) {
  ts1 <- wavelet_3_df_date_w[, pairs[1, i]]
  ts2 <- wavelet_3_df_date_w[, pairs[2, i]]
  
  # Analyze coherency between the two time series
  coherency_result <- analyze.coherency(my.data = cbind(ts1, ts2),
                                        loess.span = 0, 
                                        dt = 1, # number of observations per time unit
                                        dj = 1/12, 
                                        #type = "xwavelet",
                                        lowerPeriod = 2, 
                                        upperPeriod = 32, 
                                        make.pval = TRUE, method = "white.noise", params = NULL,
                                        n.sim = 100, 
                                        date.format = '%Y-%M-%d', date.tz = NULL, 
                                        verbose = TRUE)
  
  # Store the result
  coherency_results <- c(coherency_results, list(coherency_result))
}

## Now I extract the results from the coherency_results list ----
date_col <- wavelet_3_df_date_w$date |>
  as_tibble_col(column_name = 'date') |>
  dplyr::mutate(sample_id_num = 1:nrow(wavelet_3_df_date_w))

coherency_result[[1]]
coherency_result$Coherency

pairs_ed <- pairs |>
  as_tibble() |>
  t() |>
  as_tibble() |>
  dplyr::mutate(coherence_asvs = paste0(V1,'.', V2))

coherency_result_tibble <- coherency_result$Coherency |> 
  as_tibble()

colnames(coherency_result_tibble) <- pairs_ed$coherence_asvs

coherency_result_tibble <- coherency_result_tibble |>
  #dplyr::mutate(values = 'coherency')|>
  pivot_longer(cols= starts_with('asv'), values_to = 'coherency') 

coherency_result_pval_tibble <- coherency_result$Coherence.pval |>
  as_tibble() 

colnames(coherency_result_pval_tibble) <- pairs_ed$coherence_asvs

coherency_result_pval_tibble <- coherency_result_pval_tibble |>
  #dplyr::mutate(values = 'pvalue') |>
  pivot_longer(cols = starts_with('asv'), values_to = 'pvalue')

coherency_results_pvalue_tb_3 <- coherency_result_tibble |>
  bind_cols(coherency_result_pval_tibble) 

period_tb <- coherency_results$Period

significant_3_coherency <- coherency_results_pvalue_tb_3 |>
  dplyr::filter(pvalue < 0.05) |>
  dplyr::filter(!str_detect(coherency, 'NaN')) |>
  separate(coherency, into = c('coherency_real', 'radiant'),  sep = "[1-9](?=[+-][0-1.])", remove = FALSE)

significant_3_coherency |> 
  dplyr::mutate(radiant = str_replace(radiant, 'i', '')) |>
  dplyr::filter(as.numeric(radiant) < 3 &
                  as.numeric(radiant) > 0)

clusters_3 <- significant_3_coherency |>
  distinct(name...1)

## plot the results from FL -----
significant_02_coherency

asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2') |>
  dplyr::filter(asv_num %in% c('asv27', 'asv555', 'asv237', 'asv563', 'asv114', 'asv11')) |>
  dplyr::filter(abundance_type == 'rclr') |>
  dplyr::mutate(cluster = case_when(asv_num %in% c('asv27', 'asv555', 'asv237', 'asv563') ~ 'cluster1',
                                    asv_num %in% c('asv114', 'asv11') ~ 'cluster2')) |>
  group_by(date, asv_num_f, fraction, cluster, order_f, family_f) |>
  dplyr::reframe(abund_max = sum(abundance_value)) |>
  ggplot(aes(date, abund_max))+
  #scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette_family_assigned_bloo)+
  labs(y = 'rCLR', x = 'Time', fill = 'Family')+
  geom_area(aes(group = asv_num_f, fill = family_f), position = 'stack')+
  #facet_wrap(vars(harbor_group), ncol = 1, scales = 'free_y')+
  facet_wrap(vars(cluster),  scales = 'free_y', ncol = 2)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom',
        aspect.ratio = 4/11,
        panel.grid = element_blank(), text = element_text(size = 12),
        strip.text = element_text(margin = margin(2, 2, 2, 2)),
        plot.margin = unit(c(0.2, 5, 0.5, 0.5), "cm"),
        legend.key.size = unit(3, 'mm'))
