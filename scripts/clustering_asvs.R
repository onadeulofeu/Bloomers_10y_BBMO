library(tidyverse)
library(fields) ##add the legend
library(magrittr)
library(ggdendro)
library(gridExtra) ## combine the dendogram and the heatmap
library(factoextra) ## visualize hierarchical clusters

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
# 
# wavelets_result_tibble_tax_3_biased <- read.csv2('data/wavelets_result_tibble_tax_3_biased.csv')
# wavelets_result_tibble_tax_02_biased <- read.csv2('data/wavelets_result_tibble_tax_02_biased.csv')

# I use the dataset that has the results with wavelets analysis without removing all samples that could be affected by the margin effects but some were
## trying to reduce the bias but not completely

wavelets_result_tibble_tax_3_biased <- read.csv('../data/wavelets_analysis/wavelets_result_ed_tibble_tax_3_biased_red.csv')
wavelets_result_tibble_tax_02_biased <- read.csv('../data/wavelets_analysis/wavelets_result_ed_tibble_tax_02_biased_red.csv')

asv_tab_all_bloo_z_tax <- read.csv2('data/detect_bloo/asv_tab_all_bloo_z_tax_new_assign_checked.csv') ##using dada2 classifier assign tax with silva 138.1

tax_bbmo_10y_new <- asv_tab_all_bloo_z_tax |>
  dplyr::select(asv_num, seq, domain, phylum, class, order, family, genus) |>
  distinct()

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
pdf("results/figures/hierarchical_clustering_scales/clust_d1_PA.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
pdf("results/figures/hierarchical_clustering_scales/clust_d2_PA.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
pdf("results/figures/hierarchical_clustering_scales/clust_d3_PA.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
pdf("results/figures/hierarchical_clustering_scales/clust_d4_PA.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
pdf("results/figures/hierarchical_clustering_scales/clust_s4_PA.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
sub_grp_3 <- cutree(hc3, k = 6)
sub_grp_4 <- cutree(hc4, k = 6)
sub_grp_5 <- cutree(hc4, k = 6)

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

palette_clustering <- c("#ffe37e",
"#ff7763",
"#e5a0be",
"#c23939",
"#4cb76a",
"#518535",
"#804c90",
"#6c97c3",
"#2b347a")

fviz_cluster(list(data = df_1, cluster = sub_grp_1))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_2 , cluster = sub_grp_2))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_3 , cluster = sub_grp_3))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_4 , cluster = sub_grp_4))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_5 , cluster = sub_grp_5))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))


###table with the label corresponding to their signals-------
### for each transformation used each bloomer has a different label only if it gives a clear strong signal 
#### this table will contain the asv num + for each transformation which label did it got

###The most important signals for our dataset have been the d1, d3, and s5 they are the ones that we will use for 
### labeling the bloomers.
df_5 |>
  str()

sub_grp_5 |>
  str()


## With this analysis we will pick a representative of each group to perform the random tree analysis-----






### The same for the FL fraction----
time_series_1 <- wavelets_result_tibble_tax_02_biased |>
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

time_series_2 <- wavelets_result_tibble_tax_02_biased |>
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

time_series_3 <- wavelets_result_tibble_tax_02_biased |>
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

time_series_4 <- wavelets_result_tibble_tax_02_biased |>
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

time_series_5 <- wavelets_result_tibble_tax_02_biased |>
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
pdf("results/figures/hierarchical_clustering_scales/clust_d1_FL.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
pdf("results/figures/hierarchical_clustering_scales/clust_d2_FL.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
pdf("results/figures/hierarchical_clustering_scales/clust_d3_FL.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
pdf("results/figures/hierarchical_clustering_scales/clust_d4_FL.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
pdf("results/figures/hierarchical_clustering_scales/clust_s4_FL.pdf", width = 8, height = 4)  # Adjust width and height as needed

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
sub_grp_3 <- cutree(hc3, k = 6)
sub_grp_4 <- cutree(hc4, k = 6)
sub_grp_5 <- cutree(hc4, k = 6)

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

fviz_cluster(list(data = df_1, cluster = sub_grp_1))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_2 , cluster = sub_grp_2))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_3 , cluster = sub_grp_3))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_4 , cluster = sub_grp_4))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))

fviz_cluster(list(data = df_5 , cluster = sub_grp_5))+
  scale_color_manual(values = palette_clustering)+
  scale_fill_manual(values = palette_clustering)+
  theme_bw()+
  labs(title = '')+
  theme(text = element_text(size = 6))


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
  dplyr::select(decimal_date, wavelets_result_ed asv_num, wavelets_transformation) |>
  dplyr::filter(wavelets_transformation == 'd1' ) |>
  dplyr::select(-wavelets_transformation) |>
  pivot_wider(names_from = asv_num, values_from = wavelets_result)

time_series_2 <- wavelets_result_tibble_tax_3_biased |>
  dplyr::select(decimal_date, wavelets_result_ed asv_num, wavelets_transformation) |>
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
