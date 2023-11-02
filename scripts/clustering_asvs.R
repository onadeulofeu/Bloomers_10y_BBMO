library(tidyverse)
## Fuzzy C-Means
### FCM is a soft clustering alogrithm proposed by Bezdek. Unlike K-means algorithm in which each data object is the
### the member of only one cluster, a data object is the member of all clusters with varyinh degrees of fuzzy membership
### between 0 and 1 in FCM. Hence, the data objects closer to the centers of clusters have higher degrees of membership 
### than objects scattered in the borders of clusters

library(ppclust)

## a hierarchical clustering method is one which works by partitioning the data into groups with
## increasing similar features.

## prior to clustering we need to decide the adequate number of clusters we need for our dataset:
## for this purpose we not to implement a clustering algorithm using a range of possible numbers 
## of cluster, and then comparison of these indices will indicate which number has a high degree
## of fit without over-fitting.
## Although there are many internal indexes that have originally been proposed
## for working with hard membership degrees preoduced by the K-means and its variants,
## most of these indexes cannot be used for fuzzy clustering results.


asv_tab_all_bloo_z_tax |>
  colnames()

# I want to select only the bloomers that bloom en each particular fraction
bloo_3 # vector with ASV number of potential bloomers in PA fraction
bloo_02 # vector with ASV number of potential bloomers in FL fraction

bloo_3_w <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '3' &
                  asv_num %in% bloo_3) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(asv_num, sample_id, abundance_value) |>
  pivot_wider(id_cols = sample_id, values_from = abundance_value, names_from = asv_num) |>
  as.data.frame(row.names = sample_id)

bloo_02_w <- asv_tab_all_bloo_z_tax |>
  dplyr::filter(fraction == '0.2' &
                  asv_num %in% bloo_02) |>
  dplyr::filter(abundance_type == 'relative_abundance') |>
  dplyr::select(asv_num, sample_id, abundance_value) |>
  pivot_wider(id_cols = sample_id, values_from = abundance_value, names_from = asv_num) |>
  as.data.frame(row.names = sample_id)

# test |>
#   class()
# 
# test |>
#   glimpse()

par(mar = c(1, 1, 1, 1)) # we use this to fix the error with margins too large. 

cluster_3_15c <- fcm(bloo_3_w[,-1], centers = 15, m = 4)
cluster_3_10c <- fcm(bloo_3_w[,-1], centers = 10)
cluster_3_5c <- fcm(bloo_3_w[,-1], centers = 5)
cluster_3_4c <- fcm(bloo_3_w[,-1], centers = 4)
cluster_3_3c <- fcm(bloo_3_w[,-1], centers = 3)

summary(cluster_3_15c)
summary(cluster_3_10c)
summary(cluster_3_5c)
summary(cluster_3_4c)
summary(cluster_3_3c)

plotcluster(cluster_3_15c, cp=5)
plotcluster(cluster_3_10c, cp=1)
plotcluster(cluster_3_5c, cp=1)
plotcluster(cluster_3_4c, cp=1)
plotcluster(cluster_3_3c, cp=1)

res.fcm2_15c <- ppclust2(cluster_3_15c, "kmeans")
res.fcm2_10c <- ppclust2(cluster_3_10c, "kmeans")
res.fcm2_5c <- ppclust2(cluster_3_5c, "kmeans")
res.fcm2_4c <- ppclust2(cluster_3_4c, "kmeans")
res.fcm2_3c <- ppclust2(cluster_3_3c, "kmeans")

factoextra::fviz_cluster(res.fcm2_15c, data = bloo_3_w[,-1], 
                         ellipse.type = "convex",
                         palette = "jco",
                         repel = TRUE)

factoextra::fviz_cluster(res.fcm2_10c, data = bloo_3_w[,-1], 
                         ellipse.type = "convex",
                         palette = "jco",
                         repel = TRUE)

factoextra::fviz_cluster(res.fcm2_5c, data = bloo_3_w[,-1], 
                         ellipse.type = "convex",
                         palette = "jco",
                         repel = TRUE)

factoextra::fviz_cluster(res.fcm2_4c, data = bloo_3_w[,-1], 
                         ellipse.type = "convex",
                         palette = "jco",
                         repel = TRUE)

factoextra::fviz_cluster(res.fcm2_3c, data = bloo_3_w[,-1], 
                         ellipse.type = "convex",
                         palette = "jco",
                         repel = TRUE)

library(cluster)
res.fcm3_15 <- ppclust2(cluster_3_15c, "fanny")
res.fcm3_10 <- ppclust2(cluster_3_10c, "fanny")
res.fcm3_5 <- ppclust2(cluster_3_5c, "fanny")
res.fcm3_4 <- ppclust2(cluster_3_4c, "fanny")
res.fcm3_3 <- ppclust2(cluster_3_3c, "fanny")

# cluster::clusplot(scale(test[,-1]), res.fcm3$cluster,  
#                   main = "Cluster plot of Bloomers in the PA fraction",
#                   color=TRUE, labels = 2, lines = 2, cex=1)

cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_15$cluster,  
                  main = "Cluster plot of Bloomers in the PA fraction",
                  color=TRUE, labels = 2, lines = 2, cex=1)

cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_10$cluster,  
                  main = "Cluster plot of Bloomers in the PA fraction",
                  color=TRUE, labels = 2, lines = 2, cex=1)

cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_5$cluster,  
                  main = "Cluster plot of Bloomers in the PA fraction",
                  color=TRUE, labels = 2, lines = 2, cex=1)

cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_4$cluster,  
                  main = "Cluster plot of Bloomers in the PA fraction",
                  color=TRUE, labels = 2, lines = 2, cex=1)

cluster::clusplot(scale(bloo_3_w[,-1]), res.fcm3_3$cluster,  
                  main = "Cluster plot of Bloomers in the PA fraction",
                  color=TRUE, labels = 2, lines = 2, cex=1)

## Validation of the clustering results find the best clustering option
library(fclust)

## 15 clusters
res.fcm4 <- ppclust2(cluster_3_15c, "fclust")
idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
idxpe <- PE(res.fcm4$U)
idxpc <- PC(res.fcm4$U)
idxmpc <- MPC(res.fcm4$U)

cat("Partition Entropy: ", idxpe)
cat("Partition Coefficient: ", idxpc)
cat("Modified Partition Coefficient: ", idxmpc)
### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
### the index takes the maximum value.
cat("Fuzzy Silhouette Index: ", idxsf)

## 10 clusters
res.fcm4 <- ppclust2(cluster_3_10c, "fclust")
idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
idxpe <- PE(res.fcm4$U)
idxpc <- PC(res.fcm4$U)
idxmpc <- MPC(res.fcm4$U)

cat("Partition Entropy: ", idxpe)
cat("Partition Coefficient: ", idxpc)
cat("Modified Partition Coefficient: ", idxmpc)
### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
### the index takes the maximum value.
cat("Fuzzy Silhouette Index: ", idxsf)

## 5 clusters
res.fcm4 <- ppclust2(cluster_3_5c, "fclust")
idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
idxpe <- PE(res.fcm4$U)
idxpc <- PC(res.fcm4$U)
idxmpc <- MPC(res.fcm4$U)

cat("Partition Entropy: ", idxpe)
cat("Partition Coefficient: ", idxpc)
cat("Modified Partition Coefficient: ", idxmpc)
### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
### the index takes the maximum value.
cat("Fuzzy Silhouette Index: ", idxsf)

## 4 clusters
res.fcm4 <- ppclust2(cluster_3_4c, "fclust")
idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
idxpe <- PE(res.fcm4$U)
idxpc <- PC(res.fcm4$U)
idxmpc <- MPC(res.fcm4$U)

cat("Partition Entropy: ", idxpe)
cat("Partition Coefficient: ", idxpc)
cat("Modified Partition Coefficient: ", idxmpc)
### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
### the index takes the maximum value.
cat("Fuzzy Silhouette Index: ", idxsf)

## 3 clusters
res.fcm4 <- ppclust2(cluster_3_3c, "fclust")
idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
idxpe <- PE(res.fcm4$U)
idxpc <- PC(res.fcm4$U)
idxmpc <- MPC(res.fcm4$U)

cat("Partition Entropy: ", idxpe)
cat("Partition Coefficient: ", idxpc)
cat("Modified Partition Coefficient: ", idxmpc)
### The fuzzy silhouete Index (FSI) the optimal number of clusters (k) is such that the 
### the index takes the maximum value.
cat("Fuzzy Silhouette Index: ", idxsf)

##group by clustering ASVs and then plot using streamchart or area-----
### 3 clusters
cluster_membership <- cluster_3_4c$
cluster_3_3c$u
cluster_membership <- cluster_3_3c$v |>
  as.data.frame() |>
  rownames_to_column(var = 'cluster_id') |>
  t() 

cluster_id <- cluster_membership[1,]

cluster_membership |>
  colnames() <- cluster_id 

cluster_membership <- cluster_membership[-1,]

cluster_membership <- cluster_membership |>
  as.data.frame() |>
  rownames_to_column(var = 'asv_num') |>
  as_tibble() |>
  dplyr::mutate(across(-asv_num, as.numeric)) |>
  dplyr::mutate(cluster_belonging = case_when("Cluster 1" > "Cluster 2" & "Cluster 1" > "Cluster 3" ~ 'cluster_1', 
                                              "Cluster 2" > "Cluster 1" &  "Cluster 2" > "Cluster 3" ~ 'cluster_2', 
                                              "Cluster 3" > "Cluster 1" & "Cluster 3" > "Cluster 2" ~ 'cluster_3'))
## problem when they are equal then it decides for cluster 3...


  # pivot_longer(cols = cluster_id, values_to = 'cluster_member_ship_degree', names_to = 'cluster_id') |>
  # #group_by(asv_num) |>
  # dplyr::filter(case_when(cluster_member_ship_degree) )
### 4 clusters
cluster_3_4c$u
cluster_membership <- cluster_3_4c$v |>
  as.data.frame() |>
  rownames_to_column(var = 'cluster_id') |>
  t() 

cluster_id <- cluster_membership[1,]

cluster_membership |>
  colnames() <- cluster_id 

cluster_membership <- cluster_membership[-1,]

cluster_membership <- cluster_membership |>
  as.data.frame() |>
  rownames_to_column(var = 'asv_num') |>
  as_tibble() |>
  dplyr::mutate(across(-asv_num, as.numeric)) |>
  dplyr::mutate(cluster_belonging = case_when("Cluster 1" > "Cluster 2" & "Cluster 1" > "Cluster 3" ~ 'cluster_1', 
                                              "Cluster 2" > "Cluster 1" &  "Cluster 2" > "Cluster 3" ~ 'cluster_2', 
                                              "Cluster 3" > "Cluster 1" & "Cluster 3" > "Cluster 2" ~ 'cluster_3',
                                              "Cluster 4" > "Clutser 1" & "Cluster 4" > "Clutser 2" & "Cluster 4" > "Clutser 3" ~ 'cluster_4'))
