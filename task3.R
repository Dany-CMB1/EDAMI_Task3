
#Lab3 - Clustering

#authors:
# --- Daniel BANNISTER
# --- Achille NOTRE DAME

library(dbscan)
library(fpc)
library(cluster)
library(factoextra)
library(ggplot2)

######## Cardiotocography ########
#---------- Data import ----------
#http://archive.ics.uci.edu/ml/datasets/Cardiotocography

download.file('https://staff.elka.pw.edu.pl/~rbembeni/dane/cardioto_all_corr.csv','cardioto_all_corr.csv')
ctg_all <- read.csv("cardioto_all_corr.csv",row.names = 1)

download.file('https://staff.elka.pw.edu.pl/~rbembeni/dane/cardioto_noClass_corr.csv','cardioto_noClass_corr.csv')
ctg_noClass <- read.csv("cardioto_noClass_corr.csv",row.names = 1)

# Default params Kmeans on unpreprocessed data
ctg_noClass.kmeans_basic <- eclust(ctg_noClass, "kmeans", k=max(ctg_all$CLASS), graph=TRUE)

# Silhouette evaluation
silinfo <- fviz_silhouette(ctg_noClass.kmeans_basic, palette="jco")
basic_avg_width <- mean(silinfo$data$sil_width)
basic_avg_width
#Average silhouette width: 0.2001795

# Rand index
clusters_stats <- cluster.stats(dist(ctg_noClass), ctg_all$CLASS, ctg_noClass.kmeans_basic$cluster)
basic_rand_index <- clusters_stats$corrected.rand
basic_rand_index
# Rand Index: 0.11

#---------- DATA PREPROCESSING ----------
colnames(ctg_noClass)
nrow(ctg_noClass)
colnames(ctg_all)
nrow(ctg_all)

# We look for and remove rows containing null values
colSums(is.na(ctg_noClass))
data <- na.omit(ctg_noClass)
colSums(is.na(ctg_all))
data <- na.omit(ctg_all)

# We search for and remove duplicated lines
sum(duplicated(ctg_noClass))
data <- data[!duplicated(ctg_noClass), ]
sum(duplicated(ctg_all))
data <- data[!duplicated(ctg_all), ]
# We now have a clean data set to work with

#---------- Data statistics and analysis ----------

# For ctg_all
# Summary of the columns
summary(ctg_all)
# Structure of the colmuns and their data types
str(ctg_all)
# Histograms for numeric columns
hist(ctg_all$Mean, main = "Histogram of Mean", xlab = "Mean")
# Boxplots to detect potential outliers
boxplot(ctg_all$Variance, main = "Boxplot of Variance", horizontal = TRUE)
# We count unique values
table(ctg_all$CLASS)
# Proportion of each category
prop.table(table(ctg_all$NSP))

#---------- Clustering ----------
MAX_NUM_GROUPS <- 15

#---------- Kmeans ----------
kmeans.rand_index <- vector(mode = "double" ,length = MAX_NUM_GROUPS)


for (i in 1:MAX_NUM_GROUPS) {
  kmean <- eclust(ctg_noClass, "kmeans", k=i, graph=FALSE)
  # Rand index
  clusters_stats <- cluster.stats(dist(ctg_noClass), ctg_all$CLASS, kmean$cluster)
  kmeans.rand_index[i] <- clusters_stats$corrected.rand
}

# Average silouhette width per number of groups => find optimal number of groups
fviz_nbclust(ctg_noClass, kmeans, method = "silhouette") +
  theme_classic() +
  ggtitle("Average Silhouette - Kmeans") +
  geom_hline(yintercept = basic_avg_width, linetype = "dashed", color = "red")

# ==> 2 is the optimal number of groups according to the average silhouette criteria
# The average silhouette seem to be quite slim in general, indicating poor clustering performance,
# but still better than the basic kmeans in most case (not in the case of k=10, reference grouping)

# Plot Rand Index per number of groups
plot(1:MAX_NUM_GROUPS, kmeans.rand_index, type = "b",
     xlab = "Number of groups",
     ylab = "Rand Index",
     main = "Rand Index per Number of Groups - Kmeans", # Use 'main' for the title
     pch = 16, col = "blue") # Add points and color for better visualization

# Add horizontal line for basic Rand Index
abline(h = basic_rand_index, col = "red", lty = 2) # Dashed red line for baseline

# ==> 6 is the optimal number of groups according to the Rand Index criteria
# Rand index is in average above the value found with basic kmeans, even for
# k=10 (basic grouping)

#---------- Partitioning around medoids ----------
pam.rand_index <- vector(mode = "double" ,length = MAX_NUM_GROUPS)

for (i in 1:MAX_NUM_GROUPS) {
  pam <- eclust(ctg_noClass, "pam", k=i, graph=FALSE)
  
  # Rand index
  clusters_stats <- cluster.stats(dist(ctg_noClass), ctg_all$CLASS, pam$cluster)
  pam.rand_index[i] <- clusters_stats$corrected.rand
}

# Average silouhette width per number of groups => find optimal number of groups
fviz_nbclust(ctg_noClass, pam, method = "silhouette") +
  theme_classic() +
  ggtitle("Average Silhouette - PAM") +
  geom_hline(yintercept = basic_avg_width, linetype = "dashed", color = "red")

# ==> Agan, 2 clusters seem to be optimal according to the Avergae silhouette criteria

# Plot Rand Index per number of groups
plot(1:MAX_NUM_GROUPS, pam.rand_index, type = "b",
     xlab = "Number of groups",
     ylab = "Rand Index",
     main = "Rand Index per Number of Groups - PAM", # Use 'main' for the title
     pch = 16, col = "blue") # Add points and color for better visualization

# Add horizontal line for basic Rand Index
abline(h = basic_rand_index, col = "red", lty = 2) # Dashed red line for baseline

# ==> 6 clusters is optimal according to the Rand Index criteria
# Same general conclusions as for the kmeans algorithm
# PAM does not solve our poor clustering problem

#---------- DBSCAN ----------
EPSILON <- seq(1, 10, by = 1.5)
MINPTS <- seq(5,10, by =1)

# Initialize results storage
results <- data.frame(eps = numeric(), minPts = numeric(),
                      avg_sil_width = numeric(), rand_index = numeric())

# Nested loops to iterate over eps and MinPts
for (eps in EPSILON) {
  for (minPts in MINPTS) {
    
    # Run DBSCAN
    dbscan_result <- dbscan(ctg_noClass, eps = eps,MinPts = minPts)
    
    # Exclude noise points (cluster 0)
    clusters <- dbscan_result$cluster
    valid_points <- clusters != 0
    
    # Silhouette calculation (only for valid points)
    if (sum(valid_points) > 1) { # Ensure enough points for silhouette calculation
      sil <- silhouette(clusters, dist(ctg_noClass))
      avg_sil_width <- mean(sil[, 3])
    } else {
      avg_sil_width <- NA # Not enough points for silhouette
    }
    
    # Rand index calculation
    clusters_stats <- cluster.stats(dist(ctg_noClass), ctg_all$CLASS, clusters)
    rand_index <- clusters_stats$corrected.rand
    
    # Store results
    results <- rbind(results, data.frame(eps = eps, minPts = minPts,
                                         avg_sil_width = avg_sil_width,
                                         rand_index = rand_index))
  }
}

# Print results
print(results)

# Optional: Find the best parameter combination based on silhouette or Rand index
best_sil <- results[which.max(results$avg_sil_width), ]
best_rand <- results[which.max(results$rand_index), ]

print("Best parameters based on silhouette width:")
print(best_sil)

print("Best parameters based on Rand index:")
print(best_rand)

# DBSCAN seems to perform terribaly on this database. This might be because there
# are thre points (2127,2128 and 2129) that are "far" from the other points
# cf. l.23 - ctg_noClass.kmeans_basic <- eclust(ctg_noClass, "kmeans", k=max(ctg_all$CLASS), graph=TRUE)




######## Wine ########
#---------- Data import ----------
download.file('http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv', 'wine_red.csv');
download.file('http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv', 'wine_white.csv');
wineRed_ds = read.table("wine_red.csv", header = TRUE, sep=";", na.strings= "*")
wineRed_dsC <- wineRed_ds[,-12]
View(wineRed_ds)
summary(wineRed_ds)


#---------- DATA PREPROCESSING ----------
colnames(wineRed_ds)
nrow(wineRed_ds)

# We look for and remove rows containing null values
colSums(is.na(wineRed_ds))
data <- na.omit(wineRed_ds)

# We search for and remove duplicated lines
sum(duplicated(wineRed_ds))
data <- data[!duplicated(wineRed_ds), ]

# We now have a clean data set to work with

#---------- Data statistics and analysis ----------

# Summary of the columns
summary(wineRed_ds)
# Structure of the colmuns and their data types
str(wineRed_ds)
# Histogram for alcohol in each wine
hist(wineRed_ds$alcohol, main = "Alcohol Content Distribution", xlab = "Alcohol")
# Count of wine quality ratings
table(wineRed_ds$quality)
# Bar plot for quality
barplot(table(wineRed_ds$quality), main = "Distribution of Wine Quality")
# Boxplot for fixed acidity in each wine
boxplot(wineRed_ds$fixed.acidity, main = "Fixed Acidity", ylab = "Fixed Acidity")

# Link betwen different columns
plot(wineRed_ds$fixed.acidity, wineRed_ds$quality,
     main = "Fixed Acidity vs Quality", xlab = "Fixed Acidity", ylab = "Quality")
plot(wineRed_ds$alcohol, wineRed_ds$quality,
     main = "Alcohol vs Quality", xlab = "Alcohol", ylab = "Quality")


#---------- Clustering ----------
#example - wines
distC = dist(wineRed_dsC)
card.kmeans = kmeans(distC,6)
res3 = table(wineRed_ds$quality,card.kmeans$cluster )
res3

#---------- K-means ----------

# We will try different clustering methods on our wine dataset
# First, we remove quality field because we don't want to use it here, maybe to verify our results later
wineRed_ds_modified <- subset(wineRed_ds, select = -quality)

# We normalize our data
wineRed_ds_scaled <- scale(wineRed_ds_modified)

# We apply k-means algorithm
kmeans_result <- kmeans(wineRed_ds_scaled, centers = 6, nstart = 25)

# We add a field in our dataset indicating in which cluster is the wine
wineRed_ds$cluster <- kmeans_result$cluster

# We visualize our results, using a 2D projection
# We apply a PCA to have a 2D projection
pca_result <- prcomp(wineRed_ds_scaled, center = TRUE, scale. = TRUE)
pca_data <- data.frame(pca_result$x[, 1:2], Cluster = as.factor(wineRed_ds$cluster))

# We plot our clusters depending on those two components
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "Clustering K-means",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

# We compute the silhouette index in order to evaluate our clustering
silhouette_values <- silhouette(kmeans_result$cluster, dist(wineRed_ds_scaled))
cat("Mean Silhouette Index : ", mean(silhouette_values[, 3]), "\n")
# Mean Silhouette Index : 0.194333, we have a respectable clustering but it's not perfect at all

#---------- Cluster linkage ----------

wineRed_ds_modified <- subset(wineRed_ds, select = -quality)
wineRed_ds_scaled <- scale(wineRed_ds_modified)

# We compute distances matrix
distance_matrix <- dist(wineRed_ds_scaled)

# We apply the clustering method
hclust_result <- hclust(distance_matrix, method = "complete")

# We cut our tree into 6 sub trees
clusters <- cutree(hclust_result, k = 6)

# We compute silouette indexes
silhouette_values <- silhouette(clusters, distance_matrix)
cat("Mean Silhouette Index : ", mean(silhouette_values[, 3]), "\n")
# Mean Silhouette Index : 0.177395, here again it's acceptable but not perfect

# We try to visualize our tree
plot(hclust_result, labels = FALSE, main = "Cluster linkage")
rect.hclust(hclust_result, k = 6, border = 2:7) # Ajout des boîtes autour des clusters

#---------- DBSCAN ----------

wineRed_ds_modified <- subset(wineRed_ds, select = -quality)
wineRed_ds_scaled <- scale(wineRed_ds_modified)

EPSILON <- seq(0.1, 2, by = 0.15)
MINPTS <- seq(3,10, by =1)

# Initialize results storage
results <- data.frame(eps = numeric(), minPts = numeric(),
                      avg_sil_width = numeric(), rand_index = numeric())

# Nested loops to iterate over eps and MinPts
for (eps in EPSILON) {
  for (minPts in MINPTS) {
    
    # Run DBSCAN
    dbscan_result <- dbscan(wineRed_ds_scaled, eps = eps,MinPts = minPts)
    
    # Exclude noise points (cluster 0)
    clusters <- dbscan_result$cluster
    valid_points <- clusters != 0
    
    # Silhouette calculation (only for valid points)
    if (sum(valid_points) > 1) { # Ensure enough points for silhouette calculation
      sil <- silhouette(clusters, dist(wineRed_ds_scaled))
      avg_sil_width <- mean(sil[, 3])
    } else {
      avg_sil_width <- NA # Not enough points for silhouette
    }
    
    # Store results
    results <- rbind(results, data.frame(eps = eps, minPts = minPts,
                                         avg_sil_width = avg_sil_width))
  }
}

# Print results
print(results)

# Optional: Find the best parameter combination based on silhouette or Rand index
best_sil <- results[which.max(results$avg_sil_width), ]

print("Best parameters based on silhouette width:")
print(best_sil)


# Manually, we see that we got the best results for eps = 1.75 and minPts = 10
dbscan_result <- dbscan(wineRed_ds_scaled, eps = 1.75, MinPts = 10)

# Étape 4 : Résultats du clustering
print(dbscan_result)

# We add clusters field to our dataset
wineRed_ds$cluster <- dbscan_result$cluster

# We compute silhouette index
# We exclude outliers having cluster = -1
valid_points <- wineRed_ds_scaled[dbscan_result$cluster != -1, ]
valid_clusters <- dbscan_result$cluster[dbscan_result$cluster != -1]

# We calculate silhouette index
silhouette_values <- silhouette(valid_clusters, dist(valid_points))
cat("Mean silhouette index : ", mean(silhouette_values[, 3]), "\n")
# Mean silhouette index : 0.2725711

# We visualize our clustering with a 2D projection thanks to a PCA
pca_result <- prcomp(wineRed_ds_scaled, center = TRUE, scale. = TRUE)
pca_data <- data.frame(pca_result$x[, 1:2], Cluster = as.factor(dbscan_result$cluster))
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "DBSCAN",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

# Conclusion WHINE : looking only at the silhouette indexes, we could think that DBSCAN gave us the better clustering
# but in fact, it gives us only two clusters, which is a bad result. It seems that the clusters are too compact to be
# detected automatically by DBSCAN algorithm. In this situation, it seems that k-means algorithm is better because we
# can specify the number of clusters we want to have at the end.
