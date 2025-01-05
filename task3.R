#Lab3 - Clustering

#authors:
# --- Daniel BANNISTER
# --- Achille NOTRE DAME

library(dbscan)
library(fpc)
library(cluster)
library(factoextra)

######## Cardiotocography ########
  #---------- Data import ----------
#http://archive.ics.uci.edu/ml/datasets/Cardiotocography

download.file('https://staff.elka.pw.edu.pl/~rbembeni/dane/cardioto_all_corr.csv','cardioto_all_corr.csv')
ctg_all <- read.csv("cardioto_all_corr.csv",row.names = 1)

download.file('https://staff.elka.pw.edu.pl/~rbembeni/dane/cardioto_noClass_corr.csv','cardioto_noClass_corr.csv')
ctg_noClass <- read.csv("cardioto_noClass_corr.csv",row.names = 1)

ctg_noClass.kmeans_basic <- eclust(ctg_noClass, "kmeans", k=max(ctg_all$CLASS), graph=TRUE)

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

# # 1. For ctg_noClass
# # Summary of the columns
# summary(ctg_noClass)
# # Structure of the colmuns and their data types
# str(ctg_noClass)
# # Histograms for numeric columns
# hist(ctg_noClass$Mean, main = "Histogram of Mean", xlab = "Mean")
# # Boxplots to detect potential outliers
# boxplot(ctg_noClass$Variance, main = "Boxplot of Variance", horizontal = TRUE)

#2. For ctg_all
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
kmeans.wss <- vector(mode = "integer" ,length = MAX_NUM_GROUPS)
kmeans.rand_index <- vector(mode = "double" ,length = MAX_NUM_GROUPS)


for (i in 1:MAX_NUM_GROUPS) {
  kmean <- eclust(ctg_noClass, "kmeans", k=i, graph=FALSE)

  # Rand index
  clusters_stats <- cluster.stats(dist(ctg_noClass), ctg_all$CLASS, kmean$cluster)
  kmeans.rand_index[i] <- clusters_stats$corrected.rand
  
  # total within-cluster sum of squares
  kmeans.wss[i] <- kmean$tot.withinss
}

# Average silouhette width per number of groups => find optimal number of groups
fviz_nbclust(ctg_noClass, kmeans, method = "silhouette")+theme_classic()

# total within-cluster sum of squares per number of groups
plot(1:MAX_NUM_GROUPS, kmeans.wss, type = "b", 
     xlab = "Number of groups", 
     ylab = "Total within-cluster sum of squares")

# Rand Index per number of groups
plot(1:MAX_NUM_GROUPS, kmeans.rand_index, type="b",
     xlab = "Number of groups",
     ylab = "Rand index")

    #---------- Partitioning around medoids ----------
pam.wss <- vector(mode = "integer" ,length = MAX_NUM_GROUPS)
pam.rand_index <- vector(mode = "double" ,length = MAX_NUM_GROUPS)


for (i in 1:MAX_NUM_GROUPS) {
  pam <- eclust(ctg_noClass, "pam", k=i, graph=FALSE)
  
  # Rand index
  clusters_stats <- cluster.stats(dist(ctg_noClass), ctg_all$CLASS, pam$cluster)
  pam.rand_index[i] <- clusters_stats$corrected.rand
  
  # total within-cluster sum of squares
  pam.wss[i] <- pam$tot.withinss
}

# Average silouhette width per number of groups
fviz_nbclust(ctg_noClass, pam, method = "silhouette")+theme_classic()

# total within-cluster sum of squares per number of groups
plot(1:MAX_NUM_GROUPS, kmeans.wss, type = "b", 
     xlab = "Number of groups", 
     ylab = "Total within-cluster sum of squares")

# Rand Index per number of groups
plot(1:MAX_NUM_GROUPS, kmeans.rand_index, type="b",
     xlab = "Number of groups",
     ylab = "Rand index")

    #---------- DBSCAN ----------
dbscan.eps <- vector(mode = "double" ,length = MAX_NUM_GROUPS)
dbscan.rand_index <- vector(mode = "double" ,length = MAX_NUM_GROUPS)


for (i in 1:MAX_NUM_GROUPS) {
  dbscan <- dbscan(ctg_noClass, i)
  
  # Rand index
  clusters_stats <- cluster.stats(dist(ctg_noClass), ctg_all$CLASS, dbscan$cluster)
  dbscan.rand_index[i] <- clusters_stats$corrected.rand
  
  dbscan$eps
}

# total within-cluster sum of squares per number of groups
plot(1:MAX_NUM_GROUPS, kmeans.wss, type = "b", 
     xlab = "Number of groups", 
     ylab = "Total within-cluster sum of squares")

# Rand Index per number of groups
plot(1:MAX_NUM_GROUPS, dbscan.rand_index, type="b",
     xlab = "Number of groups",
     ylab = "Rand index")

  
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
