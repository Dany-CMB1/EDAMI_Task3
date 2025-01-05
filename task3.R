
#Lab3 - Clustering

#authors:
# --- Daniel BANNISTER
# --- Achille NOTRE DAME

#data set for the laboratory task
#http://archive.ics.uci.edu/ml/datasets/Cardiotocography

download.file('https://staff.elka.pw.edu.pl/~rbembeni/dane/cardioto_noClass_corr.csv','cardioto_noClass_corr.csv')
ctg_noClass <- read.csv("cardioto_noClass_corr.csv",row.names = 1)

download.file('https://staff.elka.pw.edu.pl/~rbembeni/dane/cardioto_all_corr.csv','cardioto_all_corr.csv')
ctg_all <- read.csv("cardioto_all_corr.csv",row.names = 1)



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

# 1. For ctg_noClass
# Summary of the columns
summary(ctg_noClass)
# Structure of the colmuns and their data types
str(ctg_noClass)
# Histograms for numeric columns
hist(ctg_noClass$Mean, main = "Histogram of Mean", xlab = "Mean")
# Boxplots to detect potential outliers
boxplot(ctg_noClass$Variance, main = "Boxplot of Variance", horizontal = TRUE)

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
#example
distC = dist(ctg_noClass)
card.kmeans = kmeans(distC,10)
res3 = table(ctg_all$CLASS,card.kmeans$cluster )
res3


# ------------------------------------------------------------------------------


#wines
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
