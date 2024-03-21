## ---------------------------
##
## Script Purpose: Initial exploratory PCA for analysis: 
##
## -----------------------------
##  This R script runs a PCA on the data sa a benchmark test. It does the following:
##.  1. Checks for missing values in proteomic dataset
##.  2. Imputes missing data: knn imputation in proteomic dataset
##.  3. Scales protein levels
##.  4. Performs PCA
##.  5. Generates PDFs of screeplots for supplement
##
##
## -----------------------------
##
## -----------------------------

setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics/UKB_PPP_Rscripts")

## -------LOAD DATASETS---------
MASTER_DATA <- read.csv("MASTER_DATA.csv") #contains all main data
ST3_UKB_proteins <- read.csv("ST3_UKB_proteins.csv") #contains data about proteomics data inc. protein panels

# Set 'eid' column as row names for the 'wm' dataset
rownames(MASTER_DATA) <- MASTER_DATA$eid

########################################### QUICK TEST RUNS on proteomic and FA data
# Generate a list of all the proteins of relevance, e.g.
Panel_info <- ST3_UKB_proteins %>% select(Assay.Target, Protein.panel)

all_proteins_vector <- as.character(unique(ST3_UKB_proteins$Assay.Target))

# Modify cases like 'HLA-DRA' and 'HLA-E' to 'HLA.DRA' and 'HLA.E'
all_proteins_vector <- gsub("-", ".", all_proteins_vector)

# n.b this is just the inflammatory ones, n = 368
all_protein_data <- MASTER_DATA %>% select(165:1627)

all_protein_data_non_imputed <- MASTER_DATA %>% select(eid, all_of(all_proteins_vector))

write.csv(all_protein_data_non_imputed, "all_protein_data_non_imputed.csv")


# remove any missing values dim(all_protein_data) = [1] 4644   368; dim(cleaned_all_protein_data) = 3423  368
# something is up here - reduces our dataset down to 173 ??
#all_protein_data_no_NA <- na.omit(all_protein_data)

# Calculate the number of missing values per column
missing_values <- sapply(imputed_data, function(x) sum(is.na(x)))

# Filter out columns with no missing values
missing_values <- missing_values[missing_values > 0]

# Check if there are any missing values at all
if (length(missing_values) > 0) {
  print("There are missing values in the dataframe.")
  print("Missing values per protein/column:")
  print(missing_values)
} else {
  print("There are no missing values in the dataframe.")
}

dim(all_protein_data)


###################################################################################

## Imputations of missing data
# KNN imputation

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("impute")
#https://www.rdocumentation.org/packages/impute/versions/1.46.0/topics/impute.knn

library(tidyverse)
library(impute)

raw_proteins <- read.csv('all_protein_data_non_imputed.csv')

## -------TIDY DATASETS---------
# Remove unnecessary columns like 'X' column from the data frame
# Set 'eid' column as row names for the 'wm' dataset and remove X

rownames(raw_proteins) <- raw_proteins$eid
raw_proteins <- raw_proteins[, !(names(raw_proteins) %in% c("X"))]
dim(raw_proteins)

IDs <- raw_proteins[c(1)]
data <- raw_proteins[c(2:1464)]
data <- as.matrix(data)
rownames(data) <- IDs$pseudo_ind_id
data <- t(data)
print(dim(data))
imputed <- impute.knn(data)

# Warning message:
#   In knnimp(x, k, maxmiss = rowmax, maxp = maxp) :
#   2 rows with more than 50 % entries missing;
# mean imputation used for these rows

imputed_data <- as.data.frame(t(imputed$data))
identical(as.character(rownames(imputed_data)), as.character(IDs$pseudo_ind_id))
imputed_data <- cbind(IDs, imputed_data)
write.csv(imputed_data, 'knn_imputed_olink_proteins.csv', row.names = F)
print('done')

#####################################################################
 # Run test

# Calculate the number of missing values per column
missing_values <- sapply(imputed_data, function(x) sum(is.na(x)))

# Filter out columns with no missing values
missing_values <- missing_values[missing_values > 0]

# Check if there are any missing values at all
if (length(missing_values) > 0) {
  print("There are missing values in the dataframe.")
  print("Missing values per protein/column:")
  print(missing_values)
} else {
  print("There are no missing values in the dataframe.")
}

# [1] "There are no missing values in the dataframe."

#####################################################################
# STEP 2: scaling

## PRE-SCALE Check histograms of non-scaled and scaled protein data
hist(imputed_data$GZMA)

# Rank-based inverse normalisation and scaling of each protein
for(i in colnames(imputed_data)[c(2:1464)]){
  imputed_data[,i] <- qnorm((rank(imputed_data[,i], na.last='keep')-0.5)/sum(!is.na(imputed_data[,i])))
}

## Scale protein data
imputed_data[,2:1464] <- apply(imputed_data[,2:1464], 2, scale)
mean(imputed_data[,200], na.rm = T)
sd(imputed_data[,200], na.rm = T)
write.csv(imputed_data, 'knn_imputed_olink_proteins_normalised_scaled.csv', row.names = F)

## POST-SCALE Check histograms of non-scaled and scaled protein data
hist(imputed_data$GZMA)

#####################################################################
# STEP 3: PCA


# Isolate just the protein columns of interest for PCA
prot <- imputed_data[c(2:1464)]
# Try prcomp
library(factoextra)

# Time approx 40s:
res.pca <- prcomp(prot, scale = TRUE)
names(res.pca)
# [1] "sdev"     "rotation" "center"   "scale"    "x"
loadings <- res.pca$rotation # this provides the loadings
loadings[1:5,1:4]

# the matrix x has the principal component score vectors
dim(res.pca$x)
# [1] 4644 1463

# Compute variance explained by componenets using standard deviations
std_dev <- res.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
head(prop_varex)

#[1] 0.18301069 0.08077721 0.03479630 0.02332409 0.01801900 0.01735832

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val
write.csv(eig.val, "eig_values_proteins_PCA.csv")

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation

### PLOT EIGENVALUES AND CUMULATIVE VARIANCE WITH CORR PLOT
var <- eig.val
var$num <- 1:1463
names(var)[1] <- "mes"
# var[1,1] <- 300
var$col <- ifelse(var$mes >= 1,  "#CF9FFF","lightgrey")
dim(var[which(var$mes >= 1),]) # 159 eigenvalues greater than or equal to 1

pdf("proteins_eigen_vals.pdf")
ggplot(var, 
       aes(num, mes, col = col)
       ) +
  geom_point(size = 1, 
             alpha = 0.5) + 
  labs(title = "Principal components analyses for the 1,463 protein analytes", 
       subtitle = "Eigenvalues for PCs",
       x = "PCs", 
       y = "Eigenvalue") +
  theme_classic() +
  theme(plot.title = element_text(size=11, face= "bold"),
        plot.subtitle = element_text(size=10),
        axis.text.y = element_text(size = 10),
        #color = Panel),
        axis.text.x = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.title.x = element_text(size = 10,face="bold"),
        legend.position = "none") +
 # xlab("Principal component") + 
 # ylab("Eigenvalue") + 
  geom_hline(yintercept=1, 
             color = "#D70040", 
             size=0.5,
             linetype='dotted') +
  scale_color_manual(values=c("#CF9FFF", "lightgrey"))
  
dev.off()

names(var)[3] <- "cum"

pdf("proteins_cum_var.pdf")

ggplot(var, aes(num, cum)) +
  geom_bar(stat = "identity", 
           fill = "#55bcc2",
           alpha = 0.6) + 
  labs(title = "Principal components analyses for the 1,463 protein analytes", 
       subtitle = "Cumulative proportion",
       x = "PCs", 
       y = "Cumulative proportion") +
  theme_classic() +
  theme(plot.title = element_text(size=11, face= "bold"),
        plot.subtitle = element_text(size=10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.title.x = element_text(size = 10,face="bold"),
        legend.position = "none") +
    geom_hline(yintercept=80, 
               color = "#D70040", 
               size=0.5,
               linetype='dotted') 
  #xlab("Principal component") + ylab("Cumulative proportion")
  
dev.off()

# ggcorrplot
library(ggcorrplot)
corr <- cor(prot)
pdf("proteins_corr_plot.pdf", width = 50, height = 50)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           outline.col = "white") + theme(text = element_text(size = 0.2))
dev.off()


devtools::session_info()
