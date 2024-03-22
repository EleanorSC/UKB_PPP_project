

## ---------------------------
##
## Script Purpose: Another PLS attempt: 
##
## -----------------------------
##  This R script runs a PLS on the data as a benchmark test. It does the following:
##.  1. 
##.  2. 
##.  3. 
##.  4. 
##.  5. 
##
## This script assumes you have loaded your data into two data frames: protein_data and gm_volume
## X (protein_data): A matrix of dimensions 4644 (individuals) x 1463 (proteins)
## Y (gm_volume): A vector or single-column data frame of grey matter volume for 4644 individuals
## -----------------------------
##
## -----------------------------

setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics/UKB_PPP_Rscripts")

## -------PACKAGES---------

install.packages("pls")

## -------ANALYSIS---------

#Step 1: Preparing the Data

# Load necessary libraries
library(pls)
#library(plsdepot)  # For alternative PLS methods if needed



## -------LOAD DATASETS---------
MASTER_DATA <- read.csv("MASTER_DATA.csv") #contains all main data
# Set 'eid' column as row names for the 'wm' dataset
rownames(MASTER_DATA) <- MASTER_DATA$eid

# use imputed_data or download from the .csv file saved (approx 30s)
imputed_data <- read.csv("knn_imputed_olink_proteins_normalised_scaled.csv")


## ------- ASSIGNING ---------
# Define GM + eid columns
gm_eid <- MASTER_DATA[, c("eid", "GMicv")]

## ------- SCALE ---------

scaled_gm <- scale(gm_eid[c(2)])

# Standardizing the data (Z-score normalization)
X <- as.matrix(imputed_data[c(2:1464)])
Y <- as.matrix(scaled_gm)


# Step 2: Running PLS Analysis

# Running PLS
pls_result <- plsr(Y ~ X, ncomp = 10, method = "kernelpls")

# View summary
summary(pls_result)

# Data: 	X dimension: 4644 1463 
# Y dimension: 4644 1
# Fit method: kernelpls
# Number of components considered: 10
# TRAINING: % variance explained
# 1 comps  2 comps  3 comps  4 comps  5 comps  6 comps  7 comps  8 comps  9 comps  10 comps
# X    12.38    24.05    27.99    30.60    32.19    33.35    34.82    35.84    36.57     37.51
# Y    13.94    23.79    38.64    44.25    48.30    52.19    54.63    57.02    58.99     60.30


#Step 3: Permutation Test (Simplified Spin Test)
#This part simplifies the spatial autocorrelation-preserving permutation test into a basic permutation test due to complexity.


# Simplified permutation test for significance of PLS components
# Assuming a basic permutation without spatial considerations for demonstration
# NOTE this takes ~ 2 mins to run
# The higher the no. permutations, the higher the value of p?

set.seed(123)
#nperm <- 1000  # Number of permutations- 1000 will take a long time
nperm <- 10 #  - 200 takes approx 2 minutes
perm_stats <- replicate(nperm, {
  perm_Y <- sample(Y)  # Permuting Y to break the association
  perm_pls <- plsr(perm_Y ~ X, ncomp = 10, method = "kernelpls")
  explained_variance <- explvar(perm_pls)  # Explained variance for all components
  explained_variance[1]  # Extract explained variance for the first component
})

# Extract the real explained variance for the first component
real_explvar <- explvar(pls_result)[1]

# Compare the real explained variance to the distribution of permuted explained variances
p_value <- mean(perm_stats >= real_explvar)
p_value

# Interpretation:
#"In assessing the significance of the explained variance by the first component of our PLSR model relating protein data to grey matter volume, a permutation test (n=20 permutations) was conducted. The p-value obtained from this test was 0.4, indicating that the explained variance observed could not be distinguished from chance (p > 0.05). Therefore, we find insufficient evidence to assert that the relationship captured by the first component is statistically significant, suggesting that further investigation is required to elucidate the nature of the association between protein profiles and grey matter volume.

#Step 4: Cross-validation
#Cross-validation can be performed using the validationplot function or by manually splitting the data.
# Basic cross-validation visualization
validationplot(pls_result, val.type = "MSEP")

#Step 5: Visualizing Results
#Visualization can vary depending on what aspect of the PLS results you are interested in. 
#Here's how to visualize the loadings:

# Visualizing loadings for the first component
loadings <- coef(pls_result, ncomp = 1)
barplot(loadings)