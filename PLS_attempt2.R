

## ---------------------------
##
## Script Purpose: Another PLS attempt: 
##
## -----------------------------
##  This R script runs a PLS on the data as a benchmark test. It does the following:
##.  1. Loads matrices of interet: predictors (proteins) and outcomes (neuroimaging); scales neuroimaging outcome variable (e.g. grey matter volume)
##.  2. Runs simple pls model: plsr(Y ~ X, ncomp = 10, method = "kernelpls")
##.  3. Runs a permutation test (Simplified Spin Test) to evaluate p-value (variance explained by PC1)
##.  4. Performs cross-validation; MSEP plot to determine how many components to extract
##.  5. Generates a loadings plot for PC1-PC5 
##
## This script assumes you have loaded your data into two data frames: protein_data and gm_volume
## X (protein_data): A matrix of dimensions 4644 (individuals) x 1463 (proteins)
## Y (gm_volume): A vector or single-column data frame of grey matter volume for 4644 individuals
## -----------------------------
##
## -----------------------------

setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics/UKB_PPP_Rscripts")

sessionInfo()

## -------PACKAGES---------

#install.packages("pls")
## -------ANALYSIS---------

#Step 1: Preparing the Data

# Load necessary libraries
library(pls)
library(tidyverse)
library(magrittr)

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
# Cross-validation can be performed using the validationplot function or by manually splitting the data.
# Basic cross-validation visualization
validationplot(pls_result, val.type = "MSEP")

# MSEP plot indicates that 3 components should be selected (elbow)

#Step 5: Visualizing Results
#Visualization can vary depending on what aspect of the PLS results you are interested in. 
#Here's how to visualize the loadings:

# Visualizing loadings for the first 5 components
loadings_x <- as.data.frame(coef(pls_result, ncomp = 1:5))
# Rename the column using dplyr's rename function
loadings_x <- loadings_x %>%
  rename(PC1_loadings = `Y.1 comps`,
         PC2_loadings = `Y.2 comps`,
         PC3_loadings = `Y.3 comps`,
         PC4_loadings = `Y.4 comps`,
         PC5_loadings = `Y.5 comps`
         )

# Convert row names to a column called 'proteins'
loadings_x <- loadings_x %>% 
  rownames_to_column(var = "proteins")

#Not sure why colour isn't changing?
ggplot(loadings_x, 
       aes(x = reorder(proteins, PC2_loadings), 
           y = PC2_loadings,
          # fill = "#55bcc2", 
           fill = "#CF9FFF",
           alpha = 0.6)
       ) +
  geom_bar(stat = "identity") +
  labs(title = "PLS Loadings for Grey Matter", 
       subtitle = "Dimension 2 (13.9%)",
       x = "Protein", y = "Loadings for Second Component") +
  theme_classic() +
  theme(plot.title = element_text(size=11, face= "bold"),
        plot.subtitle = element_text(size=10),
        axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 2),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.title.x = element_text(size = 10,face="bold"),
        legend.position = "none") +
  coord_flip()

#########

# Extracting loadings for the first component
loadings <- as.data.frame(loadings(pls_result)[, 1])

# Load the dplyr package
library(dplyr)

# Rename the column using dplyr's rename function
loadings <- loadings %>%
  rename(loadings = `loadings(pls_result)[, 1]`)

# Convert row names to a column called 'proteins'
loadings <- loadings %>% 
  rownames_to_column(var = "proteins")

ggplot(loadings, 
       aes(x = reorder(proteins, loadings), 
                       y = loadings,
                       fill = "#D70040", 
                       alpha = 0.6)
       ) +
  geom_bar(stat = "identity") +
  labs(title = "PLS Loadings for Grey Matter", 
       subtitle = "Dimension 1 (13.9%)",
       x = "Protein", y = "Loadings for First Component") +
  theme_classic() +
  theme(plot.title = element_text(size=11, face= "bold"),
        plot.subtitle = element_text(size=10),
        axis.text.y = element_text(size = 2),
        #color = Panel),
        axis.text.x = element_text(size = 2),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.title.x = element_text(size = 10,face="bold"),
        legend.position = "none") +
  coord_flip()

# Since loadings are typically in a complex structure, ensure it's a numeric vector before plotting
loadings_vector <- as.numeric(loadings)

# Create a bar plot of the loadings for the first component
barplot(loadings_vector)



## -------PLS2 but for FA across 27 TRACTS---------

########################## Run PLS2 but on tract FA data

## -------LOAD DATASETS---------
MASTER_DATA <- read.csv("MASTER_DATA.csv") #contains all main data
# Set 'eid' column as row names for the 'wm' dataset
rownames(MASTER_DATA) <- MASTER_DATA$eid

# use imputed_data or download from the .csv file saved (approx 30s)
imputed_data <- read.csv("knn_imputed_olink_proteins_normalised_scaled.csv")


## ------- ASSIGNING ---------
# Define white matter tracts columns
wm_cols <- c("FMaj_FA", "FMin_FA", "lAR_FA", "lATR_FA", "lCingG_FA", "lCingPH_FA", 
             "lCST_FA", "lIFOF_FA", "lILF_FA", "lML_FA", "lPTR_FA", "lSLF_FA", 
             "lSTR_FA", "lUnc_FA", "MCP_FA", "rAR_FA", "rATR_FA", "rCingG_FA", 
             "rCingPH_FA", "rCST_FA", "rIFOF_FA", "rILF_FA", "rML_FA", "rPTR_FA", 
             "rSLF_FA", "rSTR_FA", "rUnc_FA")

# Extract white matter tracts data
wm <- MASTER_DATA[, c("eid", wm_cols)]

# Extract proteins data by selecting columns not in wm_cols, also ignore the column eid
#proteins <- df_neuroimaging[, !(names(df_neuroimaging) %in% wm_cols)]

## ------- SCALE ---------

scaled_FA <- scale(wm[c(2:28)])

# Standardizing the data (Z-score normalization)
X <- as.matrix(imputed_data[c(2:1464)])
Y <- as.matrix(scaled_FA)

# Step 2: Running PLS Analysis

#library(pls)

# Running PLS
pls_result <- plsr(Y ~ X, ncomp = 10, method = "kernelpls")

# View summary
summary(pls_result)

# some plots:

# coefplot from the pls package creates plots

if (FALSE) {
  coefplot(pls_result, ncomp = 1:6)
  plot(pls_result, plottype = "coefficients", ncomp = 1:4) # Equivalent to the previous
  ## Plot with legend:
  coefplot(pls_result, 
           ncom = 1:4, 
           legendpos = "bottomright") +
}

# If we are to recreate these using ggplot we need to extract coefs from the pls_result object
# Creating a custom version of the coefplot from the pls package using ggplot2 in R:
#  (1) need to extract the regression coefficients from our PLS model object,
#   (2) then plot these coefficients using ggplot2
#   (3) finally, apply my specified theme settings.

# Step 1: Extract Regression Coefficients
# First, I need to extract the regression coefficients from the PLS model object. 
# I can do this using the coef() function from the pls package. 
# Note that if I want to consider a specific number of components, 
# I should specify this when extracting the coefficients.


# Assuming pls_model is your PLS model object

# Running PLS
pls_model <- plsr(Y ~ X, ncomp = 4, method = "kernelpls")

summary(pls_model)

# Step 1 :Extract coefficients for the specified number of components (e.g., for all components)
coefficients <- coef(pls_model, ncomp = pls_model$ncomp)

coefficients <- coef(pls_model)
# Convert to a dataframe for easier handling with ggplot

coef_df <- as.data.frame(coefficients)

# Keep information about proteins
coef_df$Protein <- rownames(coef_df)

# Step 2: Prepare Data for ggplot
#The dataframe needs to be in a long format for ggplot, so you might need to reshape it.

library(tidyr)

coef_df_long <- pivot_longer(coef_df, 
                             cols = -Protein, 
                             names_to = "Component", 
                             values_to = "Coefficient")

coef_df_long <- separate(coef_df_long, 
                         col = Component, 
                         into = c("WM_tract", "Component"), 
                         sep = "\\.", 
                         remove = FALSE) 

#### Nb. causes R to encounter a fatal error 

coef_df_long_FMaj_FA <- coef_df_long %>% filter(WM_tract == "FMaj_FA")

ggplot(coef_df_long_FMaj_FA, 
       aes(x = reorder(Protein, Coefficient),
           y = Coefficient, 
           group = Component,
           color = Component)
       ) +
  geom_line() + # or geom_point(), depending on what you're trying to visualize
  theme_classic() +
  coord_flip() +
  theme(plot.title = element_text(size=11, face= "bold"),
        plot.subtitle = element_text(size=10),
        axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 5, 
                                   #angle = 90, 
                                   hjust = 1), # Rotate x labels for clarity
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.title.x = element_text(size = 10, face="bold")) +
#  facet_grid(~WM_tract) +
  labs(x = "Protein", y = "Regression Coefficient", title = "PLS Model Coefficients")


# Step 3: Use ggplot
ggplot(coef_df_long, aes(x = Variable, y = Coefficient)) +
  geom_line() + # or geom_point(), depending on your preference
  theme_classic() +
  theme(plot.title = element_text(size=11, face= "bold"),
        plot.subtitle = element_text(size=10),
        axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 2, angle = 90, hjust = 1), # Rotate x labels for clarity
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.title.x = element_text(size = 10, face="bold")) +
  labs(x = "Variable", y = "Regression Coefficient", title = "PLS Model Coefficients")




## ------- FIGURES FROM PLS PACKAGE ---------

############
# X <- as.matrix(imputed_data[c(2:1464)])
# Y <- as.matrix(scaled_FA)
# pls_model <- plsr(Y ~ X, ncomp = 4, method = "kernelpls")


## ------- HOW MANY COMPONENTS? ---------

#Suggestions for the optimal number of components 

# Conveniently, the pcr function from the pls package has an implementation 
# for k-fold Cross Validation. We simply need to set validation = CV and segments = 20
# to perform 20-fold Cross Validation with PC regression. 
# If we don’t specify ncomp, pcr will select the maximum number of PCs that can be used for the CV.

# Note that our training data X consists of 4,000 + observations (rows).  
# If we perform 20-fold CV, that means we will split the data in 20 groups,
# so each group will consist of 4644/20 = 232 observations. 

# At each CV cycle, one group will be left out and the model will be trained on
# the remaining groups. This leaves us with 232 training observations for each CV
# cycle, so the maximal number of components that can be used in the linear regression
# is 231


# Taking 30+ minutes.. due to LOO Leave-One-Out (LOO) cross-validation, 
# which is computationally intensive because it fits the model 4,000+ times 
# (once for each observation left out). 
#This will considerably increase the runtime compared to no validation 
#or a simpler validation method like k-fold cross-validation.

start_time <- Sys.time()
yarn.pls <- plsr(Y ~ X, 
                 scale = TRUE,
                 ncomp = 8, 
                 validation = "CV")
end_time <- Sys.time()

# Calculate the time taken
time_taken <- end_time - start_time
print(time_taken)

# Plots incremental 
plot(yarn.pls, plottype = "validation")

optimal_ncomp <- selectNcomp(yarn.pls, 
                             method = "onesigma", 
                             plot = TRUE)

plot(yarn.pls, plottype = "validation")

selectNcomp(yarn.pls, 
            "onesigma", 
            plot = TRUE, 
            ylim = c(0, 3)
            )

selectNcomp(yarn.pls, "randomization", plot = TRUE)
selectNcomp(yarn.pls, "randomization", plot = TRUE,
            ncomp = 10, ylim = c(0, 3))

## ------- BIPLOTS ---------
biplot(pls_model, which = "x")
biplot(pls_model, which = "y")
biplot(pls_model, which = "scores")
biplot(pls_model, which = "loadings")

## ------- VALIDATION PLOTS ---------
#validationplot(pls_model, estimate = "all")
plot(RMSEP(pls_model, estimate = "all"))

## Plot R2:
plot(pls_model, "validation", val.type = "R2")

## Plot R2, with a legend:
plot(pls_model, 
     "validation", 
     val.type = "MSEP", 
     legendpos = "top") # R >= 2.1.0

## ------- PREDICTION PLOTS ---------
mod <- plsr(Y ~ X, ncomp = 6, validation = "CV")

plot(mod, ncomp = 1:6) # Equivalent to the previous

# ## Both cross-validated and test set predictions:
#  predplot(mod, 
#           ncomp = 4:6, 
#           which = c("validation", "test"),
#            newdata = yarn[!yarn$train,]
# 
#           
# mod.sens <- plsr(sensory ~ chemical, ncomp = 4, data = oliveoil)
# if (FALSE) plot(mod.sens, ncomp = 2:4) # Several responses gives several plots

