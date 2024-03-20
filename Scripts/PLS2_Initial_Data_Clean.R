## ---------------------------
##
## Script Purpose: Initial data clean of PLS for analysis: 
##
## -----------------------------
##  This R script tidies and cleans the proteomic and neuroimaging data in preparation for PLS
##
##
##
## -----------------------------
##
## -----------------------------

setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics/UKB_PPP_Rscripts")


## -------PACKAGES---------

install.packages("tidyverse")
install.packages("magrittr")

library(tidyverse)
library(magrittr)

####### 1. using pls package #######

## -------LOAD DATASETS---------
MASTER_DATA <- read.csv("MASTER_DATA.csv") #contains all main data
ST3_UKB_proteins <- read.csv("ST3_UKB_proteins.csv") #contains data about proteomics data inc. protein panels

# Generate a list of all the proteins of relevance, e.g.
ST3_UKB_inflam_proteins <- ST3_UKB_proteins %>% filter(Protein.panel == "Inflammation")

# Assuming 'data_frame' is your data frame and 'proteins' is the column of interest
inflam_proteins_list <- as.list(ST3_UKB_inflam_proteins$Assay.Target)
inflam_proteins_vector <- as.character(unique(ST3_UKB_inflam_proteins$Assay.Target))


## -------TIDY DATASETS---------
# Remove unnecessary columns like 'X' column from the data frame
MASTER_DATA <- MASTER_DATA[, !(names(MASTER_DATA) %in% c("X"))]

# Set 'eid' column as row names for the 'wm' dataset
rownames(MASTER_DATA) <- MASTER_DATA$eid

# Examine the histograms of some of these variables
hist(MASTER_DATA$GMicv)

## ------- ASSIGNING ---------
# Define WM + eid columns
grey_matter_eid <- MASTER_DATA[, c("eid", "GMicv")]

# Undefined columns selected, so probably an issue with _ or . in protein data;
# Make sure that all - are converted into . or vice versa
inflam_proteins_eid <- MASTER_DATA[, c("eid", inflam_proteins_vector)]

# To fix, create a function that loops through these inflammatory proteins:
subset_inflammatory_proteins<- function(df, protein_list) {
  # Initialize a vector to hold the names of proteins not found in the dataframe
  not_found_proteins <- c()
  
  # Always include 'eid' in the subset
  subset_cols <- c("eid")
  
  # Loop through the protein list and check if each protein is a column in the dataframe
  for (protein in protein_list) {
    if (protein %in% colnames(df)) {
      subset_cols <- c(subset_cols, protein)
    } else {
      not_found_proteins <- c(not_found_proteins, protein)
    }
  }
  
  # Warn about proteins not found in the dataframe
  if (length(not_found_proteins) > 0) {
    warning("The following proteins are not in the dataframe: ", paste(not_found_proteins, collapse = ", "))
  }
  
  # Subset the dataframe to include only the columns for 'eid' and the found proteins
  subset_df <- df[, subset_cols, drop = FALSE]
  
  return(subset_df)
}

# Assuming 'main_df' is your large dataframe and 'inflam_proteins_vector' is your list of proteins
# Call the function with your dataframe and list of proteins
inflam_proteins_eid <- subset_inflammatory_proteins(MASTER_DATA, inflam_proteins_vector)

# Warning message:
#In subset_inflammatory_proteins(MASTER_DATA, inflam_proteins_vector) :
#  The following proteins are not in the dataframe: HLA-DRA, HLA-E

# Modify 'HLA-DRA' and 'HLA-E' to 'HLA.DRA' and 'HLA.E'
inflam_proteins_vector <- gsub("HLA-DRA", "HLA.DRA", inflam_proteins_vector)
inflam_proteins_vector <- gsub("HLA-E", "HLA.E", inflam_proteins_vector)

inflam_proteins_eid <- MASTER_DATA[, c("eid", inflam_proteins_vector)]

## ------- CLEANING ---------
# We have some missing values in the protein and wm data

# Step 1: Handle missing values
# Remove rows with any missing values in either data frame
# This will align the rows by removing any cases with missing data in either set
inflam_proteins <- MASTER_DATA[, c(inflam_proteins_vector)]
grey_matter <- MASTER_DATA[, c("GMicv")]

# 4,644 x 370
combined_data <- merge(grey_matter_eid, inflam_proteins_eid, by = "eid")

# 3,423 x 370
combined_data_clean <- na.omit(combined_data)

# Set 'eid' column as row names for the 'wm' dataset
rownames(combined_data_clean) <- combined_data_clean$eid

# Remove the 'eid' column from the dataset, now that it's set as row names
combined_data_clean <- combined_data_clean[, !(names(combined_data_clean) %in% "eid")]

# Step 2: Separate the cleaned data back into 'gm' and 'proteins'
gm <- combined_data_clean[1]  # Assuming the first column is grey matter
proteins <- combined_data_clean[2:ncol(combined_data_clean)]  # The rest are 'proteins'

# Step 3: Rank-based inverse normal transformation and scaling for 'proteins'
proteins_scaled <- as.data.frame(lapply(proteins, function(x) {
  ranks <- rank(x, na.last = "keep")
  trans <- qnorm((ranks - 0.5) / length(ranks))
  scale(trans)  # standardize to have mean 0 and SD 1
}))

hist(proteins$ACTN4)
hist(proteins_scaled$ACTN4)

# Standardize each variable in 'wm' to have mean 0 and standard deviation 1
gm_scaled <- as.data.frame(lapply(gm, scale))

hist(gm$GMicv)
hist(gm_scaled$GMicv)


## ------- MODEL PLS ---------
X <- as.matrix(proteins_scaled)
Y <- as.matrix(gm_scaled)

# takes approx 1-2 minutes to run
pls_loo <- plsr(Y ~ X, ncomp = 6, data = data.frame(X, Y), validation = "LOO")
summary(pls_loo)

validationplot(pls_loo, val.type = "MSEP")
plot(RMSEP(pls_loo))
plot(pls_loo, ncomp = 5, asp = 1, line = TRUE)
plot(pls_loo, plottype = "scores", comps = 1:5)
explvar(pls_loo)


# Extract loadings for the first two components
loadings_x <- coef(pls_loo, ncomp = 1:5)

# Convert to a data frame for easier manipulation
loadings_df <- as.data.frame(loadings_x)
names(loadings_df) <- c("Component1", "Component2", "Component3", "Component4", "Component5")
loadings_df$Protein <- rownames(loadings_df)


ggplot(top_loadings_df, aes(x = reorder(Protein, Component1), 
                            y = Component1)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  # Use actual color values from the Color column
  geom_text(aes(x = reorder(Protein, Component1), 
                y = Component1 + 0.02,  # Slightly offset text for visibility
                label = Protein), 
            hjust = 0.5, vjust = 0, check_overlap = TRUE, size = 3) +
  labs(title = "PLS Loadings for Top Contributor Proteins for Grey Matter", x = "Protein", y = "Loading") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +  # Flip coordinates for horizontal bars
  ylim(-0.001, 0.001)  # Set y-axis limits, which correspond to original x-axis limits
