## ---------------------------
##
## Script Purpose: Initial overview of data for power analyses
##
## -----------------------------
##  General notes for neuroimaging data [structural MRI]
## -----------------------------
##
##
## -----------------------------
##  General notes for UKB proteomic data
## -----------------------------
##          Field ID	Description
##          30900	Number of proteins measured
##          30901	Plate used for sample run
##          30902	Well used for sample run
## ---------------------------

setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics")


## ----------------------------#
# install relevant packages

install.packages("tidyverse")
install.packages("magrittr")
install.packages("tidyr")

if (!requireNamespace("styler", quietly = TRUE)) install.packages("styler", dependencies = TRUE)
# if (!requireNamespace("lintr", quietly = TRUE)) install.packages("lintr", dependencies = TRUE)

library(tidyverse)
library(magrittr)
library(tidyr)
library(styler)

## ----------------------------#
# use a linter

styler::style_file("UKB_PPP_inital_discovery.R")



# find out which proteins are in the Olink inflammatory panel, from Supplement of 'Genetic regulation of the human plasma proteome in 54306 UK Biobank participants'

ST3_UKB_proteins <- read.csv("ST3_UKB_proteins.csv")

ST3_UKB_inflam_proteins <- ST3_UKB_proteins %>% filter(Protein.panel == "Inflammation")

## ----------------------------#
# find out how many participants have both proteomic and connectome data
# n = 37,283 participants with connectomes - first visit (first imaging visit)
UKB_connectome_IDs <- read.csv("subjects_37284.csv")

# n = 52,705 with proteomic at baseline
UKB_proteomics_IDs <- read.table("UKBB_proteomics_baseline_IDs", comment = "", header = TRUE)

# rename the ID columns
names(UKB_connectome_IDs)[names(UKB_connectome_IDs) == "X1000019"] <- "eid"
names(UKB_proteomics_IDs)[names(UKB_proteomics_IDs) == "eid"] <- "eid"

# Merge proteomics ID with connectome ID
# n = 4,869
UKB_connectomes_and_proteomics <- merge(UKB_connectome_IDs,
  UKB_proteomics_IDs,
  by = "eid"
)

UKBiobank_imaging <- readxl::read_excel("UKBiobank_imaging.xlsx")

# n = 7,080 neuroimaging data + proteomics
UKB_neuroimaging_and_proteomics <- merge(UKBiobank_imaging,
  UKB_proteomics_IDs,
  by = "eid"
)

# n = 7,080 neuroimaging data + proteomics
UKB_neuroimaging_connectome_and_proteomics <- merge(UKB_neuroimaging_and_proteomics,
  UKB_connectome_IDs,
  by = "eid"
)

# NOTE [TAKES 3 minutes to run]
# read in proteomics data
d1 <- read.table("olink_data.txt", header = TRUE)

# keep ins_index=0 (baseline data)
d2 <- d1[d1$ins_index == "0", ]

# remove ins_index
d3 <- subset(d2, select = -c(ins_index))

#######

# convert long format to wide format [TAKES A WHILE TO RUN ~1 hour]
# Record the start time
start_time <- Sys.time()

UKB_proteomic_data <- reshape(d3,
  idvar = "eid",
  timevar = "protein_id",
  direction = "wide"
)

# Record the end time
end_time <- Sys.time()

# Calculate the elapsed time
elapsed_time <- end_time - start_time

# Print a warning or message
warning(paste("The code took approximately", round(as.numeric(elapsed_time), 2), "minutes to run."))

#### The code took approximately 40.1 minutes to run.

# ID information - a meta data file containing age at baseline, sex and the Field IDs
meta_data <- read.table("UKB_PROTEOMICS_FOR_SARAH_01AUGUST2023_GD.txt", comment = "", header = TRUE)

# Olink Data-Coding 143
# Name: UniProt meaning for OLINK Protein ID
# Description: Relates the integer Protein ID presented in the OLINK dataset to the UniProt meaning.
# This is a flat (unstructured) list which uses integers to represent categories or special values.

ID_codes <- readr::read_tsv("coding143.tsv")

ID_codes <- ID_codes %>% tidyr::separate(meaning, c("Protein", "Protein_full_name"), ";")

# Keep a copy of original UKB_proteomic_data as reshape took ages to run
UKB_proteomic_data_raw <- UKB_proteomic_data

names(UKB_proteomic_data_raw) <- sub("^result.", "", names(UKB_proteomic_data_raw))
names(UKB_proteomic_data_raw) <- ID_codes$Protein[match(names(UKB_proteomic_data_raw), ID_codes$coding)]
names(UKB_proteomic_data_raw)[1] <- "eid"

# n = 4,644 connectome + neuroimaging data + proteomics
# [UKB_neuroimaging_connectome_and_proteomics]

MASTER_DATA <- merge(UKB_neuroimaging_connectome_and_proteomics,
  UKB_proteomic_data_raw,
  by = "eid"
)


write.csv(MASTER_DATA, "MASTER_DATA.csv") # 66.4MB





skimr::skim(MASTER_DATA$GMicv)


MASTER_DATA$GMicv <- as.numeric(MASTER_DATA$GMicv)

## First let's model this on just the first 100 participants
# test_sample_GM <- MASTER_DATA %>% slice(1:100) %>% select(GMicv)


GM_proteomics <- MASTER_DATA[, c(18, 164:1626)]

library(caret)


# Split the data into training and test set
set.seed(123)
training.samples <- GM_proteomics$GMicv %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- GM_proteomics[training.samples, ]
test.data <- GM_proteomics[-training.samples, ]


# Build the model on training set
set.seed(123)
model <- train(
  GMicv ~ .,
  data = train.data, method = "pls",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Plot model RMSE vs different values of components
plot(model)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model$bestTune

# Summarize the final model
summary(model$finalModel)

# Make predictions
predictions <- model %>% predict(test.data)
# Model performance metrics
data.frame(
  RMSE = caret::RMSE(predictions, test.data$medv),
  Rsquare = caret::R2(predictions, test.data$medv)
)


## make this example reproducible
# set.seed(1)
#
## fit PCR model
# model <- plsr(GMicv~IL6+IL10+ICAM1+ICAM5+VEGFA, data=MASTER_DATA, scale=TRUE, validation="CV")
#
## view summary of model fitting
## Once we’ve fit the model, we need to determine the number of PLS components worth keeping.
## The way to do so is by looking at the test root mean squared error (test RMSE) calculated by the k-fold cross-validation:
# summary(model)
#
# validationplot(model)
# validationplot(model, val.type="MSEP")

### TSBE
# Replace NA values with the mean of each protein
GM_proteomics_impute <- apply(GM_proteomics[, 2:1464], 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

GM_proteomics_noNAs <- na.omit(GM_proteomics_impute)

GM_proteomics_noNAs <- as.data.frame(GM_proteomics_noNAs)
# Separate grey matter volume and protein data
grey_matter_volume <- GM_proteomics[, 1]
protein_data <- GM_proteomics_noNAs[, 1:1463] # Assuming protein data is in columns 2 to 1001


# Apply t-SNE to reduce the dimensionality of protein data
set.seed(123) # for reproducibility
tsne <- tsne(protein_data, initial_dims = 2)

tsne <- data.frame(tsne)
pdb <- cbind(tsne, grey_matter_volume)

ggplot(pdb, aes(x = X1, y = X2, color = grey_matter_volume)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") + # Adjust color scale as needed
  labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2", color = "Grey Matter Volume") +
  ggtitle("t-SNE Plot of Protein Levels vs. Grey Matter Volume")
