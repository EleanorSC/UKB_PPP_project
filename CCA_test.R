## ---------------------------
##
## Script Purpose: Initial testing of CCA for analysis
##
## -----------------------------
##  General notes for neuroimaging data [structural MRI]
## -----------------------------
##
##
## -----------------------------
setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics")

# Install and load necessary packages

install.packages("tidyverse")
install.packages("magrittr")
install.packages("tidyr")
install.packages("dplyr")

library(tidyverse)
library(magrittr)
library(tidyr)
library(dplyr)

## -----------------------------

if (!requireNamespace("styler", quietly = TRUE)) install.packages("styler", dependencies = TRUE)
if (!requireNamespace("lintr", quietly = TRUE)) install.packages("lintr", dependencies = TRUE)

library(styler)
library(lintr)

## -----------------------------

# Style and lint the script
styler::style_file("CCA_test.R")
lintr::lint("CCA_test.R")

## MASTER_DATA is created using the UKB_PPP_initial_discovery.R script
## First let's model this on just the first 100 participants
# test_sample_GM <- MASTER_DATA %>% slice(1:100) %>% select(GMicv)


# Let's select all FA first:

#FA_data <- MASTER_DATA[, 
#                       c(1,
#                           grep("_FA", names(MASTER_DATA),
#                                value = TRUE)
#                           ), 
#                       drop = FALSE]
#
FA_data <- MASTER_DATA %>%
  select(1, matches("_FA"))

proteomics_data <- MASTER_DATA[,
                               c(1,
                                 164:1626), drop = FALSE]

# Check for missing data
proteomics_data_no_NAs <- proteomics_data[complete.cases(proteomics_data), ]
FA_data_no_NAs <- FA_data[complete.cases(FA_data), ]

# Ok, there are a BUNCH of NAs in proteomic data, our n reduces from n = 4644 to n=172

# Impute NAs with the mean of each column
proteomic_data_imputed <- proteomics_data

# Impute NAs with the mean of each column
for (col in colnames(proteomic_data_imputed)) {
  proteomic_data_imputed[[col]][is.na(proteomic_data_imputed[[col]])] <- mean(proteomic_data_imputed[[col]], na.rm = TRUE)
}

# Check for missing data
proteomic_data_imputed_no_NAs <- proteomic_data_imputed[complete.cases(proteomic_data_imputed), ]

# Perform an inner join based on the 'eid' column
merged_data <- inner_join(proteomic_data_imputed_no_NAs, 
                          FA_data_no_NAs,
                          by = "eid")

### 4644 obs of 1491 variables! 
#skimr::skim(merged_data)

# 1.1 Standardization
#merged_data <- merged_data[complete.cases(merged_data), ]

merged_data <- merged_data %>%
  mutate_all(as.numeric)

FA_data_2 <- merged_data %>%
  select(1, matches("_FA"))

proteomic_data_2 <- merged_data %>%
  select(1, -matches("_FA"))


# Remove rows with missing values
proteomic_data_3 <- na.omit(proteomic_data_2)
FA_data_3 <- na.omit(FA_data_2)

#Scale data
proteomic_data_3 <- scale(proteomic_data_2)
FA_data_3 <- scale(FA_data_2)

proteomic_data_4 <- as.data.frame(proteomic_data_3)
FA_data_4 <- as.data.frame(FA_data_3)

glimpse(FA_data_4)

proteomic_data_5 <- as.numeric(proteomic_data_4)
FA_data_5 <- as.numeric(FA_data_4)

glimpse(FA_data_4)

# 2.1 Apply CCA
cca_result <- cancor(
  proteomic_data_4,
  FA_data_4
)

any(is.nan(proteomic_data_4))
any(is.nan(FA_data_4))

numeric_proteomic_data <- as.numeric(proteomic_data_4[[1]])
numeric_FA_data <- as.numeric(FA_data_4[[1]])


numeric_FA_data <- as.numeric(as.data.frame((FA_data_4)))
numeric_FA_data <- as.numeric(unlist(FA_data_4))

# Check for NaN in the numeric data
any(is.nan(numeric_proteomic_data))
# Check for NaN in the numeric data
any(is.nan(numeric_FA_data))


FA_data_TEST <- merged_data %>%
  select(matches("_FA"))

proteomic_data_TEST <- merged_data %>%
  select(-matches("_FA"))

glimpse(FA_data_TEST)
numeric_FA_data <- as.numeric(unlist(FA_data_TEST))

numeric_proteomic_data <- as.numeric(unlist(proteomic_data_TEST))

### Merged data issue, uneven number of rows

# 2.1 Apply CCA
cca_result <- cancor(
  numeric_proteomic_data,
  numeric_FA_data
)


# 3.1 Canonical Correlation Coefficients
canonical_correlation_coefficients <- cor(cca_result$x, cca_result$y)

# 3.2 Canonical Variables
canonical_variables_proteomic <- cca_result$xcoef
canonical_variables_neuroimaging <- cca_result$ycoef

# 4.1 Biplot
biplot(cca_result)




subset_size <- 100  # Adjust as needed
cancor_result <- cancor(numeric_proteomic_data[1:subset_size, ], numeric_FA_data[1:subset_size, ])

## This takes 20:25 

# Some NAs introduced by coercian, check data:
# Check unique values in the column
unique_values <- unique(FA_data$lATR_FA)
skimr::skim(unique_values) # one value is missing, so we will remove

FA_data <- FA_data[complete.cases(FA_data), ]
skimr::skim(FA_data)

proteomic_data_standardized <- scale(proteomics_data)
FA_data_standardized <- scale(FA_data)

# Perform an inner join based on the 'eid' column
merged_data <- inner_join(as.data.frame(proteomic_data_standardized), 
                          as.data.frame(FA_data_standardized),
                          by = "eid")

# 2.1 Apply CCA
cca_result <- cancor(
  proteomic_data_standardized,
  FA_data_standardized
)

# unequal number of rows....



# 2.2 Statistical Significance (Optional)
# Conduct permutation tests or other significance tests
# Check documentation for permutation tests in the 'cancor' function or use external packages

# 3.1 Canonical Correlation Coefficients
canonical_correlation_coefficients <- cor(cca_result$x, cca_result$y)

# 3.2 Canonical Variables
canonical_variables_proteomic <- cca_result$x
canonical_variables_neuroimaging <- cca_result$y

# 4.1 Biplot
biplot(cca_result)
