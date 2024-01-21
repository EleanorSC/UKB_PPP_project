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

install.packages("lintr")
library(lintr)

# To check the entire script, use the following command:
# lint()
# lint(filename = "CCA_test.R")
# lint("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics/CCA_test.R", fix = TRUE)

# Install and load necessary packages
if (!requireNamespace("styler", quietly = TRUE)) install.packages("styler", dependencies = TRUE)
# if (!requireNamespace("lintr", quietly = TRUE)) install.packages("lintr", dependencies = TRUE)

library(styler)
# library(lintr)

# Style and lint the script
styler::style_file("CCA_test.R")
lintr::lint("CCA_test.R")

## MASTER_DATA is created using the UKB_PPP_initial_discovery.R script
## First let's model this on just the first 100 participants
# test_sample_GM <- MASTER_DATA %>% slice(1:100) %>% select(GMicv)


# Let's select all FA first:

FA_data <- MASTER_DATA[, c(
  "eid",
  grep("_FA", names(MASTER_DATA),
    value = TRUE
  )
), drop = FALSE]


GM_proteomics <- MASTER_DATA[, c(18, 164:1626)]


# 2.1 Apply CCA
cca_result <- cancor(
  proteomic_data_standardized,
  neuroimaging_data_standardized
)

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
