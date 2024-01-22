## ---------------------------
##
## Script Purpose: Initial testing of PLS for analysis
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

## ----------------------------- #### PLS attempt

#write.csv(X, "proteomic_subsample.csv")
#write.csv(Y, "FA_subsample.csv")

X <- read.csv("proteomic_subsample.csv")
Y <- read.csv("FA_subsample.csv")

## -----------------------------
install.packages("pls")
library(pls)
## -----------------------------

X <- as.matrix(proteomic_data_TEST)  # the PA tests
Y <- as.matrix(FA_data_TEST) 

pls_result <- pls::mvr(X ~ Y, ncomp = min(nrow(Y), ncol(Y)) - 1, method = "kernelpls")

plot(pls_result)

validationplot(pls_result)
