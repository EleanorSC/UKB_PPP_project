setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics")

## ----------------------------# 
# install relevant packages

install.packages("tidyverse")
install.packages("magrittr")

## ----------------------------#
# find out which proteins are in the Olink inflammatory panel, from Supplement of 'Genetic regulation of the human plasma proteome in 54306 UK Biobank participants'

ST3_UKB_proteins <-read.csv("ST3_UKB_proteins.csv")

ST3_UKB_inflam_proteins <- ST3_UKB_proteins %>% filter(Protein.panel == "Inflammation")

## ----------------------------#
# find out how many participants have both proteomic and connectome data
# n = 37,283 participants with connectomes - first visit (first imaging visit)
UKB_connectome_IDs <- read.csv("subjects_37284.csv")

# n = 52,705 with proteomic at baseline
UKB_proteomics_IDs <- read.table("UKBB_proteomics_baseline_IDs", comment="", header=TRUE)

# rename the ID columns 
names(UKB_connectome_IDs)[names(UKB_connectome_IDs) == 'X1000019'] <- 'eid'
names(UKB_proteomics_IDs)[names(UKB_proteomics_IDs) == 'eid'] <- 'eid'

# Merge proteomics ID with connectome ID
# n = 4,869
UKB_connectomes_and_proteomics <- merge(UKB_connectome_IDs, 
                       UKB_proteomics_IDs, 
                       by = "eid")

#UKBiobank_imaging <- readxl::read_excel("UKBiobank_imaging.xlsx")

# n = 7,080 neuroimaging data + proteomics
UKB_neuroimaging_and_proteomics <- merge(UKBiobank_imaging, 
                                        UKB_proteomics_IDs, 
                                        by = "eid")


