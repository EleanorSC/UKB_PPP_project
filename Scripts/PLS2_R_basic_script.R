## ---------------------------
##
## Script Purpose: Initial testing of PLS for analysis
##
## -----------------------------
##  This R script defines functions for performing Partial Least Squares Regression (PLS2), 
##  orthogonalized Covariance PLS2 (oCPLS2), and Kernel PLS2 (KPLS2), 
##  aimed at examining associations between predictors (e.g., proteomic data) 
##  and responses (e.g., neuroimaging data). 
## -----------------------------
##
##
## -----------------------------


## -----------------------------
##  In part 1, this R script uses the package pls and additionally
##  - loads a subset of data 
##  - tidies data
##  - draws initial plots 
## -----------------------------

setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics/UKB_PPP_Rscripts")


## -------PACKAGES---------

install.packages("tidyverse")
install.packages("magrittr")

library(tidyverse)
library(magrittr)


## -------LOAD DATASETS---------
df_p_data <- read.csv("ST3_UKB_proteins.csv") #contains protein panel information
main_df <-read.csv("merged_data.csv") #contains proteomics + neuroimaging


## -------TIDY DATASETS---------
# Remove unnecessary columns like 'X' column from the data frame
#df_neuroimaging <- main_df[, !(names(main_df) %in% c("eid"))]
df_neuroimaging <- main_df[, !(names(df_neuroimaging) %in% c("X"))]

# Display the first few rows of the data frame
head(df_neuroimaging, n = 10)

## ------- ASSIGNING ---------
# Define white matter tracts columns
wm_cols <- c("FMaj_FA", "FMin_FA", "lAR_FA", "lATR_FA", "lCingG_FA", "lCingPH_FA", 
             "lCST_FA", "lIFOF_FA", "lILF_FA", "lML_FA", "lPTR_FA", "lSLF_FA", 
             "lSTR_FA", "lUnc_FA", "MCP_FA", "rAR_FA", "rATR_FA", "rCingG_FA", 
             "rCingPH_FA", "rCST_FA", "rIFOF_FA", "rILF_FA", "rML_FA", "rPTR_FA", 
             "rSLF_FA", "rSTR_FA", "rUnc_FA")

# Extract white matter tracts data
wm <- df_neuroimaging[, c("eid", wm_cols)]

# Extract proteins data by selecting columns not in wm_cols, also ignore the column eid
#proteins <- df_neuroimaging[, !(names(df_neuroimaging) %in% wm_cols)]

# First, create a vector that combines 'eid' with the wm_cols
cols_to_exclude <- c("X", wm_cols)

# Then, select columns from df_neuroimaging that are not in cols_to_exclude
proteins <- df_neuroimaging[, !(names(df_neuroimaging) %in% cols_to_exclude)]

# Predefined colors vector in R
colors <- list(Inflammation = "tab:red", 
               Neurology = "tab:blue", 
               Cardiometabolic = "tab:orange", 
               Oncology = "tab:green")

# Initialize an empty list for the lookup
lookup <- list()

# Iterate over the column names of the proteins dataframe
for (index in seq_along(colnames(proteins))) {
  name <- colnames(proteins)[index]
  panel <- df_p_data[df_p_data$`Assay.Target` == name, ]$`Protein.panel`
  
  # Handle the case where panel might be empty or have multiple entries
  if (length(panel) == 0) {
    panel <- "" # Default value if no match is found
    color <- "tab:grey"
  } else {
    panel <- panel[1] # Take the first match if there are multiple
    color <- colors[[panel]] # Use the panel name to lookup the color
  }
  
  # Add the entry to the lookup list
  lookup[[index]] <- list(name = name, panel = panel, color = color)
}

# Optionally, to print or inspect the lookup list
print(lookup)

# Define the function in R
dictionary_head <- function(dictionary, num_items = 5) {
  items_to_display <- head(dictionary, num_items)
  for (index in seq_along(items_to_display)) {
    cat(paste0(names(items_to_display)[index], ": ", toString(items_to_display[[index]]), "\n"))
  }
}

# Assuming 'lookup' is a list created from previous steps
# Run the function to display the first num_items from lookup
dictionary_head(lookup)

# Test our dictionary :
lookup[2]

# Select 'Assay Target' and 'Protein panel' columns and then filter rows where 'Assay Target' is 'THBS4'
df2 <- subset(df_p_data, select = c(Assay.Target, Protein.panel))
matched_items <- df2[df2$`Assay.Target` == 'THBS4', ]

# Display the matched items
print(matched_items)

# Consider setting eid as index
# Merge wm and proteins


proteins_and_wmtracts <- merge(wm, proteins, by = "eid")


# Set 'eid' column as row names for the 'wm' dataset
rownames(wm) <- wm$eid

# Remove the 'eid' column from the dataset, now that it's set as row names
wm <- wm[, !(names(wm) %in% "eid")]

# Repeat the process for the 'proteins' dataset
rownames(proteins) <- proteins$eid
proteins <- proteins[, !(names(proteins) %in% "eid")]



## -----------------------------
##
## PART 1b: RUN INITIAL TEST PLS model
## -----------------------------

#str(proteins)
#pls_m1 <- plsr(wm ~ proteins, ncomp = 10, data = data.frame(X, Y), validation = "LOO")

X <- as.matrix(proteins)
Y <- as.matrix(wm)

# Fit the PLS model
pls_model <- plsr(Y ~ X, ncomp = 2, scale = TRUE)

# Summary of the model
summary(pls_model)

# Data: 	X dimension: 4644 1463 
# Y dimension: 4644 27
# Fit method: kernelpls
# Number of components considered: 2
# TRAINING: % variance explained
# 1 comps  2 comps
# X            0.3130    1.162
# FMaj_FA      5.8041    7.409
# FMin_FA     13.7581   17.154
# lAR_FA       5.9194    8.094
# lATR_FA     13.1910   16.682
# lCingG_FA    7.3206    9.083
# lCingPH_FA   2.5109    3.286
# lCST_FA      5.2031    6.513
# lIFOF_FA    14.9698   18.502
# lILF_FA     14.5338   17.920
# lML_FA       0.6607    1.020
# lPTR_FA      9.8447   12.256
# lSLF_FA     14.2753   17.352
# lSTR_FA      6.8537    8.654
# lUnc_FA      8.6862   10.849
# MCP_FA       2.4966    3.720
# rAR_FA       5.4761    7.903
# rATR_FA     12.5069   15.289
# rCingG_FA    7.3711    8.675
# rCingPH_FA   1.6682    2.326
# rCST_FA      4.9039    6.382
# rIFOF_FA    15.1197   18.690
# rILF_FA     14.9330   18.573
# rML_FA       0.9066    1.230
# rPTR_FA     10.8904   13.703
# rSLF_FA     13.3671   16.523
# rSTR_FA      6.4716    8.239
# rUnc_FA      9.0083   11.181


# Extract loadings of the model
# Plot the loadings for the predictors (X)
# Plotting predictor (X) loadings
plot(pls_model, comps = 1:2, 
    plottype = "loadings")

plot(pls_model, comps = 1:2, 
     plottype = "scores")

##############

# Extract loadings for the first two components
loadings_x <- coef(pls_model, ncomp = 1:2)

# Convert to a data frame for easier manipulation
loadings_df <- as.data.frame(loadings_x)
names(loadings_df) <- c("Component1", "Component2")
loadings_df$Protein <- rownames(loadings_df)

# Assume you have a threshold or specific criteria for top contributors
# For demonstration, let's select top contributors based on absolute values in Component 1
threshold <- quantile(abs(loadings_df$Component1), 0.9) # Top 10% as an example
top_contributors <- loadings_df[abs(loadings_df$Component1) > threshold,]


# Function to retrieve color from lookup
getColorForProtein <- function(proteinName, lookup) {
  # Default color if not found
  defaultColor <- "grey"
  
  # Search through lookup
  for(item in lookup) {
    if(item$name == proteinName) {
      return(item$color)
    }
  }
  
  return(defaultColor)
}

# Apply the function to assign colors
loadings_df$Color <- sapply(loadings_df$Protein, getColorForProtein, lookup=lookup)

# Replace 'tab:blue' etc. with qppropriate hex codes:
loadings_df$Color[loadings_df$Color == "tab:blue"] <- "#0000FF"
loadings_df$Color[loadings_df$Color == "tab:orange"] <- "#FFCC99"
loadings_df$Color[loadings_df$Color == "tab:green"] <- "#40E0D0"
loadings_df$Color[loadings_df$Color == "tab:red"] <-  "#8B0000"
loadings_df$Color[loadings_df$Color == "tab:grey"] <-  "#F0F0F0"
#loadings_df$Color

threshold <- quantile(abs(loadings_df$Component1), 0.9) # Adjust threshold as needed
loadings_df$TopContributor <- abs(loadings_df$Component1) > threshold


library(ggplot2)

ggplot(top_loadings_df, aes(x = reorder(Protein, Component1), 
                            y = Component1, 
                            fill = Color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  # Use actual color values from the Color column
  geom_text(aes(x = reorder(Protein, Component1), 
                y = Component1 + 0.02,  # Slightly offset text for visibility
                label = Protein), 
            hjust = 0.5, vjust = 0, check_overlap = TRUE, size = 3) +
  labs(title = "PLS Loadings for Top Contributor Proteins", x = "Protein", y = "Loading") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +  # Flip coordinates for horizontal bars
  ylim(-0.001, 0.001)  # Set y-axis limits, which correspond to original x-axis limits



## -----------------------------
##
## PART 2a: TEST OTHER R PACKAGE PLSREG2
## -----------------------------
install.packages("plsdepot")
library(plsdepot)

# N.b missing data to consider
PLS2_object<- plsreg2(X, Y, comps = 2, crosval = FALSE)




# X.scores
# The scatter plot shows how well the LVs from the protein data align or correlate with the LVs from the DTI data. 
# If there are clear patterns or clusters, it indicates a potential relationship between protein expression and DTI metrics.

# Extract scores
scores_x <- PLS2_object$x.scores
scores_y <- PLS2_object$y.scores

# Convert to a data frame
scores_df <- data.frame(Score_X1 = scores_x[,1], Score_Y1 = scores_y[,1], 
                        Score_X2 = scores_x[,2], Score_Y2 = scores_y[,2])

# Plotting
ggplot(scores_df, aes(x = Score_X1, y = Score_Y1)) +
  geom_point(aes(color = "Predictor Scores", alpha = 0.6)) +
  geom_point(aes(x = Score_X2, 
                 y = Score_Y2, 
                 color = "Response Scores", 
                 alpha = 0.01)) +
  scale_color_manual(values = c("Predictor Scores" = "blue", 
                                "Response Scores" = "red")) +
  theme_minimal() +
  labs(title = "PLS Model Scores for the First Two Components",
       x = "First Component", y = "Second Component") +
  theme(legend.title = element_blank())
#######

#LOADINGS
# Assuming you have a way to extract or calculate loadings for plotting

# Extract loadings for the first component (hypothetical example)
loadings_x1 <- PLS2_object$x.loads[,1]

# Convert to a data frame
loadings_df <- data.frame(Protein = names(loadings_x1), Loading = loadings_x1)

#Match up loadings_df to information regarding protein panel

# Function to retrieve color from lookup
getColorForProtein <- function(proteinName, lookup) {
  # Default color if not found
  defaultColor <- "grey"
  
  # Search through lookup
  for(item in lookup) {
    if(item$name == proteinName) {
      return(item$color)
    }
  }
  
  return(defaultColor)
}

# Apply the function to assign colors
loadings_df$Color <- sapply(loadings_df$Protein, getColorForProtein, lookup=lookup)


# Predefined colors vector in R
# colors <- list(Inflammation = "tab:red", 
#                Neurology = "tab:blue", 
#                Cardiometabolic = "tab:orange", 
#                Oncology = "tab:green")

# Replace 'tab:blue' etc. with original panels
loadings_df$Panel[loadings_df$Color == "tab:blue"] <- "neurology"
loadings_df$Panel[loadings_df$Color == "tab:orange"] <- "cardiometabolic"
loadings_df$Panel[loadings_df$Color == "tab:green"] <- "oncology"
loadings_df$Panel[loadings_df$Color == "tab:red"] <-  "inflammatory"
loadings_df$Panel[loadings_df$Color == "tab:grey"] <-  "misc"


loadings_df$Panel <- as.factor(loadings_df$Panel)

#n.b we have two instances where a protein is assigned "misc" so ensure to cross-check here

loadings_df %<>% filter(Panel == "neurology" |
                         Panel == "cardiometabolic" |
                         Panel == "inflammatory" |
                         Panel == "oncology")

# Plotting loadings
ggplot(loadings_df, aes(x = reorder(Protein, Loading), y = Loading, fill = Panel)) +
  geom_bar(stat = "identity",
           #fill =Panel, 
           alpha = 0.4) +
  coord_flip() +
  theme_classic() +
  facet_grid(~Panel) +
  labs(title = "PLS2 Loadings for Predictor Variables - Component 1",
       x = "Proteins", y = "Loading") +
  theme(plot.title = element_text(size=11, face= "bold"),
        axis.text.y = element_text(size = 1),
        axis.text.x = element_text(size = 5),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.title.x = element_text(size = 10,face="bold"),
        legend.position="bottom",
        legend.title = element_text(size=10, 
                                  face="bold"))

###### Plot smaller smample

#library(dplyr)

# Assuming loadings_df is your dataframe
# Select a random 10% of the rows
set.seed(42)  # For reproducibility
random_subset_df <- loadings_df %>% sample_n(size = floor(0.1 * nrow(loadings_df)))

# View the first few rows of the subset
head(random_subset_df)

# Plotting loadings
ggplot(random_subset_df, aes(x = reorder(Protein, Loading), y = Loading, fill = Panel)) +
  geom_bar(stat = "identity",
           #fill =Panel, 
           alpha = 0.4) +
  coord_flip() +
  theme_classic() +
  facet_grid(~Panel) +
  labs(title = "PLS2 Loadings for Predictor Variables - Component 1",
       x = "Proteins", y = "Loading") +
  theme(plot.title = element_text(size=11, face= "bold"),
        axis.text.y = element_text(size = 4),
                                   #color = Panel),
        axis.text.x = element_text(size = 5),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.title.x = element_text(size = 10,face="bold"),
        legend.position="bottom",
        legend.title = element_text(size=10, 
                                    face="bold"))





plot(PLS2_object[["x.scores"]]) # explanation
plot(PLS2_object[["y.scores"]]) # explanation
plot(PLS2_object[["x.loads"]]) # explanation

# plot variables (circle of correlations)
plot(PLS2_object, what="variables")

# plot observations (as points)
plot(PLS2_object, what="observations")

# plot observations with labels
plot(PLS2_object, what="observations", show.names=TRUE)


## -----------------------------
##
## PART 2: TRAN AND TEST SPLIT - CAN'T GET TO WORK YET
## -----------------------------

# In this section, we will do a PLSR on the  data to illustrate the use of pls. 
# We first divide the data set into train and test data sets: 

# Assuming 'proteins' and 'wm' are your dataframes
set.seed(42) # Ensure reproducibility

# Calculate the number of observations for the training set
# For example, to use 75% of the data for training
training_size <- floor(0.75 * nrow(proteins))

# Randomly sample indices for the training data
training_indices <- sample(seq_len(nrow(proteins)), size = training_size)

# Split the predictors into training and testing sets
X_train <- proteins[training_indices, ]
X_test <- proteins[-training_indices, ]

# Split the outcomes into training and testing sets
Y_train <- wm[training_indices, ]
Y_test <- wm[-training_indices, ]

# Now, X_train, Y_train are your training predictors and outcomes, 
# and X_test, Y_test are your testing predictors and outcomes

## -----------------------------
##
## PART 3: RUNNING PLS MODEL
## -----------------------------

# Fit the PLS model
# Assuming Y_train is a dataframe of response variables (white matter tracts)
# and X_train is a dataframe of predictor variables (proteins)

# If X_train or Y_train are lists or not in the expected format, convert them to data frames or matrices
X_train_df <- as.data.frame(X_train)
Y_train_df <- as.data.frame(Y_train)

pls_model_1<- plsr(Y_train_df ~ X_train_df, 
                   ncomp = 5, 
                   method = "kernelpls", 
                   scale = TRUE, 
                   validation = "CV")

# Convert to matrices if they're not already (this step is crucial for numerical analysis)
X_train_mat <- as.matrix(X_train_df)
Y_train_mat <- as.matrix(Y_train_df)

# Fit the PLS model using matrices
pls_model <- plsr(Y_train_mat ~ X_train_mat, ncomp = 5, method = "kernelpls", scale = TRUE, validation = "CV")



# View the summary of the model
summary(pls_model)

# Plotting the variance explained by each component
plot(pls_model, "variance")

# Plotting the predicted vs observed response for the first component
# This can be useful to assess model performance visually
plot(pls_model, "prediction", ncomp = 1)

# To assess model performance using cross-validation
validationplot(pls_model, val.type = "MSEP")


pls_model_1 <- plsr(octane ~ NIR, 
                    ncomp = 10, 
                    data = gasTrain, 
                    validation = "LOO") 

#This fits a model with 10 components, 
#and includes leave-one-out (LOO) cross-validated predictions [12]. 
#We can get an overview of the fit and validation results with the summary method:

## -----------------------------
##
## PART 1: Simulation
## -----------------------------









setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics")

## -----------------------------
##
## PART 1: Simulation
## -----------------------------

#Simulated data example


install.packages("pls")
install.packages("MASS") # For simulating correlated data
install.packages("ggplot2")

library(pls)
library(MASS) # For simulating correlated data
library(ggplot2)


## -----------------------------
##
## PART 1: Simulation
## -----------------------------

# Start measuring time
start_time <- Sys.time()

# Set seed for reproducibility
set.seed(123)

# Number of participants
n <- 4644

# Simulate highly collinear predictor variables (proteins)
Sigma_x <- matrix(0.8, 2000, 2000) + diag(0.2, 2000) # Covariance matrix for predictors
X <- mvrnorm(n, rep(0, 2000), Sigma_x) # Simulated predictors

# Simulate highly collinear response variables (white matter fractional anisotropy)
Sigma_y <- matrix(0.8, 27, 27) + diag(0.2, 27) # Covariance matrix for responses
Y <- mvrnorm(n, rep(0, 27), Sigma_y) # Simulated responses

# Perform PLS2
pls2_model <- plsr(Y ~ X, ncomp = 10, method = "kernelpls", scale = TRUE, validation = "CV")

# Summary of the model
summary(pls2_model)

# Plotting
# Scores plot for the first two components
scores <- scores(pls2_model)[,1:2] # Extract scores for the first two components

# Convert to a data frame for ggplot
scores_df <- data.frame(Component1 = scores[,1], Component2 = scores[,2])

ggplot(scores_df, aes(x = Component1, y = Component2)) +
  geom_point(alpha = 0.5) +
  ggtitle("Scores Plot for the First Two PLS Components") +
  xlab("Component 1") +
  ylab("Component 2")

# Stop measuring time
end_time <- Sys.time()

# Calculate the elapsed time
elapsed_time <- end_time - start_time

# Save the execution time to a log file
log_file <- "execution_log.txt"
cat("Linear regression execution time:", elapsed_time, "seconds", file = log_file)

# Print a message to the console
cat("Linear regression completed in", elapsed_time, "seconds. Log saved to", log_file, "\n")


## -----------------------------
##
## END Of PART 1: Simulation
## -----------------------------





#PLS2 and ptPLS2

# 
#  PLS2 Function
#  Input: Predictor matrix X, response matrix Y, and the number of latent variables A.
#  Output: Score matrix T, post-transformed score matrix Tpt, weight matrix W, W* (a transformed weight matrix), and regression coefficient matrix B.
#  Process:
#  Initializes matrices W, P (loading matrix), T (score matrix), and E (residual matrix) with zeros.

#  In a loop for each latent variable:
#  Calculates the weight vectors by decomposing the covariance matrix of Y and the current residuals E using singular value decomposition (SVD).
#  Projects E onto the weights to get scores.
#  Computes loadings P and deflates E by removing the explained variance.
#  Calculates W* and B for regression coefficients.
#  Determines if post-transformation is needed based on the number of significant singular values.
#  If post-transformation is needed, it calls step.3 to compute the post-transformed score matrix Tpt.
#  Returns a list containing T, Tpt, W, W*, and B.
#  step.2 and step.3 Functions
#  Used within the PLS2 function to compute matrices necessary for the post-transformation of the scores, depending on the significance of the singular values from the decomposition of YXW.
#  oCPLS2 Function
#  Input: Similar to PLS2 but includes an additional matrix Z for constraints.
#  Output: Similar to PLS2.
#  Process:
#  Initial steps similar to PLS2, but incorporates the constraint matrix Z to modify the predictor matrix before the main algorithm begins.
#  The rest of the process is similar to PLS2, including the optional post-transformation step.
#  KPLS2 Function
#  Input: Predictor matrix X, response matrix Y, number of latent variables A, kernel function parameters (a, b, c), and a test set Xtest.
#  Output: Normalized score matrix Tn for the training set, normalized score matrix Tnpred for the test set, calculated and predicted responses.
#  Process:
#  Computes a kernel matrix K using the polynomial kernel function with the given parameters.
#  In a loop for each latent variable:
#  Extracts scores and updates F (similar to Y but deflated in each iteration).
#  Updates the kernel matrix K and deflates F.
#  Calculates the predicted response for the test set using the kernel matrix of the test set.
#  Returns normalized score matrices for both training and test sets, along with calculated and predicted responses.
#  K.matrix and K.matrixpred Functions
#  Used within the KPLS2 function to compute the kernel matrix for the training set and the prediction kernel matrix for the test set, respectively, using the specified polynomial kernel function.
#  This script is structured to offer flexibility in analyzing associations between various types of data through PLS2, with modifications to incorporate constraints (oCPLS2) and non-linear relationships through kernel methods (KPLS2).
#  

# Define the PLS2 function
PLS2 <- function(X, Y, A) {
  # Initialize matrices for weights (W), loadings (P), scores (T), and residuals (E)
  W <- matrix(rep(0, A * ncol(X)), nrow = ncol(X), ncol = A)
  P <- matrix(rep(0, A * ncol(X)), nrow = ncol(X), ncol = A)
  T <- matrix(rep(0, A * nrow(X)), nrow = nrow(X), ncol = A)
  E <- X # Initial residuals are the predictors themselves
  
  # Iteratively calculate weights, scores, and loadings for A latent variables
  for (i in 1:A) {
    # Calculate weight vector for the current latent variable
    W[, i] <- svd(t(Y) %*% E)$v[, 1]
    # Calculate scores
    T[, i] <- E %*% W[, i]
    # Calculate loadings
    P[, i] <- t(E) %*% T[, i] / sum(T[, i] * T[, i])
    # Update residuals
    E <- E - T[, i] %*% t(P[, i])
  }
  
  # Calculate the transformed weight matrix W*
  Ws <- W %*% solve(t(P) %*% W)
  # Calculate the regression coefficient matrix B
  B <- Ws %*% solve(t(T) %*% T) %*% t(T) %*% Y
  
  # Determine necessity for post-transformation based on singular values
  d1 <- svd(t(Y) %*% X %*% W)
  Np <- length(which(d1$d > 10 ^ -8))
  G <- step.2(X, Y, W) # Helper function for post-transformation
  if (A <= Np)
    Tpt <- T
  if (A > Np)
    Tpt <- step.3(X, Y, W, G)$Tpt # Perform post-transformation if needed
  
  # Compile results into a list
  pls2 <- list(T = T, Tpt = Tpt, W = W, Wstar = Ws, B = B)
  return(pls2)
}

# Helper function for step 2 in post-transformation
step.2 <- function(X, Y, W) {
  # Initialize matrix G for transformation
  G <- matrix(rep(0, ncol(W) * ncol(W)), ncol = ncol(W))
  I <- diag(rep(1, ncol(W))) # Identity matrix
  
  # Singular value decomposition to determine significant singular values
  d1 <- svd(t(Y) %*% X %*% W)
  N <- length(which(d1$d > 10 ^ -8))
  
  # Compute matrix for post-transformation
  V <- Re(d1$v[, 1:N])
  d2 <- eigen((I - V %*% t(V)) %*% (t(W) %*% t(X) %*% X %*% W))
  M <- length(which(Re(d2$values) > 10 ^ -8))
  Go <- Re(d2$vectors[, 1:M])
  Gp <- Re(eigen(t(I - Go %*% t(Go)) %*% (I - Go %*% t(Go)))$vectors[, 1:(ncol(W) - M)])
  G <- cbind(Go, Gp)
  return(G)
}

# Helper function for step 3 in post-transformation
step.3 <- function(X, Y, W, G) {
  # Initialize transformed score matrix (Tpt) and residuals (E, F)
  Tpt <- matrix(rep(0, nrow(X) * ncol(W)), ncol = ncol(W))
  E <- X
  F <- Y
  Wpt <- W %*% G # Transformed weight matrix
  
  # Calculate transformed scores and update residuals
  for (i in 1:ncol(W)) {
    Tpt[, i] <- E %*% Wpt[, i]
    Q <- Tpt[, i] %*% t(Tpt[, i]) / sum(Tpt[, i] ^ 2)
    E <- E - Q %*% E
    F <- F - Q %*% F
  }
  
  # Calculate additional matrices for the post-transformed model
  Ppt <- t(X) %*% Tpt %*








#Input: matrix X of the predictors, matrix Y of the responses, number A of latent variables of the
#model.
#Output: score matrix T, score matrix of the post-transformed model Tpt, weight matrix W, matrix

Wstar (W * ), matrix of the regression coefficient B.
PLS2 <- function(X, Y, A) {
  W <- matrix(rep(0, A * ncol(X)), nrow = ncol(X), ncol = A)
  P <- matrix(rep(0, A * ncol(X)), nrow = ncol(X), ncol = A)
  T <- matrix(rep(0, A * nrow(X)), nrow = nrow(X), ncol = A)
  E <- X
  for (i in 1:A) {
    W[, i] <- svd(t(Y) %*% E)$v[, 1]
    T[, i] <- E %*% W[, i]
    P[, i] <- t(E) %*% T[, i] / sum(T[, i] * T[, i])
    E <- E - T[, i] %*% t(P[, i])
  }
  Ws <- W %*% solve(t(P) %*% W)
  B <- Ws %*% solve(t(T) %*% T) %*% t(T) %*% Y
  d1 <- svd(t(Y) %*% X %*% W)
  Np <- length(which(d1$d > 10 ^ -8))
  G <- step.2(X, Y, W)
  if (A <= Np)
    Tpt <- T
  if (A > Np)
    Tpt <- step.3(X, Y, W, G)$Tpt
  pls2 <- list(
    T = T,
    Tpt = Tpt,
    W = W,
    Wstar = Ws ,
    B = B
  )
  return(pls2)
}
step.2 <- function(X, Y, W) {
  G <- matrix(rep(0, ncol(W) * ncol(W)), ncol = ncol(W))
  I <- diag(rep(1, ncol(W)))
  d1 <- svd(t(Y) %*% X %*% W)
  N <- length(which(d1$d > 10 ^ -8))
  V <- Re(d1$v[, 1:N])
  Metabolites 2019, 9, 51
  doi:10.3390 / metabo9030051 S2 of S4
  d2 <- eigen((I - V %*% t(V)) %*% (t(W) %*% t(X) %*% X %*% W))
  M <- length(which(Re(d2$values) > 10 ^ -8))
  Go <- Re(d2$vectors[, 1:M])
  Gp <-
    Re(eigen(t(I - Go %*% t(Go)) %*% (I - Go %*% t(Go)))$vectors[, 1:(ncol(W) -
                                                                        M)])
  G <- cbind(Go, Gp)
  return(G)
}
step.3 <- function(X, Y, W, G) {
  Tpt <- matrix(rep(0, nrow(X) * ncol(W)), ncol = ncol(W))
  E <- X
  F <- Y
  Wpt <- W %*% G
  for (i in 1:ncol(W)) {
    Tpt[, i] <- E %*% Wpt[, i]
    Q <- Tpt[, i] %*% t(Tpt[, i]) / sum(Tpt[, i] ^ 2)
    E <- E - Q %*% E
    F <- F - Q %*% F
  }
  Ppt <- t(X) %*% Tpt %*% solve(t(Tpt) %*% Tpt)
  Wspt <- Wpt %*% solve(t(Ppt) %*% Wpt)
  ptmodel <- list(Tpt = Tpt, Wpt = Wpt, Wspt = Wspt)
  return(ptmodel)
}


# oCPLS2
#Input:matrix X of the predictors, matrix Y of the responses, matrix Z of the constraints, number A of
#latent variables of the model.
#Output:score matrix T, score matrix of the post - transformed model Tpt, weight matrix W, matrix
#Wstar (W * ), matrix of the regression coefficient B.

oCPLS2 <- function(X, Y, Z, A) {
  W <- matrix(rep(0, A * ncol(X)), nrow = ncol(X), ncol = A)
  P <- matrix(rep(0, A * ncol(X)), nrow = ncol(X), ncol = A)
  T <- matrix(rep(0, A * nrow(X)), nrow = nrow(X), ncol = A)
  B <- t(Z) %*% X
  h <- svd(B)
  R <- length(which(h$d > 10 ^ -8))
  V <- h$v[, 1:R]
  Q <- diag(rep(1, ncol(X))) - V %*% t(V)
  E <- X
  F <- Y
  for (i in 1:A) {
    W[, i] <- svd(t(F) %*% E %*% Q)$v[, 1]
    T[, i] <- E %*% W[, i]
    P[, i] <- t(E) %*% T[, i] / sum(T[, i] * T[, i])
    E <- E - T[, i] %*% t(P[, i])
    F <- F - T[, i] %*% t(T[, i]) %*% F / sum(T[, i] * T[, i])
  }
  Ws <- W %*% solve(t(P) %*% W)
  B <- Ws %*% solve(t(T) %*% T) %*% t(T) %*% Y
  d1 <- svd(t(Y) %*% X %*% W)
  Np <- length(which(d1$d > 10 ^ -8))
  G <- step.2(X, Y, W)
  if (A <= Np)
    Tpt <- T
  if (A > Np)
    Tpt <- step.3(X, Y, W, G)$Tpt
  ocpls2 <- list(
    T = T,
    Tpt = Tpt,
    W = W,
    Wstar = Ws,
    B = B
  )
  return(ocpls2)
}

# KPLS2
#Input: matrix X of the predictors, matrix Y of the responses, number A of latent variables of the
#model, (a, b, c) parameters of the polynomial kernel k(x,y)=(a(xty)+b)p, matrix Xtest of the test set to
#be predicted
#Output: normalized score matrix Tn, normalized score matrix of the test set Tnpred, calculated
#response Ycalc, predicted response Ypred

KPLS2 <- function(X, Y, A, a, b, p, Xtest) {
  Tn <- matrix(rep(0, nrow(X) * A), ncol = A)
  U <- matrix(rep(0, nrow(Y) * A), ncol = A)
  K <- K.matrix(X, a, b, p)
  F <- Y
  for (i in 1:A) {
    Tn[, i] <- Re(eigen(K %*% F %*% t(F))$vectors[, 1])
    U[, i] <- F %*% t(F) %*% Tn[, i]
    Q <- diag(1, nrow(X)) - Tn[, i] %*% t(Tn[, i])
    K <- Q %*% K %*% Q
    F <- Q %*% F
  }
  Ymod <- Y - F
  K <- K.matrix(X, a, b, p)
  Kp <- K.matrixpred(Xtest, X, a, b, p)
  Ypred <- Kp %*% U %*% solve(t(Tn) %*% K %*% U) %*% t(Tn) %*% Y
  Tnpred <- Kp %*% U %*% solve(t(Tn) %*% K %*% U)
  r <- list(
    Tn = Tn ,
    Tnpred = Tnpred,
    Ycalc = Ymod,
    Ypred = Ypred
  )
  return(r)
  
}
K.matrix <- function(X, a, b, p) {
  Km <- matrix(rep(-999, nrow(X) * nrow(X)), nrow = nrow(X))
  for (i in 1:nrow(X)) {
    for (j in 1:nrow(X)) {
      Km[i, j] <- (a * t(X[i, ]) %*% X[j, ] + b) ^ p
    }
  }
  p <- matrix(rep(1, nrow(X)), ncol = 1)
  Qc <- diag(1, nrow(X)) - (1 / nrow(X)) * p %*% t(p)
  Km <- Qc %*% Km %*% Qc
  return(Km)
}
K.matrixpred <- function(X1, X2, a, b, p) {
  Kp <- matrix(rep(-999, nrow(X1) * nrow(X2)), nrow = nrow(X1))
  K <- matrix(rep(-999, nrow(X2) * nrow(X2)), nrow = nrow(X2))
  for (i in 1:nrow(X1)) {
    for (j in 1:nrow(X2)) {
      Kp[i, j] <- (a * t(X1[i, ]) %*% X2[j, ] + b) ^ p
    }
  }
  for (i in 1:nrow(X2)) {
    for (j in 1:nrow(X2)) {
      K[i, j] <- (a * t(X2[i, ]) %*% X2[j, ] + b) ^ p
    }
  }
  p <- matrix(rep(1, nrow(X2)), ncol = 1)
  pp <- matrix(rep(1, nrow(X1)), ncol = 1)
  Qc <- diag(1, nrow(X2)) - (1 / nrow(X2)) * p %*% t(p)
  P <- (1 / nrow(X2)) * pp %*% t(p)
  Kp <- (Kp - P %*% K) %*% Qc
  return(Kp)
}