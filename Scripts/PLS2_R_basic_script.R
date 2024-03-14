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
setwd("/Users/eleanorc_worklaptop/desktop/UKB_Proteomics")


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