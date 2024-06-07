library(ggplot2)
library(MASS)
library(energy)
library(dHSIC)

dCor_own <- function(x, y, R) {
  
  x <- dist(x)
  y <- dist(y)
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  a_bar_k <- rowMeans(x)
  a_bar <- mean(x)
  A_kj <- sweep(x, 1, a_bar_k)
  A_jk <- sweep(A_kj, 2, a_bar_k)
  A <- A_jk + a_bar
 
  b_bar_k <- rowMeans(y)
  b_bar <- mean(y)
  B_kj <- sweep(y, 1, b_bar_k)
  B_jk <- sweep(B_kj, 2, b_bar_k)
  B <- B_jk + b_bar
  
  dCov <- sqrt(mean(A * B))
  dVarX <- sqrt(mean(A * A))
  dVarY <- sqrt(mean(B * B))
  
  dCor <- dCov / sqrt(dVarX * dVarY)
  
  reps <- numeric(R)
  n <- nrow(x)

  for (r in 1:R) {
    perm <- sample(1:n, n, replace = FALSE)
    permutedB <- B[perm, perm]
    dcov_perm <- mean(A * permutedB)
    reps[r] <- sqrt(dcov_perm) / (sqrt(dVarX) * sqrt(dVarY))
  }
  
  p.value <- (sum(reps >= dCor) + 1) / (R + 1)
  
  return(list(dCor = dCor, p_value = p.value))
}


rbf_matrix <- function(X, Y, sigma) {
  
  if (ncol(X)==1) {
  G <- X^2
} else {
  G <- rowSums((X^2))
}
  if (ncol(Y)==1) {
  H <- Y^2
} else {
  H <- rowSums((Y^2))
}
  
  Q <- matrix(G, nrow = nrow(X), ncol = nrow(Y), byrow = TRUE)
  R <- matrix(H, nrow = nrow(X), ncol = nrow(Y), byrow = FALSE)
  H <- Q + R - 2 * tcrossprod(X, Y)
  H <- exp(-H / (2 * sigma^2))
  return(H)
}


hsicTestBoot <- function(X, Y, bandwidth_x, bandwidth_y, R) {
  
  m <- nrow(X)
  
  K <- rbf_matrix(X, X, bandwidth_x)
  L <- rbf_matrix(Y, Y, bandwidth_y)
  
  bone <- matrix(1, nrow = m, ncol = 1)
  H <- diag(m) - 1/m * matrix(1, nrow = m, ncol = m)
  
  Kc <- H %*% K %*% H
  
  testStat <- sum(Kc * L) / (m^2)
  
  HSICarr <- numeric(R)
  
  for (whichSh in 1:R) {
    indL <- sample(m)
    HSICarr[whichSh] <- sum(Kc * L[indL, indL]) / (m^2)
  }
  
  p_value <- sum(HSICarr >= testStat) / R
  
  return(list(testStat = testStat, p_value = p_value))
}


hsicTestGamma <- function(X, Y, alpha, bandwidth_x, bandwidth_y,) {
  m <- nrow(X)
  bone <- matrix(1, nrow = m, ncol = 1)
  H <- diag(m) - 1/m * matrix(1, nrow = m, ncol = m)

  K <- rbf_matrix(X, X, bandwidth_x)
  L <- rbf_matrix(Y, Y, bandwidth_y)

  Kc <- H %*% K %*% H
  Lc <- H %*% L %*% H

  testStat <- 1/m * sum(Kc * Lc)

  varHSIC <- (1/6 * Kc * Lc)^2
  varHSIC <- (1/m/(m-1)) * (sum(varHSIC) - sum(diag(varHSIC)))
  varHSIC <- 72 * (m-4) * (m-5) / m / (m-1) / (m-2) / (m-3) * varHSIC

  K <- K - diag(diag(K))
  L <- L - diag(diag(L))

  muX <- 1/m/(m-1) * t(bone) %*% (K %*% bone)
  muY <- 1/m/(m-1) * t(bone) %*% (L %*% bone)

  mHSIC <- 1/m * (1 + muX * muY - muX - muY)
  al <- mHSIC^2 / varHSIC
  bet <- varHSIC * m / mHSIC

  thresh <- qgamma(1 - alpha, shape = al, scale = bet)

  return(list(thresh = thresh, testStat = testStat, al = al, bet = bet))
}
