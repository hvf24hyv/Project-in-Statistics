library(ggplot2)
library(MASS)
library(energy)
library(dHSIC)

set.seed(123)
ds <- 2
rhos <- c(0.01, 0.05, 0.1, 0.15, 0.3)
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500

create_cov_matrix <- function(rho, p, d) {
  Sigma <- diag(1, d)
  Sigma[1:(d/2), ((d/2)+1):d] <- rho
  Sigma[((d/2)+1):d, 1:(d/2)] <- rho
  return(Sigma)
}

adjust_for_positivity <- function(Sigma) {
  epsilon <- 0.01
  while(!is_positive_definite(Sigma)) {
    Sigma <- Sigma + diag(rep(epsilon, nrow(Sigma)))
    epsilon <- epsilon + 0.01
  }
  return(Sigma)
}

is_positive_definite <- function(matrix) {
  eigenvalues <- eigen(matrix)$values
  return(all(eigenvalues > 0))
}

all_results <- list()

for (alpha in alphas) {
  d <- ds
  p <- d / 2
  q <- d / 2
  results <- vector("list", length(rhos))
  names(results) <- as.character(rhos)

  for (i in seq_along(rhos)) {
    rho <- rhos[i]
    Sigma <- create_cov_matrix(rho, p, d)

    if (!is_positive_definite(Sigma)) {
      Sigma <- adjust_for_positivity(Sigma)
    }

    test_results <- replicate(n_simulations, {
      data <- mvrnorm(n = 200, mu = rep(0, d), Sigma = Sigma)
      dcor.test(data[, 1:p], data[, (p+1):d], R=200)$p.value < alpha
    })
    rejection_rate <- mean(test_results)
    se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
    results[[i]] <- c(rejection_rate, se)
  }

  all_results[[as.character(alpha)]] <- data.frame(
    rho = rhos,
    rejection_rate = sapply(results, function(x) x[1]),
    se = sapply(results, function(x) x[2]),
    alpha = rep(alpha, length(rhos))
  )
}

all_results_df <- do.call(rbind, all_results)

plot <- ggplot(all_results_df, aes(x = alpha, y = rejection_rate, color = factor(rho))) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se, ymax = rejection_rate + 1.96 * se, fill = factor(rho)), alpha = 0.2) +
  labs(title = "Gaussian Samples",
       x = "Significance Level (Alpha)",
       y = "Rejection Rate",
       color = "Rho") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )+
  guides(fill = FALSE)

print(plot)



set.seed(123)
d <- 200
fixed_rho <- 0.01
num_pairs <- seq(20, 100, by =20) 
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500
library(MASS) 
library(energy) 
create_cov_matrix <- function(rho, num_pairs, d) {
Sigma <- diag(1, d)

if (num_pairs > 0) {
Sigma[1:num_pairs, (d/2 +1):(d/2+num_pairs)] <- rho
Sigma[(d/2+1):(d/2+num_pairs), 1:num_pairs] <- rho
}
return(Sigma)
}
adjust_for_positivity <- function(Sigma) {
epsilon <- 0.01
while(!is_positive_definite(Sigma)) {
Sigma <- Sigma + diag(rep(epsilon, nrow(Sigma)))
epsilon <- epsilon + 0.01
}
return(Sigma)
}
is_positive_definite <- function(matrix) {
eigenvalues <- eigen(matrix)$values
return(all(eigenvalues > 0))
}
all_results <- list()
for (alpha in alphas) {
results <- vector("list", length(num_pairs))
names(results) <- as.character(num_pairs)
for (i in seq_along(num_pairs)) {
n_pairs <- num_pairs[i]
Sigma <- create_cov_matrix(fixed_rho, n_pairs, d)
if (!is_positive_definite(Sigma)) {
Sigma <- adjust_for_positivity(Sigma)
}
test_results <- replicate(n_simulations, {
data <- mvrnorm(n = 200, mu = rep(0, d), Sigma = Sigma)
dcor.test(data[, 1:n_pairs], data[, (n_pairs+1):(2*n_pairs)], R=200)$p.value < alpha
})
rejection_rate <- mean(test_results)
se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
results[[i]] <- c(rejection_rate, se)
}
all_results[[as.character(alpha)]] <- data.frame(
num_pairs = num_pairs,
rejection_rate = sapply(results, function(x) x[1]),
se = sapply(results, function(x) x[2]),
alpha = rep(alpha, length(num_pairs))
)
}
all_results_df <- do.call(rbind, all_results)
library(ggplot2)
plot <- ggplot(all_results_df, aes(x = alpha, y = rejection_rate, color = factor(num_pairs))) +
geom_line() +
geom_ribbon(aes(ymin = rejection_rate - 1.96 * se, ymax = rejection_rate + 1.96 * se, fill = factor(num_pairs)), alpha = 0.2) +
labs(title = "Gaussian Samples",
x = "Significance Level (Alpha)",
y = "Rejection Rate",
color = "Number of Correlated Pairs") +
theme_minimal()+
theme(
plot.title = element_text(size = 20, face = "bold"),
axis.title.x = element_text(size = 16, face = "bold"),
axis.title.y = element_text(size = 16, face = "bold"),
axis.text.x = element_text(size = 12),
axis.text.y = element_text(size = 12),
legend.title = element_text(size = 14),
legend.text = element_text(size = 12),
    legend.position = "bottom"
)+
  guides(fill = FALSE)
print(plot)




 set.seed(123)
ds <- 2
A <- c(1,2,10)
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500

simulate_data <- function(A, N) {
  Theta <- runif(N, 0, 2 * pi)
  eps1 <- rnorm(N, 0, 1) / 4
  eps2 <- rnorm(N, 0, 1) / 4
  X <- A * cos(Theta) + eps1
  Y <- A * sin(Theta) + eps2
  X_2 <- runif(N, 0, 1)
  Y_2 <- runif(N, 0, 1)
  return(data.frame(X, Y, X_2, Y_2))
}

all_results <- list()

for (alpha in alphas) {
  d <- ds
  p <- d / 2
  q <- d / 2
  results <- vector("list", length(A))
  names(results) <- as.character(A)

  for (i in seq_along(A)) {
    
    test_results <- replicate(n_simulations, {
      data <- simulate_data(A[i],200)
      dcor.test(data[,c('X','X_2')], data[,c('Y','Y_2')], R=100)$p.value < alpha
    })
    rejection_rate <- mean(test_results)
    se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
    results[[i]] <- c(rejection_rate, se)
  }

  all_results[[as.character(alpha)]] <- data.frame(
    A = A,
    rejection_rate = sapply(results, function(x) x[1]),
    se = sapply(results, function(x) x[2]),
    alpha = rep(alpha, length(A))
  )
}

all_results_df <- do.call(rbind, all_results)

plot <- ggplot(all_results_df, aes(x = alpha, y = rejection_rate, color = factor(A))) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se, ymax = rejection_rate + 1.96 * se, fill = factor(A)), alpha = 0.2) +
  labs(title = "M1",
       x = "Significance Level (Alpha)",
       y = "Rejection Rate",
       color = "A") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )+
  guides(fill = FALSE)

print(plot)




set.seed(123)
dimensions <- c(2,10,50,100,200)
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500
A <- 3

simulate_data <- function(A, N, d) {
  Theta <- runif(N, 0, 2 * pi)
  eps1 <- rnorm(N, 0, 1) / 4
  eps2 <- rnorm(N, 0, 1) / 4
  X <- A * cos(Theta) + eps1
  Y <- A * sin(Theta) + eps2
  
  df <- data.frame(X1 = X, Y1 = Y)
  
    for (i in 2:d) {
      df[[paste0("X", i)]] <- runif(N, 0, 1)
      df[[paste0("Y", i)]] <- runif(N, 0, 1)
    }
  
  return(df)
}

all_results <- list()

for (alpha in alphas) {
  for (d in dimensions) {
    p <- d / 2
    q <- d / 2
    
    test_results <- replicate(n_simulations, {
      data <- simulate_data(A, 30, d)
    x_cols <- paste("X", 1:d, sep='')
    y_cols <- paste("Y", 1:d, sep='')
      dcor.test(data[, c(x_cols)], data[, c(y_cols)], R=100)$p.value < alpha
    })
    rejection_rate <- mean(test_results)
    se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
    
    all_results[[paste(alpha, d)]] <- data.frame(
      Dimension = d,
      Rejection_rate = rejection_rate,
      SE = se,
      Alpha = alpha
    )
  }
}

all_results_df <- do.call(rbind, all_results)

library(ggplot2)
plot <- ggplot(all_results_df, aes(x = Alpha, y = Rejection_rate, color = as.factor(Dimension))) +
  geom_line() +
  geom_ribbon(aes(ymin = Rejection_rate - 1.96 * SE, ymax = Rejection_rate + 1.96 * SE, fill = as.factor(Dimension)), alpha = 0.2) +
  labs(title = "M1",
       x = "Significance Level (Alpha)",
       y = "Rejection Rate",
       color = "Dimension") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  guides(fill = FALSE)

print(plot)





set.seed(123)
ds <- 2
rhos <- c(0.1,0.5,1)
#rhos <- seq(0.1, 1, by=0.2)
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500

simulate_data <- function(n, rho) {
  X <- runif(n, -1, 1)  
  X_2 <- runif(n, 0, 1)
  epsilon <- rnorm(n)   
  Y <- abs(X)^rho * epsilon
  Y_2 <- runif(n, 0, 1) 
  return(data.frame(X, Y, X_2, Y_2))
}

all_results <- list()

for (alpha in alphas) {
  results <- vector("list", length(rhos))
  names(results) <- as.character(rhos)

  for (i in seq_along(rhos)) {
    
    test_results <- replicate(n_simulations, {
      data <- simulate_data(200,rhos[i])
      dcor.test(data[,c('X','X_2')], data[,c('Y','Y_2')], R=150)$p.value < alpha
    })
    rejection_rate <- mean(test_results)
    se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
    results[[i]] <- c(rejection_rate, se)
  }

  all_results[[as.character(alpha)]] <- data.frame(
    rhos = rhos,
    rejection_rate = sapply(results, function(x) x[1]),
    se = sapply(results, function(x) x[2]),
    alpha = rep(alpha, length(rhos))
  )
}

all_results_df <- do.call(rbind, all_results)

plot <- ggplot(all_results_df, aes(x = alpha, y = rejection_rate, color = factor(rhos))) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se, ymax = rejection_rate + 1.96 * se, fill = factor(rhos)), alpha = 0.2) +
  labs(title = "M2",
       x = "Significance Level (Alpha)",
       y = "Rejection Rate",
       color = "Rho") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )+
  guides(fill = FALSE)

print(plot)






set.seed(123)
ds <- c(2,10,50,100,200)
rho <- 1
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500

simulate_data <- function(n, rho, d) {
  X <- runif(n, -1, 1)  
  epsilon <- rnorm(n)   
  Y <- abs(X)^rho * epsilon
  df <- data.frame(X1 = X, Y1 = Y)
  
    for (i in 2:d) {
      df[[paste0("X", i)]] <- runif(n, 0, 1)
      df[[paste0("Y", i)]] <- runif(n, 0, 1)
    }
  
  return(df)
}

all_results <- list()

for (alpha in alphas) {
  results <- vector("list", length(ds))
  names(results) <- as.character(ds)

  for (i in seq_along(ds)) {
    
    test_results <- replicate(n_simulations, {
      data <- simulate_data(200,rho,ds[i])
      x_cols <- paste("X", 1:ds[i], sep='')
      y_cols <- paste("Y", 1:ds[i], sep='')
      dcor.test(data[, c(x_cols)], data[, c(y_cols)], R=100)$p.value < alpha
    })
    rejection_rate <- mean(test_results)
    se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
    results[[i]] <- c(rejection_rate, se)
  }

  all_results[[as.character(alpha)]] <- data.frame(
    ds = ds,
    rejection_rate = sapply(results, function(x) x[1]),
    se = sapply(results, function(x) x[2]),
    alpha = rep(alpha, length(ds))
  )
}

all_results_df <- do.call(rbind, all_results)

plot <- ggplot(all_results_df, aes(x = alpha, y = rejection_rate, color = factor(ds))) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se, ymax = rejection_rate + 1.96 * se, fill = factor(ds)), alpha = 0.2) +
  labs(title = "M2",
       x = "Significance Level (Alpha)",
       y = "Rejection Rate",
       color = "Dimension") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )+
  guides(fill = FALSE)

print(plot)









set.seed(123)
correlated_pairs <- c(2,100)
rho <- 1
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500
n <- 200

simulate_data <- function(n, rho, num_correlated_pairs) {
  X <- runif(200, -1, 1)  
  epsilon <- rnorm(200)   
  Y <- abs(X)^rho * epsilon
  df <- data.frame(X1 = X, Y1 = Y)
  
  for (i in 2:num_correlated_pairs) {
    X <- runif(200, -1, 1)
    epsilon <- rnorm(200)
    Y <- abs(X)^rho * epsilon
    df[[paste0("X", i)]] <- X
    df[[paste0("Y", i)]] <- Y
  }
  
  for (i in (num_correlated_pairs + 1):n) {
    df[[paste0("X", i)]] <- runif(200, 0, 1)
    df[[paste0("Y", i)]] <- runif(200, 0, 1)
  }
  
  return(df)
}

all_results <- list()

for (alpha in alphas) {
  results <- vector("list", length(correlated_pairs))
  names(results) <- as.character(correlated_pairs)

  for (i in seq_along(correlated_pairs)) {
    num_correlated_pairs <- correlated_pairs[i]
    
    test_results <- replicate(n_simulations, {
      data <- simulate_data(n, rho, num_correlated_pairs)
      x_cols <- paste("X", 1:n, sep='')
      y_cols <- paste("Y", 1:n, sep='')
      dcor.test(data[, c(x_cols)], data[, c(y_cols)], R = 100)$p.value < alpha
    })
    rejection_rate <- mean(test_results)
    se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
    results[[i]] <- c(rejection_rate, se)
  }

  all_results[[as.character(alpha)]] <- data.frame(
    correlated_pairs = correlated_pairs,
    rejection_rate = sapply(results, function(x) x[1]),
    se = sapply(results, function(x) x[2]),
    alpha = rep(alpha, length(correlated_pairs))
  )
}

all_results_df <- do.call(rbind, all_results)

plot <- ggplot(all_results_df, aes(x = alpha, y = rejection_rate, color = factor(correlated_pairs))) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se, ymax = rejection_rate + 1.96 * se, fill = factor(correlated_pairs)), alpha = 0.2) +
  labs(title = "M2",
       x = "Significance Level (Alpha)",
       y = "Rejection Rate",
       color = "Number of Correlated Pairs") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  guides(fill = "none")

print(plot) 








set.seed(123)

ds <- 2
A <- c(1, 2, 10)
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500
n_samples <- 200  

fy_given_x <- function(y, x, a) {
  (1 / (2 * pi)) * (1 + sin(a * x) * sin(a * y))
}

simulate_data <- function(a, N) {
  X <- rep(NA, N)
  Y <- rep(NA, N)
  X_2 <- runif(N, 0, 1)
  Y_2 <- runif(N, 0, 1)
  
  for (i in 1:N) {
    X[i] <- runif(1, -pi, pi)
    
    while (TRUE) {
      y_proposed <- runif(1, -pi, pi)
      f_proposal <- 1 / (2 * pi)
      f_target <- fy_given_x(y_proposed, X[i], a)
      acceptance_prob <- f_target / f_proposal
      
      if (runif(1) < acceptance_prob) {
        Y[i] <- y_proposed
        break
      }
    }
  }
  return(data.frame(X, Y, X_2, Y_2))
}

all_results <- list()

for (alpha in alphas) {
  d <- ds
  p <- d / 2
  q <- d / 2
  results <- vector("list", length(A))
  names(results) <- as.character(A)
  
  for (i in seq_along(A)) {
    test_results <- replicate(n_simulations, {
      data <- simulate_data(A[i], n_samples)
      dcor.test(data[,c('X','X_2')], data[,c('Y','Y_2')], R=100)$p.value < alpha
    })
    rejection_rate <- mean(test_results)
    se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
    results[[i]] <- c(rejection_rate, se)
  }
  
  all_results[[as.character(alpha)]] <- data.frame(
    A = A,
    rejection_rate = sapply(results, function(x) x[1]),
    se = sapply(results, function(x) x[2]),
    alpha = rep(alpha, length(A))
  )
}

all_results_df <- do.call(rbind, all_results)

plot <- ggplot(all_results_df, aes(x = alpha, y = rejection_rate, color = factor(A))) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se, ymax = rejection_rate + 1.96 * se, fill = factor(A)), alpha = 0.2) +
  labs(title = "M3",
       x = "Significance Level (Alpha)",
       y = "Rejection Rate",
       color = "A") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  guides(fill = FALSE)

print(plot)







set.seed(123)
ds <- c(2,200)
A <- 1
alphas <- seq(0.01, 1, by = 0.02)
n_simulations <- 500

fy_given_x <- function(y, x, a) {
  (1 / (2 * pi)) * (1 + sin(a * x) * sin(a * y))
}

simulate_data <- function(a, N, d) {
  X <- rep(NA, N)
  Y <- rep(NA, N)
  for (i in 1:N) {
    X[i] <- runif(1, -pi, pi)
    
    while (TRUE) {
      y_proposed <- runif(1, -pi, pi)
      f_proposal <- 1 / (2 * pi)
      f_target <- fy_given_x(y_proposed, X[i], a)
      acceptance_prob <- f_target / f_proposal
      
      if (runif(1) < acceptance_prob) {
        Y[i] <- y_proposed
        break
      }}}
  df <- data.frame(X1 = X, Y1 = Y)
  
  for (i in 2:d) {
    df[[paste0("X", i)]] <- runif(N, 0, 1)
    df[[paste0("Y", i)]] <- runif(N, 0, 1)
    }
  return(df)
}

all_results <- list()

for (alpha in alphas) {
  results <- vector("list", length(ds))
  names(results) <- as.character(ds)

  for (i in seq_along(ds)) {
    
    test_results <- replicate(n_simulations, {
      data <- simulate_data(A,200,ds[i])
      x_cols <- paste("X", 1:ds[i], sep='')
      y_cols <- paste("Y", 1:ds[i], sep='')
      dcor.test(data[, c(x_cols)], data[, c(y_cols)], R=100)$p.value < alpha
    })
    rejection_rate <- mean(test_results)
    se <- sqrt((rejection_rate * (1 - rejection_rate)) / n_simulations)
    results[[i]] <- c(rejection_rate, se)
  }

  all_results[[as.character(alpha)]] <- data.frame(
    ds = ds,
    rejection_rate = sapply(results, function(x) x[1]),
    se = sapply(results, function(x) x[2]),
    alpha = rep(alpha, length(ds))
  )
}

all_results_df <- do.call(rbind, all_results)

plot <- ggplot(all_results_df, aes(x = alpha, y = rejection_rate, color = factor(ds))) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se, ymax = rejection_rate + 1.96 * se, fill = factor(ds)), alpha = 0.2) +
  labs(title = "M3",
       x = "Significance Level (Alpha)",
       y = "Rejection Rate",
       color = "Dimension") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )+
  guides(fill = FALSE)

print(plot)
