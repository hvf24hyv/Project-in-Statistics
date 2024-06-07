library(ggplot2)
library(MASS)
library(energy)
library(dHSIC)

set.seed(1)
n <- 100
x <- rnorm(n, 0, 1)
y <- rnorm(n, 0, 1)
num_resamples_list <- c(10, 50, 100, 200,250, 500,1000)

pvalues <- vector()
times <- vector()
for (num_resamples in num_resamples_list) {
    start_time <- Sys.time()
    test_result <- dCor_own(x, y, R = num_resamples)
    end_time <- Sys.time()
    
    pvalues <- c(pvalues, test_result$p_value)
    times <- c(times, as.numeric(end_time - start_time, units="secs"))
}

data <- data.frame(num_resamples = num_resamples_list, pvalues = pvalues, times = times)

p1 <- ggplot(data, aes(x = num_resamples, y = pvalues)) +
    geom_line() +geom_line(color = "blue")+
    labs(y = "p-value")+
  theme_minimal()+theme(axis.title.x = element_blank())+labs(title = "Calibration")+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

p2 <- ggplot(data, aes(x = num_resamples, y = times)) +
    geom_line(color = "blue") +
    labs(y = "Time (in seconds)", x = "Number of permutations")+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

gridExtra::grid.arrange(p1, p2, nrow = 2)


n <- 200
n_simulations <- 500
alpha_levels <- seq(0.1, 1, by = 0.1)

bandwidths <- list("median", 0.1, 2, 3)

run_simulation <- function(bandwidth) {
  sapply(alpha_levels, function(alpha) {
    rejections <- replicate(n_simulations, {
      x <- matrix(rnorm(n), nrow = n)
      y <- matrix(rnorm(n), nrow = n)
      
      if (bandwidth == "median") {
        bandwidth_x <- median_bandwidth(x)
        bandwidth_y <- median_bandwidth(y)
      } else {
        bandwidth_x <- bandwidth
        bandwidth_y <- bandwidth
      }

      test_result <- hsicTestBoot(x,y,bandwidth_x,bandwidth_y,100)
      test_result$p_value < alpha
    })
    
    mean_rejections <- mean(rejections)
    se_rejections <- sqrt(mean_rejections * (1 - mean_rejections) / n_simulations)
    lower_bound <- mean_rejections - 1.96 * se_rejections
    upper_bound <- mean_rejections + 1.96 * se_rejections
    
    return(c(mean_rejections, lower_bound, upper_bound))
  })
}

data <- lapply(bandwidths, run_simulation)
names(data) <- bandwidths

plot_data <- do.call(rbind, lapply(1:length(bandwidths), function(i) {
  kernel_data <- data.frame(
    alpha_levels = alpha_levels,
    rejection_rates = data[[i]][1, ],
    lower_bound = data[[i]][2, ],
    upper_bound = data[[i]][3, ],
    bandwidth = as.character(bandwidths[i])
  )
  return(kernel_data)
}))

ggplot(plot_data, aes(x = alpha_levels, y = rejection_rates, color = bandwidth, fill = bandwidth)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Significance Level (Alpha)", y = "Rejection Rate", title = "Rejection Rate vs. Alpha") +
  scale_color_manual(values = c("blue", "green", "red", "purple")) +
  scale_fill_manual(values = c("blue", "green", "red", "purple")) +
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


linear <- function(x, y) {
  t(x)%*%y
}
n <- 200
n_simulations <- 500
alpha_levels <- seq(0.1, 1, by = 0.1)
kernels <- c("gaussian", "linear")

run_simulation <- function(kernel) {
  sapply(alpha_levels, function(alpha) {
    rejections <- replicate(n_simulations, {
      x <- matrix(rnorm(n * 1), nrow = n, ncol = 1)
      y <- matrix(rnorm(n * 1), nrow = n, ncol = 1)

      test_result <- hsicTestBoot(x,y,bandwidth_x,bandwidth_y,100, kernel)
      test_result$p.value < alpha
    })
    
    mean_rejections <- mean(rejections)
    se_rejections <- sqrt(mean_rejections * (1 - mean_rejections) / n_simulations)
    lower_bound <- mean_rejections - 1.96 * se_rejections
    upper_bound <- mean_rejections + 1.96 * se_rejections
    
    return(c(mean_rejections, lower_bound, upper_bound))
  })
}

data <- lapply(kernels, run_simulation)
names(data) <- kernels

plot_data <- do.call(rbind, lapply(1:length(kernels), function(i) {
  kernel_data <- data.frame(
    alpha_levels = alpha_levels,
    rejection_rates = data[[i]][1, ],
    lower_bound = data[[i]][2, ],
    upper_bound = data[[i]][3, ],
    kernel = kernels[i]
  )
  return(kernel_data)
}))

ggplot(plot_data, aes(x = alpha_levels, y = rejection_rates, color = kernel, fill = kernel)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Significance Level (Alpha)", y = "Rejection Rate", title = "Rejection Rate vs. Alpha") +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = FALSE)



set.seed(1)
n <- 200
x <- matrix(rnorm(100 * n, 0, 1), ncol = 100)
y <- matrix(rnorm(100 * n, 0, 1), ncol = 100)
num_resamples_list <- c(10, 50, 100, 200,250, 500,1000)

pvalues <- vector()
times <- vector()
num_eval <- 10
for (num_resamples in num_resamples_list) {
    start_time <- Sys.time()
    test_result <- hsicTestBoot(x,y,bandwidth_x,bandwidth_y,100)
    end_time <- Sys.time()
    
    pvalues <- c(pvalues, test_result$p.value)
    times <- c(times, as.numeric(end_time - start_time, units="secs"))
}

data <- data.frame(num_resamples = num_resamples_list, pvalues = pvalues, times = times)

p1 <- ggplot(data, aes(x = num_resamples, y = pvalues)) +
    geom_line() +geom_line(color = "blue")+
    labs(y = "p-value")+
  theme_minimal()+theme(axis.title.x = element_blank())+labs(title = "Calibration")+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

p2 <- ggplot(data, aes(x = num_resamples, y = times)) +
    geom_line(color = "blue") +
    labs(y = "Time (in seconds)", x = "Number of permutations")+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

gridExtra::grid.arrange(p1, p2, nrow = 2)
