# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)
library(bayesplot)

# Set seed for reproducibility
set.seed(123)

# Function to simulate data from Hoppe urn model
simulate_hoppe_urn <- function(alpha, base_distribution, total_draws) {
  urn <- rep("mutator", alpha)
  color_count <- data.frame(color = base_distribution, count = 0, stringsAsFactors = FALSE)

  for (n in 1:total_draws) {
    if (runif(1) < alpha / (alpha + n)) {
      new_color <- sample(base_distribution, 1)
      urn <- c(urn, new_color)
      color_count$count[color_count$color == new_color] <- color_count$count[color_count$color == new_color] + 1
    } else {
      drawn_ball <- sample(urn, 1)
      urn <- c(urn, drawn_ball)
      if (drawn_ball != "mutator") {
        color_count$count[color_count$color == drawn_ball] <- color_count$count[color_count$color == drawn_ball] + 1
      }
    }
  }

  return(color_count)
}

# Define base distribution and alpha values
base_distribution <- c("red", "blue", "green", "yellow", "orange")
alpha_values <- c(1, 2, 3, 5)
total_draws <- 500

# Set a non-uniform prior for the Dirichlet distribution parameters
prior_alpha <- c(2, 3, 1, 4, 5)

# Simulate data for each alpha
color_counts <- list()
for (alpha in alpha_values) {
  color_counts[[paste0("alpha_", alpha)]] <- simulate_hoppe_urn(alpha, base_distribution, total_draws)
}

# Stan model code as a string
stan_model_code <- "
data {
  int<lower=0> N;               // Number of observations
  int<lower=0> K;               // Number of colors
  int<lower=0> counts[K];       // Counts of each color
  vector<lower=0>[K] alpha;     // Concentration parameter
}
parameters {
  simplex[K] theta;             // Proportions for each color
}
model {
  theta ~ dirichlet(alpha);     // Prior
  counts ~ multinomial(theta);  // Likelihood
}
"

# Compile the Stan model from the string
stan_model <- stan_model(model_code = stan_model_code)

# Function to fit the model using MCMC
fit_model <- function(counts, alpha) {
  data_list <- list(N = sum(counts), K = length(counts), counts = counts, alpha = alpha)
  fit <- sampling(stan_model, data = data_list, iter = 2000, chains = 4)
  return(fit)
}

# Fit the model for each alpha
fits <- list()
for (alpha in alpha_values) {
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  fits[[paste0("alpha_", alpha)]] <- fit_model(counts, prior_alpha)
}

# Function to plot prior, observed, and posterior density
plot_results <- function(fit, counts, alpha, base_distribution) {
  posterior_samples <- extract(fit)$theta
  posterior_density <- apply(posterior_samples, 2, mean)

  prior_dist <- data.frame(color = base_distribution, count = prior_alpha)
  observed_dist <- data.frame(color = base_distribution, count = counts)
  posterior_dist <- data.frame(color = base_distribution, density = posterior_density)

  # Plot prior distribution
  prior_plot <- ggplot(data = prior_dist, aes(x = color, y = count, fill = color)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Prior Distribution (alpha =", alpha, ")")) +
    xlab("Color") +
    ylab("Count") +
    theme_minimal() +
    scale_fill_manual(values = c("red", "blue", "green", "yellow", "orange"))
  ggsave(paste0("prior_distribution_alpha_", alpha, ".png"), plot = prior_plot)

  # Plot observed distribution
  observed_plot <- ggplot(data = observed_dist, aes(x = color, y = count, fill = color)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Observed Distribution (alpha =", alpha, ")")) +
    xlab("Color") +
    ylab("Count") +
    theme_minimal() +
    scale_fill_manual(values = c("red", "blue", "green", "yellow", "orange"))
  ggsave(paste0("observed_distribution_alpha_", alpha, ".png"), plot = observed_plot)

  # Plot posterior density overlaying observed distribution
  combined_plot <- ggplot() +
    geom_bar(data = observed_dist, aes(x = color, y = count, fill = color), stat = "identity", alpha = 0.6) +
    geom_line(data = posterior_dist, aes(x = color, y = density * sum(counts), group = 1, color = color), size = 1.5) +
    ggtitle(paste("Observed and Posterior Density (alpha =", alpha, ")")) +
    xlab("Color") +
    ylab("Count / Density") +
    theme_minimal() +
    scale_fill_manual(values = c("red", "blue", "green", "yellow", "orange")) +
    scale_color_manual(values = c("red", "blue", "green", "yellow", "orange"))
  ggsave(paste0("combined_posterior_density_alpha_", alpha, ".png"), plot = combined_plot)

  # Plot MCMC convergence diagnostics using bayesplot
  mcmc_trace_plot <- mcmc_trace(as.array(fit), pars = "theta[1]") + ggtitle(paste("MCMC Trace Plot (alpha =", alpha, ")"))
  ggsave(paste0("mcmc_trace_plot_alpha_", alpha, ".png"), plot = mcmc_trace_plot)

  mcmc_density_plot <- mcmc_dens_overlay(as.array(fit), pars = "theta[1]") + ggtitle(paste("MCMC Density Overlay (alpha =", alpha, ")"))
  ggsave(paste0("mcmc_density_plot_alpha_", alpha, ".png"), plot = mcmc_density_plot)
}

# Plot results for each alpha
for (alpha in alpha_values) {
  fit <- fits[[paste0("alpha_", alpha)]]
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  plot_results(fit, counts, alpha, base_distribution)
}
