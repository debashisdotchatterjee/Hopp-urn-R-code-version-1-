# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

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

# Simulate data for each alpha
color_counts <- list()
for (alpha in alpha_values) {
  color_counts[[paste0("alpha_", alpha)]] <- simulate_hoppe_urn(alpha, base_distribution, total_draws)
}


####################################


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


####################################


# Compile the Stan model
stan_model <- stan_model("dirichlet_multinomial.stan")

# Function to fit the model using MCMC
fit_model <- function(counts, alpha) {
  data_list <- list(N = sum(counts), K = length(counts), counts = counts, alpha = rep(alpha, length(counts)))
  fit <- sampling(stan_model, data = data_list, iter = 2000, chains = 4)
  return(fit)
}

# Fit the model for each alpha
fits <- list()
for (alpha in alpha_values) {
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  fits[[paste0("alpha_", alpha)]] <- fit_model(counts, alpha)
}


#########################################


# Function to plot prior, observed, and posterior density
plot_results <- function(fit, counts, alpha) {
  posterior_samples <- extract(fit)$theta
  posterior_density <- apply(posterior_samples, 2, mean)

  prior_dist <- data.frame(color = base_distribution, count = rep(alpha, length(base_distribution)))
  observed_dist <- data.frame(color = base_distribution, count = counts)
  posterior_dist <- data.frame(color = base_distribution, density = posterior_density)

  # Plot prior distribution
  prior_plot <- ggplot(data = prior_dist, aes(x = color, y = count, fill = color)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Prior Distribution (alpha =", alpha, ")")) +
    xlab("Color") +
    ylab("Count") +
    theme_minimal()
  ggsave(paste0("prior_distribution_alpha_", alpha, ".png"), plot = prior_plot)

  # Plot observed distribution
  observed_plot <- ggplot(data = observed_dist, aes(x = color, y = count, fill = color)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Observed Distribution (alpha =", alpha, ")")) +
    xlab("Color") +
    ylab("Count") +
    theme_minimal()
  ggsave(paste0("observed_distribution_alpha_", alpha, ".png"), plot = observed_plot)

  # Plot posterior density
  posterior_plot <- ggplot(data = posterior_dist, aes(x = color, y = density, fill = color)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Posterior Density (alpha =", alpha, ")")) +
    xlab("Color") +
    ylab("Density") +
    theme_minimal()
  ggsave(paste0("posterior_density_alpha_", alpha, ".png"), plot = posterior_plot)
}

# Plot results for each alpha
for (alpha in alpha_values) {
  fit <- fits[[paste0("alpha_", alpha)]]
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  plot_results(fit, counts, alpha)
}
