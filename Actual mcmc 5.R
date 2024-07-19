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
  real<lower=0> prior_alpha;    // Prior on alpha
}
parameters {
  simplex[K] theta;             // Proportions for each color
  real<lower=0> alpha;          // Concentration parameter
}
model {
  alpha ~ normal(prior_alpha, 1);  // Prior on alpha
  theta ~ dirichlet(rep_vector(alpha, K)); // Prior on theta
  counts ~ multinomial(theta);  // Likelihood
}
"

# Compile the Stan model from the string
stan_model <- stan_model(model_code = stan_model_code)

# Function to fit the model using MCMC
fit_model <- function(counts, prior_alpha) {
  data_list <- list(N = sum(counts), K = length(counts), counts = counts, prior_alpha = prior_alpha)
  fit <- sampling(stan_model, data = data_list, iter = 4000, warmup = 2000, chains = 4, control = list(adapt_delta = 0.95))
  return(fit)
}

# Fit the model for each alpha
fits <- list()
for (alpha in alpha_values) {
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  fits[[paste0("alpha_", alpha)]] <- fit_model(counts, prior_alpha = alpha)
}

# Function to plot prior, observed, and posterior density
plot_results <- function(fit, counts, alpha, base_distribution) {
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
  mcmc_trace_plot <- mcmc_trace(as.array(fit), pars = "alpha") + ggtitle(paste("MCMC Trace Plot (alpha =", alpha, ")"))
  ggsave(paste0("mcmc_trace_plot_alpha_", alpha, ".png"), plot = mcmc_trace_plot)

  mcmc_density_plot <- mcmc_dens_overlay(as.array(fit), pars = "alpha") + ggtitle(paste("MCMC Density Overlay (alpha =", alpha, ")"))
  ggsave(paste0("mcmc_density_plot_alpha_", alpha, ".png"), plot = mcmc_density_plot)
}

# Plot results for each alpha
for (alpha in alpha_values) {
  fit <- fits[[paste0("alpha_", alpha)]]
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  plot_results(fit, counts, alpha, base_distribution)
}


############################

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(rstan)
library(bayesplot)

# Set seed for reproducibility
set.seed(123)

# Define the alpha values to be tested
alpha_values <- c(1, 2, 3, 5)

# Define the base distribution (colors)
base_distribution <- c("red", "blue", "green", "yellow", "orange", "purple", "pink", "brown", "black", "grey")

# Define a color palette for consistency
color_palette <- c("red", "blue", "green", "yellow", "orange", "purple", "pink", "brown", "black", "grey")

# Function to simulate Hoppe urn process
simulate_hoppe_urn <- function(alpha, num_draws = 5000) {
  urn <- rep("mutator", alpha)
  color_counts <- data.frame(color = base_distribution, count = 0, stringsAsFactors = FALSE)

  for (n in 1:num_draws) {
    if (runif(1) < alpha / (alpha + n)) {
      new_color <- sample(base_distribution, 1)
      urn <- c(urn, new_color)
      color_counts$count[color_counts$color == new_color] <- color_counts$count[color_counts$color == new_color] + 1
    } else {
      drawn_ball <- sample(urn, 1)
      urn <- c(urn, drawn_ball)
      if (drawn_ball != "mutator") {
        color_counts$count[color_counts$color == drawn_ball] <- color_counts$count[color_counts$color == drawn_ball] + 1
      }
    }
  }

  return(color_counts)
}

# Function to run MCMC for posterior density
run_mcmc <- function(color_counts, alpha_prior) {
  # Define the Stan model
  stan_model_code <- "
  data {
    int<lower=1> K;  // number of colors
    int<lower=0> y[K];  // observed counts of each color
    real<lower=0> alpha_prior;  // prior for alpha
  }
  parameters {
    simplex[K] theta;  // proportion of each color
    real<lower=0> alpha;  // concentration parameter
  }
  model {
    alpha ~ gamma(2, 1 / alpha_prior);  // prior on alpha
    y ~ multinomial(theta);  // likelihood
  }
  "

  # Compile the Stan model
  stan_model <- stan_model(model_code = stan_model_code)

  # Data for Stan model
  stan_data <- list(
    K = nrow(color_counts),
    y = color_counts$count,
    alpha_prior = alpha_prior
  )

  # Run Stan model
  fit <- sampling(stan_model, data = stan_data, iter = 2000, chains = 4)

  return(fit)
}

# Function to plot results
plot_results <- function(fit, counts, alpha, base_distribution) {
  # Extract posterior samples
  posterior_samples <- extract(fit)

  # Convert samples to data frame for plotting
  posterior_data <- data.frame(posterior_samples$theta)
  colnames(posterior_data) <- base_distribution

  # Plot observed counts
  observed_plot <- ggplot(data = data.frame(color = base_distribution, count = counts), aes(x = color, y = count, fill = color)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Observed Distribution (alpha =", alpha, ")")) +
    xlab("Color") +
    ylab("Count") +
    theme_minimal() +
    scale_fill_manual(values = color_palette)

  ggsave(paste0("observed_distribution_alpha_", alpha, ".png"), plot = observed_plot)

  # Plot posterior density
  posterior_plot <- ggplot() +
    geom_density(data = melt(posterior_data), aes(x = value, fill = variable), alpha = 0.5) +
    #geom_vline(xintercept = alpha, color = "red", linetype = "dashed", size = 1) +
    ggtitle(paste("Posterior Density color bin (alpha =", alpha, ")")) +
    xlab("Color Count") +
    ylab("Density") +
    theme_minimal() +
    scale_fill_manual(values = color_palette)

  ggsave(paste0("posterior_density_alpha_", alpha, ".png"), plot = posterior_plot)

  # Plot MCMC trace
  trace_plot <- mcmc_trace(as.array(fit), pars = "alpha") +
    ggtitle(paste("MCMC Trace Plot (alpha =", alpha, ")")) +
    theme_minimal()

  ggsave(paste0("mcmc_trace_plot_alpha_", alpha, ".png"), plot = trace_plot)

  # Plot MCMC density
  density_plot <- mcmc_dens(as.array(fit), pars = "alpha") +
    geom_vline(xintercept = alpha, color = "red", linetype = "dashed", size = 1) +
    ggtitle(paste("MCMC Density Overlay (alpha =", alpha, ")")) +
    theme_minimal()

  ggsave(paste0("mcmc_density_plot_alpha_", alpha, ".png"), plot = density_plot)
}

# Main simulation and plotting process
fits <- list()
color_counts <- list()
for (alpha in alpha_values) {
  cat("Simulating for alpha =", alpha, "\n")
  color_counts[[paste0("alpha_", alpha)]] <- simulate_hoppe_urn(alpha)
  fits[[paste0("alpha_", alpha)]] <- run_mcmc(color_counts[[paste0("alpha_", alpha)]], alpha)
}

# Plot results for each alpha
for (alpha in alpha_values) {
  fit <- fits[[paste0("alpha_", alpha)]]
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  plot_results(fit, counts, alpha, base_distribution)
}


####################

library(gridExtra)

# Function to plot prior vs posterior densities for each color bin
plot_prior_vs_posterior <- function(fit, prior_counts, alpha, base_distribution) {
  # Extract posterior samples
  posterior_samples <- extract(fit)

  # Convert samples to data frame for plotting
  posterior_data <- data.frame(posterior_samples$theta)
  colnames(posterior_data) <- base_distribution

  # Create a data frame for prior counts
  prior_data <- data.frame(color = base_distribution, count = prior_counts)
  prior_data <- prior_data %>% mutate(proportion = count / sum(count))

  # Create plots for each color bin
  plot_list <- lapply(base_distribution, function(color) {
    ggplot() +
      geom_density(data = posterior_data, aes_string(x = color, fill = color), alpha = 0.5) +
      geom_vline(xintercept = prior_data$proportion[prior_data$color == color], color = "red", linetype = "dashed", size = 1) +
      ggtitle(paste("Prior vs Posterior Density for", color, "(alpha =", alpha, ")")) +
      xlab("Proportion") +
      ylab("Density") +
      theme_minimal() +
      scale_fill_manual(values = setNames(rep(color_palette, each = length(base_distribution)), base_distribution))
  })

  # Arrange plots in a grid
  grid_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))

  # Save the grid plot
  ggsave(paste0("prior_vs_posterior_density_alpha_", alpha, ".png"), plot = grid_plot, width = 14, height = 10)
}

# Generate and save plots for each alpha
for (alpha in alpha_values) {
  fit <- fits[[paste0("alpha_", alpha)]]
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  plot_prior_vs_posterior(fit, rep(alpha, length(base_distribution)), alpha, base_distribution)
}


###################


library(ggplot2)
library(gridExtra)

# Function to plot prior vs posterior densities for each color bin
plot_prior_vs_posterior <- function(fit, prior_counts, alpha, base_distribution) {
  # Extract posterior samples
  posterior_samples <- extract(fit)

  # Convert samples to data frame for plotting
  posterior_data <- data.frame(posterior_samples$theta)
  colnames(posterior_data) <- base_distribution

  # Create a data frame for prior counts
  prior_data <- data.frame(color = base_distribution, count = prior_counts)
  prior_data <- prior_data %>% mutate(proportion = count / sum(count))

  # Create plots for each color bin
  plot_list <- lapply(base_distribution, function(color) {
    ggplot() +
      geom_density(data = posterior_data, aes_string(x = color), fill = color, alpha = 0.5) +
      geom_vline(xintercept = prior_data$proportion[prior_data$color == color], color = "red", linetype = "dashed", size = 1) +
      ggtitle(paste("Prior vs Posterior Density for", color, "(alpha =", alpha, ")")) +
      xlab("Proportion") +
      ylab("Density") +
      theme_minimal() +
      scale_fill_manual(values = setNames(rep(color_palette, each = length(base_distribution)), base_distribution))
  })

  # Arrange plots in a grid
  grid_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))

  # Save the grid plot
  ggsave(paste0("prior_vs_posterior_density_alpha_", alpha, ".png"), plot = grid_plot, width = 14, height = 10)
}

# Generate and save plots for each alpha
for (alpha in alpha_values) {
  fit <- fits[[paste0("alpha_", alpha)]]
  counts <- color_counts[[paste0("alpha_", alpha)]]$count
  plot_prior_vs_posterior(fit, rep(alpha, length(base_distribution)), alpha, base_distribution)
}
