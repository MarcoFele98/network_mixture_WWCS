library(data.table)
library(ggplot2)
library(ggraph)
library(tidygraph)

options(scipen = 10)
theme_set(cowplot::theme_cowplot())

# Read data file
data_synchrony <- read.csv("data/baboon_data.csv") |>
  as.data.table()

ggplot(data_synchrony,
       aes(dominance_1 + dominance_2, observed_hours)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~day)

ggplot(data_synchrony[, .(n_obs = .N), by = pair_id]) +
  geom_point(aes(pair_id, n_obs))

# Prepare STAN data
stan_data <- list(N = nrow(data_synchrony),
                  N_ind = length(unique(c(data_synchrony$id_1, 
                                        data_synchrony$id_2))),
                  n_hours_same = data_synchrony$observed_hours,
                  id_1 = data_synchrony$id_1,
                  id_2 = data_synchrony$id_2,
                  dominance_1 = data_synchrony$dominance_1,
                  dominance_2 = data_synchrony$dominance_2)

# Compile 
network_model_stan <- rstan::stan_model("scripts/network_model.stan")

# Sample 
posterior <- rstan::sampling(network_model_stan,
                             data = stan_data,
                             cores = getOption("mc.cores", 4L), # parallelize
                             chains = 4,
                             iter = 4000)

posterior_draws <- as.data.table(posterior)

# MCMC diagnostics
bayesplot::rhat(posterior) # should be close to 1
bayesplot::mcmc_trace(posterior_draws[, c("beta", "theta", "alpha[1]", "alpha[2]")]) |>
  ggsave(filename = "figures/MCMC_diagnostics.png",
         width = 7,
         height = 5,
         bg = "white")

# Visualize results
(p1 <- ggplot(posterior_draws) +
  geom_histogram(aes(x = beta,
                     y=..count../sum(after_stat(count)))) +
  geom_vline(aes(xintercept = 2.5), 
             color = "red", linewidth = 1) +
  ylab("Posterior probability") +
  xlab("Beta"))

(p2 <- ggplot(posterior_draws) +
  geom_histogram(aes(x = theta,
                     y=..count../sum(after_stat(count)))) +
  geom_vline(aes(xintercept = 0.5), 
             color = "red", linewidth = 1) +
  xlim(c(0,1)) +
  xlab("Theta") +
  ylab("Posterior probability"))

(p3 <- ggplot(posterior_draws[, c("alpha[1]")] |>
         setnames(old = c("alpha[1]"),
                  new = c("alpha_1"))) +
  geom_histogram(aes(x = alpha_1,
                     y=..count../sum(after_stat(count)))) +
  geom_vline(aes(xintercept = log(0.2/(1-0.2))), 
             color = "red", linewidth = 1) +
  xlab("Alpha_1") +
  ylab("Posterior probability"))

(p4 <- ggplot(posterior_draws[, c("alpha[2]")] |>
         setnames(old = c("alpha[2]"),
                  new = c("alpha_2"))) +
  geom_histogram(aes(x = alpha_2,
                     y=..count../sum(after_stat(count)))) +
  geom_vline(aes(xintercept = log(0.4/(1-0.4)) / 2 + mean((1:15)/15)), # divided by two bc I generated data with one pair level intercept, which is equivalent. Also convert to centralized scale 
             color = "red", linewidth = 1) +
  xlab("Alpha_2") +
  ylab("Posterior probability"))

fig_1 <- ggpubr::ggarrange(p1, p2, p3, p4)
ggsave("figures/posterior.png",
       fig_1,
       width = 7,
       height = 5,
       bg = "white")

# see accuracy of location prediction
check_preds <- copy(data_full_m)[, ":="(prob_same = posterior_draws[, grepl('prob_same', names(posterior_draws)), with = F
                                                                    ][, unlist(lapply(.SD, median))],
                                        sd_prob_same = posterior_draws[, grepl('prob_same', names(posterior_draws)), with = F
                                        ][, unlist(lapply(.SD, sd))],
                                        same = location_1 == location_2)]

ggplot(check_preds) +
  geom_boxplot(aes(same, prob_same)) +
  geom_jitter(aes(same, prob_same, color = sd_prob_same)) +
  scale_color_viridis_c(name = "Uncertainty",
                        direction = -1) +
  xlab("Same location") +
  ylab("Posterior probability same location")

ggsave("figures/location_accuracy.png",
       width = 7,
       height = 5,
       bg = "white")

# plot random effects
check_individual <- data.table(
  pred_effect = posterior_draws[, grepl('transformed_id', names(posterior_draws)), with = F
][, unlist(lapply(.SD, median))],
sd_pred_effect = posterior_draws[, grepl('transformed_id', names(posterior_draws)), with = F
][, unlist(lapply(.SD, sd))],
real_effect = individual_effects,
id = 1:n_ind,
n_obs = as.numeric(table(c(data_synchrony$id_1, data_synchrony$id_2)))
)

ggplot(check_individual) +
  geom_point(aes(id, real_effect, size = n_obs)) +
  geom_pointrange(aes(x = id, y = pred_effect,
                      ymin = pred_effect - sd_pred_effect,
                      ymax = pred_effect + sd_pred_effect), 
                  color = "red") +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  xlab("Individual") +
  ylab("Effect") +
  guides(size = "none")

ggsave("figures/random_effects.png",
       width = 5,
       height = 5,
       bg = "white")

# plot edges 
edges_median <- posterior_draws[, grepl('edges', names(posterior_draws)), 
                                with = F
                                ][ , unlist(lapply(.SD, median))]
edges_sd <- posterior_draws[, grepl('edges', names(posterior_draws)), 
                            with = F
                            ][ , unlist(lapply(.SD, sd))]

matrix_edges <- matrix(ncol = n_ind, nrow = n_ind)
matrix_edges_sd <- matrix(ncol = n_ind, nrow = n_ind)
index <- 1
index_2 <- 1
for(i in 1:length(edges_median)) {
  matrix_edges[index, index_2] <- edges_median[i]
  matrix_edges_sd[index, index_2] <- edges_sd[i]
  index <- index + 1
  if(index == 16) {
    index_2 <- index_2 + 1
    index <- 1
  }
}

net <- as_tbl_graph(matrix_edges, 
                    directed = F) |>
  activate(nodes) |>
  mutate(id = 1:15,
         dominance = 15:1/15) |>
  activate(edges) |>
  left_join(as_tbl_graph(matrix_edges_sd, 
                         directed = F) |>
              activate(edges) |>
              rename(sd = weight) |>
              as.data.frame(edges))

ggraph(net, layout = "stress") +
  geom_edge_link(aes(color = weight, 
                     alpha = sd),
                 linewidth = 1) +
  geom_node_point(aes(size = dominance)) +
  scale_edge_color_viridis(option = "viridis") +
  guides(size = "none",
         alpha = "none")

ggsave("figures/network.png",
       width = 5,
       height = 5,
       bg = "white")

# Plot centrality
data_centrality <- data.table(centrality = posterior_draws[, grepl('centrality', names(posterior_draws)),  with = F
                                                           ][ , unlist(lapply(.SD, median))],
                              centrality_sd = posterior_draws[, grepl('centrality', names(posterior_draws)),  with = F
                                                              ][ , unlist(lapply(.SD, sd))],
                              id = 1:15,
                              dominance = 15:1/15)


ggplot(data_centrality) +
  geom_pointrange(aes(x = dominance, y = centrality,
                      ymin = centrality - centrality_sd, ymax = centrality + centrality_sd))

ggsave("figures/centrality.png",
       width = 5,
       height = 5,
       bg = "white")

