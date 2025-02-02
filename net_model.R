library(data.table)
library(ggplot2)


# Create network. This part is really trash, I did it fast
n_days <- 20
n_ind <- 15
minimum_observed_per_day <- 5
# Parameters
individual_effects <- rnorm(15, 0, 2)
alpha <- 0.2 # baseline probability of being synchronized
beta <- 2.5 # effect (in log odds) of increasing dominance of one individual by one
  
# Simulate data
individuals_present <- c()
n_ind_obs <- c()
for(day in 1:n_days) {
  ind_per_day <- sample(minimum_observed_per_day:n_ind, 1)
  individuals_present <- c(individuals_present, 
                           sample(n_ind, ind_per_day))
  n_ind_obs <- c(n_ind_obs, ind_per_day)
}
day_id <- c()
for(day in 1:n_days) {
  day_id <- c(day_id, 
              rep(day, n_ind_obs[day]))
}
who_when_where <- data.table(id = individuals_present,
                             day = day_id
)[order(day, id)
][, ":="(location = sample(1:2, 
                           .N, 
                           replace = T), # 0.5 probability on average of being in same location
         dominance = 1/id)]
data_full <- data.table(day = numeric(), 
                        id_1 = numeric(),  
                        id_2 = numeric(), 
                        location_1 = numeric(), 
                        location_2 = numeric(), 
                        dominance_1 = numeric(), 
                        dominance_2 = numeric())
for(day_i in 1:n_days) {
  print(day_i)
  
  day_net <- as.data.table(
    expand.grid(id_1 = who_when_where[day == day_i]$id,
                id_2 = who_when_where[day == day_i]$id) |>
      dplyr::left_join(who_when_where[day == day_i] |>
                  setnames(old = c("id", "location", "dominance"),
                           new = c("id_1", "location_1", "dominance_1"))) |>
      dplyr::left_join(who_when_where[day == day_i, 
                               c("id", "location", "dominance")] |>
                  setnames(old = c("id", "location", "dominance"),
                           new = c("id_2", "location_2", "dominance_2")))
  )[, c( "day", "id_1", "id_2", "location_1", "location_2", "dominance_1", "dominance_2")]
  
  data_full <- rbind(data_full, day_net)
}
data_full_m <- data_full[id_1 != id_2
][, ":="(synchrony = ifelse(location_1 != location_2, 
                            log(0.2/(1-0.2)),
                            log(0.4/(1-0.4)) + individual_effects[id_1] + individual_effects[id_2] + beta * (dominance_1 + dominance_2)))
][, ":="(p = exp(synchrony)/(1 + exp(synchrony))) # inverse logit function
][, ":="(observed_hours = rbinom(.N, 24, prob = p))
][, ":="(edge_id = apply(.SD[, c("id_1", "id_2")], 
                         1, 
                         FUN = function(x) {
                           paste(sort(x), collapse = "")
                         }))
][order(edge_id)
][, ":="(pair_id = rleid(edge_id),
         by = edge_id)
][, .SD[1], by = list(day, pair_id)
][, ":="(same_location = ifelse(location_1 == location_2,
                                T,
                                F))
  ][, c("day", "id_1", "id_2", "pair_id", "location_1", "location_2", "same_location", "p", "synchrony", "dominance_1", "dominance_2", "observed_hours")]

write.csv(data_full_m[, c("day", "id_1", "id_2", "dominance_1", "dominance_2", "pair_id", "observed_hours", "same_location")],
          "C:/Users/marco/Desktop/mio/programming_projects/WWCS/data/observation_data.csv",
          row.names = F)
