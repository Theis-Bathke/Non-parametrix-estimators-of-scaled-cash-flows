monte_carlo_2d<-function(sim){

n <- length(sim)
# number of states minus 1 ???
p <- max(unique(unlist(lapply(sim, function(Z) Z$states)))) - 1
# Time-points of last jump of every individual (censoring times?)
r_times <- unlist(lapply(sim, FUN = function(Z) tail(Z$times, 1)))
# all jump-times without the jump into the first state
t_pool <- unlist(lapply(sim, FUN = function(Z) Z$times[-1]))
# gives a vector that finds the individuals corresponding to every jump
individuals <- c()
for (i in 1:n) {
  individuals <- c(individuals, rep(i, length(sim[[i]]$times[-1])))
}
# List of all the jumps with initial state latter state
jumps_pool <- matrix(NA, 0, 2)
for (i in 1:n) {
  v <- sim[[i]]$states
  jumps_pool <- rbind(jumps_pool, cbind(
    rev(rev(v)[-1]),
    v[-1]
  ))
}
# indexes for the sorting of all the jump-times
order_of_times <- order(t_pool)
# ordered jump-times with added 0 at the start
ordered_times <- c(40, t_pool[order_of_times])
# gives the individual that corresponds to the ordered jumps
ordered_individuals <- c(NA, individuals[order_of_times])
# gives the corresponding initial and last state from the ordered jumps
ordered_jumps <- jumps_pool[order_of_times, ]

index <- which(ordered_times > 65)[1]

rho <- function(t) {
  if (t >= 65) {
    return(1)
  } else if (t < 40) {
    return(0)
  } else {
    return((V_1_plus(t) - V_1_minus(t)) / V_1_plus(t))
  }
}



reserve_mc_cont<-rep(0, length(ordered_times))
reserve_mc_cont[1] <- -b_0
for (i in 1:n) {
  x <- sim[[i]]
  if (length(x$times) == 2) {
    if (x$states[2] == 3) {
      timei <- findInterval(x$times[2], ordered_times)
      for (j in 2:(timei )) {
        reserve_mc_cont[j] <- reserve_mc_cont[j] - 1 / n *
          premium * interest_f(ordered_times[j - 1], ordered_times[j])
      }
      reserve_mc_cont[timei] <- reserve_mc_cont[timei] + 1 / n * (V_1_plus(x$times[2])
      - V_1_minus(x$times[2])) * 1 / (1 + interest / 100)^(x$times[2] - 40)
    }
    if (x$states[2] == 4) {
      timei <- findInterval(x$times[2], ordered_times)
      for (j in 2:min(index, (timei ))) {
        reserve_mc_cont[j] <- reserve_mc_cont[j] - 1 / n *
          premium * interest_f(ordered_times[j - 1], min(65, ordered_times[j]))
      }
      for (j in index:max(index, (timei ))) {
        reserve_mc_cont[j] <- reserve_mc_cont[j] + 1 / n *
          pension * interest_f(max(65, ordered_times[j - 1]), ordered_times[j])
      }
    }
  }
  if (length(x$times) == 3) {
    if (x$states[3] == 3) {
      timei <- findInterval(x$times[2], ordered_times)
      for (j in 2:(timei )) {
        reserve_mc_cont[j] <- reserve_mc_cont[j] - 1 / n *
          premium * interest_f(ordered_times[j - 1], ordered_times[j])
      }
      times <- findInterval(x$times[3], ordered_times)
      reserve_mc_cont[times ] <- reserve_mc_cont[times ] + 1 / n * rho(x$times[2]) *
        V_1_plus(x$times[3]) * 1 / (1 + interest / 100)^(x$times[3] - 40)
    }


    if (x$states[3] == 4) {
      timei <- findInterval(x$times[2], ordered_times)
      for (j in 2:(timei )) {
        reserve_mc_cont[j] <- reserve_mc_cont[j] - 1 / n *
          premium * interest_f(ordered_times[j - 1], ordered_times[j])
      }
      timed<-findInterval(x$times[3], ordered_times)
      for (j in index:max(index, (timed ))) {
        reserve_mc_cont[j] <- reserve_mc_cont[j] + 1 / n * rho(x$times[2]) *
          pension * interest_f(max(65, ordered_times[j - 1]), ordered_times[j])
      }
    }
  }

}
reserve_mc_recover <- cumsum(reserve_mc_cont)
return(list(reserve_cont=reserve_mc_recover,ordered_times=ordered_times))
}