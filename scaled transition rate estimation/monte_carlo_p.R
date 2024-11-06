monte_carlo_p<-function(sim_c,rho)
{
  #addition of 3' and 4' to sim
  for ( i in 1:length(sim_c)){
    if(length(sim_c[[i]]$states)==3){
      if(sim_c[[i]]$states[3]!=2){
        sim_c[[i]]$states[3]<-sim_c[[i]]$states[3]+2
      }
    }
  }
  
  n <- length(sim_c)
  # number of states minus 1 
  p <- max(unique(unlist(lapply(sim_c, function(Z) Z$states)))) - 1
  # Time-points of last jump of every individual (censoring times?)
  r_times <- unlist(lapply(sim_c, FUN = function(Z) tail(Z$times, 1)))
  # all jump-times without the jump into the first state
  t_pool <- unlist(lapply(sim_c, FUN = function(Z) Z$times[-1]))
  # gives a vector that finds the individuals corresponding to every jump
  individuals <- c()
  for (i in 1:n) {
    individuals <- c(individuals, rep(i, length(sim_c[[i]]$times[-1])))
  }
  #calculating the scaling factor for every individual
  scalingfactor <- c()
  for (i in 1:n) {
    scalingfactor <- c(scalingfactor,rep(ifelse(sim_c[[i]]$states[2]==2,rho(sim_c[[i]]$times[2]),1),length(sim_c[[i]]$times)-1))
  }
  # List of all the jumps with initial state latter state
  jumps_pool <- matrix(NA, 0, 2)
  for (i in 1:n) {
    v <- sim_c[[i]]$states
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
  #gives the scaling factors corresponding to the ordered jumps
  ordered_scalings <- c(NA,scalingfactor[order_of_times])
  # gives the corresponding initial and last state from the ordered jumps
  ordered_jumps <- jumps_pool[order_of_times, ]
  #catches the censoring times:
  R_times <- unlist(lapply(sim_c, FUN = function(Z) tail(Z$times, 1)))
  #safes the time-points that are a censoring time or "finish a path"
  decisions <- ordered_times %in% R_times
  # gives a Vector that shows the change in the counting processes
  ordered_n <- rep(list(matrix(0, p + 1, p + 1)), nrow(ordered_jumps))
  for (i in 1:nrow(ordered_jumps)) {
    ordered_n[[i]][ordered_jumps[i, ][1], ordered_jumps[i, ][2]] <- 1
  }
  #constructs a list of scaled jumps
  ordered_n_scaled <- rep(list(matrix(0, p + 1, p + 1)), nrow(ordered_jumps))
  for (i in 1:nrow(ordered_jumps)) {
    ordered_n_scaled[[i]][ordered_jumps[i, ][1], ordered_jumps[i, ][2]] <- ordered_scalings[i+1]
  }
  #catches jumps that do not happen because of censoring
  ordered_n <- lapply(ordered_n, FUN = function(Z) Z - diag(diag(Z)))
  ordered_n_scaled <- lapply(ordered_n_scaled, FUN = function(Z) Z - diag(diag(Z)))
  # the difference of the indicator-function, this can be used to
  # calculate the amount of people in the 4 states, by
  # adding the corresponding vectors together
  colsums_of_n <- lapply(ordered_n, function(N) {
    colSums(N - t(N))
  })
  #here we need to change the formula for the differences in the indicator functions because the formula
  #in theory changes. In general we would need a vector which indicates which states are free policy 
  #and which are not. We hard-code this here so there is no need for that.
  colsums_of_n_scaled<- list()

  for (t in 1:(length(ordered_times)-1)){
    colsums_of_n_scaled[[t]]<-colSums(ordered_n_scaled[[t]] - 
                                        t(ordered_n_scaled[[t]]))
    if (colsums_of_n[[t]][1]==-1){
      colsums_of_n_scaled[[t]][1]<- -1}
  }
  
  out <- out2 <- list()
  # cumulated jumps from initial states to final states
  out[[1]] <- ordered_n[[1]] * n^(-1)
  for (tm in 2:(length(ordered_times) - 1)) {
    out[[tm]] <- out[[tm - 1]] + ordered_n[[tm]] * n^(-1)
  }
  out_scaled <- out2_scaled <- list()
  # cumulated scaled jumps from initial states to final states
  out_scaled[[1]] <- ordered_n_scaled[[1]] * n^(-1)
  for (tm in 2:(length(ordered_times) - 1)) {
    out_scaled[[tm]] <- out_scaled[[tm - 1]] + ordered_n_scaled[[tm]] * n^(-1)
  }
  # vector of differences of the scaled estimators for every jump-time.
  # Censoring is included
  out2[[1]] <- colsums_of_n[[1]] * n^(-1)
  if (decisions[2]) {
    wch <- ordered_individuals[2]
    end_state <- as.numeric(1:(p + 1) == tail(sim_c[[wch]]$states, 1))
    out2[[1]] <- out2[[1]] - end_state * n^(-1)
  }
  for (tm in 2:(length(ordered_times) - 1)) {
    out2[[tm]] <- out2[[tm - 1]] + colsums_of_n[[tm]] * n^(-1)
    if (decisions[tm+1]) {
      wch <- ordered_individuals[tm+1]
      end_state <- as.numeric(1:(p + 1) == tail(sim_c[[wch]]$states, 1))
      out2[[tm]] <- out2[[tm]] - end_state * n^{
        -1
      }
    }
    
  }
  out2_scaled[[1]] <- colsums_of_n_scaled[[1]] * n^(-1)
  if (decisions[2]) {
    wch <- ordered_individuals[2]
    end_state <- as.numeric(1:(p + 1) == tail(sim_c[[wch]]$states, 
                                              1))
    out2_scaled[[1]] <- out2_scaled[[1]] - end_state * n^{-1} * ordered_scalings[2]
  }
  for (tm in 2:(length(ordered_times) - 1)) {
    out2_scaled[[tm]] <- out2_scaled[[tm - 1]] + colsums_of_n_scaled[[tm]] * n^(-1)
    if (decisions[tm+1]) {
      wch <- ordered_individuals[tm+1]
      end_state <- as.numeric(1:(p + 1) == tail(sim_c[[wch]]$states, 
                                                1))
      out2_scaled[[tm]] <- out2_scaled[[tm]] - end_state * n^{-1} * ordered_scalings[tm + 1]
    }
  }
  # Vector of Starting values of the indicator vector for every insured
  # is only really needed if not all insured start in the correct state
  I_initial <- lapply(sim_c, FUN = function(Z) {
    as.numeric(1:(p + 1) == head(Z$states, 1))
  })
  # Starting value for the estimator (connected to the previous definition)
  I0 <- Reduce("+", I_initial) * n^(-1)
  # Constructing the estimator for the censored state occupation probabilities
  It <- lapply(out2, FUN = function(N) {
    return(I0 + N)
  })
  # Constructing the scaled estimator for the censored state occupation probabilities
  It_scaled <- lapply(out2_scaled, FUN = function(N) {
    return(c(I0) + N)
  })


  return(list(n=n,It=It_scaled,ordered_times=ordered_times))
}

