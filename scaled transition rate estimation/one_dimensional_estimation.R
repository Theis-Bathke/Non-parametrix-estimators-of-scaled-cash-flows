aalen_johansen_1d<-function(sim_c,rho)
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
# increments of the jump-processes
increments <- list()
increments[[1]] <- out[[1]]
for (i in 2:length(out)) {
  increments[[i]] <- out[[i]] - out[[i - 1]]
}
# scaled increments of the jump-processes
increments_scaled <- list()
increments_scaled[[1]] <- out_scaled[[1]]
for (i in 2:length(out)) {
  increments_scaled[[i]] <- out_scaled[[i]] - out_scaled[[i - 1]]
}
# construction of the increments for Lambda
contribution_first_scaled <- increments_scaled[[1]] / I0
# handling 0/0
contribution_first_scaled[is.nan(contribution_first_scaled)] <- 0
# singular summands in the sum for the definition of the estimator
# of Lambda
contributions_scaled <- mapply(FUN = function(a, b) {
  res <- b / a
  res[is.nan(res)] <- 0
  return(res)
}, It_scaled[-length(It_scaled)], increments_scaled[-1], SIMPLIFY = FALSE)
# cummulating these singular summands of Lambda leads to the whole development
# of Lambda
cumsums_scaled <- list()
cumsums_scaled[[1]] <- contribution_first_scaled
for (i in 2:length(out_scaled)) {
  cumsums_scaled[[i]] <- contributions_scaled[[i - 1]] + cumsums_scaled[[i - 1]]
}

contribution_first <- increments[[1]] / I0
# exception handling
contribution_first[is.nan(contribution_first)] <- 0
# singular summands in the sum for the definition of the estimator
# of Lambda
contributions <- mapply(FUN = function(a, b) {
  res <- b / a
  res[is.nan(res)] <- 0
  return(res)
}, It[-length(It)], increments[-1], SIMPLIFY = FALSE)
# cummulating these singular summands of Lambda leads to the whole development
# of Lambda
cumsums <- list()
cumsums[[1]] <- contribution_first
for (i in 2:length(out)) {
  cumsums[[i]] <- contributions[[i - 1]] + cumsums[[i - 1]]
}

# solving the integral-equation by discretization
sc_rowSums <- function(D_1,D_2){
  x=sum(D_1[1,c(2,3,4)])
  x=c(x,sum(D_2[2,c(5,6)]))
  x=c(x,0,0,0,0)
  return(x)
}
aj_scaled <- list()
aj_scaled[[1]] <- I0
Delta_1 <- cumsums[[1]]
Delta_2 <- cumsums_scaled[[1]]
aj_scaled[[2]] <- aj_scaled[[1]] + as.vector(aj_scaled[[1]] %*% Delta_2) - aj_scaled[[1]] *
  sc_rowSums(Delta_1,Delta_2)
for (i in 2:length(out_scaled)) {
  Delta_1 <- contributions[[i - 1]]
  Delta_2 <- contributions_scaled[[i - 1]]
  aj_scaled[[i + 1]] <- aj_scaled[[i]] + as.vector(aj_scaled[[i]] %*% Delta_2) -
    aj_scaled[[i]] * sc_rowSums(Delta_1,Delta_2)
}
aj <- list()
aj[[1]] <- I0
Delta <- cumsums[[1]]
aj[[2]] <- aj[[1]] + as.vector(aj[[1]] %*% Delta) - aj[[1]] *
  rowSums(Delta)
for (i in 2:length(out)) {
  Delta <- contributions[[i - 1]]
  aj[[i + 1]] <- aj[[i]] + as.vector(aj[[i]] %*% Delta) -
    aj[[i]] * rowSums((Delta))
}

return(list(n=n,It=It_scaled,aj_scaled=aj_scaled,ordered_times=ordered_times,
            contribution_first_scaled=contribution_first_scaled,
            contributions_scaled=contributions_scaled,
            It_scaled=It_scaled,increments_scaled=increments_scaled))
}

