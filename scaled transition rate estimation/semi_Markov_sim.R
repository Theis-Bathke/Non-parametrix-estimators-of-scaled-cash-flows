
sim_example<-function(n,censoring){
#n as the number of simulated individuals
#censoring is a vector of length n which 
#censoring<-runif(n,50,130) as an example for censoring

if (n!=length(censoring)){
  print("Censoring vector does not have the correct length")
  return()
}
mu_12 <- function(t, u) {
  ifelse(t <= 65, 0.1, 0)
}
mu_14_24 <- function(t, u) {
  0.0005 + 10^(5.728 - 10 + 0.038 * t)
}


mu_13 <- function(t, u) {
  ifelse(t <= 65, 0.05, 0)
}
mu_23 <- function(t, u) {
  ifelse(t <= 65, 0.05 + ifelse((1 / 2 <= u) * (u < 5 / 2), 0.2, 0), 0)
}

jump_rate <- function(i, t, u) {
  if (i == 1) {
    mu_14_24(t, u) + mu_12(t, u) + mu_13(t, u)
  } else if (i == 2) {
    mu_14_24(t, u) + mu_23(t, u)
  } else {
    0
  }
}

mark_dist <- function(i, s, v) {
  if (i == 1) {
    c(
      0, mu_12(s, v) / jump_rate(i, s, v),
      mu_13(s, v) / jump_rate(i, s, v), mu_14_24(s, v) / jump_rate(i, s, v)
    )
  } else if (i == 2) {
    c(
      0, 0, mu_23(s, v) / jump_rate(i, s, v),
      mu_14_24(s, v) / jump_rate(i, s, v)
    )
  } else {
    c(0, 0, 0, 0)
  }
}

bs_r <- c(jump_rate(1, 110, 0), jump_rate(2, 110, 0), 0, 0) # bound on rates
tn_r <- 110 # Time horizon
abs_r <- c(FALSE, FALSE, TRUE, TRUE) # absorbing states
t_r <- 40 # initial time
u_r <- 40 # initial duration
i_r <- 1 # initial state




sim <- list()
for (i in 1:n) {
  sim[[i]] <- AalenJohansen::sim_path(i_r,
                       rates = jump_rate, dists = mark_dist, t = t_r,
                       u = u_r, tn = censoring[i], abs = abs_r,
                       bs = bs_r
  )
}
return(sim)
}