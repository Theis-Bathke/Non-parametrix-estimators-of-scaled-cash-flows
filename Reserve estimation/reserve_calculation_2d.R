reserve_2d<-function(estimation,b_0,premium,interest,interest_f,pension,rho,V_1_plus,V_1_minus){
reserve <- numeric(length(estimation$ordered_times))
reserve[1] <- -b_0
rho_calc <- sapply(estimation$ordered_times, rho)
# Reserve calculation with no free policy option
for (t in 2:length(estimation$ordered_times)) {
  reserve[t] <- reserve[t - 1] + estimation$p[[t - 1]][1] *
    interest_f(estimation$ordered_times[t - 1], estimation$ordered_times[t]) *
    ifelse(estimation$ordered_times[t] <= 65, -premium, pension) * (1 + interest / 100)^(40)
}


# Calculation of the reserve including the free policy option
index <- which(estimation$ordered_times > 65)[1]

reserve_comp <- numeric(length(estimation$ordered_times))
reserve_comp[1] <- -b_0
# calculation of the reserve before free policy has
# an impact (transition payments regard surrender still need to be addressed)
reserve_comp[1] <- -b_0
for (t in 2:(index - 1)) {
  reserve_comp[t] <- reserve_comp[t - 1] + estimation$p[[t - 1]][1] *
    interest_f(estimation$ordered_times[t - 1], estimation$ordered_times[t]) * -premium
}
# transition from before pension to after pension ( no free policy option)
reserve_comp[index] <- reserve_comp[index - 1] + estimation$p[[index - 1]][1] *
  interest_f(estimation$ordered_times[index - 1], 65) * -premium + estimation$p[[index - 1]][1] *
  interest_f(65, estimation$ordered_times[index]) * pension
# reserves for pensioners with no free policy option
for (t in (index + 1):(length(estimation$ordered_times))) {
  reserve_comp[t] <- reserve_comp[t - 1] + estimation$p[[t - 1]][1] *
    interest_f(estimation$ordered_times[t - 1], estimation$ordered_times[t]) * pension
}
# second summand: free policy option and pension payments
# where the payments are  not in the integral
# Integral_t constructs the integral  over [40,65]times[40, t_k-1]
Integral_t <- 0
for (t1 in (2:(index - 1))) {
  Integral_temp <- 0
  for (t2 in (2:(index - 1))) {
    if (estimation$I_t_2d[1 + (t1 - 2) * 4, 1 + (t2 - 2) * 4]
    > 10^(-10)) {
      Integral_temp <- Integral_temp + estimation$p_11[t1 - 1, t2 - 1] *
        estimation$ordered_n_2d_matrix_1212[t1 - 1, t2 - 1] *
        1 / estimation$I_t_2d[1 + (t1 - 2) * 4, 1 + (t2 - 2) * 4] * 1 / estimation$n
    }
    if (estimation$I_t_2d[1 + (t1 - 2) * 4, 2 + (t2 - 2) * 4]
    > 10^(-10)) {
      Integral_temp <- Integral_temp + estimation$p_12[t1 - 1, t2 - 1] *
        1 / estimation$I_t_2d[1 + (t1 - 2) * 4, 2 + (t2 - 2) * 4] *
        (-estimation$ordered_n_2d_matrix_1223[t1 - 1, t2 - 1] -
          estimation$ordered_n_2d_matrix_1224[t1 - 1, t2 - 1]) * 1 / estimation$n
    }
  }
  Integral_t <- Integral_t + Integral_temp * rho_calc[t1] # the factor rho(u_1) is constant for all u_2, thus we can factor estimation$It
}

# calculation of the first free-policy option summand
# d_C_2 constructs the payments that stem from pension payments after the free policy option in the second state
d_C_2 <- numeric(c(length(estimation$ordered_times)))
d_C_2[index] <- Integral_t * pension * interest_f(65, estimation$ordered_times[index])

# reserve_comp[index] <- reserve_comp[index] + Integral_t * pension *
#  interest_f(65, estimation$ordered_times[index])
# calculation of the integral over the intervals [0,65]\times (65,index]
# d_I_t calculate the new integral that is added in every iterative step.
d_I_t <- numeric(max(c(length(estimation$ordered_times) - index + 1, 1)))
for (t1 in (2:(index - 1))) {
  if (estimation$I_t_2d[1 + (t1 - 2) * 4, 1 + (index - 2) * 4]
  > 10^(-10)) {
    d_I_t[1] <- d_I_t[1] + rho_calc[t1] * estimation$p_11[t1 - 1, index - 1] *
      estimation$ordered_n_2d_matrix_1212[t1 - 1, index - 1] *
      1 / estimation$I_t_2d[1 + (t1 - 2) * 4, 1 + (index - 2) * 4] * 1 / estimation$n
  }
  if (estimation$I_t_2d[1 + (t1 - 2) * 4, 2 + (index - 2) * 4]
  > 10^(-10)) {
    d_I_t[1] <- d_I_t[1] + rho_calc[t1] * estimation$p_12[t1 - 1, index - 1] *
      1 / estimation$I_t_2d[1 + (t1 - 2) * 4, 2 + (index - 2) * 4] *
      (-estimation$ordered_n_2d_matrix_1223[t1 - 1, index - 1] -
        estimation$ordered_n_2d_matrix_1224[t1 - 1, index - 1]) * 1 / estimation$n
  }
}
# calculation of the Integral for the next step
Integral_t <- Integral_t + d_I_t[1]
# free policy payments after 65:
for (t2 in ((index + 1):(length(estimation$ordered_times)))) {
  # we only add one summand to d_C_2 because the second one is always zero
  # because the integral over the payment functions is zero
  # the fact that we integrate until s- does not effect these calculations because
  # the payments are continuous
  # print(Integral_t)
  d_C_2[t2] <- d_C_2[t2 - 1] + Integral_t * pension *
    interest_f(estimation$ordered_times[t2 - 1], estimation$ordered_times[t2])

  for (t1 in (2:(index - 1))) {
    if (estimation$I_t_2d[1 + (t1 - 2) * 4, 1 + (t2 - 2) * 4]
    > 10^(-10)) {
      d_I_t[t2 - index + 1] <- d_I_t[t2 - index + 1] + rho_calc[t1] * estimation$p_11[t1 - 1, t2 - 1] *
        estimation$ordered_n_2d_matrix_1212[t1 - 1, t2 - 1] *
        1 / estimation$I_t_2d[1 + (t1 - 2) * 4, 1 + (t2 - 2) * 4] * 1 / estimation$n
    }
    if (estimation$I_t_2d[1 + (t1 - 2) * 4, 2 + (t2 - 2) * 4]
    > 10^(-10)) {
      d_I_t[t2 - index + 1] <- d_I_t[t2 - index + 1] + rho_calc[t1] * estimation$p_12[t1 - 1, t2 - 1] *
        1 / estimation$I_t_2d[1 + (t1 - 2) * 4, 2 + (t2 - 2) * 4] *
        (-estimation$ordered_n_2d_matrix_1223[t1 - 1, t2 - 1] -
          estimation$ordered_n_2d_matrix_1224[t1 - 1, t2 - 1]) * 1 / estimation$n
    }
  }
  Integral_t <- Integral_t + d_I_t[t2 - index + 1]
}
reserve_comp_3 <- reserve_comp + d_C_2



# transition payments:
# one-dimensional jumps:
# recursive calculation of the transition payments
# the deleted censoring in the code from Christian and Martin could be problematic
# vectorization could speed this up
jump_t <- numeric(length(estimation$ordered_times))
V_tech <- sapply(estimation$ordered_times, V_1_plus) - sapply(estimation$ordered_times, V_1_minus)
jump_t[2] <- V_tech[2] * estimation$increments[[1]][1, 3] * 1/(1 + interest / 100)^(estimation$ordered_times[2] - 40)
for (t in 3:(length(estimation$ordered_times))) {
  jump_t[t] <- jump_t[t - 1]
  if (estimation$It[[t - 2]][1] > 10^(-10)) {
    jump_t[t] <- jump_t[t] + estimation$p[[t - 1]][1] * V_tech[t] * 1 / estimation$It[[t - 2]][1] *
      estimation$increments[[t - 1]][1, 3] * ifelse(estimation$ordered_times[t] <= 65, 1, 1) *
      1 / (1 + interest / 100)^(estimation$ordered_times[t] - 40)
  }
}

reserve_comp_1 <- reserve_comp_3 + jump_t
# reserve_comp_1 misses only the free policy jumps from 2 to 3

# two-dimensional jumps from 1 to 2 and 2 to 3:
# assume that the transition payments are V_1(u_1)*(1+r)^(u_2-u_1)
# construction of the matrices
I_t_2d_12 <- estimation$I_t_2d[1 + ((2:(length(estimation$ordered_times) + 1)) - 2) * 4, 2 + ((2:(length(estimation$ordered_times) + 1)) - 2) * 4]

delta_Lambda_1 <- matrix(0, nrow = length(estimation$ordered_times) - 1, ncol = length(estimation$ordered_times) - 1)
for (t1 in (1:(length(estimation$ordered_times) - 1))) {
  for (t2 in (1:(length(estimation$ordered_times) - 1))) {
    if (I_t_2d_12[t1, t2] > 10^(-10)) {
      delta_Lambda_1[t1, t2] <- estimation$ordered_n_2d_matrix_1223[t1, t2] *
        estimation$p_12[t1, t2] * 1 / I_t_2d_12[t1, t2] * 1 / estimation$n
    }
  }
}
# integral_2d is a matrix that calculates all the pairs of times where there is a
# weighted jump. We can simplify this because rho does not change in u_2
# and the payments do not change in u_1
integral_2d <- matrix(0, nrow = length(estimation$ordered_times) - 1, ncol = length(estimation$ordered_times) - 1)
V_plus_calc <- sapply(estimation$ordered_times, V_1_plus)
for (t1 in (1:(length(estimation$ordered_times) - 1))) {
  integral_2d[t1, ] <- delta_Lambda_1[t1, ] * rho_calc[t1 + 1]
}
for (t1 in (1:(length(estimation$ordered_times) - 1))) {
  integral_2d[, t1] <- integral_2d[, t1] * V_plus_calc[t1 + 1] *
    1/(1 + interest / 100)^(estimation$ordered_times[t1 + 1] - 40)
}
# we only need to add in one dimension, because adding over [40,65]times[40,t] is eaqul to [40,t] times[ 40,t]
# because we can only jump from 2 to 3 after we already jumped from 1 to 2.
jump_t_2 <- numeric(length(estimation$ordered_times))
for (t1 in (2:length(estimation$ordered_times))) {
  jump_t_2[t1] <- jump_t_2[t1 - 1] + sum(integral_2d[1:(index - 2), t1 - 1])
}


reserve_comp_2 <- reserve_comp_1 + jump_t_2
# reserve_comp_2 contains all the payments.


return(reserve_comp_2)
}
