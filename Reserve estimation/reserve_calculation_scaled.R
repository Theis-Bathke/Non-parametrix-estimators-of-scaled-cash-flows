reserve_scaled_calc<-function(est,b_0,premium,interest,interest_f,pension,V_1_plus,V_1_minus){
index <- which(est$ordered_times > 65)[1]
reserve_scaled <- numeric(length(est$ordered_times))
reserve_scaled[1] <- -b_0
V_plus_calc <- sapply(est$ordered_times, V_1_plus)
V_tech <- sapply(est$ordered_times, V_1_plus) - sapply(est$ordered_times, V_1_minus)
# reserve_scaled calculation in state 1
for (t in 2:(index - 1)) {
  reserve_scaled[t] <- reserve_scaled[t - 1] + est$aj_scaled[[t - 1]][1] *
    interest_f(est$ordered_times[t - 1], est$ordered_times[t]) * -premium
}
# transition from before pension to after pension ( no free policy option)
reserve_scaled[index] <- reserve_scaled[index - 1] + est$aj_scaled[[index - 1]][1] *
  interest_f(est$ordered_times[index - 1], 65) * -premium + est$aj_scaled[[index - 1]][1] *
  interest_f(65, est$ordered_times[index]) * pension
# reserves for pensioners with no free policy option
for (t in (index + 1):(length(est$ordered_times))) {
  reserve_scaled[t] <- reserve_scaled[t - 1] + est$aj_scaled[[t - 1]][1] *
    interest_f(est$ordered_times[t - 1], est$ordered_times[t]) * pension
}
#reserve_scaled calculation in state 2
reserve_scaled_1_sum <- numeric(length(est$ordered_times))
for (t in 2:length(est$ordered_times)) {
  reserve_scaled_1_sum[t] <- reserve_scaled_1_sum[t-1] + est$aj_scaled[[t - 1]][2] *
    interest_f(max(est$ordered_times[t - 1],65), est$ordered_times[t]) *
    ifelse(est$ordered_times[t] <= 65, 0, pension)
}
reserve_scaled_1<-reserve_scaled + reserve_scaled_1_sum
#reserve_scaled calculation for the jump 2-3
b_ij<-function(t){
  b<-matrix(0,nrow=6,ncol=6)
  b[1,3]<-V_tech[t]
  b[2,5]<- V_plus_calc[t]
  return(b)
}
reserve_scaled_2_sum <- numeric(length(est$ordered_times))
reserve_scaled_2_sum[2]<- reserve_scaled_2_sum[1] + sum(est$aj_scaled[[1]] * est$contribution_first_scaled*b_ij(2))*
  1/(1 + interest / 100)^(est$ordered_times[2] - 40)
for (t in 3:length(est$ordered_times)) {
  reserve_scaled_2_sum[t] <- reserve_scaled_2_sum[t-1] + sum(est$aj_scaled[[t-1]] * est$contributions_scaled[[t-2]]*b_ij(t))*
    1/(1 + interest / 100)^(est$ordered_times[t] - 40)
}
reserve_scaled_2<-reserve_scaled_1 + reserve_scaled_2_sum
return(reserve_scaled_2)
}