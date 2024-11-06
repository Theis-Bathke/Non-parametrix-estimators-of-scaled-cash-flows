source("1d_univ.R")
source("monte_carlo_p.R")
source("one_dimensional_estimation.R")
source("semi_Markov_sim.R")

# data cannot be loaded, since the estimation of the scaled transition rates
# depend on the seed. 
set.seed(4)
n<-2000
n_mc<-10000
c<-runif(n,60,120)

#simulation
sim<-sim_example(n,c)
sim_m<-sim_example(n_mc,rep(120,n_mc))
#Setting:
#initial payment
b_0 <- 100000
# b_0 <- 0
#Premium payment
premium <- 10000
#pension (normally calculated via the technical basis)
p_11_technical <- function(x) {
  exp(-(10^((19 * x) / (500) - 159 / (125))) / (38 * log(10))
      - x / (2000)
      + 1 / (38 * 10^(159 / (125)) * log(10)))
}
pension <- (b_0 + premium * integrate(p_11_technical, 40, 65)$value) /
  (integrate(p_11_technical, 65, Inf)$value)
# pension<-10
## rho (normally calculated via the technical basis)
V_1_plus <- function(t) {
  if (t > 65) {
    result <- integrate(p_11_technical, t, Inf)$value *
      pension
  } else {
    result <- integrate(p_11_technical, 65, Inf)$value *
      pension
  }
  return(result)
}
# V_1_plus <- function(t){return(0)}
V_1_minus <- function(t) {
  if (t > 65) {
    result <- 0
  } else if (t > 40) {
    result <- integrate(p_11_technical, t, 65)$value * premium
  } else {
    result <- b_0 + integrate(p_11_technical, 40, 65)$value * premium
  }
  return(result)
}
rho <- function(t) {
  if (t >= 65) {
    return(1)
  } else if (t < 40) {
    return(0)
  } else if (V_1_plus(t)>10^(-12)){
    return((V_1_plus(t) - V_1_minus(t)) / V_1_plus(t))
  } else{return(0)}
}



est_cons<-aalen_johansen_1d(sim,rho)
p_2_scaled<-lapply(est_cons$aj_scaled,function(x)return(x[2]))
p_2_scaled_censored<-lapply(est_cons$It,function(x)return(x[2]))
est_var1<-aalen_johansen_var(sim,rho)
est_var2<-aalen_johansen_var(sim,rho)
est_var3<-aalen_johansen_var(sim,rho)

est_monte_carlo<-monte_carlo_p(sim_m,rho)

p_2_var1<-lapply(est_var1$aj,function(x)return(x[2]))
p_2_var2<-lapply(est_var2$aj,function(x)return(x[2]))
p_2_var3<-lapply(est_var3$aj,function(x)return(x[2]))

p_2_monte<-lapply(est_monte_carlo$It,function(x)return(x[2]))


index_mc_time<-which(est_monte_carlo$ordered_times[-1]>100)[1]
index_var_1_time<-which(est_var1$ordered_times[-1]>100)[1]
index_var_2_time<-which(est_var2$ordered_times[-1]>100)[1]
index_var_3_time<-which(est_var3$ordered_times[-1]>100)[1]
index_cons_time<-which(est_cons$ordered_times[-1]>100)[1]

plot(est_cons$ordered_times[1:index_cons_time-1],p_2_scaled[1:index_cons_time-1],col = "red",lty = 2,type = "l",xlab="Age",ylab="",cex.axis=1.2)
lines(est_var1$ordered_times[1:index_var_1_time-1],p_2_var1[1:index_var_1_time-1],col= "green",lty = 3)
lines(est_var2$ordered_times[1:index_var_2_time-1],p_2_var2[1:index_var_2_time-1],col= "green",lty = 3)
lines(est_var3$ordered_times[1:index_var_3_time-1],p_2_var3[1:index_var_3_time-1],col= "green",lty = 3)
lines(est_monte_carlo$ordered_times[-1][1:index_mc_time-1],p_2_monte[1:index_mc_time-1],col= "black")
# lines(est_cons$ordered_times[-1],p_2_scaled_censored,col = "blue")
legend("topright",                
       legend = c("CMAJ","SAJ", "True"),  
       col = c("green", "red","black"),   
       lty = c(3, 2,1),             
       lwd = 1,
       cex = 1) 

#calculation of the scaled transition probabilities with the change of measure
#technique for the histogram
p_2_var<-rep(0,100)
for (i in 1:100){
  current<-aalen_johansen_var(sim,rho)
  index=which(current$ordered_times>=60)[1]-1
  p_2_var[i]<-current$aj[[index]][2]
}
hist(p_2_var,breaks=15,main = ""  , xlab = "",cex.axis=1.2,cex.lab=1.2,border="green",freq = TRUE)

index_mc=which(est_monte_carlo$ordered_times>=60)[1]-1
p_2_monte_carlo_60=p_2_monte[index_mc]
index=which(est_cons$ordered_times>=60)[1]-1
p_2_cons_60=p_2_scaled[index]

abline(v = c(p_2_monte_carlo_60,p_2_cons_60,mean(p_2_var)), col = c("black","red","green"), lwd = 2, lty = c(1,2,3))
legend(
  "topright",                    
  legend = c("CMAJ", "SAJ","True"), 
  lty = c(3, 2,1),    
  col = c("green", "red","black"),
  lwd = 1,
  y.intersp = 1, # Increase vertical spacing
  x.intersp = 1, # Increase horizontal spacing
  cex = 1
)

