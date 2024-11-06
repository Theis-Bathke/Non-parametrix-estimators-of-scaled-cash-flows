source("reserve_calculation_2d.R")
source("reserve_calculation_scaled.R")
source("two_dimensional_estimation.R")
source("Monte_carlo_calculation.R")
source("one_dimensional_estimation.R")
load("simulated_Data_reserve_2000.Rdata")
#load("simulated_Data_reserve_5000.Rdata")
load("Monte_Carlo_Data.Rdata")
#Setting:
#initial payment
b_0 <- 100000
#Premium payment
premium <- 10000
#interest function

interest_f <- function(t1, t2) {
  return(t2-t1)
}
interest = 0
#pension (normally calculated via the technical basis)
p_11_technical <- function(x) {
  exp(-(10^((19 * x) / (500) - 159 / (125))) / (38 * log(10))
      - x / (2000)
      + 1 / (38 * 10^(159 / (125)) * log(10)))
}
pension <- (b_0 + premium * integrate(p_11_technical, 40, 65)$value) /
  (integrate(p_11_technical, 65, Inf)$value)
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
#estimation 2-D
estimation<-aalen_johansen_2d(sim)
#estimation 1-D (dependent on rho):
est<-aalen_johansen_1d(sim,rho)
#reserve calculation
reserve_2<-reserve_2d(estimation,b_0,premium,
                      interest,interest_f,pension,rho,V_1_plus,V_1_minus)
monte_carl<-monte_carlo_2d(sim_mc)
reserve_sc<-reserve_scaled_calc(est,b_0,premium,
                                interest,interest_f,pension,V_1_plus,V_1_minus)
#graphical analysis
index_2d <- which(estimation$ordered_times > 100)[1]
index_mc <- which(monte_carl$ordered_times > estimation$ordered_times[index_2d-1])[1]
if (is.na(index_2d)){index_2d=length(estimation$ordered_times)+1}
if (is.na(index_mc)){index_mc=length(monte_carl$ordered_times)+1}
options(scipen = 999)  
plot(estimation$ordered_times[1:(index_2d-1)],reserve_2[1:(index_2d-1)],ylim=c(min(reserve_2[1:(index_2d-1)]),max(reserve_2[(index_2d-1)],monte_carl$reserve_cont[(index_mc-1)],reserve_sc[(index_2d-1)])) ,col = "blue",type = "l", lty = 5,xlab="Age",ylab="",cex.axis=1.2,cex.lab=1.2)
lines(monte_carl$ordered_times[1:(index_mc-1)],monte_carl$reserve_cont[1:(index_mc-1)],col= "black",lty=1)
lines(est$ordered_times[1:(index_2d-1)],reserve_sc[1:(index_2d-1)],col= "red", lty=2)
legend("bottomright",                
       legend = c( "2dAJ","SAJ","True"),  
       col = c( "blue","red","black"),   
       lty = c(5, 2,1),             
       lwd = 1,
       cex = 1,
       y.intersp = 1.3, # Increase vertical spacing
       x.intersp = 1 # Increase horizontal spacing
)

