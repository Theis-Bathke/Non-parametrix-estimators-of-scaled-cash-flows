aalen_johansen_2d<-function(sim)
{


n <- length(sim)
# number of states minus 1
p <- max(unique(unlist(lapply(sim, function(Z) Z$states)))) - 1
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
# gives a Vector that shows the change in the counting processes
ordered_n <- rep(list(matrix(0, p + 1, p + 1)), nrow(ordered_jumps))
for (i in 1:nrow(ordered_jumps)) {
  ordered_n[[i]][ordered_jumps[i, ][1], ordered_jumps[i, ][2]] <- 1
}
ordered_n <- lapply(ordered_n, FUN = function(Z) Z - diag(diag(Z)))
#catches the censoring times:
R_times <- unlist(lapply(sim, FUN = function(Z) tail(Z$times, 1)))
#safes the time-points that are a censoring time or "finish a path"
decisions <- ordered_times %in% R_times
# the difference of the indicator-function, this can be used to
# calculate the amount of people in the 4 states, by
# adding the corresponding vectors together
colsums_of_n <- lapply(ordered_n, function(N) {
  colSums(N - t(N))
})
out <- out2 <- list()
# accumulated scaled jumps from initial states to final states
out[[1]] <- ordered_n[[1]] * n^(-1)
for (tm in 2:(length(ordered_times) - 1)) {
  out[[tm]] <- out[[tm - 1]] + ordered_n[[tm]] * n^(-1)
}

# vector of differences of the estimators for every jump-time.
# Censoring is included
out2[[1]] <- colsums_of_n[[1]] * n^(-1)
if (decisions[2]) {
  wch <- ordered_individuals[2]
  end_state <- as.numeric(1:(p + 1) == tail(sim[[wch]]$states, 1))
  out2[[1]] <- out2[[1]] - end_state * n^(-1)
}
for (tm in 2:(length(ordered_times) - 1)) {
  out2[[tm]] <- out2[[tm - 1]] + colsums_of_n[[tm]] * n^(-1)
    if (decisions[tm+1]) {
       wch <- ordered_individuals[tm+1]
      end_state <- as.numeric(1:(p + 1) == tail(sim[[wch]]$states, 1))
      out2[[tm]] <- out2[[tm]] - end_state * n^(-1)
   }
}
# Vector of Starting values of the indicator vector for every insured
# is only really needed if not all insured start in the correct state
I_initial <- lapply(sim, FUN = function(Z) {
  as.numeric(1:(p + 1) == head(Z$states, 1))
})
# Starting value for the estimator (connected to the previous definition)
I0 <- Reduce("+", I_initial) * n^(-1)
# Constructing the estimator for the censored state occupation probabilities
It <- lapply(out2, FUN = function(N) {
  return(I0 + N)
})
# increments of the jump-processes
increments <- list()
increments[[1]] <- out[[1]]
for (i in 2:length(out)) {
  increments[[i]] <- out[[i]] - out[[i - 1]]
}
# construction of the increments for Lambda
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
# accumulating these singular summands of Lambda leads to the whole development
# of Lambda
cumsums <- list()
cumsums[[1]] <- contribution_first
for (i in 2:length(out)) {
  cumsums[[i]] <- contributions[[i - 1]] + cumsums[[i - 1]]
}
# solving the integral-equation by recursively
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

# adding the diagonal to Lambda
cumsums <- append(list(matrix(0, p + 1, p + 1)), lapply(cumsums,
  FUN = function(M) {
    M_out <- M
    diag(M_out) <- -rowSums(M)
    M_out
  }
))

# two-dimensional estimation

ordered_n_2d_sparse <- slam::simple_sparse_array(
  matrix(c(1, 1, 1, 1, 1, 1),
    ncol = 6, byrow = TRUE
  ), 0,
  dim = c(
    nrow(ordered_jumps),
    nrow(ordered_jumps), p + 1, p + 1, p + 1, p + 1
  )
)
ordered_n_2d_sparse[1, 1, 1, 1, 1, 1] <- 0
for (tm1 in 1:(length(ordered_times) - 1)) {
  for (tm2 in 1:(length(ordered_times) - 1)) {
    if (ordered_individuals[tm1 + 1] == ordered_individuals[tm2 + 1]) {
      temp <- c(
        tm1, tm2,
        which(ordered_n[[tm1]] == 1, arr.ind = TRUE),
        which(ordered_n[[tm2]] == 1, arr.ind = TRUE)
      )
      if (length(temp)==6){ #this if clause removes jumps like 2->2 which happens because of censoring
      ordered_n_2d_sparse$i <- rbind(
        ordered_n_2d_sparse$i,
        temp)
      ordered_n_2d_sparse$v <- c(ordered_n_2d_sparse$v, 1)
      }
    }
  }
}


colsums_n_2d <- matrix(0,
  nrow = nrow(ordered_jumps) * (p + 1),
  ncol = nrow(ordered_jumps) * (p + 1)
)


x1 <- slam::rollup(ordered_n_2d_sparse[, , , , , ], c(3, 5),
  FUN = sum, DROP = TRUE
)
x2 <- slam::rollup(ordered_n_2d_sparse[, , , , , ], c(4, 6),
  FUN = sum, DROP = TRUE
)
x3 <- slam::rollup(ordered_n_2d_sparse[, , , , , ], c(3, 6),
  FUN = sum, DROP = TRUE
)
x4 <- slam::rollup(ordered_n_2d_sparse[, , , , , ], c(4, 5),
  FUN = sum, DROP = TRUE
)

indexes <- cbind((x1$i[, 1] - 1) * (p + 1) + x1$i[, 3], (x1$i[, 2] - 1) * (p + 1) + x1$i[, 4])
indexes <- rbind(indexes, cbind((x2$i[, 1] - 1) * (p + 1) + x2$i[, 3], (x2$i[, 2] - 1) * (p + 1) + x2$i[, 4]))
indexes <- rbind(indexes, cbind((x3$i[, 1] - 1) * (p + 1) + x3$i[, 3], (x3$i[, 2] - 1) * (p + 1) + x3$i[, 4]))
indexes <- rbind(indexes, cbind((x4$i[, 1] - 1) * (p + 1) + x4$i[, 3], (x4$i[, 2] - 1) * (p + 1) + x4$i[, 4]))


values <- c(
  x1$v, x2$v, -1 * x3$v, -1 * x4$v
)
colsums_n_2d[indexes] <- values




colsums_n_2d_censored <- matrix(0,
                       nrow = nrow(ordered_jumps) * (p + 1),
                       ncol = nrow(ordered_jumps) * (p + 1)
)
for (tm1 in 2:(length(ordered_times))) {
  for (tm2 in 2:(length(ordered_times))) {
    if (ordered_individuals[tm1]==ordered_individuals[tm2]){
      colsums_n_2d_censored[((tm1 - 2) * (p + 1) + 1):((tm1 - 1) * (p + 1)-2),
      ((tm2 - 2) * (p + 1) + 1):((tm2 - 1) * (p + 1)-2)]<-
        colsums_n_2d[((tm1 - 2) * (p + 1) + 1):((tm1 - 1) * (p + 1)-2),
                              ((tm2 - 2) * (p + 1) + 1):((tm2 - 1) * (p + 1)-2)]
    }
  }
}
I_t_2d <- matrix(0,
  nrow = length(ordered_times) * (p + 1),
  ncol = length(ordered_times) * (p + 1)
)


I_t_2d[1:(p + 1), 1:(p + 1)] <- outer(I0, I0, "*")

for (tm in 2:(length(ordered_times))) {
  temp1 <- which(It[[tm - 1]] > 10^(-14))
  I_t_2d[(tm-1)*(p+1)+ temp1,1]<- It[[tm-1]][temp1]
  I_t_2d[1,(tm-1)*(p+1)+ temp1]<- It[[tm-1]][temp1]
}
#gives the decisions that are actually censored and not only jumps into 3 or 4
decisions_part<-logical(length(decisions))
for (i in 2:length(decisions)){
  if(decisions[i]==TRUE){
    wch<-ordered_individuals[i]
    if (tail(sim[[wch]]$states,1)<=2){
      decisions_part[i]<-decisions[i]
    }
  }
}


for (tm1 in 2:(length(ordered_times))) {
  for (tm2 in 2:(length(ordered_times))) {
    I_t_2d[
      ((tm1 - 1) * (p + 1) + 1):(tm1 * (p + 1)),
      ((tm2 - 1) * (p + 1) + 1):(tm2 * (p + 1))
    ] <- -I_t_2d[
      ((tm1 - 2) * (p + 1) + 1):((tm1 - 1) * (p + 1)),
      ((tm2 - 2) * (p + 1) + 1):((tm2 - 1) * (p + 1))
    ] + I_t_2d[
      ((tm1 - 1) * (p + 1) + 1):((tm1) * (p + 1)),
      ((tm2 - 2) * (p + 1) + 1):((tm2 - 1) * (p + 1))
    ] + I_t_2d[
      ((tm1 - 2) * (p + 1) + 1):((tm1 - 1) * (p + 1)),
      ((tm2 - 1) * (p + 1) + 1):((tm2) * (p + 1))
    ] + colsums_n_2d_censored[
      ((tm1 - 2) * (p + 1) + 1):((tm1 - 1) * (p + 1)),
      ((tm2 - 2) * (p + 1) + 1):((tm2 - 1) * (p + 1))
    ] * 1 / n
    if ((decisions_part[tm1 ] | decisions_part[tm2 ]) &
        ordered_individuals[tm1 ] == ordered_individuals[tm2 ]) {
      wch <- ordered_individuals[tm1]
      if (tail(sim[[wch]]$states,1)==2){
        if (tm1>tm2){
          I_t_2d[(tm1 - 1) * (p + 1) + 2,(tm2 - 1) * (p + 1) + 1]<-I_t_2d[(tm1 - 1) * (p + 1) +2,(tm2 - 1) * (p + 1) + 1]+1/n
          I_t_2d[(tm1 - 1) * (p + 1) + 2,(tm2 - 1) * (p + 1) + 2]<-I_t_2d[(tm1 - 1) * (p + 1) +2,(tm2 - 1) * (p + 1) + 2]-1/n
        }else if (tm1<tm2){
          I_t_2d[(tm1 - 1) * (p + 1) + 1,(tm2 - 1) * (p + 1) + 2]<-I_t_2d[(tm1 - 1) * (p + 1) +1,(tm2 - 1) * (p + 1) + 2]+1/n
          I_t_2d[(tm1 - 1) * (p + 1) + 2,(tm2 - 1) * (p + 1) + 2]<-I_t_2d[(tm1 - 1) * (p + 1) +2,(tm2 - 2) * (p + 1) + 2]-1/n
        } else{
          I_t_2d[(tm1 - 1) * (p + 1) + 2,(tm2 - 1) * (p + 1) + 2]<-I_t_2d[(tm1 - 1) * (p + 1) +2,(tm2 - 1) * (p + 1) + 2]+1/n
        }
        
      }
      if (tail(sim[[wch]]$states,1)==1){
        I_t_2d[(tm1 - 1) * (p + 1) + 1,(tm2 - 1) * (p + 1) + 1]<-I_t_2d[(tm1 - 1) * (p + 1) +1,(tm2 - 1) * (p + 1) + 1]+1/n
      }
    }
  }
}
from_16d_array_to_matrix <- function(x, small_rep) {
  # small_rep is normally equal to (p+1)
  index_i <- seq(0, 0, length.out = nrow(x$i))
  index_j <- seq(0, 0, length.out = nrow(x$i))
  index_v <- seq(0, 0, length.out = nrow(x$i))
  for (k in 1:nrow(x$i)) {
    index_i[k] <- (x$i[k, 1] - 1) *
      small_rep^2 + (x$i[k, 3] - 1) * small_rep + x$i[k, 4]
    index_j[k] <- (x$i[k, 2] - 1) *
      small_rep^2 + (x$i[k, 5] - 1) * small_rep + x$i[k, 6]
    index_v[k] <- x$v[k]
  }
  result <- Matrix::sparseMatrix(index_i, index_j, x = index_v, repr = "T")
  result@Dim <- as.integer(c(x$dim[1] * x$dim[3]^2, x$dim[2] * x$dim[4]^2))
  return(result)
}
ordered_n_2d <- from_16d_array_to_matrix(ordered_n_2d_sparse, 4)





#function to calculate Lambda 
calcualation_lambda <- function(It, N_t, index, times, states) {
  Lambda <- array(0, dim = c(times, times))
  Nt_matrix <- as.matrix(N_t[
    ((index[1] - 1) * states + index[2]) + (2:(times) - 2) * states^2,
    ((index[3] - 1) * states + index[4]) + (2:(times) - 2) * states^2
  ])
  for (tm1 in 2:(times)) {
    #print(tm1 / times)
    for (tm2 in 2:(times)) {
      if (It[index[1] + (tm1 - 2) * states, index[3] + (tm2 - 2) * states]
          > 10^(-10)) {
        Lambda[tm1, tm2] <-
          1 / n *
          Nt_matrix[tm1 - 1, tm2 - 1] *
          1 / It[
            index[1] + (tm1 - 2) * states,
            index[3] + (tm2 - 2) * states
          ] +
          Lambda[tm1 - 1, tm2] +
          Lambda[tm1, tm2 - 1] -
          Lambda[tm1 - 1, tm2 - 1]
      } else {
        Lambda[tm1, tm2] <-
          Lambda[tm1 - 1, tm2] +
          Lambda[tm1, tm2 - 1] -
          Lambda[tm1 - 1, tm2 - 1]
      }
    }
  }
  return(Lambda)
}



# construction of the two-dimensional transition probabilities.
p_11 <- array(0, dim = c(
  nrow(ordered_jumps) + 1, nrow(ordered_jumps) + 1
))
p_12 <- array(0, dim = c(
  nrow(ordered_jumps) + 1, nrow(ordered_jumps) + 1
))
p_13 <- array(0, dim = c(
  nrow(ordered_jumps) + 1, nrow(ordered_jumps) + 1
))
p_11[1, ] <- sapply(aj, function(x) if(x[1]>10^(-10)){return(x[1])}else{return(0)})
p_11[, 1] <- sapply(aj, function(x) if(x[1]>10^(-10)){return(x[1])}else{return(0)})
p_12[1, ] <- sapply(aj, function(x) if(x[2]>10^(-10)){return(x[2])}else{return(0)})
p_13[1, ] <- sapply(aj, function(x) if(x[3]>10^(-10)){return(x[3])}else{return(0)})
ordered_n_2d_matrix_1212 <- as.matrix(ordered_n_2d[
  ((1 - 1) * (p + 1) + 2) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2,
  ((1 - 1) * (p + 1) + 2) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2
])
ordered_n_2d_matrix_1313 <- as.matrix(ordered_n_2d[
  ((1 - 1) * (p + 1) + 3) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2,
  ((1 - 1) * (p + 1) + 3) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2
])
ordered_n_2d_matrix_1414 <- as.matrix(ordered_n_2d[
  ((1 - 1) * (p + 1) + 4) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2,
  ((1 - 1) * (p + 1) + 4) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2
])
for (tm1 in 2:(length(ordered_times))) {
  #print(tm1 / length(ordered_times))
  for (tm2 in 2:(length(ordered_times))) {
    if (I_t_2d[1 + (tm1 - 2) * 4, 1 + (tm2 - 2) * 4]
        > 10^(-14)) {
      p_11[tm1, tm2] <- (-1) * p_11[tm1 - 1, tm2 - 1] +
        p_11[tm1, tm2 - 1] +
        p_11[tm1 - 1, tm2] +
        p_11[tm1 - 1, tm2 - 1] * 1 / n *
        (ordered_n_2d_matrix_1212[tm1 - 1, tm2 - 1] +
           ordered_n_2d_matrix_1414[tm1 - 1, tm2 - 1] +
           ordered_n_2d_matrix_1313[tm1 - 1, tm2 - 1]) *
        1 / I_t_2d[1 + (tm1 - 2) * 4, 1 + (tm2 - 2) * 4]
    } else {
      p_11[tm1, tm2] <- (-1) * p_11[tm1 - 1, tm2 - 1] +
        p_11[tm1, tm2 - 1] +
        p_11[tm1 - 1, tm2]
    }
    if (p_11[tm1, tm2]<10^(-15)) {
      p_11[tm1, tm2] <- 0
    }
  }
}


# I_t_2d[1 + ((2:11) - 2) * 4, 2 + ((2:11) - 2) * 4]
ordered_n_2d_matrix_1223 <- as.matrix(ordered_n_2d[
  ((1 - 1) * (p + 1) + 2) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2,
  ((2 - 1) * (p + 1) + 3) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2
])
ordered_n_2d_matrix_1224 <- as.matrix(ordered_n_2d[
  ((1 - 1) * (p + 1) + 2) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2,
  ((2 - 1) * (p + 1) + 4) + (2:(nrow(ordered_jumps) + 1) - 2) * (p + 1)^2
])
for (tm2 in 2:(length(ordered_times))) {
  #print(tm1 / length(ordered_times))
  for (tm1 in 2:(length(ordered_times))) {
    p_12[tm1, tm2] <- (-1) * p_12[tm1 - 1, tm2 - 1] +
      p_12[tm1, tm2 - 1] +
      p_12[tm1 - 1, tm2]
    if (I_t_2d[1 + (tm1 - 2) * 4, 1 + (tm2 - 2) * 4]
        > 10^(-14)) {
      p_12[tm1, tm2] <- p_12[tm1, tm2] + p_11[tm1 - 1, tm2 - 1] * 1 / n *
        (-1) * ordered_n_2d_matrix_1212[tm1 - 1, tm2 - 1] *
        1 / I_t_2d[1 + (tm1 - 2) * 4, 1 + (tm2 - 2) * 4]
    } 
    if (I_t_2d[1 + (tm1 - 2) * 4, 2 + (tm2 - 2) * 4]
        > 10^(-14)) {
      p_12[tm1, tm2] <- p_12[tm1, tm2] + p_12[tm1 - 1, tm2 - 1] * 1 / n *
        (ordered_n_2d_matrix_1223[tm1 - 1, tm2 - 1] +
           ordered_n_2d_matrix_1224[tm1 - 1, tm2 - 1]) *
        1 / I_t_2d[1 + (tm1 - 2) * 4, 2 + (tm2 - 2) * 4]
    }

  }

}




return(list(n=n,p_11=p_11,p_12=p_12,p=aj,
            ordered_times=ordered_times,increments=increments,
            It=It,I_t_2d=I_t_2d,
            ordered_n_2d_matrix_1212=ordered_n_2d_matrix_1212,
            ordered_n_2d_matrix_1223=ordered_n_2d_matrix_1223,
            ordered_n_2d_matrix_1224=ordered_n_2d_matrix_1224))
}

