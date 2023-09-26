

#' # Title: Transmission & Exposure distributions - Mark 2 
#' * Author: Beth Savagar (RVC)
#' * Date: 30 March 2021
#' 

#' 
#' #########################################################################################################
#' ### Inputs
#' #########################################################################################################
#'
#' * T_het <- "heterogeneous"/"homogeneous"
#' * E_het <- heterogeneous"/"homogeneous""
#' * N <- population size
#' * mR0 <- mean of negative binomial distribution (close to target R0)
#' * k_estimate <- dispersion parameter of negative binomial
#' * R0_target <- target R0 value for population
#'

#' #########################################################################################################
#' ### Function
#' #########################################################################################################
#'


dist <- function(
  T_het, # homo or hetero
  E_het, # hom or het
  N,
  mR0,
  k_estimate,
  R0_target,
  correlation
){
  
  
  if (T_het == "homogeneous" &
      E_het == "homogeneous") {
    T_dist <- rep(mR0, N)
    E_dist <- rep(1 / N, N)
    R0est <- round(R0calc(T_dist, E_dist), 1)
      while (R0est != R0_target) {
        T_dist <- rep(mR0, N)
        E_dist <- rep(1 / N, N)
        R0est <- round(R0calc(T_dist, E_dist), 1)
      }
    
  } else if (T_het == "heterogeneous" &
             E_het == "homogeneous") {
    T_dist <- rnbinom(N, size = k_estimate, mu = mR0)
    E_dist <- rep(1 / N, N)
    R0est <- round(R0calc(T_dist, E_dist), 1)
    
    while (R0est != R0_target) {
      T_dist <- rnbinom(N, size = k_estimate, mu = mR0)
      E_dist <- rep(1 / N, N)
      R0est <- round(R0calc(T_dist, E_dist), 1)
    }
    
  } else if (T_het == "homogeneous" &
             E_het == "heterogeneous") {
    T_dist <- rep(mR0, N)
    E_dist <- rnbinom(N, size = k_estimate, mu = mR0)
    E_dist <- E_dist / sum(E_dist) #normalised
    R0est <- round(R0calc(T_dist, E_dist), 1)
    
    while (R0est != R0_target) {
      T_dist <- rep(mR0, N)
      E_dist <- rnbinom(N, size = k_estimate, mu = mR0)
      E_dist <- E_dist / sum(E_dist) #normalised
      R0est <- round(R0calc(T_dist, E_dist), 1)
    }
    
  } else if (T_het == "heterogeneous" &
             E_het == "heterogeneous") {
    # define target R0 and target correlation
    
    E_dist <- c()
    if(correlation == "positive"){
      E_dist <- T_dist
    }else if(correlation == "negative"){
      E_dist[order(T_dist)] <- rev(sort(T_dist))
    }else if(correlation == "none"){
      E_dist <- sample(T_dist)
    }
    
    R0est <- round(R0calc(T_dist, E_dist/sum(E_dist)), 1)
    
  }
return(list(T_dist, E_dist))
}

#' #########################################################################################################
#' ### Example
#' #########################################################################################################

dist_transmission <- function(N,m,k){
  
  if(k == "hom"){
    T_dist <- rep(m, N)
    E_dist <- rep(1 / N, N)
    R0est <- round(R0calc(T_dist, E_dist), 1)
    # while (R0est != m) {
    #   T_dist <- rep(m, N)
    #   E_dist <- rep(1 / N, N)
    #   R0est <- round(R0calc(T_dist, E_dist), 1)
    # }
  }else{
    k <- as.numeric(k)
    T_dist <- rnbinom(N, size = k, mu = m)
    E_dist <- rep(1 / N, N)
    R0est <- round(R0calc(T_dist, E_dist), 1)
    # 
    # while (R0est != m) {
    #   T_dist <- rnbinom(N, size = k, mu = m)
    #   E_dist <- rep(1 / N, N)
    #   R0est <- round(R0calc(T_dist, E_dist), 1)
    # }
  }
  
  return(list(T_dist, E_dist))
}



  
