# Code to calculate proportion of simulations which become endemic

Vendemicity <- function(Tmat, TimeStop, Tvac){ # input function is matrix of numbers infected across generations (columns) for multiple simulations (rows)
  
  denom <- sum(Tmat[,Tvac-1]!=0) / nrow(Tmat) 
  
  num <- sum(Tmat[,TimeStop]!=0) / nrow(Tmat) # sum the number of simulations which result in 0 infections in the final generation (extinction), divided by the total number of simulaitons (number of row of matrix)
  
  P_endemic <- num/denom
  
  return(P_endemic) # return the proportion of simulations which become endemic
}
