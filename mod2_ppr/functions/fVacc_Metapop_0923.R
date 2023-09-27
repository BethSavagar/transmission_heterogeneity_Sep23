#' # Title: Vaccination Function for conTime model
#' * Author: Beth Savagar (RVC)
#' * Date: 16 Jun 2021
#' 

#' 
#' #########################################################################################################
#' ### Inputs
#' #########################################################################################################
#'
#' * `vT` <- Transmission potential distribtion
#' * `vS` <- Exposure potential distribution
#' * `V_prop` <- proportion to be vaccinated
#' * `V_exp` <- vaccine exponent: determines extent of selectivity 
#'     *  # `V_exp` == 0 is random vaccination, v_exp >0 is positive targeting, v_exp < 0 is negative targeting (convenience)
#'

#' #########################################################################################################
#' ### Function
#' #########################################################################################################
#'
vaccination <- function(T_dist, 
                        E_dist, 
                        V_prop, 
                        V_exp){
  
#V_exp == 0 is random vaccination, v_exp >0 is positive targeting, v_exp < 0 is negative targeting (convenience)
    N <- length(T_dist)
 
    vTI <- T_dist +1
    vEI <- E_dist +1
    
    v <- vTI /sum(vTI) * vEI/sum(vEI)
    v <- v/sum(v)
    
    pV <- v^V_exp # probability to be vaccinated
    pV <- pV/ sum(pV) # normalised
    v_ID <- sample(N, V_prop*N, replace = FALSE, prob = pV) # vID = id of vaccinated units

  
  return(v_ID)
}
