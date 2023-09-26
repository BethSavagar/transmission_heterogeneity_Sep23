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

# All lines commented out below to avoid running when function is called.
# 
# library(ggplot2)
# library(ggpubr)
# #
# N <- 1e4
# k <- Inf
# m <- 2
# 
# T_dist <- rnbinom(N, size = k, mu = m)
# E_dist <- rep(1/N, N)
# #vE[order(vT)] <- rev(sort(vT))
# V_prop <- 0.1
# V_exp <- 5
# heterogeneity <- "fixed"
# 
# #vID <- vaccination(T_dist, E_dist, V_prop, V_exp)
# 
# 
# vTI <- T_dist +1
# vEI <- E_dist +1
# 
# v <- vTI /sum(vTI) + vEI/sum(vEI)
# 
# pV <- v^V_exp # probability to be vaccinated
# pV <- pV/ sum(pV)
# vID <- sample(N, V_prop*N, replace = FALSE, prob = pV)
# 
# pVdf <- as.data.frame(cbind(T_dist, pV))
# ggplot(data = pVdf, aes(x = T_dist, y = pV))+
#   geom_point()+
#   theme_bw()+
#   labs(x = "transmission potential", y = "proabability of vaccination")
# 
# 
# 
# ## Loop with different exponents and heterogeneity
# N <- 1e4
# m <- 2
# ks <- c(0.1, Inf)
# exps <- c(-5,0,5)
# pars <- expand.grid(ks,exps)
# 
# pVdf <- c()
# 
# 
# for(i in 1:nrow(pars)){
#   k <- pars[i,1]
#   V_exp <- pars[i,2]
#   
#   T_dist <- rnbinom(N, size = k, mu = m)
#   E_dist <- rep(1/N, N)
#   
#   vTI <- T_dist +1
#   vEI <- E_dist +1
#   
#   v <- vTI /sum(vTI) + vEI/sum(vEI)
#   
#   pV <- v^V_exp # probability to be vaccinated
#   pV <- pV/ sum(pV)
#   
#   new <- cbind("k" = k,"V_exp"=V_exp,"T_dist"=T_dist,"risk_all" =v, "pV" =pV)
#   
#   pVdf <- as.data.frame(rbind(pVdf, new))
#   
# }
# 
# #colnames(pVdf) <- c("k", "V_exp", "T_dist", "risk_all", "pV")
# expr1 <- bquote("Convenience," ~ V[exp] == -5 )
# expr2 <- bquote("Random," ~ V[exp] == 0 )
# expr3 <- bquote("Targeted," ~ V[exp] == 5 )
# 
# exp_labs <- c(expr1,expr2,expr3)
# 
# Tdist_pV_plot <- ggplot(pVdf, aes(x = T_dist, y = pV, group = as.factor(V_exp), col = as.factor(V_exp)))+
#   geom_line(size = 0.75)+
#   facet_wrap(~k, scales = "free", labeller = "label_both")+
#   labs(title = "Probability of vaccination by transmission potential", x = "Transmission Potential", y = "Probability of Vaccination", col = "Vaccination Strategy \n (Vaccine Exponent)")+
#   scale_color_discrete(labels = exp_labs)+
#   theme_bw()
# 
# Tdist_pV_LOGplot <- ggplot(pVdf, aes(x = T_dist, y = log(pV), group = as.factor(V_exp), col = as.factor(V_exp)))+
#   geom_line(size = 0.75)+
#   facet_wrap(~k, scales = "free_x", labeller = "label_both")+
#   labs(title = "LOG-SCALE : Probability of vaccination by transmission potential",x = "Transmission Potential", y = "Probability of vaccination, log-scale", col = "Vaccination Strategy \n (Vaccine Exponent)")+
#   scale_color_discrete(labels = c(exp_labs))+
#   theme_bw()
# 
# Tdist_pV_stacked <- ggarrange(Tdist_pV_plot, Tdist_pV_LOGplot, ncol = 1)
# 
# #ggsave("Tdist_pV_stacked.png",Tdist_pV_stacked,width = 4,height = 3, unit = "in" )
# 
# ?ggsave
# 
# ## highlighting targeted units
# points(T_dist, pV)
# 
# plot(T_dist)
# points(vID,T_dist[vID], col = "red")
