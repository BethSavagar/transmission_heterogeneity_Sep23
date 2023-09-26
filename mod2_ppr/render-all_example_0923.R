# Example code for Transmission Heterogeneity Model I - Theoretical Agent-Based Model #

filepath <- "/storage/users/bsavagar/transmission_heterogeneity/"

source(paste0(filepath,"scripts/setup.R"))
source(paste0(filepath,"functions/V_Continuous_IBM_function_0823.R"))
source(paste0(filepath,"functions/V_function_ct.R"))
source(paste0(filepath,"functions/TEdists_CT_070621.R"))
source(paste0(filepath,"functions/pEndemicity_function_update.R"))
source(paste0(filepath,"functions/R0calc.R"))
#
###############################################################################################################
# Set Up for Parallel Computing
# 
cores=detectCores()
cl <- makeCluster(50) #not to overload your computer
# cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
########
# File saving:
#
tdate <- Sys.Date()
filename <- paste0("vaccination_output_", tdate, ".RData")

source(paste0(filepath,"scripts/params_vaccination_update0723.R")) 
#on local machine see CLEAN/sims_update_Nov22 & sims_update_Jul23

N <- 1e4 # population
TimeStop <- 7*365
output <- "I_counts" # changed to counts output from I_counts 2/7/21 # back to i_counts 120721 to look at longer term pattern

## Number of runs
runs <- 5e3 # up to 1000 once running effectively
 
params <- params %>% 
  mutate(prop_vaccinated = as.character(prop_vaccinated)) %>%
  filter(exponent == 0) %>%
  # for jul23 rerun simulation filter for values where extinction is not 0 or 1
  filter(
    (k==0.1 & prop_vaccinated %in% seq(0.04,0.5,0.02)) |
     (k==0.3 & prop_vaccinated %in% seq(0.24,0.56,0.02)) |
      (k==1 & prop_vaccinated %in% seq(0.34,0.6,0.02)) |
      (k==Inf & prop_vaccinated %in% seq(0.34,0.6,0.02)) |
      (k=="hom" & prop_vaccinated %in% seq(0.4,0.6,0.02))
    ) %>%
  mutate(prop_vaccinated = as.numeric(prop_vaccinated)) 


# k=0.1	0.04-0.5
# k=0.3	0.24-0.56
# k=1	0.34-0.6
# K=Inf	0.34-0.6
# K=Hom?	0.4-0.6


epi_List <- vector(mode = "list", length = nrow(params))

system.time(
for(i in 1:nrow(params)){
  
  m <- params[i, "m"] %>% pull()
  k <- params[i, "k"] %>% pull()
  I_period <- params[i, "I_period"] %>% pull()
  R_period <- params[i, "R_period"] %>% pull()
  IR_period <- I_period+R_period
  V_prop <- params[i, "prop_vaccinated"] %>% pull()
  V_exp <- params[i, "exponent"] %>% pull()
  # distributions  
  # code below to be incorporated in IBM function
  # dists <- dist_transmission(N,m,k)
  # T_dist <- dists[[1]]
  # E_dist <- dists[[2]]
  print(i)
  
  Runs_tracker <-
    foreach (z = 1:runs, .combine = "rbind") %dopar% {
      
      V_start <- round(365*3,digits=0)
      
      continuousIBM(
        TimeStop,
        DigitsRound,
        ts,
        N,
        m,
        k,
        states,
        # seeds,# n, m, k, added above to generate distributions independently in each run.
        # T_dist,
        # E_dist,
        I_period,
        R_period,
        IR_period,
        V_period,
        I_time,
        IS_time,
        # time of return to susceptibility for infected units, set to 0 for units which are never infected
        VS_time,
        # time of return to susceptibility of vaccinated units, set to 0 for units which are never vaccinated
        V_schedule,
        V_prop,
        V_exp,
        output
      )
    }
  
  
  epi_List[[i]] <- Vendemicity(Runs_tracker, TimeStop, 1093) # calc endemicity at 2.5 years & end of sim
  params[i,"pEnd"] <- epi_List[[i]]
  write.csv(params, paste0(filepath, "outputs/", "rand_overwrite.csv"))
}

)

tdate <- Sys.Date()
filename <- paste0("vaccination_output_Rand_", tdate, ".csv")

write.csv(params, paste0(filepath, "outputs/", filename))
