# Example script rendering Model I
# Author: Beth Savagar
# Date: 26.09.23

mod_filepath <- "mod1_IBM"


##################
## SETUP ##
##################

# load libraries
source(paste0(mod_filepath,"/scripts/setup.R"))

# load functions
source(paste0(mod_filepath, "/functions/continuousIBM.R")) # continuousIBM function
source(paste0(mod_filepath, "/functions/vaccination.R")) # vaccination function
source(paste0(mod_filepath, "/functions/TEdists.R")) # dist and dist_transmission function
source(paste0(mod_filepath, "/functions/R0calc.R")) # R0calc
source(paste0(mod_filepath, "/functions/pEndemicity.R")) # Vendemicity


##################
## PARAMETERS ##
##################

source(paste0(mod_filepath,"/data-parameters/fixed-pars.R")) 
#on local machine see CLEAN/sims_update_Nov22 & sims_update_Jul23

##################
## MODEL SETUP ##
##################

# Model Output: 
output <- "I_counts" # "Counts" for SIRV counts, I_counts for only Infecteds 

# epiList:
# List of matrices, 1 matrix per scenario in Params dataframe
# Matrices contain output of each simulation ("runs")
epi_List <- vector(mode = "list", length = nrow(params)) # List to store output fo

# Simulation: 
runs <- 5e3 # up to 1000 once running effectively

# Loop through each scenario: 
for(i in 1:nrow(params)){
  
  m <- params[i, "m"] %>% pull() # m = R0
  k <- params[i, "k"] %>% pull() # k = dispersion parameter for NBD
  I_period <- params[i, "I_period"] %>% pull() # Infectious period
  R_period <- params[i, "R_period"] %>% pull() # Immune period
  IR_period <- I_period+R_period # total infectious and immune period
  V_prop <- params[i, "prop_vaccinated"] %>% pull() # proportion vaccinated each round
  V_exp <- params[i, "exponent"] %>% pull() # vaccine strategy (see explanation in `fixed_pars.R`)
 
  print(i)
  
  
  # Run simulations (=runs) for given scenario (i) and store output in Runs_tracker matrix
  Runs_tracker <-
    foreach (z = 1:runs, .combine = "rbind") %dopar% {
      
      continuousIBM(
        TimeStop, # duration of simulation
        DigitsRound, # digits for timstep
        ts, # timestep
        N, # population size
        m, # = R0, mean transmission potential 
        k, # dispersion parameter
        states, # model states
     
        I_period, # infectious period
        R_period, # immune period
        IR_period, # infectious + immune period
        V_period, # vaccination duration
        I_time, # time of infective contact 
        IS_time, # time of return to susceptibility for infected units, set to 0 for units which are never infected
        VS_time, # time of return to susceptibility of vaccinated units, set to 0 for units which are never vaccinated
        V_schedule, # vaccination schedule (start time, frequency, number of rounds)
        V_prop, # proportion vaccinated at each round
        V_exp, # vaccine exponent, defines vaccination strategy
        output # model output
      )
    }

  # Proportion Endmeicity: 
    # - calculate the proportion of simulations that are endemic at time of vaccination (V_schedule-2), but extinct at end of simulation. 
  pEndemic <- Vendemicity(Runs_tracker, TimeStop, V_schedule - 2) 
  
  # Store output in params dataframe
  params[i,"pEnd"] <- pEndemic
  
  # re-save after each iteration?
  # write.csv(params, paste0(mod_filepath, "/output/output_mod1_",Sys.Date(), ".csv"))
}



tdate <- Sys.Date()
filename <- paste0("output_mod1_", tdate, ".csv")

write.csv(params, paste0(mod_filepath, "/outputs/", filename))
