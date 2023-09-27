# Example script rendering Model II: PPR Metapopulation model
# Author: Beth Savagar
# Date: 27.09.23

# Vaccination over 4 timesteps

rm(list=ls())
mod_filepath <- "mod2_ppr"

##################
## SETUP ##
##################

# load libraries
library(doParallel)
library(foreach)
library(tidyverse)
library(reshape2)
library(here)

# load functions
source(paste0(mod_filepath, "/functions/vfrunPPR_0923.R")) # run model with transmission within and between units
source(paste0(mod_filepath, "/functions/fPPR_OnePop_0923.R")) # simulates transmission within a unit (sub-population)
source(paste0(mod_filepath, "/functions/vfPPR_Metapop_0923.R")) # simulate transmission between units (meta-pop), with vaccination.
source(paste0(mod_filepath, "/functions/fVacc_Metapop_0923.R")) # vaccination
source(paste0(mod_filepath, "/functions/pEndemicity_ppr.R")) # calc. endemicity


##################
## PARAMETERS ##
##################

# posterior parameters
Post <- read.csv(paste0(mod_filepath, "/data-parameters/Posterior.csv"))
# fixed parameters script: 
source(paste0(mod_filepath, "/data-parameters/fixed-pars_0923.R")) # vacciantion over 4 timesteps


# Heterogeneous Transmission:
Variable_Transmission <- c("yes") # "yes": negative binomial; "no": all populations have the same transmission potential

# Variable Parameters: 
# Heterogeneous transmission and vaccination scenarios: 
k_vals <- c(0.1, 0.3, 0.5, 1, Inf) # dispersion parameter
v_units <- seq(0,1,0.05) # proportion vaccinated during each campaign
v_exp <- c(-4,-0.5,0,0.5,4) # vaccine exponent, defines vaccination strategy

# Create parameter dataframe storing parameter combinations for simulation: 
params <- list(
  hetero = Variable_Transmission, # heterogeneous transmission
  k = k_vals, # dispersion parameter
  v_units_vals = v_units, # proportion vaccinated
  v_exp = v_exp # vaccine exponent (strategy)
) %>%
  cross_df() %>%
  mutate(k = as.numeric(k),
         v_units_vals = as.numeric(v_units_vals))


##################
## MODEL SETUP ##
##################

## Metapopulation timesteps (1 timesteps = 10 days: infectious period of individual animal)
Timesteps <- 900 

# Model Output: 
output <- "I_counts"

# Simulation
runs <- 1000 # number of iterations per scenario

# Set-up for running simulations in parallel: 
#
#
#

##################
## RUN MOD ##
##################

# Loop through each scenario 
for(i in 1:nrow(params)){
  
  k_estimate <- params[i,"k"] %>% pull() # dispersion parameter for negative binomial distribution
  v_units <- params[i,"v_units_vals"] %>% pull() # proportion vaccinated during campaigns
  v_strat <- params[i,"v_exp"] %>% pull() # vaccine exponent (defines strategy)
  print(i)
  
  Runs_tracker <-
    foreach (z = 1:runs, .combine = "rbind") %dopar% {
      runPPR( # output of runPPR is an incidence within the metapopulation
        
        Post, # posterior parameters:
        n, # number of units (sub-pops) in metapopulation
        NInfSeededArea, # number of infected animals in seed unit (subpopulation)
        MinNoI, # minimum number infected within unit
        MinNoR, # minimum number recovered within unit
        TimeSteps_OnePop, # number of timeteps for within population transmission
        
        # within-unit demographic params
        Birth, # birth rate within subpop(equal throughout metapop)
        Exit_1_dt, # survival rate in group 1
        Exit_2_dt, # survival rate in gorup 2
        Age_r, # maturity rate 
        Mort_PPR, # ppr mortality
        PropAge_1, # proportion of animals in age group 1 (within sub pop)
        PropAge_2, # proportion of animals in age group 2 (within sub pop
        
        Npop, # number of animals in sub-pop (unit)
        Beta_w, # transmission rate within unit (sub-pop)
        nSeeds, # number of seed units (initial)
        Timesteps, # timesteps for metapopulation tranmission (where 1 TS = animal infection period, 10 days)
        Freq_ReIntr, # frequency of reseeding infection in metapopulation
        TimeLimit_ReIntr, # upper timelimit for reseeding events
        nReIntr, # number of units infected at reseeding event
        k_estimate, # k value 
        output, # formt of output:
        
        v_units, # vaccination coverage of units in metapopulation
        v_strat, # c(-4,0,4) targeting of units in metapopualtion
        vaccine_params # parameters for implementing vaccination within units
        )
    }
  
  # Proportion Endemicity: 
  # - calculate the proportion of simulations that are endemic at time of vaccination (V_start[1]-2), but extinct at end of simulation. 
  ppr_endemicity <-  Vendemicity(Runs_tracker, Timesteps, v_start[1]-2)
  
  # Store output in params dataframe
  params[i,"END"] <- ppr_endemicity
  
  # re-save after each iteration?
  # write.csv(params, paste0(mod_filepath, "/output/output_mod2_",Sys.Date(), ".csv"))
}

# Save completed output
tdate <- Sys.Date()
filename <- paste0("output_mod2_", tdate, ".csv")

write.csv(params, paste0(mod_filepath, "/output/", filename))

