
# Model I parameters (26.09.23)
# Explore impact of transmission heterogeneity on endemicity of directly transmitted pathogen, under vaccination. 

debug <- T # set to F if want to save parameter csv

#########################
## VARAIBLE PARAMETERS ##
#########################

# Values as per manuscript examples:

I_period_values <- c(20) # infectious period = 20days
ratio_values <- c(3) # alpha (infections - immune period ratio)
m_values <- 2 # m=R0
k_values <- c(0.1,0.3,1,Inf,"hom") # k = dispersion parameter for NBD
proportion_vaccinated <- seq(0,1,0.02) # proportion vaccinated

# vaccine exponent:
V_exp <- c(-4,-0.5, 0, 0.5, 4)  # Determines vacc strategy (V_exp<0 Convenience, V_exp = 0 Random, V_exp >0 Targeted)

# Create parameter dataframe sotring parameter combinations for simulation: 
params <- list(
  m = m_values,
  k = k_values,
  I_period = I_period_values,
  ratio = ratio_values,
  prop_vaccinated = proportion_vaccinated,
  exponent = V_exp
) %>% 
  cross_df() %>%
  mutate(R_period = I_period*ratio)

# Save parameter dataframe if needed:

if(debug==F){
  tdate<- Sys.Date()
  write.csv(params, paste0("output/params_",Sys.Date(), ".csv"))
}

#######################
## STATIC PARAMETERS ##
#######################

# Model Setup 

TimeStop <- 7*365 # length of a simulation (in days)
DigitsRound <- 0 # number of digits for rounding ts
ts <- 0 # first timestep
N <- 1e4 # population size
states <- c("S","I","R","V")

# Seed Infections

seeds <- 0.01*N # number of initial infections
I0 <- sample(1:N, seeds) # ids of first infected units

# Set vectors to store start state transitions:
I_time <- rep(Inf, length = N) # (S-->I) time of infective contact, i.e. beginning of infectious period. Set to Inf for units which are not infected.
I_time[I0] <- ts # set contact time for seed infections to 0 (ts)
IS_time <- rep(0, length = N) # (R-->S) time of return to susceptibility for infected units, set to 0 for units which are never infected
VS_time <- rep(0, length = N) # (V-->S) time of return to susceptibility for vaccinated units, set to 0 if never vaccinated


# Vaccination Parameters

V_frequency <- 1 # interval between vaccination rounds (if rounds >1)
V_rounds <- 1 # number of vaccination rounds
V_start <- round(365*3,digits=0) # start time of vaccination
# V_prop <- 0 # proportion vaccinated if not defined above
# V_exp <- 0 # vaccine strategy if not defined above

# store vaccination schedule:
V_schedule <- seq(from = V_start, # vaccination start generation
                  to = V_start + ((V_rounds - 1) * V_frequency), # start generation + number of rounds* frequency of rounds
                  by = V_frequency) # frequency of vaccination

# Duration of vaccination: 
V_period <- Inf # set to duration of simulation. 
if (V_period == Inf) {
  V_period <- TimeStop + 1
}