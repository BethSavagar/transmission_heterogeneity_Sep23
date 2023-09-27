
##############################################################
############# FIXED PARAMETERS: Do Not Modify ################
##############################################################

#-------------------------------------------------------------
# Model parameters: timesteps, seed infections, etc
#-------------------------------------------------------------

## Number of units (sub-populations)
n <- 10000

## Number of animals in metapopulation (number of animals within metapopulation?)
POPRegion_L <- 17366944

# Number of animals per unit (subpopulation)
Npop <- round(POPRegion_L/n)

## number of seed units (sub-populations)
nSeeds <- 0.01*n

## Number of animals infected in seed units (sub-pop)
NInfSeededArea <- 1

## Minimum number of infected/recovered animals within a unit (sub-population)
# how are these defined?
MinNoI <- 0.2
MinNoR <- 0.01

## Frequency of reseeding events in timesteps (*10 days)
Freq_ReIntr <- 50

## Upper limit for reseeding events in timesteps (*10 days)
TimeLimit_ReIntr <- 305

## Number of units infected during reseeding events
nReIntr <- 5



TSP_INFP <- 10 # infectious period for animal in unit / timestep in days

Tsp <- TSP_INFP # timestep in days

Day_L <- 365 # length of a year, in days
## Timesteps for within populations transmission model (nb. 1 timestep = 10 days)

TimeSteps_OnePop <- 200 # *10 days


#-------------------------------------------------------------
# Demographic parameters: ageing, survival etc ##
#-------------------------------------------------------------

## Proportion of adults and young within units (subpopulation)
P_A <- 0.6 
P_Y <- 1-P_A

## Proportion surviving from age group 1->2 in sub-populations per year?
PropSurv <- 0.5

## PPR mortality risk
Mort_PPR <- 0.5

AA <- 1 / (Day_L/Tsp) * PropSurv # survival rate per year
Birth <- P_Y/P_A / (Day_L/Tsp) #brith rate per year
Exit_2_dt <- AA*P_Y/P_A # rate at which animals leave group 2 per year, equal to rate at which animals move from group 1 - 2 * proportion in group 1 / proportin in group 2
Exit_1_dt <- 1/(Day_L/Tsp) - AA
Ageing_L <- 1 / ( Exit_2_dt/(1-Exit_1_dt)*P_A/P_Y ) * Tsp
AA - (1/(Ageing_L/Tsp)*(1-Exit_1_dt))
AA/(AA+Exit_1_dt)
Birth ; Exit_1_dt ; Exit_2_dt ; Ageing_L
AA/(AA+Exit_1_dt) ; 1/(Exit_2_dt*Day_L/Tsp) ; Birth*(Day_L/Tsp)/0.8

PropAge_1 <- P_Y
PropAge_2 <- P_A
Age_r <- Tsp / Ageing_L


# vaccination parameters script: 
# with vaccination spread over 4 weeks

v_units <- 1 # seq(0,1,0.1)

v_rounds <- 4
v_groups <- c("A","B","C","D")
v_adults <- c(1,1,0,0)
v_young <- c(1,1,1,1)
v_start <- seq(15,length.out = v_rounds, by=1)*36.5
v_start <- round(v_start, digits = 0) 
v_start_new <- c(v_start, v_start+1, v_start+2,v_start+3) %>% sort()
v_period <- Inf ; if(v_period == Inf){v_period <- 900}

v_strat <- 0 # c(-4,0,4)

vaccine_params <- data.frame(rounds = rep(1:v_rounds, each = 4),
                             v_group = rep(v_groups,4),
                             v_start = v_start_new,
                             v_end = v_start_new+v_period,
                             pA = rep(v_adults, each = 4),
                             pY = rep(v_young, each = 4)
)

