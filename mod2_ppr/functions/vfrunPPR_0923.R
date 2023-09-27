
## Run PPR metapopulation model with transmission within (homogeneous) and between (heterogeneous) units

#####################################################################################

runPPR <- function(Post, # posterior parameters:
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
){
  
  #---------------------------------------------------------
  # Identify seed unit + within-unit transmission params? ##
  #---------------------------------------------------------
  
  ## Sample parameter values from the joint posterior distribution
  # Identify seed unit in metapopulation and associated transmission potential and R0
  iS <- sample(nrow(Post),1,prob=Post$weight,replace=TRUE)
  Beta_w <- Post[iS,"Beta_w"] # within-population transmission rate
  r_b    <- Post[iS,"r_b"] # between-population reproduction number (posterior predictive value)
  
  
  #---------------------------------------------------------
  # Within-unit transmission model (homogeneous) ##
  #---------------------------------------------------------
  # simulate transmission within the seed unit: iS
  
  ## Return the sum of I/N (used to compute the force of infection and estimate Beta_b
  Res <- fPPR_OnePop(
                     n, # Number of units (sub-populations) in metapop
                     NInfSeededArea, # number of infected animals within seed unit
                     MinNoI, # min number infected animals per unit(sub-pop)
                     MinNoR, # min number recovered animals per unit(sub-pop)
                     TimeSteps_OnePop, # duration of within-unit transmission simulation (in timesteps)
                     
                     # Demographic parameters for animals within-unit:
                     Birth, # birth rate (per timestep?)
                     Exit_1_dt, # outflow g1
                     Exit_2_dt, # outflow g2
                     Age_r, # maturity rate (ageing)
                     Mort_PPR, # PPR mortality rate
                     PropAge_1, # proportion age g1
                     PropAge_2, # proportion age g2
                     
                     Npop, # number of animals within seed unit
                     Beta_w # transmission rate within seed unit
                     )
  
  ## Estimate Beta_b for each population based on r_b
  # beta_b is transmission rate between populations
  ## vExp: exposure
  if(Variable_Transmission=="no"){
    vr_b <- rep(r_b,n)
    vExp <- rep(1/n,n)
    vBeta_b <- -1/(Npop*Res*vExp[1]) * log(1-vr_b/(n-1))
  }
  
  if(Variable_Transmission=="yes"){
      vr_b <- rnbinom(n,size=k_estimate,mu=r_b)
      vExp <- rep(1/n,n)
      vBeta_b <- -1/(Npop*Res*vExp[1]) * log(1-vr_b/(n-1))
  }
  
  # k<-1 ; sum(1-exp(-vBeta_b[k]*Npop*Res*vExp[-k]))
  
  # mRes is number of animals in each state in metapopulation over time (where timestep = infectious period of unit)
  
  mRes <- PPR_Metapop_Heterogeneity(
    n, # Number of populations
    nSeeds, # number of populations where the infection is seeded
    NInfSeededArea, # number of new infections in a population (when incursion)
    MinNoI, # minimum number of infected animals, for each population 
    MinNoR, # minimum number of recovered animals, for each population 
    Birth, # birth rate (per PERIOD/timestep)
    Exit_1_dt , # rate at which animals leave the "0-1 year old" compartment (per PERIOD/timestep)
    Exit_2_dt , # rate at which animals leave the ">1 year old" compartment (per PERIOD/timestep)
    Age_r , # Ageing: rate at which units move from group 1 - 2 in a given timestep
    Mort_PPR , # mortality risk from PPR infection
    PropAge_1, # proportion of animals in 1st age category
    PropAge_2, # proportion of animals in 2nd age category
    Npop, # size of a population (number of animals in a given sub-population)
    Beta_w, # transmission rate (within a sub-population)
    vBeta_b, # transmission rate (between sub-populations)
    vExp, # Exposure probability (between populations)
    Timesteps, # Number of timesteps to simulate transmission within metapopulation
    Freq_ReIntr, # Frequency of re-introduction (of seed infections in metapopulation)
    TimeLimit_ReIntr, # Time limit for re-introduction (of seed infections in metapop)
    nReIntr, # Number of populations within which infected animals are re-introduced
    
    v_units, # vaccination coverage of units in metapopulation
    v_strat, # c(-4,0,4) targeting of units in metapopualtion
    vaccine_params # parameters for implementing vaccination within units
    )
  
  # Output 
  
  if(output == "counts"){
    return(mRes)
  }else if(output == "I_counts"){
    return(mRes[,"I"]/rowSums(mRes))
  }
  
}

# plot(mRes[,"I"]/rowSums(mRes) , type="l" )
# plot(mRes[,"R"]/rowSums(mRes) , type="l" )




