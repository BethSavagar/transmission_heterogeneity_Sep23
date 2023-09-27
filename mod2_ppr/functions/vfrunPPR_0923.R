
## Function to run PPR metapopulation model: tranmission between units (heterogeneous)

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
  Res <- fPPR_OnePop(n,
                     NInfSeededArea, # number of animals infected within seed unit
                     MinNoI,
                     MinNoR,
                     TimeSteps_OnePop,# number of timesteps for within unit transmission
                     
                     # within-unit demographic parameters:
                     Birth,
                     Exit_1_dt,
                     Exit_2_dt,
                     Age_r,
                     Mort_PPR,
                     PropAge_1,
                     PropAge_2,
                     
                     Npop, # number of animals within unit
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
  
  #mRes is number of animals in each state in metapopulation over time (where timestep = infectious period of unit)
  
  mRes <- PPR_Metapop_Heterogeneity(n, 
                                    nSeeds, 
                                    NInfSeededArea, 
                                    MinNoI, 
                                    MinNoR, 
                                    Birth,
                                    Exit_1_dt, 
                                    Exit_2_dt, 
                                    Age_r, 
                                    Mort_PPR, 
                                    PropAge_1, 
                                    PropAge_2, 
                                    Npop, 
                                    Beta_w, 
                                    vBeta_b,
                                    vExp, 
                                    Timesteps, 
                                    Freq_ReIntr, 
                                    TimeLimit_ReIntr, 
                                    nReIntr,
                                    
                                    v_units, # vaccination coverage of units in metapopulation
                                    v_strat, # c(-4,0,4) targeting of units in metapopualtion
                                    vaccine_params # parameters for implementing vaccination within units
                                    )
  
  if(output == "counts"){
    return(mRes)
  }else if(output == "I_counts"){
    return(mRes[,"I"]/rowSums(mRes))
  }
  
}

# plot(mRes[,"I"]/rowSums(mRes) , type="l" )
# plot(mRes[,"R"]/rowSums(mRes) , type="l" )




