# Continuous IBM function for Model I (26.09.23)
# An agent-based model simulating transmission of an unspecified, directly transmitted virus in a heterogeneous population

continuousIBM <- function(
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
    
){
  
  ################################
  ## Transmission Distributions ##
  ################################
  
  # Generate transmission distribution and assign transmission potentials to units in the population:
  dists <- dist_transmission(N,m,k) # generate transmission potential distribution following the negative binomial
  T_dist <- dists[[1]] # transmission distribution
  E_dist <- dists[[2]] # exposure distribution (uniform)
  
  ############
  ## OUTPUT ##
  ############
  
  # Set up model output:
  
  if(output == "full"){
    # Matrix to store the state of each unit in the population over time
    mat <- as.data.frame(matrix(0,nrow = N, ncol = TimeStop+1)) # first col is t=0
  }else{
    # Matrix to store the number of SIRV over time
    counts <- as.data.frame(matrix(0,nrow = TimeStop+1, ncol = length(states)))
  }
  
  while(ts <= TimeStop){
    
    ############################################################################################################
    ## 1. Identify units which are infected at ts
    ############################################################################################################
    
    # Identify units which are contacted & become infected at time = ts
    I0 <- which(I_time == ts) # I0 contains unit IDs
    
    # Update IS_time vector to store end of infectious+immune period for current (+ previously) infected units. 
    # IS_time = 0 if unit is never infected
    IS_time[I0] <- ts+IR_period 
    
    # Reset I_time to Inf for units which are currently infected
    I_time[I0] <- Inf
    
    ############################################################################################################
    ## 2. Generate effective contact IDs:
    ############################################################################################################
    xInt <- 10^(-DigitsRound) 
    
    # number of effective contacts for current infected units:
    C_num <- T_dist[I0]
    
    if(sum (C_num)>0){ # if there are >0 effective contacts...
      
      # sample effective contacts IDs from population with probability E_dist
      ## N.B. for small populations, self-self contact should be excluded
      
      # for each value of C_num sample the population to generate that number of effective contact IDs, store in vector C_id
      C_id <- unlist( # as vector
        lapply(C_num, # iterate over values in C_num vector
               function(x) 
                 sample(
                   1:N, # sample entire population (* not valid for small pops)
                   x, # number of effective contacts per infected unit
                   replace = FALSE, # a given infected unit cannot contact the same unit twice
                   prob = E_dist))) # E_dist defines probability of exposure for each unit
      
      # C_id contains ids of all units contacted by current infections (including duplicates, and non-S states)
      
      
      ############################################################################################################
      ## 3. Generate time of contact IDs:
      ############################################################################################################
      
      # contact should occur during infectious period of current infected unit (ts+I_period)
      # runif samples infectious period to generate uniformly distributed contact times
      xInt <- 10^(-DigitsRound) # 1 unit of time, new infections can only begin 1 unit of time after primary infection
      C_time <- runif(length(C_id), # number of units to sample
                      ts + xInt, # start of infectious period 
                      ts + I_period # end of infectious period # -xInt added on 03/06 following meeting with GF
      )
      C_time <- round(C_time, digits = DigitsRound) # rounds the contact time to the number specified by DIgits Round, i.e. if DR = 1 then time will be specific to e.g. 0.3 of a day
      
      ############################################################################################################  
      ## 4. Clean Effective Contacts
      ############################################################################################################   
      
      # C_id and C_time contain the id and time of contact for all units which are effectively contacted by a unit which is infected at time ts
      # C_id and C_time need to be cleaned to contain only suscepitble units and unique contacts
      
      ## Identify Susceptible Effective Contacts
      
      # If unit is previously infected then C_time must occur after end of IR_period, stored in IS_time i.e. C_time > IS_time 
      ## N.B. IS_time = 0 for units which are not yet infected)
      # If unit is previously vaccinaged then C_time must occur after end of V_period, stored in VS_time i.e. C_time > VS_time 
      ## N.B. VS_time = 0 for units which are not yet vaccinated)
      # If unit is infected in future then C_time must occur before the recorded effective contact, i.e. C_time < I_time
      ## N.B. I_time = Inf from all other units
      # if unit is never infected then IS_time = 0 and I_time = Inf:  IS_time < C_time < I_time
      
      w <- which( C_time > IS_time[C_id] # unit is never infected (IS_time = 0) or is infected previously and has returned to S_state (IS_time < C_time)
                  &
                    C_time > VS_time[C_id] # unit is never vaccinated (VS_time = 0) or is vaccinated previously and has returned to S_state (VS_time < C_time)
                  &
                    C_time < I_time[C_id] # unit is never infected (I_time = Inf) or is infected in future but C_time occurs sooner ( C_time < I_time)
      )
      
      # Update C_time and C_id vectors to contains only susceptible units, and contacts which occur soonest.
      C_time <- C_time[w] 
      C_id <- C_id[w]
      
      
      ############################################################################################################
      ## 5. Remove duplicate contacts & update I_time vector
      ############################################################################################################
      
      # Sort contacts in order of increasing time, discard any duplicate contacts in order of increasing C_time
      
      # if there are any new effective contacts 
      if(length(C_time)>0){
        
        # sort out by time of infection
        C_idSort <- C_id[order(C_time)]
        C_timeSort <- sort(C_time)
        
        # remove duplicate contacts (keeping soonest contact time)
        C_idClean <- C_idSort[!duplicated(C_idSort)]
        C_timeClean <- C_timeSort[!duplicated(C_idSort)]
        
        # Update I_time with new effective contact times
        I_time[C_idClean] <- C_timeClean
        # I_time contains contact time of new infected units, all other units = Inf
        
      }
      
    }

    ############################################################################################################
    ## 6. Identify next infection
    ############################################################################################################     
    # Identify the timne of the next infection
    # times in I_time should always be higher than ts (up to Inf)
    
    tsNext <- min(I_time)
    
    t1 <- ceiling(ts)  # round up ts
    t2 <- floor(tsNext - xInt) # round down
    if (t2 == Inf) { t2 <-TimeStop } # If t2 is Inf, then no new infection - end simulation
    
    tv <- min(V_schedule[V_schedule >= t1]) # tv = the next timestep for vaccination
    
    ##################
    ## STATE UPDATE ##
    ##################
    
    while (tv <= t2) {
      # update the state of the system between t1 (current infection) to tv (vaccination)
      for (t_i in t1:tv) {
        
        # for timesteps between t1:tv
        if(t_i < tv){
          if (output == "counts" |
              output == "I_counts") {
          
            ## Infected:
            # if t_i is within infectious period of unit
            I_count <-
              sum(t_i >= IS_time - IR_period # from beginning of I_period (I_time = IS_time - IR_period)
                  &
                    t_i < IS_time - R_period) # to end of I_period (R_time = IS_time - R_period)
            
            ## Recovered:
            # If t_i is within recovered period of unit
            R_count <-
              sum(t_i >= IS_time - R_period # from beginning of R_period (R_time = IS_time - R_period)
                  &
                    t_i < IS_time) # to end of IR_period (IS_time = I_time + IR_period)
            
            ## Vaccinated:
            # If t_i is within vaccinated period of unit
            V_count <-
              sum(t_i >= VS_time - V_period # from beginning of V_period
                  & t_i < VS_time) # to end of V_period
            
            
            ## Susceptible:
            S_count <-
              sum(t_i >= IS_time # IS_time = 0 for units which are not infected or IS_time = I_time+IR_period for current & previous infections
                  &
                    t_i >= VS_time) # VS_time = 0 for units which are not vaccinated or VS_time = V_time + V_period for current & previous vaccinations
            
            # Update counts with number of units in each state at ts = 0
            ## t_i+1, since first ts = 0
            counts[t_i + 1, ] <- c(S_count, I_count, R_count, V_count)
          }
          
          if (output == "full") {
            I_index <-
              which(t_i >= IS_time - IR_period # from beginning of I_period (I_time = IS_time - IR_period)
                    &
                      t_i < IS_time - R_period) # to end of I_period (R_time = IS_time - R_period)
            
            R_index <-
              which(t_i >= IS_time - R_period # from beginning of R_period (R_time = IS_time - R_period)
                    &
                      t_i < IS_time) # to end of IR_period (IS_time = I_time + IR_period)
            
            V_index <-
              which(t_i >= VS_time - V_period # from beginning of V_period
                    & t_i < VS_time) # to end of V_period)
            
            S_index <-
              which(t_i >= IS_time # IS_time = 0 for units which are not infected or IS_time = I_time+IR_period for current & previous infections
                    &
                      t_i >= VS_time) # VS_time = 0 for units which are not vaccinated or VS_time = V_time + V_period for current & previous vaccinations
            
            
            mat[I_index, t_i + 1] <- "I"
            mat[R_index, t_i + 1] <- "R"
            mat[S_index, t_i + 1] <- "S"
            mat[V_index, t_i + 1] <- "V"
          }
        }
        if(t_i == tv){
          
          # at t_i == tv vaccinate units...
          
          #########################################################################################################
          ## 7. Vaccination
          #########################################################################################################
          
          # select units for vaccination with targetting according to value of V_exp
          newV_index <- vaccination(T_dist, E_dist, V_prop, V_exp)
          
          # update VS_time (end of V_period) for vaccinated units
          VS_time[newV_index] <- tv + V_period
          
          # reset IS_time to 0 for vaccinated units (IS_time contains return to S for previously infected units)
          IS_time[newV_index] <- 0
          
          # reset I_time ffor units which have an infectious contact during vaccination period
          w <- which(I_time[newV_index] >= tv
                     & I_time[newV_index] < tv + V_period)
          
          I_time[w] <- Inf
          
        }
        
      }
      
      #########################################################################################################
      ## 8. Update Model Parameters
      #########################################################################################################
      
      t1 <- tv # next timestep for system update
      tsNext <- min(I_time) # next infection (updated after vaccination)
      t2 <- floor(tsNext - xInt) # round down time of next infection to avoid double counting
      if (t2 == Inf) { t2 <- TimeStop} # If t2 is Inf, then no new infection - end simulation
      
      tv <- min(V_schedule[V_schedule > tv]) # select the next vaccination time
      # if no more vaccinations then tv <- Inf, tv >> t2 so will move on to next for loop
      
    }
    
    
    #########################################################################################################
    ## 9. Update State of System
    #########################################################################################################
    
    # for timesteps between t1:t2 or tv:t2 depending on vaccination
    
    if (t1 <= t2) {
      for (t_i in t1:t2) {
        
        if (output == "counts" |
            output == "I_counts") {
          ## Infected:
          # if t_i is within infectious period of unit
          I_count <-
            sum(t_i >= IS_time - IR_period # from beginning of I_period (I_time = IS_time - IR_period)
                & t_i < IS_time - R_period) # to end of I_period (R_time = IS_time - R_period)
          
          ## Recovered:
          # If t_i is within recovered period of unit
          R_count <-
            sum(t_i >= IS_time - R_period # from beginning of R_period (R_time = IS_time - R_period)
                & t_i < IS_time) # to end of IR_period (IS_time = I_time + IR_period)
          
          ## Vaccinated:
          # If t_i is within vaccinated period of unit
          V_count <-
            sum(t_i >= VS_time - V_period # from beginning of V_period
                & t_i < VS_time) # to end of V_period
          
          
          ## Susceptible:
          S_count <-
            sum(t_i >= IS_time # IS_time = 0 for units which are not infected or IS_time = I_time+IR_period for current & previous infections
                & t_i >= VS_time) # VS_time = 0 for units which are not vaccinated or VS_time = V_time + V_period for current & previous vaccinations
          
          # Update counts with number of units in each state at ts = 0
          ## t_i+1, since first ts = 0
          counts[t_i + 1, ] <- c(S_count, I_count, R_count, V_count)
        }
        
        if (output == "full") {
          I_index <-
            which(t_i >= IS_time - IR_period # from beginning of I_period (I_time = IS_time - IR_period)
                  &
                    t_i < IS_time - R_period) # to end of I_period (R_time = IS_time - R_period)
          
          R_index <-
            which(t_i >= IS_time - R_period # from beginning of R_period (R_time = IS_time - R_period)
                  & t_i < IS_time) # to end of IR_period (IS_time = I_time + IR_period)
          
          V_index <-
            which(t_i >= VS_time - V_period # from beginning of V_period
                  & t_i < VS_time) # to end of V_period)
          
          S_index <-
            which(t_i >= IS_time # IS_time = 0 for units which are not infected or IS_time = I_time+IR_period for current & previous infections
                  & t_i >= VS_time) # VS_time = 0 for units which are not vaccinated or VS_time = V_time + V_period for current & previous vaccinations
          
          
          mat[I_index, t_i + 1] <- "I"
          mat[R_index, t_i + 1] <- "R"
          mat[S_index, t_i + 1] <- "S"
          mat[V_index, t_i + 1] <- "V"
        }
        
      }
      
    }
    
    ############################################################################################################
    ## 10. Update timestep
    ############################################################################################################ 
    
    ts <- tsNext
    
  }
  
  
  ############################################################################################################
  ## Define Output
  ############################################################################################################ 
  
  if(output == "full"){
    # state of each individual over time
    mat <- mat[, 1:TimeStop]
    return(mat)
  }else if(output == "I_counts"){
    # number of units in each state over time
    counts <- counts[1:TimeStop, ]
    return(counts[,2])
  }else{
    # number of units in each state over time
    counts <- counts[1:TimeStop, ]
    return(counts)
  }
  
}




