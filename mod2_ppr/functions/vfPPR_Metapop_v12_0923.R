## Simulate PPRV spread between units in metapopulation

PPR_Metapop_Heterogeneity <- function(
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
{
	
Remain_1_dt <- 1-Exit_1_dt #“retention rate”, 1 – rate of exit from age group 1 (0-1y)
Remain_2_dt <- 1-Exit_2_dt #“retention rate”, 1 – rate of exit from age group 2 (>1y)

Surv_PPR_1 <- 1 - Mort_PPR # PPR survival (same for groups 1 & 2)
Surv_PPR_2 <- 1 - Mort_PPR

# v_counter <- 1

I_VTOT <- rep(0,n) # total infected animals in units in metapop?

Decrease <- rep(0,n)
Increase <- rep(0,n)

## Allocating animals
Npop_1 <- Npop * PropAge_1 # number of animals in age group 1
Npop_2 <- Npop * PropAge_2 # number of animals in age group 2

## Number of births
nBirths <- Birth * Npop_2 # number of births per timestep (or period?), birth rate * number mature animals

# state vectors contain number of animals in each state, within each unit in the metapop
# all units are susceptible at time = 0
S_1 <- rep(Npop_1,n) 
S_2 <- rep(Npop_2,n)
I_1 <- rep(0,n)
I_2 <- rep(0,n)
R_1 <- rep(0,n)
R_2 <- rep(0,n)
#V <- rep(0,n)
V_1 <- rep(0,n)
V_2 <- rep(0,n)

## Seeding the infection
v_i <- sample(n,nSeeds,replace=FALSE) # generate ids of seed units in metapop
# for seeded sub-populations, update update S2 and I2 (seeds only within adult units)
I_2[v_i] <- NInfSeededArea # set number of infections in each seeded subpop to number specified by NInfSeededArea
S_2[v_i] <- S_2[v_i]-NInfSeededArea # update number of susceptibles in seeded sub-pops (minus number of seed infections)
S_2[v_i][ S_2[v_i] < 0 ] <- 0 # if susceptible population <0 reset to 0

## Storing results
# construct matric to store number of animals in each state in entire metapopulation
# X_1 + X_2 gives sum of animals in X state in all sub-populations for age group 1 and 2
mRes <- matrix(nrow=Timesteps+1,ncol=4)
colnames(mRes) <- c("S","I","R", "V")
mRes[1,"S"] <- sum(S_1+S_2) 
mRes[1,"I"] <- sum(I_1+I_2)
mRes[1,"R"] <- sum(R_1+R_2)
mRes[1,"V"] <- sum(V_1+V_2)
#---------------------------------------------------------
# Betwreen-unit transmission model () ##
#---------------------------------------------------------
## START: loop
for(Timestep in 1:Timesteps){ # for current timestep in metapopulation simulation

  # total number of infected/all units in each  meta-population
	I_TOT <- 0 # set total infected to 0 
	N_TOT <- 0 # set pop to 0?
	
	#---------------------------------------------------------
	# Vaccination ##
	#---------------------------------------------------------
	
	# Sample units for vaccination 
	# If vaccination timestep -->
	
	if(Timestep %in% vaccine_params[,"v_start"]){
	  v_round <- vaccine_params %>% filter(v_start == Timestep) %>% pull(rounds)
	  v_group <- vaccine_params %>% filter(v_start == Timestep) %>% pull(v_group)
	  # first vaccination timestep --> Sample units for vaccinating 
	  if(v_round ==1 & Timestep==min(vaccine_params[,"v_start"])){
	    wVacc_all <- vaccination(vBeta_b, vExp, v_units, v_strat) 
	    wVacc_split <- sample(wVacc_all, replace = F) %>% # sample the list of vaccine IDs to randomise
	      split(.,1:12)
	      
	    # split vaccine ids into groups for vaccination over 12 timesteps
	    wVacc_A <- wVacc_split$`1` 
	    wVacc_B <- wVacc_split$`2` 
	    wVacc_C <- wVacc_split$`3`
	    wVacc_D <- wVacc_split$`4`
	    
	    wVacc_E <- wVacc_split$`5` 
	    wVacc_F <- wVacc_split$`6` 
	    wVacc_G <- wVacc_split$`7`
	    wVacc_H <- wVacc_split$`8`
	    
	    wVacc_I <- wVacc_split$`9` 
	    wVacc_J <- wVacc_split$`10` 
	    wVacc_K <- wVacc_split$`11`
	    wVacc_L <- wVacc_split$`12`
	  }
	  ## END vaccination sampling
	  
	  ## Proportion vaccinated: 
	  
	  pA <- vaccine_params %>% filter(v_start == Timestep) %>% pull(pA)
	  pY <- vaccine_params %>% filter(v_start == Timestep) %>% pull(pY)
	  
	  if(v_group == "A"){
	    wVacc <- wVacc_A
	  }else if(v_group == "B"){
	    wVacc <- wVacc_B
	  }else if(v_group == "C"){
	    wVacc <- wVacc_C
	  }else if(v_group == "D"){
	    wVacc <- wVacc_D
	  }else if(v_group == "E"){
	    wVacc <- wVacc_E
	  }else if(v_group == "F"){
	    wVacc <- wVacc_F
	  }else if(v_group == "G"){
	    wVacc <- wVacc_G
	  }else if(v_group == "H"){
	    wVacc <- wVacc_H
	  }else if(v_group == "I"){
	    wVacc <- wVacc_I
	  }else if(v_group == "J"){
	    wVacc <- wVacc_J
	  }else if(v_group == "K"){
	    wVacc <- wVacc_K
	  }else if(v_group == "L"){
	    wVacc <- wVacc_L
	  }
	  
	  # V <- S+I+R
	  V_1[wVacc] <- V_1[wVacc]+(pY*(S_1[wVacc]+I_1[wVacc]+R_1[wVacc]))
	  V_2[wVacc] <- V_2[wVacc]+(pA*(S_2[wVacc]+I_2[wVacc]+R_2[wVacc]))
	  
	  S_1[wVacc] <- ((1-pY)*(S_1[wVacc]))
	  S_2[wVacc] <- ((1-pA)*(S_2[wVacc]))
	  
	  I_1[wVacc] <- ((1-pY)*(I_1[wVacc]))
	  I_2[wVacc] <- ((1-pA)*(I_2[wVacc]))
	  
	  R_1[wVacc] <- ((1-pY)*(R_1[wVacc]))
	  R_2[wVacc] <- ((1-pA)*(R_2[wVacc]))
	}

	#---------------------------------------------------------
	# Update numbers of animals in each state in Metapop ##
	#---------------------------------------------------------
	
	#pX_1, previous (before aging) S,I,R
	# X_2, X_2 contain number of units in X state in each subpopulation (length = no. sub pops in metapop)
	pS_1 <- S_1 
	pS_2 <- S_2
	pI_1 <- I_1
	pI_2 <- I_2
	pR_1 <- R_1
	pR_2 <- R_2
	pV_1 <- V_1
	pV_2 <- V_2
		
	# nX_temp = total number of units in X state (both age groups) in each sub-pop (length = no. sub pops)
	nS_Temp <- pS_1 + pS_2
	nI_Temp <- pI_1 + pI_2
	nR_Temp <- pR_1 + pR_2
	nV_Temp <- pV_1 + pV_2
	nN_Temp <- nS_Temp + nI_Temp + nR_Temp + nV_Temp # population size of each sub-population

	#identify which subpopulations are infected
	wInfe <- nI_Temp>0 
	wSusc <- nI_Temp==0

	#update total number of infected animals in all units in metapopulation
	I_TOT  <- sum(nI_Temp)
	#update total number of animals in all units in metapopulation
	N_TOT  <- sum(nN_Temp)

	# Replace nR with 0 if number of animals in R state is <MinNoR in a given unit
	pR_1[nR_Temp<MinNoR] <- 0
	pR_2[nR_Temp<MinNoR] <- 0
		
	# Logical vector of unit infection status. the same as wInfe, line 96.
	Inf_area <- nI_Temp>0
  
	#---------------------------------------------------------
	# Demographic processes ##
	#---------------------------------------------------------
	
	# Survival rate including ageing by disease state and age group, for nimals within each unit in metapopulation
	aS_1 <- Remain_1_dt * pS_1 # survival in age group 1
	Age_S_1 <- aS_1 * Age_r # (number of) units aging from age group 1 - 2
	aS_1 <- aS_1 - Age_S_1 # update number of units in S_1 after survival & aging
	aS_2 <- Remain_2_dt * pS_2 + Age_S_1 # S_2 = survival of units in S_2 + units moving into S_2 from S1
	
	# as above for infected units
	aI_1 <- Remain_1_dt * pI_1
	Age_I_1 <- aI_1 * Age_r
	aI_1 <- aI_1 - Age_I_1
	aI_2 <- Remain_2_dt * pI_2 + Age_I_1
	
	# as above for recovered units
	aR_1 <- Remain_1_dt * pR_1
	Age_R_1 <- aR_1 * Age_r
	aR_1 <- aR_1 - Age_R_1
	aR_2 <- Remain_2_dt * pR_2 + Age_R_1
	
	# as above for vaccinated units
	aV_1 <- Remain_1_dt * pV_1
	Age_V_1 <- aV_1 * Age_r
	aV_1 <- aV_1 - Age_V_1
	aV_2 <- Remain_2_dt * pV_2 + Age_V_1
	
	#	Number of reproductive animals in each unit. (age G2 and not infected)
	pN_Reprod <- pS_2 + pR_2 + pV_2

	#---------------------------------------------------------
	# Transmission processes within units ##
	#---------------------------------------------------------
	
	# Infection rate for animals within units
	ProbaInf <- 1 - exp(-Beta_w*nI_Temp/nN_Temp)

	# numbre of newly infected animals in each unit in metapopulation
	NewInf_1 <- ProbaInf * aS_1
	NewInf_2 <- ProbaInf * aS_2

	# Update number of susceptibles in each unit
	aS_1 <- aS_1 - NewInf_1
	aS_2 <- aS_2 - NewInf_2

	## BIRTHS
	# For susceptible units (subpops) all mature animals can reproduce (birth_rate*Npop_1)
	aS_1[wSusc] <- aS_1[wSusc] + nBirths
	# For infected units births are reduced as only healthy animals can reproduce (birth_rate*(S_2+R_2+V_2))
	aS_1[wInfe] <- aS_1[wInfe] + nBirths * pN_Reprod[wInfe] / Npop_2
	
	## Setting the number of infected animals to 0 at the end of the epidemic
	## When below a given threshold, or following the epidemic, at the trough
	PrevI <- nI_Temp
	NewI  <- NewInf_1 + NewInf_2 # total number of new infections in each sub-population

	# which populations have infections < than threshold infection level, set to 0
	w <- NewI<MinNoI 
	Decrease[w] <- 0
	Increase[w] <- 0
	NewInf_1[w] <- 0
	NewInf_2[w] <- 0

	Decrease[PrevI>NewI & Increase==1] <- 1
	Increase[PrevI<NewI] <- 1

	w <- Increase==1 & Decrease==1 & PrevI<NewI
	Decrease[w] <- 0
	Increase[w] <- 0
	NewInf_1[w] <- 0
	NewInf_2[w] <- 0
	
  # Recovery
	# Update number of recovered animals in each unit (from infected animals*PPR survival rate)
	aR_1 <- aR_1 + Surv_PPR_1 * aI_1
	aR_2 <- aR_2 + Surv_PPR_2 * aI_2

	# update number of infections 
	aI_1 <- NewInf_1
	aI_2 <- NewInf_2
			
	# Update number of animals in each state for all unit in metapopulation
	S_1 <- aS_1
	S_2 <- aS_2

	I_1 <- aI_1
	I_2 <- aI_2

	R_1 <- aR_1
	R_2 <- aR_2
	
	V_1 <- aV_1
	V_2 <- aV_2

	#---------------------------------------------------------
	# Transmission processes BETWEEN units ##
	#---------------------------------------------------------
	
	# Identify susceptible units 
	vSusc <- S_1 + S_2 # Total susceptible animals in each unit
	wSusc <- which((I_1+I_2)==0) # Units in susceptible state (no infections)
	nSusc <- length(wSusc) # Number of susceptible units
	
	if(nSusc>0){
	  # total force of infection for all units
		foi_TOT <- (sum(vBeta_b*nI_Temp)-vBeta_b*nI_Temp) * vExp / (N_TOT-nN_Temp)
		
		# Infection risk of susceptible units
		ProbaInf_b <- 1 - exp( - vSusc[wSusc] * foi_TOT[wSusc] )
		
		# random number allocated to each susceptible sub-pop to determine new infections
		RndNumber <- runif(nSusc,0,1)
		# generate new infections within susceptible subpopulations depending on probability of infection
		# if random number < infection risk then seed new infections.
		I_2[wSusc][RndNumber<ProbaInf_b] <- NInfSeededArea 	# update number of infected animals in new infected units
		S_2[wSusc][RndNumber<ProbaInf_b] <- S_2[wSusc][RndNumber<ProbaInf_b]-NInfSeededArea # update number of susceptible animals in new infected units
		
		S_2[wSusc][RndNumber<ProbaInf_b][ S_2[wSusc][RndNumber<ProbaInf_b] < 0 ] <- 0 #????????
	}

	#---------------------------------------------------------
	# Reseeding in metapopulation ##
	#---------------------------------------------------------
	## Re-seeding of infected units at metapopulation level
	# if timestep = timestep of seed events and < reseed threshold
	if(Timestep%%Freq_ReIntr==0 & Timestep<TimeLimit_ReIntr){
	  # which populations are susceptible?
		wSusc <- which((I_1+I_2)==0)
		# number of introductions limited by number of susceptible populations
		nSusc_Intr <- min(length(wSusc),nReIntr)
		# if there is more than 1 susceptible population in the metapop then continue with reseeding.
		if(nSusc_Intr>1){
		  
		  #generate ids of seed subpopulations
			wIntr <- sample(wSusc,nSusc_Intr,replace=FALSE)
			
			#for seed subpops, update I and S with number of new infections
			I_2[wIntr] <- NInfSeededArea
			S_2[wIntr] <- S_2[wIntr]-NInfSeededArea
			S_2[wIntr][ S_2[wIntr] < 0 ] <- 0
		}
	}

	#update mRes matrix tracking total number of units in each state in the metapopulation
	mRes[Timestep+1,"S"] <- sum(S_1+S_2)
	mRes[Timestep+1,"I"] <- sum(I_1+I_2)
	mRes[Timestep+1,"R"] <- sum(R_1+R_2)
	mRes[Timestep+1,"V"] <- sum(V_1+V_2)

}
## End: day loop

mRes

}
## end function

# plot(mRes[,2]/rowSums(mRes) , type="l" )
# plot(mRes[,3]/rowSums(mRes) , type="l" )


