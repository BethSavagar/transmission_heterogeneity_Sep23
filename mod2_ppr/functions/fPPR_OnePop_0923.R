## Within-unit (sub-population) spread of PPR

fPPR_OnePop <- function(
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
{

Timestep_Intr <- 1 # ?? Time of introduction

# -----------------------------------------------------------------------------
## Demographic Rates Calculations
# -----------------------------------------------------------------------------

Remain_1_dt <- 1-Exit_1_dt # survival rate age g1
Remain_2_dt <- 1-Exit_2_dt # survival rate age g2

# PPR survival rate. equal for g1 and g2
Surv_PPR_1 <- 1 - Mort_PPR
Surv_PPR_2 <- 1 - Mort_PPR

# Population size per age group
Npop_1 <- Npop * PropAge_1
Npop_2 <- Npop * PropAge_2

# Number of births per timestep (line 114?)
nBirths <- Birth * Npop_2 # absolute number of births, calculated from birth rate * number of mature animals (>1y)


# -----------------------------------------------------------------------------
## Disease state of population, by age group
# -----------------------------------------------------------------------------

# Susceptible (all animals)
S_1 <- Npop_1
S_2 <- Npop_2

# Infected
I_1 <- 0 ; I_2 <- 0

# Recovered animals
R_1 <- 0 ; R_2 <- 0

# -----------------------------------------------------------------------------

Increase <- 0 ; # ?? What is this (BS - 280721)
Decrease <- 0 ; # ?? What is this (BS - 280721)


Res <- 0 ; # ?? What is this (BS - 280721)

#-------------------------------------------------------------
## Transmission Simulation ##
#-------------------------------------------------------------

for(Timestep in 1:TimeSteps_OnePop){

	#-------------------------------------------------------------
	## Demographic processes ##
  #-------------------------------------------------------------
  
  # pX = temporary number of animals in each state/age group
	pS_1 <- S_1
	pS_2 <- S_2
	pI_1 <- I_1
	pI_2 <- I_2
	pR_1 <- R_1
	pR_2 <- R_2

	nS_Temp <- pS_1+pS_2 # total number of susceptibles
	nI_Temp <- pI_1+pI_2 # total number of infecteds
	nR_Temp <- pR_1+pR_2 # total number of recovereds
	nN_Temp <- nS_Temp+nI_Temp+nR_Temp # total population size
			
	if(nR_Temp<MinNoR){ pR_1 <- 0 ; pR_2 <- 0 } # if the number of recovered animals in the unit is < than the specified minimum, then set to 0 for each age group
  
	# Survival rate including ageing by disease state and age group
	
	# Susceptibles:
	aS_1 <- Remain_1_dt * pS_1 # g1 survival rate 
	Age_S_1 <- aS_1 * Age_r # g1 ageing
	aS_1 <- aS_1 - Age_S_1 # g1 loss due to ageing
	aS_2 <- Remain_2_dt * pS_2 + Age_S_1 # g2 survival rate and gain due to ageing (from g1-g2)
	
	# infected:
	aI_1 <- Remain_1_dt * pI_1
	Age_I_1 <- aI_1 * Age_r
	aI_1 <- aI_1 - Age_I_1
	aI_2 <- Remain_2_dt * pI_2 + Age_I_1
  
	# Recovered
	aR_1 <- Remain_1_dt * pR_1
	Age_R_1 <- aR_1 * Age_r
	aR_1 <- aR_1 - Age_R_1
	aR_2 <- Remain_2_dt * pR_2 + Age_R_1
	
	
	# Reproductive population (Susceptible & Recovered in age g2)
	pN_Reprod <- pS_2 + pR_2 # Birth added after the Infection process 

	
	#-------------------------------------------------------------
	## Infection processes within seed unit ##
	#-------------------------------------------------------------
	
	# Infection rate within-unit
	ProbaInf <- 1 - exp(-Beta_w * nI_Temp/nN_Temp)
	
	# calculate the NUMBER of new infections in each age category (from susceptible population)
	NewInf_1 <- ProbaInf * aS_1
	NewInf_2 <- ProbaInf * aS_2
	
	# update number of susceptible units in each age category, minus new infections
	aS_1 <- aS_1 - NewInf_1
	aS_2 <- aS_2 - NewInf_2

	## BIRTHS
	## could just do birth rate * pN_Reprod?
	aS_1 <- aS_1 + nBirths * pN_Reprod / Npop_2

	## Setting the number of infected animals to 0 at the end of the epidemic
	## When below a given threshold, or following the epidemic, at the trough
	PrevI <- nI_Temp
	NewI <- NewInf_1 + NewInf_2
	
	# If number of new infections is under specified threshold (minNoI) then set to 0
	## What are incease and decrease??
	if(NewI<MinNoI){
		Decrease <- 0
		Increase <- 0
		NewInf_1 <- 0
		NewInf_2 <- 0
	}

	if(PrevI > NewI & Increase==1 ){ Decrease <- 1 }
	if(PrevI < NewI){ Increase <- 1 }

	## Switch, increase, decrease and then 0 at the next increase
	if(Increase==1 & Decrease==1 & PrevI<NewI){
		Decrease <- 0
		Increase <- 0
		NewInf_1 <- 0
		NewInf_2 <- 0
	}
	
	## UPDATE SYSTEM
	# update number of units in each state after infection adn births are accounted for

	aR_1 <- aR_1 + Surv_PPR_1 * aI_1 # number recoveredm = number in recovered category + 
	aR_2 <- aR_2 + Surv_PPR_2 * aI_2

	aI_1 <- NewInf_1
	aI_2 <- NewInf_2
	
	S_1 <- aS_1
	S_2 <- aS_2

	I_1 <- aI_1
	I_2 <- aI_2

	R_1 <- aR_1
	R_2 <- aR_2
	
	if(Timestep==Timestep_Intr){
		I_2 <- NInfSeededArea
		S_2 <- max(0,S_2-NInfSeededArea)
	}

	## End: Ageing/Infection process
	
	Res <- Res + (I_1+I_2) / ( S_1+S_2+I_1+I_2+R_1+R_2 + Npop*(n-2) )
	# BS 03 Aug - this Res gives the proportion of infected individuals in the metapopulation
	# cumulative for the epidemic period (hence Res + I/N)
	# + Npop*(n-2), number of units in the metapopulation (i.e. sum of all units in populations excluding index infected population and current infected population?)

}
## End: Timestep loop

Res # returns I/N where I is all infected individuals in the unit over unit-level timestep (200days) and N is the size of the metapopulation?

}
## end function



