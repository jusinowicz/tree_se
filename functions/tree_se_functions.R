#=============================================================================
# Functions to run the storage effect code. This includes population models, 
# invasion growth rate scenarios, data parsing. 
#=============================================================================
#libraries: 
library(MASS)

#=============================================================================
# Random reproduction rates. This is just to randomly generate some 
# reproduction with randomly generated parameters (means, variance, 
# correlation)
#=============================================================================

random_repro = function(nspp, time, mus = NULL, cov = NULL) {
	if(is.null(mus) == TRUE ) {
		#Randomly generate species repro
		mu1 = abs(rnorm(nspp,0.1 ) ) #Mean
	}else{ 
		mu1 = mus
	}

	if(is.null(cov) == TRUE) {
		#Create a random covariance matrix: 
		# Random correlation matrix
		cm1 = matrix(runif(nspp^2), ncol = nspp)
		cm1 =  cm1 %*% t(cm1)
		diag(cm1) = 1  # Diagonal elements should be 1

		# Standard deviations of the time series
		sd1 = runif(nspp)

		# Create the covariance matrix
		sig1 = diag(sd1) %*% cm1 %*% diag(sd1)
	}else{ 
		sig1 = cov
	}

	#use mvrnorm to make time series (lognormal)
	sl = exp(mvrnorm(time, mu1, sig1))

	return(sl)

}

#=============================================================================
#Population models
#=============================================================================
#
#=============================================================================
# Usinowicz forest model: 
# Two-stage seedling, sapling and adult model from Usinowicz et al. (2012,2017)
# Returns a single step of the population dynamics as a list. List items are 
# Adults, Saplings, Adult competition, Sapling competition. 
# The competition is returned because it is potentially useful for calculating
# various coexistence metrics.
#=============================================================================
ufm_step = function ( ni, si, repro, aij, fi, del1_s ) { 

	#Run each species simultaneously
	si_comp = 1+sum(aij%*%(ni*repro) ) #Sapling competition
	ni_comp = sum(repro) #Adult competition 
	si_new = (si*fi+ni*repro/si_comp) # Seedlings
	ni_new =  ni*del1_s+(1-sum(ni*del1_s)) * si/ni_comp #Adults

	#Put these into one variable to return
	vr = list(ni = ni_new,
			  si = si_new,
			  comp_ni = ni_comp,
			  comp_si = si_comp)
	return(vr)

}

#=============================================================================
# Stump tree model: from Stump and Vasseur 2023. 
# Returns a single step of the population dynamics as a list. List items 
# are Adults, Seedling competition, Attack rate. 
# The competition is returned because it is potentially useful for calculating
# various coexistence metrics.
#=============================================================================
stm_step = function (ni, repro, att_e, fi, del1) { 

	#Run each species simultaneously
	pij = exp ( - att_e *ni ) #Herbivore attack rate
	ni_comp = sum(ni * repro* pij ) #Total apparent competition
	ni_new = ni*(1-del1+del1*(repro*pij/ni_comp) ) 
	
	#Put these into one variable to return
	vr = list( ni = ni_new, 
			   comp_ni = ni_comp, 
			   pij = pij )
	return(vr)

}

#=============================================================================
# simulate_model
# Function to simulate a number of time steps of the model from random 
# initial conditions. 
#
# Must be told which model is being run: 
# 
# model 			Which population model. Currently either the Usinowicz 
#					forest model "ufm" or the Stump tree model "stm."
#
# parms				This is a list that must include the appropriate 
#					parameters for the model. 
# 					For ufm: 
#						aij		nspp x nspp matrix of competition coefficients
#						fi 		nspp vector of sapling survival
#						del1_s  nspp vector of adult survival 
#
#					For stm:
#						
#=============================================================================

simulate_model = function ( model, repro, nspp, parms, time,
							ics_yes=FALSE, ics_use=NULL ){

	#This is written with the option to pass ICs to facilitate 
	#the invasion scenarios. 
	if(ics_yes == TRUE) { 
		if( is.null(ics_use) == FALSE ){ 
			#Use the supplied ICs 
			ni = matrix(0,time,nspp)
			ni[1,] = ics_use

		} else {
			stop("You must supply the initial conditions")
		}

	} else { 

		#Initialize population matrix with a random number
		ni = matrix(runif(1), time, nspp)

	}

	#Which model is being run? 
	if( model == "ufm"){ 
		#Create the sapling stage as well 
		si = ni

		#Iterate over time steps
		for (t in 1:(time-1)){ 
			new_step = ufm_step(ni[t,], si[t,], repro[t,], 
						parms$aij, parms$fi, parms$del1_s)
			ni[t+1,] = new_step$ni #Next step in adult pop
			si[t+1,] = new_step$si #Next step in sapling pop

		}

		return(data.frame(ni = ni, si = si))

	} else if (model == "stm"){

		#Iterate over time steps
		for (t in 1:(time-1)){ 
			new_step = stm_step(ni[t,], repro[t,], 
						parms$att_e, parms$fi, parms$del1)
			ni[t+1,] = new_step$ni #Next step in adult pop
	
		}

		return(data.frame(ni = ni))

	} else { print("Error: model not recognized")}


}

#=============================================================================
# run_invasions
# Function to simulate population dynamics under invasion conditions for each 
# species against the community. 
#
# Must be told which model is being run: 
# 
# model 			Which population model. Currently either the Usinowicz 
#					forest model "ufm" or the Stump tree model "stm."
#
# parms				This is a list that must include the appropriate 
#					parameters for the model. 
# 					For ufm: 
#						aij		nspp x nspp matrix of competition coefficients
#						fi 		nspp vector of sapling survival
#						del1_s  nspp vector of adult survival 
#
#					For stm:
#
# pop_init			The equilibrium density of the community (use simulate_model)
#						
#=============================================================================

run_invasions = function ( model, repro, nspp, parms, time=500){

	#Initialize this list for each invasion scenario
	ni_all = vector("list", nspp)
	#Which model is being run? 
	if( model == "ufm"){ 

		#Take each species in turn as invader: 
		for (s in 1:nspp){
			#First, get equilibrium density of community in the 
			#absence of each species for ICs
			ni_inv = matrix(0,time, nspp )
 			si_inv = ni_inv

			ics_use = runif(nspp)
			ics_use[s] = 0
			ufm_sim = simulate_model ("ufm", repro, nspp, parms, time,
										ics_yes=TRUE, ics_use= ics_use )

			ni_inv[1,] = as.matrix(ufm_sim[time,1:nspp] )
			si_inv[1,] = as.matrix(ufm_sim[time,(nspp+1):(nspp*2)] )

			ni_inv[1,s] = .00000001
			si_inv[1,s] = .00000001
		
			#Iterate over time steps
			for (t in 1:(time-1)){ 
				new_step = ufm_step(ni_inv[t,], si_inv[t,], repro[t,], 
							parms$aij, parms$fi, parms$del1_s)
				ni_inv[t+1,] = new_step$ni #Next step in adult pop
				si_inv[t+1,] = new_step$si #Next step in sapling pop

			}
 
 			#Add each invasion to the list
			ni_all[[s]] = data.frame(ni = ni_inv, si = si_inv)

		}

		return(ni_all)

	} else if (model == "stm"){
		
		#Take each species in turn as invader: 
		for (s in 1:nspp){
			#First, get equilibrium density of community in the 
			#absence of each species for ICs
			ni_inv = matrix(0,time, nspp )

			ics_use = runif(nspp)
			ics_use[s] = 0
			stm_sim = simulate_model ("stm", repro, nspp, parms, time,
										ics_yes=TRUE, ics_use= ics_use )

			ni_inv[1,] = as.matrix(stm_sim[time,1:nspp] )
			ni_inv[1,s] = .00000001

			#Iterate over time steps
			for (t in 1:(time-1)){ 
				new_step = stm_step(ni_inv[t,], repro[t,], 
							parms$att_e, parms$fi, parms$del1)
				ni_inv[t+1,] = new_step$ni #Next step in adult pop
	
			}
			#Add each invasion to the list
			ni_all[[s]] = data.frame(ni = ni_inv)
		}

		return(ni_all)

	} else { print("Error: model not recognized")}


}

#=============================================================================
#get_ldgrs
# Function to measure the low-density growth rates. It finds all of the growth
# rate data below a certain arbtrary cut-off point, takes the log, finds the 
# slope. 
#=============================================================================
get_ldgrs = function (pop, invader, thresh = 0.001) { 
	pop_low = pop[pop[,invader]<thresh, ] #Filter
	nx = 1:dim(pop_low)[1]
	com_pop = rowSums(as.matrix(pop_low[,-invader]))#Sum residents
	log_inv_low = log(pop_low[,invader]) #Log
	log_com_pop = log(com_pop)
	ldgr_inv = lm(log_inv_low~nx)$coefficients[2] #Invader: Fit and get slope
	ldgr_res = lm(log_com_pop~nx)$coefficients[2] #Resident: Fit and get slope

	ldgr = NULL # This is the variable to return
	ldgr$inv = ldgr_inv 
	ldgr$res = ldgr_res
	return(ldgr)
}
