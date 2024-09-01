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
	si_comp = 1+aij%*%(ni*repro)  #Sapling competition
	ni_comp = sum(si) #Adult competition 
	si_new = (si*fi+ni*repro/si_comp) # Seedlings
	ni_new =  ni*del1_s+(1-sum(ni*del1_s)) * si/ni_comp #Adults

	# si_comp = 1+aij%*%(ni[t,]*repro[t,])  #Sapling competition
	# ni_comp = sum(si[t,]) #Adult competition 
	# si_new = (si[t,]*fi+ni[t,]*repro[t,]/si_comp) # Seedlings
	# ni_new =  ni[t,]*del1_s+(1-sum(ni[t,]*del1_s)) * si[t,]/ni_comp #Adults


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
 			si_comp = ni_inv

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
				si_comp[t, ] = new_step$comp_si
			}
 
 			#Add each invasion to the list
			ni_all[[s]] = data.frame(ni = ni_inv, si = si_inv, si_comp = si_comp)

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
	#print(paste("thresh ", thresh, "nx ", length(nx)))
	com_pop = rowSums(as.matrix(pop_low[,-invader]))#Sum residents
	log_inv_low = log(pop_low[,invader]) #Log
	log_com_pop = log(com_pop)
	ldgr_inv = exp(lm(log_inv_low~nx)$coefficients[2]) #Invader: Fit and get slope
	ldgr_res = exp(lm(log_com_pop~nx)$coefficients[2] )#Resident: Fit and get slope

	ldgr = NULL # This is the variable to return
	ldgr$inv = ldgr_inv 
	ldgr$res = ldgr_res
	return(ldgr)
}


#=============================================================================
#get_dij
# Version 1: LDG only 
# Function to measure the low-density growth rates, in the form of Dij from 
# Usinowicz et al. 2012
# 
# Assume the data are pariwise interactions so xi, si, ri are each a matrix 
# with dimensions #timesteps by 2.  
#=============================================================================

get_dij = function ( xi, si, ri, di, fi,aij, invader=2){

	nspp = dim(xi)[2]
	aijm = matrix(aij[invader,],dim(xi)[1],nspp)
	si1 = fi[invader] * si[,invader] + ri[,invader]*xi[,invader]/(1+rowSums(aijm*xi*ri) )  
	Dij = mean(log(di[invader] + (1-di[invader])*(si1/(xi[,invader]*rowSums(si) ) ) ) ) 

	return(Dij)

}

#=============================================================================
#get_dij_approx1
# Version 2: Can we numerically approximate di -> 1? Not really
# Function to measure the low-density growth rates, in the form of Dij from 
# Usinowicz et al. 2012
# 
# Assume the data are pariwise interactions so xi, si, ri are each a matrix 
# with dimensions #timesteps by 2. 
#=============================================================================

get_dij_approx1 = function ( xi, si, ri, fi,aij, invader){

	#di=0.99999
	nspp = dim(xi)[2]
	aijm = matrix(aij[invader,],dim(xi)[1],nspp)
	si1 = fi[invader] * si[,invader] + ri[,invader]*xi[,invader]/(1+rowSums(aijm*xi*ri) )  
	#Dij = mean(log(di + (1-di)*(si1/(xi[,invader]*rowSums(si) ) ) ) ) 
	
	#This version matches the full approximation below. It is essentially
	#assuming that species achieve their MAX growth rates relative to the 
	#adult mortality rate, di. It is actually the correct way to do the 
	# di -> 1 approximation numerically:

	Dij = mean(log(1 + (si1/(xi[,invader]*rowSums(si) ) ) ) ) 


	return(Dij)

}

#=============================================================================
#get_dij_approx2
# Version 2: Calculate Dij based on the step after di -> 1. 
# Function to measure the low-density growth rates, in the form of Dij from 
# Usinowicz et al. 2012
#
# Assume the data are pariwise interactions so xi, si, ri are each a matrix 
# with dimensions #timesteps by 2. 
#=============================================================================


get_dij_approx2 = function( xi, si, ri, fi,aij, invader){

	nspp = dim(xi)[2]
	aijm = matrix(aij[invader,],dim(xi)[1],nspp)
	si1 = fi[invader] * si[,invader] + ri[,invader]*xi[,invader]/(1+rowSums(aijm*xi*ri) )  
	Dij = mean(log(si1/(xi[,invader]*rowSums(si) ) ))
	
	# The original version of the approximation was: 
	# Dij = mean((si1/(xi[,invader]*rowSums(si) )-1 ))
	# Which returns a larger value. Taking the log seems to be closer 
	# to what we get when measuring the LDG directly (i.e. get_dij).

	return(Dij)

}

#=============================================================================
#get_dij_approx3
# Version 3: The full approximation
# Function to measure the low-density growth rates, in the form of Dij from 
# Usinowicz et al. 2012
#
# Assume the data are pariwise interactions so xi, si, ri are each a matrix 
# with dimensions #timesteps by 2.  
#
# Notes: 
# The approximation x(t-tau)/x(t) = 1 at the boundary seems to hold well
# Approximating the numerator and the denominator through the sums seems to 
# hold well.
# Breaking the approximation down to si/sj at the boundary seems to hold well. 
# But using the original approximation doesn't seem to match with the first
# purely numerical definition used above!
# After some tinkering, this is what matches: 
# The below sum version of the approximation  = 
#				mean(log(1 + (si1/(xi[,invader]*rowSums(si) ) ) ) ) 
#=============================================================================

get_dij_approx3 = function(ri, fi, aij, tau=15, invader = 2){

	ngens = dim(ri)[1]
	nspp = dim(ri)[2]
	res = 1:nspp
	res = res[-invader]

	#fi^tau
	fi.tau = fi[invader]^seq(from=0, to=tau, by=1)
	fi.tau = fi.tau[length(fi.tau):1]

	#invader competition coefficients (remove invader)
	aijI = aij[invader, ]
	aijI = aijI[-invader]

	#resident competition coefficients (remove invader)
	aijR = aij[res,]
	aijR = aijR[-invader]

	sumDij = 0
	num1_k = matrix(0,ngens-tau,1)
	den1_k = matrix(0,ngens-tau,1)
	for (t in 2:(ngens-tau)){

		#Numerator of the approximation, sum from t-tau to t
		num1 = xi[t:(t+tau),invader]*ri[t:(t+tau),invader]/(1+aijI*ri[t:(t+tau),res])
		num1_k[t] = num1%*%fi.tau
		num1 = num1%*%fi.tau

		#Denominator of the approximation, sum from t-tau to t
		den1 = xi[t:(t+tau),invader]*ri[t:(t+tau),res]/(1+aijR*ri[t:(t+tau),res])
		den1_k[t] = den1%*%fi.tau
		den1 = den1%*%fi.tau

		#Keep a running sum
		sumDij = sumDij+(num1/den1)
	}

	#Get the expectation by dividing the length of the sum
	Dij = sumDij/(ngens-tau) -1 
	#Dij = log(Dij)

	return(Dij)

}