#=============================================================================
# Main file for the simulation of population dynamics, invasion growth rates
# and measurement of storage effect. 
# 
# There are two different models at the moment:
# 
#=============================================================================
#load libraries 
source("./functions/tree_se_functions.R")

#=============================================================================
# Global parameters, including the reproduction 
#=============================================================================
ntime = 5000
nspp = 2

#Get random repro series for two species. 
repros_sp = random_repro (nspp,ntime)

#=============================================================================
# Initialize variables
#=============================================================================
#Model 1 adapted from Usinowicz et al. 2012
#parameters
parms_ufm = NULL
parms_ufm$aij = matrix(1,nspp,nspp) #Seedling competition
parms_ufm$fij = matrix(0.6,nspp,1) #Seedling survival
del1 = c(0.1,0.1) #Adult death rate
parms_ufm$del1_s = 1 - del1 #Adult survival rate

#Model 2 adapted from Stump and Vasseur 2023
#parameters
parms_stm = NULL
parms_stm$del1 = c(0.1,0.1) #Adult death rate
parms_stm$att_e = c(1,1) #Enemy attack rate ????<====Any expectation from lit?   

#=============================================================================
# Run simulations
#=============================================================================
ufm_sim = simulate_model ( model = "ufm", repro = repros_sp, nspp=nspp, 
			parms=parms_ufm, time= ntime)
stm_sim = simulate_model ( model = "stm", repro = repros_sp, nspp=nspp, 
			parms=parms_stm, time= ntime)

#Check simulations with basic plots
par(mfrow = c(2,1))
plot(ufm_sim[,1],t="l", ylim=c(0,1), main = "UFM")
for(n in 2:nspp){
	lines( ufm_sim[,n],col="blue")
}


plot(stm_sim[,1],t="l", ylim=c(0,1), main = "STM")
for(n in 2:nspp){
	lines( stm_sim[,n],col="blue")
}

#=============================================================================
# Run the invasions to measure the storage effect. 
#=============================================================================
ufm_inv_full = run_invasions( model = "ufm", repro = repros_sp, nspp=nspp, 
			parms=parms_ufm, time= ntime)
stm_inv_full = run_invasions( model = "stm", repro = repros_sp, nspp=nspp, 
			parms=parms_stm, time= ntime)

#Check the invasions with basic plots
par(mfrow = c(nspp,2))
for(n in 1:nspp){
	plot(ufm_inv_full[[n]][,n],t="l", ylim=c(0,1), main = "UFM")
	spp_com = 1:nspp
	spp_com = spp_com[-n]
	lines( rowSums(as.matrix(ufm_inv_full[[n]][,spp_com])),col="blue")
}

for(n in 1:nspp){
	plot(stm_inv_full[[n]][,n],t="l", ylim=c(0,1), main = "STM")
	spp_com = 1:nspp
	spp_com = spp_com[-n]
	lines( rowSums(as.matrix(ufm_inv_full[[n]][,spp_com])),col="blue")
}