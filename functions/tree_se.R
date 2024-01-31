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
#repros_sp = random_repro (nspp,ntime)

#Specify means and a covariance matrix: 
# Means:
mus = c(1.01,1.1)
# Random correlation matrix
cm1 = matrix(c(1,-0.1,-0.1,1), ncol = nspp)
cm1 =  cm1 %*% t(cm1)

# Standard deviations of the time series
sd1 = rnorm(nspp,mean=1,sd=0.1)

# Create the covariance matrix
sig1 = diag(sd1) %*% cm1 %*% diag(sd1)

repros_sp = random_repro (nspp,ntime, cov = sig1)
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
# Run the invasions to check coexistence. 
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

#=============================================================================
# Measure the storage effect. Compare the average low-density growth rates 
# across two scenarios: when environmental variation is present, and when 
# the environmental response is held constant at its average.
#=============================================================================
#Run invasions with constant environment: 
repros_con = matrix(colMeans(repros_sp),ntime, nspp, byrow=T)
ufm_inv_con= run_invasions( model = "ufm", repro = repros_con, nspp=nspp, 
			parms=parms_ufm, time= ntime)
stm_inv_con = run_invasions( model = "stm", repro = repros_con, nspp=nspp, 
			parms=parms_stm, time= ntime)

#Check the invasions with basic plots
par(mfrow = c(nspp,2))
for(n in 1:nspp){
	plot(ufm_inv_con[[n]][,n],t="l", ylim=c(0,1), main = "UFM")
	spp_com = 1:nspp
	spp_com = spp_com[-n]
	lines( rowSums(as.matrix(ufm_inv_con[[n]][,spp_com])),col="blue")
}

for(n in 1:nspp){
	plot(stm_inv_con[[n]][,n],t="l", ylim=c(0,1), main = "STM")
	spp_com = 1:nspp
	spp_com = spp_com[-n]
	lines( rowSums(as.matrix(ufm_inv_con[[n]][,spp_com])),col="blue")
}

#Get the IGRs as the slope of log-transformed growth rate while it is below 
#a (somewhat arbitrary) low density

ufm_se = vector("list", nspp) 
stm_se = vector("list", nspp)
s_ind = 1:nspp

for(s in 1:nspp){ 
	#UFM
	#Get the growth rates for invader and resident community
	#with variation. 
	ufm_se[[n]]$full = get_ldgrs(ufm_inv_full[[s]], invader = s)
	#Get the growth rates for invader and resident community
	#no variation. 
	ufm_se[[n]]$con = get_ldgrs(ufm_inv_con[[s]], invader = s)
	#The storage effect is the difference between these growth rates
	ufm_se[[n]]$se = ufm_se[[n]]$full$inv - ufm_se[[n]]$full$res

	#STM
	#Get the growth rates for invader and resident community
	#with variation. 
	stm_se[[n]]$full = get_ldgrs(stm_inv_full[[s]], invader = s)
	#Get the growth rates for invader and resident community
	#no variation. 
	stm_se[[n]]$con = get_ldgrs(stm_inv_con[[s]], invader = s)
	#The storage effect is the difference between these growth rates
	stm_se[[n]]$se = stm_se[[n]]$full$inv - stm_se[[n]]$full$res

}



