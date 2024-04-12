#=============================================================================
# Main file for measuring the storage effect from real data. 
# 
# There are two different models at the moment:
# 
#=============================================================================
#load libraries 
source("./functions/tree_se_functions.R")


#=============================================================================
# Load the yearly counts
#=============================================================================
#Load data:
repros_sp = read.csv(file="./data/yearcount66.csv")

#=============================================================================
# Global parameters
#=============================================================================
nspp = dim(repros_sp)[2]
ntime = dim(repros_sp)[1]

# Define the number of iterations and number of rows to subsample
n_bootstrap = 10
n_subsample = 200

# Initialize variables to store results
#Each species competes against each other species and itself. 
all_ufm_se = array( 0, dim = c(nspp, nspp, n_bootstrap))
all_stm_se = array( 0, dim = c(nspp, nspp, n_bootstrap))

# Get the covariance matrix from data
sig1 = cov(repros_sp)

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
# ==== How do we equate attack rates to aij? 

#=============================================================================
# Run simulations
#=============================================================================
# *************Implement this for species pairs to get visual checks
# ufm_sim = simulate_model ( model = "ufm", repro = repros_sp, nspp=nspp, 
# 			parms=parms_ufm, time= ntime)
# stm_sim = simulate_model ( model = "stm", repro = repros_sp, nspp=nspp, 
# 			parms=parms_stm, time= ntime)

# #Check simulations with basic plots
# par(mfrow = c(2,1))
# plot(ufm_sim[,1],t="l", ylim=c(0,1), main = "UFM")
# for(n in 2:nspp){
# 	lines( ufm_sim[,n],col="blue")
# }


# plot(stm_sim[,1],t="l", ylim=c(0,1), main = "STM")
# for(n in 2:nspp){
# 	lines( stm_sim[,n],col="blue")
# }

#=============================================================================
# Loop through each bootstrap iteration
#=============================================================================

for (bootstrap_iter in 1:n_bootstrap) {
	ntime = n_subsample
  # Perform bootstrap resampling on repros_sp preserving correlations
  b_indices = sample(1:nrow(repros_sp), size = n_subsample, replace = TRUE)
  repros_sub= repros_sp [b_indices, ]
  
	#=============================================================================
	# Run the invasions to check coexistence. 
  	# This is currently implemented to check the invasion on a pairwise
  	# basis only! That is, it is not 1 spp vs community. 
	#=============================================================================

  	#Invader
  	for (s1 in 1:nspp) { 
  		#Cycle through residents
  		for(s2 in 1:nspp){
  			repros_sub_p = repros_sub[,c(s1,s2)]


			ufm_inv_full = run_invasions( model = "ufm", repro = repros_sub_p, nspp=nspp, 
						parms=parms_ufm, time= ntime)
			stm_inv_full = run_invasions( model = "stm", repro = repros_sub_p, nspp=nspp, 
						parms=parms_stm, time= ntime)

			# #Check the invasions with basic plots
			# par(mfrow = c(nspp,2))
			# for(n in 1:nspp){
			# 	plot(ufm_inv_full[[n]][,n],t="l", ylim=c(0,1), main = "UFM")
			# 	spp_com = 1:nspp
			# 	spp_com = spp_com[-n]
			# 	lines( rowSums(as.matrix(ufm_inv_full[[n]][,spp_com])),col="blue")
			# }

			# for(n in 1:nspp){
			# 	plot(stm_inv_full[[n]][,n],t="l", ylim=c(0,1), main = "STM")
			# 	spp_com = 1:nspp
			# 	spp_com = spp_com[-n]
			# 	lines( rowSums(as.matrix(ufm_inv_full[[n]][,spp_com])),col="blue")
			# }

			#=============================================================================
			# Measure the storage effect. Compare the average low-density growth rates 
			# across two scenarios: when environmental variation is present, and when 
			# the environmental response is held constant at its average.
			#=============================================================================
			#Run invasions with constant environment: 
			repros_con = matrix(colMeans(repros_sub_p), ntime, nspp, byrow=T)
			ufm_inv_con= run_invasions( model = "ufm", repro = repros_con, nspp=nspp, 
						parms=parms_ufm, time= ntime)
			stm_inv_con = run_invasions( model = "stm", repro = repros_con, nspp=nspp, 
						parms=parms_stm, time= ntime)

			# #Check the invasions with basic plots
			# par(mfrow = c(nspp,2))
			# for(n in 1:nspp){
			# 	plot(ufm_inv_con[[n]][,n],t="l", ylim=c(0,1), main = "UFM")
			# 	spp_com = 1:nspp
			# 	spp_com = spp_com[-n]
			# 	lines( rowSums(as.matrix(ufm_inv_con[[n]][,spp_com])),col="blue")
			# }

			# for(n in 1:nspp){
			# 	plot(stm_inv_con[[n]][,n],t="l", ylim=c(0,1), main = "STM")
			# 	spp_com = 1:nspp
			# 	spp_com = spp_com[-n]
			# 	lines( rowSums(as.matrix(ufm_inv_con[[n]][,spp_com])),col="blue")
			# }

			#Get the IGRs as the slope of log-transformed growth rate while it is below 
			#a (somewhat arbitrary) low density

			ufm_se = NULL 
			stm_se = NULL

			#UFM
			#Get the growth rates for invader and resident community
			#with variation. 
			ufm_se$full = get_ldgrs(ufm_inv_full[[1]], invader = 1)
			#Get the growth rates for invader and resident community
			#no variation. 
			ufm_se$con = get_ldgrs(ufm_inv_con[[1]], invader = 1)
			#The storage effect is the difference between these growth rates
			ufm_se$se = ufm_se$full$inv - ufm_se$full$res

			#STM
			#Get the growth rates for invader and resident community
			#with variation. 
			stm_se$full = get_ldgrs(stm_inv_full[[1]], invader = 1)
			#Get the growth rates for invader and resident community
			#no variation. 
			stm_se$con = get_ldgrs(stm_inv_con[[1]], invader = 1)
			#The storage effect is the difference between these growth rates
			stm_se$se = stm_se$full$inv - stm_se$full$res

			all_ufm_se[s1,s2,bootstrap_iter] = ufm_se$se
			all_stm_se[s1,s2,bootstrap_iter] = stm_se$se
		}
	}

}



mean_ufm = colMeans(all_ufm_se)
sd_ufm = sqrt(apply(all_ufm_se, 2, var))
mean_stm = colMeans(all_stm_se)
sd_stm = sqrt(apply(all_stm_se, 2, var))

# mean_ufm = sapply(all_ufm_se, function(ufm_se_iter) sapply(ufm_se_iter, function(se) mean(unlist(se))))
# sd_ufm = sapply(all_ufm_se, function(ufm_se_iter) sapply(ufm_se_iter, function(se) sd(unlist(se))))
# mean_stm = sapply(all_stm_se, function(stm_se_iter) sapply(stm_se_iter, function(se) mean(unlist(se))))
# sd_stm = sapply(all_stm_se, function(stm_se_iter) sapply(stm_se_iter, function(se) sd(unlist(se))))
