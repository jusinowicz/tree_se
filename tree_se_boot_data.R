#=============================================================================
# Main file for measuring the storage effect from real data. 
# 
# There are two different models at the moment:
# 
#=============================================================================
#load libraries 
source("./functions/tree_se_functions.R")


#=============================================================================
# Load the yearly counts.
# Remove species that do not meet basic sampling tresholds. 
# Scale the remaining species to standardize the mean and 
# variance effects. 
#=============================================================================
#Load data:
repros_sp = read.csv(file="./data/yearcount66.csv")

#Filter data according to e.g. total number of samples 
min.samp=100 #Used in Usinowicz et al. 2018
repros_sp.use=repros_sp [,which(colSums(repros_sp )>=min.samp)]

#Include this filter to make sure the number of YEARS sampled is also above 
#some threshold
yr.samp=3 #Used in Usinowicz et al. 2018
repros_sp.use =repros_sp.use [,which(colSums(repros_sp.use >0)>=yr.samp)]

#This scales each species' repro by first converting to logscale, 
#scaling to m=0, sd = 1, then transforming back with exp()
repros_sp_scale= exp( scale( log (repros_sp.use +1) ) )

#=============================================================================
# Global parameters
#=============================================================================
nspp = dim(repros_sp_scale)[2]
ntime = dim(repros_sp_scale)[1]

# Define the number of iterations and number of rows to subsample
n_bootstrap = 10
n_subsample = 200

# Initialize variables to store results
#Each species competes against each other species and itself. 
all_ufm_se = array( 0, dim = c(nspp, nspp, n_bootstrap))
all_stm_se = array( 0, dim = c(nspp, nspp, n_bootstrap))

# Get the covariance matrix from data
sig1 = cov(repros_sp_scale)

#=============================================================================
# Initialize variables for the models
#=============================================================================
#Model 1 adapted from Usinowicz et al. 2012
#parameters
parms_ufm = NULL
parms_ufm$aij = matrix(1,nspp,nspp) #Seedling competition
parms_ufm$fij = matrix(0.6,nspp,1) #Seedling survival
del1 = matrix(0.1,nspp,1)  #Adult death rate
parms_ufm$del1_s = 1 - del1 #Adult survival rate

#Model 2 adapted from Stump and Vasseur 2023
#parameters
parms_stm = NULL
parms_stm$del1 = matrix(0.1,nspp,1) #Adult death rate
parms_stm$att_e = matrix(1,nspp,1) #Enemy attack rate ????<====Any expectation from lit?   
# ==== How do we equate attack rates to aij? 

#=============================================================================
# Run simulations
#=============================================================================
# *************Implement this for species pairs to get visual checks
# ufm_sim = simulate_model ( model = "ufm", repro = repros_sp_scale, nspp=nspp, 
# 			parms=parms_ufm, time= ntime)
# stm_sim = simulate_model ( model = "stm", repro = repros_sp_scale, nspp=nspp, 
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
# Calculate each pairwise measurement of the storage effect. 
# Do this for multiple bootstrap iterations
#=============================================================================

for (bootstrap_iter in 1:n_bootstrap) {
	ntime = n_subsample
	# Perform bootstrap resampling on repros_sp_scale preserving correlations
	b_indices = sample(1:nrow(repros_sp_scale), size = n_subsample, replace = TRUE)
	repros_sub= repros_sp_scale [b_indices, ]
  
	#=============================================================================
	# Run the invasions to check coexistence. 
  	# This is currently implemented to check the invasion on a pairwise
  	# basis only! That is, it is not 1 spp vs community. 
	#=============================================================================

  	#Invader
  	for (s1 in 1:nspp) { 
  		#Cycle through residents
  		for(s2 in 1:nspp){

  			#Subset the reproduction and all parameters:
  			repros_sub_p = as.matrix(repros_sub[,c(s1,s2)])

  			#UFM
  			parms_ufm_sub = NULL
			parms_ufm_sub$aij = matrix(c(parms_ufm$aij[s1,s2], parms_ufm$aij[s1,s1],
										parms_ufm$aij[s2,s1], parms_ufm$aij[s2,s2]),
										2,2)
			parms_ufm_sub$fij = as.matrix(c(parms_ufm$fij[s1],parms_ufm$fij[s2]))
			parms_ufm_sub$del1_s = as.matrix(c(parms_ufm$del1_s[s1], 
												parms_ufm$del1_s[s2]))

			#STM
			parms_stm_sub = NULL
			parms_stm_sub$del1 = as.matrix(c(parms_stm$del1[s1],parms_stm$del1[s2])) #Adult death rate
			parms_stm_sub$att_e = as.matrix(c(parms_stm$att_e[s1],parms_stm$att_e[s2])) #Enemy attack rate ????<====Any expectation from lit?   

			#Run invasions 
			ufm_inv_full = run_invasions( model = "ufm", repro = repros_sub_p, nspp=2, 
						parms=parms_ufm_sub, time= ntime)
			stm_inv_full = run_invasions( model = "stm", repro = repros_sub_p, nspp=2, 
						parms=parms_stm_sub, time= ntime)

			# #Check the invasions with basic plots
			# par(mfrow = c(2,2))
			# for(n in 1:2){
			# 	plot(ufm_inv_full[[n]][,n],t="l",ylim=c(0,1), main = "UFM")
			# 	spp_com = 1:2
			# 	spp_com = spp_com[-n]
			# 	lines( rowSums(as.matrix(ufm_inv_full[[n]][,spp_com])),col="blue")
			# }

			# for(n in 1:2){
			# 	plot(stm_inv_full[[n]][,n],t="l", ylim=c(0,1), main = "STM")
			# 	spp_com = 1:2
			# 	spp_com = spp_com[-n]
			# 	lines( rowSums(as.matrix(ufm_inv_full[[n]][,spp_com])),col="blue")
			# }

			#=============================================================================
			# Measure the storage effect. Compare the average low-density growth rates 
			# across two scenarios: when environmental variation is present, and when 
			# the environmental response is held constant at its average.
			#=============================================================================
			#Run invasions with constant environment: 
			repros_con = matrix(colMeans(repros_sub_p), ntime, 2, byrow=T)
			ufm_inv_con= run_invasions( model = "ufm", repro = repros_con, nspp=2, 
						parms=parms_ufm_sub, time= ntime)
			stm_inv_con = run_invasions( model = "stm", repro = repros_con, nspp=2, 
						parms=parms_stm_sub, time= ntime)

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

			#Store the SE calculation for each boostrap iteration
			all_ufm_se[s1,s2,bootstrap_iter] = ufm_se$se
			all_stm_se[s1,s2,bootstrap_iter] = stm_se$se
		}
	}

}


#Get the average and sd across all the bootstraps for each pairwise interaction.
#Using apply is always faster than a for loop.
mean_ufm = apply(all_ufm_se, c(1, 2), function(x) mean(x)) 
sd_ufm = apply(all_ufm_se, c(1, 2), function(x) sd(x)) 
mean_stm = apply(all_stm_se, c(1, 2), function(x) mean(x)) 
sd_stm = apply(all_stm_se, c(1, 2), function(x) sd(x)) 

#=============================================================================
#Plot a histogram of the SE measurements.
#=============================================================================
mufm_plot = mean_ufm[upper.tri(mean_ufm)]
mstm_plot = mean_stm[upper.tri(mean_stm)]

par(mfrow = c(2,1))
hist(mufm_plot, main = "UFM model",xaxt="n",xlim=c(0,1))
# Add a vertical line for the overall median
abline(v = median(mufm_plot), col = "red", lwd = 2)
# Label the overall median on the x-axis
axis(1, at = c(0,median(mufm_plot), max(mufm_plot),1), labels = c(0,paste(round(median(mufm_plot),digits=3)), 
	paste(round(max(mufm_plot),digits=3)),1))

hist(mstm_plot, main = "STM model",xaxt="n",xlim=c(0,1))
# Add a vertical line for the overall median
abline(v = median(mstm_plot), col = "red", lwd = 2)
# Label the overall median on the x-axis
axis(1, at = c(0,median(mstm_plot), max(mstm_plot),1), labels = c(0,paste(round(median(mstm_plot),digits=3)), 
	paste(round(max(mstm_plot),digits=3)),1))