#=============================================================================
# Main file for measuring the storage effect from real data. 
#
# This version follows the approach in Usinowicz et al. 2012, 2017, which is 
# also originally presented conceptually in Warner and Chesson 1981. The key
# difference here is the approach to measuring the low-density growth rates
# (LDG): the LDG is measured as di (adult survival) goes to 1. This 
# standardizes the effect of di on LDGs and the storage effect. 
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
ntime = 100

# Define the number of iterations and number of rows to subsample
# n_subsample determines how long the pop sim will run. It is important
# to set this long enough so that it includes the full invasion, but 
# longer runs take more computation time. Unfortunately, I have not 
# come up with a good way to automate this so it requires some tinkering. 
n_bootstrap = 10
n_subsample = 1000

# Initialize variables to store results
#Each species competes against each other species and itself.
#Storage effect 
all_ufm_se = array( 0, dim = c(nspp, nspp, n_bootstrap))
all_stm_se = array( 0, dim = c(nspp, nspp, n_bootstrap))
sd_ufm_se = array( 0, dim = c(nspp, nspp, n_bootstrap))
sd_stm_se = array( 0, dim = c(nspp, nspp, n_bootstrap))

#The LDG (or LGR), low-density growth rate
all_ufm_lgr = array( 0, dim = c(nspp, nspp, n_bootstrap))
all_stm_lgr = array( 0, dim = c(nspp, nspp, n_bootstrap))
sd_ufm_lgr = array( 0, dim = c(nspp, nspp, n_bootstrap))
sd_stm_lgr = array( 0, dim = c(nspp, nspp, n_bootstrap))

#The LDG in constant environment
all_ufm_con = array( 0, dim = c(nspp, nspp, n_bootstrap))
all_stm_con = array( 0, dim = c(nspp, nspp, n_bootstrap))
#The Aij, competition coefficients
all_ufm_Aij = array( 0, dim = c(nspp, nspp, n_bootstrap))
all_stm_Aij = array( 0, dim = c(nspp, nspp, n_bootstrap))

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
#
# This is currently implemented to check the invasion on a pairwise
# basis only! That is, it is not 1 spp vs community. 
#=============================================================================

#Invader
for (s1 in 1:nspp) { 
  	#Cycle through residents
  	print(paste("Species 1 is ",s1,".",sep =""))

  	for(s2 in 1:nspp){

  	#UFM
  	parms_ufm_sub = NULL
	parms_ufm_sub$aij = matrix(c(parms_ufm$aij[s1,s1],parms_ufm$aij[s1,s2],
										parms_ufm$aij[s2,s1], parms_ufm$aij[s2,s2]),
										2,2)
	parms_ufm_sub$fij = as.matrix(c(parms_ufm$fij[s1],parms_ufm$fij[s2]))
	parms_ufm_sub$del1_s = as.matrix(c(parms_ufm$del1_s[s1], 
												parms_ufm$del1_s[s2]))

	#STM
	parms_stm_sub = NULL
	parms_stm_sub$del1 = as.matrix(c(parms_stm$del1[s1],parms_stm$del1[s2])) #Adult death rate
	parms_stm_sub$att_e = as.matrix(c(parms_stm$att_e[s1],parms_stm$att_e[s2])) #Enemy attack rate ????<====Any expectation from lit?   

	#=============================================================================
	# Run the invasions to check coexistence. 
	# Do this for multiple bootstrap iterations
	#=============================================================================
	ufm_inv_boot = vector("list", nspp)
	stm_inv_boot = vector("list", nspp)
	ufm_con_boot = vector("list", nspp)
	stm_con_boot = vector("list", nspp)

	for (bootstrap_iter in 1:n_bootstrap) {

		# Perform bootstrap resampling on repros_sp_scale preserving correlations
		b_indices = sample(1:nrow(repros_sp_scale), size = n_subsample, replace = TRUE)
		repros_sub= repros_sp_scale [b_indices, ]

		#Subset the reproduction and all parameters:
		repros_sub_p = as.matrix(repros_sub[,c(s1,s2)])

		#Run invasions 
		ufm_inv_full = run_invasions( model = "ufm", repro = repros_sub_p, nspp=2, 
					parms=parms_ufm_sub, time= ntime)
		stm_inv_full = run_invasions( model = "stm", repro = repros_sub_p, nspp=2, 
					parms=parms_stm_sub, time= ntime)

		for (n in 1:2){
			ufm_inv_boot[[n]] = rbind(ufm_inv_boot[[n]], unlist(c(ufm_inv_full[[n]][2,], 
										repros_sub_p[2,])) )
			stm_inv_boot[[n]] = rbind(stm_inv_boot[[n]], unlist(c(stm_inv_full[[n]][2,], 
							repros_sub_p[2,])) )
		}

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

		for (n in 1:2){
			ufm_con_boot[[n]] = rbind(ufm_con_boot[[n]], unlist(c(ufm_inv_con[[n]][2,], 
										repros_con[2,])) )
			stm_con_boot[[n]] = rbind(stm_con_boot[[n]], unlist(c(stm_inv_con[[n]][2,], 
							repros_con[2,])) )
		}


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

	}

	#Now do the calculations for LDG, SE, and Aij using the samples from 
	#each of the iterations of the model. 
	for( sp in 1:2){

			# all_ufm_lgr[s1,s2,1] = 
			# get_ldg_ufm ( xi = ufm_inv_boot[[sp]][,1:(2)],
			# 	si = ufm_inv_boot[[sp]][,(2+1):(2*2)],
			# 	ri = ufm_inv_boot[[sp]][,(3*2+1):(4*2)],
			# 	fi = parms_ufm_sub$fij,
			# 	di = parms_ufm_sub$del1_s,
			# 	aij = parms_ufm_sub$aij,
			# 	invader = sp
			# )

			#Approximate LDG for UFM as di -> 1
			ufm_lgr = 
					get_dij_ufm1( xi = ufm_inv_boot[[sp]][,1:(2)],
						si = ufm_inv_boot[[sp]][,(2+1):(2*2)],
						ri = ufm_inv_boot[[sp]][,(3*2+1):(4*2)],
						fi = parms_ufm_sub$fij,
						aij = parms_ufm_sub$aij,
						invader = sp,
						return = "full"
					)

			#Approximate LDG for UFM with constant environment as di -> 1
			ufm_con = 
					get_dij_ufm1( xi = ufm_con_boot[[sp]][,1:(2)],
						si = ufm_con_boot[[sp]][,(2+1):(2*2)],
						ri = ufm_con_boot[[sp]][,(3*2+1):(4*2)],
						fi = parms_ufm_sub$fij,
						aij = parms_ufm_sub$aij,
						invader = sp,
						return = "full"
					)

			#Take the averages
			all_ufm_lgr[s1,s2,1] = mean(ufm_lgr)
			sd_ufm_lgr[s1,s2,1]  = sd(ufm_lgr)
 
			all_ufm_con[s1,s2,1] = mean(ufm_con)

			#Storage effect is the average difference between the LDG with and 
			#without fluctuations (constant environment)
			all_ufm_se[s1,s2,1] = mean(ufm_lgr-ufm_con)
			sd_ufm_se[s1,s2,1] = sd(ufm_lgr-ufm_con)


			#The LDG can be easily converted to Aij now
			all_ufm_Aij[s1,s2,1]  = 1/exp(mean(ufm_lgr) )

			# all_stm_lgr[s1,s2,1] = 
			# 		get_ldg_stm( xi = stm_inv_boot[[sp]][,1:(2)],
			# 			ri = stm_inv_boot[[sp]][,(2+1):(2*2)],
			# 			del1 = parms_stm_sub$del1,
			# 			att_e = parms_stm_sub$att_e,
			# 			invader = sp
			# 		)

			#Approximate LDG for STM as di -> 1
			stm_lgr = 
					get_dij_stm1( xi = stm_inv_boot[[sp]][,1:(2)],
						ri = stm_inv_boot[[sp]][,(2+1):(2*2)],
						att_e = parms_stm_sub$att_e,
						invader = sp,
						return = "full"
					)

			#Approximate LDG for STM with constant environment as di -> 1
			stm_con = 
					get_dij_stm1( xi = stm_con_boot[[sp]][,1:(2)],
						ri = stm_con_boot[[sp]][,(2+1):(2*2)],
						att_e = parms_stm_sub$att_e,
						invader = sp,
						return = "full"
					)

			#Take the averages
			all_stm_lgr[s1,s2,1] = mean(stm_lgr)
			sd_stm_lgr[s1,s2,1]  = sd(stm_lgr)

			all_stm_con[s1,s2,1] = mean(stm_con)

			#Storage effect is the average difference between the LDG with and 
			#without fluctuations (constant environment)
			all_stm_se[s1,s2,1] = mean(stm_lgr - stm_con)
			sd_stm_se[s1,s2,1]  = sd(stm_lgr - stm_con)
			
			#The LDG can be easily converted to Aij now
			all_stm_Aij[s1,s2,1]  = 1/exp(mean(stm_lgr) )
		}

	}

}



#=============================================================================
#Plot a histogram of the SE measurements.
#=============================================================================
# maiju_plot = 1/(mean_ulgr *t(mean_ulgr))
# maijs_plot = 1/(mean_slgr *t(mean_slgr))
# maiju_plot = maiju_plot[upper.tri(maiju_plot)]
# maijs_plot = maijs_plot[upper.tri(maijs_plot)]
mean_ufm = all_ufm_lgr[,,1]
mean_stm = all_stm_lgr[,,1]

mean_ufm_se = all_ufm_se[,,1]
mean_stm_se = all_stm_se[,,1]

mean_ufm_aij = all_ufm_Aij[,,1]
mean_stm_aij = all_stm_Aij[,,1]



mufm_plot = mean_ufm[upper.tri(mean_ufm)]
mstm_plot = mean_stm[upper.tri(mean_stm)]

mufm_se_plot = mean_ufm_se[upper.tri(mean_ufm)]
mstm_se_plot = mean_stm_se[upper.tri(mean_stm)]

mufm_aij_plot = mean_ufm_aij[upper.tri(mean_ufm)]
mstm_aij_plot = mean_stm_aij[upper.tri(mean_stm)]

# mufm_plot = mean_ufm_se[upper.tri(mean_ufm_se)]
# mstm_plot = mean_stm_se[upper.tri(mean_stm_se)]

#Plot of the LGRs
pdf(file="logLGR_models1.pdf", height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)


par(mfrow = c(1,2))
hist(mufm_plot, main = "UFM model",xaxt="n")
# Add a vertical line for the overall median
abline(v = median(mufm_plot), col = "red", lwd = 2)
# Label the overall median on the x-axis
axis(1, at = c(0,median(mufm_plot), max(mufm_plot),1), labels = c(0,paste(round(median(mufm_plot),digits=3)), 
	paste(round(max(mufm_plot),digits=3)),1))

hist(mstm_plot, main = "STM model",xaxt="n")
# Add a vertical line for the overall median
abline(v = median(mstm_plot), col = "red", lwd = 2)
# Label the overall median on the x-axis
axis(1, at = c(0,median(mstm_plot), max(mstm_plot),1), labels = c(0,paste(round(median(mstm_plot),digits=3)), 
	paste(round(max(mstm_plot),digits=3)),1))

dev.off()


#Plot of the Storage Effect
pdf(file="SE_models1.pdf", height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow = c(1,2))
hist(mufm_se_plot, main = "UFM model",xaxt="n")
# Add a vertical line for the overall median
abline(v = median(mufm_se_plot), col = "red", lwd = 2)
# Label the overall median on the x-axis
axis(1, at = c(0,median(mufm_se_plot), max(mufm_se_plot),1), labels = c(0,paste(round(median(mufm_se_plot),digits=3)), 
	paste(round(max(mufm_se_plot),digits=3)),1))

hist(mstm_se_plot, main = "STM model",xaxt="n")
# Add a vertical line for the overall median
abline(v = median(mstm_se_plot), col = "red", lwd = 2)
# Label the overall median on the x-axis
axis(1, at = c(0,median(mstm_se_plot), max(mstm_se_plot),1), labels = c(0,paste(round(median(mstm_se_plot),digits=3)), 
	paste(round(max(mstm_se_plot),digits=3)),1))

dev.off()


#Plot of the Aij (one-sided)
pdf(file="aij_models1.pdf", height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow = c(1,2))
hist(mufm_aij_plot, main = "UFM model",xaxt="n")
# Add a vertical line for the overall median
abline(v = median(mufm_aij_plot), col = "red", lwd = 2)
# Label the overall median on the x-axis
axis(1, at = c(0,median(mufm_aij_plot), max(mufm_aij_plot),1), labels = c(0,paste(round(median(mufm_aij_plot),digits=3)), 
	paste(round(max(mufm_aij_plot),digits=3)),1))

hist(mstm_aij_plot, main = "STM model",xaxt="n")
# Add a vertical line for the overall median
abline(v = median(mstm_aij_plot), col = "red", lwd = 2)
# Label the overall median on the x-axis
axis(1, at = c(0,median(mstm_aij_plot), max(mstm_aij_plot),1), labels = c(0,paste(round(median(mstm_aij_plot),digits=3)), 
	paste(round(max(mstm_aij_plot),digits=3)),1))

dev.off()


# hist(maiju_plot, main = "UFM model",xaxt="n",xlim=c(0,1))
# # Add a vertical line for the overall median
# abline(v = median(maiju_plot), col = "red", lwd = 2)
# # Label the overall median on the x-axis
# axis(1, at = c(0,median(maiju_plot), max(maiju_plot),1), labels = c(0,paste(round(median(maiju_plot),digits=3)), 
# 	paste(round(max(maiju_plot),digits=3)),1))

# hist(maijs_plot, main = "STM model",xaxt="n",xlim=c(0,1))
# # Add a vertical line for the overall median
# abline(v = median(maijs_plot), col = "red", lwd = 2)
# # Label the overall median on the x-axis
# axis(1, at = c(0,median(maijs_plot), max(maijs_plot),1), labels = c(0,paste(round(median(maijs_plot),digits=3)), 
# 	paste(round(max(maijs_plot),digits=3)),1))

save(file = "test_se_r1.var",all_ufm_se,all_stm_se,all_ufm_lgr,all_stm_lgr,all_ufm_con,all_stm_con )