###############################################################################
# The code below calculates AijAji, which quantify the potential of recruitment 
# fluctuations to promote coexistence via the storage effect. See Usinowicz 
# et al. 2012 and Usinowicz et al. 2017. 
#  
# These calculations are designed for a simple CSV or data file with annual 
# totals of reproduction (e.g. seed production or seedling recruitment), 
# where: 
# 	columns = species ID
#	rows = year
# 
# The example code is given with a fake example data set (see below for file name).
# Not included here is the code to calculate annual counts from raw data. Although
# I am happy to furnish these upon request. 
#
# The calculation of Aij is based on finding analytical invasion growth rates
# for a particular model of competition. The model has two stages, 
# where the key component supporting the storage effect 
# is lottery competition among saplings in a persistent
# sapling bank to recruit into an adult stage. However, population dynamics 
# are not simulated here and the details of the model can be found in the above
# references. 
#
# In this version of the code, reproduction is just the actual data set.
#
#=============================================================================
# Data input and filtering
#=============================================================================
#This file should be organized with species as row, and years as columns organized 
#oldest=1 to most current=nrow

filenames = list.files(path="./", pattern='*.csv')

all_aij = list()
all_aijlist=list()

for (a in 1:length(filenames)){

	print(paste("site", a))
	#Annual seedling counts for full community
	sl.full=read.csv(paste(filenames[a]), header=T)[,-1]

	#Filter data according to e.g. total number of samples 
	min.samp=25
	sl.use=sl.full[which(rowSums(sl.full)>=min.samp),]

	#Include this filter to make sure the number of YEARS sampled is also above some threshold
	yr.samp=3
	sl.use=sl.use[which(rowSums(sl.use>0)>=yr.samp), ]

	#=============================================================================
	# OR, if data are ready for calculations 
	#=============================================================================

	#sl=read.table(file="seedlings_share_pasoh.csv", sep=",")

	#=============================================================================
	# Main body -- Calculation of Aij for all pairwise combinations
	#=============================================================================
	#sl=t(sl.use)
	sl=sl.use

	#Number of spp, number of years
	nsp=nrow(sl)
	nt=ncol(sl)
	sl.dat=as.matrix(sl)

	#Number of simulations/generations
	ngens=10^3

	#Sapling survival rate
	f=0.6

	#Length of tau, lag over which to sum for approximation
	l.tau=15

	#f^tau for the sequence length l.tay
	f.tau=f^seq(from=0, to=l.tau, by=1)

	#Pre-allocate the matrix of competitive interactions
	simAij=matrix(0,nrow=nsp,ncol=nsp)
	#Apparent competition coefficients (lowercase alpha). These are all set to 1. 
	alpha.ij=matrix(1,2,2)

		for (i in 1:nsp) {
			print(paste("species i", i))

			for (j in 1:nsp){

				switchd=0
				#Don't waste time on species interaction with itself
				if (i !=j){ 

						# Generate sequences of reproduction. We pick out a pair of species, 
						# then simulate an arbitrary number of reproduction events. 
						# This is the jackkifed version. Pairs of events are chosen from the 
						# original distribution to preserve correlations.  

						# Pick out species i and j. Simulate distribution using bootstrapping: 
						# sample pairs from the data set in random order
						sl.temp=matrix(c(sl.dat[i,]+1,sl.dat[j,]+1),ncol=2)		
						samp=1+round((nt-1)*runif(ngens))
						rsim=sl.temp[samp,]
					 
						#Scale fluctuations relative to the mean of each species
						rsiml=log(rsim)
						rsim=exp((rsiml-t(matrix(colMeans(rsiml),2,ngens)))%*%diag(diag(1/sqrt(var(rsiml)))))


						######################
						# Calculate Aij. 
						######################
			
						#This for loop performs the sums in equation 4. 
						sumS=0; #matrix(0,(ngens-l.tau))
						x2s = matrix(0,ngens,2)
						for (t in 2:(ngens-l.tau)){ 
								x2=rsim[t:(t+l.tau),2]/(1+alpha.ij[2,1]*rsim[t:(t+l.tau),1]);
								s2=f.tau%*%x2;
								x2s[t,2] = x2[1]

								x1=rsim[t:(t+l.tau),1]/(1+alpha.ij[1,1]*rsim[t:(t+l.tau),1]);
								s1=f.tau%*%x1;
				               	x2s[t,1] = x1[1]
								sumS=sumS+s2/s1;

								#if(is.infinite(sumS) || is.nan(sumS)){ sumS=0}
						}

					#The expectation of the sum
					EsumS=sumS/(ngens-l.tau)
					Dij=EsumS-1
					#Convert Dij to Aij
					simAij[i,j]=1/(1+Dij)




				} #End if
			} #End j, resident species

		} #End i, invader species

		
	
	#Once Aij is calculated for every pairwise combination, calculate the AijAji for every pair. 
	simAijAji=simAij*t(simAij)
	#For plotting
	simAijAjilist=simAijAji
	sim.lt=simAijAjilist[lower.tri(simAijAjilist)]
	simAijAjilist=as.matrix(sim.lt)

	all_aij[[a]]=simAij
	all_aijlist[[a]] = simAijAjilist
	

}


#The correct order of the sites: 
#For doubled sites:
#co = c(27, 25, 11,21, 19,17,15, 23, 13 )
#co2 = co+1

co = c(14, 13, 6,11, 10,9,8, 12, 7 )

#For Bonanza:
a1.tmp=NULL
#for (r in seq(1,9,2)){ a1.tmp=rbind (a1.tmp, unlist(all_aijlist[[r]]))}
for (r in seq(1,5,1)){ a1.tmp=rbind (a1.tmp, unlist(all_aijlist[[r]]))}
var(c(a1.tmp))

# a2.tmp=NULL
# for (r in seq(2,10,2)){ a2.tmp=rbind (a2.tmp, unlist(all_aijlist[[r]]))}
# median(a2.tmp)
# var(c(a2.tmp))

cs1 = matrix(0,10,1)
csvar1 =matrix(0,10,1)

#cs2 = matrix(0,10,1)
#csvar2 =matrix(0,10,1)
for (r in 1:9) { 
	
	aij.tmp=c(unlist(all_aijlist[[co[r]]]))
	aij.tmp = aij.tmp[aij.tmp>0]
	cs1[r] = median(aij.tmp,na.rm=T)

	# aij2.tmp=c(unlist(all_aijlist[[co2[r]]]))
	# aij2.tmp = aij2.tmp[aij2.tmp>0]
	# cs2[r] = median(aij2.tmp,na.rm=T)
	
	csvar1[r] = var(aij.tmp,na.rm=T) 
	#csvar2[r] = var(aij2.tmp,na.rm=T)
}

cs1[10] = median(c(a1.tmp,na.rm=T))
csvar1[10] = var(c(a1.tmp,na.rm=T))

#cs2[10] = median(c(a2.tmp,na.rm=T))
#csvar2[10] = var(c(a2.tmp,na.rm=T))

lat1=c(-0.68, 2.97, 9.15, 18.32, 22.059, 24.76, 35, 36.93,42.38, 64.86)
lm1=lm(cs1~lat1)
summary(lm1)


#save this as a file	
save(all_aij, all_aijlist, file="simaijall_2.var")


######################
#Plotting
######################

#The histograms of AijAji, for simulated data
# par(mfrow = c(ceiling(length(filenames))/4,4 ))
# for (a in 1:length(filenames)){

# 	hist(unlist(all_aijlist[[a]]),xlim=c(0,1.2), breaks=40, main=print(filenames[a]))

# }

pdf(file="simaijall_null.pdf", height=5.5, width=5.5, family='Helvetica', pointsize=12)
plot(lat1,cs1,ylim=c(0,1), ylab = "median A_ijA_ji",xlab="latitude")
points(lat1,cs2,col="red")
dev.off()

pdf(file="simaijall_real_jk_leverage.pdf", height=5.5, width=5.5, family='Helvetica', pointsize=12)
plot(lm1,which=5)
dev.off()