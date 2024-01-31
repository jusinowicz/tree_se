library(MASS)

ngens = 5000
nspp = 2

#Randomly generate species repro
mu1 = abs(rnorm(nspp,1,0.1 ) ) #Mean
#Create a random covariance matrix: 

# Random correlation matrix
cm1 = matrix(runif(nspp^2), ncol = nspp)
cm1 =  cm1 %*% t(cm1)
diag(cm1) = 1  # Diagonal elements should be 1

# Standard deviations of the time series
sd1 = runif(nspp)

# Create the covariance matrix
sig1 = diag(sd1) %*% cm1 %*% diag(sd1)

#use mvrnorm to make time series (lognormal)
sl = exp(mvrnorm(ngens, mu1, sig1))

#Model 1 adapted from Usinowicz et al. 2012
si1 = matrix(0.1, ngens, nspp) #Saplings
ni1 = matrix(0.1, ngens, nspp) #Adults
csi1 = matrix(0, ngens, nspp) #Sapling competition
cni1 = matrix(0, ngens, nspp) #Adult competition

#parameters
aij = matrix(1,nspp,nspp) #Seedling competition
fij = matrix(0.6,nspp,1) #Seedling survival
del1 = c(0.1,0.1) #Adult death rate
del1_s = 1 - del1 #Adult survival rate

#Simulate over all time
for(t in 1:(ngens-1){
	#Run each species simultaneously
	csi1[t,] = 1+sum(aij%*%(ni1[t,]*sl[t,]) ) #Sapling competition
	cni1[t,] = sum(sl[t,]) #Adult competition 
	si1[t+1] = (si1[t,]*fi+ni1[t,]*sl[t,]/csi1[t,]) # Seedlings
	ni1[t+1,] = ni2[t,]*del1_s+(1-sum(ni1[t,]*del1_s)) *
				 si1[t,]/cni1[t,] #Adults
}

#Model 2 adapted from Stump and Vasseur 2023
ni2 = matrix(0.1, ngens, nspp) #Adults 
pi2 = matrix(0, ngens, nspp) #Seed death
ci2 = matrix(0, ngens, nspp) #Competition

#parameters: 
att_e = c(1,1) #Enemy attack rate ????<====Any expectation from lit?   

#Simulate over all time
for(t in 1:(ngens-1) ){
	#Run each species simultaneously
	pi2[t,] = exp ( - att_e *ni2[t,] )
	ci2[t,] = sum(ni2[t,] * sl[t,] * pi2[t,] )
	ni2[t+1,] = ni2[t,]*(1-del1+del1*(sl[t,]*pi2[t,]/ci2[t,] ) ) 
}
