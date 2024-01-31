#=============================================================================
# Functions to run the storage effect code. This includes population models, 
# invasion growth rate scenarios, data parsing. 
#=============================================================================
#libraries: 

#=============================================================================
#Population models
#=============================================================================
#
#=============================================================================
# Two-stage seedling, sapling and adult model from Usinowicz et al. (2012,2017)
#=============================================================================
utm_step = function (si1, ni1, repro, aij, fi, del1_s ) { 

	#Run each species simultaneously
	csi1 = 1+sum(aij%*%(ni1[t,]*sl[t,]) ) #Sapling competition
	cni1 = sum(sl[t,]) #Adult competition 
	si1[t+1] = (si1[t,]*fi+ni1[t,]*sl[t,]/csi1[t,]) # Seedlings
	ni1[t+1,] = ni2[t,]*del1_s+(1-sum(ni1[t,]*del1_s)) *
				 si1[t,]/cni1[t,] #Adults


}