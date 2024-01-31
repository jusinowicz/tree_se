Code to run two different models of tree competition/coexistence and measure the storage effect from numerical iterations. The two models are: 
Usinowicz forest model: Two-stage seedling, sapling and adult model from Usinowicz et al. (2012,2017)
Stump tree model: from Stump and Vasseur 2023. 

The code runs reciprocal invasions for all species against the community. Then measures the storage effect by comparing the LDGRs in the presence and absence of environmental variation. 

There are two core files at the moment: 
tree_se.R: The main file which generates random growth rates then goes through each step in the calculations. 
tree_se_functions.R: All of the functions where the hard work happens. Population models, basic simulations, simulation of invasion scenarios, calculation of LDGRs. 


