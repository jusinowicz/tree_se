################################## First, define functions 
####### This function is used below to compute moving sums over a pre-specified window, k
winsum = function(x, k) { 
	n = length(x) 
	# this is a difference from the previous cell 
	y = x[ k:n ] - x[ c(1,1:(n-k)) ] 
	# find the first sum 
	y[1] = sum(x[1:k]); 
	# apply precomputed differences 
	y = cumsum(y) 
	return(y)
}

##############################################################################
#Code for parsing weekly seed counts starts here
##############################################################################

wanted=read.csv("spp_key.csv", header=FALSE, stringsAsFactors=FALSE)
seedXY = read.csv('BCI_TRAP200_20111216.csv', header=TRUE)
seedconv=read.csv("seed_conversion.csv", header=TRUE)
#seedXY$DATE.D.1=as.Date(seedXY$DATE.D.1)
seedXY$MONTH=as.numeric(format(as.Date(seedXY$DATE.D.1), "%m"))
seedXY$YEAR=as.numeric(format(as.Date(seedXY$DATE.D.1), "%Y"))

###This first part pulls out the specified species in "wanted" 

#This may need to be changed (i.e. between V2 and V3) based on the format of the original csv file above
speciescom=wanted$V2

tags=speciescom
ntags=length(tags)

seed.rec= matrix ( data=0, 1, ncol(seedXY))
#This line will set variable names which may need to be adjusted below
colnames(seed.rec)=colnames(seedXY)
index=1;
for (i in 1:ntags){
		#Variables refer to columns that may be file-specific
		seed.temp=subset(seedXY, seedXY$SPECIES.C4 == tags[i])
		#Include only instances of mature fruit (part = 1) or diaspore (part =2)
		seed.temp=subset(seed.temp,seed.temp$PART.N20==1 | seed.temp$PART.N20==2)
		plots=unique(seed.temp$TRAP.N40)
		nplots=length(plots)
	for (j in 1:nplots){
			
			seed.plot.temp=subset(seed.temp, seed.temp$TRAP.N40 == plots[j])
			seed.rec=rbind(seed.rec, seed.plot.temp)
			}
		}
#Variable refers to name of column that may be file-specific
seed.rec=subset(seed.rec, seed.rec$YEAR != 0)

#To export this
write.table(seed.rec, "seed_rain108.csv", row.names=FALSE, col.names=TRUE)

###############################################################################
## Create monthly counts per species from weekly counts. 
###############################################################################


nspp=length(speciescom)
months=seq(1,12,1)
nmnths=length(months)
years=unique(seed.rec$YEAR)
years=sort(years)
nyrs=length(years)

month.count = matrix ( data=0, nyrs*nmnths, nspp)

for (i in 1:nspp){
	#Variables refer to columns that may be file-specific
	conv.tmp=subset(seedconv, seedconv$SP4 == speciescom[i])[4]
	if (nrow(conv.tmp) ==0){conv.tmp=matrix(0,1,1)}
 	if(is.na(conv.tmp)==1){conv.tmp=matrix(0,1,1)}
	for (j in 1:nyrs) {
		annual.tmp=subset(seed.rec, seed.rec$SPECIES.C4 ==speciescom[i] & seed.rec$YEAR==years[j])
		for(m in 1:nmnths){
		mnth.tmp=subset(annual.tmp, annual.tmp$MONTH==m)
	
			#Pick out seeds vs. mature fruit and perform conversion in the case of fruit.
			conv.mat=as.matrix(as.matrix(conv.tmp)*as.numeric(mnth.tmp$PART.N20==1))
			seed.mat=as.matrix(as.numeric(mnth.tmp$PART.N20==2))
			month.count[m+(j-1)*12,i] = t(seed.mat)%*%as.matrix(mnth.tmp$QUANTITY.N61)+t(conv.mat)%*%as.matrix(mnth.tmp$QUANTITY.N61)

		}
				
	}
}

colnames(month.count)=speciescom

###<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
###<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#########################################################################################
# Turn monthly counts into annual counts. 
# It deals with seed production in overlapping calendar 
# years by (1) identifying peaks, then (2) placing all seed production associated with that peak (to the left and right)
# into the year of that peak. Years may contain an arbitrary number of peaks. It does assume that peaks are not wider
# than a single year. 
#########################################################################################
#Load the monthly data:
month.count = read.csv(file="monthcount108.csv")
speciescom = colnames(month.count)

nspp=length(speciescom)
months=seq(1,12,1)
nmnths=length(months)
nyrs = floor(dim(month.count)[1]/12) #years=unique(seed.rec$YEAR)
years=1:nyrs
#nyrs=length(years)

year.count = matrix ( data=0, nyrs, nspp)
k1=11

for (i in 1:nspp){

	sp.cnt=month.count[,i]

	# Use these as initial counts. If fruiting 
	# events are about annual, this should be fine. 
	# Otherwise, code below corrects for overlap in peak 
	# events (i.e. events that happen closer than 11 months apart) 
	sp.sum=winsum(sp.cnt, k1)
	
	#Use rle and diff to pick out the plateaus and match them 
	#with correct year (should be last year of the repeated sequence)
	sp.sum.rle=rle(sp.sum)
	
	#Find the peaks
	sp.sum.pks=as.numeric(diff(sign(diff(c(0,sp.sum.rle$values,0)))) == -2)
	
	#Identify the right edge of each plateau 
	sp.sum.month=cumsum(sp.sum.rle$lengths)
	
	#Put each edge into a year
	yin=ceiling(sp.sum.month/12)
	y.index=sp.sum.pks*yin
	sum.index=sp.sum.pks*sp.sum.rle$values
	for (y in 1:length(y.index)){year.count[y.index[y],i]=sum.index[y]}

	# Now the error checking, based on finding peaks that 
	# violate the assumption of being space more than 11 months. 
	# This takes second derivative as before to identify changes 
	# in slope, then picks out the index of each peak.
	# Most things end up here. 

	sp.cnt.rle=rle(sp.cnt)
	sp.peaks=which(rep(diff(sign(diff(c(0,sp.cnt.rle$values,0)))) == -2, 
			times= sp.cnt.rle$lengths))

	#Distance between peaks
	pks.dist=diff(sp.peaks)

	#Are peaks closer than a year? 
	less12=c(as.numeric(pks.dist<=(k1)),0)

	if (sum(less12)>0) {
	
		#Label the year that each peak occurs in.
		less12.off=c(0, less12[1:(length(less12)-1)]) 
		pks12=sp.peaks*as.numeric( (less12.off+less12)>0)
		pks12=pks12[pks12>0]
		yrs.pk = ceiling(pks12/12)
		#Reset problem years
		year.count[yrs.pk,i]=0

		for (y in 1:length(pks12) ){
			#Trace out peaks, add up seeds, and put them in 
			#appropriate year based on location of the peak
			
			if (yrs.pk[y]==1){
				#Pad with zeros to avoid accessing non-existent elements at beginning
				sp.cntB=c(matrix(0,5,1), sp.cnt)
				tmp.pk=sp.cntB[(pks12[y]):(pks12[y]+10)]
			
			}else if (yrs.pk[y]==max(nyrs)){
				#Pad with zeros to avoid accessing non-existent elements at the end
				sp.cntE=c(sp.cnt,matrix(0,5,1))
				tmp.pk=sp.cntE[(pks12[y]-5):(pks12[y]+5)]

			}else {  
				#Center peak with 5 months either side
				tmp.pk=sp.cnt[(pks12[y]-5):(pks12[y]+5)]}
		
			
			#Use to find edges of peaks, which could be a change 
			#in direction (-2) or just a flattening out (-1)
			tmp.pk.rle=rle(-tmp.pk)	
			find.pk=as.numeric(rep(diff(sign(diff(c(0,tmp.pk.rle$values,0)))) < 0, 
				times= tmp.pk.rle$lengths))	

			#By starting ind.left= 0 and ind.right = 1, when peaks 
			#share a non-zero minima it will go to the earlier 
			#(left) peak. 
			lft.start=5
			ind.left=0
			chk.pk=find.pk[lft.start]

			#cycle left until the left edge is found. Use ind.left to keep track of distance(months)
			while (chk.pk != 1){
				ind.left=ind.left+1
				lft.start=lft.start-1
				chk.pk=find.pk[lft.start] 
				if (length(chk.pk)<1) { chk.pk=1; ind.left=0}	
				if (is.na(chk.pk)) { chk.pk=1}
			}

			#Now the other side.
			right.start=7
			ind.right=1
			chk.pk=find.pk[right.start]
			
			#cycle right until the bottom of the first peak is found. 
			#Use ind.right to keep track of distance(months)
			while (chk.pk != 1){
				ind.right=ind.right+1
				right.start=right.start+1
				chk.pk=find.pk[right.start] 
				if (length(chk.pk)<1){ chk.pk=1; ind.right=0}
				if (is.na(chk.pk)) { chk.pk=1}
			}
			if ( is.na(sp.cnt[pks12[y]+ind.right])){ind.right=0}
			#Add this to the year the peak resides in
			year.count[yrs.pk[y],i] = year.count[yrs.pk[y],i] + 
					sum (sp.cnt[(pks12[y]-ind.left):(pks12[y]+ind.right)])
			
		}
    }
}

colnames(year.count)=speciescom

#Be careful about overwriting data: 
# write.table(year.count, file="yearcount108.csv", col.names=TRUE,sep=",")
# write.table(month.count, file="monthcount108.csv", col.names=TRUE,sep=",")


