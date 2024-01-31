#R code to simulate species seed production at monthly intervals and then 
#count the seeds and bin them into a year based on the presence of the peak. 
#
#Simulation happens in "get_sp"
#Peak identification happens in "get_peaks"
#Summation and binning happens in "get_seed_sums"
#
#The simulation only produces unimodal distributions, but this should work 
#regardless of the number of peaks. It will just always center the summation
#based on the largest peak. 

###Define functions
#Function to simulate a species-specific seasonal pattern of peaks
get_sp = function(nyears, mseeds = 100, a_rr = 2 ){

  #initialize matrix to store data
  sp_data = matrix(0,nyears*12,1)

  #Place the total seeds at each peak month
  peak_month = sample (1:12,1) #Which month is peak? 
  place_index = c(peak_month+(0:(nyears-1))*12)
  sp_data[place_index] = rnorm(length(place_index), mseeds, (mseeds/10) )

  # Define a range for the dispersal kernel
  lims = floor((nyears*12)/2)
  x_range = seq(-lims+1, lims, by =1)

  # Calculate the kernel values (this is an exponential kernel)
  # The only parrameter is a_rr, which controls the spread
  # (standard deviation). 
  kernel_values = a_rr/2*exp(-a_rr*abs(x_range))
  kernel_values = kernel_values/sum(kernel_values)

  # FFT of the data and kernel
  fft_data = fft(sp_data)
  fft_kernel = fft(kernel_values)

  # Convolution in frequency domain
  n = nyears*12
  convoluted_fft = fft_data * fft_kernel

  # Inverse FFT to obtain the convoluted data
  convoluted_data = Re(fft(convoluted_fft, inverse = TRUE)) / n

  return(convoluted_data )

}

# Function to identify peak months for a species
get_peaks = function(seed_data,nyears) {
  # Find peaks using diff function twice to detect changes 
  # in the the 1st derivative (i.e. get the 2nd derivative)
  peaks = which(diff(sign(diff(seed_data))) == -2) + 1
  
  #The approach with diff will miss a peak if it is in the first
  #month, so do this to fix that: 
  if(length(peaks)<nyears){ peaks = c(1,peaks)}

  return(peaks)
}

#Use the identified peak months to perform a sum around that
#peak.This function assumes a data frame in long format: 
# rows = months
# columns = species
get_seed_sums = function(seed_data, peaks) { 

  nmonths = length(seed_data)[1]
  nyears = floor(nmonths/12)
  npeaks = length(peaks)
  sums_peaks = matrix(0,nmonths,1)
  sums_years = matrix(0,nyears,1)
  #Go to the location of each peak and sum over a 12 month window
  #centered on the peak
  for(p in 1:npeaks){

    #Get the current peak
    ip = peaks[p]

    #Test if we are in year 1 and where the first peak
    #falls. Then sum over an annual window: 
    if (p == 1){ 
      if( ip >=7){ 
        sums_peaks[ip] = sum(seed_data[(ip-6):(ip+5)])
      }else {
        sums_peaks[ip] = sum(seed_data[1:(ip+5)])
      }
    }else{
        sums_peaks[ip] = sum(seed_data[(ip-6):(ip+5)])
      }

    #Put the total in our output vector  
    sums_years[(ip/12+1)] = sums_peaks[ip]

  }

  return( sums_years )

}

######Run the simulated counts for nspecies over nyears: 
nyears = 20
nspecies = 10

#Initialize seeds, peaks, annual counts
seeds = matrix(0,nyears*12,nspecies)
peaks = matrix(0,nyears,nspecies)
annual_counts = matrix(0,nyears,nspecies)

#Generate seeds for each species
for(s in 1:nspecies){ 
  seeds[,s] = get_sp(nyears)
}

#Find peaks and then get annual counts. 
for(s in 1:nspecies){ 
  peaks[,s] = get_peaks(seeds[,s],nyears)
  annual_counts[,s] = get_seed_sums(seeds[,s], peaks[,s])
}

