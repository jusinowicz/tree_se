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

# Function to simulate seasonal pattern
get_sp = function(period, peak_month, n) {
  t = 1:n
  seasonal_pattern = sin(2 * pi * (t - peak_month) / period)
  return(seasonal_pattern)
}

# Simulate time series data for 10 plant species
nspecies = 10
full_seeds = list()

for (species in 1:nspecies) {
  seed_data = matrix(0, nrow = years * months, ncol = 1)
  
  for (year in 1:nyears) {
    start_month = (year - 1) * months + 1
    end_month = year * months
    
    # Simulate a peak month for each species
    peak_month = sample(1:months, 1)
    
    # Generate seasonal pattern for the species
    seasonal_pattern = get_sp(period = months, 
      peak_month = peak_month, n = months)
    
    # Add noise to the seasonal pattern
    noise = rnorm(months, mean = 0, sd = 1)
    seed_season = seasonal_pattern + noise
    
    # Assign the data to the respective months
    seed_data[start_month:end_month, 1] = seed_season
  }

  #Could add further transforms here, e.g. to mean. 
  seed_data [seed_data < 0 ] = 0 #Make negatives 0
  full_seeds[[species]] = ts(seed_data, start = c(1, 1), frequency = months)

}

#Convert to data frame (in long format by default)
seeds = data.frame(matrix(unlist(full_seeds), nyears*nmonths, nspecies ))

# Plotting example for a species
plot(seeds[1:36,2], t="l" )

# Convert the data to a time series object (assuming 'Month' is the row index)
library(xts)
seed_counts_xts <- as.xts(seed_counts)

# Function to identify peak months for a specific species
identify_peak_month = function(species_data) {
  # Find peaks using diff function twice to detect changes 
  # in the the 1st derivative (i.e. get the 2nd derivative)
  peaks = which(diff(sign(diff(species_data))) == -2) + 1
  return(peaks)
}

# Apply the function for each species column in the dataset
peak_months = apply(seeds, 2, identify_peak_month)

# Display the peak months for each species
names(peak_months) <- colnames(seed_counts_xts)
print(peak_months)

# Set seed for reproducibility
set.seed(123)

# Number of years
nyears = 20

# Number of plant species
nspecies = 10

# Create an empty data frame to store the simulated data
# This is in long format
seed_data = data.frame(
  Year = rep(1:nyears, each = 12),
  Month = rep(1:12, times = nyears),
  Species = rep(1:nspecies, each = nyears * 12),
  Seed_Count = 0
)

# Simulate seed counts for each month and species
for (i in 1:nspecies) {
  for (j in 1:12) {
    # Generating normal distribution parameters for each species and month
    peak_month = sample(1:12, 1)  # Randomly select the peak month
    mean_count = 100  # Adjust this value for the mean seed count
    sd_count = 20  # Adjust this value for the standard deviation
    
    # Generate seed counts based on a normal distribution
    counts = rnorm(nyears * 12, mean = 0, sd = sd_count)
    
    # Shift peak month to the center and taper on either side
    counts = counts + ifelse((1:12) == peak_month, mean_count, 0)
    
    # Store the generated counts in the data frame
    seed_data$Seed_Count[(j + (i - 1) * 12) + ((1:nyears - 1) * 12)] = counts
  }
}

# Convert seed_data from long to wide format
seed_data_wide = as.data.frame(pivot_wider(
  data = seed_data,
  id_cols = c(Year, Month),
  names_from = Species,
  values_from = Seed_Count
) )


# Example data frame structure (replace this with your actual dataset)
seed_counts = data.frame(
  Month = 1:12,
  Species1 = sample(100:500, 12, replace = TRUE),
  Species2 = sample(100:500, 12, replace = TRUE)
  # Add data for Species3 to Species10 similarly...
)

# Function to find peak month for each species
find_peak_month <- function(species_data) {
  peak_month <- which.max(species_data)
  return(peak_month)
}

# Function to sum seed counts for a specified number of months before and after the peak month
sum_counts_around_peak <- function(species_data, peak_month, months_before, months_after) {
  start_month <- max(peak_month - months_before, 1)
  end_month <- min(peak_month + months_after, length(species_data))
  counts_around_peak <- sum(species_data[start_month:end_month])
  return(counts_around_peak)
}

months_before_peak <- 6
months_after_peak <- 6

# Loop through each species column, find peak month, and sum counts around the peak
for (i in 2:ncol(seed_counts)) { # Assuming Month is the first column
  species <- names(seed_counts)[i]
  peak_month <- find_peak_month(seed_counts[[i]])
  counts_around_peak <- sum_counts_around_peak(seed_counts[[i]], peak_month, months_before_peak, months_after_peak)
  
  cat("For", species, "the peak month is:", peak_month, "\n")
  cat("Total seed counts around the peak (", months_before_peak, "months before and", months_after_peak, "months after) is:", counts_around_peak, "\n\n")
}


# Load necessary libraries
library(dplyr)
library(tidyr)

# Simulate data for 20 years and 10 species
set.seed(123) # Setting seed for reproducibility

#This function will generate a time series of reproduction for each
#species that we can test functions on.

simulate_species = function(start_month, peak_month, duration, years) {
  months = rep(0, 12 * years) # Initialize monthly counts
  
  # Simulate seed counts for the reproduction event
  event_counts = c(0, sample(1:50, duration - 2, replace = TRUE), 0)
  peak_index = (start_month - 1) + which.max(event_counts)
  months[peak_index:(peak_index + length(event_counts) - 1)] <- event_counts
  
  return(months)
}


# Simulate seed counts for each species
nspecies = 10 #Number of species
nyears = 20 #Number of years

#For species phenologies: 
start_month = matrix(0,nspecies,1)
peak_month = matrix(0,nspecies,1)
duration = matrix(0,nspecies,1)

species_data = matrix(0,nyears*12,nspecies)

for (n in 1:nspecies){ 

  #Randomly define each species' phenology: 
  start_month = sample(1:12, 1) # Random start month
  peak_month = sample(1:12, 1) # Random peak month
  duration = sample(3:10, 1) # Random duration of reproduction

  #Generate the time series for each species
  species_data[,n] = simulate_species(start_month, peak_month, duration, years = 20)

}

simulate_species = function(start_month, peak_month, duration, years) {
  months = rep(0, 12 * years) # Initialize monthly counts
  
  # Simulate seed counts for the reproduction event
  event_counts <- c(0, sample(1:50, duration - 2, replace = TRUE), 0)
  peak_index <- (start_month - 1) + which.max(event_counts)
  months[peak_index:(peak_index + length(event_counts) - 1)] <- event_counts
  
  return(months)
}

# Simulate data for 10 species over 20 years
species_data <- lapply(1:10, function(i) {
  start_month <- sample(1:12, 1) # Random start month
  peak_month <- sample(1:12, 1) # Random peak month
  duration <- sample(3:10, 1) # Random duration of reproduction
  
  # Simulate seed counts for each species
  counts <- simulate_species(start_month, peak_month, duration, years = 20)
  return(counts)
} )

# Combine the simulated data into a data frame
simulated_data = data.frame(matrix(unlist(species_data), nrow = 12 * 20, byrow = TRUE))
colnames(simulated_data) = paste0("Species_", 1:10)
rownames(simulated_data) = seq(as.Date("2000-01-01"), by = "month", length.out = 12 * 20)

# Aggregate counts into yearly totals based on peak month
simulated_data_yearly <- simulated_data %>%
  mutate(year = lubridate::year(row.names(.)),
         month = lubridate::month(row.names(.))) %>%
  group_by(year) %>%
  summarize(across(starts_with("Species"), sum)) %>%
  ungroup() %>%
  select(-month)

# Convert the data to wide format with months as rows and species as columns
simulated_data_wide <- spread(simulated_data_yearly, key = "year", value = starts_with("Species"))

# Display the resulting dataset
print(simulated_data_wide)



library(tidyr)

# Function to aggregate counts across years considering wrapping
aggregate_counts = function(wide_data) {

  # Create an empty dataframe to store aggregated yearly counts
  aggregated_counts = data.frame(Year = unique(wide_data$Year), 
                                      stringsAsFactors = FALSE)
  
  # Loop through each row (year) in the wide-format data
  for (i in 1:nrow(wide_data)) {
    year <- wide_data$Year[i]
    # Initialize count for the current year
    year_count <- rep(0, ncol(wide_data) - 1) # Exclude the 'Year' column
    
    # Loop through each species
    for (j in 2:ncol(wide_data)) { # Starting from the 2nd column
      species_count <- wide_data[i, j]
      # Check if the species count should wrap to the previous year
      if (month <- which.max(wide_data[i, 2:ncol(wide_data)])) {
        if (month > 1) {
          year_count[j - 1] <- year_count[j - 1] + species_count
          # Distribute the count to the previous year
          year_count[j - 1] <- year_count[j - 1] + wide_data[i, month + 1] * (1 - ((month - 1) / 12))
        } else {
          year_count[j - 1] <- year_count[j - 1] + species_count
        }
      } else {
        year_count[j - 1] <- year_count[j - 1] + species_count
      }
    }
    
    # Store the aggregated counts for the year
    aggregated_counts[i, -1] <- year_count
  }
  
  return(aggregated_counts)
}
# Simulated data
set.seed(123)
species = rep(letters[1:3], each = 240)  # 20 years * 12 months = 240
year = rep(rep(2000:2019, each = 12), times = 3)
month = rep(1:12, times = 60)
count = round(runif(720, min = 0, max = 100))

# Combine into the dataset. Long format
long_data = data.frame(Species = species, Year = year, Month = month, Count = count)

# Make the long data set wide (species columns, years rows)
wide_data = pivot_wider(long_data, names_from = Species, values_from = Count)


