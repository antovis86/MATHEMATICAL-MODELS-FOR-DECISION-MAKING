# Clear the environment by removing all variables
rm(list=ls())

# Get the directory of the current script for dynamic file referencing
scriptdir = dirname(rstudioapi::getActiveDocumentContext()$path)

# Load the 'dataALL.csv' from the "DATA" folder located in the current working directory
# The .Platform$file.sep ensures cross-platform compatibility for file paths
dataALL = read.csv(paste0(scriptdir, .Platform$file.sep, 'DATA', .Platform$file.sep, 'dataALL.csv'), sep = ",")

# Extract the unique subject IDs from the dataset
subjects = unique(dataALL$ID)

# Initialize matrices to store d' (sensitivity) and c (bias) values for each subject and deltaNB
dPrime = matrix(NA, nrow = length(subjects), ncol = length(deltaNB))
c = matrix(NA, nrow = length(subjects), ncol = length(deltaNB))

# Initialize a vector to store group information for each subject
group = c()

# TODO: Create a loop for processing each subject (currently only processing the first subject)

# Process the first subject (idx = 1)
idx = 1
data = dataALL[dataALL$ID == subjects[idx],]  # Extract data for the current subject
data = data[data$RT > .150,]  # Filter out trials with reaction time less than 0.150 seconds
group[idx] = unique(data$Group)  # Store the group information for the current subject

# Extract the response data and unique deltaNB values (dNum)
resp = data$resp
deltaNB = sort(unique(data$dNum))

# Loop through each unique deltaNB value and calculate signal detection measures
for (dnx in 1:length(deltaNB)) {
  # Calculate hit rate (hr) and false alarm rate (far) for the current deltaNB value
  hr = (sum(resp[data$dNum == deltaNB[dnx] & data$num_rx > data$num_lx]) + 0.5) / 
    (sum(data$dNum == deltaNB[dnx] & data$num_rx > data$num_lx) + 1)
  far = (sum(resp[data$dNum == deltaNB[dnx] & data$num_rx < data$num_lx]) + 0.5) / 
    (sum(data$dNum == deltaNB[dnx] & data$num_rx < data$num_lx) + 1)
  
  # Convert hit rate and false alarm rate to z-scores using the normal distribution
  zhits <- qnorm(hr)
  zFA <- qnorm(far)
  
  # Calculate the signal detection theory measures: d' (sensitivity) and c (bias)
  dPrime[idx, dnx] <- zhits - zFA  # d' is the difference between z-scores
  c[idx, dnx] <- -(zhits + zFA) / 2  # c is the average of the z-scores, indicating bias
}

# TODO: Add a step to test for group differences in d' and c measures

# ROC Curve Analysis ------------------------------------------------------

# TODO: Plot ROC curves for the first subject (idx = 1)
# TODO: Plot ROC curves for all groups (if group information is available)

# Update dPrime to a new set of values for further analysis
dPrime = dPrime[idx, dnx]

# Define a sequence of beta (decision threshold) values ranging from 10 to -10 in increments of -0.1
betas <- seq(10, -10, by = -0.1)

# Initialize vectors to store hit and false alarm rates for each beta value
hits <- rep(NA, length(betas))
false_alarms <- rep(NA, length(betas))

# Loop through each beta value to calculate hit and false alarm rates
for (bx in 1:length(betas)) {
  # Calculate hit rate using the cumulative normal distribution for signal-present trials
  hits[bx] <- pnorm(dPrime - betas[bx])  
  
  # Calculate false alarm rate using the cumulative normal distribution for signal-absent trials
  false_alarms[bx] <- pnorm(-betas[bx])  
}

# Plot the ROC curve based on the calculated hit and false alarm rates
# The x-axis represents false alarm rate and the y-axis represents hit rate
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "False alarms", ylab = "Hits", main = "Model-Based ROC Curve")
lines(false_alarms, hits, col = "red")  # Add the ROC curve to the plot (red color)
