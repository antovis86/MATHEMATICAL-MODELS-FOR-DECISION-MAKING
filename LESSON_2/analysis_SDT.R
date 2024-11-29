rm(list=ls())

# Get the directory of the current script for dynamic file referencing
scriptdir = dirname(rstudioapi::getActiveDocumentContext()$path)

# Read the CSV dataset from the "DATA" folder located in the current working directory
data = read.csv(paste0(scriptdir, .Platform$file.sep, 'DATA', .Platform$file.sep, 'Dataset_MR.csv'), sep = ",")

# Remove rows where the 'CompExp.thisRepN' column has missing values (NA) (i.e. take task rows only)
data = data[!is.na(data$CompExp.thisRepN),]

# Extract and sort unique values in the 'N1' column (i.e., standards)
Standards = sort(unique(data$N1))

# Define the short and long probes for each standard based on 'N2' values
S_Probes = sort(unique(data$N2[data$N1 == Standards[1]]))  # Probes related to the first standard (short standard)
L_Probes = sort(unique(data$N2[data$N1 == Standards[2]]))  # Probes related to the second standard (long standard)



# Short Standard Near Probes-------------------------------------------------------------------------
# Calculate the hit rate (probability of correctly identifying a target) for the first short standard (S1)
hits_S1 = (sum(data$targetResp.keys[data$N1 == Standards[1] & data$N2 == S_Probes[3]] == "l") + 0.5) /
  (length(data$targetResp.keys[data$N1 == Standards[1] & data$N2 == S_Probes[3]]) + 1)

# Calculate the false alarm rate (FA) for the first short standard (S1)
FA_S1 = (sum(data$targetResp.keys[data$N1 == Standards[1] & data$N2 == S_Probes[2]] == "l") + 0.5) /
  (length(data$targetResp.keys[data$N1 == Standards[1] & data$N2 == S_Probes[2]]) + 1)

# Compute z-scores for hits and false alarms to calculate d' and c (bias)
zhits_S1 <- qnorm(hits_S1)
zFA_S1 <- qnorm(FA_S1)

# Calculate the signal detection theory measures: d' (sensitivity) and c (bias)
dPrimeS1 <- zhits_S1 - zFA_S1
cS1 <- -(zhits_S1 + zFA_S1) / 2

# Print results for S1
print(paste("S1 d':", dPrimeS1))
print(paste("S1 Bias (c):", cS1))



# Short Standard Far Probes (S2)-------------------------------------------------------------------------

# Repeat the same calculations for the second short standard (S2)
hits_S2 = (sum(data$targetResp.keys[data$N1 == Standards[1] & data$N2 == S_Probes[4]] == "l") + 0.5) /
  (length(data$targetResp.keys[data$N1 == Standards[1] & data$N2 == S_Probes[4]]) + 1)
FA_S2 = (sum(data$targetResp.keys[data$N1 == Standards[1] & data$N2 == S_Probes[1]] == "l") + 0.5) /
  (length(data$targetResp.keys[data$N1 == Standards[1] & data$N2 == S_Probes[1]]) + 1)

zhits_S2 <- qnorm(hits_S2)
zFA_S2 <- qnorm(FA_S2)

# Calculate d' and c for S2
dPrimeS2 <- zhits_S2 - zFA_S2
cS2 <- -(zhits_S2 + zFA_S2) / 2

# Print results for S2
print(paste("S2 d':", dPrimeS2))
print(paste("S2 Bias (c):", cS2))


# Long Standard Near Probes (L1)---------------------------------------------------------------

# Calculate hit rate and FA for the first long standard (L1)
hits_L1 = (sum(data$targetResp.keys[data$N1 == Standards[2] & data$N2 == L_Probes[3]] == "l") + 0.5) /
  (length(data$targetResp.keys[data$N1 == Standards[2] & data$N2 == L_Probes[3]]) + 1)

FA_L1 = (sum(data$targetResp.keys[data$N1 == Standards[2] & data$N2 == L_Probes[2]] == "l") + 0.5) /
  (length(data$targetResp.keys[data$N1 == Standards[2] & data$N2 == L_Probes[2]]) + 1)

zhits_L1 <- qnorm(hits_L1)
zFA_L1 <- qnorm(FA_L1)

# Calculate d' and c for L1
dPrimeL1 <- zhits_L1 - zFA_L1
cL1 <- -(zhits_L1 + zFA_L1) / 2

# Print results for L1
print(paste("L1 d':", dPrimeL1))
print(paste("L1 Bias (c):", cL1))


# Long Standard Far Probes (L2) -------------------------------------------
hits_L2 = (sum(data$targetResp.keys[data$N1 == Standards[2] & data$N2 == L_Probes[4]] == "l") + 0.5) /
  (length(data$targetResp.keys[data$N1 == Standards[2] & data$N2 == L_Probes[4]]) + 1)

FA_L2 = (sum(data$targetResp.keys[data$N1 == Standards[2] & data$N2 == L_Probes[1]] == "l") + 0.5) /
  (length(data$targetResp.keys[data$N1 == Standards[2] & data$N2 == L_Probes[1]]) + 1)

zhits_L2 <- qnorm(hits_L2)
zFA_L2 <- qnorm(FA_L2)

# Calculate d' and c for L2
dPrimeL2 <- zhits_L2 - zFA_L2
cL2 <- -(zhits_L2 + zFA_L2) / 2

# Print results for L2
print(paste("L2 d':", dPrimeL2))
print(paste("L2 Bias (c):", cL2))



# ROC Curve ---------------------------------------------------------------

dPrime = dPrimeL1

# Define a sequence of beta (decision threshold) values from 10 to -10
betas <- seq(10, -10, by = -0.1)

# Initialize vectors for storing hit and false alarm rates
hits <- rep(NA, length(betas))
false_alarms <- rep(NA, length(betas))

# Loop through different beta values to calculate hit and false alarm rates for each
for (bx in 1:length(betas)) {
  # Calculate hit rate using the cumulative normal distribution for signal-present trials
  hits[bx] <- pnorm(dPrime - betas[bx])  
  
  # Calculate false alarm rate using the cumulative normal distribution for signal-absent trials
  false_alarms[bx] <- pnorm(-betas[bx])  
}

# Plot the ROC curve based on the hit and false alarm rates
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "False alarms", ylab = "Hits", main = "Model-Based ROC Curve")
lines(false_alarms, hits, col = "red")
