rm(list=ls())
library(WienR)
# Get the directory of the current script for dynamic file referencing
scriptdir = dirname(rstudioapi::getActiveDocumentContext()$path)

# Read the CSV dataset from the "DATA" folder located in the current working directory
dataALL = read.csv(paste0(scriptdir, .Platform$file.sep, 'DATA', .Platform$file.sep, 'dataALL.csv'), sep = ",")

# Extract the unique subject IDs from the dataset
subjects = 

  
# Initialize results list  
p_ddm_all <- list()

# TODO: Create a loop for processing each subject (currently only processing the first subject)
  

# Define the function to be minimized
fit_DDM_log <- function(params, RT, RESP,N1,N2) {
  # Transform parameters to their native space
  a <- exp(params[1])
  v1 <- exp(params[2])
  w <- 1 / (1 + exp(-params[3]))
  Ter = min(RT) / (1 + exp(-params[4])) # non decision time can be between 0 and min(RT)
  RT <- pmax(.Machine$double.eps, RT - Ter)
  v = v1*log(N2/N1)
  
  # Initialize llik as a vector of NA values
  llik <- rep(NA, length(RT))
  # Loop over each trial
  for (ntrial in 1:length(RT)) {
    P <- tryCatch({
      if (RESP[ntrial] == 1) {
        WienerPDF(t = RT[ntrial], response = "upper", a = a, v = v[ntrial], w = w)$value
      } else {
        WienerPDF(t = RT[ntrial], response = "lower", a = a, v = v[ntrial], w = w)$value
      }
    }, error = function(e) {
      message(paste("Error in trial", ntrial, ":", e$message))
      return(0)  # Return 0 if there's an error
    })
    
    # If P > 0, store log(P), otherwise store log(eps)
    if (!is.na(P) && P > 0) {
      llik[ntrial] <- log(P)
    } else {
      llik[ntrial] <- log(.Machine$double.eps)
    }
  }
  # Sum response likelihoods
  loglik <- sum(llik)  # Remove NA values to prevent errors
  
  return(-loglik)  # Negative log-likelihood for minimization
}


# Process the first subject (idx = 1)
idx = 1  
# Remove rows where the 'CompExp.thisRepN' column has missing values (NA) (i.e. take task rows only)
data = dataALL[dataALL$ID==subjects[idx],]
data=data[data$RT>.150,]
rt = data$RT
resp = data$resp





# Define initial points
x0 <- c(0, 0, 0,0)

result <- optim(
  par = x0, 
  fn = fit_DDM_log, 
  RT = rt, 
  RESP = resp,
  N1 = data$num_lx,
  N2 = data$num_rx,
  method = "BFGS",  # Similar to 'quasi-newton' in MATLAB
)


# Extract results and transform parameters back to their native space
p_ddm <- list(
  a = exp(result$par[1]),
  v = exp(result$par[2]),
  w = 1 / (1 + exp(-result$par[3])),
  Ter = min(rt) / (1 + exp(-result$par[4])) 
)

# Display results
cat("Results:\n\n")
cat("Parameter estimates:\n")
print(p_ddm)

# Store the p_ddm list in the p_ddm_list
p_ddm_all[[idx]] <- p_ddm

