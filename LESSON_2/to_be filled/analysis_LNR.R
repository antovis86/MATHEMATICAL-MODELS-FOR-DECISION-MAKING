rm(list=ls())
library(WienR)
# Get the directory of the current script for dynamic file referencing
scriptdir = dirname(rstudioapi::getActiveDocumentContext()$path)

# Read the CSV dataset from the "DATA" folder located in the current working directory
dataALL = read.csv(paste0(scriptdir, .Platform$file.sep, 'DATA', .Platform$file.sep, 'dataALL.csv'), sep = ",")

subj_ID = 1
# Remove rows where the 'CompExp.thisRepN' column has missing values (NA) (i.e. take task rows only)
data = dataALL[dataALL$ID==subj_ID,]
data=data[data$RT>.150,]
rt = data$RT
resp = data$resp

# Define the function to be minimized
fit_LNR_log <- function(params, RT, RESP,N1,N2) {
  # Transform parameters to their native space
  a <- params[1]
  v <- params[2]
  sigma <- exp(-params[3])
  Ter = min(RT) / (1 + exp(-params[4])) # non decision time can be between 0 and min(RT)
  RT <- pmax(.Machine$double.eps, RT - Ter)
  mu1 = a - v*log(N2/N1)
  mu2 = a - v*log(N2/N1)
  
  # Initialize llik as a vector of NA values
  llik <- rep(NA, length(RT))
  # Loop over each trial
  for (ntrial in 1:length(RT)) {
    if (RESP[ntrial]==1) {
      mu_pdf = mu1[ntrial]
      mu_cdf = mu2[ntrial]
    } else {
      mu_pdf = mu2[ntrial]
      mu_cdf = mu1[ntrial]
    }
    
    pdf_value <- dlnorm(RT[ntrial], mu_pdf, sigma)
    survival_value <- 1 - plnorm(RT[ntrial], mu_cdf, sigma)
    P = pdf_value*survival_value
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

# Define initial points
x0 <- c(0, 0, 0,0)

result <- optim(
  par = x0, 
  fn = fit_LNR_log, 
  RT = rt, 
  RESP = resp,
  N1 = data$num_lx,
  N2 = data$num_rx,
  method = "BFGS",  # Similar to 'quasi-newton' in MATLAB
)


# Extract results and transform parameters back to their native space
p_lnr <- list(
  a = result$par[1],
  v = result$par[2],
  sigma = exp(result$par[3]),
  Ter = min(rt) / (1 + exp(-result$par[4])) 
)

# Display results
cat("Results:\n\n")
cat("Parameter estimates:\n")
print(p_lnr)
