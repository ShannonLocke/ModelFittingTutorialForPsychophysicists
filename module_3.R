# Module 3
#
# Created by SML Dec 2020

# Preamble:
library("dplyr")
library("rstan")
library("shinystan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set working directory:
dataPath <- "/Users/shannonlocke/Documents/Research/Library/Teaching/GradPsychophysics/2020/module_3"
if (!dir.exists(dataPath)) {dataPath <- getwd}
setwd(dataPath)

#######################################
# PART 1: FITTING INDIVIDUAL OBSERVER #
#######################################

# Read raw response data:
fetchFile <- paste(dataPath, '/simIndivObs_module3.txt', sep = "")
indata <- read.table(fetchFile, header = T)
nTrials <- length(indata$stimidx)

# MCMC sampling instructions:
iter <- 2000
chains <- 4

# MCMC sampling boundaries:
dprime_range <- c(0, 4)
k_range <- c(-2,2)

# Create list for rstan script:
data <- list(nTrials = nTrials,
             stimidx = indata$stimidx,
             resp = indata$resp,
             dprime_range = dprime_range,
             k_range = k_range)

# Fit stan model:
fit <- stan("model_SDT_basic.stan", data=data, iter=iter, chains=chains)

# See fit summary and output samples:
print(fit)
plot(fit)
traceplot(fit)
pairs(fit)
launch_shinystan(fit)

#######################################
# PART 2: FITTING MULTIPLE OBSERVERS #
#######################################

# Read raw response data:
fetchFile <- paste(dataPath, '/simMultObs_module3.txt', sep = "")
indata <- read.table(fetchFile, header = T)
nTrials <- length(indata$stimidx)
resp <- as.matrix(indata)
resp <- resp[,-1]
nObs <- ncol(resp)
  
# MCMC sampling instructions:
iter <- 2000
chains <- 4

# MCMC sampling boundaries:
dprime_range <- c(0, 4)
k_range <- c(-2,2)

# Create list for rstan script:
data <- list(nObs = nObs,
             nTrials = nTrials,
             stimidx = indata$stimidx,
             resp = resp,
             dprime_range = dprime_range,
             k_range = k_range)

# Fit stan model:
fitMult <- stan("model_SDT_hierarchical.stan", data=data, iter=iter, chains=chains)

# See fit summary and output samples:
print(fitMult)
launch_shinystan(fitMult)
sims <- rstan::extract(fitMult)

# Store results in data frame:
res <- data.frame(subjID = seq(1:nObs))
res$dHierarchical <- colMeans(sims$dprime)
res$kHierarchical <- colMeans(sims$k)

###############################################
# PART 3: COMPARE FITS: BASIC VS HIERARCHICAL #
###############################################

# Fit each subject individually:
for (sel_subj in seq(1,nObs)){
  data <- list(nTrials = nTrials,
               stimidx = indata$stimidx,
               resp = resp[,sel_subj],
               dprime_range = dprime_range,
               k_range = k_range)
  fit <- stan("model_SDT_basic.stan", data=data, iter=iter, chains=chains)
  sims <- rstan::extract(fit)
  res$dBasic[sel_subj] <- mean(sims$dprime)
  res$kBasic[sel_subj] <- mean(sims$k)
}

# Export data:
fname <- "/data_compareFits.txt"
newFile <- paste(dataPath, fname, sep = "")
write.table(res, newFile, quote = FALSE, row.names = FALSE, sep=" ")
