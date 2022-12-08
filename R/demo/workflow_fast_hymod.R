# This script provides an application example of the 
# Fourier Amplitude Sensitivity Test (FAST).
# FAST uses the Fourier decomposition of the model output
# to approximate the variance-based first-order sensitivity indices.
# See help of 'FAST_indices' for more details and references.
#
# In this workflow, FAST is applied to the rainfall-runoff Hymod model
# (see help of 'hymod_sim' for more details) 
# applied to the Leaf catchment in Mississipi, US (Sorooshian et al.,1983).
# The inputs subject to SA are the 5 model parameters, and the scalar 
# output for SA is a performance metric.
# 
# FAST estimates are compared to those obtained by the 'conventional'
# resampling approach used in Variance-Based SA
# [see help of 'vbsa_indices'].
#
# REFERENCES:
#
# Sorooshian, S., Gupta, V., Fulton, J. (1983). Evaluation of maximum 
# likelihood parameter estimation techniques for conceptual rainfall-runoff
# models: Influence of calibration data variability and length on model 
# credibility. Water Resour. Res., 19, 251-259.

# This script prepared by Francesca Pianosi and Fanny Sarrazin (for Matlab version)
# Isabella Gollini for the SAFER package,
# University of Bristol, 2014
# mail to: isabella.gollini@bristol.ac.uk

####################################
# Step 1 (load the package)
####################################

library(SAFER)

####################################
# Step 2: setup the model and define input ranges
####################################

# Load data:
data(LeafCatch)
dat <- LeafCatch[1:365,]

# Number of uncertain parameters subject to SA:
M <- 5

# Parameter ranges (from literature):
DistrPar <- list( c(0, 400), c(0, 2), c(0, 1), c(0, 0.1), c(0.1, 1))

# Parameter distributions:
DistrFun  <- "unif"

# Name of parameters (will be used to costumize plots):
X_labels <- c("Sm", "beta", "alfa", "Rs", "Rf") 

# Define output:
myfun <- "hymod_nse"


######################################
# Step 3: approximate first-order sensitivity indices by FAST
######################################

# FAST sampling:

Fsamp <- FAST_sampling(DistrFun, DistrPar, M)
X <- Fsamp$X
s <- Fsamp$s

# Run the model and compute model output at sampled parameter sets:
Y <- model_execution(myfun, X, dat = dat)
 
Si_fast <- FAST_indices(Y, M)$Si

#####################
# Step 4: Convergence analysis
#####################

# The FAST_sampling function used above automatically sat the sample size
# to the minimum possible given the number M of inputs (see help of
# FAST_sampling to learn about this). In our case this is:

N_fast <- length(Y)

# We can now assess whether FAST estimates would change if using a larger
# number of samples:

NNfast <- seq(N_fast, N_fast + 5000, 2500)

X_NNfast <- lapply(NNfast[-1], function(n) FAST_sampling(DistrFun, DistrPar, M, N = n)$X)

##
# WARNING: Takes long to run! (~ 5min)
##

# uncomment if you have time to wait!

Y_NNfast <- lapply(X_NNfast, function(xn) model_execution(myfun, xn, dat = dat))

Si_NNfast <- sapply(Y_NNfast, function(yn) FAST_indices(yn, M)$Si)

Si_fast_conv <- t(cbind(Si_fast, Si_NNfast))

plot_convergence(NNfast, Si_fast_conv,  xlab = "model evals", ylab = "1st-order sensitivity", main = "FAST", labels = X_labels)

#####################
# Step 5: Comparison with VBSA
#####################

# Here we compare FAST estimates with those obtained by the "conventional"
# resampling approach used in Variance-Based SA
# [see help of "vbsa_indices"].

# Set the base sample size for VBSA in such a way that
# the total number of model executions be the same as FAST:
Nvbsa <- ceiling(NNfast / (M + 2))

# VBSA sampling:

SampStrategy <- "lhs"
X <- lapply(Nvbsa, function(n) AAT_sampling(SampStrategy, M, DistrFun, DistrPar, 2 * n))
XABC <- lapply(X, function(x) vbsa_resampling(x))

# Run the model and compute model output at sampled parameter sets:
	YA <- lapply(XABC, function(x) model_execution(myfun, x$XA, dat = dat)) 
	YB <- lapply(XABC, function(x) model_execution(myfun, x$XB, dat = dat)) 
	YC <- lapply(XABC, function(x) model_execution(myfun, x$XC, dat = dat))
	
# Use output samples to estimate the indices at different sub-sample sizes:
	NNvbsa <- floor(NNfast / (M + 2))
	Si_vbsa_conv <- sapply(1:length(NNvbsa), function(n) vbsa_convergence(c(YA[[n]], YB[[n]], YC[[n]]), M, NNvbsa[n])$Si)

Si_vbsa_conv <- t(Si_vbsa_conv)

# Compare
xlim <- range(NNfast, NNvbsa * (M + 2))

ylim <- range(0, 1, Si_vbsa_conv, Si_fast_conv)

dev.new()
par(mfrow = c(2,1), mar = c(4, 4, 3, 3))

plot_convergence(NNfast, Si_fast_conv, xlim = xlim, ylim = ylim, xlab = "model evals", ylab = "1st-order sensitivity", main = "FAST", labels = X_labels)
plot_convergence(NNvbsa * (M + 2), Si_vbsa_conv, xlim = xlim, ylim = ylim, xlab = "model evals", ylab = "1st-order sensitivity", main = "vbsa", labels = X_labels)
