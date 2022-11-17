# This script provides an application example of the 
# Fourier Amplitude Sensitivity Test (FAST).
# FAST uses the Fourier decomposition of the model output
# to approximates the variance-based first-order sensitivity indices.
# See help of 'FAST_indices' for more details and references.
#
# In this workflow, FAST is applied to the Sobol g-function
# with varying number of inputs and different parameterizations for the
# function parameters not subject to SA
# [see help of 'sobol_g_function'].
# 
# FAST estimates are compared to those obtained by the 'conventional'
# resampling approach used in Variance-Based SA
# [see help of 'vbsa_indices'].
#
# This script can be used to reproduce Figure 2 in:
# Saltelli and Bolado (1998), An alternative way to compute Fourier
# amplitude sensitivity test (FAST), Computational Statistics & Data
# Analysis, 26, 445-460.

# This script prepared by Francesca Pianosi and Fanny Sarrazin (for Matlab version)
# Isabella Gollini for the SAFER package,
# University of Bristol, 2014
# mail to: isabella.gollini@bristol.ac.uk

## Step 1 (load the package)

library(SAFER)

# Define output function:
myfun  <- "sobol_g_function"
# Define input distribution and ranges:
DistrFun <- "unif" 
DistrPar <- c(0, 1)
    
aa <- c(0, 1, 9, 99) # options for the (fixed) parameters
MM <- 5:11 # options for the number of inputs

aMmat <- expand.grid(a = aa, M = MM)

aM <- apply(aMmat, 1, as.list)

Si_ex <- sapply(aM, function(h) attributes(sobol_g_function(runif(h$M), h$a))$Si_ex )

Fsamp <- lapply(1:length(aM), function(n) FAST_sampling(DistrFun, DistrPar, aM[[n]]$M))

Y <- lapply(1:length(aM), function(n) model_execution(myfun, Fsamp[[n]]$X, a = aM[[n]]$a)) # size (Ns,1)

Si_fast <- lapply(1:length(aM), function(n) FAST_indices(Y[[n]], aM[[n]]$M)$Si)

# VBSA

SampStrategy <- "lhs" 
Nfast <- sapply(Y, length)   

# Option 1: set the base sample size for VBSA in such a way that
# the total number of model executions be the same as FAST:

Nvbsa <- ceiling(Nfast / (aMmat[,2] + 2))

# Option 2: set the base sample size to 4096 independently by the
# sample size used for FAST (as done in Saltelli and Bolado, 1998)
# [THIS OPTION MAY BE TIME-CONSUMING!]
#    Nvbsa <- rep(4096, length(Nfast)) 

X <- lapply(1:length(aM), function(n) AAT_sampling(SampStrategy, aM[[n]]$M, DistrFun, DistrPar, 2 * Nvbsa[aM[[n]]$M]))

XABC <- lapply(X, vbsa_resampling)
YA <- lapply(1:length(aM), function(n) model_execution(myfun, XABC[[n]]$XA, a = aM[[n]]$a)) # size (N,1)
YB <- lapply(1:length(aM), function(n) model_execution(myfun, XABC[[n]]$XB, a = aM[[n]]$a)) # size (N,1)
YC <- lapply(1:length(aM), function(n) model_execution(myfun, XABC[[n]]$XC, a = aM[[n]]$a)) # size (N*M,1)

Si_vbsa <- lapply(1:length(aM), function(n) vbsa_indices(YA[[n]],YB[[n]],YC[[n]])[1,])    

SM <- cbind(sapply(Si_ex, mean), sapply(Si_vbsa, mean), sapply(Si_fast, mean))


par(mfrow = c(2,2), mar = c(4, 4, 3, 3))

for(a in aa){

matplot(MM, SM[aMmat[,1] == a,], ylim = c(-.1,1), xlab = "number of inputs M", ylab = "1st-order sensitivity", main = paste("a(i) = ", a, sep =""), pch = 1:3)
legend("topright", c("Analytic", "VBSA", "FAST"), pch = 1:3, col = 1:3, bg = "white")

}


# Plot number of model executions against number of inputs:
dev.new()

plot(MM, unique(Nvbsa) * (MM + 2), col = 2, pch = 2, xlab = "Number of inputs M", ylab = "Number of model executions N")
points(MM, unique(Nfast), col = 3, pch = 3)

legend("topleft", c("VBSA", "FAST"), pch = 2:3, col = 2:3, bg = "white")
