# This script applies Variance-Based Sensitivity Analysis to the
# Ishigami-Homma function.
# This function is commonly used to test approximation procedures
# of variance-based indices because its output variance, first-order and 
# total-order indices (or main and total effects) can be analytically
# computed. Therefore this script mainly aims at analyzing the accuracy
# and convergence of the function vbsa_indices. 

# This script prepared by Francesca Pianosi and Fanny Sarrazin (for Matlab version)
# Isabella Gollini for the SAFER package,
# University of Bristol, 2014
# mail to: isabella.gollini@bristol.ac.uk

library(SAFER)

# Setup the model and define input ranges

fun_test  <- "ishigami_homma_function"
M <- 3
distr_fun <- "unif"
distrpar <-  c(-pi, pi)

# Compute the exact values of the output variance (V) and of 
# the first-order (Si_ex) and total-order (STi_ex) 
# variance-based sensitivity indices (this is possible in this
# very specific case because V, Si_ex and STi_ex can be computed
# analytically)

ihfun <- ishigami_homma_function(runif(M))

Si_ex <- attributes(ihfun)$Si_ex
STi_ex <- attributes(ihfun)$STi_ex

# Sample parameter space:

SampStrategy <- "lhs"
N <- 3000
X <- AAT_sampling(SampStrategy, M, distr_fun, distrpar, 2 * N)

# Apply resampling strategy for the efficient approximation of the indices:
XABC <- vbsa_resampling(X)

# Run the model and compute selected model output at sampled parameter
# sets:

YA <- model_execution(fun_test, XABC$XA) # size (N,1)
YB <- model_execution(fun_test, XABC$XB) # size (N,1)
YC <- model_execution(fun_test, XABC$XC) # size (N*M,1)

# Compute main (first-order) and total effects:

ind <- vbsa_indices(YA, YB, YC)

Si <- ind[1,]
STi <- ind[2,] 

# Plot main and total effects and compare the values estimated 
# by the function "vbsa_indices" with the exact values

dev.new()

par(mfcol = c(1, 2))

boxplot2(Si, Si_ex, leg = c("estimated", "exact"), main = "Si")
boxplot2(STi, STi_ex, leg = c("estimated", "exact"), main = "STi")

# Analyze convergence of sensitivity indices:
NN <- seq(N / 5, N, by = N/5)

conv <- vbsa_convergence(c(YA, YB, YC), M, NN)

Si <- conv$Si
STi <- conv$STi

dev.new()
par(mfrow = c(1, 2))

plot_convergence(NN * (M + 2), Si, ex = Si_ex, xlab = "model evals", ylab = "main effect")
plot_convergence(NN * (M + 2), STi, ex = STi_ex, xlab = "model evals", ylab = "total effect")