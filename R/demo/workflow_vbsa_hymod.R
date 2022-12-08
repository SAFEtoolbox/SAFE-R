# This script provides an application example of Variance Based Sensitivity
# Analysis (VBSA)
#
# METHODS
#
# We use two well established variance-based sensitivity indices: 
# - the first-order sensitivity index (or main effects)
# - the total-order sensitivity index (or total effects)
# (see help of 'vbsa_indices' for more details and references)
#
# MODEL AND STUDY AREA
#
# The model under study is the rainfall-runoff model Hymod
# (see help of function hymod_sim for more details) 
# applied to the Leaf catchment in Mississipi, USA
# (Sorooshian et al., 1983).
# The inputs subject to SA are the 5 model parameters, and the scalar 
# output for SA is (one or multiple) performance metric.
#
# INDEX
#
# Steps:
# 1. Add paths to required directories
# 2. Load data, set-up the Hymod model and define input ranges
# 3. Compute first-order (main effects) and total-order (total effects)
#    variance-based indices.
# 4: Example of how to repeat computions after adding up new 
#    input/output samples.
# 5. Example of how to compute indices when dealing with multiple outputs. 
#
# REFERENCES
#
# Sorooshian, S., Gupta, V., Fulton, J. (1983). Evaluation of maximum 
# likelihood parameter estimation techniques for conceptual rainfall-runoff
# models: Influence of calibration data variability and length on model 
# credibility. Water Resour. Res., 19, 251-259.

# This script prepared by Francesca Pianosi and Fanny Sarrazin (for Matlab version)
# Isabella Gollini for the SAFER package,
# University of Bristol, 2014
# mail to: isabella.gollini@bristol.ac.uk

## Step 1: load the package

library(SAFER)
library(ggplot2)
library(cowplot)

## Step 2: setup the model and define input ranges

# Load data:
data(LeafCatch)

dat <- LeafCatch[1:365,]

# Define input distribution and ranges:

M  <- 5 # number of uncertain parameters (Sm beta alfa Rs Rf )
DistrFun  <- "unif" # Parameter distribution
DistrPar  <- list( c(0, 400), c(0, 2), c(0, 1), c(0, 0.1), c(0.1, 1)) #Parameter ranges (from literature)
x_labels <- c("Sm", "beta", "alfa", "Rs", "Rf")

## Step 3: Compute first-order and total-order variance-based indices

myfun <- "hymod_nse"

# Sample parameter space using the resampling strategy proposed by 
# (Saltelli, 2008; for reference and more details, see help of functions
# vbsa_resampling and vbsa_indices) 

SampStrategy <- "lhs"
N <- 3000 # Base sample size.

# Comment: the base sample size N is not the actual number of input 
# samples that will be executed In fact, because of the resampling
# strategy, the total number of model executions to compute the two
# variance-based indices is equal to N*(M+2) 

X <- AAT_sampling(SampStrategy, M, DistrFun, DistrPar, 2 * N)
XABC <- vbsa_resampling(X)

# Run the model and compute selected model output at sampled parameter
# sets:

YA <- model_execution(myfun, XABC$XA, dat = dat) # size (N,1)
YB <- model_execution(myfun, XABC$XB, dat = dat) # size (N,1)
YC <- model_execution(myfun, XABC$XC, dat = dat) # size (N*M,1)

# Compute main (first-order) and total effects:

ind <- vbsa_indices(YA, YB, YC)

Si <- ind[1,]
STi <- ind[2,] 


names(Si) <- x_labels
names(STi) <- x_labels

# Plot results:

# plot main and total separately
dev.new()

p1 <- boxplot1(Si, prnam = x_labels) + ggtitle("Si")
p2 <- boxplot1(STi, prnam = x_labels) + ggtitle("STi")
plot_grid(p1, p2, nrow = 2)

# plot both in one plot
dev.new()

boxplot2(Si, STi, leg = c("main effects", "total effects"), labels = x_labels)

# Check the model output distribution (if multi-modal or highly skewed, the
# variance-based approach may not be adequate):
Y <- c(YA, YC)

dev.new()
par(mfrow = c(1, 2))
plot(ecdf(Y), xlab = "NSE", ylab = "CDF", main ="")
plot(density(Y), xlab = "NSE", ylab = "PDF", main ="")

# Compute confidence bounds:
Nboot <- 500
ind500 <- vbsa_indices(YA,YB,YC, Nboot)

Si <- ind500[1,]
Si_lb <- ind500[3,]
Si_ub <- ind500[4,]

STi <- ind500[5,]
STi_lb <- ind500[7,]
STi_ub <- ind500[8,]

dev.new()

# plot main and total separately
p3 <- boxplot1(mu=Si, lb=Si_lb, ub=Si_ub, prnam = x_labels) + ggtitle("Si")
p4 <- boxplot1(mu=STi, lb=STi_lb, ub=STi_ub, prnam = x_labels) + ggtitle("STi")
plot_grid(p3, p4, nrow = 2)

# plot both in one plot
dev.new()

boxplot2(Si, STi, Si_lb, Si_ub, STi_lb, STi_ub, leg = c("main effects", "total effects"), labels = x_labels)

# Analyze convergence of sensitivity indices:
NN <- seq(N / 5, N, by = N / 5)

vbsaconv <- vbsa_convergence(c(YA, YB, YC), M, NN)

Sic <- vbsaconv$Si
STic <- vbsaconv$STi

dev.new()
par(mfrow = c(1, 2))
plot_convergence(NN * (M + 2), Sic, xlab = "model evals", ylab = "main effect", labels = x_labels)
plot_convergence(NN * (M + 2), STic, xlab = "model evals", ylab = "total effect", labels = x_labels)

# With confidence bounds:

vbsaconv500 <- vbsa_convergence(c(YA, YB, YC), M, NN, Nboot = Nboot)

Sic <- vbsaconv500$Si
STic <- vbsaconv500$STi

Si_lbc <- vbsaconv500$Si_lb
Si_ubc <- vbsaconv500$Si_ub

STi_lbc <- vbsaconv500$STi_lb
STi_ubc <- vbsaconv500$STi_ub

dev.new()
par(mfrow = c(1, 2))
plot_convergence(NN * (M + 2), Sic, Si_lbc, Si_ubc, xlab = "model evals", ylab = "main effect", labels = x_labels)
plot_convergence(NN * (M + 2), STic, STi_lbc, STi_ubc, xlab = "model evals", ylab = "total effect", labels = x_labels)

## Step 4: Adding up new samples

N2 <- 500 # increase of base sample size
# (that means: N2*(M+2) new samples that will need to be executed)

Xext <- AAT_sampling_extend(X, DistrFun, DistrPar, 2 * (N + N2)) # extended sample 

# (it includes the already executed samples X and the new ones)

Xnew <- Xext[-(1:(2 * N)),] # extract the new input samples that need to be executed

 # # Resampling strategy:
XABC2 <- vbsa_resampling(Xnew)

# Execute model against new samples:

YA2 <- model_execution(myfun, XABC2$XA, dat = dat) # size (N2,1)
YB2 <- model_execution(myfun, XABC2$XB, dat = dat) # size (N2,1)
YC2 <- model_execution(myfun, XABC2$XC, dat = dat) # size (N2*M,1)


# Put new and old results toghether:
YAn <- c(YA, YA2) # should have length (N+N2)
YBn <- c(YB, YB2) # should have length (N+N2)
YCn <- rbind(matrix(YC, N, M), matrix(YC2, N2, M)) #  should have size (N+N2,M)
YCn <- c(YCn) # should have length ((N+N2)*M)

# Recompute indices:
Nboot <- 1000 

ind1000 <- vbsa_indices(YAn,YBn,YCn, Nboot)

Sin <- ind1000[1,]
Si_lbn <- ind1000[3,]
Si_ubn <- ind1000[4,]

STin <- ind1000[5,]
STi_lbn <- ind1000[7,]
STi_ubn <- ind1000[8,]

dev.new()

par(mfrow = c(1, 2))

boxplot2(Si, STi, Si_lb, Si_ub, STi_lb, STi_ub, main = paste(N*(M+2), "model eval."),  leg = c("main effects", "total effects"), labels = x_labels)

boxplot2(Sin, STin, Si_lbn, Si_ubn, STi_lbn, STi_ubn, main = paste((N+N2)*(M+2), "model eval."), leg = c("main effects", "total effects"), labels = x_labels)

# Step 5: case of multiple outputs 
# (In this example: RMSE and AME)

myfun <- "hymod_MulObj"
YA <- model_execution(myfun, XABC$XA, dat = dat) # size (N,P)
YB <- model_execution(myfun, XABC$XB, dat = dat) # size (N,P)
YC <- model_execution(myfun, XABC$XC, dat = dat) # size (N*M,P)

# select the j-th model output:
j <- 1 
ind1 <- vbsa_indices(YA[, j], YB[, j],YC[, j])

Si1 <- ind1[1,]
STi1<- ind1[2,]

j <- 2
ind2 <- vbsa_indices(YA[, j], YB[, j],YC[, j])

Si2 <- ind2[1,]
STi2<- ind2[2,]
 
dev.new()
par(mfrow = c(1, 2))
boxplot2(Si1, STi1, leg = c("main effects", "total effects"), main = "RMSE", labels = x_labels)
boxplot2(Si2, STi2, leg = c("main effects", "total effects"), main = "BIAS", labels = x_labels)

# If you want to add samples in this case:
N2 <- 500 # increase of base sample size (see previous Step)
Xext <- AAT_sampling_extend(X, DistrFun, DistrPar, 2*(N+N2)) # extended sample 
Xnew <- Xext[-(1:(2 * N)),] # extract the new input samples that need to be executed
# Resampling strategy:
XABC2 <- vbsa_resampling(Xnew)
# Execute the model against new samples:
YA2 <- model_execution(myfun, XABC2$XA, dat = dat) # size (N2,2)
YB2 <- model_execution(myfun, XABC2$XB, dat = dat) # size (N2,2)
YC2 <- model_execution(myfun, XABC2$XC, dat = dat) # size (N2*M,2)

# Select the j-th model output:
j <- 1 

# Put new and old results toghether:
YAn <- c(YA[,j], YA2[,j]) # should have length (N+N2)
YBn <- c(YB[,j], YB2[,j]) # should have length (N+N2)
YCn <- rbind(matrix(YC[,j], N, M), matrix(YC2[,j], N2, M)) #  should have size (N+N2,M)
YCn <- c(YCn) # should have length ((N+N2)*M)


ind1n <- vbsa_indices(YAn, YBn, YCn)

Si1n <- ind1n[1,]
STi1n<- ind1n[2,]

dev.new()
par(mfrow = c(1, 2))
boxplot2(Si1, STi1, main = paste(N*(M+2), "model eval."),  leg = c("main effects", "total effects"), labels = x_labels)
boxplot2(Si1n, STi1n, main = paste((N+N2)*(M+2), "model eval."), leg = c("main effects", "total effects"), labels = x_labels)





