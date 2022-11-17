# This script provides an application example of Regional Sensitivity
# Analysis (RSA)
#
# METHODS
#
# This script provides an example of application of Regional Sensitivity
# Analysis (RSA). RSA is applied in two different ways:
#
# - according to the original formulation by Spear and Hornberger (1980)
# where the input samples are split into two datasets depending on whether
# the corresponding output satisfies a threshold condition (for an
# application example see for instance Sieber and Uhlenbrook, 2005)
# functions: RSA_thres, RSA_plot_thres, RSA_convergence_thres
#
# - splitting the input samples into 'ngroup' datasets corresponding to 
# equally spaced ranges of the output (see for instance Tang et al., 2007)
# functions: RSA_groups, RSA_plot_groups
#
# MODEL AND STUDY AREA
#
# The model under study is the rainfall-runoff model Hymod (Boyle, 2001; 
# Wagener et al., 2001) applied to the Leaf catchment in Mississipi, USA
# (Sorooshian et al., 1983).
# The inputs subject to SA are the 5 model parameters, and the scalar 
# output for SA is (one or multiple) performance metric.
#
# INDEX
#
# Steps:
# 1. Load the package
# 2. Load data and set-up the Hymod model
# 3. Sample inputs space
# 4. Run the model against input samples
# 5a. Perform RSA with thresholds 
# or
# 5b. Perform RSA with groups
#
# If an input/output sample is already available, skip steps 1-4 and go to
# step 5a/b (see comment marked as (**)).
#
# REFERENCES
# 
# Boyle, D. (2001). Multicriteria calibration of hydrological models. 
# PhD thesis, Dep. of Hydrol. and Water Resour., Univ. of Ariz., Tucson.
#
# Sieber, A., Uhlenbrook, S. (2005). Sensitivity analyses of a distributed
# catchment model to verify the model structure, Journal of Hydrology, 
# 310(1?4), 216-235.
#
# Sorooshian, S., Gupta, V., Fulton, J. (1983). Evaluation of maximum 
# likelihood parameter estimation techniques for conceptual rainfall-runoff
# models: Influence of calibration data variability and length on model 
# credibility. Water Resour. Res., 19, 251?259.
#
# Spear, R.C. and Hornberger, G.M. (1980). Eutrophication in peel inlet,
# II, identification of critical uncertianties via generalized sensitivity
# analysis, Water Resour. Res., 14, 43?49.
#
# Tang, Y., Reed, P., Wagener, T., and van Werkhoven, K. (2007). Comparing 
# sensitivity analysis methods to advance lumped watershed model 
# identification and evaluation, Hydrol. Earth Syst. Sci., 11, 793-817.

# Wagener, T., Boyle, D., Lees, M., Wheater, H., Gupta, H., and Sorooshian, 
# S. (2001). A framework for development and application of hydrological 
# models. Hydrol. Earth Syst. Sci., 5, 13?26.

# This script prepared by Francesca Pianosi and Fanny Sarrazin (for Matlab version)
# Isabella Gollini for the SAFER package, edited by Valentina Noacco
# University of Bristol, 2020
# mail to: valentina.noacco@bristol.ac.uk

## Step 1 (load the package)

library(calibrater)
library(SAFER)
library(ggplot2)
library(matrixStats)

## Step 2 (setup the Hymod model)

# Load data:
data(LeafCatch)

dat <- LeafCatch[1:365,]

# Define inputs:
DistrFun  <- "unif" # Parameter distribution
DistrPar  <- list( c(0, 400), c(0, 2), c(0, 1), c(0, 0.1), c(0.1, 1)) #Parameter ranges (from literature)
x_labels <- c("Sm", "beta", "alfa", "Rs", "Rf")

# Define output:
myfun <- "hymod_MulObj"

## Step 3 (sample inputs space)

SampStrategy <- "lhs" # Latin Hypercube
N <- 3000 # Number of samples
M <- length(DistrPar) # Number of inputs
X <- AAT_sampling(SampStrategy, M, DistrFun, DistrPar, N)
colnames(X) <- x_labels

## Step 4 (run the model) 
Y <- model_execution(myfun, X, dat = dat)  # size (N,2)
colnames(Y) <- c("rmse", "bias")

## Step 5a (Regional Sensitivity Analysis with threshold)

# (**) Note: if you want to use input/output samples generated in
# another programme, them here and save them in two matrix 
# input : X = (N x M)
# output: Y = (N x P)
# [N=number of samples; M=number of inputs; P=number of outputs]

# Visualize input/output samples (this may help finding a reasonable value
# for the output threshold):

scatter_plots(X,Y[,1]) + ylab("rmse")

# use dev.new() if you want to open the plot in a new window

dev.new()

scatter_plots(X,Y[,2]) + ylab("bias")

# Set output threshold:
rmse_thres <- 3    #  threshold for the first obj. fun.
bias_thres <- 0.5  # behavioural threshold for the second obj. fun.

# RSA (find behavioural parameterizations):
threshold <- c(rmse_thres, bias_thres)

rsatr <- RSA_indices_thres(X, Y, threshold) 
mvd <- rsatr$stat
idxb <- rsatr$idxb

# Highlight the behavioural parameterizations in the scatter plots:

dev.new()
scatter_plots_tr(X, Y[,1], prnam = x_labels, idxb) + ylab("rmse")

dev.new()
scatter_plots_tr(X, Y[,2], prnam = x_labels, idxb) + ylab("bias")

# Plot parameter CDFs:
dev.new()
RSA_plot_thres(X, idxb, prnam = x_labels, threshold = threshold, str_legend =c("behav", "non-behav")) # add legend


# Check the ranges of behavioural parameterizations by
# Parallel coordinate plot:

mycol <- idxb
mycol[idxb == FALSE] <- gray(.7, alpha = .7)
mycol[idxb == TRUE] <-  gray(0, alpha = .7)

dev.new()
parcoord(X, col = mycol, plotorder = idxb)


# Plot the sensitivity indices (maximum vertical distance between
# parameters CDFs):

# border sets the colors of the boxplots, boxwex sets the width of the boxplot, axes = FALSE and add = TRUE allow to draw the boxplot over the parcoord plot.

dev.new()
boxplot1(mu = mvd, prnam = x_labels) + ylab("mvd") 

# Compute sensitivity indices with confidence intervals using bootstrapping

Nboot <- 1000
rsatr_b <- RSA_indices_thres(X, Y, threshold, Nboot = Nboot) 

mvd <- rsatr_b$stat
idxb <- rsatr_b$idxb

mvd_lb <- rsatr_b$stat_lb
mvd_ub <- rsatr_b$stat_ub

# Plot results:

dev.new()
boxplot1(mu = mvd, lb = mvd_lb, ub = mvd_ub, prnam = x_labels) + ylab("mvd") 

# Repeat computations using an increasing number of samples so to assess
# convergence:
NN <- seq(N / 5, N, by = N / 5 )
mvd <- RSA_convergence_thres(X, Y[,1], NN, threshold = rmse_thres) 

mvd_st <- mvd$stat

# Plot the sensitivity measures (maximum vertical distance between
# parameters CDFs) as a function of the number of samples:

dev.new()
plot_convergence(NN, mvd_st, xlab = "no of samples", ylab = "mvd", labels = x_labels)

# Repeat convergence analysis using bootstrapping to derive
# confidence bounds:

Nboot <- 1000
rsatr_b <- RSA_convergence_thres(X, Y[,1], NN,  threshold = rmse_thres, Nboot = Nboot) 

mvd <- rsatr_b$stat
idxb <- rsatr_b$idxb

mvd_lb <- rsatr_b$stat_lb
mvd_ub <- rsatr_b$stat_ub

dev.new()
plot_convergence(NN, mvd, mvd_lb, mvd_ub, xlab = "no of samples", ylab = "mvd", labels = x_labels)


## Step 5b (Regional Sensitivity Analysis with groups)

# RSA (find behavioural parameterizations):

rsa_gr <- RSA_indices_groups(X, Y[,1])

mvd_median <- rsa_gr$mvd_median
mvd_mean <- rsa_gr$mvd_mean
mvd_max <- rsa_gr$mvd_max
spread_median <- rsa_gr$spread_median
spread_mean <- rsa_gr$spread_mean
spread_max <- rsa_gr$spread_max
idx <- rsa_gr$idx
Yk <- rsa_gr$Yk

# Plot parameter CDFs:
dev.new()
RSA_plot_groups(X, idx, Yk, prnam = x_labels) + ylab("rmse")

# Compute sensitivity indices with confidence intervals using bootstrapping
Nboot <- 1000
ngroup <- 10
rsa_gr_b <- RSA_indices_groups(X, Y[,1], ngroup, Nboot)

# Statistics across all bootstrap resamples
mvd_median <- rsa_gr_b$mvd_median
mvd_mean <- rsa_gr_b$mvd_mean
mvd_max <- rsa_gr_b$mvd_max
spread_median <- rsa_gr_b$spread_median
spread_mean <- rsa_gr_b$spread_mean
spread_max <- rsa_gr_b$spread_max
idx <- rsa_gr_b$idx
Yk <- rsa_gr_b$Yk

# Compute mean and confidence intervals of the sensitivity indices across the
# bootstrap resamples:
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

mvd_median_m <- colMeans(mvd_median) # median
mvd_median_lb <-  colQuantiles(mvd_median,probs=alfa/2) # Lower bound
mvd_median_ub <- colQuantiles(mvd_median,probs=1-alfa/2) # Upper bound

mvd_mean_m <- colMeans(mvd_mean) # mean
mvd_mean_lb <-  colQuantiles(mvd_mean,probs=alfa/2) # Lower bound
mvd_mean_ub <- colQuantiles(mvd_mean,probs=1-alfa/2) # Upper bound

mvd_max_m <- colMeans(mvd_max) # max
mvd_max_lb <-  colQuantiles(mvd_max,probs=alfa/2) # Lower bound
mvd_max_ub <- colQuantiles(mvd_max,probs=1-alfa/2) # Upper bound

# Plot results:

dev.new()
boxplot1(mu = mvd_median_m, lb = mvd_median_lb, ub = mvd_median_ub, prnam = x_labels) + ylab("mvd median") 
