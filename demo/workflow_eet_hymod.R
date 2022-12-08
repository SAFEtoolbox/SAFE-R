# This script provides a basic application example
# of the Elementary Effects Test. Useful to get started with the EET.
#
# METHOD
#
# This script provides an example of application of the Elementary Effects
# Test (EET) or "method of Morris" (Morris, 1991; Saltelli et al., 2008).
#
# The EET is a One-At-the-Time method for global Sensitivity Analysis.
# It computes two indices for each input:
# i) the mean (mi) of the EEs, which measures the total effect of an input
# over the output;
# ii) the standard deviation (sigma) of the EEs, which measures the degree
# of interactions with the other inputs.
# Both sensitivity indices are relative measures, i.e. their value does not
# have any specific meaning per se but it can only be used in pair-wise
# comparison (e.g. if input x(1) has higher mean EEs than input x(3) than
# x(1) is more influential than x(3)).
#
# For an application example in the environmental domain, see for instance
# Nguyen and de Kok (2007).
#
# MODEL AND STUDY AREA
#
# The model under study is the rainfall-runoff model Hymod
# (see help of function hymod_sim.m for more details) 
# applied to the Leaf catchment in Mississipi, USA
# (Sorooshian et al., 1983).
# The inputs subject to SA are the 5 model parameters, and the scalar 
# output for SA is a metric of model performance.
#
# INDEX
#
# Steps:
# 1. Add paths to required directories
# 2. Load data and set-up the HBV model
# 3. Sample inputs space
# 4. Run the model against input samples 
# 5. Compute the elementary effects
#
# REFERENCES
# 
# Morris, M.D. (1991), Factorial sampling plans for preliminary          
# computational experiments, Technometrics, 33(2).
#
# Nguyen, T.G. and de Kok, J.L. (2007). Systematic testing of an integrated
# systems model for coastal zone management using sensitivity and
# uncertainty analyses. Env. Mod. & Soft., 22, 1572-1587. 
#
# Saltelli, A., et al. (2008) Global Sensitivity Analysis, The Primer,
# Wiley.

# This script prepared by Francesca Pianosi and Fanny Sarrazin (for Matlab version)
# Isabella Gollini for the SAFER package,
# University of Bristol, 2014
# mail to: isabella.gollini@bristol.ac.uk

## Step 1 (load the package)

library(SAFER)

## Step 2 (setup the Hymod model)

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
X_labels <- c("Sm","beta","alfa","Rs","Rf") 

# Define output:
myfun <- "hymod_nse"


## Step 3 (sample inputs space)

r <- 100 # Number of Elementary Effects
# [notice that the final number of model executions will be equal to
# r * (M + 1)]

# # option 1: use the sampling method originally proposed by Morris (1991):
# L <- 6  # number of levels in the uniform grid
# design_type  <- "trajectory" # (note used here but required later)
# X <- Morris_sampling(r, xmin, xmax, L) # (r * (M + 1), M)

# option 2: Latin Hypercube sampling strategy
SampStrategy <- "lhs" # Latin Hypercube
design_type <- "radial"
# other options for design type:
# design_type  = "trajectory"
X <- OAT_sampling(r, M, DistrFun, DistrPar, SampStrategy, design_type)

# Step 4 (run the model) 
Y <- model_execution(myfun, X, dat = dat) # size (r*(M+1),1)

## Step 5 (Computation of the Elementary effects)

# Compute Elementary Effects:
EETind <- EET_indices(r, DistrPar, X, Y, design_type)

EE <- EETind$EE
mi <- EETind$mi
sigma <- EETind$sigma 

# Plot results in the plane (mean(EE),std(EE)):

dev.new()
EET_plot(mi, sigma,  xlab = "Mean of EEs", ylab = "Sd of EEs",  labels = X_labels)

# Use bootstrapping to derive confidence bounds:

Nboot <-100

EETind100 <- EET_indices(r, DistrPar, X, Y, design_type, Nboot)

EE <- EETind100$EE
mi <- EETind100$mi
sigma <- EETind100$sigma
mi_lb <- EETind100$mi_lb
mi_ub <- EETind100$mi_ub
sigma_lb <- EETind100$sigma_lb
sigma_ub <- EETind100$sigma_ub

# Plot bootstrapping results in the plane (mean(EE),std(EE)):
#EET_plot

dev.new()
EET_plot(mi, sigma, mi_lb, mi_ub, sigma_lb, sigma_ub, labels = X_labels)

# Repeat computations using a decreasing number of samples so as to assess
# if convergence was reached within the available dataset:
rr <- seq(r / 5, r, by = r / 5)
EETconv<- EET_convergence(EE, rr)

m_r <- EETconv$m_r

# Plot the sensitivity measure (mean of elementary effects) as a function 
# of model executions:

dev.new()
plot_convergence(rr * (M + 1), m_r, labels = X_labels, xlab = "no of model executions", ylab = "mean of EEs")

# Repeat convergence analysis using bootstrapping:
Nboot <- 100
rr <- seq(r / 5, r, by = r / 5)
EETconv100<- EET_convergence(EE, rr, Nboot)

m_r <- EETconv100$m_r
s_r <- EETconv100$s_r
m_lb_r <- EETconv100$m_lb_r
m_ub_r <- EETconv100$m_ub_r

dev.new()

plot_convergence(rr * (M + 1), m_r, m_lb_r, m_ub_r, xlab = "no of model executions", ylab = "mean of EEs", labels = X_labels)
