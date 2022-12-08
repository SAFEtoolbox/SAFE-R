# This script provides an application example of the 
# Elementary Effects Test to the HBV rainfall-runoff model.
# Useful to learn about how to use the EET when one of the parameter takes
# discrete values.
#
# METHOD
#
# see description in "workflow_eet_hymod"
#
# MODEL AND STUDY AREA
#
# The model under study is the HBV rainfall-runoff model,
# the inputs subject to SA are the 13 model parameters,
# and the outputs for SA are a set of performance metrics.
# See help of function "hbv_snow_objfun" for more details
# and references about the HBV model and the objective functions.
#
# The case study area is is the Nezinscot River at Turner center, Maine, 
# USA (USGS 01055500, see http://waterdata.usgs.gov/nwis/nwismap)

# This script prepared by Francesca Pianosi and Fanny Sarrazin (for Matlab version)
# Isabella Gollini for the SAFER package,
# University of Bristol, 2014
# mail to: isabella.gollini@bristol.ac.uk


## Step 1 (load the package)

library(SAFER)

## Step 2 (setup the HBV model))

data(hbvdata)

Case  <-  1 # Case = 1: interflow is dominant / Case = 2: percolation is dominant
warmup <- 365 # Model warmup period (days)


dateRange <- c("1948-10-01", "1953-09-30")

t_start <- which(hbvdata$date == dateRange[1])
t_end <- which(hbvdata$date == dateRange[2])
time <- t_start:t_end

dat <- hbvdata[time, 1:5]

X_labels <- c("TS","CFMAX","CFR","CWH","BETA","LP","FC","PERC","K0","K1","K2","UZL","MAXBAS") # uncertain parameters

M <- length(X_labels) # number of uncertain parameters

# Parameter ranges (from Kollat et al.(2012)) 

xmin <- c(-3, 0, 0, 0, 0, 0.3, 1, 0, 0.05, 0.01, 0.05, 0, 1)
xmax <- c(3, 20, 1, 0.8, 7, 1, 2000, 100, 2, 1, 0.1, 100, 6)

# Parameter distributions:
DistrFun  <- c(rep("unif", M-1), "unifd") # discrete uniform for MAXBAS

# Define output:
myfun <- "hbv_snow_objfun"

## Step 3 (sample inputs space)

r <- 100 # Number of samples

# # Option 1: use the sampling method originally proposed by Morris (1991).
# #
# # This requires specifying the number of levels in the uniform grid (L).
# # In this specific case, because one of the parameters (MAXBAS) is 
# # discrete-valued, the value of "L" must coincide with the number of
# # feasible values for the discrete parameter 
# # (6 in this case, ranging from 1 to 6).
 # L <- 6 
# design_type <- "trajectory" # (note used here but required later)
# X <- Morris_sampling(r, xmin, xmax, L) # (r*(M+1),M)

# # Option 2: Latin Hypercube sampling strategy.
# #
# # In this case, we do not need to make any specific choice for the tuning
# # parameters of the sampling strategy, however it may happen that during
# # the sampling we will get several "warning" messages due to the fact that
# # the "OAT_sampling" function checks that each consecutive sampled input
# # vectors differ at least by one component, and if they don't, randomly
# # change one of the two.

SampStrategy <- "lhs" # Latin Hypercube
design_type <- "radial" # other option is "trajectory"

# Convert parameter ranges to the format needed by the OAT_sampling function:
DistrPar <- vector("list", M)
names(DistrPar) <- X_labels
for(i in 1:M){
	DistrPar[[i]] <- c(xmin[i], xmax[i])
}

X <- OAT_sampling(r, M, DistrFun, DistrPar, SampStrategy, design_type)

## From now on just the same as in "workflow_eet_hymod"...

# Step 4 (run the model) 

# WARNING: it may take ~5 minutes to run the model

Y <- model_execution(myfun, X, dat = dat, warmup = warmup, Case = Case)

# Step 5 (Computation of the Elementary effects)

# Choose one among multiple outputs for subsequent analysis:
Yi <- Y[,2] 

#  Compute indices:

EETind <- EET_indices(r, DistrPar, X, Yi, design_type)

EE <- EETind$EE
mi <- EETind$mi
sigma <- EETind$sigma 

dev.new()
EET_plot(mi, sigma, labels = X_labels)

# Use bootstrapping to derive confidence bounds:

Nboot <-100

EETind100 <- EET_indices(r, DistrPar, X, Yi, design_type, Nboot)

EE <- EETind100$EE
mi <- EETind100$mi
sigma <- EETind100$sigma
mi_lb <- EETind100$mi_lb
mi_ub <- EETind100$mi_ub
sigma_lb <- EETind100$sigma_lb
sigma_ub <- EETind100$sigma_ub

# Plot bootstrapping results in the plane (mean(EE),std(EE)):

dev.new()
EET_plot(mi, sigma, mi_lb, mi_ub, sigma_lb, sigma_ub, labels = X_labels)

# Convergence analysis:

rr <- seq(r / 5, r, by = r / 5)
EETconv<- EET_convergence(EE, rr)

m_r <- EETconv$m_r

# Plot the sensitivity measure (mean of elementary effects) as a function 
# of model execution:

dev.new()

plot_convergence(rr * (M + 1), m_r, xlab = "no of model executions", ylab = "mean of EEs", labels = X_labels)

#  Convergence analysis with bootstrapping:
Nboot <- 100
rr <- seq(r / 5, r, by = r / 5)
EETconv100<- EET_convergence(EE, rr, Nboot)

m_r <- EETconv100$m_r
s_r <- EETconv100$s_r
m_lb_r <- EETconv100$m_lb_r
m_ub_r <- EETconv100$m_ub_r

colnames(m_r) <- X_labels

dev.new()

plot_convergence(rr * (M + 1), m_r, m_lb_r, m_ub_r, xlab = "no of model executions", ylab = "mean of EEs", labels = X_labels)
