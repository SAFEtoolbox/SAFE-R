# This script applies the PAWN method (Pianosi and Wagener, 2015,2018) to 
# the Ishigami-Homma function.
#
# REFERENCES
#
# Pianosi, F. and Wagener, T. (2018), Distribution-based sensitivity
# analysis from a generic input-output sample, Env. Mod. & Soft., 108, 197-207.
#
# Pianosi, F. and Wagener, T. (2015), A simple and efficient method
# for global sensitivity analysis based on cumulative distribution
# functions, Env. Mod. & Soft., 67, 1-11.
#
# This script was prepared by Valentina Noacco for the SAFER package
# (by Francesca Pianosi and Fanny Sarrazin for the Python version)
# University of Bristol, 2020
# mail to: valentina.noacco@bristol.ac.uk

## Step 1 (load the packages)

library(caTools)
library(calibrater)
library(SAFER)
library(ggplot2)
library(gridExtra)
library(matrixStats)

## Step 2 (setup the model and define input ranges)

fun_test  <- "ishigami_homma_function"
M <- 3
distr_fun <- "unif"
distrpar <-  c(-pi, pi)

## Step 3 (sample inputs space)

SampStrategy <- "lhs"
N <- 2000
X <- AAT_sampling(SampStrategy, M, distr_fun, distrpar, N)

## Step 4 (run the model) 

Y <- model_execution(fun_test, X)

## Step 5 (apply PAWN)

n <- 10 # number of conditioning intervals

x_labels = c('x(1)', 'x(2)', 'x(3)')

# Compute and plot conditional and unconditional CDFs:
dev.new()
pawn_cdf <- pawn_plot_CDF(X, Y, n=10, n_col=3, y_label='output y', labelinput=x_labels)
YF <- pawn_cdf$YF
FU <- pawn_cdf$FU
FC <- pawn_cdf$FC
xc <- pawn_cdf$xc

# Compute and plot KS statistics for each conditioning interval:
KS <- pawn_ks(YF, FU, FC)
dev.new()
KS_all <- pawn_plot_ks(YF, FU, FC, xc, n_col=3, x_labels = x_labels)

# Compute PAWN sensitivity indices:
pawn_ind <- pawn_indices(X, Y, n)

KS_median <- pawn_ind$KS_median
KS_mean <- pawn_ind$KS_mean
KS_max <- pawn_ind$KS_max

# Plot results:
dev.new()
p1 <- boxplot1(as.vector(KS_median), prnam = x_labels) + ylab("KS (median)") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- boxplot1(as.vector(KS_mean), prnam = x_labels) + ylab("KS (mean)") +
  theme(plot.title = element_text(hjust = 0.5))
p3 <- boxplot1(as.vector(KS_max), prnam = x_labels) + ylab("KS (max)") +
  theme(plot.title = element_text(hjust = 0.5))
grid.arrange(grobs = list(p1, p2, p3), ncol = 3)

# Use bootstrapping to derive confidence bounds:
Nboot <- 1000

# Compute sensitivity indices for Nboot bootstrap resamples
# (Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
pawn_ind <- pawn_indices(X, Y, n, Nboot)
KS_median <- pawn_ind$KS_median
KS_mean <- pawn_ind$KS_mean
KS_max <- pawn_ind$KS_max

# KS_median and KS_mean and KS_max have shape (Nboot, M)
# Compute mean and confidence intervals of the sensitivity indices across the
# bootstrap resamples:
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

KS_stat <- KS_median
KS_median_m <- colMeans(KS_stat) # mean
KS_median_lb <-  colQuantiles(KS_stat,probs=alfa/2) # Lower bound
KS_median_ub <- colQuantiles(KS_stat,probs=1-alfa/2) # Upper bound

KS_stat <- KS_mean
KS_mean_m <- colMeans(KS_stat) # mean
KS_mean_lb <-  colQuantiles(KS_stat,probs=alfa/2) # Lower bound
KS_mean_ub <- colQuantiles(KS_stat,probs=1-alfa/2) # Upper bound

KS_stat <- KS_max
KS_max_m <- colMeans(KS_stat) # mean
KS_max_lb <- colQuantiles(KS_stat,probs=alfa/2) # Lower bound
KS_max_ub <- colQuantiles(KS_stat,probs=1-alfa/2) # Upper bound

# Plot bootstrapping results:
dev.new()
p1 <- boxplot1(mu = KS_median_m, lb = KS_median_lb, ub = KS_median_ub, prnam = x_labels) + ylab("KS (median)") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- boxplot1(mu = KS_mean_m, lb = KS_mean_lb, ub = KS_mean_ub, prnam = x_labels) + ylab("KS (mean)") +
  theme(plot.title = element_text(hjust = 0.5))
p3 <- boxplot1(mu = KS_max_m, lb = KS_max_lb, ub = KS_max_ub, prnam = x_labels) + ylab("KS (max)") +
  theme(plot.title = element_text(hjust = 0.5))
grid.arrange(grobs = list(p1, p2, p3), ncol = 3)


# Analyze convergence of sensitivity indices:
NN <- seq(N / 5, N, by = N / 5 )

#( Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
pawn_conv <- pawn_convergence(X, Y, n, NN)
KS_median_c <- pawn_conv$KS_median
KS_mean_c <- pawn_conv$KS_mean
KS_max_c <- pawn_conv$KS_max

KS_median_c <- do.call("rbind",KS_median_c)
KS_mean_c <- do.call("rbind",KS_mean_c)
KS_max_c <- do.call("rbind",KS_max_c)

# Plot convergence
dev.new()
par(mfrow=c(3,1))
plot_convergence(NN, KS_median_c, xlab = "no of model executions", ylab = "KS (median)", labels = x_labels, panel.first = grid())
plot_convergence(NN, KS_mean_c, xlab = "no of model executions", ylab = "KS (mean)", labels = x_labels, panel.first = grid())
plot_convergence(NN, KS_max_c, xlab = "no of model executions", ylab = "KS (max)", labels = x_labels, panel.first = grid())

# Analyze convergence using bootstrapping to derive confidence intervals
#( Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
pawn_conv <- pawn_convergence(X, Y, n, NN, Nboot)

KS_median_c <- pawn_conv$KS_median
KS_mean_c <- pawn_conv$KS_mean
KS_max_c <- pawn_conv$KS_max

# Compute mean and confidence intervals of the sensitivity indices across the
# bootstrap resamples:
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

KS_stat <- KS_max_c
KS_max_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_max_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_max_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

KS_stat <- KS_mean_c
KS_mean_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_mean_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_mean_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

KS_stat <- KS_median_c
KS_median_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_median_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_median_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

# Plot convergence results:
dev.new()
par(mfrow=c(3,1))
plot_convergence(NN, KS_median_c_m, KS_median_c_lb, KS_median_c_ub, xlab = "no of model executions", ylab = "KS (median)", labels = x_labels, panel.first = grid())
plot_convergence(NN, KS_mean_c_m, KS_mean_c_lb, KS_mean_c_ub, xlab = "no of model executions", ylab = "KS (mean)", labels = x_labels, panel.first = grid())
plot_convergence(NN, KS_max_c_m, KS_max_c_lb, KS_max_c_ub, xlab = "no of model executions", ylab = "KS (max)", labels = x_labels, panel.first = grid())


## Step 6 (identification of influential and non-influential inputs)
# This is done by adding an articial 'dummy' input to the list of the model inputs. 
# The sensitivity indices for the dummy parameter estimate the approximation error of the
# sensitivity indices. For reference and more details, see help of the function pawn_indices

# Sensitivity indices using bootstrapping for the model inputs and the dummy input:
# Use bootstrapping to derive confidence bounds:
Nboot <- 1000
# Compute sensitivity indices for Nboot bootstrap resamples. We analyse KS_max
# only (and not KS_median and KS_mean) for screening purposes.
# (Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
pawn_ind <- pawn_indices(X, Y, n, Nboot, dummy = TRUE)

KS_max <- pawn_ind$KS_max # KS_max has dim (Nboot, M)
KS_dummy <- pawn_ind$KS_dummy # KS_dummy has dim (Nboot, 1)

# Compute mean and confidence intervals of the sensitivity indices across the
# bootstrap resamples:
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

KS_stat <- KS_max
KS_max_m <- colMeans(KS_stat) # mean
KS_max_lb <- colQuantiles(KS_stat,probs=alfa/2) # Lower bound
KS_max_ub <- colQuantiles(KS_stat,probs=1-alfa/2) # Upper bound

KS_stat <- KS_dummy
KS_dummy_m <- mean(KS_stat) # mean
KS_dummy_lb <-  quantile(KS_dummy,alfa/2) # Lower bound
KS_dummy_ub <- quantile(KS_dummy,1-alfa/2) # Upper bound

# Combine KS max for all inputs and for dummy to plot
KS_max_d_m <- c(KS_max_m,KS_dummy_m) 
KS_max_d_lb <- c(KS_max_lb,KS_dummy_lb)
KS_max_d_ub <- c(KS_max_ub,KS_dummy_ub)

# Plot bootstrapping results:
dev.new()
boxplot1_dummy(mu = KS_max_d_m, lb = KS_max_d_lb, ub = KS_max_d_ub, prnam = x_labels) + ylab("KS") +
  theme(plot.title = element_text(hjust = 0.5))

# Analyze convergence using bootstrapping to derive confidence intervals
#( Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
NN <- seq(N / 5, N, by = N / 5 )
pawn_conv <- pawn_convergence(X, Y, n, NN, Nboot, dummy = TRUE)

KS_median_c <- pawn_conv$KS_median
KS_mean_c <- pawn_conv$KS_mean
KS_max_c <- pawn_conv$KS_max
KS_dummy_c <- pawn_conv$KS_dummy

# Calculate statistics across bootstrap resamples (mean, lower and upper bounds of sensitivity indices):
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

KS_stat <- KS_median_c
KS_median_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_median_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_median_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

KS_stat <- KS_mean_c
KS_mean_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_mean_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_mean_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

KS_stat <- KS_max_c
KS_max_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_max_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_max_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

KS_stat <- KS_dummy_c
KS_dummy_c_m <- sapply(KS_stat,mean) # mean
KS_dummy_c_lb <-  sapply(KS_dummy_c,quantile,alfa/2) # Lower bound
KS_dummy_c_ub <- sapply(KS_dummy_c,quantile,1-alfa/2) # Upper bound

# Combine KS max for all inputs and for dummy to plot
KS_max_d_c_m <- unname(cbind(KS_max_c_m,KS_dummy_c_m))
KS_max_d_c_lb <- unname(cbind(KS_max_c_lb,KS_dummy_c_lb))
KS_max_d_c_ub <- unname(cbind(KS_max_c_ub,KS_dummy_c_ub))

x_labels_dummy = c('x(1)', 'x(2)', 'x(3)','dummy')

# Plot convergence results:
dev.new()
plot_convergence(NN, KS_max_d_c_m, KS_max_d_c_lb, KS_max_d_c_ub, xlab = "no of model executions", ylab = "KS", labels = x_labels_dummy, panel.first = grid())

