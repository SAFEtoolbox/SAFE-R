#' hbv snow objfun
#'
#' This function simulates the snow accumulation/melting process (via the internal function \code{snow_routine} and the rainfall-runoff process (via the HBV model by Seibert (1997)) and returns 6 objective functions (see Kollat et al, 2002).
#'
#' @param param vector (13) of model parameters. 
#'  Snow routine parameters:
#' \itemize{
#'   \item \code{Ts} = threshold temperature [C]
#'   \item \code{CFMAX} = degree day factor [mm/C]
#'   \item \code{CFR} = refreezing factor [-]
#'   \item \code{CWH} = Water holding capacity of snow [-]
#'   \item \code{temp} = temperature [C]
#' } 
#' HBV parameters: 
#' \itemize{
#'   \item \code{BETA} = Exponential parameter in soil routine [-]
#'   \item \code{LP} = evapotranspiration limit [-]
#'   \item \code{FC} = field capacity [mm]
#'   \item \code{PERC} = maximum flux from Upper to Lower Zone [mm/Dt]
#'   \item \code{K0} = Near surface flow coefficient (ratio) [1/Dt]
#'   \item \code{K1} = Upper Zone outflow coefficient (ratio) [1/Dt]
#'   \item \code{K2} = Lower Zone outflow coefficient (ratio) [1/Dt]
#'   \item \code{UZL} = Near surface flow threshold [mm]
#' }
#' @param dat dataset containing time series of precipitation (\code{prec}), evapotranspiration (\code{ept}), observed flow (\code{flow}) and temperature (\code{temp})

#' @param warmup scalar, warmup time
#' @param Case scalar, 1 or 2, indicates the preferred path. \code{Case = 1}: runoff from the upper zone, \code{Case = 2}: percolation.
#'
#' @return List containing: 
#' \itemize{
#'   \item \code{f} = vector (6) of objective functions 1: AME, 2: NSE, 3: BIAS, 4: TRMSE, 5: SFDCE, 6: RMSE
#'   \item  \code{Q_sim} = vector (\code{N}) time series of simulated flow
#'   \item \code{STATES} = matrix \code{(N, 5)} time series of simulated storages (all in mm) 1: water content of snowpack (snow component) 2: water content of snowpack (liquid component) 3: water content of soil (soil moisture) 4: water content of upper reservoir of flow routing routine 5: water content of lower reservoir of flow routing routine
#'   \item \code{FLUXES} = matrix \code{(N, 8)} time series of simulated fluxes (all in mm/Dt) 1: refreezing 2: snowmelt 3: actual evapotranspiration 4: recharge (water flux from soil moisture accounting module to flow routing module) 5: percolation (water flux from upper to lower reservoir of the flow routing module) 6: runoff from upper reservoir 7: runoff from lower reservoir
#'}
#'
#' @references Seibert,J.(1997)."Estimation of Parameter Uncertainty in the HBV Model". Nordic Hydrology.28(4/5).247-262.
#'
#' Kollat,J.B.,Reed,P.M.,Wagener,T.(2002)."When are multiobjective calibration trade-offs in hydrologic models meaningful?". Water resources research, VOL.48, W03520.

#' @export

hbv_snow_objfun <- function(param, dat, warmup, Case){

prec <- dat$prec
temp <- dat$temp
flow <- dat$flow
ept <- dat$ept

# Comments:
# * Model components: snow routine (optional)- soil moisture, upper zone
# and lower routine - flow routing routine

N <- length(prec)
STATES <- matrix(nrow = N + 1, ncol = 5)
FLUXES<- matrix(nrow = N, ncol = 7)

snr <- snow_routine(temp, prec, param[1:4])

P <- snr$P
STATES[,1:2] <- snr$STATES
FLUXES[,1:2] <- snr$FLUXES

hbvsres <- hbv_sim(P, ept, param[5:13], Case)

Q_sim <- hbvsres$Q_sim
STATES[,3:5] <- hbvsres$STATES
FLUXES[,3:7] <- hbvsres$FLUXES

f <- numeric(6)

Qs <- Q_sim[-(1:warmup)]
Qo <- flow[-(1:warmup)]

N <- length(Qs)

N_67 <- floor(N * 0.67)
N_33 <- floor(N * 0.33)
Qs_sort <- sort(Qs)
Qo_sort <- sort(Qo)

lambda <- 0.3
Zs<- ((1 + Qs)^lambda - 1) / lambda
Zo <- ((1 + Qo)^lambda - 1) / lambda

f[1] <- mean(abs(Qs - Qo)) #AME
f[2] <- 1 - sum((Qs - Qo)^2) / sum((mean(Qo) - Qo)^2) #NSE
f[3] <- abs(mean(Qs - Qo)) #BIAS
f[4] <- sqrt(mean((Zs - Zo)^2))# TRMSE
f[5] <- abs((Qs_sort[N_33] - Qs_sort[N_67])/(Qo_sort[N_33] - Qo_sort[N_67]) - 1) * 100 #SFDCE
f[6] <- sqrt(mean((Qs - Qo)^2)) #RMSE

robj <- f

attributes(robj) <- list(Q_sim = Q_sim, STATES = STATES, FLUXES = FLUXES)

return(robj)

}