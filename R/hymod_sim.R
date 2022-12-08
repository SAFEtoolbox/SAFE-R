#' Hymod rainfall-runoff model
#'
#' This function simulates the Hymod rainfall-runoff model
#'
#' flow_sim = time series of simulated flows
#' This function requires that the time series of rainfall (\code{rain})
#' and potential evaporation (\code{evap}) be defined as global variables
#'
#' @param rain vector \code{(T)} time series of rainfall
#' @param evap vector \code{(T)} time series of potential evaporation
#' @param param 5 elements vector of model parameters \code{(Smax, beta, alfa, Rs, Rf)}
#' 
#' @references  Boyle, D. (2001). Multicriteria calibration of hydrological models. PhD thesis, Dep. of Hydrol. and Water Resour., Univ. of Ariz., Tucson.
#'
#' Wagener, T., Boyle, D., Lees, M., Wheater, H., Gupta, H., and Sorooshian, S. (2001). A framework for development and application of hydrological models. Hydrol. Earth Syst. Sci., 5, 13-26.
#' 
#' @export

hymod_sim <- function(rain, evap, param){ 

## --------------------------
## Recover model parameters:
## --------------------------
Sm <- max(.Machine$double.eps, param[1]) # Maximum Soil Moisture (cannot be zero!)
beta <- param[2] # Exponential parameter in soil routine [-]
alfa <- param[3] # Partitioning factor [-]
Rs   <- param[4] # Slow reservoir outflow coefficient (ratio) [1/Dt]  
Rf   <- param[5] # Fast reservoir outflow coefficient (ratio) [1/Dt] 

N_step <- length(rain) # number of time steps in the simulation horizon

## -----------------------
## Initialize variables:
## ---------------------- 
Pe  <- numeric(N_step) # Recharge from the soil [mm/Dt]
Ea  <- numeric(N_step) # Actual Evapotranspiration [mm/Dt]
sm <- numeric(N_step + 1) # Soil Moisture [mm]
sL <- numeric(N_step + 1) # Slow reservoir moisture [mm]
sF1 <- numeric(N_step + 1) # Fast reservoir 1 moisture [mm]
sF2 <- numeric(N_step + 1) # Fast reservoir 2 moisture [mm]
sF3 <- numeric(N_step + 1) # Fast reservoir 3 moisture [mm]

QsL <- numeric(N_step)  # Slow flow [mm/Dt]
QsF <- numeric(N_step)  # Fast flow [mm/Dt]

# simulation:
for (t in 1:N_step){
	
	# --------------------------
	#   Soil Moisture Dynamics:
	# --------------------------

	FF  <- 1 - ( 1 - sm[t] / Sm )^beta
	Pe[t] <- FF * rain[t] # Compute the value of the outflow 
	# (we assumed that this process is faster than evaporation)
	
	sm_temp  = max(min(sm[t] + rain[t] - Pe[t], Sm), 0) # Compute the water 
    # balance with the value of the outflow 
    Pe[t] = Pe[t] + max(sm[t] + rain[t] - Pe[t] - Sm, 0) + min(sm[t] + rain[t] - Pe[t], 0)
    # adjust Pe by an amount equal to the possible negative sm amount or 
    # to the possible sm amount above Sm.
    
    W = min(abs(sm[t] / Sm ),1) # Correction factor for evaporation
    Ea[t] = W * evap[t] # Compute the evaporation
    sm[t+1] = max(min(sm_temp - Ea[t], Sm), 0) # Compute the water balance 
    Ea[t]= Ea[t] + max(sm_temp - Ea[t] - Sm, 0) + min(sm_temp - Ea[t], 0) # adjust Ea 
    # by an amount equal to the possible negative sm amount or to the 
    # possible sm amount above Sm 
   
   # -------------------------
  #   Groundwater Dynamics:
  # -------------------------
   
    # slow flow
    QsL[t] <- Rs * sL[t]
	sL[t+1] <- sL[t] + (1 - alfa) * Pe[t] - QsL[t]
	# fast flow
	sF1[t+1] <- sF1[t] +  alfa * Pe[t] - Rf * sF1[t]
	sF2[t+1] <- sF2[t] +  Rf * sF1[t] - Rf * sF2[t]
	QsF[t]  <- Rf * sF3[t]
	sF3[t+1] <- sF3[t] +  Rf * sF2[t] - QsF[t]
	
}

Qsim <- QsL + QsF

STATES <- list(sm = sm,sL = sL, sF1 = sF1, sF2 = sF2, sf3 = sF3)
FLUXES <- list(Pe = Pe, Ea = Ea, QsL = QsL, QsF = QsF)

attributes(Qsim) <- list(STATES = STATES, FLUXES = FLUXES)

return(Qsim)
}