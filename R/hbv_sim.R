hbv_sim <- function(P, ept, param, Case){

stopifnot(Case %in% 1:2)

# References: 
#
# Seibert,J.(1997)."Estimation of Parameter Uncertainty in the HBV Model".
# Nordic Hydrology.28(4/5).247-262.
#
# Comments:
# * The Capillary flux (from upper tank to soil moisture accounting module)
# is not considered
# * The recharge from the soil to the upper zone is considered to be a
# faster process than evapotranspiration.
# * The preferential path from the upper zone can be modified 
#           - Case 1: interflow is dominant
#           - Case 2: percolation is dominant


# --------------------------
# Recover model parameters:
# --------------------------

BETA <- param[1] # Exponential parameter in soil routine [-]
LP <- param[2] # evapotranspiration limit [-]
FC <- max(.Machine$double.eps, param[3]) # field capacity [mm] cannot be zero

PERC <- param[4] # maximum flux from Upper to Lower Zone [mm/Dt]
K0 <- param[5] # Near surface flow coefficient (ratio) [1/Dt]  
K1 <- param[6] # Upper Zone outflow coefficient (ratio) [1/Dt]  
K2 <- param[7] # Lower Zone outflow coefficient (ratio) [1/Dt]  
UZL <- param[8] # Near surface flow threshold [mm]

MAXBAS <- max(1, round(param[9])) # Flow routing coefficient [Dt]


N <- length(ept) # number of time samples


# ---------------------
# SOIL MOISTURE ROUTINE
# ---------------------

EA <- numeric(N) # Actual Evapotranspiration [mm/Dt]
SM <- numeric(N + 1) #  Soil Moisture [mm]
R  <- numeric(N) # Recharge (water flow from Soil to Upper Zone) [mm/Dt]
UZ <- numeric(N + 1) #  Upper Zone moisture [mm]
LZ <- numeric(N + 1) #  Lower Zone moisture [mm]
RL <- numeric(N) #  Recharge to the lower zone [mm]
Q0 <- numeric(N) # Outflow from Upper Zone [mm/Dt]
Q1 <- numeric(N) #  Outflow from Lower Zone [mm/Dt]

for(t in 1:N){
    
# --------------------------
#    Soil Moisture Dynamics:
# --------------------------

    R[t] <- P[t] * (SM[t] / FC)^BETA # Compute the value of the recharge to the 
    #upper zone (we assumed that this process is faster than evaporation)
    SM_dummy <- max(min(SM[t] + P[t] - R[t], FC), 0) # Compute the water balance 
    # with the value of the recharge  
    R[t] <- R[t] + max(SM[t] + P[t] - R[t] - FC, 0) + min(SM[t] + P[t] - R[t], 0) #adjust R 
    # by an amount equal to the possible negative SM amount or to the 
    # possible SM amount above FC
    
    EA[t] <- ept[t] * min(SM_dummy / (FC * LP), 1) # Compute the evaporation
    SM[t+1] <- max(min(SM_dummy - EA[t], FC), 0) # Compute the water balance 
    
    EA[t] <- EA[t] + max(SM_dummy - EA[t] - FC, 0) + min(SM_dummy - EA[t], 0) # adjust EA
    # by an amount equal to the possible negative SM amount or to the 
    # possible SM amount above FC
     
# --------------------
# Upper Zone dynamics:
# --------------------

    if (Case == 1){
    # Case 1: Preferred path = runoff from the upper zone 
  	    Q0[t] <- max(min(K1 * UZ[t] + K0 * max(UZ[t] - UZL, 0), UZ[t]), 0)     
        RL[t] <- max(min(UZ[t] - Q0[t], PERC), 0)
	} else { #if Case==2
    # Case 2: Preferred path = percolation
        RL[t] <- max(min(PERC, UZ[t]), 0)
        Q0[t] <- max(min(K1 * UZ[t] + K0 * max(UZ[t] - UZL, 0), UZ[t] - RL[t]),0)
    }

  UZ[t+1] <- UZ[t] + R[t] - Q0[t] - RL[t]

# --------------------
# Lower Zone dynamics: 
# --------------------

    Q1[t] <- max(min(K2 * LZ[t], LZ[t]), 0)
    LZ[t + 1] <- LZ[t] + RL[t] - Q1[t]
   
}
 
Q <-  Q0 + Q1 # total outflow (mm/Dt)

# --------------------
# FLOW ROUTING ROUTINE
# --------------------

c <- trimf(1:MAXBAS, 0, (MAXBAS + 1)/2, MAXBAS + 1) #(Seibert,1997)
c <- c / sum(c) # vector of normalized coefficients - (1,MAXBAS)
Q_sim <- numeric(N)
for(t in MAXBAS:N){
    Q_sim[t] <- sum(c * Q[(t - MAXBAS + 1):t])  #(Seibert,1997)
    }


STATES <- cbind(SM, UZ, LZ)
FLUXES <- cbind(EA, R, RL, Q0, Q1)

robj <- list(Q_sim = Q_sim, STATES = STATES, FLUXES = FLUXES)

return(robj)
}

