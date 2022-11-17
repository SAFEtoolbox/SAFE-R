snow_routine <- function(temp, prec, param){


###########################
# Recover model parameters:
###########################

Ts    <- param[1] # threshold temperature [C]
CFMAX <- param[2] # degree day factor [mm/C]
CFR   <- param[3] # refreezing factor [-]
CWH   <- param[4] # Water holding capacity of snow [-]

N <- length(prec) #number of time samples

###########################
# SNOWPACK ROUTINE
###########################

P <- numeric(N) # snowmelt leaving the snowpack/recharge to the soil [mm/Dt]
  
rain <- prec
rain[temp < Ts] <- 0  # [mm/Dt]
snow <- prec
snow[temp >= Ts] <- 0 # [mm/Dt]
Ta <- temp - Ts
Ta[temp < Ts] <- 0 # Active Temperature for snowmelt
Tn <- Ts-temp
Tn[temp >= Ts] <- 0 # Active Temperature for refreezing
m  <- numeric(N) # snowmelt [mm/Dt]
rfz <- numeric(N) # refreezing [mm/Dt]
v  <- matrix(0, N + 1, 1) # snowpack depth [mm]: solid component
vl <- matrix(0, N + 1, 1) # snowpack depth [mm]: liquid component   

for(t in 1:N){
        
    m[t]  <- min(CFMAX * Ta[t], v[t])                        
    rfz[t] <- min(CFR * CFMAX * Tn[t], vl[t])                   
    #   snowpack dynamics: solid component
    v[t + 1] <- v[t] - m[t] + snow[t] + rfz[t] 
    #   snowpack dynamics: liquid component
    vl[t + 1] <- vl[t] + m[t] + rain[t] - rfz[t] 
    if (vl[t + 1] > CWH * v[t + 1]){ # if the liquid component exceed the snow pack
        # holding capacity
        P[t] <- vl[t + 1] - CWH * v[t + 1]
        vl[t+1] <- CWH * v[t + 1]
        }
}

STATES <- cbind(v, vl)
FLUXES <- cbind(rfz, m)    

robj <- list(P = P, STATES = STATES, FLUXES = FLUXES)

return(robj)

}