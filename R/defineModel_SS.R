#' defineModel_SS
#'
#' @details Defines the occupancy model in Nimble for one species. Currently this implements a simple occupancy model.
#' Next steps are to add the random walk prior and the 3 list length options.
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclStateRE should there be a site-level random effect in the state model?
#' @return a set of code
#' @import nimble
#' @export

###########################################################################

defineModel_SS <- function(inclPhenology = TRUE,
                           inclStateRE = FALSE){

  modelcode <- nimbleCode({
    ######################### state model
    for(j in 1:nsite){
        for(t in 1:nyear){
          if(inclStateRE){ # currently only estimates linear trend
            logit(psi[j,t]) <- lam.0 + Trend * (t-1) + eta[j]
          } else {
            logit(psi[j,t]) <- lam.0 + Trend * (t-1)
          }
          z[j,t] ~ dbern(psi[j,t]) # True occupancy status
    }}

    ######################### state model priors
    if(inclStateRE){
      for(j in 1:nsite) {eta[j] ~ dnorm(0, sd=sd.eta)} # site-level random effect
      sd.eta ~ T(dt(0, 1, 1), 0, 10) # constrained
    }
    Trend ~ dnorm(0, tau=0.001)
    lam.0 ~ dnorm(0, tau=3) # highly constrained

    ######################### Obs model
    for(k in 1:nvisit) {
      y[k] ~ dbern(prob = Py[k]) # Observed data
      Py[k] <- z[site[k], year[k]] * p1[k]
      if(inclPhenology){
        logit(p1[k]) <- alpha.0 + alpha.1 * (f_JD[JulDate[k]] - max(f_JD[1:365]))
      } else {
        logit(p1[k]) <- alpha.0
      }
    }

    ######################### Seasonality shared effect
    if(inclPhenology){
      for (d in 1:365){
        f_JD[d] <- 1/((2*3.14159265359)^0.5 * beta2) * exp(-((d - (beta1))^2 / (2* beta2^2)))
      }
    }
    ######################### Obs model priors
    alpha.0 ~ dnorm(-2, tau = 1/2.72) # logit detection probability per pan trap at peak phenology (or mean across year).
    if(inclPhenology){
      alpha.1 ~ T(dt(0, 1, 1), 0, Inf) # constrained to be positive
    } # scaling parameter for detection

    if(inclPhenology){
      beta1 ~ dunif(100, 250) # peak detectability/activity. Not constrained to fall within the field season (c(100, 250))
      beta2 ~ T(dt(0, 1, 1), 10, 200) # Half Cauchy. Stdev of phenology. At sd=500 the curve is entirely flat
    }
    #########################  derived parameters
    for(t in 1:nyear){
      psi.fs[t] <- mean(z[1:nsite,t])
    }
  })
  return(modelcode)
}

