#' defineModel_SS
#'
#' @details Defines the occupancy model in Nimble for one species. Currently this implements a simple occupancy model.
#' This is very close to the random walk formulation used in sparta
#' There are some minor differences: priors are generally more constrained and `alpha.0` is simplified (`alpha.p` in sparta)
#' @param inclPhenology should the model account for seasonal variation?
#' @param ListLen should list length be modelled as continuous (`cont`), categorical (`cat`) or ignored (`NULL`)
#' @param inclStateRE should there be a site-level random effect in the state model?
#' @return a set of code
#' @import nimble
#' @export

###########################################################################

defineModel_SS <- function(inclPhenology = TRUE,
                           ListLen = NULL,
                           inclStateRE = FALSE){

  modelcode <- nimbleCode({
    ######################### state model
    for(j in 1:nsite){
        for(t in 1:nyear){
          if(inclStateRE){ # currently only estimates linear trend
            logit(muZ[j,t]) <- psi[t] + eta[j]
          } else {
            logit(muZ[j,t]) <- psi[t]
          }
          z[j,t] ~ dbern(muZ[j,t]) # True occupancy status
    }}

    ######################### state model priors
    if(inclStateRE){
      for(j in 1:nsite) {eta[j] ~ dnorm(0, sd = sd.eta)} # site-level random effect
      sd.eta ~ T(dt(0, 1, 1), 0, 10) # constrained
    }
    lam.0 ~ dnorm(0, tau=1) # highly constrained
    sd.psi ~ T(dt(0, 1, 1), 0, Inf) # half Cauchy

    # random walk
    psi[1] <- lam.0
    for(t in 2:nyear){
      psi[t] ~ dnorm(psi[t-1], sd = sd.psi)
    }

    ######################### Obs model
    for(k in 1:nvisit) {
      y[k] ~ dbern(prob = Py[k]) # Observed data
      Py[k] <- z[site[k], year[k]] * p1[k]
      if(inclPhenology){
        phenoComp[k] <- alpha.1 * (f_JD[JulDate[k]] - max(f_JD[1:365]))
      } else {
        phenoComp[k] <- 0
      }
      if(!is.null(ListLen)){
        if(ListLen == "cat") {
          ListLenComp[k] <- gamma.1 * DT2[k] + gamma.2 * DT3[k]
        } else if(ListLen == "cont") {
          ListLenComp[k] <- gamma.1 * log(L[k])
        } else {
          stop("invalid List Length option")
        }
      } else {
        ListLenComp[k] <- 0
      }
      logit(p1[k]) <- alpha.0[year[k]] + phenoComp[k] + ListLenComp[k]
    }

    ######################### Seasonality shared effect
    if(inclPhenology){
      for (d in 1:365){
        f_JD[d] <- 1/((2*3.14159265359)^0.5 * beta2) * exp(-((d - (beta1))^2 / (2* beta2^2)))
      }
    }

    ######################### Obs model priors
    if(inclPhenology){
      alpha.1 ~ T(dt(0, 1, 1), 0, Inf) # scaling parameter for detection: constrained to be positive
      beta1 ~ dunif(100, 250) # peak detectability/activity. Not constrained to fall within the field season (c(100, 250))
      beta2 ~ T(dt(0, 1, 1), 10, 200) # Half Cauchy. Stdev of phenology. At sd=500 the curve is entirely flat
      mu.alpha ~ dnorm(0, tau = 1/2.72) # logit detection probability per pan trap at peak phenology (or mean across year).
    } else {
      mu.alpha ~ dnorm(-2, tau = 1)  # logit detection probability throughout the year
    }

    for(t in 1:nyear){
      alpha.0[t] ~ dnorm(mu.alpha, sd = sd.alpha)
    }
    sd.alpha ~ T(dt(0, 1, 1), 0, 10) # constrained

    if(!is.null(ListLen)){
      if(ListLen == "cat") {
        gamma.1 ~ dnorm(0, tau = 1/2.72)
        gamma.2 ~ dnorm(0, tau = 1/2.72)
      } else if(ListLen == "cont") {
        gamma.1 ~ T(dt(0, 1, 1), 0, Inf) # constrained to be positive
      }
    }

    #########################  derived parameters
    for(t in 1:nyear){
      psi.fs[t] <- mean(z[1:nsite,t])
    }
  })
  return(modelcode)
}

