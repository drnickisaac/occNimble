#' runModel
#'
#' @details Runs an occupancy model for multiple species.
#' @param dataConstants dataframe produced by formatData()
#' @param obsData dataframe produced by formatData()
#' @param dataSumm$stats dataframe produced by formatData()
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclStateRE should there be a site-level random effect in the state model?
#' @param multiSp should the model be run as a multispecies model, or many single-species models?
#' @param parallelize should the chains be run as separate processes on different cores?
#' @param allPars if `TRUE` then all model parameters are monitored. If `FALSE`, just `lam.0bda` and `Trend`.
#' @param n.iter number of iterations for the Nimble model. Default is 1000.
#' @param n.burn number of iterations for the burn-in. If `NULL` (the default), then it will be set to `n.iter/2`.
#' @param n.thin thinning for the MCMC chains. Defaults to 5
#' @param n.chain number of MCMC chains. Defaults to 3
#' @param maxSp maximum number of species to be modelled
#' @return a set of year effects
#' @export
#' @import nimble
#' @import pbmcapply
#' @import parallel
#' @import coda
#' @import reshape2
#' @import lme4

#' #' @examples
#' \dontrun{
#'
#' set.seed(123)
#'
#' # Create data
#' n <- 15000 #size of dataset
#' nyear <- 20 # number of years in data
#' nSamples <- 100 # set number of dates
#' nSites <- 50 # set number of sites
#'
#' # Create somes dates
#' first <- as.Date(strptime("2010/01/01", format="%Y/%m/%d"))
#' last <- as.Date(strptime(paste(2010+(nyear-1),"/12/31", sep=''), format="%Y/%m/%d"))
#' dt <- last-first
#' rDates <- first + (runif(nSamples)*dt)
#'
#' # taxa are set as random letters
#' taxa <- sample(letters, size = n, TRUE)
#'
#' # taxa are set as random letters: weight so that some species are more commonly observed
#' taxa <- sample(letters, size = n, TRUE, prob = seq(from = 0.01, to = 0.4, length.out = 26))
#'
#' # sites are visited randomly
#' site <- sample(paste('site_', 1:nSites, sep=''), size = n, TRUE)
#'
#' # Format the data
#' inData <- data.frame(species = taxa,
#'                      siteID = site,
#'                      survey = survey,
#'                      year = as.numeric(format(as.Date(survey, format="%d/%m/%Y"),"%Y")))
#'
#' # prepare the data for the model (includes removing species found on few sites)
#' formattedData <- formatData(inData)
#'
#' # run the model with these data for one species (very small number of iterations)
#' results <- runModel(taxa_name = taxa[1],
#'                       n_iterations = 50,
#'                       burnin = 15,
#'                       occDetdata = visitData$occDetdata,
#'                       spp_vis = visitData$spp_vis,
#'                       write_results = FALSE,
#'                       provenance  = "sparta test dataset")
#' }
#'
###############################################################

runModel <- function(dataConstants,
                     obsData,
                     dataSumm,
                     useNimble = TRUE,
                     inclPhenology = TRUE,
                     inclStateRE = FALSE,
                     multiSp = FALSE,
                     parallelize = FALSE,
                     allPars = FALSE,
                     n.iter = 1000,
                     n.burn = NULL,
                     n.thin = 5,
                     n.chain = 3,
                     maxSp = 9999){

  ###################################################################

  if(grepl("mingw32", sessionInfo()$platform) & parallelize == TRUE){
    warning("Parallelization not yet implemented for Windows machines: setting to FALSE")
    parallelize <- FALSE
  }

  ###################################################################
  # kludge required: if maxSp is set to 1 (or zero) then problems arise later.
  if(maxSp < 2) maxSp <- 2

  # truncate the dataset if there are too many species
  if(dim(obsData$y)[1] > maxSp){
    obsData <- lapply(obsData, function(x) x[1:maxSp,])
    dataSumm$occMatrix <- dataSumm$occMatrix[1:maxSp,,]
    dataSumm$stats <- dataSumm$stats[1:maxSp,]
    dataConstants$nsp <- maxSp
    print(paste('Warning: only the first', maxSp, 'will be used in modelling: others will be ignored'))
  }

    ###################################################################
  if(useNimble) {
    if(is.null(n.burn)) n.burn = n.iter/2

    if(multiSp == TRUE){ # Multispecies option

      # step 1 define the model code
      modelcode <- defineModel_MS(inclPhenology = inclPhenology,
                                  inclStateRE = inclStateRE)

      init.vals <- list(z = dataSumm$occMatrix,
                        lam.0 = cloglog(dataSumm$stats$naiveOcc),
                        gamma.0 = rep(ilogit(0.2), times=maxSp),
                        Trend = rnorm(n=1),
                        spTr = rnorm(n=maxSp),
                        tau.trend = 1)
      if(inclPhenology){
        init.vals$beta1 <- rep(180, times=maxSp)
        init.vals$beta2 <- rep(50, times=maxSp)
        init.vals$gamma.1 <- rep(1, times=maxSp)
      }
      if(inclStateRE){
        init.vals$sd.eta <- 2
        init.vals$eta <- rnorm(n=dataConstants$nsite, mean=0, sd=2)
      }

      # step 2 create an operational from from NIMBLE/BUGS code
      model <- nimbleModel(code = modelcode,
                           constants = dataConstants,
                           data = obsData,
                           inits = init.vals)

      params <- c("mu.lambda","Trend")
      if(allPars) {
        params <- c(params, 'lam.0','gamma.0', 'psi.fs', 'tau.trend')
        if(inclPhenology) params <- c(params, "beta1", "beta2", 'gamma.1')
        if(inclStateRE) params <- c(params, "sd.eta")
      }

      # step 3 build an MCMC object using buildMCMC(). we can add some customization here
      occMCMC <- buildMCMC(model,
                           monitors = params,
                           thin = n.thin,
                           useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible

      # step 3 before compiling the MCMC object we need to compile the model first
      Cmodel <- compileNimble(model)

      # test whether the model is fully initialised
      if(is.na(Cmodel$calculate())) {stop("model not fully initialized")}
      Cmodel$initializeInfo()

      # now the MCMC (project = NIMBLE model already associated with a project)
      CoccMCMC <- compileNimble(occMCMC, project = model)

      # and now we can use either $run or runMCMC() on the compiled model object.
      if(parallelize){
        av_cores <- parallel::detectCores() - 1
        runMCMC_samples <- pbmcapply::pbmclapply(1:n.chain, function(i)
          runMCMC(
            mcmc = CoccMCMC,
            nburnin = n.burn,
            niter = n.iter,
            nchains = 1, samplesAsCodaMCMC = T),
          mc.cores = av_cores)

      } else {
        runMCMC_samples <- runMCMC(CoccMCMC,
                                   nburnin = n.burn,
                                   niter = n.iter,
                                   nchains = n.chain, samplesAsCodaMCMC = T)
      }
      yearEff <- runMCMC_samples

      ############################################ end multispecies

    } else { # sequential single-species option
      # step 1 define the model code
      modelcode <- defineModel_SS(inclPhenology = inclPhenology,
                                  inclStateRE = inclStateRE)

      init.vals <- list(z = dataSumm$occMatrix[1,,], # value for species 1
                        lam.0 = ilogit(dataSumm$stats$naiveOcc)[1], # value for species 1
                        gamma.0 = ilogit(0.2),
                        Trend = rnorm(n=1))

      if(inclPhenology){
        init.vals$beta1 <- 180
        init.vals$beta2 <- 50
        init.vals$gamma.1 <- 1
      }
      if(inclStateRE){
        init.vals$sd.eta <- 2
        init.vals$eta <- rnorm(n=dataConstants$nsite, mean=0, sd=2)
      }

      # step 2 create an operational from from NIMBLE/BUGS code
      model <- nimbleModel(code = modelcode,
                           constants = dataConstants[!names(dataConstants) %in% "nsp"],
                           data = lapply(obsData, function(x) x[1,]), # values for species 1
                           inits = init.vals)

      # step 3 build an MCMC object using buildMCMC(). we can add some customization here
      params <- c("mu.lambda","Trend")
      if(allPars) {
        params <- c(params, 'lam.0','gamma.0', 'psi.fs')
        if(inclPhenology) params <- c(params, "beta1", "beta2", 'gamma.1', "alpha.1")
        if(inclStateRE) params <- c(params, "sd.eta")
      }

      occMCMC <- buildMCMC(model,
                           monitors = params,
                           useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible

      # step 3 before compiling the MCMC object we need to compile the model first
      Cmodel <- compileNimble(model)

      # now the MCMC (project = NIMBLE model already associated with a project)
      CoccMCMC <- compileNimble(occMCMC, project = model)

      ####################################################################################

      single_species_model <- function(sp, spDat, dataSumm,
                                       n.iter, n.burn, n.thin, n.chain,
                                       Cmodel, CoccMCMC){

        # apparent occupancy for this species
        Z <- dataSumm$occMatrix[sp,,]

        # write an informative message about this species' data
        nS <- sum(rowSums(Z>0, na.rm=T)>0)
        spName <- dataSumm$stats$species
        print(paste0("Now running ", spName[sp], ", which is present on ", nS, " sites"))

        # add the data for the species of interest
        Cmodel$setData(spDat)

        # finish initialization
        spInits <- list(z = Z,
                        lam.0 = cloglog(dataSumm$stats$naiveOcc)[sp])
        if(inclPanTrap) spInits$alpha.0 = ilogit(dataSumm$stats$reportingRate_z1[sp]) # replace with reportingRate_1 when I can calculate it
        Cmodel$setInits(spInits)

        # test whether the model is fully initialised
        if(is.na(Cmodel$calculate())) {stop("model not fully initialized")}
        Cmodel$initializeInfo()

        # and now we can use $run on the compiled model object.
        samplesList <- list()
        for(i in 1:n.chain){
          CoccMCMC$run(niter = n.iter,
                       nburnin = n.burn,
                       chain = i,
                       thin = n.thin,
                       reset = TRUE)
          samplesList[[i]] <- as.matrix(CoccMCMC$mvSamples)
        }
        samplesList <- coda::as.mcmc.list(lapply(samplesList, as.mcmc))
      }

      ####################################################################################

      ####### run the model for each species

      if(parallelize){
        #av_cores <- parallel::detectCores() - 1
        yearEff <- pbmcapply::pbmclapply(1:maxSp, function(i){
          single_species_model(sp=i,
                               spDat=lapply(obsData, function(x) x[i,]),
                               dataSumm=dataSumm,
                               n.iter = n.iter,
                               n.burn = n.burn,
                               n.thin = n.thin,
                               n.chain = n.chain,
                               Cmodel, CoccMCMC)
        },
        mc.cores = getOption("mc.cores", 7L)  #av_cores
        )
      } else {
        yearEff <- lapply(1:maxSp, function(i){
          single_species_model(sp=i,
                               spDat=lapply(obsData, function(x) x[i,]),
                               dataSumm = dataSumm,
                               n.iter = n.iter,
                               n.burn = n.burn,
                               n.thin = n.thin,
                               n.chain = n.chain,
                               Cmodel, CoccMCMC)
        }
        )
      }
      names(yearEff) <- dimnames(dataSumm$occMatrix)[[1]][1:maxSp]
      attr(yearEff, "modelCode") <- model$getCode()
    }

    #####################################################################

  }
  else {
    # for simplicity, let's just report the annual total count across all data types
    totalObs <- sapply(1:maxSp, function(i) rowSums(sapply(obsData, function(x) x[i,])))

    mData <- with(dataConstants, data.frame(site=site, year=year))
    mData <- melt(cbind(mData, totalObs), id=1:2)
    names(mData)[3] <- "species"
    # NB mData has one row per round

    yearEff <- t(sapply(unique(mData$species), function(sp){
      spDat <- subset(mData, species == sp)
      #subset to the sites with >0 observations
      occSites <- which(acast(spDat, site~., value.var="value", sum) > 0)
      mod <- glm(value ~ year + factor(site),
                 data = subset(spDat, site %in% occSites),
                 family = "poisson")
      return(as.numeric(summary(mod)$coefficients["year",1:2]))
    }))
    dimnames(yearEff)[[2]] <- c("Estimate", "Std. Error")
    dimnames(yearEff)[[1]] <- paste0("species", dimnames(obsData$y2)[[1]])

  }
  return(yearEff)
}
