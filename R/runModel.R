#' runModel
#'
#' @details Runs an occupancy model for multiple species.
#' @param dataConstants dataframe produced by formatData()
#' @param obsData dataframe produced by formatData()
#' @param dataSumm$stats dataframe produced by summariseData()
#' @param format Either "Nimble" (default) or "spOcc"
#' @param ListLen should be included as continuous ("Cont"), categorical ("Cat") or excluded (NULL)
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclStateRE should there be a site-level random effect in the state model?
#' @param multiSp should the model be run as a multispecies model, or many single-species models?
#' @param parallelize should the chains be run as separate processes on different cores?
#' @param allPars Either a list of parameters to monitor or logical statement: if `TRUE` then all model parameters are monitored. If `FALSE`the just `Trend`.
#' @param n.iter number of iterations for the Nimble model. Default is 1000.
#' @param n.burn number of iterations for the burn-in. If `NULL` (the default), then it will be set to `n.iter/2`.
#' @param n.thin thinning for the MCMC chains. Defaults to 5
#' @param n.chain number of MCMC chains. Defaults to 3
#' @param maxSp maximum number of species to be modelled
#' @return a set of year effects
#' @import nimble
#' @import pbmcapply
#' @import parallel
#' @import coda
#' @import reshape2
#' @import spOccupancy

#' @examples
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
#' # Summarise the data
#' dataSumm <- with(formattedData, summariseData(obsData, dataConstants))
#'
#' # run the model with these data for the first three species in the dataset
#' results <- runModel(formattedData$dataConstants,
#'                    formattedData$obsData,
#'                    dataSumm,
#'                    maxSp = 3)
#' }
#' @export
#'
###############################################################

runModel <- function(dataConstants,
                     obsData,
                     dataSumm,
                     format = "Nimble",
                     inclPhenology = FALSE,
                     ListLen = NULL,
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

  ####### the number of species to be modelled
  nSpMod <- min(dim(obsData$y)[1], maxSp)

  if(is.null(n.burn)) n.burn <- round(n.iter/2)

  ###################################################################
  # truncate the dataset if there are too many species
  if(dim(obsData$y)[1] > nSpMod){
    obsData$y <- obsData$y[1:nSpMod,]
    dataSumm$occMatrix <- dataSumm$occMatrix[1:nSpMod,,]
    dataSumm$stats <- dataSumm$stats[1:nSpMod,]
    dataConstants$nsp <- nSpMod
    print(paste('Warning: only the first', nSpMod, 'will be used in modelling: others will be ignored'))
  }

  # put the species names in a vector
  spNames <- dimnames(formattedData$obsData$y)[[1]][1:nSpMod]

  ###################################################################

  if(format == "Nimble") {

    #####################################
    # define the model parameters and which should be monitored
    modPars <- c('lam.0', 'psi.fs', "sd.psi", "alpha.0", "mu.alpha", "sd.alpha")
    # work out what other parameters there are
    if(inclPhenology) modPars <- c(modPars, "alpha.1", "beta1", "beta2")
    if(inclStateRE) modPars <- c(modPars, "sd.eta")
    if(!is.null(ListLen)){
      modPars <- c(modPars, "gamma.1")
      if(ListLen == "cat") modPars <- c(modPars, "gamma.1")
    }

    if(all(is.logical(allPars))){
      if(allPars == TRUE) {
        params <- modPars
      } else {
        params <- c("lam.0")
      }
    } else {
      # check that the manually supplied set of parameters is valid
      params <- allPars
      if(!all(params %in% modPars)){
        badPars <- setdiff(params, modPars)
        warning(paste(badPars, "not recognised"))
        if(all(params %in% badPars)) {
          stop("no valid parameters listed")
        } else {
          params <- setdiff(params, badPars)
        }
      }
    }
    print(paste("Monitoring:", params))

    #####################################

    if(multiSp == TRUE){ # Multispecies option - not edited for simple occupancy

      stop("Multispecies not yet implemented")

    } else { # sequential single-species option

      # step 1 define the model code
      modelcode <- defineModel_SS(inclPhenology = inclPhenology,
                                  inclStateRE = inclStateRE)
      #print("model definition read in")

      init.vals <- list(z = dataSumm$occMatrix[1,,], # value for species 1
                        lam.0 = logit(dataSumm$stats$naiveOcc[1] * 0.99), # value for species 1
                        sd.psi = rnorm(n = 1, sd = 0.02),
                        psi = dnorm(n=dataConstants$nyear, sd = 0.02),
                        mu.alpha = -1,
                        sd.alpha = 2,
                        alpha.0 = rnorm(n=dataConstants$nyear, mean=0, sd=2)
                        )

      if(inclPhenology){
        init.vals$beta1 <- 180
        init.vals$beta2 <- 50
        init.vals$alpha.0 <- 0
      } else {
        init.vals$alpha.0 <- -2
      }
      if(inclStateRE){
        init.vals$sd.eta <- 2
        init.vals$eta <- rnorm(n=dataConstants$nsite, mean=0, sd=2)
      }
      if(!is.null(ListLen)) {
        if(ListLen == "cont"){
          init.vals$gamma.1 <- 0.1
        } else if(ListLen == "cat"){
          init.vals$gamma.1 <- 0.1
          init.vals$gamma.2 <- 0.1
        } else {
          stop("invalid List Length option")
        }
      }
      #print("initial values set")

      # step 2 create an operational from from NIMBLE/BUGS code
      model <- nimbleModel(code = modelcode,
                           constants = dataConstants[!names(dataConstants) %in% "nsp"],
                           data = lapply(obsData, function(x) x[1,]), # values for species 1
                           inits = init.vals)
      #print("step 2 complete")

      # step 3 build an MCMC object using buildMCMC(). we can add some customization here
      occMCMC <- buildMCMC(model,
                           monitors = params,
                           useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible
      #print("model build complete")

      # step 4 before compiling the MCMC object we need to compile the model first
      Cmodel <- compileNimble(model)
      #print("model compilation complete")

      # now the MCMC (project = NIMBLE model already associated with a project)
      CoccMCMC <- compileNimble(occMCMC, project = model)
      #print("compilation of MCMC object complete")
      ####################################################################################

      single_species_model <- function(sp, spDat, dataSumm,
                                       n.iter, n.burn, n.thin, n.chain,
                                       Cmodel, CoccMCMC){

        # apparent occupancy for this species
        Z <- dataSumm$occMatrix[sp,,]

        # write an informative message about this species' data
        nS <- sum(rowSums(Z, na.rm=T)>0)
        spName <- dataSumm$stats$species[sp]
        print(paste0("Now running ", spName, ", which is present on ", nS, " sites"))

        # add the data for the species of interest
        Cmodel$setData(spDat)

        # finish initialization
        spInits <- list(z = Z,
                        lam.0 = logit(dataSumm$stats$naiveOcc[sp] * .99)) # to avoid numeric problem
        Cmodel$setInits(spInits)

        # test whether the model is fully initialised
        if(is.na(Cmodel$calculate())) {
          print(Cmodel$initializeInfo())
          stop("model not fully initialized")
          }

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

      if(parallelize){
        #av_cores <- parallel::detectCores() - 1
        yearEff <- pbmcapply::pbmclapply(1:nSpMod, function(i){
          single_species_model(sp = i,
                               spDat = lapply(obsData, function(x) x[i,]),
                               dataSumm = dataSumm,
                               n.iter = n.iter,
                               n.burn = n.burn,
                               n.thin = n.thin,
                               n.chain = n.chain,
                               Cmodel, CoccMCMC)
        },
        mc.cores = getOption("mc.cores", 7L)  #av_cores
        )
      } else {
        yearEff <- lapply(1:nSpMod, function(i){
          single_species_model(sp = i,
                               spDat = lapply(obsData, function(x) x[i,]),
                               dataSumm = dataSumm,
                               n.iter = n.iter,
                               n.burn = n.burn,
                               n.thin = n.thin,
                               n.chain = n.chain,
                               Cmodel, CoccMCMC)
        }
        )
      }
      #names(yearEff) <- dimnames(dataSumm$occMatrix)[[1]][1:nSpMod] # should be same as below
      names(yearEff) <- dimnames(formattedData$obsData$y)[[1]][1:nSpMod]
      attr(yearEff, "modelCode") <- model$getCode()
    }

    #####################################################################

  }
  else if(format == "spOcc") {

    # set priors and initial values
    priors <- list(alpha.normal = list(mean = 0, var = 5),
                   beta.normal = list(mean = 0, var = 5))

    inits <- list(alpha = 0,
                  beta = 0)

    # model specification for the occupancy and detectability
    occ.formula <- ~ as.factor(Site) + as.factor(Year)
    if(!is.null(ListLen)){
      if(grepl("cont", ListLen, ignore.case = TRUE)){
        det.formula <- ~ as.factor(year) + logL
      } else if(grepl("cat", ListLen, ignore.case = TRUE)){
        det.formula <- ~ as.factor(year) + as.factor(DT2) + as.factor(DT3)
      }
    } else {
      det.formula <- ~ as.factor(year)
    }

    if(inclPhenology){warning("phenology option not yet implemented")}

    single_species_spOcc <- function(sp, y, spName=NULL,
                                     dat,
                                     inits,
                                     n.iter, n.burn, n.thin, n.chain,
                                     occ.formula,
                                     det.formula){
      dat$y <- y
      Z <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
      inits$z <- Z

      # write an informative message about this species' data
      nS <- sum(rowSums(Z, na.rm=T)>0)
      if(!is.null(spName))
        print(paste0("Now running ", spName, ", which is present on ", nS, " sites"))

      # MCMC settings
      # I find it easier to think in terms of samples (to be comparable to sparta)
      # so I create that first and work back to batch size
      n.batch <- 200
      batch.length <- n.iter/n.batch

      # now fit
      out <- tPGOcc(occ.formula = occ.formula,
                    det.formula = det.formula,
                    data = dat,
                    inits = inits,
                    priors = priors,
                    n.batch = n.batch,
                    batch.length = batch.length,
                    n.omp.threads = n.chain,
                    verbose = TRUE,
                    n.report = 1000,
                    n.burn = n.burn,
                    n.thin = n.thin,
                    n.chain = n.chain)
    }

    if(parallelize){
      #av_cores <- parallel::detectCores() - 1
      yearEff <- pbmcapply::pbmclapply(1:nSpMod, function(i){
        obsData$scaff$focal <- obsData$scaff$variable == spNames[i]
        obsData$scaff <- unique(obsData$scaff[, c("siteID", "year", "Replicate", "focal")])
        spDat <- acast(obsData$scaff, siteID ~ year ~ Replicate, value.var = "focal",
                       fun = function(x) ifelse(length(x) > 0, max(x), -999))
        spDat[spDat < 0] <- NA
        single_species_spOcc(sp = i,
                             y = spDat,
                             spName = spNames[i],
                             dat = dataConstants,
                             inits = inits,
                             n.iter = n.iter,
                             n.burn = n.burn,
                             n.thin = n.thin,
                             n.chain = n.chain,
                             occ.formula = occ.formula,
                             det.formula = det.formula)
      },
      mc.cores = getOption("mc.cores", 7L)  #av_cores
      )
    } else {
      yearEff <- lapply(1:nSpMod, function(i){
        obsData$scaff$focal <- obsData$scaff$variable == spNames[i]
        obsData$scaff <- unique(obsData$scaff[, c("siteID", "year", "Replicate", "focal")])
        spDat <- acast(obsData$scaff, siteID ~ year ~ Replicate, value.var = "focal",
                       fun = function(x) ifelse(length(x) > 0, max(x), -999))
        spDat[spDat < 0] <- NA
        single_species_spOcc(sp = i,
                             y = spDat,
                             spName = spNames[i],
                             dat = dataConstants,
                             inits = inits,
                             n.iter = n.iter,
                             n.burn = n.burn,
                             n.thin = n.thin,
                             n.chain = n.chain,
                             occ.formula = occ.formula,
                             det.formula = det.formula)
      }
      )
    }
    names(yearEff) <- spNames
    attr(yearEff, "modelCode") <- list(occ = occ.formula, det = det.formula)

  #####################################################################

  } else {
    stop("format not known")
  }
  return(yearEff)
}
