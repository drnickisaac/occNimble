#' formatData
#'
#' @details Formats data ready for Nimble
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param inData A dataset produced by the simulations
#' @param inclPhenology should the model account for seasonal variation?
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param maxSite defines a limit on the number of sites in the database
#' @return list of two data frames
#' @import reshape2
#' @export

formatData <- function(inData,
                       inclPhenology = TRUE,
                       minSite = 1, maxSite = 999){

  # check that every year has data
  yrs <- min(inData$year):max(inData$year)
  if(any(!yrs %in% inData$year)){
    missingYear <- setdiff(1:attr(inData, "years"), inData$year)
    stop(paste0(missingYear, " has no records in the input data"))
  }

  ### now limit the data to MaxSite whilst preserving the attributes that will be needed later.
  if(maxSite < length(unique(inData$siteID))){
    if(maxSite < 10) {maxSite <- 10}
    print(paste("Subsetting the dataset to", maxSite,"sites"))
    temp <- list(attr(inData, "trend"), attr(inData, "sp_pool"))
    inData <- subset(inData, siteID %in% paste0("site_",1:maxSite))
    attr(inData, "trend") <- temp[[1]]
    attr(inData, "sp_pool") <- temp[[2]]
  }

  castDat <- dcast(inData, year + siteID + survey ~ "nsp",
                   value.var = "species", fun = length, fill = 0)

  # Julian date
  castDat$jday <- as.numeric(format(as.POSIXlt(castDat$survey, format = "%d%b%y"), "%j"))

  # renumber years starting at 1
  castDat$Year <- castDat$year # save a copy of the original
  castDat$year <- 1+ castDat$Year - min(castDat$Year)

  # set the minimum number of sites, as there are some species that were never observed.
  if(minSite < 1) minSite <- 1

  # now restrict the data to species that occur on at least `minSite` sites
  # apparent occupancy matrix across all species:site combos for all data types
  sp_site <- (acast(inData, species~siteID, value.var = "year", function(x) max(x) > 0, fill = 0))

  sp_n_Site <- rowSums(sp_site)
  sp2incl <- which(sp_n_Site > minSite)
  nExcl <- length(sp_n_Site) - length(sp2incl)
  print(paste('Note:',nExcl,'species out of', length(sp_n_Site), 'have been excluded because they occur on', minSite, 'sites or fewer'))
  if(length(sp2incl) > 0) {
    print(paste('We proceed to modelling with', length(sp2incl), 'species'))
  } else {
    stop(paste0("There are no species with enough sites model"))
  }

  # create metadata object
  md <- formatMetadata(inData)

  md$datastr$sp_n_Site <- data.frame(species = names(sp_n_Site), nSite = as.numeric(sp_n_Site))
  md$settings <- list(sp_modelled = length(sp2incl),
                      minSite = minSite,
                      minYr = min(castDat$Year),
                      maxYr = max(castDat$Year))

  dataConstants <- list(nsp = as.numeric(md$settings["sp_modelled"]),
                        nsite = md$dataPars$nSites,
                        nvisit = nrow(castDat),
                        nyear = md$dataPars$nYears,
                        year = castDat$year,
                        site = as.numeric(gsub(castDat$siteID, patt="site_", repl=""))
                        )

  if(inclPhenology){dataConstants$JulDate <- castDat$jday}

  # extract the observations and populate the obsData list
  obsData <- list()

  obsMat <- acast(inData, year + survey + siteID ~ species,
                    value.var = "year", fun = function(x) max(x) > 0, fill=0)
  obsData$y = t(obsMat)[sp2incl,]



  return(list(dataConstants = dataConstants,
              obsData = obsData,
              md = md))
}
