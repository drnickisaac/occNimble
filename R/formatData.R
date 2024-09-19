#' formatData
#'
#' @details Formats data ready for Nimble
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param inData A dataset with columns "species", "siteID", "survey" and "year"
#' @param format Either "Nimble" (default) or "spOcc"
#' @param ListLen should be included as continuous ("Cont"), categorical ("Cat") or excluded (NULL)
#' @param inclPhenology should the model account for seasonal variation?
#' @param minYrPerSite the minimum number of years with data for a site to be included (defaults to 2, as in sparta)
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param minRecs the threshold minimum number of records for a species to be considered (after applying `minSite`)
#' @return list of two data frames
#' @import reshape2
#' @import magrittr
#' @import dplyr
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
#' }
#'
#' @export

formatData <- function(inData,
                       format = "Nimble",
                       ListLen = NULL,
                       inclPhenology = TRUE,
                       minYrPerSite = 2,
                       minSite = 1,
                       minRecs = 1,
                       verbose = TRUE){

  # check that every year has data
  yrs <- min(inData$year):max(inData$year)
  if(any(!yrs %in% inData$year)){
    missingYear <- setdiff(yrs, inData$year)
    stop(paste0(missingYear, " has no records in the input data"))
  }

  ################################################
  ### subset the data to sites occurring in N years
  #(NB sparta does this at a later stage)

  if(minYrPerSite > 1){
    temp <- inData %>%
      distinct(siteID, year) %>%
      group_by(siteID) %>%
      count()

  sites_to_include <- subset(temp, n >= minYrPerSite)$siteID

  if(verbose){
    print(paste(length(unique(inData$siteID)) - length(sites_to_include),
              "sites out of",
              length(unique(inData$siteID)),
              "have visits from fewer than",
              minYrPerSite,
              "years and are being discarded."))
    }

    inData <- subset(inData, siteID %in% sites_to_include)
  }

  ########################################################
  # summarise the raw data at species level

  spRec <- inData %>%
    distinct(species, survey, siteID, year) %>%
    group_by(species) %>%
    count()

  spSite <- inData %>%
    distinct(species, siteID) %>%
    group_by(species) %>%
    count()

  spSumm <- as.data.frame(spRec)
  names(spSumm)[2] <- "recs"
  spSumm$sites <- spSite[,2]

  ########################################################
  # Filter species (part 1)
  # We don't remove them from the dataset at this stage, but just identify the ones to remove later

  # set the minimum number of sites & records, as there are some species that were observed too rarely to model.
  # first check there are no impossible numbers
  if(minSite < 1) minSite <- 1
  if(minRecs < 1) minRecs <- 1

  if(minRecs > 1 | minSite > 1){
    # now identify species that occur on at least `minSite` sites and minRecs records
    sp2incl <- subset(spSumm, recs >= minRecs & sites >= minSite)$species

    nOrig <- nrow(spSumm)
    nExcl <- nOrig - length(sp2incl)

    if(verbose & nExcl > 0) {
      print(paste('Note:',nExcl,'species out of', nOrig, 'have been excluded because they occur on fewer than', minSite, 'sites or have fewer than', minRecs, 'records'))
      if(length(sp2incl) > 0) {
        print(paste('We proceed to modelling with', length(sp2incl), 'species'))
      } else {
        stop(paste0("There are no species with enough data to model"))
      }
    }
  } else {
    sp2incl <- spSumm$species
  }

  ##################### format data with one row per visit
  castDat <- dcast(inData, year + siteID + survey ~ "nsp",
                   value.var = "species", fun = length, fill = 0)

  # Julian date
  castDat$jday <- as.numeric(format(as.POSIXlt(castDat$survey, format = "%d%b%y"), "%j"))

  # renumber years starting at 1
  castDat$Year <- castDat$year # save a copy of the original
  castDat$year <- 1+ castDat$Year - min(castDat$Year)

  ########################################################

  # create metadata object
  md <- formatMetadata(inData)

  md$datastr$spSumm <- spSumm
  md$settings <- list(sp_modelled = length(sp2incl),
                      minSite = minSite,
                      minYr = min(castDat$Year),
                      maxYr = max(castDat$Year))

   ########################################################

  # extract the observations and populate the obsData list
  # this is the Nimble format. It is used in summariseData
  # it is not the format used by spOccupancy: we have to create that later.

  obsData <- list()

  obsMat <- acast(inData, year + survey + siteID ~ species,
                  value.var = "year", fun = function(x) max(x) > 0, fill=0)
  obsData$y = t(obsMat)[sp2incl,]

  ########################################################

  if(format == "Nimble"){
    dataConstants <- list(nsp = as.numeric(md$settings["sp_modelled"]),
                          nsite = md$dataPars$nSites,
                          nvisit = nrow(castDat),
                          nyear = md$dataPars$nYears,
                          year = castDat$year,
                          site = as.numeric(factor(castDat$siteID))
                          )
    # need some way to link site number with site identity

    if(!is.null(ListLen)){
      if(grepl("cont", ListLen, ignore.case = TRUE)){
        dataConstants$L <- castDat$nsp
      } else if(grepl("cat", ListLen, ignore.case = TRUE)){
        dataConstants$DT2 <- as.numeric(castDat$nsp %in% 2:3)
        dataConstants$DT3 <- as.numeric(castDat$nsp > 3)
      }
    }
    if(inclPhenology){dataConstants$JulDate <- castDat$jday}

    if(verbose){
      print(with(dataConstants, paste("Formatted data contains",
                                      nvisit, "visits to",
                                      nsite, "sites in",
                                      nyear, "years."
                                      )))
    }

  ########################################################

  } else if(format == "spOcc") {
    ##### for spOccupancy
    bd <- castDat %>%
      group_by(siteID, year) %>%
      mutate(Replicate = 1:n())

    # dimensions, site by year by reps
    nsite <- md$dataPars$nSites
    nyear <- md$dataPars$nYears
    nreps <- max(bd$Replicate)

    dataConstants <- list(
      # occupancy covariates are indexed by site and year
      occ.covs = list(
        Year = as.factor(matrix(rep(1:nyear, nsite),
                      nsite, nyear, byrow=TRUE)),
        Site = as.factor(matrix(rep(1:nsite, nyear),
                        nsite, nyear, byrow=FALSE))
      ),
      # detection covariates are indexed by site, year and replicate
      det.covs = list(
        year = as.factor(acast(bd, siteID ~ year ~ Replicate, fill = NA, value.var = "year"))
      )
    )

    if(!is.null(ListLen)){
      if(grepl("cont", ListLen, ignore.case = TRUE)){
        dataConstants$det.covs$logL = log(acast(bd, siteID ~ year ~ Replicate, fill = NA, value.var = "nsp"))
      } else if(grepl("cat", ListLen, ignore.case = TRUE)){
        bd$DT2 <- as.numeric(bd$nsp %in% 2:3)
        bd$DT3 <- as.numeric(bd$nsp > 3)
        dataConstants$det.covs$DT2 = as.factor(acast(bd, siteID ~ year ~ Replicate, fill = NA, value.var = "DT2"))
        dataConstants$det.covs$DT3 = as.factor(acast(bd, siteID ~ year ~ Replicate, fill = NA, value.var = "DT3"))
      }
    }
    if(inclPhenology){
      dataConstants$det.covs$jday = acast(bd, siteID ~ year ~ Replicate, fill = NA, value.var = "jday")
    }

    # create dataset with the same number of rows but with species in columns
    # this is what sparta would normally need. Nimble does not but spOccupancy does
    spObs <- acast(inData, year + siteID + survey ~ species,
                   value.var = "year", fun = function(x) length(x) > 0, fill = 0)

    if(nrow(spObs) != nrow(bd)) stop("Mismatching number of rows")

    temp <- melt(cbind(bd, spObs), id = 1:ncol(bd))

    #obsData <- list(y=acast(temp, siteID ~ year ~ Replicate ~ variable))
    obsData$scaff <- temp

  ########################################################

  } else
    stop("format not known")

  ########################################################

  return(list(dataConstants = dataConstants,
              obsData = obsData,
              md = md))
}
