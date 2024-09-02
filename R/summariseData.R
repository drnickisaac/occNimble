#' summariseData
#'
#' @details Summary statistics on the input data
#' @param obsData dataframe produced by formatData()
#' @param dataConstants dataframe produced by formatData()
#' @return summary statistics
#' @import reshape2
#' @export

summariseData <- function(obsData, dataConstants){

  # calculate for each visit whether the species was observed across 1-3 data types (depending on what exists)
  obs <- obsData$y

  # which sites are occupied?
  temp <- data.frame(cbind(site=dataConstants$site, year=dataConstants$year, t(obs)))
  temp <- melt(temp, id=1:2) %>%
    mutate(value = value>0) %>%
    group_by(site, year, variable) %>%
    summarise(occ = max(value))
  occMatrix <- acast(as.data.frame(temp), variable ~ site ~ year, value.var = "occ", fill = 0)

  occSpSite <- apply(occMatrix,c(1,2),max)

  # melt the occupancy per species:site combo
  mOSS <- melt(occSpSite)
  names(mOSS) <- c("species", "site", "occ")

  ############## Reporting Rate
  RR <- if(is.null(obsData$y1)) 0.000001 else rowMeans(obsData$y1/5)

  RR1 <- if(is.null(obsData$y1)) NULL else {
    temp <- melt(as.data.frame(cbind(site=dataConstants$site, t(obsData$y1/5))), id=1)
    names(temp)[c(2,3)] <- c("species", "TrapRate")
    as.data.frame(merge(mOSS, temp, all=TRUE) %>%
                    group_by(species) %>%
                    summarise(mean(TrapRate[occ == 1])))
  }

  ##############

  stats <- data.frame(
    species = dimnames(occMatrix)[[1]],
    naiveOcc = apply(occSpSite, 1, mean)
  )

  if(!is.null(obsData$y1)){
    stats$reportingRate <- RR # per pan trap
    stats$reportingRate_z1 <- as.numeric(RR1[,2]) # per visit on occupied sites
  }

  return(list(
    occMatrix = occMatrix,
    stats = stats))
}
