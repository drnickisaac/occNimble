#' Format Metadata
#'
#' @param indat A dataset
#' @return A one line dataframe
#' @export


formatMetadata <- function(indat){
  list(
    dataPars = data.frame(nDates = length(unique(indat$survey)),
             nSites = length(unique(indat$siteID)),
             nYears = length(unique(indat$year)),
             sp_pool = length(unique(indat$species))
  ))
}
