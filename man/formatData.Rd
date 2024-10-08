% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatData.R
\name{formatData}
\alias{formatData}
\title{formatData}
\usage{
formatData(
  inData,
  format = "Nimble",
  ListLen = NULL,
  inclPhenology = TRUE,
  minYrPerSite = 2,
  minSite = 1,
  minRecs = 1,
  verbose = TRUE
)
}
\arguments{
\item{inData}{A dataset with columns "species", "siteID", "survey" and "year"}

\item{format}{Either "Nimble" (default) or "spOcc"}

\item{ListLen}{should be included as continuous ("Cont"), categorical ("Cat") or excluded (NULL)}

\item{inclPhenology}{should the model account for seasonal variation?}

\item{minYrPerSite}{the minimum number of years with data for a site to be included (defaults to 2, as in sparta)}

\item{minSite}{the threshold minimum number of sites for a species to be considered for modelling}

\item{minRecs}{the threshold minimum number of records for a species to be considered (after applying \code{minSite})}
}
\value{
list of two data frames
}
\description{
formatData
}
\details{
Formats data ready for Nimble
}
\examples{
\dontrun{

set.seed(123)

# Create data
n <- 15000 #size of dataset
nyear <- 20 # number of years in data
nSamples <- 100 # set number of dates
nSites <- 50 # set number of sites

# Create somes dates
first <- as.Date(strptime("2010/01/01", format="\%Y/\%m/\%d"))
last <- as.Date(strptime(paste(2010+(nyear-1),"/12/31", sep=''), format="\%Y/\%m/\%d"))
dt <- last-first
rDates <- first + (runif(nSamples)*dt)

# taxa are set as random letters
taxa <- sample(letters, size = n, TRUE)

# taxa are set as random letters: weight so that some species are more commonly observed
taxa <- sample(letters, size = n, TRUE, prob = seq(from = 0.01, to = 0.4, length.out = 26))

# sites are visited randomly
site <- sample(paste('site_', 1:nSites, sep=''), size = n, TRUE)

# Format the data
inData <- data.frame(species = taxa,
                     siteID = site,
                     survey = survey,
                     year = as.numeric(format(as.Date(survey, format="\%d/\%m/\%Y"),"\%Y")))

# prepare the data for the model (includes removing species found on few sites)
formattedData <- formatData(inData)
}

}
