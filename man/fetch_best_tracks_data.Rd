% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load-parse-hurdat-data.R
\name{fetch_best_tracks_data}
\alias{fetch_best_tracks_data}
\title{Fetch Best Tracks data from NOAA (HURDAT2)}
\usage{
fetch_best_tracks_data(
  path,
  src = "https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2020-052921.txt"
)
}
\arguments{
\item{path}{character: the path to output downloaded data (*.csv)}

\item{src}{character: the path to most recent HURDAT2 Best tracks data. May need updating each year.}
}
\description{
This function is used to download, parse, and format the current tropical storm database
(HURDAT2) from NOAA. The \code{src} parameter may need to be updated each year as the path
to the input data is updated. Visit the NOAA website for full \href{https://www.nhc.noaa.gov/data/#hurdat}{data description}
and \href{https://www.nhc.noaa.gov/data/hurdat/hurdat2-format-nov2019.pdf}{metadata} for the HURDAT2 database
}
\examples{
# Download recent data from HURDAT2 (NOAA)
path = 'hurdat_data.csv'
fetch_best_tracks_data(path)

# load data for trackID AL142018 (Hurricane Michael)
track = load_hurdat_track(path, trackID = 'AL142018')

# Retrieve a US shapefile and project to UTM16
library(maptools)
data("wrld_simpl")
us = wrld_simpl['USA',]
us = spTransform(us, '+init=epsg:32616')

output_raster = hurrecon_run(track, land=us, max_rad_km = 100, res_m = 500, max_interp_dist_km = 1)
raster::plot(output_raster)

}
