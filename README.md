# Model Hurricane Wind Fields Using HURRECON

<img src=img/katrina.png width=200 align=right>

The purpose of this package is to implement the HURRECON model described in Boose et al. 2001. The model uses observations of hurricane location, maximum wind speeds, and radius at which particular wind speeds occur from the NOAA HURDAT2 database. The model provides functions to gap fill between 6-hour observations, fit wind profiles from HURDAT observations, and generate rasters of maximum wind speed and direction of maximum wind speed continuously along the hurricane track. 

# Install `hurrecon`

If you haven't installed packages from Github before or used the `devtools` package, you'll need to first ensure that Rtools is properly installed. It is installed separately from typical packages.

* Install the correct version of Rtools using this link: https://cran.r-project.org/bin/windows/Rtools/rtools40.html
* Be sure to follow the directions at the bottom regarding *Putting Rtools on the PATH*. This is critical for a proper installation. You'll need to create a file in your Documents folder named '.Renviron' that contains the following line: `PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"` and restart R.

Next, you should be able to get the latest released version of `hurrecon` from github

```
install.packages('devtools')
devtools::install_github('jbcannon/hurrecon')
```
# Example Usage

`hurrecon` can be used both to download data from the NOAA best tracks database, and process a hurricane track. Below is example usage.

```
# Download recent data from HURDAT2 (NOAA)
library(hurrecon)
path = 'hurdat_data.csv'  #path to output downloaded data
fetch_best_tracks_data(path)

# load data for trackID AL142018 (Hurricane Michael)
track = load_hurdat_track(path, trackID = 'AL142018')

# A shapefile representing "land" is needed as an input to reduce surface wind speeds over land
# Retrieve a US shapefile and project to UTM16
library(maptools)
library(sp)
data("wrld_simpl")
us = wrld_simpl['USA',]
us = spTransform(us, '+init=epsg:32616')
sp::plot(track)
sp::plot(us, add=TRUE)
 
output_raster = hurrecon_run(track, land=us, max_rad_km = 100, res_m = 500, max_interp_dist_km = 1)
library(raster)
plot(output_raster$Vs)
```
