#' Run HURRECON model for single Hurricane Track
#' 
#' This function takes a single hurricane track from `load_hurdat_track` and 
#' estimates a wind field producing the maximum sustained wind (`Vs` in ms^-1)
#' and direction of maximum wind speed (`d` in radians)
#' @examples
#' # Download recent data from HURDAT2 (NOAA)
#' path = 'hurdat_data.csv'
#' fetch_best_tracks_data(path)
#' 
#' # load data for trackID AL142018 (Hurricane Michael)
#' track = load_hurdat_track(path, trackID = 'AL142018')
#' 
#' # Retrieve a US shapefile and project to UTM16
#' library(maptools)
#' data("wrld_simpl")
#' us = wrld_simpl['USA',]
#' us = spTransform(us, '+init=epsg:32616')
#' 
#' output_raster = hurrecon_run(track, land=us, max_rad_km = 100, res_m = 500, max_interp_dist_km = 1)
#' raster::plot(output_raster)
#' 
#' @param trk SpatialPointsDataFrame of single track data returned from `load_hurdat_track` 
#' @param max_rad_km integer: maximum distance away from eye to model windspeed
#' @param res_m integer: spatial resolution of output raster
#' @param max_interp_dist_km integer: with more than one observation, `max_interp_dist_km`
#' is the spatial resolution of the linear interpolation between concurrent hurrecon observation.
#' When `max_interp_dist_km` = 1, an interpolated observation is placed every 1 km between two
#' observations further apart than 1 km. Values < 1 km are not recommend to avoid gaps between
#' wind fields
#' @param proj string: default projecgtion for analyses. Recommend keeping at UTM16 (epsg:32616)
#' @param land SpatialPolygons: This is spatial layer representing land. This is used to identify when the hurricane is on land to incorporate appropriate adjustment factors
#' @param aoi area of interest, track locations outside of aoi are ignored. This can be used
#' to limit results for a specific spatial region (county, state, country, etc.). Must be
#' in same projection as `trk`.
#' @export
hurrecon_run = function(trk, max_rad_km = 100, res_m = 500, max_interp_dist_km = 1, proj = '32616',aoi = NULL) {
  # check that track is valid
  val = all(class(trk) == c('sf', 'tbl_df', 'tbl', 'data.frame'))
  if(!val) stop('trk must be a valid sf object from load_hurdat_track')
  
  val = all(c('track_id', 'record', 'time', 'lon', '34kt_se') %in% colnames(trk))
  if(!val) stop('trk does not appear to contain valid columns. Load track with load_hurdat_track')
  
  if(!'land' %in% ls()) {
    data('geographic')
  }
  
  if(!'radius_models' %in% ls()) {
    data('radius_models')
  }
  
  # Check if there are missing values in size prediction
  needs_pred = apply(dplyr::select(sf::st_drop_geometry(trk), contains('34') | contains('50') | contains('64')), 1, function(x) any(x<0)) & trk$max_speed >= 34
  if(any(needs_pred)) stop('Missing size data in some tracks, us gap fill size information containng value -999')
  vars = grep('kt_', colnames(trk), value=TRUE)
  for(v in vars) trk[dplyr::pull(trk, v) == -999,v] = 0

  #Add movement columns
  movement = get_Ha_Hv(trk)
  trk = dplyr::bind_cols(trk, movement)
  
  dens_fact = find_densification(trk, max_interp_dist_km);
  track_densif = densify(trk, dens_fact, land)
  if(is.null(aoi)) {
    track_densif_aoi = track_densif
    } else {
      track_densif_aoi = sf::st_intersection(track_densif, sf::st_transform(aoi, sf::st_crs(trk)))
    }
  hurrecon_out = hurrecon(track_densif_aoi, max_radius_km = max_rad_km, res_m = res_m)
  return(hurrecon_out)
}