#'  Fetch Best Tracks data from NOAA (HURDAT2)
#' 
#' This function is used to download, parse, and format the current tropical storm database
#' (HURDAT2) from NOAA. The `src` parameter may need to be updated each year as the path
#' to the input data is updated. Visit the NOAA website for full [data description](https://www.nhc.noaa.gov/data/#hurdat)
#' and [metadata](https://www.nhc.noaa.gov/data/hurdat/hurdat2-format-nov2019.pdf) for the HURDAT2 database
#' @examples
#' # Download recent data from HURDAT2 (NOAA)
#' library(terra)
#' library(hurrecon)
#' data("geographic")
#' 
#' path = 'hurdat_data.csv'
#' fetch_best_tracks_data(path)
#' 
#' # load data for trackID AL142018 (Hurricane Michael)
#' track = load_hurdat_track(path, trackID = 'AL142018')
#'  
#' # Quick one for example
#' output = hurrecon_run(track, land=land, max_rad_km = 100, res_m = 500, max_interp_dist_km = 50)
#' plot(land)
#' plot(output, add = TRUE)
#' plot(land)
#' 
#' @param path character: the path to output downloaded data (*.csv)
#' @param src character: the path to most recent HURDAT2 Best tracks data. May need updating each year.
#' @export
fetch_best_tracks_data = function(path, src = 'https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2020-052921.txt') {
  if(!dir.exists(dirname(path))) {
    stop(paste0('directory', dirname(path), 'does not exist'))
  }
  cat('downloading best track data from\n', src)
  data = readLines(src)
  cat(' parsing track data')
  parsed_data = list()
  i=1
  prog_yr = 1850
  while(i <= length(data)) {
    entry = data[i]
    header = trimws(strsplit(entry, ',')[[1]])
    track_id = header[1]
    track_name = header[2]
    rows_to_follow = as.numeric(header[3])
    body = data[(i+1):(i+rows_to_follow)] 
    body = strsplit(body, ',')
    body = lapply(body, trimws)
    body = lapply(body, function(x) x[1:20]) # Fix database error with extra comma
    body = do.call(rbind,body)
    colnames(body) = c('date','time',
                       'record', 'status',
                       'lat', 'lon',
                       'max_speed', 'min_press',
                       '34kt_ne', '34kt_se', '34kt_sw', '34kt_nw',
                       '50kt_ne', '50kt_se', '50kt_sw', '50kt_nw',
                       '64kt_ne', '64kt_se', '64kt_sw', '64kt_nw')
    body = tibble::as_tibble(body)
    body$lat = ifelse(grepl('N', body$lat), gsub('N', '', body$lat), paste0('-',gsub('S', '', body$lat)))
    body$lon = ifelse(grepl('E', body$lon), gsub('E', '', body$lon), paste0('-',gsub('W', '', body$lon)))
    hd = tibble::tibble(track_id = track_id, track_name = track_name)
    track_info = dplyr::bind_cols(hd,body)
    parsed_data[[length(parsed_data)+1]] = track_info
    curr_yr = substr(track_info$date[1],1,4)
    i = i + rows_to_follow + 1
    if(curr_yr != prog_yr) {
      cat('year', prog_yr, 'complete\n')
      prog_yr = curr_yr
    }
  }
  cat('writing best tracks to disk\n', path,'\n')
  out_data = do.call(rbind, parsed_data)
  readr::write_csv(out_data, path)
  return(out_data)
}

#' Load Full Track and Return as POINT simple feature
#' 
#' This function returns an attributed feature containing tropical cyclone positions
#' and attributes from HURDAT2 database.
#' @examples
#' # Download recent data from HURDAT2 (NOAA)
#' library(terra)
#' library(hurrecon)
#' data("geographic")
#' 
#' path = 'hurdat_data.csv'
#' fetch_best_tracks_data(path)
#' 
#' # load data for trackID AL142018 (Hurricane Michael)
#' track = load_hurdat_track(path, trackID = 'AL142018')
#'  
#' # Quick one for example
#' output = hurrecon_run(track, land=land, max_rad_km = 100, res_m = 500, max_interp_dist_km = 50)
#' plot(land)
#' plot(output, add = TRUE)
#' plot(land)
#' 
#' @param path character: path to parsed HURDAT2 database downloaded using `fetch_best_tracks_data()`.
#' @param trackID string: trackID from HURDAT database (e.g., AL122018 indicates Hurricane Michael, the twelfth Atlantic hurricane of the 2018 season)
#' @param proj output CRS for feature. Must be a projection coordinate system (e.g., UTM) for distance-based calculations
#' @export

load_hurdat_track = function(path, trackID, proj=32616) {
  db = suppressMessages(readr::read_csv(path))
  
  if(length(trackID) > 1) stop('only one track_id allowed at a time')
  valid = all(c('track_id', 'record', 'date', 'max_speed', 'status', 'lat') %in% colnames(db)) 
  
  if(!valid) stop('path does not appear to be valid HURDAT2 dataset. Use fetch_best_tracks_data() to download valid dataset')
  missing = !trackID %in% db$track_id
  
  if(!trackID %in% db$track_id) stop('track', trackID, '\nmissing or invalid')
  
  crs_valid = !is.na(sf::st_crs(proj)$input)
  if(!crs_valid) stop('proj = ', as.character(proj), ' invalid. Use EPSG code (e.g., 32616)')
  
  db = dplyr::filter(db,  track_id %in% trackID)
  db[, c('lon2', 'lat2')] = db[, c('lon', 'lat')]
  output = sf::st_as_sf(db, coords = c('lon2', 'lat2'), crs=4326)
  output = sf::st_transform(output, crs = proj)
  output[, c('utmx', 'utmy')] = sf::st_coordinates(output)
  
  # Check if there needs to be some gap filling
  yr = max(as.numeric(stringr::str_sub(output$date,1,4)))
  if(yr < 2004) {
    cat('pre-2004 track requires gap filling\n')
    output = size_pred(output)
    } else {
      # Do any rows need size prediction?
      needs_pred = (apply(sf::st_drop_geometry(dplyr::select(output, dplyr::contains('34'))),1,function(i) all(i==-999)) & output$max_speed>34)
      if(any(needs_pred)) {
        cat('some size info missing, gap filling', sum(needs_pred), 'rows\n')
        output = size_pred(output)
      }
    }
  return(output)
}
