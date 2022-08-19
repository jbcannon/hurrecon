#' @export
get_Ha_Hv = function(trk){
  trk$dt = paste0(trk$date, trk$time)
  trk$dt = as.POSIXct(trk$dt, format= '%Y%m%d%H%M', tz = 'GMT')
  trk$Vh = 0
  trk$Ah = 0
  for(i in 2:nrow(trk)){
    # calc distance between consecutive points
    tmp = trk[(i-1):i, ]
    dist_m = max(sf::st_distance(tmp, by_element = FALSE))
    #Horizontal angle
    coords = sf::st_coordinates(tmp)
    HA = apply(coords, 2, diff)
    HA = atan2(HA[1], HA[2])
    trk$Ah[i] = HA
    # Horizontal velocity
    time_diff = as.numeric(diff(tmp$dt)*(60*60))
    Vh = dist_m / time_diff
    trk$Vh[i] = Vh
  }
  return(sf::st_drop_geometry(trk[, c('Ah', 'Vh')]))
}

#' @export
hurrecon = function(track_densif, res_m, max_radius_km){
  ext = suppressWarnings(cumulative_extent(track_densif, res_m = res_m, max_radius_km = max_radius_km))
  Vs_out = ext
  Vs_out[] = 0
  track_densif = subset(track_densif, max_speed >= 34)
  if(nrow(track_densif) < 1) return(list(Vs=Vs_out, D=D_out)) else {
    len = nrow(track_densif)
    for(i in 1:len){
      x = track_densif[i,]
      if(x$`50kt_ne` == 0 & x$`50kt_se` == 0 & x$`50kt_sw` == 0 & x$`50kt_nw` == 0) next

      Vs_tmp = suppressWarnings(Vs(track_densif[i,], max_radius_km = max_radius_km, template=ext))
      v = terra::extend(Vs_tmp$Vs, ext)
      v[is.na(v)]= 0
      ismax = v > Vs_out
      Vs_out = (Vs_out * !ismax) + (v * ismax)

      cat(i, ' of ', len, ' (', round(i/len*100,1), '% complete)\n', sep = '')
    }
    return(Vs_out)
  }
}

getfun<-function(x) {
  if(length(grep("::", x))>0) {
    parts<-strsplit(x, "::")[[1]]
    getExportedValue(parts[1], parts[2])
  } else {
    x
  }
}

find_densification = function(track_pts, max_dist_km){
  # Get maximum distance between points
  d = sf::st_distance(track_pts) / 1000

  #get distances between adjacent points
  if(nrow(d)>2) {x = diag(d[-1,-ncol(d)])} else {x = d[-1,-ncol(d)]}
  densification_factor = ceiling(x/max_dist_km)
  return(densification_factor)
}

densify = function(track_pts, factor, land){
  constant_vars = c('track_id', 'track_name', 'date', 'time', 'record', 'status', 'lat', 'lon')
  interpolate_vars = c('utmx', 'utmy', 'max_speed', 'min_press', grep('kt_', colnames(track_pts), value=TRUE), 'Ah', 'Vh')
  out_pts = list()

  for(i in 1:(nrow(track_pts)-1)){
    local_factor = as.numeric(factor[i])
    pair = sf::st_drop_geometry(track_pts[i:(i+1),])
    out_df = list()
    for(v in constant_vars){
      x = data.frame(rep(pair[1,v], local_factor+1))
      x=t(x)
      colnames(x) = v
      out_df[[length(out_df)+1]] = x
    }
    for(v in interpolate_vars){
      if(any(is.na(pair[,v]))) x= tibble::tibble(rep(NA, local_factor+1)) else {
        x = data.frame(approx(pair[,v], n = local_factor+1)$y)
      }
      colnames(x) = v
      out_df[[length(out_df)+1]] = x
    }

    out_df = do.call(getfun('dplyr::bind_cols'), out_df)
    out_pts[[length(out_pts)+1]] = out_df
  }

  # merge all dataframes together
  if(length(out_pts)==1) out_pts = out_pts[[1]] else {
    for(i in 2:length(out_pts)){
      out_pts[[i]] = out_pts[[i]][-1,]
    }
    out_pts = do.call(getfun('dplyr::bind_rows'), out_pts)
  }

  #convert to sf
  coords = out_pts[, c('utmx', 'utmy')]
  out_pts[, c('spatialdatX', 'spatialdatY')] = coords
  out_pts = sf::st_as_sf(out_pts, coords = c('spatialdatX', 'spatialdatY'), crs = sf::st_crs(track_pts))

  on_land = sf::st_intersects(out_pts, sf::st_transform(land, sf::st_crs(track_pts)), sparse=FALSE)
  out_pts$cov = ifelse(on_land, 'L', 'W')

  return(out_pts)
}

# Extract radius of each speed threshold and get max or mean, output in metric
wind_profile_data = function(t_spdf, fun = 'max'){
  #convert values to numeric
  for(i in grep('kt_|max_speed', colnames(t_spdf), value = TRUE)){
    t_spdf[, i] = as.numeric(as.character(t_spdf[, i]))
  }
  Rm = list()
  for(k in c(34, 50, 64)){
    Rm[[length(Rm)+1]] = apply(t_spdf[grep(k, colnames(t_spdf))], 1, FUN = fun, na.rm = TRUE)
  }; Rm = do.call('cbind', Rm);
  Rm = Rm*1.852 # convert nautical miles to km
  max_speed = t_spdf$max_speed * 0.514444 # convert knots to m/s
  Rm = cbind(max_speed, Rm)
  colnames(Rm) = c('Vm_ms', paste0('Rkm_', c(34,50,64), 'kt'))

  return(Rm)
}

# Function to fit B coefficient from relative wind profile from Boose et al. 2004
get_B_coeff = function(profile){
  max_speed = profile[1]
  profile = data.frame(V = c(34, 50, 64)/max_speed , R = c(profile[2:4]))
  if(profile[2,2] == 0) {return(data.frame(B=NA, Rm_km = NA, Vm_ms = max_speed))}
  if(profile[3,2] == 0) profile = profile[1:2,]
  profile = rbind(data.frame(V=0, R=0.1), profile)
  model_fit_found = FALSE
  j = 0
  pars = list(B=0.5, Rm=10)
  while(!model_fit_found & j <= 5) {
    mod = try(nls(V ~ wind_profile_fxn(R, Rm, B),
                  data = profile,
                  start = list(B = pars$B, Rm = pars$Rm),
                  algorithm = 'port', lower = c(0.1,1), upper = c(4,200),
                  control = nls.control(maxiter=1000)), silent = TRUE)
    if(!'try-error' %in% class(mod)) {model_fit_found = TRUE} else {
      if(j == 0) pars = list(B=0.5, Rm=1) # Try a few different starting points
      if(j == 1) pars = list(B=1, Rm=10)
      if(j == 2) pars = list(B=1, Rm=1)
      if(j == 3) pars = list(B=1, Rm=20)
      if(j == 4) pars = list(B=1, Rm=30)
      if(j == 5) pars = list(B=0.5, Rm=20)
      j = j + 1
    }
  }
  if(model_fit_found) {coeff = coef(mod)} else {coeff = c(NA, NA)}
  x = data.frame(B=coeff[1], Rm_km = coeff[2], Vm_ms = max_speed*0.51444)
  return(list(coeff=x, profile=profile))
}

#'  Return wind profile data from HURDAT2 observation
#'
#' This function takes a single HURDAT2 observation retreived from `load_hurdat_track()`
#' and returns a wind profileis used to download, parse, and format the current tropical storm database
#' @examples
#' # Download recent data from HURDAT2 (NOAA)
#` hurdat_database = 'path/to/save/database.txt'
#' fetch_best_tracks_data(hurdat_database)
#' track = load_hurdat_track('AL142018') # Load data from hurricane Micahel
#' observation = track[17,] #select a single observation to analyze
#' get_wind_profiles(observation)
#' @param obs SpatialPointsDataFrame: a single hurricane observation represented by
#' a single row of an object returned from `load_hurdat_track()`
#' @export
get_wind_profiles = function(obs){
  #convert values to numeric then to metric
  obs = sf::st_drop_geometry(obs)
  for(i in grep('kt_|max_speed', colnames(obs), value = TRUE)){
    obs[,i] = as.numeric(as.character(obs[,i]))
  }
  obs$max_speed_ms = obs$max_speed * 0.514444 # convert knots to m/s
  for(i in grep('kt_', colnames(obs), value = TRUE)) obs[,i] = obs[,i] *1.852   # convert nautical mi to km

  # Get profile for each direction
  Rm_ne = unlist(obs[, c('max_speed', grep('_ne', colnames(obs), value=TRUE))])
  Rm_se = unlist(obs[, c('max_speed', grep('_se', colnames(obs), value=TRUE))])
  Rm_nw = unlist(obs[, c('max_speed', grep('_nw', colnames(obs), value=TRUE))])
  Rm_sw = unlist(obs[, c('max_speed', grep('_sw', colnames(obs), value=TRUE))])
  x = rbind(Rm_ne, Rm_se, Rm_nw, Rm_sw)

  # If there are 0 readings in some direction for 50kt, then use average readings
  if(any(x[,3] = 0)) {x[,2] = mean(x[,2]); x[,3] = mean(x[,3])}

  # if there are any negative readings in the third column, zero them out.
  x[, x[,3] < 0] = 0

  ne = get_B_coeff(x[1,])
  se = get_B_coeff(x[2,])
  nw = get_B_coeff(x[3,])
  sw = get_B_coeff(x[4,])

  ne$profile$dir='ne'
  se$profile$dir='se'
  nw$profile$dir='nw'
  sw$profile$dir='sw'

  x = rbind(ne$coeff,se$coeff,nw$coeff,sw$coeff)
  x= cbind(data.frame(dir=c('ne','se','nw','sw')), x)
  x = tibble::as_tibble(x)
  out2 = rbind(ne$profile, se$profile, nw$profile, sw$profile)
  out2 = tibble::as_tibble(out2)
  return(list(summary=x, profiles=out2))
}

#' @export
Vs = function(obs, template, max_radius_km){
  # Create blank raster centered over eye
  eye = obs
  eye_ext = sf::st_buffer(eye, dist = max_radius_km*1000)
  ex = template
  ex = terra::crop(ex, terra::ext(eye_ext))
  ex[] = 0

  #calculate radius and direction raster
  R = terra::distance(ex, terra::vect(eye))
  pts = tibble::as_tibble(terra::crds(R, df=TRUE))
  pts$x = pts$x-eye$utmx
  pts$y = pts$y-eye$utmy
  A = ex
  A[] = atan2(pts$y, pts$x)
  A = -A +2.5*pi
  A = terra::ifel(A>=2*pi, A-2*pi, A)
  T = ex
  T[]=eye$Ah
  T = T-A
  T = terra::ifel(T >= 2*pi, T - 2*pi, T)
  T = terra::ifel(T < 0, T + 2*pi, T)
  wp = get_wind_profiles(obs)$summary

  F = as.numeric(ifelse(obs$cov=='W', 1, 0.8))

  Vs = list()
  for(d in c('ne', 'se', 'sw', 'nw')){
    tmp_wp = subset(wp, dir==d)
    Vm = tmp_wp$Vm_ms
    Vh = eye$Vh
    B = tmp_wp$B
    Rm = tmp_wp$Rm_km * 1000
    S = 1 #per Boose
    Vs[[d]] = F * (Vm - S*(1-sin(T))*Vh/2) * sqrt((Rm/R)^B*exp(1-(Rm/R)^B))
  }

  weights = list()
  for(d in 0:3){
    dir = pi/4 + d*(pi/2)
    diff = acos(cos(abs(dir-A)))
    diff[diff>pi/2] = NA
    diff = 1 - diff/max(diff[], na.rm=TRUE)
    dir = c('ne', 'se', 'sw', 'nw')
    weights[[dir[d+1]]] = diff
  }

  Vs_out = list()
  for(i in 1:4){
    x = weights[[i]] * Vs[[i]]
    x[is.na(x)] = 0
    Vs_out[[i]] =x
  };
  Vs_out = Vs_out[[1]] + Vs_out[[2]] + Vs_out[[3]] + Vs_out[[4]]
  Vs = Vs_out

  #Calculate wind direction
  I = as.numeric(ifelse(eye$cov=='W', 20, 40)*pi/180)
  Az = A - pi
  Az = terra::ifel(Az < 0, Az + 2*pi, Az)
  D = Az - pi - I
  D = terra::ifel(D < 0, D + 2*pi, D)
  return(list(Vs=Vs,D=D))
}

#' @export
cumulative_extent = function(track, res_m, max_radius_km){
  track = sf::st_buffer(track, dist=max_radius_km*1000)
  trk_ext = terra::ext(track)
  trk_ext[c(1,3)] = floor(trk_ext[c(1,3)]/res_m)*res_m
  trk_ext[c(2,4)] = ceiling(trk_ext[c(2,4)]/res_m)*res_m
  #ex =
  ex = terra::rast(trk_ext, resolution = res_m,
                crs = sf::st_crs(track, parameters=TRUE)$srid)
  return(ex)
}

size_pred = function(dat) {
  if(!exists('radius_models')) {stop('must first run \`data("radius_models\")\`')}
  if(!exists('tab_asymmetry')) {stop('must first run \`data("radius_models\")\`')}

  cat('Predicting hurricane size. See Cannon et al. for details\n')

  check_pred = function(v) {
    df = sf::st_drop_geometry(dat)
    x = dat$max_speed >= v
    y = apply(sf::st_drop_geometry(dat[,paste0(v, 'kt_', c('ne', 'se', 'sw', 'nw'))]),1,mean) == -999
    return(x & y)
  }

  press_known = ifelse(dat$min_press > -999, 'p', 'np')
  needs_pred = sapply(34, function(v) check_pred(v))

  df = dat[, c('max_speed', 'min_press', 'lat')]
  df = sf::st_drop_geometry(df)
  for(i in 1:3) df[, i] = as.numeric(unlist(df[,i]))
  #predict inner-most windspeed
  r34 = sapply(1:nrow(df), function(i) ifelse(df$max_speed[i]>34, predict(radius_models[[1]][[press_known[i]]], df[i,]),0))
  r50 = sapply(1:nrow(df), function(i) ifelse(df$max_speed[i]>50, predict(radius_models[[2]][[press_known[i]]], df[i,]),0))
  r64 = sapply(1:nrow(df), function(i) ifelse(df$max_speed[i]>64, predict(radius_models[[3]][[press_known[i]]], df[i,]),0))
  size = list(r34=r34, r50=r50, r64=r64)

  #apply assymtry factors
  assym = tab_asymmetry$ratio_mn
  out = list()
  for(i in 1:3){
      x = size[[i]]
      x = do.call(rbind, lapply(x, function(i) i*assym))
      x = tibble::as_tibble(x)
      colnams = paste0(substr(names(size)[i],2,3), 'kt_', tab_asymmetry$direction)
      colnames(x) = colnams
      out[[length(out)+1]] = x
  }
  out = tibble::as_tibble(do.call(cbind, out))

  x = expand.grid(x=c(34,50,64), y=c('ne','nw','se','sw'));x= x[order(x$x),]
  col_replace = with(x, paste0(x,'kt_',y))
  input = sf::st_drop_geometry(dplyr::select(dat, col_replace))

  dat = dplyr::select(dat, -col_replace)
  for(i in which(needs_pred)) input[i,] = out[i,]
  dat = dplyr::bind_cols(dat, input)

  return(dat)
}

#' Function to fit HURRECON wind profile as a function
#' This function will map relative wind profile based on 2 or 3 radius with known speed
#' @examples
#' radius = c(seq(0,1,by=0.01), seq(2,150,1))
#' wind_profile_fxn(radius, Rm = 4, B=1.3)
#' @param R numeric: vector of radii at which to return wind profile
#' @param Rm numeric: Radius at which maximum wind speed occurs
#' @param B numeric: Decay coefficient for wind profile. See Boose et al. 2004
#' @export
wind_profile_fxn = function(R, Rm, B) sqrt((Rm/R)^B *exp(1- (Rm/R)^B)) # From Boose et al. 2004

