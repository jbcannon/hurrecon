get_Ha_Hv = function(t_spdf){
  t_spdf$dt = paste0(t_spdf$date, t_spdf$time)
  t_spdf$dt = as.POSIXct(t_spdf$dt, format= '%Y%m%d%H%M', tz = 'GMT')
  t_spdf$Vh = 0
  t_spdf$Ah = 0
  for(i in 2:length(t_spdf)){
    # calc distance between consecutive points
    tmp = t_spdf[(i-1):i, ]
    dist_m = max(raster::pointDistance(tmp))
    #Horizontal angle
    coords = sp::coordinates(tmp)
    HA = apply(coords, 2, diff)
    HA = atan2(HA[1], HA[2])
    t_spdf$Ah[i] = HA
    # Horizontal velocity
    tmp = tmp@data
    time_diff = as.numeric(diff(tmp$dt)*(60*60))
    Vh = dist_m / time_diff
    t_spdf$Vh[i] = Vh
  }
  return(t_spdf@data[, c('Ah', 'Vh')])
}

hurrecon = function(track_densif, res_m, max_radius_km){
  ext = suppressWarnings(cumulative_extent(track_densif, res_m = res_m, max_radius_km = max_radius_km))
  Vs_out = ext
  D_out = ext
  track_densif = subset(track_densif, max_speed >= 34)
  if(nrow(track_densif) < 1) return(list(Vs=Vs_out, D=D_out)) else {
    len = nrow(track_densif)
    Vs_out[]=0
    for(i in 1:len){
      x = track_densif[i,]
      if(x$`50kt_ne` == 0 & x$`50kt_se` == 0 & x$`50kt_sw` == 0 & x$`50kt_nw` == 0) next
      Vs_tmp = suppressWarnings(Vs(track_densif[i,], proj=proj, res_m = res_m, max_radius_km = max_radius_km, template=ext))
      v = raster::extend(Vs_tmp$Vs, ext)
      d = raster::extend(Vs_tmp$D, ext)
      ismax = v > Vs_out
      #update Vout and D with new max
      Vs_out[ismax] = v[ismax]
      D_out[ismax] = d[ismax]
      cat(i, ' of ', len, ' (', round(i/len*100,1), '% complete)\n', sep = '')
    }
    return(list(Vs=Vs_out, D=D_out))  
  }
  
}


find_densification = function(track_pts, max_dist_km){
  # Get maximum distance between points
  d = as.matrix(dist(sp::coordinates(track_pts)))
  d = d/1000
  d = diag(d[-1,-ncol(d)]) #get distances between adjacent points
  densification_factor = ceiling(d/max_dist_km)
  return(densification_factor)
}

densify = function(track_pts, factor, land, proj){
  constant_vars = c('track_id', 'track_name', 'date', 'time', 'record', 'status', 'lat', 'lon')
  interpolate_vars = c('utmx', 'utmy', 'max_speed', 'min_press', grep('kt_', colnames(track_pts@data), value=TRUE), 'Ah', 'Vh')
  out_pts = list()
  
  for(i in 1:(length(track_pts)-1)){
    local_factor = factor[i]
    pair = track_pts[i:(i+1),]
    out_df = list()
    for(v in constant_vars){
      x = data.frame(rep(pair@data[1,v], local_factor+1))
      colnames(x) = v
      out_df[[length(out_df)+1]] = x
    }
    for(v in interpolate_vars){
      if(any(is.na(pair@data[,v]))) x= data.frame(rep(NA, local_factor+1)) else {
        x = data.frame(approx(pair@data[,v], n = local_factor+1)$y)
      }
      colnames(x) = v
      out_df[[length(out_df)+1]] = x
    }
    out_df = do.call('cbind', out_df)
    out_pts[[length(out_pts)+1]] = out_df
  }
  # merge all dataframes together
  for(i in 2:length(out_pts)){
    out_pts[[i]] = out_pts[[i]][-1,]
  }
  out_pts = do.call(rbind, out_pts)
  #out_pts = out_pts[, colnames(track_pts@data)]
  
  #convert to spdf
  coords = out_pts[, c('utmx', 'utmy')]
  coords = sp::SpatialPoints(coords)
  out_pts = sp::SpatialPointsDataFrame(coords, out_pts)
  sp::proj4string(out_pts) = proj
  
  out_pts$cov = 'W'
  land_pts = out_pts[land,]
  out_pts[row.names(land_pts), 'cov'] = 'L'
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
  mod = try(nls(V ~ wind_profile_fxn(R, Rm, B),
                data = profile,
                start = list(B = 1.3, Rm = 40),
                algorithm = 'port', lower = c(1.0,1), upper = c(1.7,200),
                control = nls.control(maxiter=1000)), silent = TRUE)
  if('try-error' %in% class(mod)) {
    coeff = c(NA, NA)
  } else {
    coeff = coef(mod)  
  }
  x = data.frame(B=coeff[1], Rm_km = coeff[2], Vm_ms = max_speed*0.51444)
  return(list(coeff=x, profile=profile))
}

# Feed a single observation and get wind profiles in each quadrant
get_wind_profiles = function(obs){
  #convert values to numeric then to metric
  for(i in grep('kt_|max_speed', colnames(obs@data), value = TRUE)){
    obs@data[,i] = as.numeric(as.character(obs@data[,i]))
  } 
  obs$max_speed_ms = obs$max_speed * 0.514444 # convert knots to m/s
  for(i in grep('kt_', colnames(obs@data), value = TRUE)) obs@data[,i] = obs@data[,i] *1.852   # convert nautical mi to km
  
  # Get profile for each direction
  Rm_ne = unlist(obs[, c('max_speed', grep('_ne', colnames(obs@data), value=TRUE))]@data)
  Rm_se = unlist(obs[, c('max_speed', grep('_se', colnames(obs@data), value=TRUE))]@data)
  Rm_nw = unlist(obs[, c('max_speed', grep('_nw', colnames(obs@data), value=TRUE))]@data)
  Rm_sw = unlist(obs[, c('max_speed', grep('_sw', colnames(obs@data), value=TRUE))]@data)
  x = rbind(Rm_ne, Rm_se, Rm_nw, Rm_sw)
  
  # If there are 0 readings in some direction for 50kt, then use average readings
  if(any(x[,3] == 0)) {x[,2] = mean(x[,2]); x[,3] = mean(x[,3])}
  
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
  
  out2 = rbind(ne$profile, se$profile, nw$profile, sw$profile)
  
  return(list(summary=x, profiles=out2))
}

Vs = function(obs, proj, template, max_radius_km, res_m){
  # Create blank raster centered over eye
  eye = obs
  eye_ext = rgeos::gBuffer(eye, width = max_radius_km*1000)
  ex = template
  ex = raster::crop(ex, raster::extent(eye_ext))
  ex[] = 0
  
  #calculate radius and direction raster
  R = raster::distanceFromPoints(ex,eye)
  pts = as.data.frame(raster::rasterToPoints(ex))
  pts$x = pts$x-eye$utmx
  pts$y = pts$y-eye$utmy
  pts = atan2(pts[,'y'], pts[,'x'])
  A = ex
  A[] = pts
  A = -A +2.5*pi
  A[A >= 2*pi] = A[A>=2*pi] - 2*pi
  T = ex
  T[]=eye$Ah
  T = T-A
  T[T >= 2*pi] = T[T>=2*pi] - 2*pi
  T[T < 0] = T[T<0] + 2*pi
  
  wp = get_wind_profiles(obs)$summary
  
  F = ifelse(obs$cov=='W', 1, 0.8)
  
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
  I = ifelse(eye$cov=='W', 20, 40)*pi/180
  Az = A - pi
  Az[Az < 0] = Az[Az< 0] + 2*pi
  D = Az - pi - I
  D[D<0] = D[D<0]+2*pi
  
  return(list(Vs=Vs,D=D))
}


cumulative_extent = function(track, res_m, max_radius_km){
  track = rgeos::gBuffer(track, width=max_radius_km*1000)
  trk_ext = raster::extent(track)
  trk_ext[c(1,3)] = floor(trk_ext[c(1,3)]/res_m)*res_m
  trk_ext[c(2,4)] = ceiling(trk_ext[c(2,4)]/res_m)*res_m
  ex = raster::raster( resolution=res_m, ext=trk_ext)
  sp::proj4string(ex) = suppressWarnings(sp::proj4string(track))
  suppressWarnings(ex[] <- NA)
  return(ex)
}

size_pred = function(dat, mods) {
  press_known = ifelse(as.numeric(dat$min_press) > 0,'p','np')
  
  df = dat[, c('max_speed', 'min_press', 'lat')]
  df[, 1:3] = as.numeric(df[,1:3])
  #predict inner-most windspeed
  r34 = ifelse(df$max_speed>34, predict(mods[[1]][[press_known]], df),0)
  r50 = ifelse(df$max_speed>50, predict(mods[[2]][[press_known]], df),0)
  r64 = ifelse(df$max_speed>64, predict(mods[[3]][[press_known]], df),0)
  
  # predict outwind speeds from inner if availabile
  if(r64 > 0) r50 = 1.521*r64+15.834
  if(r50 > 0) r34 = 1.707*r50+31.473
  size = list(r34=r34, r50=r50, r64=r64)
  
  #apply assymtry factors
  assym = c(1.28, 1.11, 0.69, 0.91)
  out = c()
  colnams = c()
  for(i in 1:3){
    x = size[[i]] * assym
    out = c(out, x)
    colnams = c(colnams, paste0(substr(names(size)[i],2,3), 'kt_', c('ne', 'se','sw', 'nw')))
  }
  
  #format output
  out = as.data.frame(t(out))
  colnames(out) = colnams
  return(out)
}

# Function to map relative wind profile based on 2 or 3 radius with known speed
wind_profile_fxn = function(R, Rm, B) sqrt((Rm/R)^B *exp(1- (Rm/R)^B)) # From Boose et al. 2004

