renice_session <- function() system(paste0("renice -n 15 -p ", Sys.getpid()))

get.XAX_ann_zonal <- function(X, M, start.year = 1850, cmip) {
  

  ncol <- sqrt(length(X[1][[1]][,1]) /2)
  nrow <- 2*ncol
  
  # OVERALL DIMENSION
  s = sum(sapply(X = X, FUN=function(x) dim(x)[2]))
  
  MAX <- data.frame(matrix(nrow=s,ncol=9))
  names(MAX) = c("vari", "res", "file.name", "cmip", "mod", "modcl", "scen", "ens.mem", "year")
  YAX <- matrix(ncol= ncol, nrow=s)
  names(YAX) <- "ZM"
  
  XAX <- matrix(nrow=s,ncol= nrow * ncol)
  
  cc <- 0
  for (k in 1:length(X)){
    

    cat("\r ",k)
    Xsc <- t(X[[k]])
    # Ysc <- Y[[k]][seq(mon, dim(X[[k]])[3], 12),]
    
    for (pc in 1:(dim(X[[k]])[2])){
      cc <- cc+1
      MAX[cc,] <- as.character(c(M$vari[k], M$res[k], M$file.name[k], cmip, M$mod[k], M$modcl[k], M$scen[k], M$ens.mem[k], start.year-1+as.numeric(pc)))
      #YAX[cc,] <- Xsc[pc,] %*% areaw
      
      YAX[cc,] <- colMeans(matrix(Xsc[pc,], nrow, ncol))
      
      XAX[cc,]  <- Xsc[pc,]
    }
    
    print(MAX[cc,])
  }
  Y <- list()
  Y$ZM <- YAX
  return(list(X = XAX, Y=Y, M=MAX))
}  

#========================================================

extract.forced.response.zonal <- function(XAX_scen, nmem=1) {
  
  nbands <- length(XAX_scen$Y$ZM[1,])
  scen.un = unique(XAX_scen$M$scen)
  l = dim(XAX_scen$M)[1]
  
  XAX_scen$Y$fraw = matrix(nrow = l, ncol = nbands)
  XAX_scen$Y$fr_rm   = matrix(nrow = l, ncol = nbands)

  
  for (s in 1:length(scen.un)) {
    print(s)
    cur.scen = scen.un[s]
    
    # Determine number of members per model:
    no.mem = sapply(X = unique(XAX_scen$M$mod), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod & XAX_scen$M$scen == cur.scen)])))
    cur.mod.un = names(which(no.mem >= nmem))
    fraw <- list()
    for (mod in 1:length(cur.mod.un)) {
      for (lat in 1:nbands) {
        
      ix = which(XAX_scen$M$mod == cur.mod.un[mod] & XAX_scen$M$scen == cur.scen)
      

    
      ens.mat = sapply(X = unique(XAX_scen$M$ens.mem[ix]), FUN=function(cur.ens) XAX_scen$Y$ZM[which(XAX_scen$M$mod == cur.mod.un[mod] & XAX_scen$M$scen == cur.scen & XAX_scen$M$ens.mem == cur.ens),lat])
      fraw = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat))
      # fraw_fIV = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat[,-k]))
      
      # fr.loess_0.75 = loess(fr.raw ~ c(1:231), span = 0.75, degree = 2)$fittedn=dim(ens.mat)[1]
      # 231 * 0.25 / n
      n = dim(ens.mat)[1]
      fr_rm <- rollmean(fraw[,1], k = 30, fill = NA) %>% rep(no.mem[mod])
      #fl = apply(X = fraw, MARGIN = 2, FUN=function(x) loess(x ~ c(1:n), span = 231 * 0.25 / n, degree = 2)$fitted)
      XAX_scen$Y$fraw[ix,lat] = c(fraw)
      XAX_scen$Y$fr_rm[ix,lat] = fr_rm
      }
      }
    
   }
  return(XAX_scen) 
} 


#=============================================== 
LENS_pr_adjust_units <- function(LENS_XAX) {
  LENS_XAX$Y$fraw <- LENS_XAX$Y$fraw*24*3600*30
  LENS_XAX$Y$AGMT <- LENS_XAX$Y$AGMT*24*3600*30
  LENS_XAX$Y$fr_rm   <- LENS_XAX$Y$fr_rm*24*3600*30
  LENS_XAX$X      <- LENS_XAX$X *24*3600*30
  LENS_XAX$Y$ZM <-   LENS_XAX$Y$ZM *24*3600*30
  LENS_XAX$Y$CZM <-   LENS_XAX$Y$CZM *24*3600*30
  return(LENS_XAX)
}

#==============================================
LENS_psl_adjust_units <- function(LENS_XAX){
  
  GFDL.ix = which(LENS_XAX$M$mod %in% c("GFDL-CM3","GFDL-ESM2M"))
  LENS_XAX$Y$AGMT[GFDL.ix] <- LENS_XAX$Y$AGMT[GFDL.ix] * 100
  LENS_XAX$Y$fraw[GFDL.ix] <- LENS_XAX$Y$fraw[GFDL.ix] * 100
  LENS_XAX$Y$fl[GFDL.ix] <- LENS_XAX$Y$fl[GFDL.ix] * 100
  LENS_XAX$X[GFDL.ix,] <-  LENS_XAX$X[GFDL.ix,] *100
  return(LENS_XAX)
}


#=============================================
fit_ridge_global <- function(pr_na, LENS_pr_XAX, LENS_X_sec, train.ix, multi) {
  ## (1.1) Train CESM:
  # create crossclass vector
  crossclass <- tibble(mod = as.factor(LENS_pr_XAX$M$mod[train.ix])) %>% 
    mutate(crossclass = as.numeric(mod)) %>% .$crossclass
  
  crossclass <- rep(1:10, length(train.ix)/10) %>% sort()
  
  #create weights according to ensemble size  
  models <- LENS_pr_XAX$M$mod[train.ix] %>% unique()
  weights_list <- list()
  
  
  for( y in 1:length(models)){
    n_ens <- LENS_pr_XAX$M$ens.mem[train.ix[which(LENS_pr_XAX$M$mod[train.ix]== models[y])]] %>% unique() %>% length()
    len_mod <- LENS_pr_XAX$M$ens.mem[which(LENS_pr_XAX$M$mod[train.ix] == models[y])] %>% length()
    weights_list[[y]] <- rep(1/n_ens, len_mod)
  }
  weights <- unlist(weights_list)
  
  if(multi == T){
    X_masked <-  cbind(LENS_pr_XAX$X[train.ix, -pr_na], LENS_X_sec)
  }else{         X_masked <-  LENS_pr_XAX$X[train.ix, -pr_na]}
  
  #parallel computation of the five models for each extrema
  
  
  mod <-   cv.glmnet(x = X_masked, y = LENS_pr_XAX$Y$GM[train.ix], 
              alpha = 0, foldid = crossclass,  parallel = T, weights = weights)
  

  
  return(mod) 
}



#===========================================================
fit_ridge_weighted_zonal <- function(pr_na, LENS_pr_XAX, LENS_X_sec, train.ix, multi, target, d0){
    ### TRAIN ON GLMNET
    registerDoParallel(cores=24)
    glmnet.control(fdev=0, devmax=1, mnlam = 100)
    
    ## (1.1) Train CESM:
    # create crossclass vector
    crossclass <- tibble(mod = as.factor(LENS_pr_XAX$M$mod[train.ix])) %>% 
      mutate(crossclass = as.numeric(mod)) %>% .$crossclass

    #create weights according to ensemble size  
    models <- LENS_pr_XAX$M$mod[train.ix] %>% unique()
    weights_list <- list()
  
    
   for( y in 1:length(models)){
    n_ens <- LENS_pr_XAX$M$ens.mem[train.ix[which(LENS_pr_XAX$M$mod[train.ix]== models[y])]] %>% unique() %>% length()
    len_mod <- LENS_pr_XAX$M$ens.mem[which(LENS_pr_XAX$M$mod[train.ix] == models[y])] %>% length()
    weights_list[[y]] <- rep(1/n_ens, len_mod)
    }
    weights <- unlist(weights_list)

    if(multi == T){
      X_masked <-  cbind(LENS_pr_XAX$X[train.ix, -pr_na], LENS_X_sec)
    }else{         X_masked <-  LENS_pr_XAX$X[train.ix, -pr_na]}
    
    #parallel computation of the models at each latitude
    lats <- which(!is.na(target[16,]))

    if(missing(d0)){ #no penalty
      mod <- foreach(lat = lats) %dopar% {
        cv.glmnet(x = X_masked, y = target[train.ix, lat],
                  alpha = 0, foldid = crossclass,  parallel = T, weights = weights)}
      
#      }else if(length(d0)==  1){ #exponentially increasing distancy penalty
#    mod <- foreach(lat = lats) %dopar% {
#      penalty <- get_penalty(d0, lat)
#      cv.glmnet(x = X_masked, y = target[train.ix, lat],
#                alpha = 0, foldid = crossclass,  parallel = T, weights = weights, penalty.factor = penalty[-pr_na])}
    
      }else if(length(d0 > 1)){ #when custom penalty vector is supplied
        mod <- foreach(lat = lats) %dopar% {
          cv.glmnet(x = X_masked, y = target[train.ix, lat],
                    alpha = 0, foldid = crossclass,  parallel = T, weights = weights, penalty.factor = d0[-pr_na])}
      }
  return(mod) 
}



#=========================================================
map_betas_norm <- function(mask, mod, sds, lambda = "min", break_sec){
  
  if (lambda == "1se"){
    ix = which(mod$lambda == mod$lambda.1se)
  }else{
    ix = which(mod$lambda == mod$lambda.min)}
  
  if (missing(break_sec)){
  betas <- mod$glmnet.fit$beta[,ix] * sds
  }else{
    betas <- mod$glmnet.fit$beta[break_sec,ix] * sds
  }
  
  beta.shape = raster(nrows = 36, ncols = 72, xmn = 0, xmx = 360)
  
  
  dummy <- rep(T, 2592)
  dummy[mask] <- F
  dummy_martix = matrix(dummy, 72, 36)
  dummy_martix[which(dummy_martix)] = betas
  values(beta.shape) = c(matrix(dummy_martix, 72, 36)[,36:1])
  
  spdf <- as(beta.shape, "SpatialPixelsDataFrame")
  spdf <- as.data.frame(spdf)  %>% 
    mutate(x = case_when(x > 180 ~ x-360, x < 180 ~ x),
           layer = case_when(layer != 0 ~ layer, layer == 0 ~ NA_real_)) %>% 
    rename("Betas" = layer) %>% 
    mutate(NormBetas = Betas/ max(abs(Betas), na.rm = T))

  
  return(spdf)
}



#===================================================


map_betas_zonalband <- function(lambda, mask, beta, ix.na, lat.ix){
  beta.ix <- seq(1:(2592 - length(ix.na)))
  
  betas=cbind(beta[beta.ix,lambda])
  beta.shape = raster(nrows = 36, ncols = 72, xmn = 0, xmx = 360)
  
  vec1 = matrix(mask, 72, 36)
  vec1[which(vec1)] = betas
  values(beta.shape) = c(matrix(vec1, 72, 36)[,36:1])
  
  zone <- matrix(rep(NA, 2592), 72, 36)
  zone[,lat.ix] <- 1
  zone <- c(matrix(zone, 72, 36)[,36:1])
  
  spdf <- as(beta.shape, "SpatialPixelsDataFrame")
  spdf <- as.data.frame(spdf)  %>% mutate(x = case_when(x > 180 ~ x-360, x < 180 ~ x),
                                          layer = case_when(layer != 0 ~ layer, layer == 0 ~ NA_real_),
                                          zone = zone) %>% 
    rename("Betas" = layer)
  
  return(spdf)
}



#===========================================================


calc_naive_ADZMP <- function(LENS_pr_XAX, LENS_pr_norm_XAX, ix.na){
require(zoo)
tot_len <- length(LENS_pr_XAX$M$year)
LENS_models <- unique(LENS_pr_XAX$M$mod)

#initialize various storage containters
AZMP_naive <-         matrix(NA, tot_len, 36)
AZMP_smooth_naive <-  matrix(NA, tot_len, 36)
AZMP_naive_norm <-    matrix(NA, tot_len, 36)
fraw_naive <-         matrix(NA, tot_len, 5)
ADZMP_naive <-        matrix(NA, tot_len, 5)


#loop selecting individual models
for (mod in 1:length(LENS_models)){

mod.ix <- which(LENS_pr_XAX$M$mod == LENS_models[i])
n_ens <- LENS_pr_XAX$M$ens.mem[mod.ix] %>% unique() %>% length()

#calculate zonal means only using grid points where observations are present
for (i in 1:length(mod.ix)){
  naive_dummy <- LENS_pr_XAX$X[mod.ix[i],]
  naive_dummy_norm <- LENS_pr_norm_XAX$X[mod.ix[i],]
  
  naive_dummy[ix.na] <- NA
  naive_dummy_norm[ix.na] <- NA
    
  AZMP_naive[mod.ix[i], ] <- naive_dummy %>% matrix(72,36) %>% colMeans(na.rm = T)
  AZMP_naive_norm[mod.ix[i], ] <- naive_dummy_norm %>% matrix(72,36) %>% colMeans(na.rm = T)
}

#apply 30 year running mean
k <- 30
n_zone <- 36
timesteps <- LENS_pr_XAX$M$year[mod.ix] %>% unique() %>% length()


for (i in 1:n_zone){
  AZMP_zone <-  AZMP_naive[mod.ix, i] %>% matrix(timesteps)
  AZMP_smooth_naive[mod.ix ,i] <- c(rollapplyr(AZMP_zone, k, mean, partial = TRUE))
}

#find zonal maxima

s_max <- list()
s_min <- list()
t_max <- list()
n_min <- list()
n_max <- list()

for( i in 1:length(mod.ix)){
  
  t_max[[i]] <- which.max(AZMP_smooth_naive[mod.ix[i],10:26])+ 9
  s_min[[i]] <- which.min(AZMP_smooth_naive[mod.ix[i],9:t_max[[i]]])+ 8
  n_min[[i]] <- which.min(AZMP_smooth_naive[mod.ix[i],t_max[[i]]:29])+ t_max[[i]]-1
  
  s_max[[i]] <- which.max(AZMP_smooth_naive[mod.ix[i],9:s_min[[i]]])+ 8
  n_max[[i]] <- which.max(AZMP_smooth_naive[mod.ix[i],n_min[[i]]:32])+ n_min[[i]] -1
}

ix_smax <- unlist(s_max)
ix_smin <- unlist(s_min)
ix_tmax <- unlist(t_max)
ix_nmin <- unlist(n_min)
ix_nmax <- unlist(n_max)

#find the value of the extreme values
AZMP_smax <- diag(AZMP_naive_norm[mod.ix, ix_smax])
AMPZ_smin <- diag(AZMP_naive_norm[mod.ix, ix_smin])
AZMP_tmax <- diag(AZMP_naive_norm[mod.ix, ix_tmax])
AZMP_nmin <- diag(AZMP_naive_norm[mod.ix, ix_nmin])
AZMP_nmax <- diag(AZMP_naive_norm[mod.ix, ix_nmax])

AZMP_temp <- as.matrix(cbind(AZMP_smax, AMPZ_smin, AZMP_tmax, AZMP_nmin, AZMP_nmax))

fraw_temp <- list()
for(i in 1:dim(AZMP_temp)[2]){
  fraw_temp[[i]] <- matrix(AZMP_temp[,i], timesteps) %>% rowMeans() %>% rep(n_ens)
}

#fraw_naive[mod.ix,] <- matrix(unlist(fraw_temp), length(mod.ix), 5)
ADZMP_naive[mod.ix,] <- AZMP_temp

}
return(ADZMP_naive)

}


#=================================================================
rep_row <- function(vec, n){
  mat <- matrix(vec, 1, length(vec))
  mat <- mat[rep(seq_len(nrow(mat)), n), ]
  return(mat)
}
#==========================================================

get_sum_FR <- function(LENS_XAX, span) {
  temp <- list()
  for (lat in 1:length(latitudes)){
    temp[[lat]] <- tibble(
      FR = LENS_XAX$Y$fr_rm[,lat],
      year = LENS_XAX$M$year,
      model = LENS_XAX$M$mod,
      ens =  LENS_XAX$M$ens.mem,
    ) %>% filter(year %in% span) %>% 
      group_by(model, ens) %>% 
      summarise(across(where(is.double), ~lm(. ~ year)$coefficients[2]),
                .groups = "drop") %>% 
      tibble(lat = latitudes[lat])
    
  }
  FR_sum <- bind_rows(temp) %>% group_by(lat) %>% 
    summarise(med_trend = median(FR),
              upp_trend = max(FR),
              low_trend = min(FR))
  
  return(FR_sum)
}

#==================================================================
#produces a df where the predictions for each individual ensemble member or just the the first one.
get_pred <-  function(LENS_X, 
                      LENS_X_second, 
                      mod, 
                      lambda, 
                      multi){

   if (multi == F){X_test <- LENS_X
   }else if (multi == T){
     X_test <- cbind(LENS_X, LENS_X_second)}
    
    #calculate prediction, of the form Y = X*B + A
    if (lambda == "1se"){lambda_ix <- which(mod$lambda == mod$lambda.1se)
    }else if (lambda == "min"){
      lambda_ix <- which(mod$lambda == mod$lambda.min)}
    
    pred = X_test %*% mod$glmnet.fit$beta[,lambda_ix] + rep_row(mod$glmnet.fit$a0[lambda_ix], dim(X_test)[1])

  return(pred)
}


#===========================================================
get_land_na <- function(res = 5){

library(rworldmap)
library(raster)
  raster.template = raster(res = res, xmn = 0, xmx=360, ymn = -90, ymx=90)

  y <- dim(raster.template)[1]
  x <- dim(raster.template)[2]
  
coord = coordinates(raster.template) #save coords,
corr.cor <- which(coord[,1] > 180)
coord[corr.cor,1] = coord[corr.cor,1] - 360 
test.points = SpatialPoints(coords = coord)
crs(test.points) = crs(countriesCoarse)
ctry.vec <- over(x = test.points, y = countriesCoarse)$ADM0_A3 %in% as.character(countriesCoarse$ADM0_A3)
#ctry.vec[2232:2592] <- F #removes antarticta

ctry.vec <- c(matrix(!ctry.vec, x, y)[,y:1])
land_na = which(ctry.vec)

return(land_na)
}

#===================================================

get.XAX_ann_piC <- function(X, ncol = 72*36, M, start.year = 1, center = T, scale = F, cmip) {
  
  # OVERALL DIMENSION
  s = sum(sapply(X = X, FUN=function(x) dim(x)[2]))
  
  MAX <- data.frame(matrix(nrow=s,ncol=9))
  names(MAX) = c("vari", "res", "file.name", "cmip", "mod", "modcl", "scen", "ens.mem", "year")
  YAX <- data.frame(matrix(nrow=s, ncol=5))
  names(YAX) <- c("AGMT", "GMT.10y", "GMT.20y", "GMT.30y", "GMT.50y")
  
  XAX <- matrix(nrow=s,ncol=ncol)
  
  cc <- 0
  for (k in 1:length(X)){
    cat("\r ",k)
    Xsc <- t(X[[k]])
    
    if (center == T & scale == F) {
      Xsc = scale(Xsc, T, F)
    } else if (scale == T & scale == T) {
      Xsc = scale(Xsc, T, T)
    }
    # Ysc <- Y[[k]][seq(mon, dim(X[[k]])[3], 12),]
    
    # do slow averages:
    AGMT = c(Xsc %*% areaw)
    AGMT.10y = rollmean(x = AGMT, k = 10, fill = NA)
    AGMT.20y = rollmean(x = AGMT, k = 20, fill = NA)
    AGMT.30y = rollmean(x = AGMT, k = 30, fill = NA)
    AGMT.50y = rollmean(x = AGMT, k = 50, fill = NA)
    
    for (pc in 1:(dim(X[[k]])[2])){
      cc <- cc+1
      MAX[cc,] <- as.character(c(M$vari[k], M$res[k], M$file.name[k], cmip, M$mod[k], M$modcl[k], M$scen[k], M$ens.mem[k], start.year-1+as.numeric(pc)))
      YAX[cc,1] <- AGMT[pc]
      YAX[cc,2] <- AGMT.10y[pc]
      YAX[cc,3] <- AGMT.20y[pc]
      YAX[cc,4] <- AGMT.30y[pc]
      YAX[cc,5] <- AGMT.50y[pc]
      XAX[cc,]  <- Xsc[pc,] 
    }
    print(MAX[cc,])
  }
  
  # save to .RData file:
  return(list(X = XAX, Y=YAX, M=MAX))
}


flip_burger2 <- function(X){
  X <- c(matrix(X,72, 36)[,36:1])
  return(X)}


get_obs_XAX <- function(obs_brick){
  obs_XAX <- list()
  obs_XAX$X <- values(obs_brick) %>% apply(MARGIN = 2, FUN = flip_burger2) %>% t()
  obs_XAX$Y$ZM <- apply(X = obs_XAX$X, FUN = function(X) {matrix(X,72,36) %>% colMeans(na.rm = T)}, MARGIN = 1)
  obs_XAX$M$year <- names(obs_brick) %>% substr(2,5) %>% as.integer()
  return(obs_XAX)
}

get_coast <- function(){
coast <-  map_data(rworldmap::coastsCoarse) %>% dplyr::tibble() %>% dplyr::rename("lati" = "lat") 
return(coast)}

convert.to.pacific <- function(raster.in){
  rl <- crop(raster.in, extent(c(xmin(raster.in),xmax(raster.in)-180,ymin(raster.in), ymax(raster.in))))
  # set extent from rr to rleft:
  rr <- (crop(raster.in, extent(c(xmax(raster.in)-180,xmax(raster.in),ymin(raster.in), ymax(raster.in)))))
  extent(rl) <- extent(c(xmin(rl)+360,xmax(rl)+360,ymin(raster.in), ymax(raster.in)))
  r.out = merge(rr,rl)
  names(r.out) = names(raster.in)
  return(r.out)
}
#======================
convert.to.pacific <- function(raster.in){
  rl <- crop(raster.in, extent(c(xmin(raster.in),xmax(raster.in)-180,ymin(raster.in), ymax(raster.in))))
  # set extent from rr to rleft:
  rr <- (crop(raster.in, extent(c(xmax(raster.in)-180,xmax(raster.in),ymin(raster.in), ymax(raster.in)))))
  extent(rl) <- extent(c(xmin(rl)+360,xmax(rl)+360,ymin(raster.in), ymax(raster.in)))
  r.out = merge(rr,rl)
  names(r.out) = names(raster.in)
  return(r.out)
}
#====================
get_area_weights <- function(res = 2.5){
  raster.template = raster(res = res, xmn = 0, xmx=360, ymn = -90, ymx=90)
  y <- dim(raster.template)[1]
  x <- dim(raster.template)[2]
  areaw=c(matrix(values(raster::area(raster.template)), x,y)[,y:1]) / sum(c(matrix(values(raster::area(raster.template)), x,y)[,y:1]))
  
  return(areaw)
}
#===================
get_lonlatDT <- function(res = 2.5){
  r_template = raster(res = res, xmn = -180, xmx=180, ymn = -90, ymx=90)
  r_template <- convert.to.pacific(r_template)
  lonlatDT <- data.table(
    lat = c(matrix(coordinates(r_template)[,2],144, 72)[,72:1]),
    lon = coordinates(r_template)[,1])
  lonlatDT[lon > 180, lon := lon -360]
  lonlatDT[, lonlat := paste(lon,lat)]
  return(lonlatDT)
}

#======================
pub_theme <- theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.background = element_blank(),#element_rect(fill = "white", size = 4, colour = "white"),
    legend.title = element_blank()
  ) 


#========================
mse <- function(x,y) mean((x-y)^2)
rmse <- function(x,y)  sqrt(mean((x-y)^2))
 
 #======================
 get_na <- function(Obs_XAX, years){
   na <- Obs_XAX$X[which(Obs_XAX$M$year %in%  years),] %>% 
     colMeans() %>% 
     is.na() %>% 
     which() 
   return(na)
 }
#=====================

FR_comp <- function(seas, obs_XAX, mask, mod_fr, lambda, startyear, multi = F, obs_XAX_sec){

  temp <- list()
  temp_LENS <- list()
  temp_PiC <- list()
  
  time_ix <- which(obs_XAX$M$year %in% startyear:2014)
  if (!missing(obs_XAX_sec)){
    time_ix_sec <- which(obs_XAX_sec$M$year %in% startyear:2014)
  }else{time_ix_sec <- 0}
  trend_length <- length(startyear:2014)
  for (lat in 1:36){
    
    pred <- get_pred(obs_XAX$X[time_ix, -mask],
                     obs_XAX_sec$X[time_ix_sec,], 
                     mod_fr[[lat]], 
                     lambda, 
                     multi = multi) %>% c()
    temp[[lat]] <-  data.table(pred = pred,
                               lat = latitudes[lat],
                               year = as.integer(obs_XAX$M$year[time_ix]))
    
    
    if (seas == "DJF"){
      
      GFDL_ix = which(CMIP5_pr_PIC_DJF_XAX$M$mod %in% c("GFDL-CM3","GFDL-ESM2M"))
      LENS_ix <- which(LENS_pr_DJF_5d00_XAX_norm$M$year %in% startyear:2014)
      LENS_pred <- get_pred(LENS_pr_DJF_5d00_XAX_norm$X[LENS_ix, -mask],
                            LENS_psl_DJF_5d00_XAX_norm$X[LENS_ix, ], 
                            mod_fr[[lat]], 
                            lambda, 
                            multi = multi) %>% c() 
      temp_LENS[[lat]] <- data.table(LENS_pred = LENS_pred,
                                     lat = latitudes[lat],
                                     year = as.integer(LENS_pr_DJF_5d00_XAX_norm$M$year[LENS_ix]),
                                     ens = LENS_pr_DJF_5d00_XAX_norm$M$ens.mem[LENS_ix],
                                     mod = LENS_pr_DJF_5d00_XAX_norm$M$mod[LENS_ix])
      
      PiC_pred <- get_pred(CMIP5_pr_PIC_DJF_XAX$X[-GFDL_ix,-mask],
                           CMIP5_psl_PIC_DJF_XAX$X[-GFDL_ix,], 
                           mod_fr[[lat]],
                           lambda,
                           multi = multi) %>% c()
      temp_PiC[[lat]] <-  data.table(PiC_pred = PiC_pred,
                                     lat = latitudes[lat],
                                     year = as.integer(CMIP5_pr_PIC_DJF_XAX$M$year[-GFDL_ix]),
                                     mod = CMIP5_pr_PIC_DJF_XAX$M$mod[-GFDL_ix])
    }
    if (seas == "JJA"){
      
      GFDL_ix = which(CMIP5_pr_PIC_JJA_XAX$M$mod %in% c("GFDL-CM3","GFDL-ESM2M"))
      LENS_ix <- which(LENS_pr_JJA_5d00_XAX_norm$M$year %in% startyear:2014)
      LENS_pred <- get_pred(LENS_pr_JJA_5d00_XAX_norm$X[LENS_ix, -mask],
                            LENS_psl_JJA_5d00_XAX_norm$X[LENS_ix, ], 
                            mod_fr[[lat]], 
                            lambda, 
                            multi = multi) %>% c() 
      temp_LENS[[lat]] <- data.table(LENS_pred = LENS_pred,
                                     lat = latitudes[lat],
                                     year = as.integer(LENS_pr_JJA_5d00_XAX_norm$M$year[LENS_ix]),
                                     ens = LENS_pr_JJA_5d00_XAX_norm$M$ens.mem[LENS_ix],
                                     mod = LENS_pr_JJA_5d00_XAX_norm$M$mod[LENS_ix])
      
      PiC_pred <- get_pred(CMIP5_pr_PIC_JJA_XAX$X[-GFDL_ix,-mask],
                           CMIP5_psl_PIC_JJA_XAX$X[-GFDL_ix,], 
                           mod_fr[[lat]],
                           lambda,
                           multi = multi) %>% c()
      temp_PiC[[lat]] <-  data.table(PiC_pred = PiC_pred,
                                     lat = latitudes[lat],
                                     year = as.integer(CMIP5_pr_PIC_JJA_XAX$M$year[-GFDL_ix]),
                                     mod = CMIP5_pr_PIC_JJA_XAX$M$mod[-GFDL_ix])
    }
    
  }
  
  #Obs trend calculation
  obs_FR <- rbindlist(temp)[,lm(pred ~year)$coefficients[2]*10, by = lat] %>% setnames( "V1" ,"trend")
  
  #LENS trend calculation
  LENS_FR <- rbindlist(temp_LENS)[,lm(LENS_pred ~year)$coefficients[2]*10, by = .(lat, ens, mod)] %>% setnames( "V1" ,"trend")
  
  pic_length <- length(startyear:2014)
  #PIC trend calculation
  PIC_sum <- rbindlist(temp_PiC)[, group := base::cut(year, breaks =seq(0,1500,pic_length)), by = .(mod,lat)][
    ,.(lm(PiC_pred ~year)$coefficients[2]*10, .N),  by = .(mod,lat,group)][N==pic_length] %>% setnames( "V1" ,"trend")
  
  
  return(list(obs_FR, 
              LENS_FR,
              PIC_sum)
  )
}
