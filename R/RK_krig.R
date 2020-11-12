RK_krig <- function(data,grid,data.val,rasters) {

  control <- trainControl(method = "repeatedcv",number = 10,repeats = 5)
  train.grid <- expand.grid(committees = c(1,5,9),neighbors = c(1,5,9))

  ## Prep the data
  error.temp <- data.frame(rep = seq(1,1,1),rmse = rep(NA,1),
                           mse = rep(NA,1),mae = rep(NA,1),
                           rsq = rep(NA,1),rsqmod = rep(NA,1))

  # extract overlapping points from the raster stack, and the subsample points
  DSM_data <- raster::extract(rasters,data,sp = 1,method = "simple")
  mDat <- as.data.frame(DSM_data)

  DSM_VAL <- raster::extract(rasters,data.val,sp = 1,method = "simple")
  vDat <- as.data.frame(DSM_VAL)

  #crs(DSM_data)
  #crs(rasters)

  ## Apply cubist model
  #ptm <- proc.time()
  cub.mod <- train(x = mDat[,names(rasters)],y = mDat[,l],
                   method = "cubist", trControl = control,
                   tuneGrid = train.grid) #metric = "Rsquared" to use instead of RMSE for optimal model slection
  #proc.time() - ptm
  #derive model residuals - model predictions subtracted from the residual
  mDat$residual <- mDat[,l] - predict(cub.mod,newdata=mDat,
                                      neighbors = cub.mod[["bestTune"]][["neighbors"]])
  #mean(mDat$residual)

  coordinates(mDat) <- ~Easting+Northing
  crs(mDat) <- proj4Str

  ## Krige the residuals
  # variogram
  vgm <- autofitVariogram(mDat$residual~1,mDat,model = c("Sph","Exp","Gau","Ste"))

  #cubist predictions
  Cubist.pred.V <- predict(cub.mod,newdata = vDat,neighbors = cub.mod[["bestTune"]][["neighbors"]])

  #make residuals predictions
  RK.preds.V <- as.data.frame(krige(residual ~ 1,mDat,model=vgm$var_model,newdata=DSM_VAL))

  # sum of the two components together
  RK.preds.fin <- Cubist.pred.V + RK.preds.V[,3]

  #vtakeaway rk.preds.fin from measured for rmse
  error.temp$rmse <- rmse(data.val@data[[l]],RK.preds.fin)
  error.temp$mse <- mse(data.val@data[[l]],RK.preds.fin)
  error.temp$mae <- mae(data.val@data[[l]],RK.preds.fin)
  error.temp$rsq <- rsq(data.val@data[[l]],RK.preds.fin)
  error.temp$rsqmod <- max(cub.mod[["results"]][["Rsquared"]])

  return(error.temp)

}
