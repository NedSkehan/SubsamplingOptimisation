#' A function for repetition of Regression Kriging
#'
#' This function will iterate through a vector of subsample sizes
#'
#' @param css vector of subsample lengths
#' @param depth_spdf SpatialPointsDataFrame of data
#' @param grid grid to krige across
#' @param rasters RasterStack of environmental covariate rasters
#' @param replications numeric of required replications
#' @param layers chr string of layer names for naming in output list
#'
#' @return list of error indicies for each subsample size
#' @export
#' @examples
#' rep.function.rk(vec,spdf,grid,rasters)


rep.function.rk <- function(css,depth_spdf,grid,rasters,replications,layers) {

  error.vec <- data.frame(css = css,
                          rmse = rep(NA, length(css)),
                          mse = rep(NA, length(css)),
                          mae = rep(NA, length(css)),
                          rsq = rep(NA, length(css)))

  error.dist <- list()

  mod.error <- list()

  ct <- 0

  #cl <- parallel::makeCluster(numcores)
  #doParallel::registerDoParallel(cl)
  #rep.list <-
  for (s in css) {
    # create a list of all repetitions
    set.seed(s)
    subsample.rep <- replicate(replications,depth_spdf[sample(1:nrow(depth_spdf),s+15,replace = F),],
                               simplify = F)


    error.reps <- data.frame(rep = seq(1,replications,1),
                             rmse = rep(NA,replications),
                             mse = rep(NA,replications),
                             mae = rep(NA,replications),
                             rsq = rep(NA,replications))

    mod.error.reps <- list()

    for (j in 1:length(subsample.rep)) {

      # partition subsample into calibration and validation
      #cal.ind <- subsample.rep[[j]]
      subsample.cal <- subsample.rep[[j]][16:nrow(subsample.rep[[j]]),]
      subsample.val <- subsample.rep[[j]][1:15,]

      # this produces a data frame with errors of the krige estimate, for the single layer,
      # at the single rep
      depth <- RK_krig(subsample.cal,grid,subsample.val,rasters)

      error.reps$rmse[j] <- depth$rmse
      error.reps$mse[j] <- depth$mse
      error.reps$mae[j] <- depth$mae
      error.reps$rsq[j] <- depth$rsq
      mod.error.reps[j] <- depth$rsqmod


    }

    names(mod.error.reps) <- layers

    ct <- ct + 1
    error.vec$rmse[ct] <- mean(error.reps$rmse)
    error.vec$mse[ct] <- mean(error.reps$mse)
    error.vec$mae[ct] <- mean(error.reps$mae)
    error.vec$rsq[ct] <- mean(error.reps$rsq)

    error.dist[[ct]] <- error.reps

    mod.error[[ct]] <- mod.error.reps

  }
  #parallel::stopCluster(cl)

  names(error.dist) <- css
  names(mod.error) <- css
  output <- list(error.vec,error.dist,mod.error)

  return(output)
}

#' A Regression Kriging function
#'
#' This function will perform regression kriging
#'
#' @param data SpatialPointsDataFrame of calibration data
#' @param grid grid for kriging
#' @param data.val SpatialPointsDataFrame of validation data
#' @param rasters RasterStack of environmental covariate rasters
#' @param proj4Str chr string of the projection of data
#' @param l numeric column index of variable for kriging
#'
#' @return list of error indicies for each subsample size
#' @export
#' @examples
#' rep.function.rk(vec,spdf,grid,rasters)

RK_krig <- function(data,grid,data.val,rasters,proj4Str,l) {

  control <- caret::trainControl(method = "repeatedcv",number = 10,repeats = 5)
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
  cub.mod <- caret::train(x = mDat[,names(rasters)],y = mDat[,l],
                   method = "cubist", trControl = control,
                   tuneGrid = train.grid) #metric = "Rsquared" to use instead of RMSE for optimal model slection
  #proc.time() - ptm
  #derive model residuals - model predictions subtracted from the residual
  mDat$residual <- mDat[,l] - stats::predict(cub.mod,newdata=mDat,
                                      neighbors = cub.mod[["bestTune"]][["neighbors"]])
  #mean(mDat$residual)

  sp::coordinates(mDat) <- ~Easting+Northing
  raster::crs(mDat) <- proj4Str

  ## Krige the residuals
  # variogram
  vgm <- automap::autofitVariogram(mDat$residual~1,mDat,model = c("Sph","Exp","Gau","Ste"))

  #cubist predictions
  Cubist.pred.V <- stats::predict(cub.mod,newdata = vDat,neighbors = cub.mod[["bestTune"]][["neighbors"]])

  #make residuals predictions
  RK.preds.V <- as.data.frame(gstat::krige(residual ~ 1,mDat,model=vgm$var_model,newdata=DSM_VAL))

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
