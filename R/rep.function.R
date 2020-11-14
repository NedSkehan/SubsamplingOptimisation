#' A function for repetition of Oridinary Kriging
#'
#' This function will iterate through a vector of subsample sizes
#'
#' @param css vector of subsample lengths
#' @param depth_spdf SpatialPointsDataFrame of data
#' @param grid grid to krige across
#' @param replications numeric of required replications
#'
#' @return list of error indicies for each subsample size
#' @export
#' @examples
#' rep.function(vec,spdf,grid)


rep.function <- function(css,depth_spdf,grid,replications) {

  error.vec <- data.frame(css = css,
                          rmse = rep(NA, length(css)),
                          mse = rep(NA, length(css)),
                          mae = rep(NA, length(css)),
                          rsq = rep(NA, length(css)))

  error.dist <- list()

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



    for (j in 1:length(subsample.rep)) {

      # partition subsample into calibration and validation
      #cal.ind <- subsample.rep[[j]]
      subsample.cal <- subsample.rep[[j]][16:nrow(subsample.rep[[j]]),]
      subsample.val <- subsample.rep[[j]][1:15,]

      # this produces a data frame with errors of the krige estimate, for the single layer,
      # at the single rep
      depth <- OK_krig(subsample.cal,grid,subsample.val)

      error.reps$rmse[j] <- depth$rmse
      error.reps$mse[j] <- depth$mse
      error.reps$mae[j] <- depth$mae
      error.reps$rsq[j] <- depth$rsq
    }

    ct <- ct + 1
    error.vec$rmse[ct] <- mean(error.reps$rmse)
    error.vec$mse[ct] <- mean(error.reps$mse)
    error.vec$mae[ct] <- mean(error.reps$mae)
    error.vec$rsq[ct] <- mean(error.reps$rsq)

    error.dist[[ct]] <- error.reps

  }
  #parallel::stopCluster(cl)

  names(error.dist) <- css
  output <- list(error.vec,error.dist)

  return(output)
}

#' An Ordinary Kriging function
#'
#' This function will perform orindary kriging
#'
#' @param data SpatialPointsDataFrame of calibration data
#' @param grid grid for kriging
#' @param data.val SpatialPointsDataFrame of validation data
#' @param proj4Str chr string of the projection of data
#' @param l numeric column index of variable for kriging
#'
#' @return list of error indicies for each subsample size
#' @export
#' @examples
#' rep.function.rk(vec,spdf,grid)


OK_krig <- function(data,grid,data.val,proj4Str,l) {

  data.df <- as.data.frame(data)
  data.df <- stats::na.omit(data.df)

  # initialise data frame for prediction results
  #preds <- as.data.frame(matrix(c(rep(NA,1)),nrow = nrow(data.val),ncol = 1))
  #names(preds) <- "layer"

  error.temp <- data.frame(rep = seq(1,1,1),rmse = rep(NA,1),
                           mse = rep(NA,1),mae = rep(NA,1),
                           rsq = rep(NA,1))


  # variogram
  vgm <- automap::autofitVariogram(data.df[,l]~1,data,model = c("Sph","Exp","Gau","Ste"))
  #plot(vgm)

  # ordinary kriging
  dat.krig <- gstat::krige(data.df[,l]~1,data,grid,model = vgm$var_model)
  #plot(dat.krig$var1.pred)

  # convert to SPDF to extract points on the grid
  map.df <- as.data.frame(dat.krig)
  #map.spdf <- SpatialPointsDataFrame(coords=map.df[,1:2], data = map.df, proj4string = crs(proj4Str))

  # convert to raster to plot
  kriged.map <- raster::rasterFromXYZ(map.df,crs = sp::CRS(proj4Str))
  raster::plot(kriged.map$var1.pred)
  raster::plot(data,add = T,pch = 3,col="red")
  raster::plot(data.val,add = T,pch = 1)
  graphics::legend("topright",inset=0.01,legend = c("OK Cal","Val"),col = c("red","black"),cex = 0.8,pch=c(3,1))

  output.pred <- kriged.map$var1.pred
  output.var <- kriged.map$var1.var

  #Create raster stack of kriged layers
  #if(map.idx==1){
  #output <- output.pred
  #}else{
  #  output <- stack(output,output.pred)
  #}

  #cat("Layer",layers[id],"is complete","\n")

  ##### calculate predicted values
  prediction <- as.data.frame(gstat::krige(data.df[,l]~1,data,model=vgm$var_model,newdata = data.val))
  preds <- prediction[,3]

  # temp df for storing each error
  error.temp$rmse <- rmse(data.val@data[[l]],preds)
  error.temp$mse <- mse(data.val@data[[l]],preds)
  error.temp$mae <- mae(data.val@data[[l]],preds)
  error.temp$rsq <- rsq(data.val@data[[l]],preds)

  # compile the output
  #names(output) <- layers
  #output.list <- list(output.pred,preds)

  return(error.temp)
}
