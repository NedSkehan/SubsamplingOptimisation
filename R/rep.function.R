#' A function for repetition of Oridinary Kriging
#'
#' This function will iterate through a vector of subsample sizes
#'
#' @param css vector of subsample lengths
#' @param depth_spdf SpatialPointsDataFrame of data
#' @param grid grid to krige across
#' @param replications numeric of required replications
#' @param property char of property name in 'depth_spdf' for kriging. Case sensitive
#'
#' @return list of error indicies for each subsample size
#' @export
#' @examples



rep.function <- function(css,depth_spdf,grid,replications,property) {

  error.vec <- data.frame(css = css,
                          rmse = rep(NA, length(css)),
                          mse = rep(NA, length(css)),
                          mae = rep(NA, length(css)),
                          rsq = rep(NA, length(css)),
                          rsq_adj = rep(NA, length(css)),
                          sde = rep(NA, length(css)))

  error.dist <- list()

  ct <- 0


  for (s in css) {
    # create a list of all repetitions
    set.seed(s)
    subsample.rep <- replicate(replications,depth_spdf[sample(1:nrow(depth_spdf),s+15,replace = F),],
                               simplify = F)


    error.reps <- data.frame(rep = seq(1,replications,1),
                             rmse = rep(NA,replications),
                             mse = rep(NA,replications),
                             mae = rep(NA,replications),
                             rsq = rep(NA,replications),
                             rsq_adj = rep(NA,replications),
                             sde = rep(NA,replications))



    for (j in 1:length(subsample.rep)) {

      # partition subsample into calibration and validation
      #cal.ind <- subsample.rep[[j]]
      subsample.cal <- subsample.rep[[j]][16:nrow(subsample.rep[[j]]),]
      subsample.val <- subsample.rep[[j]][1:15,]

      # this produces a data frame with errors of the krige estimate, for the single layer,
      # at the single rep
      depth <- OK_krig(subsample.cal,grid,subsample.val,property)

      error.reps$rmse[j] <- depth$rmse
      error.reps$mse[j] <- depth$mse
      error.reps$mae[j] <- depth$mae
      error.reps$rsq[j] <- depth$rsq
      error.reps$rsq_adj[j] <- depth$rsq_adj
      error.reps$sde <- depth$sde
    }

    ct <- ct + 1
    error.vec$rmse[ct] <- mean(error.reps$rmse)
    error.vec$mse[ct] <- mean(error.reps$mse)
    error.vec$mae[ct] <- mean(error.reps$mae)
    error.vec$rsq[ct] <- mean(error.reps$rsq)
    error.vec$rsq_adj[ct] <- mean(error.reps$rsq_adj)
    error.vec$sde <- mean(error.reps$sde)

    error.dist[[ct]] <- error.reps

  }

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
#' @param property char of property name in 'data' for kriging. Case sensitive
#'
#' @return list of error indicies for each subsample size
#' @export
#' @examples



OK_krig <- function(data,grid,data.val,property) {

    l <- grep("property",colnames(data))

    data.df <- as.data.frame(data)
    data.df <- stats::na.omit(data.df)

    # initialise data frame for prediction results
    error.temp <- data.frame(rep = seq(1,1,1),rmse = rep(NA,1),
                             mse = rep(NA,1),mae = rep(NA,1),
                             rsq = rep(NA,1),rsq_adj = rep(NA,1),
                             sde = rep(NA,1))


    # variogram
    vgm <- automap::autofitVariogram(data.df[,l]~1,data,model = c("Sph","Exp","Gau","Ste"))
    #plot(vgm)

    ######################## Just for a visual ############################
    # ordinary kriging
    #dat.krig <- krige(data.df[,l]~1,data,newdata = grid,model = vgm$var_model)

    # convert to SPDF to extract points on the grid
    #map.df <- as.data.frame(dat.krig)

    # convert to raster to plot
    #kriged.map <- rasterFromXYZ(map.df,crs = CRS(proj4Str))
    #plot(kriged.map$var1.pred)
    #plot(data,add = T,pch = 3,col="red")
    #plot(data.val,add = T,pch = 1)
    #legend("topright",inset=0.01,legend = c("OK Cal","Val"),col = c("red","black"),cex = 0.8,pch=c(3,1))

    #output.pred <- kriged.map$var1.pred
    #output.var <- kriged.map$var1.var

    ######################## END: Just for a visual ############################

    ##### calculate predicted values
    prediction <- as.data.frame(gstat::krige(data.df[,l]~1,data,model=vgm$var_model,newdata = data.val))
    preds <- prediction[,3]

    # temp df for storing each error
    error.temp$rmse <- rmse(data.val@data[[l]],preds)
    error.temp$mse <- mse(data.val@data[[l]],preds)
    error.temp$mae <- mae(data.val@data[[l]],preds)
    error.temp$rsq <- rsq(data.val@data[[l]],preds)
    error.temp$rsq_adj <- rsq_adj(data.val@data[[l]],preds,nrow(data.df))
    error.temp$sde <- sde(data.val@data[[l]],preds)

    return(error.temp)
  }
