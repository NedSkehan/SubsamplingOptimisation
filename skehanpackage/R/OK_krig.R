OK_krig <- function(data,grid,data.val) {

  data.df <- as.data.frame(data)
  data.df <- na.omit(data.df)

  # initialise data frame for prediction results
  #preds <- as.data.frame(matrix(c(rep(NA,1)),nrow = nrow(data.val),ncol = 1))
  #names(preds) <- "layer"

  error.temp <- data.frame(rep = seq(1,1,1),rmse = rep(NA,1),
                           mse = rep(NA,1),mae = rep(NA,1),
                           rsq = rep(NA,1))


  # variogram
  vgm <- autofitVariogram(data.df[,l]~1,data,model = c("Sph","Exp","Gau","Ste"))
  #plot(vgm)

  # ordinary kriging
  dat.krig <- krige(data.df[,l]~1,data,grid,model = vgm$var_model)
  #plot(dat.krig$var1.pred)

  # convert to SPDF to extract points on the grid
  map.df <- as.data.frame(dat.krig)
  #map.spdf <- SpatialPointsDataFrame(coords=map.df[,1:2], data = map.df, proj4string = crs(proj4Str))

  # convert to raster to plot
  kriged.map <- rasterFromXYZ(map.df,crs = CRS(proj4Str))
  plot(kriged.map$var1.pred)
  plot(data,add = T,pch = 3,col="red")
  plot(data.val,add = T,pch = 1)
  legend("topright",inset=0.01,legend = c("OK Cal","Val"),col = c("red","black"),cex = 0.8,pch=c(3,1))

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
  prediction <- as.data.frame(krige(data.df[,l]~1,data,model=vgm$var_model,newdata = data.val))
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
