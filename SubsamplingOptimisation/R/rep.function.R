rep.function <- function(css,depth_spdf) {

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
