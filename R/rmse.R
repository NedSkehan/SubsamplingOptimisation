#' A function for calculating Root Mean Square Error
#'
#' This function will calculate RMSE
#'
#' @param actual vector of validation data
#' @param predicted vector of predicted data
#'
#' @return numeric of RMSE
#' @export
#' @examples


rmse <- function(actual,predicted) {
  sqrt(mean((actual - predicted)^2))
}

#' A function for calculating Mean Average Error
#'
#' This function will calculate MAE
#'
#' @param actual vector of validation data
#' @param predicted vector of predicted data
#'
#' @return numeric of MAE
#' @export
#' @examples


mae <- function(actual,predicted) {
  sum(abs(actual-predicted))/length(actual)
}

#' A function for calculating Mean Square Error
#'
#' This function will calculate MSE
#'
#' @param actual vector of validation data
#' @param predicted vector of predicted data
#'
#' @return numeric of MSE
#' @export
#' @examples


mse <- function(actual,predicted) {
  sum((actual-predicted)^2)/length(actual)
}

#' A function for calculating Rsquared
#'
#' This function will calculate RSQ
#'
#' @param actual vector of validation data
#' @param predicted vector of predicted data
#'
#' @return numeric of RSQ
#' @export
#' @examples


rsq <- function(actual,predicted) {
  stats::cor(actual,predicted)^2
}

#' A function for calculating Standard Deviation of residuals
#'
#' This function will calculate SDE
#'
#' @param actual vector of validation data
#' @param predicted vector of predicted data
#'
#' @return numeric of SDE
#' @export
#' @examples


sde <- function(actual,predicted) {
  stats::sd(actual-predicted)
}

#' A function for calculating Adjusted Rsquared
#'
#' This function will calculate RSQ_ADJ
#'
#' @param actual vector of validation data
#' @param predicted vector of predicted data
#' @param p numeric of number of predictors used
#'
#' @return numeric of RSQ
#' @export
#' @examples


rsq_adj <- function(actual,predicted,p) {
  1-(1-(stats::cor(actual,predicted)^2))*((p-1)/(p-length(predicted)-1))
}
