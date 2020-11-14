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
#' rmse(actual,predicted)

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
#' mae(actual,predicted)

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
#' mse(actual,predicted)

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
#' rsq(actual,predicted)

rsq <- function(actual,predicted) {
  stats::cor(actual,predicted)^2
}
