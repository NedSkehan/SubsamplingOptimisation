mse <- function(actual,predicted) {
  sum((actual-predicted)^2)/length(actual)
}
