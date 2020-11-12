mae <- function(actual,predicted) {
  sum(abs(actual-predicted))/length(actual)
}
