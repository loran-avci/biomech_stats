## completeFun.R
## Author: Loran Avci
## Version: 26.08.2021

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
