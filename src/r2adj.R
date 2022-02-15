## r2adj.R
## Author: Loran Avci
## Version: 26.08.2021

r2adj <- function(R2, n, p){
  r2a <- 1-(1-R2)*((n-1)/(n-p))
  return(r2a)
}