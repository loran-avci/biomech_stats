# substitute NAs in independent vars with column median
interpol_dat <- function(dat, indep, method = 'median'){
  for(i in indep){
    if (is.numeric(dat[,i])){
      if (method == 'median'){
        dat[is.na(dat[,i]), i] <- median(dat[,i], na.rm = TRUE)
      }
      else if (method == 'mean'){
        dat[is.na(dat[,i]), i] <- mean(dat[,i], na.rm = TRUE)
      }else{message("Error: invalid method")}
    }
  }
  return(dat)
}
