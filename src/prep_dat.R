prep_dat <- function(dat){
  dat$E756del <- as.factor(dat$E756del)
  dat$Sex <- as.factor(dat$Sex)
  # change encoding of Gender Variable
  dat$Sex <- ifelse(dat$Sex == "1","m","f")
  dat$Sex <- as.factor(dat$Sex)
  # change numeric activity  variables to factor 
  #10Y
  quants <- round(as.numeric(quantile(dat$MarxScore10y, probs = c(0.3,0.6, 1))));quants
  dat$MarxScore10y_f <- ifelse(dat$MarxScore10y  < quants[1],"low", ifelse(dat$MarxScore10y > quants[2], "high", "medium"))
  dat$MarxScore10y_f <- factor(dat$MarxScore10y_f, levels = c("low", "medium", "high"))
  table(dat$MarxScore10y_f, dat$MarxScore10y)
  # 1Y
  quants <- round(as.numeric(quantile(dat$MarxScore1y, probs = c(0.3,0.6, 1))));quants
  dat$MarxScore1y_f <- ifelse(dat$MarxScore1y  < quants[1],"low", ifelse(dat$MarxScore1y > quants[2], "high", "medium"))
  dat$MarxScore1y_f <- factor(dat$MarxScore1y_f, levels = c("low", "medium", "high"))
  table(dat$MarxScore1y_f, dat$MarxScore1y)
  #PAS
  dat$PAS_f <- ifelse(dat$PAS > median(dat$PAS), "high", "low")
  dat$PAS_f <- factor(dat$PAS_f, levels = c("low", "high"))
  return(dat)
}
