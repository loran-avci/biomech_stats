## get_gene_models.R
## Author: Loran Avci, avci@zhaw.ch
## Version: 30.08.2021

get_gene_models <- function(dat, response, indep, interact = TRUE, interpol = TRUE, save_path = paste0(getwd(),'/data/'), best = "RMSE"){
  "
  description : { Saves a dataframe with performance metrics from lm and best models by AIC for each response. }
  arguments   : { dat       : dataframe containing data for responses and independent variables
                  response  : vector containing all relevant targets
                  indep     : vector containing all relevant independent variables (confounders)
                  interact  : if TRUE (default), interaction termes are possible for E756del
                  interpol  : if TRUE (default), missing data is replaced by median (other methods possible)
                  save_path : directory to save dataframe with best models
                  best      : if RMSE (default) best model is with min RMSE or max R2adj, else max R2
                }
  returns     : {Dataframe}
  "
  
  # interpolate missing data
  if(interpol){
    dat <- interpol_dat(dat, indep, method = 'median')
    file_p <- paste0(save_path, "model_data_interpol_",best,".RData")
  }else{file_p <- paste0(save_path, "model_data_no-interpol_",best,".RData")}
  
  # init df
  bestmods <- NULL
  # iterate through targets 
  for (i in seq_along(response)){
    message("\n",i ,". Target: ", response[i])
    tic() # processing time
    dat0 <- completeFun(dat, c(response[i], indep))
    if (dim(dat0)[1] > 20){
      # interaction with mutation
      if (interact ){
        f1 <- paste0(response[i],' ~ ', paste0( indep, collapse = " + "))
        f2 <- paste0( "E756del*", indep[-1], collapse = " + ")
        form <- paste0(f1, " + ", f2)
        fit0 <- lm(as.formula(form), dat0)
      }
      else {
        form <- paste0(response[i],' ~ ', paste0( indep, collapse = " + "))
        fit0 <- lm(as.formula(form), dat0)
      }
      # init df for perormance
      perf_df <- NULL
      fitted_value <- NULL
      # variable selection criteria, grid search for optimal K in stepwise AIC
      k_vec <- seq(from = 1.5 ,to = log(dim(dat0)[1])+1, by = 0.05)
      for (k in k_vec) {
        fit1 <- step(fit0, direction = "both", k = k , trace=0, scope = list(lower = ~ E756del))
        mod <- as.formula(fit1$call)
        # Leave one out Cross validation
        for(j in 1:nrow(dat0)){
          validation <- dat0[j, indep]
          training <- dat0[-j, c(response[i],indep)]
          fit1 <- lm(formula = mod , data = training )
          fitted_value[j] <- predict(fit1, newdata = validation)[1]
        }
        #performance on Leave one Out
        perf<- data.frame(
          K = k,
          R2 = R2(fitted_value,dat0[,response[i]], na.rm = TRUE),
          #R2adj = R2(fitted_value,dat0[,response[i]],formula = "traditional", na.rm = TRUE),
          R2adj = r2adj(R2 = R2(fitted_value,dat0[,response[i]], na.rm = TRUE),n = length(fitted_value), p = length(fit1$coefficients)),
          RMSE = RMSE(fitted_value, dat0[,response[i]],na.rm = TRUE),
          MAE = MAE(fitted_value,dat0[,response[i]] , na.rm = TRUE),
          MAPE = mean(abs((dat0[,response[i]]-fitted_value)/dat0[,response[i]]),na.rm = TRUE),
          N = length(fitted_value)
        )
        # rounding Performance, adding predictors and targets 
        perf <- round(perf,4)
        perf$mod <- as.character(mod)[3]
        perf$target <- as.character(mod)[2]
        perf_df <- rbind(perf_df, perf)
      }
      # best model selection by RMSE or R2
      if (best == "RMSE"){
        bestmods <- rbind(bestmods,perf_df[which.min(perf_df$RMSE),])
      }
      else if (best == "R2adj"){
        bestmods <- rbind(bestmods,perf_df[which.max(perf_df$R2adj),])
      }
      else{
        bestmods <- rbind(bestmods,perf_df[which.max(perf_df$R2),])
      }
    }
    else{print(paste0("Not enough Obs for: ", response[i]))}
    toc()
  }
  # save data
  rownames(bestmods) <- NULL
  save(bestmods, file = file_p)
  return(bestmods)
}