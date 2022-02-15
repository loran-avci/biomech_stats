## plot_gene_models.R
## Author: Loran Avci, avci@zhaw.ch
## Version: 30.08.2021

plot_gene_models <- function(dat, bestmods, dat_path = paste0(getwd(),'/data/') , save_path =  paste0(getwd(), '/reports/')){
  "
  description : { Plots performance metrics from lm and best models by AIC for each response. }
  arguments   : { dat       : dataframe containing data for responses and independent variables
                  bestmods  : 
                  dat_path  :
                  save_path : directory to save dataframe with best models
                }
  returns     : {RStudioGD}
  "
  tic()
  # interpolate missing data
  if(paste0(str_split(bestmods,"[_ .]")[[1]][3], collapse = "_") == 'interpol'){
    dat <- interpol_dat(dat, indep, method = 'median')
  }
  # setup ----
  reg_mode <- paste0(str_split(bestmods,"[_ .]")[[1]][3:4], collapse = "_")
  file_pdf <- gsub(":","", paste0("TensMechProp_",reg_mode, "_", gsub((" "), "_", Sys.time(), fixed = TRUE), '.pdf'))
  file_pdf <- paste0(save_path, file_pdf)
  
  load(paste0(dat_path,bestmods))
  
  message("Save PDF in : ", file_pdf)
  pdf(file = file_pdf, width = 28.3 , height =  11.7, paper = 'a4r', title = 'mutation_regression')
  response <- bestmods$target
  
  # init excel outputs 
  df_all <- NULL
  fitlist <- list()
  fitlist_2 <- list()
  tar_list <- NULL
  
  # plot summary df
  tt <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.3)),
    colhead = list(fg_params=list(cex = 0.3)),
    rowhead = list(fg_params=list(cex = 0.3)))
  grid.arrange(tableGrob(bestmods, theme = tt),top = textGrob('Summary',vjust = -0.1, gp = gpar(fontsize = 10)))
  # iterate through response Variables
  for (i in seq_along(response)){
    
    message("\n",i ,". Target: ", response[i])
    # extract predictor variables
    st <- str_split(bestmods$mod[i], "[ :]")[[1]]
    st <- unique(st)
    indep <- st[st != "" & st!= "+"]
    dat0 <- completeFun(dat, c(response[i], indep))
    
    if (dim(dat0)[1] > 20){
      # ols linear regression----
      form <- paste0(bestmods$target," ~ ", bestmods$mod)[i]
      fit1 <- lm(as.formula(form), dat0)
      
      # excel output
      fitlist[i] <-  list(round(fit1$coefficients,4))
      fitlist_2[i] <-  list(round(summary(fit1)$coefficients[,c(4)],4))
      tar_list[i] <- response[i]
      
      # coefficients
      df <- as.data.frame(round(summary(fit1)$coefficients[,c(1,2,4)],4))
      df$Predictor <- rownames(df)
      df <- df[,c(4,1,2,3)]
      df$Sig <- ifelse(df[,4] < 0.05, "x", "")
      rownames(df) <- NULL
      colnames(df)[4] <- "P-value"
      # model performance df
      fitted_value <- NULL
      fitted_value <- fit1$fitted.values
      perf<- data.frame(
        R2 =  round(summary(fit1)$r.squared,4),
        R2adj = round(summary(fit1)$adj.r.squared,4),
        RMSE = RMSE(fitted_value, dat0[,response[i]],na.rm = TRUE),
        MAE = MAE(fitted_value,dat0[,response[i]] , na.rm = TRUE),
        MAPE = mean(abs((dat0[,response[i]]-fitted_value)/dat0[,response[i]]),na.rm = TRUE),
        N = nobs(fit1))
      perf <- round(perf,4)
      rownames(perf) <- "Inference Model"
      # fitted vs obs. plot
      g_ols <- qplot(x = dat0[,response[i]], y= fitted_value)+geom_point()+geom_abline(intercept = 0, slope = 1, col="red") +  
        xlim(min(fitted_value, dat0[,response[i]], na.rm = TRUE), max(fitted_value, dat0[,response[i]], na.rm = TRUE))+
        ylim(min(fitted_value, dat0[,response[i]], na.rm =TRUE),max(fitted_value, dat0[,response[i]], na.rm = TRUE))+
        xlab("Observed") +ylab("Fitted") + ggtitle("OLS Inference Model")
      
      # Leave one Out CV ----
      fitted_value_loo <- NULL
      for(j in 1:nrow(dat0)){
        if (length(indep<2)){ 
          validation <-  as.data.frame(dat0[j, indep])
          colnames(validation) <- indep
        }
        else{
          validation <-  as.data.frame(dat0[j, indep])
        }
        
        training <- as.data.frame(dat0[-j, c(response[i],indep)])
        fit0 <- lm(formula = form , data = training )
        fitted_value_loo[j] <- predict(fit0, newdata = validation)[1]
      }
      
      dfloo<- data.frame(
        K = bestmods[i, 1],
        R2 = R2(fitted_value_loo, dat0[,response[i]], na.rm = TRUE),
        R2adj = R2(fitted_value_loo,dat0[,response[i]],formula = "traditional", na.rm = TRUE),
        RMSE = RMSE(fitted_value_loo, dat0[,response[i]],na.rm = TRUE),
        MAE = MAE(fitted_value_loo,dat0[,response[i]] , na.rm = TRUE),
        MAPE = mean(abs((dat0[,response[i]]-fitted_value_loo)/dat0[,response[i]]),na.rm = TRUE),
        N = nrow(dat0)
      )
      dfloo <- round(dfloo,4)
      rownames(dfloo) <- "Leave-One-Out CV Model"
      
      g_loo <- qplot(x = dat0[,response[i]], y= fitted_value_loo)+geom_point()+geom_abline(intercept = 0, slope = 1, col="red") +  
        xlim(min(fitted_value_loo, dat0[,response[i]], na.rm = TRUE), max(fitted_value_loo, dat0[,response[i]], na.rm = TRUE))+
        ylim(min(fitted_value_loo, dat0[,response[i]], na.rm =TRUE),max(fitted_value_loo, dat0[,response[i]], na.rm = TRUE))+
        xlab("Observed") +ylab("Fitted")+ ggtitle("Leave-One-Out Cross Validation")
      
      rownames(dfloo) <- "Leave-One-Out CV Model"
      
      # plot dfs and plots in pdf
      tt <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = 0.5)),
        colhead = list(fg_params=list(cex = 0.5)),
        rowhead = list(fg_params=list(cex = 0.5)))
      
      grid.arrange(tableGrob(df, theme = tt), tableGrob(perf, theme = tt),g_ols, 
                   #ggplot()+ theme_void(),
                   tableGrob(dfloo, theme = tt),g_loo,
                   ncol=3,layout_matrix =  rbind(c(1,2,3),
                                                 c(1,4,5)),top = textGrob(response[i], gp = gpar(fontsize = 20)))
      
      par(mfrow=c(2,2))
      plot.lmSim(obj = fit1, SEED = 2021, Nsim = 100)
      plot(fit1, which = 5)
    }else{
      print(paste0("Not enough Obs for: ", response[i]))
    }
  }
  
  #merge excel output
  df_coef <- plyr::ldply(fitlist, rbind)
  df_pval <- plyr::ldply(fitlist_2, rbind)
  df_coef$Target <- tar_list
  df_pval$Target <- tar_list
  n <- dim(df_coef)[2]
  df_coef <- df_coef[,c(n, 1:n-1)]
  df_pval <- df_pval[,c(n, 1:n-1)]
  # write excel
  write_xlsx(df_coef,"coef_file.xlsx")
  write_xlsx(df_pval,"pval_file.xlsx")
  
  toc()
  dev.off()
}

