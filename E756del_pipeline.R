## E756del_pipeline.R
## Author: Loran Avci
## Version: 18.08.2021

# setup env. ----
rm(list = ls()) # clear env. 
cat("\014") # clear console
try(dev.off(dev.list()["RStudioGD"]), silent = TRUE) # clear graphics

# adjust these paths if necessary
setwd("C:/Users/avci/Documents/Zivi/")
src_path <- paste0(getwd(),'/src/')
dat_path <- paste0(getwd(),'/data/')
rep_path <- paste0(getwd(), '/reports/')

# load packages ----
require(car)
require(readr)
require(robustbase)
require(caret)
require(data.table)
require(gridExtra)
require(ggfortify)
require(grid)
require(condformat)
require(ggplotify)
require(tictoc)
require(stringr)

# load data from local repo ----
dat <- as.data.frame(read_csv(paste0(dat_path,"Master.csv"), col_types = cols())) 
RespList <- as.vector(read_csv(paste0(dat_path,"RespList.csv"), col_types = cols())); response <- RespList$Resp
# define possible independent variables (confounders)
indep <- c('E756del', 'Age', 'Sex', 'PAS_f', 'MarxScore1y_f', 'Height','MarxScore10y_f', 
           'BMI', 'CSA_Ach','CSA_Pat','zeroLength', 'AchTendonLength','PennAngle_mean')

# load functions ----
source(paste0(src_path,'plot-lmSim.R'))
source(paste0(src_path,'get_gene_models.R'))
source(paste0(src_path,'plot_gene_models.R'))
source(paste0(src_path,'completeFun.R'))
source(paste0(src_path,'r2adj.R'))
source(paste0(src_path,'interpol_dat.R'))
source(paste0(src_path,'prep_dat.R'))

# data preparation ----
dat <- prep_dat(dat)

# model generation ----
get_gene_models(dat=dat, 
                response=response, 
                indep=indep, 
                interact = TRUE, 
                interpol = TRUE, 
                save_path = dat_path, 
                best = "R2adj")

# plot model ----
plot_gene_models(dat = dat,  
                 bestmods = 'model_data_interpol_R2adj.RData', 
                 dat_path = dat_path , 
                 save_path = rep_path)

