#TRAIN RIDGE MODELS TO PREDICT ZONAL MEAN PRECIPITATION FROM INCLOMPLETE COVERAGE
#THERE ARE PRETRAINED MODELS IN THE DIRECTORY

#REPLACE THIS WITH THE PATH FROM WHERE YOU LAUNCH THIS SCRIPT
basedir <- "/net/xenon/climphys_backedup/maegli/Precipitation/Scripts/pub/"
setwd(basedir)

# Read functions: 
source("./functions.R")
source("./00_global_vars.R")


#renice_session()
registerDoParallel(cores=24) #use for parallel computing during training
mod <- list()


#GHCN_na_1900_DJF <- get_na(GHCN_pr_5d00_DJF_anom_XAX, 1900:2014)
#GHCN_na_1900_ <- get_na(GHCN_pr_5d00_JJA_anom_XAX$X, 1900:2014)

# A model is trained for every latitude
# Training these models can require considerable computational time, as well as memory but there are predtrained models provided
mod <- list()
for (i in 1:length(train_ix_list_single_to2005)){
  mod[[i]] <-  fit_ridge_weighted_zonal(GPCC_na_1950_DJF, #Indecies of grid cells which are removed from the training data (options: GHCN, GPCC, DJF JAA)
                                        LENS_pr_DJF_5d00_XAX_norm, #Primary training data, also contains meta data
                                        LENS_psl_DJF_5d00_XAX_norm$X[train_ix_list_single_to2005[[i]],], #Secondary Training data, only used if multi = T
                                        train_ix_list_single_to2005[[i]], #training indecies (options: train_ix_list_single_to2005, train_ix_list_single_to2005)
                                        multi = T, # Flag if secondary training data should be included (options: T, F)
                                        target = LENS_pr_DJF_5d00_XAX_norm$Y$ZM #Matrix of the target (36 latitudes)
                                        )
}

#save model output
save(list = "mod", file  = "./models/model_as_truth/mod_DJF_prpslGPCC1950_ZM_single.RData")


