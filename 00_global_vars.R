list.of.packages <-
  c('raster',
    'fields',
    'ncdf4',
    'glmnet',
    'foreach',
    'ggplot2',
    'dplyr',
    'tidyr',
    'rworldmap',
    'doParallel',
    'data.table')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, 'Package'])]
if (length(new.packages)) install.packages(new.packages)

library(raster)
library(fields)
library(ncdf4)
library(glmnet)
library(foreach)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rworldmap)
library(doParallel)
library(data.table)

load("./data/DJF/LENS_pr_DJF_5d00_XAX_norm.RData")
load("./data/DJF/LENS_psl_DJF_5d00_XAX_norm.RData")
load("./data/JJA/LENS_pr_JJA_5d00_XAX_norm.RData")
load("./data/JJA/LENS_psl_JJA_5d00_XAX_norm.RData")


#=================================================
###load variables into global enviroment to be used later

#create lists for saving stuff

LENS_models <-  LENS_pr_DJF_5d00_XAX_norm$M$mod %>% unique()
ensmem_names <- LENS_pr_DJF_5d00_XAX_norm$M$ens.mem %>% unique()
latitudes <- seq(-87.5, 87.5, 5)

train_ix_list2005 <- list() 
test_ix_list2005 <- list()

train_ix_list_single <- list()
train_ix_list_single_to2005 <- list()
#=========================================================================
#================ define training indices =========================

#use these to exclude single climate model to later test on
for (i in 1:length(LENS_models)){
  train_ix_list2005[[i]] = which(LENS_pr_DJF_5d00_XAX_norm$M$mod  != LENS_models[i] &
                                   LENS_pr_DJF_5d00_XAX_norm$M$year > 1920 &
                                   LENS_pr_DJF_5d00_XAX_norm$M$year <= 2005)
  test_ix_list2005[[i]]  = which(LENS_pr_DJF_5d00_XAX_norm$M$mod  == LENS_models[i] & 
                                   LENS_pr_DJF_5d00_XAX_norm$M$year > 1920 &
                                   LENS_pr_DJF_5d00_XAX_norm$M$year <= 2005)
}

#use to train models which are to applied to observations
train_ix_list_single[[1]] = which(LENS_pr_DJF_5d00_XAX_norm$M$year > 1920)
train_ix_list_single_to2005[[1]] = which(LENS_pr_DJF_5d00_XAX_norm$M$year > 1920 &
                                           LENS_pr_DJF_5d00_XAX_norm$M$year <= 2005)


###=====================================================================
### SLECT DATA WHICH WILL BE REMOVED FROM THE TRAINING DATASET

#load obs data
load("./data/DJF/GHCN_pr_5d00_DJF_anom_XAX.RData")
load("./data/JJA/GPCC_pr_5d00_JJA_anom_XAX.RData")

load("./data/JJA/GHCN_pr_5d00_JJA_anom_XAX.RData")
load("./data/DJF/GPCC_pr_5d00_DJF_anom_XAX.RData")

#Define training masks
GHCN_TrueAll_1950_DJF <- !is.na(colMeans(GHCN_pr_5d00_DJF_anom_XAX$X[which(GHCN_pr_5d00_DJF_anom_XAX$M$year %in%  1950:2014),]))
lonlatDT5 <- get_lonlatDT(res = 5) %>% unique()
lonlatDT5[, GHCN := GHCN_TrueAll_1950_DJF]


GHCN_na_1950_JJA <- get_na(GHCN_pr_5d00_JJA_anom_XAX, 1950:2014)
GHCN_na_1950_DJF <- get_na(GHCN_pr_5d00_DJF_anom_XAX, 1950:2014)


GPCC_pr_5d00_DJF_anom_XAX$M$nobs[GPCC_pr_5d00_DJF_anom_XAX$M$nobs == 0] <- NA
GPCC_pr_5d00_JJA_anom_XAX$M$nobs[GPCC_pr_5d00_JJA_anom_XAX$M$nobs == 0] <- NA

GPCC_na_1950_DJF <- GPCC_pr_5d00_DJF_anom_XAX$M$nobs[GPCC_pr_5d00_DJF_anom_XAX$M$year %in% 1950:2014,] %>% 
  colMeans() %>% 
  is.na() %>% 
  which() 

GPCC_na_1950_JJA <- GPCC_pr_5d00_JJA_anom_XAX$M$nobs[GPCC_pr_5d00_JJA_anom_XAX$M$year %in% 1950:2014,] %>% 
  colMeans() %>% 
  is.na() %>% 
  which() 

#used for the simple average calculation

GHCN_TrueAll_1950_DJF <- !is.na(colMeans(GHCN_pr_5d00_DJF_anom_XAX$X[which(GHCN_pr_5d00_DJF_anom_XAX$M$year %in%  1950:2014),]))
GHCN_TrueAll_1950_JJA <- !is.na(colMeans(GHCN_pr_5d00_JJA_anom_XAX$X[which(GHCN_pr_5d00_JJA_anom_XAX$M$year %in%  1950:2014),]))

lonlatDT5 <- get_lonlatDT(res = 5) %>% unique()
lonlatDT5[, GHCN_DJF := GHCN_TrueAll_1950_DJF]
lonlatDT5[, GHCN_JJA := GHCN_TrueAll_1950_JJA]

