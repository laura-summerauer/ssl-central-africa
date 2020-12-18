
############################################################################
############################################################################
###                                                                      ###
###  PREDICT SIX SELECTED HOLD-OUT REGIONS FROM CSSL BY AFSIS SSL        ###
###                                                                      ###
############################################################################
############################################################################

# remove objects in environment
rm(list = ls())

## load required packages
pkgs <- list("qs", "prospectr", "tidyverse", "resemble", "magrittr", "future", "doSNOW")
sapply(pkgs, library, character.only = TRUE)


#################################################################
##                          Load data                          ##
#################################################################

# load afsis ssl
afsis <- qread("out/data/afsis.qs")  %>%
  select(sample_id, region_name, tc, tn, spc)


# load congo ssl
cssl_all <- qread("out/data/cssl_assigned.qs") 


# select only validation samples/prediction set defined in 
# the R script 'R/01_define_modeling-sets.R'
cssl <- cssl_all %>% filter(set_type == "validation_set") %>% 
  select(sample_id, region_name, tc, tn, spc)

dim(cssl)
# 6 x 20 spiking samples are removed

##################################################################
##                       Resample spectra                       ##
##################################################################

afsis_wavs <- afsis$spc %>% colnames() %>% as.numeric
cssl_wavs <- cssl$spc %>% colnames() %>% as.numeric

# 4000 cm^(-1) - 600 cm^(-1) with a resolution of 16 cm^(-1)
new_wavs <- seq(4000, 600, -16)

afsis$spc_resampled <- afsis$spc %>%
  resample(wav = afsis_wavs, new.wav = new_wavs)

cssl$spc_resampled <- cssl$spc %>%
  resample(wav = cssl_wavs, new.wav = new_wavs)


################################################################################
##  Pre-process spectra using Savitzky-Golay filter; MSC will be done in the  ##
##  next step                                                                 ##
################################################################################

## preprocess 
# "global" preprocessing optimization 
sg_ssl_m <- 2
sg_ssl_p <- 2
sg_ssl_w <- 17

afsis$spc_sg2 <- afsis$spc_resampled %>%
  savitzkyGolay(m = sg_ssl_m, p = sg_ssl_p, w = sg_ssl_w)

cssl$spc_sg2 <- cssl$spc_resampled %>%
  savitzkyGolay(m = sg_ssl_m, p = sg_ssl_p, w = sg_ssl_w)


#################################################################
##                 Run MBL without any spiking                 ##
#################################################################

# load prepared mbl function
source("R/Rhelper/my_custom_mbl.R")

# set size of neighbors to be tested
ks <- seq(150, 500, by = 25)


### total carbon
# tc set; remova NAs
cssl_val_tc <- cssl[!is.na(cssl$tc), ]

## prepare cores for parallel execution
n_cores <- availableCores() - 1
clust <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(clust)
getDoParWorkers()

mbl_tc_general <- my_custom_mbl(Xr = afsis$spc_sg2,
                                Yr = afsis$tc,
                                Xu = cssl_val_tc$spc_sg2,
                                Yu = cssl_val_tc$tc,
                                k = ks,
                                two_step = TRUE,
                                cor_diss = 0.15,
                                cor_k_range = c(50, 500),
                                ws = 71,
                                spike = NULL,
                                method = local_fit_wapls(5, 30),
                                diss_method = "cor",
                                group = afsis$region_name,
                                scale = TRUE,
                                msc = TRUE)

registerDoSEQ()
try(stopCluster(clust))

mbl_tc_general
plot(mbl_tc_general)
mbl_tc_general$rmse_index_optimal
plot_my_custom_mbl(mbl_tc_general, xlim = c(150, 500), main = paste("TC: CSSL by AfSIS"))


# tn
cssl_val_tn <- cssl[!is.na(cssl$tn), ]

## prepare cores for parallel execution
n_cores <- availableCores() - 1
clust <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(clust)
getDoParWorkers()

mbl_tn_general <- my_custom_mbl(Xr = afsis$spc_sg2,
                                Yr = afsis$tn,
                                Xu = cssl_val_tn$spc_sg2,
                                Yu = cssl_val_tn$tn,
                                k = ks,
                                two_step = TRUE,
                                cor_diss = 0.15,
                                cor_k_range = c(50, 500),
                                ws = 71,
                                spike = NULL,
                                method = local_fit_wapls(5, 30),
                                diss_method = "cor",
                                group = afsis$region_name,
                                scale = TRUE,
                                msc = TRUE)

registerDoSEQ()
try(stopCluster(clust))

mbl_tn_general
plot(mbl_tn_general)
mbl_tn_general$rmse_index_optimal
plot_my_custom_mbl(mbl_tn_general, xlim = c(150, 500), main = paste("TN: CSSL by AfSIS"))


