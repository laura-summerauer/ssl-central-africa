
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
cssl <- qread("out/data/cssl.qs") %>% 
  rename(region_name = province_name) %>%
  select(sample_id, region_name, tc, tn, spc)



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



#################################################################
##        Select spiking samples using k-means sampling        ##
#################################################################

# pre-process spectra
# pre-processing optimization (see script ...)
sg_ssl_m <- 2
sg_ssl_p <- 2
sg_ssl_w <- 17

# SavitzkiGolay
afsis$spc_sg2 <- afsis$spc_resampled %>%
  savitzkyGolay(m = sg_ssl_m, p = sg_ssl_p, w = sg_ssl_w)

cssl$spc_sg2 <- cssl$spc_resampled %>%
  savitzkyGolay(m = sg_ssl_m, p = sg_ssl_p, w = sg_ssl_w)

# MSC
afsis$spc_sg2_msc <- afsis$spc_sg2 %>%
  msc()

cssl$spc_sg2_msc <- cssl$spc_sg2 %>%
  msc(reference_spc = attr(afsis$spc_sg2_msc, "Reference spectrum"))


### set spiking samples
set.seed(161020)
kmeans <- naes(X = cssl$spc_sg2_msc, k = 20, pc = 0.99, iter.max = 10, .center = TRUE, .scale = FALSE)

# define sets
spiking_set <- cssl[kmeans$model, ]
spiking_set$spiking_order <- c(1:length(kmeans$model))
spiking_set$set_type <- "spiking_set"

validation_set <- cssl[-kmeans$model, ]
validation_set$spiking_order <- NA
validation_set$set_type <- "validation_set"



#################################################################
##                 Run MBL without any spiking                 ##
#################################################################

# load prepared mbl function
source("R/Rhelper/my_custom_mbl.R")

# set size of neighbors to be tested
ks <- seq(150, 500, by = 25)


### total carbon
# tc set
cssl_val_tc <- validation_set[!is.na(validation_set$tc), ]

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
cssl_val_tn <- validation_set[!is.na(validation_set$tn), ]

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
plot_my_custom_mbl(mbl_tn_general, xlim = c(150, 500), main = paste("TN:CSSL by AfSIS"))




##################################################################
##                     Run MBL with spiking                     ##
##################################################################

# load prepared mbl function
source("R/Rhelper/my_custom_mbl_spiking.R")

## prepare cores for parallel execution
n_cores <- availableCores() - 1
clust <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(clust)
getDoParWorkers()

for (i in 1:20){

  cssl_val_tc <- validation_set[!is.na(validation_set$tc), ]
  cssl_spike_tc <- spiking_set[!is.na(spiking_set$tc) & spiking_set$spiking_order <= i, ] 

  ssl <- rbind(cssl_spike_tc %>% select(sample_id, region_name, tc, spc_sg2), 
               afsis %>% select(sample_id, region_name, tc, spc_sg2))

  # set vector with spiking indexes
  spiking_idx <- c(1:i)

  ## The neighborhood sizes to test
  ks <- seq(150, 500, by = 25)


  mbl_tc <- my_custom_mbl_spike(Xr = ssl$spc_sg2,
                          Yr = ssl$tc,
                          Xu = cssl_val_tc$spc_sg2,
                          Yu = cssl_val_tc$tc,
                          k = ks,
                          two_step = TRUE,
                          cor_diss = 0.15,
                          cor_k_range = c(50, 500),
                          ws = 71,
                          spike = spiking_idx,
                          method = local_fit_wapls(5, 30),
                          diss_method = "cor",
                          group = ssl$region_name,
                          scale = TRUE,
                          msc = TRUE)

  # qsave(mbl_tc, paste("out/final_models/spiking/", "CSSLbyAfsis_TC", "_nr-spike-", i, ".qs", sep=""))

}

registerDoSEQ()
try(stopCluster(clust))


# tn
## prepare cores for parallel execution
n_cores <- availableCores() - 1
clust <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(clust)
getDoParWorkers()

for (i in 1:20){

  cssl_val_tn <- validation_set[!is.na(validation_set$tn), ]

  cssl_spike_tn <- spiking_set[!is.na(spiking_set$tn) & spiking_set$spiking_order <= i, ]

  ssl <- rbind(cssl_spike_tn %>% select(sample_id, region_name, tn, spc_sg2), 
               afsis %>% select(sample_id, region_name, tn, spc_sg2))

  # set vector with spiking indexes
  spiking_idx <- c(1:i)

  ## The neighborhood sizes to test
  ks <- seq(150, 500, by = 25)


  mbl_tn <- my_custom_mbl_spike(Xr = ssl$spc_sg2,
                                Yr = ssl$tn,
                                Xu = cssl_val_tn$spc_sg2,
                                Yu = cssl_val_tn$tn,
                                k = ks,
                                two_step = TRUE,
                                cor_diss = 0.15,
                                cor_k_range = c(50, 500),
                                ws = 71,
                                spike = spiking_idx,
                                method = local_fit_wapls(5, 30),
                                diss_method = "cor",
                                group = ssl$region_name,
                                scale = TRUE,
                                msc = TRUE)

  # qsave(mbl_tn, paste("out/final_models/spiking/", "CSSLbyAfsis_TN", "_nr-spike-", i, ".qs", sep=""))

}

registerDoSEQ()
try(stopCluster(clust))
