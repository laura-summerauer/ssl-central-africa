
################################################################################
################################################################################
###                                                                          ###
###  PREDICT SIX INDIVIDUAL HOLD-OUT REGIONS BY AFSIS SSL TOGETHER WITH THE  ###
###  FIVE REMAINING REGIONS                                                  ###
###                                                                          ###
################################################################################
################################################################################

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
  select(sample_id, region_name, tc, tn, spc) %>% 
  arrange(region_name)


# load congo ssl
cssl <- qread("out/data/cssl_assigned.qs") %>% 
  select(sample_id, region_name, holdout_region, set_type, spiking_order, tc, tn, spc) %>% 
  arrange(region_name)


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
##                 Apply Savitzky-Golay filter                 ##
#################################################################


## preprocess 
# "global" preprocessing optimization 
sg_ssl_m <- 2
sg_ssl_p <- 2
sg_ssl_w <- 17

afsis$spc_sg2 <- afsis$spc_resampled %>%
  savitzkyGolay(m = sg_ssl_m, p = sg_ssl_p, w = sg_ssl_w)

cssl$spc_sg2 <- cssl$spc_resampled %>%
  savitzkyGolay(m = sg_ssl_m, p = sg_ssl_p, w = sg_ssl_w)



##################################################################################
##  Run MBL for each hold-out regions 1) without spiking (Strategy 2) and once  ##
##  with spiking of 1–20 spiking samples (Strategy 3)                           ##
##################################################################################


# load prepared mbl function
source("R/Rhelper/my_custom_mbl.R")
source("R/Rhelper/my_custom_mbl_spiking.R")

### total carbon

## prepare cores for parallel execution
n_cores <- availableCores() - 1
clust <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(clust)
getDoParWorkers()

for (ith_region in paste0("Region", 1:6)) {

  # validation set
  region_val <- cssl[!is.na(cssl$tc) & cssl$holdout_region == ith_region & cssl$set_type == "validation_set", ]
  # spiking set
  region_spike <- cssl[!is.na(cssl$tc) & cssl$holdout_region == ith_region & cssl$set_type == "spiking_set", ]
  
  # cssl data without hold-out region
  cssl_noreg <- cssl[!cssl$holdout_region == ith_region, ]
  
  # merge both ssls (5 remaining regions of cssl with afsis ssl)
  ssl <- rbind(afsis %>% select(sample_id, region_name, tc, spc_sg2), 
               cssl_noreg %>% select(sample_id, region_name, tc, spc_sg2))
  
  
  # define numbers of neighbors to be tested (depending on region)
  if (ith_region %in% c("Region3", "Region5")) {
    ks <- seq(300, 500, by = 25)
  } else {
    ks <- seq(150, 500, by = 25)
  }
  
  
  # grouping factor for nearest neighbor validation
  ith_group_nospike <- as.factor(ssl$region_name)
  
  # run mbl without any spiking
  mbl_tc_nospike <- my_custom_mbl(Xr = ssl$spc_sg2[!is.na(ssl$tc),],
                                  Yr = ssl$tc[!is.na(ssl$tc)],
                                  Xu = region_val$spc_sg2,
                                  Yu = region_val$tc,
                                  k = ks,
                                  two_step = TRUE,
                                  cor_diss = 0.15,
                                  cor_k_range = c(50, 500),
                                  ws = 71,
                                  spike = NULL, # no spiking
                                  method = local_fit_wapls(5, 30),
                                  diss_method = "cor",
                                  group = ith_group_nospike,
                                  scale = TRUE,
                                  msc = TRUE)

  # # save object
  # qsave(mbl_tc_nospike, paste("out/final_models/nospike/", ith_region, "_TC_no-spike.qs", sep=""))


  # run mbl with 1 to 20 spiking samples 
  for (ith_spike in 1:20) {
    
    to_spike <- region_spike[region_spike$spiking_order <= ith_spike, ] # select spiking samples be order of kmeans calibration sampling order
    
    # spiked ssl
    ssl_spiked <- rbind(to_spike %>% select(sample_id, region_name, tc, spc_sg2),
                        ssl)
    
    # define indexes of spiking samples
    spiking_idx <- c(1:ith_spike) 
    
    
    # grouping factor for nearest neighbor validation
    ith_group_spike <- as.factor(ssl_spiked$region_name)
    
    
    mbl_tc <- my_custom_mbl_spike(Xr = ssl_spiked$spc_sg2[!is.na(ssl_spiked$tc),],
                              Yr = ssl_spiked$tc[!is.na(ssl_spiked$tc)],
                              Xu = region_val$spc_sg2,
                              Yu = region_val$tc,
                              k = ks, 
                              two_step = TRUE,
                              cor_diss = 0.15, 
                              cor_k_range = c(50, 500), 
                              ws = 71,
                              spike = spiking_idx, # spiking idx
                              method = local_fit_wapls(5, 30),
                              diss_method = "cor", 
                              group = ith_group_spike,
                              scale = TRUE,
                              msc = TRUE)
    
    # # save object
    # qsave(mbl_tc, paste("out/final_models/spiking/", ith_region, "_TC_nr-spike-", ith_spike, ".qs", sep=""))
      
  }
}

# stop parallel processing
registerDoSEQ()
try(stopCluster(clust))


### total nitrogen

## prepare cores for parallel execution
n_cores <- availableCores() - 1
clust <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(clust)
getDoParWorkers()

for (ith_region in paste0("Region", 1:6)) {
  
  # validation set
  region_val <- cssl[!is.na(cssl$tn) & cssl$holdout_region == ith_region & cssl$set_type == "validation_set", ]
  # spiking set
  region_spike <- cssl[!is.na(cssl$tn) & cssl$holdout_region == ith_region & cssl$set_type == "spiking_set", ]
  
  # congo ssl without hold-out region
  cssl_noreg <- cssl[!cssl$holdout_region == ith_region, ]
  
  # merge both ssls (5 remaining regions of cssl with afsis ssl)
  ssl <- rbind(afsis %>% select(sample_id, region_name, tn, spc_sg2), 
               cssl_noreg %>% select(sample_id, region_name, tn, spc_sg2))
  
  
  # grouping factor for nearest neighbor validation
  ith_group <- as.factor(ssl$region_name)
  
  # define numbers of neighbors to be tested (depending on region)
  if (ith_region %in% c("Region3", "Region5")) {
    ks <- seq(300, 500, by = 25)
  } else {
    ks <- seq(150, 500, by = 25)
  }
  
  # grouping factor for nearest neighbor validation
  ith_group_nospike <- as.factor(ssl$region_name)
  
  # run mbl without any spiking
  mbl_tn_nospike <- my_custom_mbl(Xr = ssl$spc_sg2,
                                  Yr = ssl$tn,
                                  Xu = region_val$spc_sg2,
                                  Yu = region_val$tn,
                                  k = ks, 
                                  two_step = TRUE,
                                  cor_diss = 0.15, 
                                  cor_k_range = c(50, 500), 
                                  ws = 71,
                                  spike = NULL, # no spiking
                                  method = local_fit_wapls(5, 30),
                                  diss_method = "cor", 
                                  group = ith_group_nospike,
                                  scale = TRUE,
                                  msc = TRUE)
  
  # # save object
  # qsave(mbl_tn_nospike, paste("out/final_models/nospike/", ith_region, "_TN_no-spike.qs", sep=""))
  
  # run mbl with 1 to 20 spiking samples
  for (ith_spike in 1:20) {
    
    to_spike <- region_spike[region_spike$spiking_order <= ith_spike, ] # select spiking samples be order of kmeans calibration sampling order
    
    # spiked ssl
    ssl_spiked <- rbind(to_spike %>% select(sample_id, region_name, tn, spc_sg2),
                        ssl)
    
    # define indexes of spiking samples
    spiking_idx <- c(1:ith_spike)
    
    # grouping factor for nearest neighbor validation
    ith_group_spike <- as.factor(ssl_spiked$region_name)
    
    # run mbl
    mbl_tn <- my_custom_mbl_spike(Xr = ssl_spiked$spc_sg2,
                                  Yr = ssl_spiked$tn,
                                  Xu = region_val$spc_sg2,
                                  Yu = region_val$tn,
                                  k = ks, 
                                  two_step = TRUE,
                                  cor_diss = 0.15, 
                                  cor_k_range = c(50, 500), 
                                  ws = 71,
                                  spike = spiking_idx, # spiking idx
                                  method = local_fit_wapls(5, 30),
                                  diss_method = "cor", 
                                  group = ith_group_spike,
                                  scale = TRUE,
                                  msc = TRUE)
    
    # # save object
    # qsave(mbl_tn, paste("out/final_models/spiking/", ith_region, "_TN_nr-spike-", ith_spike, ".qs", sep=""))
    
  }
}

registerDoSEQ()
try(stopCluster(clust))

