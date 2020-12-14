
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
cssl <- qread("out/data/cssl.qs") %>% 
  rename(region_name = province_name) %>%
  select(sample_id, region_name, tc, tn, spc) %>% 
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


##################################################################
##        Define validation/spiking sets for each region        ##
##################################################################

## preprocess 
# "global" preprocessing optimization 
sg_ssl_m <- 2
sg_ssl_p <- 2
sg_ssl_w <- 17

afsis$spc_sg2 <- afsis$spc_resampled %>%
  savitzkyGolay(m = sg_ssl_m, p = sg_ssl_p, w = sg_ssl_w)

cssl$spc_sg2 <- cssl$spc_resampled %>%
  savitzkyGolay(m = sg_ssl_m, p = sg_ssl_p, w = sg_ssl_w)

cssl$spc_sg2_msc <- cssl$spc_sg2 %>%
  msc()


# run a loop over every single region and select 20 validation samples
assigned_set <- list()

for (i in unique(cssl$region_name)){
  
  sub <- cssl %>% filter(!is.na(tc)) %>% filter(region_name == i) 
  no_TC <- cssl %>% filter(is.na(tc)) %>% filter(region_name == i)  # add samples with missing TC values later again
  
  # select by k-means (tendency to outlier sampling with Kennard-Stone -> see above)
  set.seed(161020)
  kmeans <- naes(X = sub$spc_sg2_msc, k = 20, pc = 0.99, iter.max = 10, .center = TRUE, .scale = FALSE)
  
  spiking_set <- sub[kmeans$model, ]
  spiking_set$spiking_order <- c(1:length(kmeans$model))
  spiking_set$set_type <- "spiking_set"
  
  validation_set <- sub[-kmeans$model, ]
  validation_set$spiking_order <- NA
  validation_set$set_type <- "validation_set"
  
  no_TC$spiking_order <- NA
  no_TC$set_type <- "validation_set"
  
  data_region <- rbind(spiking_set, validation_set, no_TC)
  assigned_set[[i]] <- data_region
  
}



cssl_assigned <- bind_rows(assigned_set) %>% 
  mutate(holdout_region = ifelse(region_name == "Haut-Katanga", "Region1", 
                                 ifelse(region_name == "Sud-Kivu", "Region2",
                                        ifelse(region_name == "Tshopo", "Region3",
                                               ifelse(region_name == "Tshuapa", "Region4",
                                                      ifelse(region_name == "Iburengerazuba", "Region5",
                                                             ifelse(region_name == "Kabarole", "Region6", "SSL"))))))) %>% 
  arrange(region_name)


# see how big the different sets are (number of samples)
cssl_assigned %>% group_by(holdout_region, set_type) %>% summarize(n())
cssl_assigned %>% colnames()

# rename column names
colnames(cssl_assigned$spc) <- colnames(cssl$spc)
colnames(cssl_assigned$spc_resampled) <- colnames(cssl$spc_resampled)
colnames(cssl_assigned$spc_sg2) <- colnames(cssl$spc_sg2)
colnames(cssl_assigned$spc_sg2_msc) <- colnames(cssl$spc_sg2_msc)

# structure of data
str(cssl_assigned)

plot(colMeans(cssl_assigned$spc))
plot(colMeans(cssl$spc))
plot(colMeans(cssl_assigned$spc_resampled))
plot(colMeans(cssl$spc_resampled))
plot(colMeans(cssl_assigned$spc_sg2_msc))
plot(colMeans(cssl$spc_sg2_msc))

# object size
format(object.size(cssl_assigned), "Mb")


##################################################################
##               Run MBL for each hold-out region               ##
##################################################################

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
  region_val <- cssl_assigned[!is.na(cssl_assigned$tc) & cssl_assigned$holdout_region == ith_region & cssl_assigned$set_type == "validation_set", ]
  # spiking set
  region_spike <- cssl_assigned[!is.na(cssl_assigned$tc) & cssl_assigned$holdout_region == ith_region & cssl_assigned$set_type == "spiking_set", ]
  
  # cssl data without hold-out region
  cssl_noreg <- cssl_assigned[!cssl_assigned$holdout_region == ith_region, ]
  
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
  region_val <- cssl_assigned[!is.na(cssl_assigned$tn) & cssl_assigned$holdout_region == ith_region & cssl_assigned$set_type == "validation_set", ]
  # spiking set
  region_spike <- cssl_assigned[!is.na(cssl_assigned$tn) & cssl_assigned$holdout_region == ith_region & cssl_assigned$set_type == "spiking_set", ]
  
  # congo ssl without hold-out region
  cssl_noreg <- cssl_assigned[!cssl_assigned$holdout_region == ith_region, ]
  
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

