###################################################################################
###################################################################################
###                                                                             ###
###  DEFINE SETS WHICH WILL BE USED FOR THE SPECTRAL MODELING: SELECT REGIONS,  ###
###  DEFINE PREDICTION SETS AND SPIKING SETS                                    ###
###                                                                             ###
###################################################################################
###################################################################################


# remove objects in environment
rm(list = ls())

## load required packages
pkgs <- list("qs", "prospectr", "tidyverse", "resemble", "magrittr", "future", "doSNOW")
sapply(pkgs, library, character.only = TRUE)


#################################################################
##                          Load data                          ##
#################################################################


# load congo ssl
cssl <- qread("out/data/cssl.qs") %>% 
  rename(region_name = province_name) %>%
  select(sample_id, region_name, tc, tn, spc) %>% 
  arrange(region_name)


##################################################################
##                       Resample spectra                       ##
##################################################################

cssl_wavs <- cssl$spc %>% colnames() %>% as.numeric

# 4000 cm^(-1) - 600 cm^(-1) with a resolution of 16 cm^(-1)
new_wavs <- seq(4000, 600, -16)

# resample 
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

# save object
qsave(cssl_assigned, "out/data/cssl_assigned.qs")

