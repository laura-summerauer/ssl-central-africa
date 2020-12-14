###########################################################################
###########################################################################
###                                                                     ###
###  USE SPECTRAL PROJECTIONS TO FIND THE OPTIMAL PRE-TREATMENTS        ###
###                                                                     ###
###########################################################################
###########################################################################

# load required packages
pkgs <- list("qs", "prospectr", "resemble", "magrittr")
sapply(pkgs, library, character.only = TRUE)



#################################################################
##                          Read data                          ##
#################################################################

afsis <- qread("out/data/afsis.qs")

cssl <- qread("out/data/cssl.qs")


##################################################################
##                       Resample spectra                       ##
##################################################################

afsis_wavs <- afsis$spc %>% colnames() %>% as.numeric
cssl_wavs <- cssl$spc %>% colnames() %>% as.numeric
wavenumber_label <- expression(paste('Wavenumbers [cm'^-1,']'))

# range in nm
10000000/range(afsis_wavs)
10000000/range(cssl_wavs)

# new wavenumbers (for comparison)
new_wavs <- seq(4000, 600, -16)

# resampling
afsis$spc_resampled <- afsis$spc %>% 
  resample(wav = afsis_wavs, new.wav = new_wavs)

cssl$spc_resampled <- cssl$spc %>% 
  resample(wav = cssl_wavs, new.wav = new_wavs)


# explore mean similarity
plot(x = new_wavs,
     colMeans(standardNormalVariate(cssl$spc_resampled)),
     col = "red", 
     type = "l", 
     xlim = c(4000, 500),
     xlab = wavenumber_label,
     ylab = "SNV(Absorbance)")
lines(x = new_wavs,
      colMeans(standardNormalVariate(afsis$spc_resampled)),
      col = "dodgerblue")
grid(lty = 1, col = "#8080804D")



#################################################################
##                   Reconstruction analysis                   ##
#################################################################

pls_sel <- "opc"
max_pls <- 20

# Set the parameters for savitzkyGolay first derivative
sg_1_m <- 1
sg_1_p <- 1
sg_1_w <- 9

# Set the parameters for savitzkyGolay second derivative
sg_2_m <- 2
sg_2_p <- 2
sg_2_w <- 15


# Set the parameters for savitzkyGolay second derivative
sg_3_m <- 2
sg_3_p <- 2
sg_3_w <- 17

# Set the parameters for savitzkyGolay second derivative
sg_4_m <- 2
sg_4_p <- 2
sg_4_w <- 11


results <- NULL



# load the required function
source("R/Rhelper//get_reconstruction_error.R")

# SG1 ---------------------------------------------------------------------

afsis$sg_1 <- afsis$spc_resampled %>% 
  savitzkyGolay(m = sg_1_m, p = sg_1_p, w = sg_1_w) 

cssl$sg_1 <- cssl$spc_resampled %>%
  savitzkyGolay(m = sg_1_m, p = sg_1_p, w = sg_1_w) 


sg_1 <- get_reconstruction_error(afsis$sg_1, 
                                 Yr = afsis[,c("tc", "tn")], 
                                 Xu = cssl$sg_1, 
                                 pc_selection = list(pls_sel, max_pls))
sg_1
results[["sg_1"]] <- unlist(sg_1)


# SG2 ---------------------------------------------------------------------

afsis$sg_2 <- afsis$spc_resampled %>% 
  savitzkyGolay(m = sg_2_m, p = sg_2_p, w = sg_2_w) 


cssl$sg_2 <- cssl$spc_resampled %>%
  savitzkyGolay(m = sg_2_m, p = sg_2_p, w = sg_2_w) 


sg_2 <- get_reconstruction_error(afsis$sg_2, 
                                 Yr = afsis[,c("tc", "tn")], 
                                 Xu = cssl$sg_2, 
                                 pc_selection = list(pls_sel, max_pls))
sg_2
results[["sg_2"]] <- unlist(sg_2)


# SG3 ---------------------------------------------------------------------

afsis$sg_3 <- afsis$spc_resampled %>% 
  savitzkyGolay(m = sg_3_m, p = sg_3_p, w = sg_3_w)

cssl$sg_3 <- cssl$spc_resampled %>% 
  savitzkyGolay(m = sg_3_m, p = sg_3_p, w = sg_3_w)

sg_3 <- get_reconstruction_error(afsis$sg_3, 
                                 Yr = afsis[,c("tc", "tn")], 
                                 Xu = cssl$sg_3, 
                                 pc_selection = list(pls_sel, max_pls))
sg_3
results[["sg_3"]] <- unlist(sg_3)


# SG4 ---------------------------------------------------------------------

# sg 4
afsis$sg_4 <- afsis$spc_resampled %>% 
  savitzkyGolay(m = sg_4_m, p = sg_4_p, w = sg_4_w)

cssl$sg_4 <- cssl$spc_resampled %>% 
  savitzkyGolay(m = sg_4_m, p = sg_4_p, w = sg_4_w)

sg_4 <- get_reconstruction_error(afsis$sg_4, 
                                 Yr = afsis[,c("tc", "tn")], 
                                 Xu = cssl$sg_4, 
                                 pc_selection = list(pls_sel, max_pls))
sg_4
results[["sg_4"]] <- unlist(sg_4)


# SG1 + MSC ---------------------------------------------------------------

afsis$sg_1_msc <- afsis$sg_1 %>%  
  msc()

cssl$sg_1_msc <- cssl$sg_1 %>%  
  msc(reference_spc = attr(afsis$sg_1_msc, "Reference spectrum"))

sg_1_msc <- get_reconstruction_error(afsis$sg_1_msc, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$sg_1_msc, 
                                     pc_selection = list(pls_sel, max_pls))
sg_1_msc
results[["sg_1_msc"]] <- unlist(sg_1_msc)



# SG2 + MSC ---------------------------------------------------------------

afsis$sg_2_msc <- afsis$sg_2 %>%  
  msc()

cssl$sg_2_msc <- cssl$sg_2 %>%  
  msc(reference_spc = attr(afsis$sg_2_msc, "Reference spectrum"))

sg_2_msc <- get_reconstruction_error(afsis$sg_2_msc, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$sg_2_msc, 
                                     pc_selection = list(pls_sel, max_pls))
sg_2_msc
results[["sg_2_msc"]] <- unlist(sg_2_msc)


# SG3 + MSC ---------------------------------------------------------------

afsis$sg_3_msc <- afsis$sg_3 %>%  
  msc()

cssl$sg_3_msc <- cssl$sg_3 %>%  
  msc(reference_spc = attr(afsis$sg_3_msc, "Reference spectrum"))

sg_3_msc <- get_reconstruction_error(afsis$sg_3_msc, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$sg_3_msc, 
                                     pc_selection = list(pls_sel, max_pls))
sg_3_msc
results[["sg_3_msc"]] <- unlist(sg_3_msc)


# SG4 + MSC ---------------------------------------------------------------

afsis$sg_4_msc <- afsis$sg_4 %>%  
  msc()

cssl$sg_4_msc <- cssl$sg_4 %>%  
  msc(reference_spc = attr(afsis$sg_4_msc, "Reference spectrum"))

sg_4_msc <- get_reconstruction_error(afsis$sg_4_msc, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$sg_4_msc, 
                                     pc_selection = list(pls_sel, max_pls))
sg_4_msc
results[["sg_4_msc"]] <- unlist(sg_4_msc)



# MSC + SG1 ---------------------------------------------------------------

afsis$msc_sg_1 <- afsis$spc_resampled %>%  
  msc() %>% 
  savitzkyGolay(m = sg_1_m, p = sg_1_p, w = sg_1_w) 


cssl$msc_sg_1 <- cssl$spc_resampled %>%  
  msc(reference_spc = colMeans(afsis$spc_resampled)) %>% 
  savitzkyGolay(m = sg_1_m, p = sg_1_p, w = sg_1_w) 

msc_sg_1 <- get_reconstruction_error(afsis$msc_sg_1 , 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$msc_sg_1 , 
                                     pc_selection = list(pls_sel, max_pls))

msc_sg_1
results[["msc_sg_1"]] <- unlist(msc_sg_1)


# MSC + SG2 ---------------------------------------------------------------

afsis$msc_sg_2 <- afsis$spc_resampled %>%  
  msc() %>% 
  savitzkyGolay(m = sg_2_m, p = sg_2_p, w = sg_2_w) 


cssl$msc_sg_2 <- cssl$spc_resampled %>%  
  msc(reference_spc = colMeans(afsis$spc_resampled)) %>% 
  savitzkyGolay(m = sg_2_m, p = sg_2_p, w = sg_2_w) 

msc_sg_2 <- get_reconstruction_error(afsis$msc_sg_2, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$msc_sg_2, 
                                     pc_selection = list(pls_sel, max_pls))

msc_sg_2
results[["msc_sg_2"]] <- unlist(msc_sg_2)


# MSC + SG3 ---------------------------------------------------------------

afsis$msc_sg_3 <- afsis$spc_resampled %>%  
  msc() %>% 
  savitzkyGolay(m = sg_3_m, p = sg_3_p, w = sg_3_w) 


cssl$msc_sg_3 <- cssl$spc_resampled %>%  
  msc(reference_spc = colMeans(afsis$spc_resampled)) %>% 
  savitzkyGolay(m = sg_3_m, p = sg_3_p, w = sg_3_w) 

msc_sg_3 <- get_reconstruction_error(afsis$msc_sg_3, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$msc_sg_3, 
                                     pc_selection = list(pls_sel, max_pls))

msc_sg_3
results[["msc_sg_3"]] <- unlist(msc_sg_3)



# MSC + SG4 ---------------------------------------------------------------

afsis$msc_sg_4 <- afsis$spc_resampled %>%  
  msc() %>% 
  savitzkyGolay(m = sg_4_m, p = sg_4_p, w = sg_4_w) 


cssl$msc_sg_4 <- cssl$spc_resampled %>%  
  msc(reference_spc = colMeans(afsis$spc_resampled)) %>% 
  savitzkyGolay(m = sg_4_m, p = sg_4_p, w = sg_4_w) 

msc_sg_4 <- get_reconstruction_error(afsis$msc_sg_4 , 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$msc_sg_4 , 
                                     pc_selection = list(pls_sel, max_pls))

msc_sg_4
results[["msc_sg_4"]] <- unlist(msc_sg_4)



# SG1 + detrend1 ----------------------------------------------------------

afsis$sg_1_detrend_1 <- afsis$sg_1 %>%  
  detrend(wav = as.numeric(colnames(afsis$sg_1)), p = 1)

cssl$sg_1_detrend_1 <- cssl$sg_1 %>%  
  detrend(wav = as.numeric(colnames(afsis$sg_1)), p = 1)


sg_1_detrend_1 <- get_reconstruction_error(afsis$sg_1_detrend_1, 
                                           Yr = afsis[,c("tc", "tn")], 
                                           Xu = cssl$sg_1_detrend_1, 
                                           pc_selection = list(pls_sel, max_pls))

sg_1_detrend_1
results[["sg_1_detrend_1"]] <- unlist(sg_1_detrend_1)



# SG1 + detrend2 ----------------------------------------------------------

afsis$sg_1_detrend_2 <- afsis$sg_1 %>%  
  detrend(wav = as.numeric(colnames(afsis$sg_1)), p = 2)

cssl$sg_1_detrend_2 <- cssl$sg_1 %>%  
  detrend(wav = as.numeric(colnames(afsis$sg_1)), p = 2)


sg_1_detrend_2 <- get_reconstruction_error(afsis$sg_1_detrend_2, 
                                           Yr = afsis[,c("tc", "tn")], 
                                           Xu = cssl$sg_1_detrend_2, 
                                           pc_selection = list(pls_sel, max_pls))

sg_1_detrend_2
results[["sg_1_detrend_2"]] <- unlist(sg_1_detrend_2)


# SG2 + detrend1 ----------------------------------------------------------

afsis$sg_2_detrend_1 <- afsis$sg_2 %>%  
  detrend(wav = as.numeric(colnames(afsis$sg_2)), p = 1)

cssl$sg_2_detrend_1 <- cssl$sg_2 %>%  
  detrend(wav = as.numeric(colnames(afsis$sg_2)), p = 1)


sg_2_detrend_1 <- get_reconstruction_error(afsis$sg_2_detrend_1, 
                                           Yr = afsis[,c("tc", "tn")], 
                                           Xu = cssl$sg_2_detrend_1, 
                                           pc_selection = list(pls_sel, max_pls))

sg_2_detrend_1
results[["sg_2_detrend_1"]] <- unlist(sg_2_detrend_1)

# --- 2.10 pretreat 10: sg 2 + detrend 2 ----
afsis$sg_2_detrend_2 <- afsis$sg_2 %>%  
  detrend(wav = as.numeric(colnames(afsis$sg_2)), p = 2)


# SG2 + detrend2 ----------------------------------------------------------

cssl$sg_2_detrend_2 <- cssl$sg_2 %>%  
  detrend(wav = as.numeric(colnames(afsis$sg_2)), p = 2)


sg_2_detrend_2 <- get_reconstruction_error(afsis$sg_2_detrend_2, 
                                           Yr = afsis[,c("tc", "tn")], 
                                           Xu = cssl$sg_2_detrend_2, 
                                           pc_selection = list(pls_sel, max_pls))

sg_2_detrend_2
results[["sg_2_detrend_2"]] <- unlist(sg_2_detrend_2)



# SG1 + snv ---------------------------------------------------------------

afsis$sg_1_snv <- afsis$sg_1 %>%  
  standardNormalVariate()

cssl$sg_1_snv <- cssl$sg_1 %>%  
  standardNormalVariate()

sg_1_snv <- get_reconstruction_error(afsis$sg_1_snv, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$sg_1_snv, 
                                     pc_selection = list(pls_sel, max_pls))
sg_1_snv
results[["sg_1_snv"]] <- unlist(sg_1_snv)



# SG2 + snv ---------------------------------------------------------------

afsis$sg_2_snv <- afsis$sg_2 %>%  
  standardNormalVariate()

cssl$sg_2_snv <- cssl$sg_2 %>%  
  standardNormalVariate()

sg_2_snv <- get_reconstruction_error(afsis$sg_2_snv, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$sg_2_snv, 
                                     pc_selection = list(pls_sel, max_pls))
sg_2_snv
results[["sg_2_snv"]] <- unlist(sg_2_snv)


# snv + SG1 ---------------------------------------------------------------

afsis$snv_sg_1 <- afsis$spc_resampled %>%  
  standardNormalVariate() %>% 
  savitzkyGolay(m = sg_1_m, p = sg_1_p, w = sg_1_w) 

cssl$snv_sg_1 <- cssl$spc_resampled %>%  
  standardNormalVariate() %>% 
  savitzkyGolay(m = sg_1_m, p = sg_1_p, w = sg_1_w) 

snv_sg_1 <- get_reconstruction_error(afsis$snv_sg_1 , 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$snv_sg_1 , 
                                     pc_selection = list(pls_sel, max_pls))

snv_sg_1
results[["snv_sg_1"]] <- unlist(snv_sg_1)

# --- 2.14 pretreat 14: snv + sg 2 ----
afsis$snv_sg_2 <- afsis$spc_resampled %>%  
  standardNormalVariate() %>% 
  savitzkyGolay(m = sg_2_m, p = sg_2_p, w = sg_2_w) 


# snv + SG2 ---------------------------------------------------------------

cssl$snv_sg_2 <- cssl$spc_resampled %>%  
  standardNormalVariate() %>% 
  savitzkyGolay(m = sg_2_m, p = sg_2_p, w = sg_2_w) 

snv_sg_2 <- get_reconstruction_error(afsis$snv_sg_2, 
                                     Yr = afsis[,c("tc", "tn")], 
                                     Xu = cssl$snv_sg_2, 
                                     pc_selection = list(pls_sel, max_pls))

snv_sg_2
results[["snv_sg_2"]] <- unlist(snv_sg_2)




#################################################################
##               Results of reconstruction error               ##
#################################################################

pre_processing_results <- do.call("rbind", results) %>% as.data.frame()

# filter for lowest scaled errors
pre_processing_results[pre_processing_results$scaled_error < 0.04,]

