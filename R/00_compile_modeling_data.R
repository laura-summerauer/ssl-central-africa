
################################################################################
################################################################################
###                                                                          ###
###  COMPILE BOTH SOIL SPECTRAL LIBRARIES (SSLS) FOR MODELING: CSSL (OUR     ###
###  CENTRAL AFRICAN SPECTRAL LIBRARY) AND THE CONTINENTAL SPECTRAL LIBRARY  ###
###  FROM ICRAF/AFSIS SSL                                                    ###
###                                                                          ###
################################################################################
################################################################################

# Script author: Laura Summerauer

# remove objects in environment
rm(list = ls())

# load required packages
pkgs <- c("tidyverse", "qs")
lapply(pkgs, library, character.only = TRUE)


##################################################################
##                        Load CSSL data                        ##
##################################################################

# spectra (averaged spectra for each replicates per soil sample)
cssl_spec <- read_csv("data/spectra_data/cssl_spectra_manuscript_subset/cssl_spectra.csv") %>% 
  arrange(sample_id)

cssl_wavs <- colnames(cssl_spec)[grep("^[0-9]", colnames(cssl_spec))]

# reference data 
cssl_ref <- read_csv("data/reference_data/cssl_refdata_manuscript_subset/cssl_refdata.csv") %>% 
  arrange(sample_id)

# metadata
cssl_meta <- read_csv("data/field_metadata/cssl_metadata_manuscript_subset/cssl_metadata.csv") %>% 
  arrange(sample_id)

# merge reference and metadata
cssl <- inner_join(cssl_meta, cssl_ref, by = "sample_id")

# define a new data frame to merge all three datasets
cssl$spc <- as.matrix(cssl_spec[,cssl_wavs])

# # save data
qsave(cssl, "out/data/cssl.qs")

#################################################################
##                       Load AfSIS data                       ##
#################################################################

# spectra (averaged spectra for each replicates per soil sample)
afsis_spec <- read_csv("data/spectra_data/afsis_spectra/afsis_sectra.csv") %>% 
  arrange(sample_id)

afsis_wavs <- colnames(afsis_spec)[grep("^[0-9]", colnames(afsis_spec))]

# reference data 
afsis_ref <- read_csv("data/reference_data/afsis_refdata/afsis_refdata.csv") %>% 
  arrange(sample_id)

# metadata
afsis_meta <- read_csv("data/field_metadata/afsis_metadata/afsis_metadata.csv") %>% 
  rename(sample_id = ssn, country_code = country, region_name = site) %>% 
  arrange(sample_id)

# merge reference and metadata
afsis <- inner_join(afsis_meta, afsis_ref, by = "sample_id")

# define a new data frame to merge all three datasets
afsis$spc <- as.matrix(afsis_spec[,afsis_wavs])

# save data
qsave(afsis, "out/data/afsis.qs")

