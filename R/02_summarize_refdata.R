
##################################################################################
##################################################################################
###                                                                            ###
###  SUMMARIZE REFERENCE DATA FOR BOTH SOIL SPECTRAL LIBRARIES CSSL AND AFSIS  ###
###  SSL                                                                       ###
###                                                                            ###
##################################################################################
##################################################################################


# remove objects in environment
rm(list = ls())

# load required packages
pkgs <- c("tidyverse", "qs")
lapply(pkgs, library, character.only = TRUE)



#################################################################
##                          Load data                          ##
#################################################################


afsis <- qread("out/data/afsis.qs")

# we want to present the entire library, therefore we load the entire dataset here
cssl_ref <- read_csv("data/reference_data/cssl_all_refdata.csv")
cssl_meta <- read_csv("data/field_metadata/cssl_all_metadata.csv")
cssl <- inner_join(cssl_ref, cssl_meta, by = "sample_id") %>% 
  mutate(SSL = "CSSL")



#################################################################
##                Summarize CSSL reference data                ##
#################################################################


# filter for core regions 
(core_regions <- cssl %>% group_by(region_name) %>% 
  summarize(n_samples_perregion = n()) %>%
  filter(n_samples_perregion > 80) %>%  # filter regions with > 80 samples 
  select(region_name) %>% unique()
  )


cssl_core_reg <- cssl %>% 
  filter(region_name %in% core_regions$region_name) %>% 
  filter(!is.na(tc) | !is.na(tn)) %>% 
  select(sample_id, country_name, region_name, tc, tn)


cssl_rem <- cssl %>% 
  filter(!region_name %in% core_regions$region_name) %>% 
  filter(!is.na(tc) | !is.na(tn)) %>% 
  select(sample_id, country_name, region_name, tc, tn) 

chem_all <- bind_rows(cssl_core_reg, cssl_rem)


(summary_tc <- chem_all %>% 
  group_by(region_name) %>% 
  filter(!is.na(tc)) %>% 
  summarize(n_tc = n(), mean_tc = mean(tc), median_tc = median(tc), min_tc = min(tc), max_tc = max(tc))
)

(summarytn <- chem_all %>% 
  group_by(region_name) %>% 
  filter(!is.na(tn)) %>% 
  summarize(n_tn = n(), mean_tn = mean(tn), median_tn = median(tn), min_tn = min(tn), max_tn = max(tn))
)




##################################################################
##                Summarize AfSIS reference data                ##
##################################################################


afsis <- qread("out/data/afsis.qs")

(afsis_summary_tc <- afsis %>% 
  mutate(SSL = "AfSIS") %>% 
  group_by(SSL) %>% 
  filter(!is.na(tc)) %>% 
  summarize(n_tc = n(), mean_tc = mean(tc), median_tc = median(tc), min_tc = min(tc), max_tc = max(tc))
)


(afsis_summary_tn <- afsis %>% 
  mutate(SSL = "AfSIS") %>% 
  group_by(SSL) %>% 
  filter(!is.na(tn)) %>% 
  summarize(n_tn = n(), mean_tn = mean(tn), median_tn = median(tn), min_tn = min(tn), max_tn = max(tn))
)

