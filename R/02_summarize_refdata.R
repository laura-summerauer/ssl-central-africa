
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
cssl <- qread("out/data/cssl_assigned.qs")


#################################################################
##                Summarize CSSL reference data                ##
#################################################################


(cssl_summary_tc <- cssl %>% 
  group_by(region_name, set_type) %>% 
  filter(!is.na(tc)) %>% 
  summarize(n_tc = n(), mean_tc = mean(tc), median_tc = median(tc), min_tc = min(tc), max_tc = max(tc))
)

(cssl_summary_tn <- cssl %>% 
  group_by(region_name, set_type) %>% 
  filter(!is.na(tn)) %>% 
  summarize(n_tn = n(), mean_tn = mean(tn), median_tn = median(tn), min_tn = min(tn), max_tn = max(tn))
)


##################################################################
##                Summarize AfSIS reference data                ##
##################################################################

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

