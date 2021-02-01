# 00_process_raw_eBird.R
# Author: Benjamin R. Goldstein
# Date: ???

### Script 4: Run parallelized batches of model fits.
# (1) Source helper files and read in data
# (2) Set up lists of target SSRs
# (3) Batch run models

# NOTE that this script was developed to accomodate the case that some SSRs
# are already fit, but some calls may error if none have been. However, if
# you hit an error and keep running subsequent lines, the script should work.

library(tidyverse)
library(parallel)

##### 1. Source helper files and read in data #####
source("03_fit_one_ssr_final.R")
# This might error if no fits have been run yet, but it won't stop the script
# from working
tryCatch({
  source("read_results_helper_file.R") 
}, error = function(e){})


ca_obs <- read_csv("intermediate/CA_obs_processed.csv")
ca_checklists <- read_csv("intermediate/CA_checklists_w_covariates.csv")
targets <- read_csv("intermediate/chosen_ssrs.csv") %>% 
           mutate(subregion = center, species = name_clean) %>% 
           select(subregion, species) %>% 
           distinct()


##### 2. Set up lists of target SSRs #####
tryCatch(
  {targets <- #sample_n(targets, nrow(targets), replace = F) %>% 
   targets %>%
    arrange(-subregion) %>%
    filter(!paste0(subregion, species) %in% paste0(ssrs_completed$sr, ssrs_completed$species))},
  error = function(e){}
)

targets_as_list <- list()
for (i in 1:nrow(targets)) {
  targets_as_list[[i]] <- list(species = targets$species[i],
                               subregion = targets$subregion[i])
}

targets_as_list_list <- list()
lastInd <- 0
for (i in 1:min(nCores+1, length(targets_as_list))) {
  nextInd <- min(lastInd + ceiling(length(targets_as_list) / nCores),
                 length(targets_as_list))
  targets_as_list_list[[i]] <- targets_as_list[(lastInd+1):nextInd]
  lastInd <- nextInd
}

for (i in 1:length(targets_as_list_list)) {
  if (is.null(targets_as_list_list[[i]][1][[1]])) {
    targets_as_list_list[[i]] <- 
      targets_as_list_list[[i]][2:length(targets_as_list_list[[i]])]
  }
}

##### 3. Batch run models. #####

# Set up cluster:
nCores <- 38
cl <- parallel::makeCluster(nCores)
capture <- clusterEvalQ(cl, {
  source("03_fit_one_ssr_final.R")
  ca_obs <- read_csv("intermediate/CA_obs_processed.csv")
  ca_checklists <- read_csv("intermediate/CA_checklists_w_covariates.csv")
})

# Run.
parLapply(cl = cl, X = targets_as_list_list, function(xx) {
  lapply(xx, function(x) {
    tryCatch({
      fit_one_ssr(subregion = x$subregion,
                  species = x$species,
                  checklists = ca_checklists,
                  obs = ca_obs, nCores = 0, verbose = FALSE, 
                  glmm_only = FALSE,
                  overwrite = FALSE,#c("GLMM_Pois", "GLMM_Nbin",
                  #  "Nmix_BP", "Nmix_BNB",
                  # "Nmix_BBP", "Nmix_BBNB"),
                  path_prefix = "output/onemodel_results/"
      )
    }, error = function(err) {as.character(err)}
    )
  })
})
