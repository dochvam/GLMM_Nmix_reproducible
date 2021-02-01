# 07_batch_posttests.R
# Author: Benjamin R. Goldstein
# Date: 2/1/2021

# This script calls goodness-of-fit checks for all ssrs whose models
# have already been fit

# (1) Read in info and prepare cores
# (2) Batch stability checks
# (3) Batch goodness-of-fit checks
# (4) Check autocorrelation


source("06_post_tests.R")
source("read_results_helper_file.R")

### 1. Read in info and prepare cores

nCores <- 38
cl <- makeCluster(nCores)
capture <- clusterEvalQ(cl, source("06_post_tests.R"))

ssr_strs_completed <- paste0("_", ssrs_completed$sr, "_", ssrs_completed$species)

stability_files_completed <- list.files("output/posttests", pattern = "stability",
                                        full.names = TRUE)
gof_files_completed <- list.files("output/posttests", pattern = "gof",
                                        full.names = TRUE)
targets <- ssrs_completed %>%
  filter(!paste0("output/posttests/stability", ssr_strs_completed, ".csv") %in%
           stability_files_completed)

targets_as_list <- list()
ct <- 1
for (i in 1:nrow(targets)) {
  this_ssr_str <- paste0("_", targets$sr[i], "_", targets$species[i])
  if (length(targets_as_list) < i %% nCores + 1) {
    targets_as_list[[i %% nCores + 1]] <- this_ssr_str
  } else {
    targets_as_list[[i %% nCores + 1]] <- c(targets_as_list[[i %% nCores + 1]],
                                            this_ssr_str)
  }
}

##### 2. Batch stability checks #####
nbin_stabil_list <- parLapply(cl,
                              targets_as_list, 
                              function(xx) {
  for (i in 1:length(xx)) {
    x <- xx[i]
    if (!file.exists(paste0("output/posttests/stability", x, ".csv"))) {
      tryCatch({
        BP_res <- stability_check(x, moddist = "Nmix_BP", onemodel_path = "output/onemodel_oneyear")
        BP_res$moddist <- "Nmix_BP"
        BNB_res  <- stability_check(x, moddist = "Nmix_BNB", onemodel_path = "output/onemodel_oneyear")
        BNB_res$moddist <- "Nmix_BNB"
        BBP_res  <- stability_check(x, moddist = "Nmix_BBP", onemodel_path = "output/onemodel_oneyear")
        BBP_res$moddist <- "Nmix_BBP"
        BBNB_res <- stability_check(x, moddist = "Nmix_BBNB", onemodel_path = "output/onemodel_oneyear")
        BBNB_res$moddist <- "Nmix_BBNB"
  
        write_csv(
          bind_rows(BP_res, BNB_res, BBP_res, BBNB_res),
          paste0("output/posttests/stability", x, ".csv")
        )
        NA
      }, error = function(err) {as.character(err)})
    } else {
      cat("Already completed")
    }
  }
})


##### 3. Batch goodness-of-fit checks #####

gof_list <- parLapply(cl,
                      targets_as_list, 
  function(xx) {
    for (i in 1:length(xx)) {
      x <- xx[i]
      if (!file.exists(paste0("output/posttests/gof", x, ".csv"))) {
        tryCatch({
          this_result <- gof_by_ssr(x, onemodel_path = "output/onemodel_oneyear/")
         
          write_csv(
            this_result,
            paste0("output/posttests/gof", x, ".csv")
          )
          NA
        }, error = function(err) {as.character(err)})
      } else {
        NA
      }
    }
  }
)


##### 4. Check autocorrelation #####

autocorr_results_list <- list()
for (i in 1:nrow(ssrs_completed)) {
  if (i %% 10 == 0) cat(i, "\t")
  if (file.exists(paste0("output/posttests/gof_", 
                         ssrs_completed$sr[i], "_", 
                         ssrs_completed$species[i], ".csv"))) {
    autocorr_results_list[[i]] <-
      check_autocorr(species = ssrs_completed$species[i], 
                     subregion = ssrs_completed$sr[i],
                     resid_path = "output/residuals", 
                     onemodel_path = "output/onemodel_oneyear", 
                     glmm_only = FALSE)
  } else {
    autocorr_results_list[[i]] <- NULL
  }
}
autocorr_df <- do.call(rbind, autocorr_results_list)
write_csv(autocorr_df, "output/posttests/autocorr_df.csv")



##### (5). Batch checks against unmarked #####
# The following code, if uncommented, will run checks against unmarked.
# This was useful for model validation during development but is very slow.

# AICs_list <- parLapply(cl, targets_as_list, 
#    function(xx) {
#      for (i in 1:length(xx)) {
#        x <- xx[i]
#        if (!file.exists(paste0("output/posttests/umcomp", x, ".csv"))) {
#          tryCatch({
#            
#            bp_res <- readRDS(paste0("output/onemodel_oneyear/Nmix_BP", x, ".RDS"))
#            bp_comp <- refit_as_unmarked(bp_res)
#            
#            bnb_res <- readRDS(paste0("output/onemodel_oneyear/Nmix_BNB", x, ".RDS"))
#            bnb_comp <- refit_as_unmarked(bnb_res)
# 
#            saveRDS(
#              list(bp_comp = bp_comp, bnb_comp = bnb_comp),
#              paste0("output/posttests/umcomp", x, ".RDS")
#            )
#            
#            "Completed."
#          }, error = function(err) {as.character(err)})
#        } else {
#          "Already done."
#        }
#      }
#    }
# ) 
# 
# 
# umcomp_files <- list.files("output/posttests", pattern = "umcomp",
#                            full.names = TRUE)
# umcomp_values <- lapply(umcomp_files, function(x) {
#   res <- readRDS(x)
#   data.frame(
#     AIC_diff = c(res$bp_comp$refit_AIC - res$bp_comp$original_AIC,
#                  res$bnb_comp$refit_AIC - res$bnb_comp$original_AIC),
#     dist = c("BP", "BNB")
#   )
# }) %>% 
#   do.call(what = rbind)

