# read_results_helper_file.R
# Author: Benjamin R. Goldstein
# Date: 2/1/2021

# Sourcing this file reads in summary information for all completed
# models. It's needed by scripts 4-8.

library(tidyverse)

logit <- function(x) log(x / (1-x))
dec <- 4 # decisiveness threshold

results_folder <- onemodel_path <- "output/onemodel_results"
num__s <- str_count(results_folder, "_")

ssr_site_info <- read_csv("intermediate/chosen_ssrs.csv") %>% 
  rename("ssr_n" = "n") %>% 
  left_join(read_csv("intermediate/chosen_subregions.csv"),
            by = c("center" = "cID")) %>% 
  mutate(report_rate = ssr_n / n_cl)
ssr_site_info$abund_type[is.na(ssr_site_info$abund_type)] <- "Global"

# Get list of files
all_onemodels <- list.files(results_folder, full.names = TRUE)
all_onemodels_df <- do.call(rbind,
    lapply(all_onemodels, function(x) {
      temp <- readRDS(x)
      # cat(x, "\n")
      split_file <- unlist(strsplit(unlist(strsplit(x, "_")), "/"))
      
      log_theta <- temp$coefficients$est[temp$coefficients$param == "Log theta"]
      if (length(log_theta) == 0) log_theta <- NA
      
      log_s <- temp$coefficients$est[temp$coefficients$param == "Log s"]
      if (length(log_s) == 0) log_s <- NA
      
      if (split_file[3 + num__s] == "GLMM" && split_file[4 + num__s] == "Nbin") {
        nb_phi <- summary(temp$fit)$sigma
      } else {
        nb_phi <- NA
      }
      
      data.frame(model = split_file[3 + num__s],
                 distribution = split_file[4 + num__s],
                 moddist = paste0(split_file[3 + num__s], "_", split_file[4 + num__s]),
                 num_inf = sum(is.infinite(temp$models_tried$AIC)),
                 num_nan = sum(is.nan(temp$coefficients$se)),
                 num_params = sum(temp$best_row[1:(ncol(temp$best_row) - 1)]),
                 time_wrote = file.info(x)$mtime,
                 sr = temp$subregion,
                 species = temp$species,
                 glmm_nb_phi = nb_phi,
                 log_theta = log_theta, 
                 log_s = log_s,
                 time_taken = as.numeric(temp$time_taken, units = "mins"),
                 AIC = temp$AIC,
                 fn = x)
    }))

 ssrs_completed <- all_onemodels_df %>% 
  count(sr, species) %>% 
  filter(n == 6) %>% 
  select(-n) %>% 
  left_join(ssr_site_info %>% 
              select(-method, -abund_type) %>% distinct(), 
            by = c("species" = "name_clean", "sr" = "center")) %>% 
  left_join(all_onemodels_df %>% 
              distinct(sr, species, glmm_nb_phi) %>% 
              filter(!is.na(glmm_nb_phi)))
 
all_onemodels_df_completed <- all_onemodels_df %>% 
  filter(paste0(species, sr) %in% paste0(ssrs_completed$species, 
                                         ssrs_completed$sr)) %>% 
  left_join(ssrs_completed, by = c("species", "sr")) %>%
  group_by(sr, species) %>% 
  arrange(sr, species) %>% 
  mutate(chosen = AIC == min(AIC), rank = as.factor(order(order(AIC))))

all_onemodels_df_incomplete <- all_onemodels_df %>% 
  filter(!(paste0(species, sr) %in% paste0(ssrs_completed$species, 
                                         ssrs_completed$sr)))


for (i in 1:nrow(ssrs_completed)) {
  thismods <- all_onemodels_df_completed %>% 
    filter(species == ssrs_completed$species[i], sr == ssrs_completed$sr[i])
  
  ssrs_completed$choice[i] <- thismods$moddist[which.min(thismods$AIC)]
  
  unchosen <- thismods %>% filter(!moddist == ssrs_completed$choice[i])
  
  ssrs_completed$nextbest[i] <- unchosen$moddist[which.min(unchosen$AIC)]
  ssrs_completed$diff_to_next[i] <- thismods$AIC[which.min(thismods$AIC)] - 
    unchosen$AIC[which.min(unchosen$AIC)]
  
  unchosen <- unchosen %>% filter(!moddist == ssrs_completed$nextbest[i])
  
  ssrs_completed$thirdbest[i] <- unchosen$moddist[which.min(unchosen$AIC)]
  
  ssrs_completed$diff_to_third[i] <- thismods$AIC[which.min(thismods$AIC)] - 
    unchosen$AIC[which.min(unchosen$AIC)]
  
  this_nmix <- thismods %>% 
    filter(model == "Nmix")
  this_nmix_nobb <- thismods %>% 
    filter(moddist %in% c("Nmix_BP", "Nmix_BNB"))
  this_glmm <- thismods %>% 
    filter(model == "GLMM")
  this_glmm_nonb <- thismods %>% 
    filter(moddist == "GLMM_Pois")
  
  this_nobb <- bind_rows(this_nmix_nobb, this_glmm)
  this_nonb <- bind_rows(this_nmix, this_glmm_nonb)
  
  ssrs_completed$best_Nmix_AIC[i] <- min(this_nmix$AIC)
  ssrs_completed$best_Nmix_noBB_AIC[i] <- min(this_nmix_nobb$AIC)
  ssrs_completed$noBB_choice[i] <- this_nobb$moddist[which.min(this_nobb$AIC)]
  ssrs_completed$best_GLMM_AIC[i] <- min(this_glmm$AIC)
  ssrs_completed$best_GLMM_noNB_AIC[i] <- min(this_nonb$AIC)
  ssrs_completed$noGNB_choice[i] <- this_nonb$moddist[which.min(this_nonb$AIC)]
  
  ssrs_completed$strength_cat[i] <- 
    if (abs(ssrs_completed$diff_to_next[i]) < dec) "Indecisive" 
  else "Decisive"
}

ssrs_completed$choicecat <- paste0(ssrs_completed$strength_cat, " ",
                                   substr(ssrs_completed$choice, 1, 4))

to_exclude <- all_onemodels_df %>% filter(glmm_nb_phi > 1e12)
ssrs_completed <- ssrs_completed %>% 
  filter(!(paste0(sr, species) %in% 
           paste0(to_exclude$sr, to_exclude$species)))
all_onemodels_df <- all_onemodels_df %>% 
  filter(!(paste0(sr, species) %in% 
             paste0(to_exclude$sr, to_exclude$species)))
all_onemodels_df_completed <- all_onemodels_df_completed %>% 
  filter(!(paste0(sr, species) %in% 
             paste0(to_exclude$sr, to_exclude$species)))


abund_types <- ssr_site_info %>%
  filter(!(abund_type == "Global")) %>%
  select(center, name_clean, abund_type, n_cl, report_rate)
# ssrs_completed <- left_join(ssrs_completed, abund_types,
#                             by = c("species" = "name_clean", "sr" = "center", 
#                                    "n_cl", "report_rate"))

# In order: GLMM negbin, GLMM pois, Nmix negbin, Nmix pois
model_fill_colors <- c(
  "Nmix_BBNB" = "#1F78B4",
  "Nmix_BBP"  = "#A6CEE3",
  "Nmix_BNB"  = "#33A02C",
  "Nmix_BP"   = "#B2DF8A",
  "GLMM_Nbin" = "#FF7F00", 
  "GLMM_Pois" = "#FDBF6F"
)

model_line_colors <- c(
  "Nmix_BBNB" = "#1F78B4",
  "Nmix_BBP"  = "#70b6db",
  "Nmix_BNB"  = "#23781e",
  "Nmix_BP"   = "#76c92a",
  "GLMM_Nbin" = "#FF7F00", 
  "GLMM_Pois" = "#FDBF6F"
)

tf_colors <- c("#f00000", "#000000")
tf_shapes <- c(1, 19)

