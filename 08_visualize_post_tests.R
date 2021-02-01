# 08_visualize_post_tests
# Author: Benjamin R. Goldstein
# Date: 2/1/2021

# This file produces figures 3 and S5-S7


##### 1. Setup #####
library(tidyverse)
source("06_post_tests.R")
source("read_results_helper_file.R")
library(grid)
library(gridExtra)

posttest_path <- "output/posttests"
color_facet_strips <- function(p) {
  g <- ggplotGrob(p)
  strips <- which(grepl('strip', g$layout$name))
  for (i in strips) {
    label_value <- 
      g$grobs[[i]]$grobs[[1]]$children[[
        which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      ]]$children[[1]]$label
    
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <-
      model_fill_colors[names(model_fill_colors) == label_value]
  }
  g
}



##### 2. Plot GOF #####

# Read in GOF
ct <- 1
gof_list <- list()
for (i in 1:nrow(ssrs_completed)) {
  this_fn <- paste0(posttest_path, "/gof_", ssrs_completed$sr[i],
                    "_", ssrs_completed$species[i], ".csv")
  if (file.exists(this_fn)) {
    gof_list[[ct]] <- read_csv(file = this_fn)
    ct <- ct + 1
  }
}
gof_df <- do.call(rbind, gof_list)
gof_df$chosenAIConly[gof_df$chosenAIConly == "GLMM_P"] <- "GLMM_Pois"
gof_df$fit_model[gof_df$fit_model == "GLMM_P"] <- "GLMM_Pois"
gof_df$chosenAIConly[gof_df$chosenAIConly == "GLMM_NB"] <- "GLMM_Nbin"
gof_df$fit_model[gof_df$fit_model == "GLMM_NB"] <- "GLMM_Nbin"


# Number of BBNB selected gof failures
gof_df %>% 
  filter(pvalue < 0.05, 
         test == "Uniformity") %>% 
  count(fit_model)

# 
gof_df %>% 
  filter(pvalue < 0.05, test == "Uniformity", chosenAIConly == "Nmix_BBNB") %>% 
  group_by(species, subregion) %>% 
  summarize(nfail = sum(pvalue < 0.05), .groups = "drop") %>% 
  count(nfail)


# Visualize GOF
gof_plot_A <- gof_df %>% 
  mutate(`Selected model` = fit_model == chosenAIConly) %>% 
  filter(test == "Uniformity",
         resid_type != "Marginal RQ") %>%
  ggplot(aes(pvalue, fill = `Selected model`)) +
  geom_histogram(bins = 20) + 
  scale_fill_manual(values = tf_colors, name = "Selected") +
  # scale_shape_manual(values = tf_shapes, name = "Selected") +
  geom_vline(xintercept = 0.05, col = "gray") +
  # xlab("P-value from KS test for normality of residuals") +
  xlab("") +
  scale_x_continuous(breaks = c(0.05, 0.5, 0.9)) + ylab("") +
  # ggtitle("A. All models") +
  theme(#axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid = element_line(color = "#ededed"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(
          size = 0.5, colour = "gray", fill = NA,
        )
  ) +
  facet_wrap(~fit_model, nrow = 2)


# gof_legend <- get_legend(gof_plot_A)
# gof_plot_A <- gof_plot_A + theme(legend.position = "none")
gof_plot_A <- color_facet_strips(gof_plot_A)
ggsave("output/plots/Fig3_GOF.jpg", gof_plot_A, device = "jpg",
       width = 8, height = 4)

# extract ssrs where BB-NB failed but was chosen
failed_nbin <- gof_df %>% 
  filter(fit_model == "GLMM_Nbin", 
         chosenAIConly == "GLMM_Nbin",
         pvalue < 0.05) %>% 
  select(species, subregion) %>% 
  mutate(ssr = paste0(subregion, species))

gof_plot_B <- gof_df %>% 
  mutate(ssr = paste0(subregion, species)) %>% 
  filter(test == "Uniformity", 
         ssr %in% failed_nbin$ssr,
         resid_type != "Marginal RQ") %>% 
  ggplot(aes(pvalue, fill = chosenAIConly==fit_model)) +
  geom_histogram(bins = 20) + 
  scale_fill_manual(values = tf_colors, name = "Selected") +
  # scale_shape_manual(values = tf_shapes, name = "Selected") +
  geom_vline(xintercept = 0.05, col = "gray") +
  # ggtitle("B. Datasets selected as GLMM Nbin, with p < 0.05 in GLMM Nbin  GoF") +
  xlab("P-value from KS test for normality of residuals") +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) + ylab("") +
  theme(#axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid = element_line(color = "#ededed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(size = 0.5, 
                                colour = "gray", 
                                fill = NA),
    legend.position = "none"
  ) +
  facet_wrap(~fit_model, nrow = 2)

gof_plot_B <- color_facet_strips(gof_plot_B)
ggsave("output/plots/FigS5_GOF.jpg", gof_plot_B, device = "jpg",
       width = 8, height = 4)




##### 3. Plot abundance value comparisons #####

accept_abd_names <- c("elevation", "precip", "tmax")


signs_agree <- list()
sigs_agree <- list()
mod_coeffs <- list()
for (i in 1:nrow(ssrs_completed)) {
  ssr_str <- paste0("_", ssrs_completed$sr[i], "_", ssrs_completed$species[i])
  
  mod_files <- list.files(onemodel_path, pattern = ssr_str, full.names = T)
  mod_abd_coeffs <- lapply(lapply(mod_files, readRDS), function(x) {
    if ("mixture" %in% names(x)) {
      x$coefficients %>% 
        filter(param %in% accept_abd_names | param == "log(p * lambda)") %>% 
        mutate(moddist = c("B-P" = "Nmix_BP",
                           "B-NB" = "Nmix_BNB",
                           "BB-P" = "Nmix_BBP",
                           "BB-NB" = "Nmix_BBNB"
        )[x$mixture])
    } else {
      x$coefficients %>% 
        filter(param %in% accept_abd_names | param == "(Intercept)") %>% 
        mutate(moddist = c("GLMM_Pois", "GLMM_Nbin")[1 + identical(x$fam, nbinom2)])
    }
  }) %>% 
    do.call(what = rbind)
  
  mod_abd_coeffs$is.sig <- sign(mod_abd_coeffs$est - 1.96 * mod_abd_coeffs$se) ==
    sign(mod_abd_coeffs$est + 1.96 * mod_abd_coeffs$se)
  
  signs_agree[[i]] <- mod_abd_coeffs %>%
    group_by(param) %>% 
    summarize(n_signs_pos = sum(sign(est) == 1), .groups = "drop") %>% 
    mutate(species = ssrs_completed$species[i], sr = ssrs_completed$sr[i])
  
  sigs_agree[[i]] <- mod_abd_coeffs %>%
    group_by(param) %>% 
    summarize(n_sigs = sum(is.sig), .groups = "drop") %>% 
    mutate(species = ssrs_completed$species[i], sr = ssrs_completed$sr[i])
  
  mod_coeffs[[i]] <- mod_abd_coeffs %>% 
    mutate(species = ssrs_completed$species[i], sr = ssrs_completed$sr[i])
}

signs_agree_df <- do.call(rbind, signs_agree)
sigs_agree_df <- do.call(rbind, sigs_agree)
mod_coeffs_df <- do.call(rbind, mod_coeffs)

mod_coeffs_df$param[mod_coeffs_df$param %in% c("(Intercept)", "log(p * lambda)")] <-
  "Intercept"

est_plots <- list()
se_plots <- list()
model_combos_list <- list()

for (j in 1:4) {
  parname <- unique(mod_coeffs_df$param)[j]
  this_coeffs <- mod_coeffs_df %>% filter(param == parname)
  
  model_combos <- as.data.frame(t(combn(unique(this_coeffs$moddist), 2))) %>% 
    rename(model1 = V1, model2 = V2) %>% 
    mutate(R2 = NA, slope = NA, slope_diff_from1 = NA, 
           int = NA, int_diff_from0 = NA)
  
  combdat <- list()
  lmlist_est <- list()
  lmlist_se <- list()
  
  model_combos_est <- model_combos %>% mutate(type = "est")
  model_combos_se <- model_combos %>% mutate(type = "se")
  
  for (i in 1:nrow(model_combos)) {
    mod1 <- model_combos[i,1]; mod2 <- model_combos[i,2]
    
    dat1 <- this_coeffs %>% 
      filter(moddist == mod1) %>% 
      select(moddist, species, sr, est, se) %>% 
      rename(model1 = moddist, est1 = est, se1 = se)
    dat2 <- this_coeffs %>% 
      filter(moddist == mod2) %>% 
      select(moddist, species, sr, est, se) %>% 
      rename(model2 = moddist, est2 = est, se2 = se)
    combdat[[i]] <- left_join(dat1, dat2) %>% 
      filter(est1 > -10 & est1 < 10 & est2 > -10 & est2 < 10) %>% 
      filter(se1 < 10 & se2 < 10)
    
    if (parname == "(Intercept)") {
      combdat[[i]] <- combdat[[i]] %>% 
        mutate(est1 = exp(est1), est2 = exp(est2))
    }
    
    # Point estimates
    lmlist_est[[i]] <- lm(est2 ~ est1, data = combdat[[i]])
    model_combos_est$R2[i] <- summary(lmlist_est[[i]])$r.squared
    model_combos_est$slope[i] <- summary(lmlist_est[[i]])$coefficients[2,1]
    model_combos_est$int[i] <- summary(lmlist_est[[i]])$coefficients[1,1]
    
    model_combos_est$slope_diff_from1[i] <- 
      sign(summary(lmlist_est[[i]])$coefficients[2,1] + 
             1.96 * summary(lmlist_est[[i]])$coefficients[2,2] - 1) ==
      sign(summary(lmlist_est[[i]])$coefficients[2,1] - 
             1.96 * summary(lmlist_est[[i]])$coefficients[2,2] - 1)
    
    model_combos_est$int_diff_from0[i] <- 
      sign(summary(lmlist_est[[i]])$coefficients[1,1] + 
             1.96 * summary(lmlist_est[[i]])$coefficients[1,2]) ==
      sign(summary(lmlist_est[[i]])$coefficients[1,1] - 
             1.96 * summary(lmlist_est[[i]])$coefficients[1,2])
    
    # Standard errors
    lmlist_se[[i]] <- lm(se2 ~ se1, data = combdat[[i]])
    model_combos_se$R2[i] <- summary(lmlist_se[[i]])$r.squared
    model_combos_se$slope[i] <- summary(lmlist_se[[i]])$coefficients[2,1]
    model_combos_se$int[i] <- summary(lmlist_se[[i]])$coefficients[1,1]
    
    model_combos_se$slope_diff_from1[i] <- 
      sign(summary(lmlist_se[[i]])$coefficients[2,1] + 
             1.96 * summary(lmlist_se[[i]])$coefficients[2,2] - 1) ==
      sign(summary(lmlist_se[[i]])$coefficients[2,1] - 
             1.96 * summary(lmlist_se[[i]])$coefficients[2,2] - 1)
    
    model_combos_se$int_diff_from0[i] <- 
      sign(summary(lmlist_se[[i]])$coefficients[1,1] + 
             1.96 * summary(lmlist_se[[i]])$coefficients[1,2]) ==
      sign(summary(lmlist_se[[i]])$coefficients[1,1] - 
             1.96 * summary(lmlist_se[[i]])$coefficients[1,2])
    
  }
  combdat_df <- do.call(rbind, combdat)
  
  est_plots[[j]] <- ggplot(combdat_df, aes(est1, est2)) + 
    geom_point(cex = 0.7) +
    facet_grid(model1~model2) +
    geom_smooth(method='lm', formula = y~x)+
    geom_abline(slope = 1, intercept = 0, col = "lightgray") +
    theme_minimal() +
    scale_x_continuous(breaks = c(-8, -4, 0, 4, 8)) +
    scale_y_continuous(breaks = c(-8, -4, 0, 4, 8)) +
    theme(panel.grid.minor = element_blank()) +
    coord_fixed() +
    xlab("point estimate") + 
    ylab("point estimate") +
    ggtitle(paste0(c("A. ", "B. ", "C. ", "D. ")[j], parname))
  
  ggsave(filename = paste0("output/plots/FigS6_est_", parname, ".jpg"), 
         plot = est_plots[[j]], width = 5, height = 5)
  
  se_plots[[j]] <- ggplot(combdat_df, aes(se1, se2)) + 
    geom_point(cex = 0.7) +
    facet_grid(model1~model2) +
    geom_smooth(method='lm', formula = y~x)+
    geom_abline(slope = 1, intercept = 0, col = "lightgray") +
    theme_minimal() +
    scale_x_continuous(breaks = c(-8, -4, 0, 4, 8)) +
    scale_y_continuous(breaks = c(-8, -4, 0, 4, 8)) +
    theme(panel.grid.minor = element_blank()) +
    coord_fixed() +
    xlab("standard error") + 
    ylab("standard error") +
    ggtitle(paste0(c("A. ", "B. ", "C. ", "D. ")[j], parname))
  
  ggsave(filename = paste0("output/plots/FigS7_se_", parname, ".jpg"), 
         plot = se_plots[[j]], width = 5, height = 5)
  
  
  model_combos_list[[j]] <- bind_rows(model_combos_est, model_combos_se) %>% 
    mutate(param = parname)
}