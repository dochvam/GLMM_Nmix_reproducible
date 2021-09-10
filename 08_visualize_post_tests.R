# 08_visualize_post_tests.R
# Author: Benjamin R. Goldstein
# Date: 2/1/2021

### Script 8: Visualize posttests
# This sccript contains code for generating figures 3-5 and supplemental figure
# 3. It also creates summaries of AIC stability.
##### 1. Setup #####
library(tidyverse)
source("06_post_tests.R")
source("read_results_helper_file.R")
library(grid)
library(gridExtra)
library(glmmTMB)

posttest_path <- "output/posttests_oneyear/"
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


# Collect files by looping over ssr info (saves string extraction of info)
stability_basic_list <- list()
stability_list <- list()
ct <- 0
for (i in 1:nrow(ssrs_completed)) {
  ssr_str <- paste0("_", ssrs_completed$sr[i], "_", ssrs_completed$species[i])
  thisfile <- paste0(posttest_path, "/stability", ssr_str, ".csv")
  if (file.exists(thisfile)) {
    ct <- ct + 1
    this_result <- read_csv(thisfile)
    this_result$ssr_str <- ssr_str
    this_result$species <- ssrs_completed$species[i]
    this_result$sr <- ssrs_completed$sr[i]
    this_result$choice <- ssrs_completed$choice[i]
    
    
    
    dlist <- list()
    basic_summary_df <- data.frame(
      species = ssrs_completed$species[i],
      sr = ssrs_completed$sr[i],
      choice = ssrs_completed$choice[i],
      moddist = names(model_fill_colors)[1:4],
      dlpl_topK = NA,
      dlpol_topK = NA,
      dAIC_topK = NA
    )
    
    
    for (j in 1:4) {
      tempdf <- this_result %>% 
        filter(moddist == names(model_fill_colors[j]))
      
      this_ddf <- tempdf
      this_ddf$est <- this_ddf$est - this_ddf$est[this_ddf$K == min(this_ddf$K)]
      
      dlist[[j]] <- this_ddf
      
      basic_summary_df$dlpl_topK[j] <- 
        this_ddf$est[this_ddf$param == "log(p * lambda)" &
                       this_ddf$K == max(this_ddf$K)] - 
        this_ddf$est[this_ddf$param == "log(p * lambda)" &
                       this_ddf$K == max(this_ddf$K[this_ddf$K != max(this_ddf$K)])]
      basic_summary_df$dlpol_topK[j] <- 
        this_ddf$est[this_ddf$param == "log(p / lambda)" &
                       this_ddf$K == max(this_ddf$K)] - 
        this_ddf$est[this_ddf$param == "log(p / lambda)" &
                       this_ddf$K == max(this_ddf$K[this_ddf$K != max(this_ddf$K)])]
      basic_summary_df$dAIC_topK[j] <- 
        this_ddf$est[this_ddf$param == "AIC" &
                       this_ddf$K == max(this_ddf$K)] - 
        this_ddf$est[this_ddf$param == "AIC" &
                       this_ddf$K == max(this_ddf$K[this_ddf$K != max(this_ddf$K)])]
      
    }
    ddf <- do.call(rbind, dlist)
    
    stability_list[[ct]] <- ddf
    stability_basic_list[[ct]] <- basic_summary_df
  }
}
stability_df <- do.call(rbind, stability_list)

stability_summary_df <- do.call(rbind, stability_basic_list) %>% 
  left_join(all_onemodels_df[, c("species", "sr", "moddist", "num_nan")])


##### Stability stats #####
stability_summary_df <- stability_summary_df %>% 
  mutate(`Stable in AIC` = abs(dAIC_topK) < 0.1) %>% 
  mutate(stable_lpl = abs(dlpl_topK) < 0.1) %>% 
  mutate(stable_lpol = abs(dlpol_topK) < 0.1)

# Rates of stability in AIC
stability_summary_df %>% 
  group_by(moddist, `Stable in AIC`) %>% 
  filter(!`Stable in AIC`) %>% 
  count() %>% 
  mutate(rate = n / nrow(ssrs_completed))

# Rates of stability in log(p*lambda)
stability_summary_df %>% 
  group_by(moddist, stable_lpl) %>% 
  filter(!stable_lpl) %>% 
  count() %>% 
  mutate(rate = n / nrow(ssrs_completed))

# Rates of stability in log(p/lambda)
stability_summary_df %>% 
  group_by(moddist, stable_lpol) %>% 
  filter(!stable_lpol) %>% 
  count() %>% 
  mutate(rate = n / nrow(ssrs_completed))



stability_summary_df %>% 
  group_by(`Stable in AIC`, stable_lpl, stable_lpol) %>% 
  count()


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
         resid_type != "Marginal RQ",
         `Selected model`) %>%
  ggplot(aes(pvalue)) +
  geom_histogram(bins = 20, fill = "black") + 
  geom_vline(xintercept = 0.05, col = "gray") +
  xlab("") +
  scale_x_continuous(breaks = c(0.05, 0.5, 0.9)) + ylab("") +
  theme(
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


gof_plot_A <- color_facet_strips(gof_plot_A)
ggsave("output/plots/Fig3.jpg", gof_plot_A, device = "jpg",
       width = 8, height = 4)

failed <- gof_df %>% 
  filter(pvalue < 0.05, 
         test == "Uniformity",
         fit_model %in% c("GLMM_Nbin", "GLMM_Pois"), 
         chosenAIConly == fit_model) %>% 
  mutate(ssr = paste0(subregion, species))

gof_plot_B <- gof_df %>% 
  mutate(ssr = paste0(subregion, species)) %>% 
  filter(test == "Uniformity", 
         ssr %in% failed$ssr,
         resid_type != "Marginal RQ") %>% 
  ggplot(aes(pvalue, fill = chosenAIConly==fit_model)) +
  geom_histogram(bins = 20) + 
  scale_fill_manual(values = tf_colors, name = "Selected") +
  geom_vline(xintercept = 0.05, col = "gray") +
  xlab("P-value from KS test for normality of residuals") +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) + ylab("") +
  theme(
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
ggsave("output/plots/FigS3.jpg", gof_plot_B, device = "jpg",
       width = 8, height = 4)




##### 3. Plot abundance value comparisons #####

accept_abd_names <- c("elevation")


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

for (j in 1:2) {
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
  }
  combdat_df <- do.call(rbind, combdat) %>% 
    mutate(cat = paste0(model1, " - ", model2))
  
  combdat_df$est_diff <- log(combdat_df$est1) - log(combdat_df$est2)
  combdat_df$se_log_ratio <- log(combdat_df$se1) - log(combdat_df$se2)
  
  if (parname == "elevation") parname <- "Elevation"

  est_plots[[j]] <- combdat_df %>% 
    ggplot(aes(reorder(cat, est_diff, mean, na.rm = T), est_diff)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 0) +
    theme_minimal() +
    scale_y_continuous(limits = c(-4.5, 4)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_rect(fill= "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    coord_flip() +
    xlab("") + 
    ylab("") +
    ggtitle(paste0(c("(a) ", "(b) ", "(c) ", "(d) ")[j], parname))
  
  model_combos_list[[j]] <- bind_rows(model_combos_est, model_combos_se) %>% 
    mutate(param = parname)
}


grid.arrange(grobs = est_plots, nrow=1, bottom = "Difference in estimates") %>%
  ggsave(filename = "output/plots/Fig4.jpg",
         width = 8, height = 7)



this_coeffs <- mod_coeffs_df %>% filter(param %in% c("Intercept", "elevation"))

ylim1 <- boxplot.stats(this_coeffs$se)$stats[c(1, 5)]

(this_coeffs %>% 
    mutate(param = ifelse(param == "Intercept", "(a) Intercept",
                          ifelse(param == "elevation", "(b) Elevation",
                                 NA))) %>% 
    ggplot() +
    geom_boxplot(aes(moddist, se, fill = moddist), outlier.shape = NA) + 
    ylim(c(0, 2)) +
    scale_fill_manual("", values = model_fill_colors) +
    xlab("") + ylab("Estimated standard error") +
    theme_minimal() +
    coord_flip() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = "white")) +
    facet_wrap(~factor(param, levels = c("(a) Intercept", "(b) Elevation")))) %>% 
  ggsave(filename = "output/plots/Fig5.jpg", 
         width = 8, height = 4)



autocorr_df <- read_csv("output/posttests_oneyear/autocorr_df.csv")
mean(autocorr_df$p.value < 0.05, na.rm = T)
