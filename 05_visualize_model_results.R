# 05_visualize_model_results.R
# Author: Benjamin R. Goldstein
# Date: 2/1/2021

# This script produces figures 1, 2, and 6 as well as supplemental figures 1-2
# for the manuscript.
# It also produces some tabular summaries reported in Results.

library(tidyverse)
library(gtable)
library(grid)
library(gridExtra)

source("read_results_helper_file.R")


### Make some tables ###
ssrs_completed %>% 
  count(choice) %>% 
  mutate(rate = n / nrow(ssrs_completed))

# How often was the best model an N-mixture model?
ssrs_completed %>% 
  mutate(model_choice = substr(choice, 1, 4)) %>% 
  count(model_choice)

# How often was the best N-mixture model a BB N-mix model?
nrow(ssrs_completed) - sum(ssrs_completed$best_Nmix_AIC == 
                             ssrs_completed$best_Nmix_noBB_AIC)
# How often was the best GLMM a Nbin?
nrow(ssrs_completed) - sum(ssrs_completed$best_GLMM_AIC == 
                             ssrs_completed$best_GLMM_noNB_AIC)


all_onemodels_df_completed %>% 
  filter(rank == 6, moddist == "Nmix_BP") %>% 
  nrow()

# How often was the choice decisive in favor of N-mixture models?
ssrs_completed %>% 
  mutate(magnitude = best_Nmix_AIC - best_GLMM_AIC) %>% 
  summarize(rate = sum(magnitude < -dec))


# How often was the choice decisive in favor of GLMMs?
ssrs_completed %>% 
  mutate(magnitude = best_Nmix_AIC - best_GLMM_AIC) %>% 
  summarize(rate = sum(magnitude > dec))

# How often was the choice decisive in favor of N-mixture BB models?
ssrs_completed %>% 
  mutate(magnitude = best_Nmix_AIC - best_GLMM_AIC) %>% 
  filter(grepl("BB", choice)) %>% 
  summarize(rate = sum(magnitude < -dec))


# EXCLUDE BB:
# How often was the best model an N-mixture model if you exclude beta-binomial?
ssrs_completed %>% 
  mutate(model_choice = substr(noBB_choice, 1, 4)) %>% 
  count(model_choice)
# How often was the choice decisive in favor of N-mixture models?
ssrs_completed %>% 
  mutate(magnitude = best_Nmix_noBB_AIC - best_GLMM_AIC) %>% 
  summarize(rate = sum(magnitude < -dec))
# How often was the choice decisive in favor of GLMMs?
ssrs_completed %>% 
  mutate(magnitude = best_Nmix_noBB_AIC - best_GLMM_AIC) %>% 
  summarize(rate = sum(magnitude > dec))

# EXCLUDE GLMM NB:
# How often was the best model an N-mixture model if you exclude beta-binomial?
ssrs_completed %>% 
  mutate(model_choice = substr(noGNB_choice, 1, 4)) %>% 
  count(model_choice)
# How often was the choice decisive in favor of N-mixture models?
ssrs_completed %>% 
  mutate(magnitude = best_Nmix_noBB_AIC - best_GLMM_AIC) %>% 
  summarize(rate = sum(magnitude < -dec))
# How often was the choice decisive in favor of GLMMs?
ssrs_completed %>% 
  mutate(magnitude = best_Nmix_noBB_AIC - best_GLMM_AIC) %>% 
  summarize(rate = sum(magnitude > dec))

#### Figure 1 ####
aodf_w_abd <- all_onemodels_df_completed %>% 
  left_join(ssr_site_info, by = c("species" = "name_clean",
                                  "sr" = "center"))
more_om_df <- aodf_w_abd %>% filter(abund_type == "More")
less_om_df <- aodf_w_abd %>% filter(abund_type == "Less")
wide_om_df <- aodf_w_abd %>% filter(abund_type == "Global")

ssrc_w_abd <- ssrs_completed %>% 
  left_join(ssr_site_info, by = c("species" = "name_clean",
                                  "sr" = "center"))
more_ssrc <- ssrc_w_abd %>% filter(abund_type == "More")
less_ssrc <- ssrc_w_abd %>% filter(abund_type == "Less")
wide_ssrc <- ssrc_w_abd %>% filter(abund_type == "Global")

ssrc_w_abd <- ssrc_w_abd %>% select(-abund_type, -method) %>% distinct()

wide_spec_ranks <- wide_om_df %>% 
  ggplot(aes(rank, fill = moddist)) + 
  geom_bar(position = position_fill()) +
  ggtitle(paste0("Widespread species (n = ", nrow(wide_om_df)/6, ")")) +
  scale_fill_manual(values = model_fill_colors) +
  scale_color_manual(values = model_line_colors) +
  xlab("Rank") + ylab("Rate") + labs(fill = "Model", col = "Model") +
  theme_minimal(base_size = 7) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))

more_abd_ranks <- more_om_df %>%
  ggplot(aes(rank, fill = moddist)) +
  geom_bar(position = position_fill()) +
  ggtitle(paste0("Highly reported species (n = ", nrow(more_om_df)/6, ")")) +
  scale_fill_manual(values = model_fill_colors) +
  scale_color_manual(values = model_line_colors) +
  ylab("") + xlab("Rank") +
  theme_minimal(base_size = 7) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

less_abd_ranks <- less_om_df %>%
  ggplot(aes(rank, fill = moddist)) +
  geom_bar(position = position_fill()) +
  ggtitle(paste0("Less reported species (n = ", nrow(less_om_df)/6, ")")) +
  scale_fill_manual(values = model_fill_colors) +
  scale_color_manual(values = model_line_colors) +
  ylab("") + xlab("Rank") +
  theme_minimal(base_size = 7) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") 

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
rank_legend <- get_legend(wide_spec_ranks)
wide_spec_ranks <- wide_spec_ranks + theme(legend.position = "none")

rank_laid_out <- grid.arrange(
  grobs = list(
    wide_spec_ranks, more_abd_ranks, less_abd_ranks, rank_legend
  ),
  layout_matrix = matrix(c(
    1,1,1,2,2,2,3,3,3,
    1,1,1,2,2,2,3,3,3,
    1,1,1,2,2,2,3,3,3,
    1,1,1,2,2,2,3,3,3,
    4,4,4,4,4,4,4,4,4
  ), byrow = T, nrow = 5)
)


# Delta AIC plots
daic_dat_to_plot <- ssrc_w_abd %>%
  mutate(dummy = -100) %>%
  mutate(glmm_minus_nmix = best_GLMM_AIC - best_Nmix_AIC) %>%
  mutate(ssr = as.factor(paste0(species, sr))) %>% 
  mutate(ssr = reorder(ssr, -glmm_minus_nmix)) %>% 
  arrange(ssr)

margin <- theme(plot.margin = unit(c(rep(0, 4)), "cm"))

daic_dat_to_plot$is_x_intercept <- FALSE

decisive <- T
for (i in 1:nrow(daic_dat_to_plot)) {
  if (abs(daic_dat_to_plot$glmm_minus_nmix[i]) >= dec) {
    if (!decisive) {
      daic_dat_to_plot$is_x_intercept[i] <- T
    }
    decisive <- T
  } else {
    if (decisive) {
      daic_dat_to_plot$is_x_intercept[i] <- T
    }
    decisive <- F
  }
}
x_intercepts <- daic_dat_to_plot$ssr[daic_dat_to_plot$is_x_intercept]


daic_dat_to_plot$group <- "Nmix"
daic_dat_to_plot$group[daic_dat_to_plot$glmm_minus_nmix < dec] <- "Indecisive"
daic_dat_to_plot$group[daic_dat_to_plot$glmm_minus_nmix < -dec] <- "GLMM"

mean_region_indices <- numeric(3)
mean_region_indices[1] <- median(which(daic_dat_to_plot$group == "Nmix"))
mean_region_indices[2] <- median(which(daic_dat_to_plot$group == "Indecisive"))
mean_region_indices[3] <- median(which(daic_dat_to_plot$group == "GLMM"))


label_dat <- data.frame(
  ssr = daic_dat_to_plot$ssr[mean_region_indices],
  glmm_minus_nmix = -100,
  value = c("Decisive N-mixture",
            "Indecisive",
            "Decisive GLMM")
)

daic_dots <- daic_dat_to_plot %>% 
  ggplot() +
  geom_rect(data = data.frame(xmin = -Inf, xmax = Inf,
                              ymin = 0, ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill=model_fill_colors["Nmix_BBP"], alpha=0.2) +
  geom_rect(data = data.frame(xmin = -Inf, xmax = Inf,
                              ymin = -Inf, ymax = 0),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill=model_fill_colors["GLMM_Pois"], alpha=0.2) +
  geom_point(aes(ssr, glmm_minus_nmix), cex = 1.5, pch = 1) +
  # geom_col(aes(ssr, glmm_minus_nmix, fill = choice), 
  #          col = NA, width = 0.5) +
  # facet_wrap(~abund_type) +
  # geom_hline(yintercept = 0) +
  geom_vline(xintercept = x_intercepts, col = "#474747", size = 1) +
  ylab("\u0394 AIC") +
  labs(col = "Model", fill = "Model") +
  # geom_hline(yintercept = c(dec, -dec), col = "black") +
  # geom_hline(yintercept = c(dec, -dec), col = "black") +
  scale_color_manual(values = model_line_colors) +
  scale_fill_manual(values = model_fill_colors) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  ) + margin +
  geom_text(mapping = aes(ssr, glmm_minus_nmix, label = value),
            data = label_dat)


daic_colors <- daic_dat_to_plot %>%
  ggplot() +
  geom_col(aes(x = ssr, y = dummy, fill = choice), 
           cex = 1.5, position = "dodge", width = 1) +
  # facet_wrap(~abund_type) +
  theme_void() +
  geom_hline(yintercept = 0) +
  xlab("Species-subregion dataset, sorted by \u0394 AIC") + 
  ylab("DELETE ME") + 
  labs(col = "Model", fill = "Model") +
  geom_vline(xintercept = x_intercepts, col = "#474747", size = 1) +
  scale_color_manual(values = model_line_colors) +
  scale_fill_manual(values = model_fill_colors) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(colour = "white"), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(colour = "white"),
        axis.line = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  ) +
  scale_y_continuous(limits = c(-100, 0)) +
  margin +
  scale_x_discrete(position = "top") 


daic_legend <- get_legend(daic_colors)
daic_colors <- daic_colors + theme(legend.position = "none")



daic_plot <- arrangeGrob(grobs = list(daic_colors, daic_dots, daic_legend),
                         layout_matrix = matrix(c(
                           1,1,1,1,1,3,
                           2,2,2,2,2,3,
                           2,2,2,2,2,3
                         ), byrow = T, nrow = 3), padding = 0, 
                         plot.margin = 0) #%>% 
# ggsave(filename = "output/plots/Fig2_dAIC_plot.jpg",
#        device  = "jpeg", width = 8, height = 3.5)

grid.arrange(textGrob("(a)", 0.1, just = "bottom"),
             rank_laid_out,
             textGrob("(b)", 0.1, just = "bottom"),
             daic_plot, 
             ncol = 1, heights = c(0.15, 1, 0.15, 1.5)
) %>% 
  ggsave(filename = "output/plots/Fig1.jpg",
         device = "jpeg", width = 8, height = 8)





##### Figure 2 #####

(aodf_w_abd %>% 
    ungroup() %>% 
    select(-abund_type) %>% 
    distinct() %>% 
    filter(chosen) %>%
    mutate(n_layer_overdisp = ifelse(
      moddist == "Nmix_BBNB", "2 (N-mixture BBNB)",
      ifelse(moddist %in% c("GLMM_Nbin", "Nmix_BBP", "Nmix_BNB"),
             "1 (GLMM neg bin, N-mixture BB-P and B-NB)", 
             "0 (GLMM Poisson, N-mixture B-P)")
    )) %>% 
    count(sr, n_layer_overdisp, n_cl.x) %>%
    group_by(sr) %>%
    mutate(countT= sum(n)) %>%
    group_by(n_layer_overdisp, .add = TRUE) %>%
    mutate(per=n/countT) %>% 
    ggplot(aes(log(n_cl.x), per, group = as.factor(n_layer_overdisp), 
               col = as.factor(n_layer_overdisp), 
               fill = as.factor(n_layer_overdisp),
               shape = as.factor(n_layer_overdisp))) + 
    geom_point() +
    geom_smooth(method = "lm", se = T, alpha = 0.1, cex = 0.5) +
    scale_color_manual("Number of overdispersion terms", values = unname(model_line_colors[c(4, 3 ,1)])) +
    scale_fill_manual("Number of overdispersion terms", values = unname(model_line_colors[c(4, 3 ,1)])) +
    scale_shape("Number of overdispersion terms") +
    xlab("Log number of checklists") + ylab("Rate at which category was selected by AIC") +
    theme_minimal(base_size = 8) +
    theme(plot.background = element_rect(fill = "white", color = "white")) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())) %>% 
  ggsave(filename = "output/plots/Fig2.jpg",
         device = "jpeg", width = 7.5, height = 4)





##### Figure 6 #####
library(rgdal)

radius <- 10000
srs <- read_csv("intermediate/chosen_subregions.csv")
good_CRS <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 
                 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m 
                 +no_defs")


sr_pts <- SpatialPointsDataFrame(coords = srs[, c("x", "y")],
                                 data = srs,
                                 proj = good_CRS)

# Create a polygon buffer around each point
center_buffers <- rgeos::gBuffer(sr_pts, width = radius, 
                                 byid = TRUE, id = sr_pts$cID) %>% 
  spTransform(CRS("+proj=longlat"))

label_df <- data.frame(
  name = 1:length(center_buffers),
  long = unlist(lapply(center_buffers@polygons, function(x) x@labpt[1])),
  lat = unlist(lapply(center_buffers@polygons, function(x) x@labpt[2]))
)

ca <- map_data(map = "state", region = "CA")

sr_map <- ggplot(mapping = aes(long, lat)) +
  geom_polygon(data = ca, aes(long, lat), fill = "white", col = "darkgray") +
  geom_polygon(data = center_buffers, aes(group=group),
               fill = "red", col = "red") +
  geom_text(data = label_df, aes(label = name), nudge_x = 0.25, nudge_y = 0.2) + 
  coord_fixed(1.2) +
  theme_minimal() +
  theme(plot.background = element_rect("white", colour = "white"))
ggsave(sr_map, filename = "output/plots/FigS1_sr_map.jpg", width = 5, height = 6)



#### Figure S1 ####

plot_AIC_comp_from_str <- function(ssr_str) {
  test_files <- list.files("output/onemodel_oneyear", pattern = ssr_str, 
                           full.names = TRUE)
  test_results <- lapply(test_files, function(x){
    res <- readRDS(x)
    # res$species <- "House_Finch"
    # res$subregion <- 20
    res
  })
  
  models_tried <- lapply(test_results, function(x) {
    
    nterms <- ncol(x$models_tried) - 1
    addtl <- nrow(x$coefficients) - sum(x$best_row[nterms])
    
    mt <- x$models_tried %>% 
      filter(!is.infinite(AIC))
    
    mt$code <- unlist(apply(mt, 1, function(mtrow) {
      paste0(mtrow[1:nterms], collapse = "")}))
    
    mt$complexity <- apply(mt, 1, function(thisrow) {
      thisrow <- as.numeric(thisrow[1:nterms])
      sum(thisrow) + addtl
    })
    mt
  })
  
  all_models_tried <- bind_rows(  
    models_tried[[1]] %>% mutate(model = "GLMM_Nbin"),
    models_tried[[2]] %>% mutate(model = "GLMM_Pois"),
    models_tried[[3]] %>% mutate(model = "Nmix_BBNB"),
    models_tried[[4]] %>% mutate(model = "Nmix_BBP"),
    models_tried[[5]] %>% mutate(model = "Nmix_BNB"),
    models_tried[[6]] %>% mutate(model = "Nmix_BP")
  )
  
  # models_tried$plot_order <- 
  #   paste0(str_pad(models_tried$complexity, width = 3, pad = "0"), models_tried$code)
  
  ggplot(all_models_tried) +
    # geom_line(aes(complexity, AIC, group = code)) +
    geom_point(aes(complexity, AIC, color = model), position = "jitter") +
    scale_color_manual(values= model_line_colors) +
    theme_minimal() +
    xlab("Number of parameters") +
    ggtitle(paste(gsub("_", " ", test_results[[1]]$species), 
                  "SR", test_results[[1]]$subregion))
}

AIC_plots <- list()
for (i in 1:nrow(ssrs_completed)) {
  AIC_plots[[i]] <- plot_AIC_comp_from_str(paste0("_", ssrs_completed$sr[i], "_",
                                                  ssrs_completed$species[i]))
}

### Create a set of random AIC plots
set.seed("8142")
random_indices <- sample(1:nrow(ssrs_completed), size = 6, replace = F)
random_ssrs <- paste0("_", ssrs_completed$sr, "_", ssrs_completed$species)[
  random_indices
]
AIC_plot_legend <- get_legend(AIC_plots[[1]] + theme(legend.position = "bottom"))
for (i in random_indices) {
  AIC_plots[[i]] <- AIC_plots[[i]] + theme(legend.position = "none")
}

arrangeGrob(grobs = list(
  AIC_plots[[random_indices[1]]], AIC_plots[[random_indices[2]]],
  AIC_plots[[random_indices[3]]], AIC_plots[[random_indices[4]]],
  AIC_plots[[random_indices[5]]], AIC_plots[[random_indices[6]]]),
  bottom = AIC_plot_legend, nrow = 2
) %>% 
  ggsave(filename = "output/plots/FigS1.jpg", width = 9, height = 6)



##### Figure S2 #####
spec_choice_plot <- wide_om_df %>% 
  filter(chosen) %>% 
  # count(species, moddist) %>% 
  # mutate(n = n * (-1) ^ as.numeric(moddist == "Nmix_BBP")) %>% 
  ggplot(aes(as.factor(species), sr, fill = moddist)) + 
  geom_tile() +
  ggtitle(paste0("Widespread species (n = ", nrow(wide_om_df)/6, ")")) +
  scale_fill_manual(values = model_fill_colors) +
  xlab("Species") +
  ylab("Subregion") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.x = element_blank(), 
        plot.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor.x = element_blank()) +
  coord_flip()
ggsave("output/plots/FigS1_spec_choice.jpg", device = "jpeg", spec_choice_plot,
       width = 6, height = 4)
