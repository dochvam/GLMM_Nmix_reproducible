# Script 5: Visualize results of Script 4 (paper figures 1 and 2)

library(tidyverse)
library(ggalluvial)
library(gtable)
library(grid)
library(gridExtra)

source("read_results_helper_file.R")

##### 1. Automate some summary statistics ######

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

##### 2. Ranks plot ######
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
#### Plot chosen models by AIC ####
  

wide_spec_ranks <- wide_om_df %>% 
  ggplot(aes(rank, fill = moddist)) + 
  geom_bar(position = position_fill()) +
  ggtitle(paste0("Widespread species (n = ", nrow(wide_om_df)/6, ")")) +
  scale_fill_manual(values = model_fill_colors) +
  scale_color_manual(values = model_line_colors) +
  xlab("Rank") + ylab("Rate") + labs(fill = "Model", col = "Model") +
  theme_minimal(base_size = 8) +
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
  theme_minimal(base_size = 8) +
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
  theme_minimal(base_size = 8) +
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
ggsave(filename = "output/plots/Fig1_AIC_ranks.jpg", plot = rank_laid_out, 
       device = "jpeg", width = 8, height = 3.5, dpi = 900)



##### 3. Delta AIC plots #####
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



arrangeGrob(grobs = list(daic_colors, daic_dots, daic_legend),
            layout_matrix = matrix(c(
              1,1,1,1,1,3,
              2,2,2,2,2,3,
              2,2,2,2,2,3
            ), byrow = T, nrow = 3), padding = 0, 
            plot.margin = 0) %>% 
  ggsave(filename = "output/plots/Fig2_dAIC_plot.jpg",
         device  = "jpeg", width = 8, height = 3.5)



##### 4. Plot subregion buffers #####
library(rgdal)

radius <- 5000
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
  name = 1:15,
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
  theme_minimal()
ggsave(sr_map, filename = "output/plots/FigS1_sr_map.jpg", width = 5, height = 6)

