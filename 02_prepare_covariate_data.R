
### Script 2: Associate eBird checklists with covariate data
# (0) Load inputs
# (1) Get data rasters
# (2) Create buffers for landcover data extraction
# (3) Assign data to checklists
# (4) Observations are aggregated to a 50m spatial grid

library(sp)
library(raster)
library(tidyverse)
library(lubridate)
library(foreign)

##### 0. Load data #####
CA_checklists <- read_csv("intermediate/CA_checklists_processed.csv")

evt1 <- raster("data/lf27569570_US_200EVT/US_200EVT/us_200evt/w001001.adf")
evt2 <- raster("data/lf48292806_US_200EVT/US_200EVT/us_200evt/w001001.adf")
evt3 <- raster("data/lf66819691_US_200EVT/US_200EVT/us_200evt/w001001.adf")
evt4 <- raster("data/lf69779486_US_200EVT/US_200EVT/us_200evt/w001001.adf")
evt_table <- 
  foreign::read.dbf("data/lf27569570_US_200EVT/US_200EVT/US_200EVT.dbf")

centers <- read_csv("intermediate/chosen_subregions.csv")
center_pts <- SpatialPointsDataFrame(coords = centers[, c("x", "y")],
                                     data = centers,
                                     proj = CRS("+proj=aea +lat_1=34 +lat_2=40.5 
                                                 +lat_0=0 +lon_0=-120 
                                                 +x_0=0 +y_0=-4000000 
                                                 +ellps=GRS80 +datum=NAD83 
                                                 +units=m +no_defs")) %>% 
              spTransform(crs(evt1))
CA_map <- maps::map(database = "state", regions = "CA", fill = T, plot = F)
CA_poly <- maptools::map2SpatialPolygons(
  CA_map, IDs = "california", proj4string = CRS("+proj=longlat +datum=WGS84")
) %>% spTransform(crs(evt1))
center_buffers <- rgeos::gBuffer(center_pts, width = 7500, 
                                 byid = TRUE, id = center_pts$cID)


r <- raster::getData("worldclim", var = "alt", res = 0.5, path = "data",
                     lon = -125, lat = 45)
r2 <- raster::getData("worldclim", var = "alt", res = 0.5, path = "data",
                      lon = -115, lat = 45)
target_raster <- raster(extent(CA_poly), resolution = 150, crs = crs(CA_poly))
r_elev <- raster::mosaic(r, r2, fun = mean) %>% 
          raster::projectRaster(target_raster)


r <- raster::getData("worldclim", var = "prec", res = 0.5, path = "data",
                     lon = -125, lat = 45)
r2 <- raster::getData("worldclim", var = "prec", res = 0.5, path = "data",
                      lon = -115, lat = 45)
r_prec <- raster::mosaic(r, r2, fun = mean) %>% 
  raster::projectRaster(target_raster)


r <- raster::getData("worldclim", var = "tmax", res = 0.5, path = "data",
                     lon = -125, lat = 45)
r2 <- raster::getData("worldclim", var = "tmax", res = 0.5, path = "data",
                      lon = -115, lat = 45)
r_tmax <- raster::mosaic(r, r2, fun = mean) %>% 
  raster::projectRaster(target_raster)


##### 1. Get rasters at each subregion #####
# Get data for each buffer

checklists_with_data <- list()

for (ctr in 1:length(center_buffers)) {
  
  ### Get the cover data needed for this raster
  needed_base_rasters <- list(evt1, evt2, evt3, evt4)[which(
    unlist(lapply(c(evt1, evt2, evt3, evt4), function (xraster) {
      !is.null(raster::intersect(extent(xraster), extent(center_buffers[ctr,])))
    })))]
  for (j in 1:length(needed_base_rasters)) {
    needed_base_rasters[[j]] <- 
      raster::crop(needed_base_rasters[[j]], extent(center_buffers[ctr,]))
  }
  if (length(needed_base_rasters) == 1) {
    cover_raster <- needed_base_rasters[[1]]
  } else if (length(needed_base_rasters) == 2) {
    cover_raster <- raster::merge(needed_base_rasters[[1]], 
                                  needed_base_rasters[[2]])
  } else if (length(needed_base_rasters) == 3) {
    cover_raster <- raster::merge(needed_base_rasters[[1]], 
                                  needed_base_rasters[[2]],
                                  needed_base_rasters[[3]])
  } else if (length(needed_base_rasters) == 4) {
    cover_raster <- raster::merge(needed_base_rasters[[1]], 
                                  needed_base_rasters[[2]],
                                  needed_base_rasters[[3]],
                                  needed_base_rasters[[4]])
  }
  
  elev_raster <- crop(r_elev, extent(center_buffers[ctr,]))
  prec_raster <- crop(r_prec, extent(center_buffers[ctr,]))
  tmax_raster <- crop(r_tmax, extent(center_buffers[ctr,]))
  values(elev_raster)[is.na(values(elev_raster))] <- 0
  values(prec_raster)[is.na(values(elev_raster))] <- 0
  values(tmax_raster)[is.na(values(elev_raster))] <- 0

  ### Get a raster of each landcover type: water, vegetation, agriculture, urban
  
  water_raster <- cover_raster
  values(water_raster) <- values(cover_raster) %in% 
    unique(evt_table$VALUE[evt_table$EVT_LF == "Water"])
  ag_raster <- cover_raster
  values(ag_raster) <- values(cover_raster) %in% 
    unique(evt_table$VALUE[evt_table$EVT_LF == "Agriculture"])
  veg_raster <- cover_raster
  values(veg_raster) <- values(cover_raster) %in% 
    unique(evt_table$VALUE[evt_table$EVT_LF %in% c("Shrub", "Herb")])
  urban_raster <- cover_raster
  values(urban_raster) <- values(cover_raster) %in% 
    unique(evt_table$VALUE[evt_table$EVT_LF == "Developed"])
  tree_raster <- cover_raster
  values(tree_raster) <- values(cover_raster) %in% 
    unique(evt_table$VALUE[evt_table$EVT_LF == "Tree"])

##### 2. Create buffers for landcover covariates #####
  
  ### Get a 600m buffer around each point of each covariate
  buffer <- 600
  buffer_mtx <- matrix(0, buffer / 30 + 1, buffer / 30 + 1)
  for (i in 1:nrow(buffer_mtx)) for (j in 1:nrow(buffer_mtx)) {
    if ((i - (nrow(buffer_mtx)/2+0.5))^2 + 
        (j - (nrow(buffer_mtx)/2+0.5))^2 <= 
        (nrow(buffer_mtx)/2 - 0.49)^2) 
      buffer_mtx[i, j] <- 1
    else buffer_mtx[i, j] <- 0
  }
  buffer_mtx <- buffer_mtx / sum(buffer_mtx)
  
  water_neigh <- focal(water_raster, buffer_mtx)
  ag_neigh <- focal(ag_raster, buffer_mtx)
  urban_neigh <- focal(urban_raster, buffer_mtx)
  tree_neigh <- focal(tree_raster, buffer_mtx)
  veg_neigh <- focal(veg_raster, buffer_mtx)
  
  ##### 3. Associate checklists w data #####
  
  checklists_this_ctr <- CA_checklists %>% 
    filter(center == center_buffers[ctr,]$cID)
  
  this_ctr_points <- SpatialPointsDataFrame(
        coords = checklists_this_ctr[, c("gx", "gy")],
        data = checklists_this_ctr,
        proj4string = CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 
                 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m 
                 +no_defs")) %>% 
    spTransform(CRSobj = crs(evt1))
  
  checklists_this_ctr$elevation <- raster::extract(elev_raster, this_ctr_points)
  checklists_this_ctr$precip <- raster::extract(prec_raster, this_ctr_points)
  checklists_this_ctr$tmax <- raster::extract(tmax_raster, this_ctr_points)
  checklists_this_ctr$pct_veg   <- raster::extract(veg_neigh, this_ctr_points)
  checklists_this_ctr$pct_water <- raster::extract(water_neigh, this_ctr_points)
  checklists_this_ctr$pct_urban <- raster::extract(urban_neigh, this_ctr_points)
  checklists_this_ctr$pct_ag    <- raster::extract(ag_neigh, this_ctr_points)
  checklists_this_ctr$pct_tree  <- raster::extract(tree_neigh, this_ctr_points)
  
  checklists_this_ctr$pct_veg_comb <- checklists_this_ctr$pct_tree +
                                      checklists_this_ctr$pct_veg
  
  checklists_with_data[[ctr]] <- checklists_this_ctr
}

 checklists_with_data_df <- do.call(bind_rows, checklists_with_data)

if (nrow(checklists_with_data_df) != nrow(CA_checklists)) stop("Lost some data...")

 
# Drop obs that are missing data
checklists_with_data_df <- checklists_with_data_df %>% 
  filter(
    !is.na(precip),
    !is.na(elevation),
    !is.na(tmax),
    !is.na(pct_veg),
    !is.na(pct_water),
    !is.na(pct_urban),
    !is.na(pct_tree),
    !is.na(pct_veg_comb),
  )
 
##### 4. Write results ######
write_csv(checklists_with_data_df, 
          "intermediate/CA_checklists_w_covariates.csv")

