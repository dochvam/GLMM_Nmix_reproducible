
### Script 1: prepare eBird data for modeling.
# (0) Load inputs
# (1) Data are cleaned/filtered
# (2) Regions are chosen
# (3) Species are chosen
# (4) Observations are aggregated to a 50m spatial grid
# (5) Input files are produced containing only obs info for selected SSRs

library(sp)
library(raster)
library(tidyverse)
library(lubridate)
library(taxadb)

##### 0. Load inputs #####

good_CRS <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 
                 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m 
                 +no_defs")

# Raw eBird data
CA_checklists_raw <-
  read_csv("../eBird_Data/subsets-2020/state_rdbs/CA_Aug2020_checklist_info.csv")
CA_obs_raw <- 
  read_csv("../eBird_Data/subsets-2020/state_rdbs/CA_Aug2020_species_counts.csv")
CA_map <- maps::map(database = "state", regions = "CA", fill = T, plot = F)
CA_poly <- maptools::map2SpatialPolygons(
  CA_map, IDs = "california", proj4string = CRS("+proj=longlat +datum=WGS84")
) %>% spTransform(good_CRS)

subregion_radius <- 5000
subregion_buffer <- 50000

##### 1. Prepare the checklist data #####

# Clean to:
# - Only complete checklists
# - Sampling protocols P21 and P22 only
# - Breeding season (April - July)
# - Remove checklists with Xs

CA_checklists <- CA_checklists_raw %>% 
  mutate(month = lubridate::month(OBSERVATION.DATE),
         year = lubridate::year(OBSERVATION.DATE)) %>% 
  filter(ALL.SPECIES.REPORTED == 1, 
         PROTOCOL.CODE %in% c("P21", "P22"),
         month %in% 4:6,
         year %in% 2020,
         EFFORT.DISTANCE.KM <= 2)

# For stationary protocol, distance = 0
CA_checklists$EFFORT.DISTANCE.KM[
  CA_checklists$PROTOCOL.CODE == "P21"
] <- 0

CA_checklists <- CA_checklists %>% 
  filter(!is.na(OBSERVATION.DATE), !is.na(DURATION.MINUTES),
         !is.na(EFFORT.DISTANCE.KM), !is.na(NUMBER.OBSERVERS), )

# Remove checklists that reported "X" values for any species
checklist_IDs_w_Xs <- 
  unique(CA_obs_raw$SAMPLING.EVENT.IDENTIFIER[CA_obs_raw$total_count == "X"])

CA_checklists_good <- CA_checklists %>% 
  filter(!(SAMPLING.EVENT.IDENTIFIER %in% checklist_IDs_w_Xs))



##### 2. Choose regions #####
# Method:
# - Create a 0.5km grid over California
# - Count number of checklists in each grid cell

CA_checklist_points <- SpatialPointsDataFrame(
  coords = CA_checklists_good[, c("LONGITUDE", "LATITUDE")],
  data = CA_checklists_good[, c("SAMPLING.EVENT.IDENTIFIER")],
  proj4string = CRS("+proj=longlat +datum=WGS84")
) %>% 
  spTransform(good_CRS)

CA_checklists_good$x <- CA_checklist_points@coords[,1]
CA_checklists_good$y <- CA_checklist_points@coords[,2]

CA_grid <- raster(ext = extent(CA_poly), resolution = 500); values(CA_grid) <- 0
CA_checklist_points$cell <- cellFromXY(CA_grid, CA_checklist_points@coords)

cell_counts <- CA_checklist_points@data %>% 
  count(cell) %>% 
  filter(!is.na(cell))

values(CA_grid)[cell_counts$cell] <- cell_counts$n

dim <- floor(1.5 * subregion_radius / 500)
dim <- dim + (1 - dim %% 2)
neighbor_sum_mtx <- matrix(1, dim, dim)
for (i in 1:dim) for (j in 1:dim) {
  if ((i - (dim/2+0.5))^2 + (j - (dim/2+0.5))^2 <= (dim/2 - 0.5)^2) 
    neighbor_sum_mtx[i, j] <- 1
  else neighbor_sum_mtx[i, j] <- 0
}
CA_grid_neighborsums <- focal(CA_grid,
                              w = neighbor_sum_mtx, 
                              fun = sum)

CA_grid_coords <- as.data.frame(xyFromCell(CA_grid_neighborsums, 
                                           1:ncell(CA_grid_neighborsums)))
CA_grid_coords$cell <- 1:ncell(CA_grid_neighborsums)
CA_grid_coords$n_cl <- values(CA_grid_neighborsums)[CA_grid_coords$cell]
CA_grid_coords_all <- CA_grid_coords %>% 
                        filter(!is.nan(n_cl)) %>% 
                        arrange(-n_cl)


# Protocol:
# - Accept the coordinate with the highest n_cl
# - Remove all coordinates within 25 km of that one

ncenter <- 15
accepted_centers <- CA_grid_coords_all[0,]
CA_grid_coords <- CA_grid_coords_all

for (i in 1:ncenter) {
  # Accept the highest (first, since they're sorted)
  accepted_centers <- bind_rows(accepted_centers, CA_grid_coords[1,])
  CA_grid_coords <- CA_grid_coords %>% 
    filter((x - accepted_centers$x[i])^2 + (y - accepted_centers$y[i])^2 >
             subregion_buffer^2) %>%
    arrange(-n_cl)
}
accepted_centers$cID <- 1:ncenter

center_pts <- SpatialPointsDataFrame(coords = accepted_centers[, c("x", "y")],
                                     data = accepted_centers,
                                     proj = good_CRS)

# Create a polygon buffer around each point
center_buffers <- rgeos::gBuffer(center_pts, width = subregion_radius, 
                                 byid = TRUE, id = center_pts$cID,)
{
  plot(CA_poly)
  plot(center_buffers, col = "red", border = "red", add = TRUE)
}

# Identify regions associated with each checklist
CA_checklists_good$center <- sp::over(CA_checklist_points, center_buffers)$cID

CA_checklists_buffered <- CA_checklists_good %>% 
  filter(!is.na(center))

##### 3. Choose species #####

# Read in the Elton Traits database to filter out pelagic birds
if (!file.exists("data/BirdFuncDat_wITIS.csv")) {
  EltonTraits <- read_tsv("data/BirdFuncDat.txt")
  EltonTraits$species_code <- 
    taxadb::get_ids(EltonTraits$Scientific, db = "itis")
  write_csv(EltonTraits, "data/BirdFuncDat_wITIS.csv")
} else {
  EltonTraits <- read_csv("data/BirdFuncDat_wITIS.csv")
}
EltonTraits <- filter(EltonTraits, !is.na(species_code))

ca_breedingbirds <- read_csv("data/all_CA_breeding_birds.txt", col_names = "species")
ca_breedingbirds$name_clean <- gsub("['.()]", "", 
                                    gsub("[ /]", "_", ca_breedingbirds$species))
  
  
  
# Now we restrict the data further.
# Get all obs associated with those checklists
CA_obs_buffered <- CA_obs_raw %>%
  filter(SAMPLING.EVENT.IDENTIFIER %in% 
         CA_checklists_buffered$SAMPLING.EVENT.IDENTIFIER) %>% 
  filter(name_clean %in% ca_breedingbirds$name_clean) %>% 
  left_join(CA_checklists_buffered[, c("SAMPLING.EVENT.IDENTIFIER", "center")])

# Get EltonTraits info for species in our dataset
eBird_species <- CA_obs_buffered %>% 
  distinct(name_clean, SCIENTIFIC.NAME)
eBird_species$species_code <- 
  taxadb::get_ids(eBird_species$SCIENTIFIC.NAME, db = "itis")
candidate_species <- eBird_species %>%  
    filter(!is.na(species_code)) %>% 
    left_join(EltonTraits, by = "species_code") %>% 
    filter(!(IOCOrder %in% c("Charadriiformes", "Procellariiformes", 
             "Pelecaniformes", "Suliformes", "Phaethontiformes", 
             "Podicipediformes"))) %>% 
    filter(PelagicSpecialist == 0)
  

# Look at all ssr possibilities by # checklists
possible_ssrs <- CA_obs_buffered %>% 
  count(center, name_clean, SCIENTIFIC.NAME) %>% 
  left_join(accepted_centers[,c("cID", "n_cl")], by = c("center" = "cID")) %>% 
  mutate(pct = n / n_cl) %>% 
  filter(pct >= 0.1)
specs_in_many <- possible_ssrs %>% count(name_clean) %>% arrange(-n)
chosen_specs_global <- specs_in_many[1:10,]
chosen_ssrs_global <- possible_ssrs %>% 
                      filter(name_clean %in% chosen_specs_global$name_clean) %>% 
                      mutate(method = "Overall", abund_type = NA) %>% 
  select(-n_cl, -pct)

chosen_species_local_list <- list()
for (i in 1:ncenter) {
  this_site_obs <- CA_obs_buffered %>% 
    filter(center == i)
  ncls <- length(unique(this_site_obs$SAMPLING.EVENT.IDENTIFIER))
  
  this_species <- CA_obs_buffered %>% 
    filter(center == i) %>%
    filter(SCIENTIFIC.NAME %in% candidate_species$SCIENTIFIC.NAME) %>% 
    count(SCIENTIFIC.NAME, name_clean) %>% 
    arrange(-n) %>% 
    filter(n >= ncls * 0.10)
  chosen_species_local_list[[i]] <- 
    this_species[c(1:10, (nrow(this_species) - 9):nrow(this_species)), ]
  chosen_species_local_list[[i]]$abund_type <- c(rep("More", 10), rep("Less", 10))
  chosen_species_local_list[[i]]$center <- i
}
chosen_ssrs_locally <- do.call(rbind, chosen_species_local_list)
chosen_ssrs_locally$method <- "Local"

chosen_ssrs <- bind_rows(chosen_ssrs_locally, chosen_ssrs_global)
duplicate_ssrs <- count(chosen_ssrs, SCIENTIFIC.NAME, center) %>% 
                  filter(n > 1)

chosen_species <- distinct(chosen_ssrs, SCIENTIFIC.NAME, center) %>% 
                  count(SCIENTIFIC.NAME) %>% 
                  left_join(candidate_species, by = "SCIENTIFIC.NAME")

# Extract SSR observations
ssr_obs_list <- list()
for (i in 1:ncenter) {
  this_sr_specs <- chosen_ssrs %>% filter(center == i)
  
  ssr_obs_list[[i]] <- CA_obs_buffered %>% 
    dplyr::select(SAMPLING.EVENT.IDENTIFIER, center, name_clean, SCIENTIFIC.NAME, 
                  total_count) %>% 
    filter(name_clean %in% this_sr_specs$name_clean &
           center == i) %>% 
    mutate(total_count = as.numeric(total_count))
    
}
ssr_obs <- do.call(bind_rows, ssr_obs_list)



##### 4. Aggregate checklist data to grid #####
gridded_dat_list <- list()
for (i in 1:ncenter) {
  this_cls <- CA_checklists_buffered %>% 
    filter(center == i)
  
  this_extent  <- extent(center_buffers[i,])
  
  grid_raster <- raster(this_extent, resolution = 50, crs = good_CRS)
  values(grid_raster) <- 1:length(grid_raster)
  
  this_pts <- SpatialPointsDataFrame(
    coords = this_cls[, c("LONGITUDE", "LATITUDE")],
    data = this_cls[, c("SAMPLING.EVENT.IDENTIFIER")],
    proj4string = CRS("+proj=longlat +datum=WGS84")
  ) %>% 
    spTransform(good_CRS)
  this_cls$cell <- raster::extract(grid_raster, this_pts, cellnumbers=TRUE)[,2]
  
  this_cls_gll <- raster::xyFromCell(grid_raster, cell = this_cls$cell)
  
  this_cls$gx <- this_cls_gll[, "x"]
  this_cls$gy  <- this_cls_gll[, "y"]
  
  gridded_dat_list[[i]] <- this_cls
}
CA_checklists_gridded <- do.call(rbind, gridded_dat_list)






##### 5. Write files #####
# Write a file with all the checklist info we want
CA_checklists_gridded %>%
  dplyr::select(SAMPLING.EVENT.IDENTIFIER, OBSERVER.ID, center, LATITUDE, 
                LONGITUDE, gx, gy, OBSERVATION.DATE, TIME.OBSERVATIONS.STARTED, 
                DURATION.MINUTES, EFFORT.DISTANCE.KM, PROTOCOL.CODE, 
                NUMBER.OBSERVERS, year) %>% 
  write_csv("intermediate/CA_checklists_processed.csv")

write_csv(ssr_obs, "intermediate/CA_obs_processed.csv")

write_csv(chosen_species, "intermediate/chosen_species.csv")
write_csv(chosen_ssrs, "intermediate/chosen_ssrs.csv")
write_csv(accepted_centers, "intermediate/chosen_subregions.csv")  
  