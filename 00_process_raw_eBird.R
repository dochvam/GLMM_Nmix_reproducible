# 00_process_raw_eBird.R
# Author: Benjamin R. Goldstein
# Date: 2/1/2021

# This file uses the package auk to transform the eBird raw dataset into a
# two-file relational database storing observation counts and checklist
# metadata. The user needs to download the eBird full release dataset as a .txt
# and indicate its filepath in the "input file" slot.

library(auk)
library(tidyverse)

input_file <- stop("Put your input file here.")
output_file <- "data/ebd_CA_auk.txt"

# With the full repository, this takes a very long time.
ebird_data <- input_file %>% 
  auk_ebd() %>% 
  auk_state(state = "US-CA") %>% 
  auk_filter(file = output_file)





# Because of checklist metadata redundancy, this file is still way too big to
# read into memory. I'll do it in chunks. First, I split the file by lines.

dir.create("data/state_splits")
system(paste("split -l 2000000", output_file, "data/state_splits/CA_"))

get_sum <- function(xs) {
  if (any(xs == "X")) "X"
  else as.character(sum(as.numeric(xs)))
}

verbose <- T
state_code <- "CA"

if (verbose) cat("Starting", state_code, "\n")

this_state_files <- list.files("data/state_splits", 
                               pattern = state_code, full.names = TRUE)

if (verbose) cat("...Reading files...\n")
suppressWarnings(
  suppressMessages(
    state_data_list <- lapply(this_state_files[1:200], function(x) {
      read_csv(x, progress = F, 
               col_types = list(EFFORT.AREA.HA = col_double(),
                                SAMPLING.EVENT.IDENTIFIER = col_character(),
                                OBSERVATION.DATE = col_date(),
                                TIME.OBSERVATIONS.STARTED = col_time(),
                                OBSERVER.ID = col_character(),
                                DURATION.MINUTES = col_double(),
                                EFFORT.DISTANCE.KM = col_double())) %>% 
        distinct(SAMPLING.EVENT.IDENTIFIER, LATITUDE, LONGITUDE, 
                 OBSERVATION.DATE, TIME.OBSERVATIONS.STARTED,
                 OBSERVER.ID, DURATION.MINUTES, EFFORT.DISTANCE.KM, 
                 PROTOCOL.CODE, EFFORT.AREA.HA, ALL.SPECIES.REPORTED,
                 NUMBER.OBSERVERS)
    })
  ))

if (verbose) cat("...Handling checklist data...\n")

state_data_df <- bind_rows(state_data_list)
rm(state_data_list)

# Get the list of checklists
state_checklists <- state_data_df %>% 
  distinct(SAMPLING.EVENT.IDENTIFIER, LATITUDE, LONGITUDE, 
           OBSERVATION.DATE, TIME.OBSERVATIONS.STARTED,
           OBSERVER.ID, DURATION.MINUTES, EFFORT.DISTANCE.KM, 
           PROTOCOL.CODE, EFFORT.AREA.HA, ALL.SPECIES.REPORTED,
           NUMBER.OBSERVERS)

# Check non-distinct cases
duplicates <- state_checklists %>% 
  count(SAMPLING.EVENT.IDENTIFIER) %>% 
  filter(n > 1)
# It looks like duplicates occur when information about effort area isn't attributed to all obs.
# This seems wrong, so the solution is to delete the row with NA effort area
state_checklists <- state_checklists %>% 
  filter(!(SAMPLING.EVENT.IDENTIFIER %in% 
             duplicates$SAMPLING.EVENT.IDENTIFIER) |
           !is.na(EFFORT.AREA.HA))
# Check non-distinct cases AGAIN
duplicates <- state_checklists %>% 
  count(SAMPLING.EVENT.IDENTIFIER) %>% 
  filter(n > 1)

if (nrow(duplicates) > 0) stop(paste0("Duplicates still found. State code ", state_code))

# Now, we have all the checklist info. We simply have to write out that and also write out
# the observation columns.
write_csv(state_checklists, 
          paste0("data/", state_code, "_Aug2020_checklist_info.csv"))

if (verbose) cat("...Finished checklists. Reading files for counts...\n")

suppressWarnings(
  suppressMessages(
    state_data_list <- lapply(this_state_files, function(x) {
      read_csv(x, progress = F,
               col_types = list(COMMON.NAME = col_character(),
                                SCIENTIFIC.NAME = col_character(),
                                OBSERVATION.COUNT = col_character())) %>% 
        select(SAMPLING.EVENT.IDENTIFIER, COMMON.NAME,
               SCIENTIFIC.NAME, OBSERVATION.COUNT)
    })
  ))
state_data_df <- bind_rows(state_data_list)
rm(state_data_list)

if (verbose) cat("...Processing counts...\n")

state_obs_info <- state_data_df %>% 
  select(SAMPLING.EVENT.IDENTIFIER, SCIENTIFIC.NAME, 
         COMMON.NAME, OBSERVATION.COUNT) %>% 
  mutate(name_clean = gsub("['.()]", "", 
                           gsub("[ /]", "_", COMMON.NAME)))

species_to_process <- unique(state_obs_info$name_clean)

data_by_species <- list()

if (verbose) cat("...Handling duplicates...\n")

for (i in 1:length(species_to_process)) { # Handle subspecies duplicates
  # writeLines(as.character(i))
  data_by_species[[i]] <- 
    state_obs_info %>% 
    filter(name_clean == species_to_process[i]) %>% 
    select(-COMMON.NAME) %>% 
    group_by(SAMPLING.EVENT.IDENTIFIER, name_clean, SCIENTIFIC.NAME) %>% 
    summarise(total_count = get_sum(OBSERVATION.COUNT), .groups = "drop")
}

all_data_unique <- bind_rows(data_by_species)

write.csv(all_data_unique, 
          paste0("data/", state_code, "_Aug2020_species_counts.csv"))

if (verbose) cat("...Done.\n")


