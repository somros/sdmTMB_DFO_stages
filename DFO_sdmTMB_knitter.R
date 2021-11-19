#This document puts the Synoptic Bottom Trawl Survey data for groundfish into the same format that we use for RACE data, with the aim of using the same sdmTMB workflow on these.
#The code then knits the markdown containing the sdmTMB workflow, plots and diagnostics
#It is based on the data as downloaded from [here](https://open.canada.ca/data/en/dataset/a278d1af-d567-4964-a109-ae1e84cbd24a). 
#Note that this will probably end up being based on biomass CPUE instead of numbers.

library(sdmTMB)
library(rbgm)
library(viridis)
library(kableExtra)
library(tidyverse)
library(sf)
library(data.table)
library(maps)
library(mapdata)
library(lubridate)

select <- dplyr::select

# Settings for sdmTMB.
cutoff <- 20

# coast for plotting
coast <- map("worldHires", regions = c("Canada", "USA"), plot = FALSE, fill = TRUE)
coast_sf <- coast %>% st_as_sf() 

# haul info
load('catch_to_CPUE_DFO/hauls_dfo.Rdata') # as accessible on AKFIN Answers with no further modifications

# atlantis bgm
atlantis_bgm <- read_bgm("data/GOA_WGS84_V4_final.bgm")

# prediciton grid, make the depth positive for consistency with RACE data. Also append coordinates, and add time dimension
load("data/atlantis_grid_depth.Rdata")
atlantis_coords <- atlantis_grid_depth %>% st_as_sf(coords = c("x", "y"), crs = atlantis_bgm$extra$projection) %>%
  st_transform(crs = "+proj=longlat +datum=WGS84") %>% dplyr::select(geometry)

atlantis_grid_template <- cbind(atlantis_grid_depth, do.call(rbind, st_geometry(atlantis_coords)) %>%
                                  as_tibble() %>% setNames(c("lon","lat")))
#######################################################################################
# load DFO data

load("catch_to_cpue_DFO/cpue_by_stage_DFO.Rdata")

# loop over groups
all_groups <- catch_all_hauls %>% select(species_code, stage) %>% distinct()

cpue_knitter <- function(this_group,this_stage){
  dfo_data <- catch_all_hauls %>% filter(species_code == this_group & stage==this_stage) %>% select(-survey)
  
  rmarkdown::render(
    'DFO_sdmTMB_template.Rmd', 
    output_file = paste0("outputs/", this_group, this_stage, "_", cutoff, '.html')
  )
}

# run for all groups, start for next group if the model will not converge
purrr::map2(all_groups$species_code, all_groups$stage, possibly(cpue_knitter,NA))
