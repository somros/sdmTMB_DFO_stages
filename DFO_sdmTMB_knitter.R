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

# read in BGM for projection
atlantis_bgm <- read_bgm("data/GOA_WGS84_V4_final.bgm")

# prediciton grid, make the depth positive for consistency with RACE data. Also append coordinates, and add time dimension
load("data/atlantis_grid_depth.Rdata")
atlantis_coords <- atlantis_grid_depth %>% st_as_sf(coords = c("x", "y"), crs = atlantis_bgm$extra$projection) %>%
  st_transform(crs = "+proj=longlat +datum=WGS84") %>% dplyr::select(geometry)

atlantis_grid_template <- cbind(atlantis_grid_depth, do.call(rbind, st_geometry(atlantis_coords)) %>%
                                  as_tibble() %>% setNames(c("lon","lat")))

# #Workflow:
#   
# 1. Concatenate catch data for all areas - starting from raw
# 2. Join with atlantis key
# 3. Concatenate tow info for all areas - starting from raw
# 4. Select species (then we will run it in loop)
# 5. Join catch and tow to add zero catches
# 6. Calculate CPUE from catch and effort.
# 7. clean and rename columns to be piped in the sdmTMB markdown

#This chunk is to run only once. Points 1-3.
# concatenate catch data for the three survey areas

path_to_surveys <- paste(getwd(),"data/surveys",sep = "/")
all_surveys <- dir(path_to_surveys)

# drop WCHG
all_surveys <- setdiff(all_surveys, "WCHG")

get_catch <- function(survey_area){
  sa <- survey_area
  sa_path <- paste(path_to_surveys,sa,sep="/")
  catch <- read.csv(paste(sa_path,paste(sa,"catch.csv",sep="_"),sep="/"))
  catch$Survey.area <- sa
  return(catch)
}
catch_all <- rbindlist(lapply(all_surveys,FUN = get_catch))

# join Atlantis group lookup information
atlantis_key <- read.csv("data/dfo_species_atlantis.csv")
atlantis_key <- atlantis_key %>% filter(Atlantis.group != "?")
atlantis_groups <- read.csv("data/GOA_Groups.csv", fileEncoding = "UTF-8-BOM")
atlantis_groups <- atlantis_groups %>% select(Code,Name,LongName)
#join into one key with long names
atlantis_key <- atlantis_key %>% left_join(atlantis_groups, by = c("Atlantis.group"="Code"))

#join catch data with Atlantis key by scientific name, which seems to be unique to each record in the DFO data
catch_all <- catch_all %>% left_join(atlantis_key %>% select(Scientific.name,Atlantis.group,Name), by = "Scientific.name")

# drop the NA in the groups - it means it is one of the undecided
catch_all <- catch_all %>% filter(!is.na(Atlantis.group))

# concatenate effort information from the three survey areas
get_effort <- function(survey_area){
  sa <- survey_area
  sa_path <- paste(path_to_surveys,sa,sep="/")
  effort <- read.csv(paste(sa_path,paste(sa,"effort.csv",sep="_"),sep="/"))
  return(effort)
}
effort_all <- rbindlist(lapply(all_surveys,FUN = get_effort)) # 4533 sets as of May 2021
# how many sets are deeper than 1000 m?
# nrow(effort_all %>% filter(Bottom.depth..m. > 1000))/nrow(effort_all)*100 # less than 1%


#This chunk will be looped for each Atlantis group. 

ag <- catch_all %>% select(Atlantis.group) %>% distinct() %>% pull()

dfo_knitter <- function(this_group) {
  this_group <- this_group
  
  species_catch <- catch_all %>% filter(Atlantis.group==this_group)
  species_catch <- species_catch %>% 
    group_by(Survey.area,Survey.Year,Trip.identifier,Set.number,Atlantis.group,Name) %>% 
    summarise(Catch.weight.kg = sum(Catch.weight..kg.), Catch.count = sum(Catch.count..pieces.)) %>%
    ungroup()
  
  species_all_hauls <- species_catch %>% full_join(effort_all) # Joining, by = c("Survey.Year", "Trip.identifier", "Set.number")
  
  # add one column with a combination of trip and set, and one for the name of the area
  species_all_hauls <- species_all_hauls %>% mutate(hauljoin = paste(Trip.identifier,Set.number,sep="_"))
  
  # calculate area swept in km2, and CPUE for weight and numbers as biom/aream/1000
  species_all_hauls <- species_all_hauls %>% mutate(area_sweptkm2 = Distance.towed..m.*Trawl.door.spread..m./1e+6,
                                                    cpue_kgkm2 = Catch.weight.kg/area_sweptkm2,
                                                    cpue_numkm2 = Catch.count/area_sweptkm2)
  
  # extract month here - RACE data is 6,7,8, but these seem to start in May. Come back to this and try with overlapping months only
  # species_all_hauls <- species_all_hauls %>% mutate(month = month(ymd(Set.date))) %>% filter(month %in% 6:8)
  
  #assume that NA catches are 0, since they come from sets with no catch
  species_all_hauls$cpue_kgkm2[which(is.na(species_all_hauls$cpue_kgkm2))] <- 0
  
  #filter out NA depth
  species_all_hauls <- species_all_hauls %>% filter(!is.na(Bottom.depth..m.))
  
  # drop some columns and rename as for the RACE bottom trawl surveys
  species_all_hauls <- species_all_hauls %>% select(Survey.Year, hauljoin, End.latitude, End.longitude, Bottom.depth..m.,Atlantis.group,Name,cpue_kgkm2,cpue_numkm2,Survey.area) %>%
    set_names(c(
      "year",
      "hauljoin",
      "lat",
      "lon",
      "depth",
      "species_code",
      "name",
      "biom_kgkm2",
      "num_km2",
      "survey"))
  
  # here it needs to call the sdmTMB document
  rmarkdown::render(
    'DFO_sdmTMB_template.Rmd', output_file = paste0("outputs/knitted_output_afterSean/", this_group, "_", cutoff, '.html')
  )
}

purrr::map(ag, possibly(dfo_knitter, NA))

#Have a look at this on a map.
# dfo_sf <- dfo_data %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)
# 
# # coast mask
# coast <- map("worldHires", regions = c("Canada", "USA"), plot = FALSE, fill = TRUE)
# coast <- coast %>% st_as_sf() #%>% st_transform(crs = atlantis_bgm$extra$projection)
# 
# ggplot()+
#   geom_sf(data = dfo_sf, aes(colour = log1p(biom_kgkm2)))+
#   scale_colour_viridis_c()+
#   geom_sf(data = coast)+
#   coord_sf(xlim = c(min(dfo_data$lon),max(dfo_data$lon)), ylim=c(min(dfo_data$lat),max(dfo_data$lat)))+
#   theme_minimal()+
#   facet_wrap(~year, ncol = 3)+
#   labs(title = paste(species,"CPUE from DFO bottom trawl survey", sep = " "))

#Plot data from all years together.
# ggplot()+
#   geom_sf(data = dfo_sf, aes(colour = log1p(biom_kgkm2)))+
#   scale_colour_viridis_c()+
#   geom_sf(data = coast)+
#   coord_sf(xlim = c(min(dfo_data$lon),max(dfo_data$lon)), ylim=c(min(dfo_data$lat),max(dfo_data$lat)))+
#   theme_minimal()+
#   labs(title = paste(species,"CPUE from DFO bottom trawl survey", sep = " "))

#Write out data.
# assign(paste("ATF","CPUE_DFO",sep="_"), dfo_data)
# save(ATF_CPUE_DFO,file=paste0("data/DFO_species_data/","ATF",".Rdata"))
