# get all species from DFO

# need to read data from all survey areas. Just read them and bind them together
library(tidyverse)
library(data.table)

setwd("C:/Users/arove/Documents/GOA/SDM/Canada")

all_files <- paste(getwd(),"data/surveys",sep = "/")

all_areas <- dir(all_files)

get_species <- function(survey_area){
  area <- survey_area
  dfo_path <- paste(all_files,area,sep="/")
  catch <- read.csv(paste(dfo_path,paste(area,"catch.csv",sep="_"),sep="/"))
  return(catch)
}

all_species <- rbindlist(lapply(all_areas,FUN = get_species))

# save all species as Rdata - it seems to be relatively small and it will be easier to pull species from that master instead of reading data in every time and subsetting

save(all_species, file = "data/DFO_species_data/all_species_catch.Rdata")

all_species <- all_species %>% dplyr::select(Species.code:English.common.name) %>% distinct()

# try and fill this with the mapping we have already done, by species name
race_species <- read.csv("data/RACE_species_goa_Atlantis.csv")

race_species <- race_species %>% dplyr::select(Atlantis.group,Common.Name,Scientific.Name)

# base matching on scientific name, as common names may vary. This will not get all of them (some of the names were missing in RACE too, some will be spelled different, etc)

match_name <- function(Scientific.name) {
  At_gr <- race_species[grep(Scientific.name,race_species$Scientific.Name,ignore.case = TRUE),]$Atlantis.group
  At_gr <- At_gr[!is.na(At_gr)]
  if(length(unique(factor(At_gr)))==1) {
    At_gr <- At_gr[1]
  } else {
    At_gr <- paste(unique(factor(At_gr)),sep = " ", collapse = " ")
  }
  return(At_gr)
}

all_species <- all_species %>% rowwise() %>% mutate(Atlantis.group=match_name(Scientific.name))

# need to finish the rest manually (also to check that everything is sensible)

write.csv(all_species,"data/dfo_species.csv",row.names=FALSE)
