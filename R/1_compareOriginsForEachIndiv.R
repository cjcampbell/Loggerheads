
# Setup, load files -------------------------------------------------------

library(tidyverse)
library(terra)
library(tidyterra)
library(sf)

# Specify an equal-area projection before doing any area analyses
myProj <- "+proj=aea +lon_0=-77 +lat_1=17 +lat_2=30 +lat_0=23.5 +datum=WGS84 +units=m +no_defs"

# List individual IDs.
ids <- list.dirs("data/GoM Loggerheads", full.names = T) %>%  # List directories that correspond to each individual.
  grep("1",. , value = T) %>% # Filter only to subdirectories.
  word(3,3,sep="/") %>% # Pull out the ID part of the file path
  word(1,1,sep=" ") # Pull out only the first (numeric) part of the ID.

# Make a list where each item in the list corresponds to an individual, and the list contains all the rasts for that indiv.
listOfRasts <- lapply(ids, function(i) {
  assigns <- list.dirs("data/GoM Loggerheads", full.names = T) %>% 
    grep(pattern = i, ., value = T) %>% 
    list.files(., recursive = T, full.names = T, pattern = ".tif$") %>% 
    rast() %>% 
    terra::project(myProj)
  return(assigns)
})
names(listOfRasts) <- ids



listOfRasts[["13327"]] %>% plot
listOfRasts[["12952"]] %>% plot


# Load sample sites -------------------------------------------------------

sampleLocations <- read.csv("data/Assignment Lat Longs.csv") %>% 
  st_as_sf(coords = c("Stranding.Longitude", "Stranding.Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs") %>% 
  st_transform(myProj) %>% 
  dplyr::mutate(lyr = factor(Turtle.ID)) %>% 
  dplyr::mutate(sample_season = case_when(Stranding.Month %in% c(12,1:3) ~ "winter", TRUE ~ "summer"))



# Load bathymetry map ------------------------------------------------------
library(marmap)

b <- marmap::getNOAA.bathy(-93, -70, 16, 32)
marmap::autoplot.bathy(b, geom=c("raster"))

bb <- rast(marmap::as.raster(b))


# Sum origins for each indiv ----------------------------------------------

summedSurfaces <- lapply(listOfRasts, function(y) {sum(y)/nlyr(y)}) %>% 
  rast()

# Visualize.
ggplot() +
  geom_spatraster(summedSurfaces, mapping = aes()) +
  geom_sf(data = sampleLocations, mapping = aes(shape = sample_season), color = "white") +
  scale_fill_viridis_c(na.value = NA) +
  # ggnewscale::new_scale_fill() +
  # geom_spatraster(bb, mapping = aes(), alpha = 0.3) +
  facet_wrap(~lyr) 
# Yellow (value = 1) means that all origins from each scute sample sample overlapped with a region,
# i.e., limited evidence of movement. If that area is nonexistant or very small, that's
# evidence of movement to me! Let's also identify the size of the areas where the 100% of the origins overlap
# For each individual.

# Make new raster showing where 100% of scute samples overlap
allOverlapping <- summedSurfaces == 1

# Calculate resolution of cells in km2
cellRes <- prod(res(allOverlapping)/1e3)

# Make frequency table and summarize.
allOverlapping %>% 
  freq() %>% 
  complete(layer, value, fill = list(count = 0)) %>% 
  dplyr::filter(value == 1) %>% 
  dplyr::mutate(
    id = ids,
    area_km2 = count*cellRes
  ) %>% 
  dplyr::select(id, area_km2)


