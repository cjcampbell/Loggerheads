
# Setup, load files -------------------------------------------------------

source("R/0_setup.r")

# Load list of tidy binary rasters.
listOfRasts <- list.files("data/surfaces/binary", full.names = T) %>% 
  lapply(rast)

names(listOfRasts) <- tools::file_path_sans_ext(list.files("out/binary assignment maps"))


# Sum origins for each indiv ----------------------------------------------

summedSurfaces <- lapply(listOfRasts, function(y) {sum(y)/nlyr(y)}) %>% 
  rast()

# Visualize.
p_overlaps <- ggplot() +
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
ggsave(p_overlaps, filename = file.path("figs/p_overlaps.png"), width = 10, height = 10, dpi = 600)

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
    #id = ids,
    area_km2 = count*cellRes
  ) %>% 
  cbind(id = names(summedSurfaces)) %>% 
  dplyr::select(id, area_km2)


