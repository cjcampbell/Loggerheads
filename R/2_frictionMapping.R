
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(leastcostpath)
library(pbapply)

myProj <- "+proj=aea +lon_0=-77 +lat_1=17 +lat_2=30 +lat_0=23.5 +datum=WGS84 +units=m +no_defs"

# Load continuous surfaces ------------------------------------------------


# List individual IDs.
ids <- list.dirs("data/GoM Loggerheads", full.names = T) %>%  # List directories that correspond to each individual.
  grep("1",. , value = T) %>% # Filter only to subdirectories.
  word(3,3,sep="/") %>% # Pull out the ID part of the file path
  word(1,1,sep=" ") # Pull out only the first (numeric) part of the ID.

# Make a list where each item in the list corresponds to an individual, and the list contains all the rasts for that indiv.
r_cont <- lapply(ids, function(i) {
  assigns <- list.dirs("data/GoM Loggerheads", full.names = T) %>% 
    grep(pattern = i, ., value = T) %>% 
    list.files(., recursive = T, full.names = T, pattern = ".asc$") %>% 
    rast() %>% 
    terra::project(myProj)
  return(assigns)
})
names(r_cont) <- ids

sampleLocations <- read.csv("data/Assignment Lat Longs.csv") %>% 
  st_as_sf(coords = c("Stranding.Longitude", "Stranding.Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs") %>% 
  st_transform(myProj) %>% 
  dplyr::mutate(lyr = factor(Turtle.ID)) %>% 
  dplyr::mutate(sample_season = case_when(Stranding.Month %in% c(12,1:3) ~ "winter", TRUE ~ "summer"))


# Get bathymetry data -----------------------------------------------------

library(marmap)
b <- marmap::getNOAA.bathy(-98, -65, 13, 33)
bb <- rast(marmap::as.raster(b)) 
bb_sea <- bb%>% 
  project(myProj) %>% 
  classify(rcl = matrix(c(0, Inf, NA), ncol = 3, byrow = TRUE), include.lowest = TRUE)
bb_cs <- bb_sea %>% 
  classify(rcl = matrix(c(0, Inf, 0, -Inf, 0, 1), ncol = 3, byrow = TRUE), include.lowest = TRUE) %>% 
  create_cs(, neighbours = 8)


## What about a path connecting multiple subsampled points?

simulatePaths <- function(my_id = "13327") {
  rr <- r_cont[[my_id]]
  startPt <- sampleLocations %>% 
    dplyr::filter(Turtle.ID == my_id) %>% 
    dplyr::select(geometry)
  
  # Determine layer order.
  lyrOrder <- as.character(1:nlyr(rr)*50) 
  mypts <- lapply( lyrOrder, function(i) {
    # i <- 1
    pts_1 <- spatSample(rr[[i]], size = 1, method = "weights", replace = TRUE, xy = T) %>% 
      dplyr::select(-3) %>% 
      st_as_sf(coords = c("x", "y"), crs = myProj)
  }) %>% do.call("rbind", .)
  # TODO make sure points are in the right order...
  mypts2 <- rbind(startPt, mypts)
  
  mypaths <- lapply(1:(nrow(mypts2)-1), function(q){
    ls <- leastcostpath::create_lcp(bb_cs, mypts2[q,], mypts2[q+1,])
  }) 
  # Check for empty linestrings and remove.
  if( any(unlist(lapply(mypaths, st_is_empty))) ) return(NULL)
  
  mypaths <- mypaths %>% 
    do.call("rbind", .) %>% 
    st_union() %>% 
    st_as_sf()
  
  return(mypaths)
  
  ggplot() +
    geom_spatraster(bb_sea, mapping = aes()) +
    geom_sf(mypaths, mapping = aes()) +
    geom_sf(startPt, mapping = aes(), color = "green") +
    geom_sf(mypts, mapping = aes(), color = "white")
}

nreps <- 5000
savepath <- paste0("bin/reps", nreps)
if(!dir.exists(savepath)) dir.create(savepath)

start_time <- Sys.time()
for(myid in as.character(sampleLocations$Turtle.ID)) {
  if(file.exists(file.path(savepath, paste0("/pathSim_", myid, ".rds")))) next
  suppressWarnings({
    suppressMessages({
      set.seed(42)
      out <- pbreplicate(nreps, simulatePaths(my_id = myid), simplify = F) %>% 
        do.call("rbind", .) 
    })
  })
  out$length <- st_length(out)
  out$id <- myid
  saveRDS(out, file = file.path(savepath, paste0("/pathSim_", myid, ".rds")))
}
end_time <- Sys.time()
end_time - start_time

out_mls <- list.files(savepath, full.names = T) %>% 
  lapply(readRDS) %>% 
  do.call("rbind", .)

# Plot example + least-cost paths.
ggplot() +
  geom_spatraster(bb_sea, mapping = aes()) +
  geom_sf(
    slice(out_mls, 1:5, .by = id), 
    mapping = aes()) +
  geom_sf(
    slice(out_mls, which.min(length), .by = id), 
    mapping = aes(), color = "yellow") +
  facet_wrap(~id)

# Plot densities of paths.
# path_densities <- lapply(as.character(sampleLocations$Turtle.ID), function(myid) {
#   if(!myid %in% unique(out_mls$id)) return(NULL)
#   r1 <- rasterize(out_mls[out_mls$id == myid,], bb_sea, field = 1, fun = "sum", na.rm = T)
#   r2 <- r1 / unlist(terra::global(r1, "sum", na.rm = T))
#   names(r2) <- myid
#   return(r2)
# }) %>% rast
# 
# ggplot() +
#   geom_spatraster(path_densities, mapping = aes()) +
#   scale_fill_viridis_c(option = "turbo", trans = "log", na.value = NA) +
#   facet_grid(~lyr)

# Plot length-weighted densities of paths.
path_weighted_densities <- lapply(as.character(sampleLocations$Turtle.ID), function(myid) {
  if(!myid %in% unique(out_mls$id)) return(NULL)
  sampled_mls <- out_mls %>% 
    dplyr::filter(id == myid) %>% 
    dplyr::mutate(
      inv_length = 1/as.numeric(length) # ,
      # scaled_length = max(length)-length
      ) %>% 
    slice_sample(n=nreps/10, weight_by = inv_length)
  r1 <- rasterize(sampled_mls, bb_sea, field = 1, fun = "sum", na.rm = T)
  r2 <- r1 / unlist(terra::global(r1, "sum", na.rm = T))
  names(r2) <- myid
  return(r2)
}) %>% rast

ggplot() +
  geom_spatraster(path_weighted_densities, mapping = aes()) +
  scale_fill_viridis_c(option = "turbo", trans = "log10", na.value = NA) +
  facet_wrap(~lyr)

# Make individual level plots of length-weighted densities.

lapply(as.character(sampleLocations$Turtle.ID), function(myid) {
  
  ggplot() +
    geom_spatraster(bb_sea, mapping = aes()) +
    scale_fill_continuous(na.value = NA) +
    ggnewscale::new_scale_fill() +
    geom_spatraster(path_weighted_densities[[myid]], mapping = aes()) +
    scale_fill_viridis_c(option = "turbo", trans = "log10", na.value = NA) +
    coord_sf(xlim = c(-15e5, 5e5), ylim = c( -8e5, 11e5))
    
  
})




# Plot number of simulations included vs. minimum track length.
lapply(1:5000, function(x) {
  out_mls %>% 
    group_by(id) %>% 
    slice(1:x) %>% 
    slice(out_mls, which.min(length), .by = id)
  
  
})


