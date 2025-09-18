
source("R/0_setup.r")

# Load continuous surfaces ------------------------------------------------

ids <- tools::file_path_sans_ext(list.files("data/surfaces/continuous"))
r_cont <- list.files("data/surfaces/continuous", pattern = "tif", full.names = T) %>% 
  lapply(FUN = rast)
names(r_cont) <- ids

# Load sample locations (not known for all individuals).
sampleLocations <- read.csv("data/Assignment Lat Longs.csv") %>% 
  st_as_sf(coords = c("Stranding.Longitude", "Stranding.Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs") %>% 
  st_transform(myProj) %>% 
  dplyr::mutate(lyr = factor(Turtle.ID)) %>% 
  dplyr::mutate(sample_season = case_when(Stranding.Month %in% c(12,1:3) ~ "winter", TRUE ~ "summer"))


# Get bathymetry data -----------------------------------------------------

b <- marmap::getNOAA.bathy(-98, -65, 13, 33)
bb <- rast(marmap::as.raster(b)) 
bb_sea <- bb%>% 
  project(myProj) %>% 
  classify(rcl = matrix(c(0, Inf, NA), ncol = 3, byrow = TRUE), include.lowest = TRUE)
bb_cs <- bb_sea %>% 
  classify(rcl = matrix(c(0, Inf, 0, -Inf, 0, 1), ncol = 3, byrow = TRUE), include.lowest = TRUE) %>% 
  create_cs(, neighbours = 8)

## Generate a path connecting multiple subsampled points.

simulatePaths <- function(my_id = "13327") {
  rr <- r_cont[[my_id]]

  # Determine layer order.
  lyrOrder <- as.character(1:nlyr(rr)*50) 
  mypts <- lapply( lyrOrder, function(i) {
    # i <- 1
    pts_1 <- spatSample(rr[[i]], size = 1, method = "weights", replace = TRUE, xy = T) %>% 
      dplyr::select(-3) %>% 
      st_as_sf(coords = c("x", "y"), crs = myProj)
  }) %>% do.call("rbind", .)
  
  # TODO make sure points are in the right order...
  
  if(my_id %in% sampleLocations$Turtle.ID) {
    startPt <- sampleLocations %>% 
      dplyr::filter(Turtle.ID == my_id) %>% 
      dplyr::select(geometry)
    mypts2 <- rbind(startPt, mypts)
  } else {
    mypts2 <- mypts
  }
  
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
  
  if(my_id %in% sampleLocations$Turtle.ID) {
  ggplot() +
    geom_spatraster(bb_sea, mapping = aes()) +
    geom_sf(mypaths, mapping = aes()) +
    geom_sf(startPt, mapping = aes(), color = "green") +
    geom_sf(mypts, mapping = aes(), color = "white")
  } else {
    ggplot() +
      geom_spatraster(bb_sea, mapping = aes()) +
      geom_sf(mypaths, mapping = aes()) +
      geom_sf(mypts, mapping = aes(), color = "white")
  }
}

nreps <- 5000
savepath <- paste0("bin/reps", nreps)
if(!dir.exists(savepath)) dir.create(savepath)

start_time <- Sys.time()
for(myid in as.character(ids)) {
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

saveRDS(out_mls, file = paste0("out/paths_", nreps, "_allSimulations.rds"))

# Plot example + least-cost paths.
p_exPaths <- ggplot() +
  geom_spatraster(bb_sea, mapping = aes()) +
  scale_fill_viridis_c("Depth (m)", option = "mako", end = 0.65) +
  ggnewscale::new_scale_fill() +
  geom_sf(
    slice(out_mls, 1:5, .by = id), 
    mapping = aes()) +
  geom_sf(
    slice(out_mls, which.min(length), .by = id), 
    mapping = aes(), color = "yellow") +
  facet_wrap(~id) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
ggsave(p_exPaths, filename = file.path("figs/example_paths.png"), width = 10, height = 10, dpi = 300)

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
path_weighted_densities <- lapply(ids, function(myid) {
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

# Plot path densities.
bb_sea_dup <- c(rep(bb_sea,length(ids)))
names(bb_sea_dup) <- ids

p_pathDesities <- ggplot() +
  geom_spatraster(bb_sea_dup, mapping = aes()) +
  scale_fill_viridis_c("Depth (m)", option = "mako", end = 0.65) +
  ggnewscale::new_scale_fill() +
  geom_spatraster(path_weighted_densities, mapping = aes()) +
  scale_fill_viridis_c("Length-weighted\npath density", option = "turbo", trans = "log10", na.value = NA, begin = 0.5) +
  facet_wrap(~lyr) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_sf(xlim = c(-15e5, 5e5), ylim = c( -8e5, 11e5))
ggsave(p_pathDesities, filename = file.path("figs/pathDensities.png"),  width = 12, height = 12, dpi = 600)


# Make individual level plots of length-weighted densities.

lapply(ids, function(myid) {
  
  ggplot() +
    geom_spatraster(bb_sea, mapping = aes()) +
    scale_fill_viridis_c("Depth (m)", option = "mako", end = 0.65, na.value = "grey50") +
    ggnewscale::new_scale_fill() +
    geom_spatraster(path_weighted_densities[[myid]], mapping = aes(), alpha = 0.9) +
    scale_fill_viridis_c("Length-weighted\npath density", option = "turbo", trans = "log10", na.value = NA, begin = 0.5) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_sf(xlim = c(-15e5, 5e5), ylim = c( -8e5, 11e5))
    
  ggsave(filename = paste0("figs/pathDensity_", myid, ".png"),  width = 6, height = 6, dpi = 600)
  
})



# Plot number of simulations included vs. minimum track length.
sens <- lapply(c(1:10, seq(10,100,by=10), seq(100, 5000, by = 100)), function(y) {
  out_mls %>% 
    as.data.frame %>% 
    slice(1:y, .by = id) %>% 
    slice( which.min(length), .by = id) %>% 
    mutate(nsims = y)
}) %>% 
  data.table::rbindlist()

ggplot(sens) +
  geom_path(aes(x = nsims, y = as.numeric(length)/1e3, group = id, color = id)) +
  scale_color_viridis_d("Turtle ID", option = "turbo") +
  
  scale_y_continuous("Minimum path length (m)", expand = c(0,0)) +
  scale_x_log10("Number of simulations") 
  # scale_x_log10("Number of simulations", expand = expansion(add = c(0, 10))) +
  # geom_text(
  #   data = dplyr::filter(sens, nsims == 5000),
  #   aes(x = 5000, y =  as.numeric(length)/1e3, label = id),
  #   hjust = 0, size = 3
  # ) +
  # theme(legend.position = "none")
  
ggsave("figs/simulation sensitivity.png", dpi = 600, width = 6, height = 6)

