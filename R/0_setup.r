library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(leastcostpath)
library(pbapply)
library(marmap)
library(data.table)

theme_set(theme_minimal())

# Specify an equal-area projection before doing any area analyses
myProj <- "+proj=aea +lon_0=-77 +lat_1=17 +lat_2=30 +lat_0=23.5 +datum=WGS84 +units=m +no_defs"
