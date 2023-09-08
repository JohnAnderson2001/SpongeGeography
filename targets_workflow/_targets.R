library(targets)
library(tarchetypes)


tar_option_set(
	packages = c("ncdf4", "tidync", "tidyverse", "magrittr", "lubridate", "gstat", "sp")
)

source("_functions.R")

list(
  tar_target(reduced_data, create_reduced_data(tidydepth)),
  tar_target(tidydepth, tidync("C:\\Users\\jjera\\Documents\\RStuff\\GEBCO_2023_sub_ice_topo.nc")),
  tar_target(tidyph, tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\pH\\allpH.nc")),
  tar_target(sponges_data_for_idw, extracts_depths_around_sponges(reduced_data, bound=0.00416667*5)),
  tar_target(inferred_depths, infer_depths(sponges_data_for_idw)),
  tar_target(depths_and_ph, extract_ph(inferred_depths, tidyph))
)

