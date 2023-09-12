library(targets)
library(tarchetypes)


tar_option_set(
	packages = c("ncdf4", "tidync", "tidyverse", "magrittr", "lubridate", "gstat", "sp", "phylolm", "geiger", "phylolm")
)

source("_functions.R")

list(
  tar_target(reduced_data, create_reduced_data(tidydepth)),
  tar_target(tidydepth, tidync("C:\\Users\\jjera\\Documents\\RStuff\\GEBCO_2023_sub_ice_topo.nc")),
  tar_target(tidyph, tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\pH\\allpH.nc")),
  tar_target(tidysilica, tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\Silica\\allsi.nc")),
  tar_target(tidytemp, tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\Temperature\\alltemp.nc")),
  tar_target(phy, read.tree("c:\\Users\\jjera\\Documents\\RStuff\\Final Bootstrap Tree\\intree.dated.final.tre")),
  tar_target(sponges_data_for_idw, extracts_depths_around_sponges(reduced_data, tidydepth, bound=0.00416667*5)),
  tar_target(inferred_depths, infer_depths(sponges_data_for_idw, reduced_data)),
  tar_target(depths_and_ph, extract_ph(inferred_depths, tidyph)),
  tar_target(depths_pH_and_silica, extract_silica(depths_and_ph, tidysilica)),
  tar_target(depths_pH_silica_and_temp, extract_temp(depths_pH_and_silica, tidytemp)),
  tar_target(pruned_tree, datatree(depths_pH_silica_and_temp, phy)),
  tar_target(model_output, phyloglm(spicules ~ idw_depths + ph + silica + temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10"))

)
