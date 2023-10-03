library(targets)
library(tarchetypes)


tar_option_set(
	packages = c("ncdf4", "tidync", "tidyverse", "magrittr", "lubridate", "gstat", "sp", "phylolm", "geiger", "phylolm", "MuMIn", "sisters")
)

source("_functions.R")

list(
  #files for spicule presence/absence analysis (some reused for skeletal complexity analysis)
  tar_target(reduced_data, create_reduced_data(tidydepth)),
  tar_target(tidydepth, tidync("C:\\Users\\jjera\\Documents\\RStuff\\GEBCO_2023_sub_ice_topo.nc")),
  tar_target(tidyph, tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\pH\\allpH.nc")),
  tar_target(tidysilica, tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\Silica\\allsi.nc")),
  tar_target(tidytemp, tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\Temperature\\alltemp.nc")),
  tar_target(phy, read.tree("c:\\Users\\jjera\\Documents\\RStuff\\Final Bootstrap Tree\\intree.dated.final.tre")),
  #environmental data extraction
  tar_target(sponges_data_for_idw, extracts_depths_around_sponges(reduced_data, tidydepth, bound=0.00416667*5)),
  tar_target(inferred_depths, infer_depths(sponges_data_for_idw, reduced_data)),
  tar_target(depths_and_ph, extract_ph(inferred_depths, tidyph)),
  tar_target(depths_pH_and_silica, extract_silica(depths_and_ph, tidysilica)),
  tar_target(depths_pH_silica_and_temp, extract_temp(depths_pH_and_silica, tidytemp)),
  tar_target(pruned_tree, datatree(depths_pH_silica_and_temp, phy)),
  #models for binary data
  tar_target(model_output_full, phyloglm(spicules ~ idw_depths + ph + silica + temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_depths_only, phyloglm(spicules ~ idw_depths, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_ph_only, phyloglm(spicules ~ ph, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_silica_only, phyloglm(spicules ~ silica, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_temperature_only, phyloglm(spicules ~ temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_only_intercept, phyloglm(spicules ~ 1, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_depth_ph, phyloglm(spicules ~ idw_depths + ph, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_depth_silica, phyloglm(spicules ~ idw_depths + silica, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_depth_temperature, phyloglm(spicules ~ idw_depths + temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_ph_silica, phyloglm(spicules ~ ph + silica, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_ph_temperature, phyloglm(spicules ~ ph + temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_silica_temperature, phyloglm(spicules ~ silica + temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_depth_ph_silica, phyloglm(spicules ~ idw_depths + ph + silica, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_depth_ph_temperature, phyloglm(spicules ~ idw_depths + ph + temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_depth_silica_temperature, phyloglm(spicules ~ idw_depths + silica + temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  tar_target(model_output_ph_silica_temperature, phyloglm(spicules ~ ph + silica + temperature, data=pruned_tree$data, phy=pruned_tree$phy, method="logistic_IG10")),
  #model comparison visualization
  tar_target(model_comparison, rbind(model_output_full$aic, model_output_depths_only$aic, model_output_ph_only$aic, model_output_silica_only$aic, model_output_temperature_only$aic, model_output_only_intercept$aic, model_output_depth_ph$aic, model_output_depth_silica$aic, model_output_depth_temperature$aic, model_output_ph_silica$aic, model_output_ph_temperature$aic, model_output_silica_temperature$aic, model_output_depth_ph_silica$aic, model_output_depth_ph_temperature$aic, model_output_depth_silica_temperature$aic, model_output_ph_silica_temperature$aic)),
  #visualizations for binary data analysis
  tar_target(spicule_presence_summaries, pruned_tree$data |> group_by(spicules) |> summarize(
    n=n(),
    median_depths=median(idw_depths),
    lower_depths=quantile(idw_depths, 0.025),
    upper_depths=quantile(idw_depths, 0.975),
    median_ph = median(ph),
    median_silica=median(silica),
    median_temperature=median(temperature)
  )),

  #files for spicule complexity analysis
  tar_target(complexityGBIF, read.delim("C:\\Users\\jjera\\Documents\\RStuff\\0163120-230224095556074.csv")),
  tar_target(complexity, read.csv("C:\\Users\\jjera\\Documents\\RStuff\\Complexity.csv")),
  tar_target(complexity_data, create_complexity_data(complexityGBIF, complexity)),
  #environmental data extraction
  tar_target(complexity_and_depths, complex_depth_extraction(complexity_data, tidydepth)),
  tar_target(complexity_depths_and_ph, complex_ph_extraction(complexity_and_depths, tidyph)),
  tar_target(complexity_depths_ph_and_silica, complex_silica_extraction(complexity_depths_and_ph, tidysilica)),
  tar_target(complexity_depths_ph_silica_and_temp, complex_temperature_extraction(complexity_depths_ph_and_silica, tidytemp)),
  tar_target(complex_pruned_tree, complex_datatree(complexity_depths_ph_silica_and_temp, phy)),
  #models for complexity data
  tar_target(complex_bm_model_output, complex_bm_model(complex_pruned_tree)),
  tar_target(complex_rr_model_output, complex_randomroot_model(complex_pruned_tree)),
  #try models with adjusted predictors
  tar_target(complex_adjusted_data, complex_data_adjustments(complexity_depths_ph_silica_and_temp)),
  tar_target(complex_adjusted_tree, complex_tree_adjustments(complex_adjusted_data, phy)),
  tar_target(adjusted_bm_model, complex_bm_model_adjustment(complex_adjusted_tree)),
  tar_target(adjusted_randomroot_model, complex_rr_model_adjustment(complex_adjusted_tree)),
  #try presence/absence models with adjusted predictors
  tar_target(adjusted_binary_data, binary_data_adjustments(depths_pH_silica_and_temp)),
  #try sister group analysis for presence/absence data
  tar_target(sister_group_model, sisters_analysis(depths_pH_silica_and_temp, phy))

)
#log mean depth, log max depth, log min depth, photic yes/no, nonphotic yes/no
#add interaction term for Calcarea (convert to calcite spicule yes/no)
#consider modeling distance from optimal pH and temp
#silica is fine untransformed
#look at regional sponge predation