#find sponges in GBIF
library(targets)
tar_load(complexityGBIF)
relevanttest <- complexityGBIF[,c("gbifID", "class", "order", "family", "genus", "species", "decimalLatitude", "decimalLongitude", "depth", "year")]
test <- relevanttest[(complexityGBIF$species %in% "Spinularia spinularia"), ] #replace string with the name of a species to see if it is in GBIF
test
