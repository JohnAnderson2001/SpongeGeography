library(ncdf4)
library(tidync)
library(tidyverse)
library(magrittr)

#extracting pH
tidyph <- tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\pH\\allpH.nc")
sponges <- read.csv("C:\\Users\\jjera\\Documents\\RStuff\\Occurrence_complexity_data_with_dates.csv")
sponges$ph <- NA

for(sponge_index in sequence(nrow(sponges))) {
    #target_latitude <- tidyph$transforms$latitude$latitude[(which.min(abs(tidyph$transforms$latitude$latitude - sponges$decimalLatitude[sponge_index])))]
    #target_longitude <- tidyph$transforms$longitude$longitude[(which.min(abs(tidyph$transforms$longitude$longitude - sponges$decimalLongitude[sponge_index])))]
    target_depth <- tidyph$transforms$depth$depth[(which.min(abs(tidyph$transforms$depth$depth - sponges$depth[sponge_index])))]
    phs <- tidyph %>% hyper_filter(
        latitude = between(latitude, sponges$decimalLatitude[sponge_index]-1, sponges$decimalLatitude[sponge_index]+1),
        longitude = between(longitude, sponges$decimalLongitude[sponge_index]-1, sponges$decimalLongitude[sponge_index]+1),
        depth = between(depth, target_depth-10, target_depth+10)
         ) %>% hyper_array()
    sponges$ph[sponge_index] <- median(phs$ph, na.rm=TRUE)
    cat("\r", sponge_index, " of ", nrow(sponges), " it found ", length(phs$ph), " values")

}

#pH extraction successful
#next step saved as spongestep1.csv

#extracting silica
tidysilica <- tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\Silica\\allsi.nc")
sponges <- read.csv("C:\\Users\\jjera\\Documents\\RStuff\\spongestep1.csv")
sponges$silica <- NA

for(sponge_index in sequence(nrow(sponges))) {
    #target_latitude <- tidyph$transforms$latitude$latitude[(which.min(abs(tidyph$transforms$latitude$latitude - sponges$decimalLatitude[sponge_index])))]
    #target_longitude <- tidyph$transforms$longitude$longitude[(which.min(abs(tidyph$transforms$longitude$longitude - sponges$decimalLongitude[sponge_index])))]
    target_depth <- tidysilica$transforms$depth$depth[(which.min(abs(tidysilica$transforms$depth$depth - sponges$depth[sponge_index])))]
    sils <- tidysilica %>% hyper_filter(
        latitude = between(latitude, sponges$decimalLatitude[sponge_index]-1, sponges$decimalLatitude[sponge_index]+1),
        longitude = between(longitude, sponges$decimalLongitude[sponge_index]-1, sponges$decimalLongitude[sponge_index]+1),
        depth = between(depth, target_depth-10, target_depth+10)
         ) %>% hyper_array()
    sponges$silica[sponge_index] <- median(sils$si, na.rm=TRUE)
    cat("\r", sponge_index, " of ", nrow(sponges), " it found ", length(sils$si), " values")

}

#extracting temp
tidytemp <- tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\Temperature\\alltemp.nc")
sponges$temperature <- NA

for(sponge_index in sequence(nrow(sponges))) {
    #target_latitude <- tidyph$transforms$latitude$latitude[(which.min(abs(tidyph$transforms$latitude$latitude - sponges$decimalLatitude[sponge_index])))]
    #target_longitude <- tidyph$transforms$longitude$longitude[(which.min(abs(tidyph$transforms$longitude$longitude - sponges$decimalLongitude[sponge_index])))]
    #target_depth <- tidysilica$transforms$depth$depth[(which.min(abs(tidysilica$transforms$depth$depth - sponges$depth[sponge_index])))]
    temps <- tidytemp %>% hyper_filter(
        latitude = between(latitude, sponges$decimalLatitude[sponge_index]-1, sponges$decimalLatitude[sponge_index]+1),
        longitude = between(longitude, sponges$decimalLongitude[sponge_index]-1, sponges$decimalLongitude[sponge_index]+1),
        #depth = between(depth, target_depth-10, target_depth+10)
         ) %>% hyper_array()
    sponges$temperature[sponge_index] <- median(temps$bottomT, na.rm=TRUE)
    cat("\r", sponge_index, " of ", nrow(sponges), " it found ", length(temps$bottomT), " values")

}
#write.csv(sponges, "spongeswithdata.csv")

#model creation
#sponges <- read.csv("spongeswithdata.csv")
cleandat <- sponges[,c("species", "SpiculeTypes", "ph", "temperature", "silica", "depth"), ]
cleandat$species <- sub(" ", "_", cleandat$species)
phdat <- aggregate(ph ~ species, FUN="mean", data=cleandat)
tempdat <- aggregate(temperature ~ species, FUN="mean", data=cleandat)
sildat <- aggregate(silica ~ species, FUN="mean", data=cleandat)
depthdat <- aggregate(depth ~ species, FUN="mean", data=cleandat)
spicdat <- aggregate(SpiculeTypes ~ species, FUN="mean", data=cleandat)
finalsponges <- merge(finalsponges, spicdat, by="species") #42 species
#write.csv(sponges, "spongeswithdata_averages.csv")
#finalsponges <- read.csv("spongeswithdata_averages.csv")
row.names(finalsponges) <- finalsponges$species

library(phylolm)
tree <- read.tree("C:\\Users\\jjera\\Documents\\RStuff\\Partition.txt.treefile")
phylomodelbm <- phylolm(SpiculeTypes ~ ph + temperature + silica + depth, data=finalsponges, phy=tree, model="BM")
summary(phylomodelbm)
phylomodelrandomroot <- phylolm(SpiculeTypes ~ ph + temperature + silica + depth, data=finalsponges, phy=tree, model="OUrandomRoot")
summary(phylomodelrandomroot)