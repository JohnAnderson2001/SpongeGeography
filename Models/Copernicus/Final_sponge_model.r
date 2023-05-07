library(ncdf4)
library(tidync)
library(tidyverse)
library(magrittr)

#cleaning GBIF data
dat <- read.delim("C:\\Users\\jjera\\Documents\\RStuff\\0163120-230224095556074.csv")
reldat <- dat[,c("gbifID", "class", "order", "family", "genus", "species", "decimalLatitude", "decimalLongitude", "depth", "eventDate", "day", "month", "year")]

complexdat <- read.csv("C:\\Users\\jjera\\Documents\\RStuff\\Complexity.csv")
colnames(complexdat)[1] ="species" #change first column so names match before merging

prefinaltestdat <- merge(reldat, complexdat, by="species", all=T) #merge with occurrence data
finaltestdat <- prefinaltestdat [complete.cases(prefinaltestdat), ] #68 species at this stage

#clean the merged dataset; remove rows with no depth and no date
finaltestdat2 <- finaltestdat[!is.na(as.numeric(finaltestdat$depth)), ] #removes all rows with depth value that is not a number
finaltestdat2$depth <- as.numeric(finaltestdat2$depth) #62 species
finaltestdat3 <- finaltestdat2[!is.na(as.numeric(finaltestdat2$decimalLatitude)), ] #60 species
sponges <- finaltestdat3[!finaltestdat3$year < 1950, ] #52 species after date cleaning
#write.csv(sponges, file="Occurrence_complexity_data_with_dates.csv")


sponges <- read.csv("C:\\Users\\jjera\\Documents\\RStuff\\Sponges_after_1950.csv")

#extracting pH
tidyph <- tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\pH\\allpH.nc")
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

#extracting silica
tidysilica <- tidync("C:\\Users\\jjera\\Documents\\RStuff\\Copernicus\\Silica\\allsi.nc")
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
#write.csv(sponges, file="anthropocene_sponges_with_data.csv")


#constructing models
sponges <- read.csv("C:\\Users\\jjera\\Documents\\RStuff\\anthropocene_sponges_with_data.csv")
cleandat <- sponges[,c("genus", "SpiculeTypes", "ph", "temperature", "silica", "depth"), ]
phdat <- aggregate(ph ~ genus, FUN="mean", data=cleandat)
tempdat <- aggregate(temperature ~ genus, FUN="mean", data=cleandat)
sildat <- aggregate(silica ~ genus, FUN="mean", data=cleandat)
depthdat <- aggregate(depth ~ genus, FUN="mean", data=cleandat)
spicdat <- aggregate(SpiculeTypes ~ genus, FUN="mean", data=cleandat)
finalsponges <- merge(phdat, tempdat, by="genus")
finalsponges <- merge(finalsponges, sildat, by="genus")
finalsponges <- merge(finalsponges, depthdat, by="genus")
finalsponges <- merge(finalsponges, spicdat, by="genus")
row.names(finalsponges) <- finalsponges$genus

library(phylolm)
phy <- read.tree("C:\\Users\\jjera\\Documents\\RStuff\\intree.dated.final.tre")
phy$tip.label <- sub("_[^_]+$", "", phy$tip.label)

RemoveDuplicateNames <- function(phy) {
    if(any(duplicated(phy$tip.label))) {
		phy$tip.label <- make.unique(phy$tip.label, sep = ".")
		cat("Duplicate names were found and have been renamed.\n")
	} else {
		cat("No duplicate names were found.\n")
	}
	return(phy)
}

library(geiger)
prunedtree <- treedata(phy=RemoveDuplicateNames(phy), data=finalsponges) #42 tips

prunedtree$data <- as.data.frame(prunedtree$data)
for (i in sequence(ncol(prunedtree$data))) {
    prunedtree$data[,i] <- as.numeric(prunedtree$data[,i])
}

phylomodelbmgenus <- phylolm(SpiculeTypes ~ ph + temperature + silica + depth, data=prunedtree$data, phy=prunedtree$phy, model="BM")
summary(phylomodelbmgenus)
phylomodelrandomrootgenus <- phylolm(SpiculeTypes ~ ph + temperature + silica + depth, data=prunedtree$data, phy=prunedtree$phy, model="OUrandomRoot")
summary(phylomodelrandomrootgenus)

library(MuMIn)
bestmodelbmgenus <- dredge(phylomodelbmgenus)
bestmodelrandomrootgenus <- dredge(phylomodelrandomrootgenus)

averagemodel <- model.avg(bestmodelrandomrootgenus, full=TRUE)
MuMIn::sw(bestmodelrandomrootgenus)