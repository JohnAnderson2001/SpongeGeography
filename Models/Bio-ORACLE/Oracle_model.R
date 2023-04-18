library(phylolm)

dat <- read.csv("Occurrence_oracle_averages.csv")

row.names(dat) <- dat$species
phylomodelbm <- phylolm(SpiculeTypes ~ pH + temp + silicate + depth, data=dat, phy=tree, model="BM")
summary(phylomodel)
phylomodelrandomroot <- phylolm(SpiculeTypes ~ pH + temp + silicate + depth, data=dat, phy=tree, model="OUrandomRoot")
summary(phylomodelrandomroot)

library(MuMIn)
all_results_BM <- MuMIn::dredge(phylomodelbm)
summary(all_results_BM)
print(all_results_BM)

all_results_OU <- MuMIn::dredge(phylomodelrandomroot)
summary(all_results_OU)
print(all_results_BM)

BM_data <- as.data.frame(print(all_results_BM))
BM_data$Model <- "BM"
OU_data <- as.data.frame(print(all_results_OU))
OU_data$Model <- "OUrandomRoot"
completedata <- rbind(BM_data, OU_data)
completedata$delta <- completedata$AICc - min(completedata$AICc)
write.csv(completedata, "Oracle_results_summary.csv")

#try old data
framedat <- read.csv("CompleteSpongeData.csv")

phdat <- aggregate(pH ~ species, FUN="mean", data=framedat)
tempdat <- aggregate(temp ~ species, FUN="mean", data=framedat)
sildat <- aggregate(silicate ~ species, FUN="mean", data=framedat)
depthdat <- aggregate(depth ~ species, FUN="mean", data=framedat)
spicdat <- aggregate(Type ~ species, FUN="mean", data=framedat)

finalolddat <- merge(finalolddat, depthdat, by="species")
#write.csv(finalolddat, "Occurrence_oracle_averages_binaryspicdata.csv")
olddat <- read.csv("Occurrence_oracle_averages_binaryspicdata.csv")
olddat$species <- sub(" ", "_", olddat$species)
row.names(olddat) <- olddat$species

phylomodel_binary <- phyloglm(Type ~ pH + temp + silicate + depth, data=olddat, phy=tree, method="logistic_IG10")
summary(phylomodel_binary)

#all_results_binary <- MuMIn::nobs(phylomodel_binary)


