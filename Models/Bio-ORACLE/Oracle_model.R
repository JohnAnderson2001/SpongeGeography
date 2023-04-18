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
