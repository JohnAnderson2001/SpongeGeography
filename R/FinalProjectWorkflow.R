#open files
finalspics <- read.csv("finalspics.csv")
finalnospics <-read.csv("finalnospics.csv")

#cleaning up depths
mean(finalspics$depth) #much deeper than nospics -- very interesting and expected!
mean(finalnospics$depth) #much shallower
t.test(finalspics$depth, finalnospics$depth) #quite significant. But anyway:

quantile(finalspics$depth)
quantile(finalnospics$depth)

finalcleanspics <- finalspics[(finalspics$depth < 100),]
finalcleannospics <- finalnospics[(finalnospics$depth < 100),]

#the rasters
library(raster)

pH <- raster("Present.Surface.pH.BOv2_2.asc")
temp <- raster("Present.Benthic.Max.Depth.Temperature.Mean.asc")
silicate <- raster("Present.Benthic.Max.Depth.Silicate.Mean.asc")

#combinations
#nospics and pH
plot(pH)
points(finalcleannospics$decimalLongitude, finalcleannospics$decimalLatitude, col="darkblue", pch=19)
finalcleannospics$pH <- extract(pH, finalcleannospics[,c("decimalLongitude","decimalLatitude")], method="simple")

#nospics and temp
plot(temp)
points(finalcleannospics$decimalLongitude, finalcleannospics$decimalLatitude, col="darkblue", pch=19)
finalcleannospics$temp <- extract(temp, finalcleannospics[,c("decimalLongitude","decimalLatitude")], method="simple")

#nospics and silicate
plot(silicate)
points(finalcleannospics$decimalLongitude, finalcleannospics$decimalLatitude, col="darkblue", pch=19)
finalcleannospics$silicate <- extract(silicate, finalcleannospics[,c("decimalLongitude","decimalLatitude")], method="simple")

#save master nospics file
#write.csv(finalcleannospics, "alldatanospics.csv", row.names=FALSE)

#spics and pH
plot(pH)
points(finalcleanspics$decimalLongitude, finalcleanspics$decimalLatitude, col="blue", pch=19)
finalcleanspics$pH <- extract(pH, finalcleanspics[,c("decimalLongitude","decimalLatitude")], method="simple")

#spics and and temp
plot(temp)
points(finalcleanspics$decimalLongitude, finalcleanspics$decimalLatitude, col="blue", pch=19)
finalcleanspics$temp <- extract(temp, finalcleanspics[,c("decimalLongitude","decimalLatitude")], method="simple")

#spics and silicate
plot(silicate)
points(finalcleanspics$decimalLongitude, finalcleanspics$decimalLatitude, col="blue", pch=19)
finalcleanspics$silicate <- extract(silicate, finalcleanspics[,c("decimalLongitude","decimalLatitude")], method="simple")

#save master spics file
#write.csv(finalcleanspics, "alldataspics.csv", row.names=FALSE)

#Analyses
#omit NAs from files first
spics$Type <- rep(1, nrow(spics))
nospics$Type <- rep(0, nrow(nospics))
alldat <- rbind(spics, nospics)
occurrence <- glm(alldat$Type ~ alldat$pH + alldat$temp + alldat$silicate + alldat$depth, family="binomial")

occurrence <- glm(alldat$Type ~ alldat$pH + alldat$temp + alldat$silicate, family="binomial") #not as good AIC
occurrence <- glm(alldat$Type ~ alldat$pH + alldat$temp + alldat$depth, family="binomial")

step(occurrence) #final model has no silicate

finalmodel <- glm(alldat$Type ~ alldat$pH + alldat$temp + alldat$depth, family="binomial")

#controlled model
library(glmmTMB)
alldat <- read.csv("CompleteSpongeData.csv")
finalmodel <- glmmTMB(Type ~ pH + temp + silicate + depth + (1 | order), family=binomial, dat=alldat)
summary(finalmodel)


#testing machine learning stuff

library(caret)

test <- createDataPartition(p=0.75, list=FALSE, y=alldat$Type)

indep_var = colnames(alldat) != "Type"
model_rf <- train(x= alldat[indep_var], y= alldat$Type, method = 'rf')

#library(caret)
#alldat <- read.csv("CompleteSpongeData.csv")
#testdat <- testdat[sample(nrow(testdat), size=1000), ]
#testdat <- alldat[sample(nrow(alldat), size=1000), ]
#newdat <- subset(testdat, select = -c(gbifID, class, family, genus, species, decimalLatitude, decimalLongitude))
#head(newdat, 10)

?predict
?RMSE
