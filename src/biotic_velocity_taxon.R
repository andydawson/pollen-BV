#setwd('C:/Users/abrow/Documents/pg-pollen')
# require(tidyr)
require(ggplot2)
# require(rasterVis)
# require(fields)
# require(rgdal)
require(raster)
# require(enmSdm)
# require(rgeos)
# require(sp)
# require(dplyr)
# require(holoSimCell)
# require(gridExtra)
# require(ggrepel)
# require(parallel)

######################################################################################################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
######################################################################################################################

# read in grid 
locs_grid <- readRDS('data/grid_4.1.RDS')

# see what it looks like
head(locs_grid)

# read in list of taxa that we have maps for
taxa = readRDS("data/taxa_4.1.RDS")

taxa

# pick one to look at
tx = 10

# read in the calculated biotic velocities
bv = readRDS(paste0('output/',taxa[tx],'_bvs_n200_v4.1.RDS'))
bv$taxon=taxa[tx] 

# see what it looks like
head(bv)

# names of the columns
colnames(bv)

# this object is the output of a function called bioticVelocity from the R package enmSdm
# probably not a good idea to try install this (at least yet); it has many dependencies
# this means that to get it working you also need to install many other packages
# read about the function on the github page for the package
# https://rdrr.io/github/adamlilith/enmSdm/man/bioticVelocity.html

# dimension of data frame
# rows by columns
dim(bv)

# why so many rows?
# if you look at the timeFrom column you can see that some times repeat
bv$timeFrom
unique(bv$timeFrom)

# this is because for each taxon, the model generates many sets of maps
# so can calculate bv for each set
# so for example there are many rows with timeFrom = -21000 
which(bv$timeFrom == -21000)
# how many?
length(which(bv$timeFrom == -21000))
# so 200 bv values for this taxon for timeFrom = -21000
# so we can learn about the spread of the values (also called the uncertainty)



# plot some of the measures
# let's plot the centroid velocity measure
# which quantifies movement in the location of the center of the blob
ggplot(bv, aes(x=as.factor(timeFrom), y=centroidVelocity)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90,hjust=1))


ggplot(bv, aes(x=centroidVelocity)) + 
  geom_histogram() 


ggplot(bv, aes(x=as.factor(timeFrom), y=centroidLat)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90,hjust=1))


ggplot(bv, aes(x=centroidLong, y=centroidLat)) + 
  geom_point()


# let's plot the  measure
# which quantifies movement in the location of
ggplot(bv, aes(x=as.factor(timeFrom), y=nsCentroidVelocity)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90,hjust=1))


ggplot(bv, aes(x=centroidVelocity)) + 
  geom_histogram() 

