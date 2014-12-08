rm(list=ls())

#################################
## Load libraries              ##
#################################

require("reshape2")
require("ggplot2")
require("grid")
require("maptools")
require("RColorBrewer")
#require("rgeos")

#################################
## Load shapefiles             ##
#################################

zone_veg <- readShapePoly("shapes/zone_veg.shp")
zone_veg@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#lakes <- readShapePoly("shapes/lakes_qc.shp")
#lakes@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
region_qc <- readShapePoly("shapes/region_qc.shp")
region_qc@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#Subset sur l'aire des polygones (prend les lacs avec une superficie supérieur à 0.005 km²)
#area <- gArea(lakes, byid=TRUE)
#area <- area[area>0.005]
#lakes <- lakes[names(area),]

# Transforme pour ggplot2
#lakes <- fortify(lakes)
zone_veg <- fortify(zone_veg)
region_qc <- fortify(region_qc)
#rivers <- fortify(rivers)
load("shapes/map_shp.rObj")

#################################
## Climate data             ##
#################################

# Climate data
climate_grid = read.csv("data/climate_grid.csv")
climate_grid = subset(climate_grid,climate_grid[,4]!="-9999")
climate_grid[,4] = climate_grid[,4]*1000
climate_grid = subset(climate_grid,climate_grid[,3]>-4)
climate_grid[,5] = climate_grid[,3] + 4
climate_grid[,6] = climate_grid[,4]

#climate_grid[,3] = (climate_grid[,4]-1.529991)/1.891697
#climate_grid[,4] = (climate_grid[,5]*1000-1044.194)/132.858
#climate_grid[,3] = (climate_grid[,4]-1.529991)/1.891697
#climate_grid[,4] = (climate_grid[,5]*1000-1044.194)/132.858

#################################
## PANEL A: SDM PROJECTIONS    ##
#################################

#climate_grid <- subset(climate_grid,lat>=rg_lat[1] & lat<=rg_lat[2])
#climate_grid <- subset(climate_grid,lon>=rg_lon[1] & lon<=rg_lon[2])

# Run the SDM
data = read.csv("fit_model/transitionsFourState.csv")
selectedVars = c("annual_mean_temp2","annual_pp2")
datSel = data[,c("state2",selectedVars)]
datSel_wo_U <- subset(datSel, state2 != "U")
datSel_wo_U$state <- droplevels(datSel_wo_U$state2)
library(nnet)
SDM1 = multinom(state2 ~ .^2 + I(annual_mean_temp2^2) + I(annual_mean_temp2^3) + I(annual_pp2^2) + I(annual_pp2^3), data = datSel_wo_U[,c("state2",selectedVars)], maxit =300)

# Projection
dataProj = data.frame(annual_mean_temp2 = climate_grid[,5], annual_pp2 = climate_grid[,6])
projSDM= predict(SDM1,new=dataProj,"prob", OOB=TRUE)

# Draw values
draw = function(p) {
  states = c("B","M","R","T")
  draw = rmultinom(n=1,size=1,prob=p)
  states[which(draw==1)]
}
States = apply(projSDM,1,draw)

# Set the figure
colors_state = c("darkcyan","palegreen3","black","orange")
gg_prob_state <- cbind(climate_grid[,c(1,2)],States)

theme_set(theme_grey(base_size=16))
ggplot_state = ggplot(gg_prob_state) +
    geom_polygon(data = subset(region_qc,hole==FALSE), aes(x = long, y = lat, group = group),fill="grey90",colour="grey60",size=0.05) +
    geom_raster(aes(lon,lat,fill=States),alpha=0.9) +
    scale_fill_manual(values=colors_state,name="States")+
    geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",color="deepskyblue2",size=0.1) +
    #geom_polygon(data = zone_veg, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.2) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    coord_equal() +
    xlab("Longitude") + ylab("Latitude")+
    theme(plot.title = element_text(lineheight=.8,face="bold"),
        panel.background = element_rect(fill = "lightskyblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90",size=0.1))


quartz(height = 5, width = 6)
ggplot_state

dev.copy2eps(file = "Figs/SDM_map_2080.eps")

#################################
## PANEL B: STM PROJECTIONS    ##
#################################

source("scripts/read_pars.R")
source("scripts/sim_ODEs.R")
source("scripts/get_eq.R")

#pars = as.list(read_pars("data/pars.txt"))
pars = as.list(read_pars("data/pars_v0.txt"))

nsteps = 80

# Initial values
prob = matrix(nr = nrow(climate_grid), nc = 4)
pred_STM = matrix(nr = nrow(climate_grid), nc=4)

# Standardize values
climate_grid_std = climate_grid
climate_grid_std[,3] = (climate_grid[,3]-1.529991)/1.891697
climate_grid_std[,4] = (climate_grid[,4]-1044.194)/132.858
climate_grid_std[,5] = climate_grid_std[,3] + 4
climate_grid_std[,6] = climate_grid_std[,4]*1.01

for(z in 1:nrow(climate_grid_std)) {
  prob[z,] = get_eq(ENV1 = climate_grid_std[z,3],ENV2 = climate_grid_std[z,4], pars)[1:4]
  pred_STM[z,] = sim_ODEs(p0 = prob[z,], climStart=as.numeric(climate_grid_std[z,c(3,4)]),climEnd = as.numeric(climate_grid_std[z,c(5,6)]), nsteps, pars)
}

# Draw values
draw = function(p) {
  states = c("B","T","M","R")
  draw = rmultinom(n=1,size=1,prob=p)
  states[which(draw==1)]
}
pred_STM[is.na(pred_STM[,4]),]=c(0,0,0,1)
pred_STM[apply(pred_STM,1,sum)==0,] = c(0,0,0,1)
States = apply(pred_STM,1,draw)

gg_prob_state <- cbind(climate_grid[,c(1,2)],States)

# Do the map
theme_set(theme_grey(base_size=16))
ggplot_state = ggplot(gg_prob_state) +
    geom_polygon(data = subset(region_qc,hole==FALSE), aes(x = long, y = lat, group = group),fill="grey90",colour="grey60",size=0.05) +
    geom_raster(aes(lon,lat,fill=States),alpha=0.9) +
    scale_fill_manual(values=colors_state,name="States")+
    geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",color="deepskyblue2",size=0.1) +
    #geom_polygon(data = zone_veg, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.2) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    coord_equal() +
    xlab("Longitude") + ylab("Latitude")+
    theme(plot.title = element_text(lineheight=.8,face="bold"),
        panel.background = element_rect(fill = "lightskyblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90",size=0.1))


quartz(height = 5, width = 6)
ggplot_state

dev.copy2eps(file = "Figs/STM_map_2080.eps")






