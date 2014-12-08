rm(list = ls())

# Read parameters
source("scripts/read_pars.R")
#pars = as.list(read_pars("data/pars.txt"))
#pars = as.list(read_pars("data/pars_v1.txt"))
#pars = as.list(read_pars("data/pars_all.txt"))
pars = as.list(read_pars("data/GenSA_test_sub25000.txt"))

#################################
## Prepare the panels          ##
#################################

#################################
## PANEL A: CLIMATE SPACE      ##
#################################

source("scripts/get_eq.R")

T = seq(-3,4,0.01)
P = seq(-2,2,0.01)
clim_space = expand.grid(T,P)

inv = matrix(nr = nrow(clim_space), nc = 2)
prob = matrix(nr = nrow(clim_space), nc = 4)

# Compute invasibility values
for(i in 1:nrow(clim_space)) {
#  inv[i,] = get_inv(ENV1 = clim_space[i,1],ENV2 = clim_space[i,2], pars)
  prob[i,] = get_eq(ENV1 = clim_space[i,1],ENV2 = clim_space[i,2], pars)[1:4]
  cat(round(i/nrow(clim_space)*100,2),'\n')
}

# Draw states
draw = function(p) {
  states = c(1,2,3,0)
  draw = rmultinom(n=1,size=1,prob=p)
  states[which(draw==1)]
}

StatesSTM = apply(prob,1,draw)

# Plot the results
Z = matrix(StatesSTM,nr = length(T), nc = length(P))
#Z = matrix(coexist,nr = length(T), nc = length(P))
quartz(height = 6, width = 6)
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("center",legend = c("Regeneration","Boreal","Temperate","Mixt"),fill = c("black","darkcyan","orange","palegreen3"),bty = "n",horiz = TRUE,cex = 1)
par(mar=c(5,5,0,2))
image(T*1.891697+1.529991 ,(P*132.858+1044.194),Z,xlab = "Annual mean temperature", ylab = "Precipitations (mm)", cex.lab = 1.5, cex.axis = 1.25, col = c("black","darkcyan","orange","palegreen3"))

dev.copy2eps(file = "Figs/STM_clim_space.eps")

#################################
## PANEL B: Map                ##
#################################

#################################
## Load libraries              ##
#################################

require("reshape2")
require("ggplot2")
require("grid")
require("maptools")
require("RColorBrewer")

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
climate_grid = subset(climate_grid,climate_grid[,3]>-4)
climate_grid[,3] = (climate_grid[,3]-1.529991)/1.891697
climate_grid[,4] = (climate_grid[,4]*1000-1044.194)/132.858

#################################
## Solve model at equilibrium  ##
#################################

eq = matrix(nr = nrow(climate_grid), nc = 4)
for(i in 1:nrow(climate_grid)) eq[i,] = get_eq(ENV1 = climate_grid[i,3],ENV2 = climate_grid[i,4], pars)[1:4]

# Clean negative values
eq[eq[,4]<0] = c(0,0,0,1)
eq[apply(eq,1,sum)==0,] = c(0,0,0,1)
StatesSTM = apply(eq,1,draw)

# Draw states
draw = function(p) {
	states = c("B","T","M","R")
	draw = rmultinom(n=1,size=1,prob=p)
	states[which(draw==1)]
}

States = apply(eq[,1:4],1,draw)

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

dev.copy2eps(file = "Figs/STM_map.eps")

