rm(list = ls())

# Read parameters
source("scripts/read_pars.R")
#pars = as.list(read_pars("data/pars.txt"))
pars = as.list(read_pars("data/pars_v0.txt"))


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

climate_grid = read.csv("data/climate_grid.csv")
climate_grid = subset(climate_grid,climate_grid[,4]!="-9999")
climate_grid = subset(climate_grid,climate_grid[,3]>-4)
climate_grid[,3] = (climate_grid[,3]-1.529991)/1.891697
climate_grid[,4] = (climate_grid[,4]*1000-1044.194)/132.858

#################################
## Compute eigenvalues            ##
#################################

p = matrix(nr = nrow(climate_grid), nc = 4)
maxEig = numeric(nrow(climate_grid))

# Compute initial conditions, final conditions and eigenvalues at final conditions
for(i in 1:nrow(climate_grid)) {

	EQ = get_eq(p0=p0,ENV1 = climate_grid[i,3],ENV2 = climate_grid[i,4], pars)
	maxEig[i] = max(EQ[5:7])

}

#################################
## Set the figure          ##
#################################
maxEig[maxEig>0] = -999
maxEig[maxEig==-999] = max(maxEig)
logEig = log(-maxEig)
gg_eig <- cbind(climate_grid[,c(1,2)],value=logEig)

plot_eig = ggplot(gg_eig) +
	geom_polygon(data = subset(region_qc,hole==FALSE), aes(x = long, y = lat, group = group),fill="grey90",colour="grey60",size=0.05) +
	geom_raster(aes(lon,lat,fill=factor(cut(value,11),rev(levels(cut(value,11)))))) +
#	facet_wrap(~variable) +
	scale_fill_manual(values=brewer.pal(11,"Spectral"),name="Resilience")+
	geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",color="deepskyblue2",size=0.1) +
	#geom_polygon(data = zone_veg, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.2) +
	scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
	coord_equal() +
	xlab("Longitude") + ylab("Latitude")+
	theme(plot.title = element_text(lineheight=.8,face="bold"),
	panel.background = element_rect(fill = "lightskyblue"),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_line(color = "grey90",size=0.1),
	strip.text.x = element_text(size = 12,face="bold" ,colour = "white"),strip.background = element_rect(colour="black", fill="black"))


quartz(height = 5, width = 6)
plot_eig

dev.copy2eps(file = "Figs/Eig_map.eps")





