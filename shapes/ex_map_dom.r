# Exemple de carte pour Dom
# By Steve Vissault
# November 24th 2014

require("ggplot2")
require("maptools")

# Importation des shapefiles 

# Shapefiles des lacs
lakes  <- readShapePoly("./water_area_qc.shp") # Lecture du shapefile
lakes@proj4string  <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") # Définit la projection géographique
lakes  <- fortify(lakes) # Transforme les polygones vers un dataframe compréhensible pour ggplot2

# Shapefiles des zones de végétations
veg_zone  <- readShapePoly("~/Documents/Spatial_Data/CANUSA_shapes/veg_zone.shp") # Lecture du shapefile
veg_zone@proj4string  <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") # Définit la projection géographique
veg_zone  <- fortify(veg_zone) # Transforme les polygones vers un dataframe compréhensible pour ggplot2

theme_set(theme_grey(base_size = 18)) # Augmente la taille de la police par défault des figures (ggplot2)

# Le fichier d'entrée "dat" doit ressembler à ca:

# 3 colonnes:
#################
# 1. x = longitude (degree, projection WGS84)
# 2. y = latitude (degree, projection WGS84)
# 3. État de la cellule (T,B,M ou R) (doit être de type facteur)

# Les colonnes 1, 2 corespondent au centroide de la cellule

# PS: J'ai une fonction qui reshape les données de sorties du modèle vers ce format. 

ggplot(dat) + # Ouvre la fenêtre ggplot2 et définit le jeu de données par défault
	geom_raster(aes(x,y,fill=States, order = rev(States))) + # Définit la couche (raster) des états (T,B,M ou R)
  	scale_fill_brewer(palette="Spectral") + # Définit la palette de couleur
  	geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",colour="dodgerblue4",size=0.1) + # Ajoute la couche des polygones ( de végétation) 
  	geom_polygon(data = veg_zone, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.3) + # Ajoute la couche des polygones (Zone de végétation) 
  	scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+ # Élimine certains espaces indésirables dans les axes
  	coord_equal() + # Fait en sorte que les carrés du raster ne soit pas rectangulaire 
  	xlab("Longitude") + ylab("Latitude")+ # Renomme les axes
  	theme(plot.title = element_text(lineheight=.8), panel.margin = unit(c(0.5,0,0,0),"in"),plot.margin = unit(c(0.1,0,0,0),"in")) # Ajuste les marges du graphiques