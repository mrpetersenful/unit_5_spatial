## 24 March 2021
## Lesson 5.2: Right whale critical habitats and mortalities

### Spatial vector data: Points, lines, polygons

## In the last class we examined rasters, or gridded spatial data. The other way we
## typically deal with spatial data is through vector data, which consists of points,
## lines, and polygons. The most common file format for passing along spatial vector
## data is the shapefile. 
#### -- points: a single x,y coordinate pair (tree, city)
#### -- lines: 2 or more points that are connected (stream, road)
#### -- polygons: 3 or more points that form a closed shape (lake, country)

## Geospatial data in vector format are often stored in a shapefile format. Shapefiles
## were developed by Esri (the ArcGIS company). Because the structure of points, 
## lines, and polygons are diferent, each individual shapefile can only contain one
## vector type (all points, all lines, or all polygons). You will not find a mixture
## of point, line and polygon objects in a single shapefile. Shapefiles are actually
## a collection of files with a common filename prefix and stored in the same 
## directory. In R, the whole set of shapefiles with a common filename prefix can be
## read in simultaneously with a single load command. 

## Objects stored in a shapefile often have a set of associated attributes that 
## describe the data. For example, a line shapefile that contains the locations of
## streams, might contain the associated stream name, and other information about
## each stream line object.

### sf package

## Let's start by loading the libraries we are going to use. The new sf package gives
## us tools to work with spatial data in data.frame and tibble formats to play nicely
## with the tidyverse workflows we have been practicing. An sf (simple features)
## object is a data.frame (or tibble) with rows of features, columns of attributes 
## and a special geometry column that contains the spatial aspects of the features. 
## The geometry columns contains a list of class sfc (simple features column), which
## is made up of individual objects of class sfg (simple features geometries) and
## a CRS (Coordinate Reference System). The sfg objects each represent the geometry
## of a single feature and contains information about the feature's coordinates, 
## dimensions and type of geometry. (Most spatial data, including everything we'll
## work with in class, is 2-dimensional with X and Y coordinates). The geometry type
## indicates the shape of a feature. Here are the most common geometry types: 

#### --POINT: a single point
#### --MULTIPOINT: multiple points
#### --LINESTRING: sequence of two or more points connected by straight lines
#### --MULTILINESTRING: multiple lines
#### --POLYGON: a closed ring with zero or more interior holes
#### --MULTIPOLYGON: multiple polygons
#### --GEOMETRYCOLLECTION: any combination of the above types

## All functions and methods in the sf package that operate on spatial data are 
## preceded by the prefix st_ which refers to spatial type. 

### North Atlantic right whales

## We will practice using spatial vector data by recreating a figure Erin made for 
## a policy-oriented paper on the critically endangered North Atlantic right whale. 
## Right whales are heavily impacted by human activities, and are frequently killed
## by ship strikes or entanglement in fishing gear. For decades, US and Canadian
## maritime agencies (especially NOAA and the DFO) implement regulations to reduce
## ship speeds, modify fishing gear and limit fishing activity in right whale hot
## spots. In 2017, right whales experienced an Unexpected Mortality Event when 17
## carcasses were discovered over the course of a single year. This was the result
## of an abrupt shift in right whale distribution in response to climate-driven 
## changes to prey availability. By mapping the boundaries of right whale critical 
## habitats in the US and Canada alongside the location of right whale carcasses, we
## can assess the efficacy of fishery and shipping regulations, and the failure of 
## these regulations in 2017.


library(sf) # simple features (spatial vector data) st_read, st_transform
# library(broom) # part of tidyverse
library(tidyverse)
library(mapdata)  # map_data
library(marmap) # getNOAA.bathy()

### Reading and organizing disparate spatial datasets

## Now we are going to load in three right whale datasets. (This data is all publicly
## available, but it required me to combine multiple, messy, sources as well as data
## requests from the North Atlantic Right Whale Consortium.) We'll load:
#### 1. shapefiles that contain polygons corresponding to the boundaries of right
####     whale critical habitats in the US in 2017
#### 2. a .csv file that contains the lat/lon boundaries of right whale critical 
####     habitats in Canada in 2017
#### 3. a .csv file that contains the location and other info about right whale 
####     carcasses found in 2017

## To read in the US critical habitat shapefiles, we use sf:st_read(). We take a look
## at the data and see that the CRS is EPSG code 4269, which is the NAD83 projection.
## [The 4-digit EPSG (European Petroleum Survey Group Geodetic Parameter Dataset) 
## is a public registry of spatial reference systems. The NAD83 projection is the
## North American Datum 1983 CRS, which is most commonly used by US federal agencies.]
## In this exercise, we are going to use CRS WGS84 [World Geodetic System 1984; 
## EPSG 4326.] so we'll use st_transform() to transform our data from one CRS to 
## another. 

## First, let's read in the carcass location data. 
carcass = read.csv('data/RW_carcasses_2017.csv')

## Now, we'll read in US critical habitat shapefiles.
## https://www.greateratlantic.fisheries.noaa.gov/educational_resources/gis/data/index.html
## This will read in the set of shapefiles.
USA_crit_hab = st_read('data/North_Atlantic_Right_Whale_Critical_Habitat','North_Atlantic_Right_Whale_Critical_Habitat') 
USA_crit_hab
## Now we're transforming our CRS to be the same. CRS = "+proj=longlat +datum=WGS84")
USA_crit_hab_sf = st_transform(USA_crit_hab, crs=4326)#crs="+proj=longlat +datum=WGS84")
class(USA_crit_hab_sf)

## To get the Canadian critical habitat boundaries, we read in the series of polygon
## points provided in a simple .csv file. These are latitude/longitude coordinates 
## with the CRS WGS84. We can convert the data frame into a simple features data
## frame using st_as_sf(). To do this, we must tell the function st_as_sf() which
## columns should be used as the X and Y coordinates, and what the CRS is. This gives
## us a simple features data frame with point geometries. 

## Next we need to turn our points into polygons (i.e., the boundaries of our critical
## habitats). To do this we can use dplyr functions that have been created for simple
## features. We'll group our POINTS by habitat and country (grouping by country
## doesn't accomplish anything in this case since all of the points are in Canada, 
## but including it in the grouping allows us to retain this column when we use
## summarize()). Then we'll use summarize to collapse the POINTS into 2 MULTIPOINT
## geometry types (one for Grand Manan Basin, one for Rosewat Basin). As you recall
## from using summarize() previously, this function drops all columns except the 
## column(s) that the data have been grouped by and the column(s) that you create in
## the summarize() function. However, when using dplyr on sf objects, the geometries
## are "sticky", meaning that the geometry column of sf objects is retained unless 
## the user deliberately removes it. So we have used summarize() to convert the POINTs
## in each habitat to a single MULTIPOINTgeometry. Whew. Finally, we used the
## parameter do_union=FALSE. This triggers the summarize() function to use st_combine()
## instead of st_union() to collapse our POINTs into MULTIPOINTs. This only matters 
## because the st_union() function does not preserve POINT order, and if our
## POINTs get out of order, then the POLYGON that we create will be wonky. Finally,
## we use st_cast() to turn our MULTIPOINT geometries into POLYGONs. 


## Load in Canadian RW critical habitat coordinates 
## http://www.dfo-mpo.gc.ca/species-especes/profiles-profils/rightwhaleNA-baleinenoireAN-eng.html
CAN_crit_hab = read.csv('data/NARW_canadian_critical_habitat_2017.csv')
head(CAN_crit_hab)
class(CAN_crit_hab)

names(CAN_crit_hab)[names(CAN_crit_hab) == "ï..lat"] = "lat"

head(CAN_crit_hab)
names(CAN_crit_hab)

head(CAN_crit_hab_sf)


## Turn data frame into sf points, then sf polygon
CAN_crit_hab_sf = CAN_crit_hab %>% 
  ## Converts to shapefile.
  st_as_sf(coords=c("lon","lat"), crs=4326) %>% 
  dplyr::group_by(habitat, country) %>% 
  ## Collapses data into multipoint; do_union=FALSE prevents reordering points.
  dplyr::summarise(do_union=FALSE) %>% 
  st_cast("POLYGON") 
## This is weird, but my column name for latitude has some weird characters in front
## of it, so I had to rename it to make this work. 
print(CAN_crit_hab_sf) # 2 simple features, with habitat and country attributes

## Now let's get our USA habitat and our Canada habitat data ready to be joined. 
## We'll give the USA data the same attributes as the Canadian data: habitat and 
## country. Then the USA and Canadian habitat data have the same attributes, the same
## geometric types (POLYGON) and the same CRS. Now we can just use rbind() to combine 
## them into a single sf object. 

## Simple: USA_crit_hab data frame to match CAN_crit_hab
## Here's our Gulf of Maine habitat.
plot(USA_crit_hab_sf$geometry[1], axes=TRUE)

## Here's our Florida and Georgia habitat.
plot(USA_crit_hab_sf$geometry[2], axes=TRUE)

USA_crit_hab_sf$habitat=c("NMFS_GOM", "SEUS")
USA_crit_hab_sf$country="USA"

head(USA_crit_hab_sf)


## This drops all other variables from shapefile.
USA_crit_hab_sf = USA_crit_hab_sf %>% dplyr::select(country, habitat, geometry) 
## Here, I had to specify that I wanted to use the dplyr function select() before 
## it would work.

## Join the USA and Canada critical habitat sf objects
crit_hab = rbind(USA_crit_hab_sf, CAN_crit_hab_sf)

### Making maps

## So that was a lot of hard work merging our USA and Canada datasets because their
## original formats were very different. Different file types (.csv vs. shapefile), 
## different geometries (coordinates vs. polygons), different attributes, and different
## CRSs. Now that we have combined these disparate datasets, we can use ggplot2 tools
## to turn this spatial vector data into maps!

## First I'm going to set GOM and GSL limits. 
lon_bounds = c(-72, -54)
lat_bounds = c(39, 53)

## And add in coastline data.
world_map = map_data("worldHires", ylim = lat_bounds, xlim = lon_bounds)

## Now I want to plot critical habitats and carcass locations.
ggplot()+
  ## Add coastline.
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "black") + 
  geom_sf(data=crit_hab, alpha = 0.5, aes(fill=country)) +
  geom_point(data = carcass, aes(x = Longitude, y = Latitude, 
                                 color = Carcass_position), size=2) + 
  ## Crop map edges.
  coord_sf(1.3, xlim = lon_bounds, ylim = lat_bounds) + 
  ylab("Latitude") + xlab("Longitude") + theme_classic()

## We can use what we learned in the previous tutorial to add bathymetry data beneath
## our critical habitat polygons and our carcass data. 


## library(ggnewscale)  
## new_scale_fill() add to ggplot to create 2nd fill scale
bath_m_raw = marmap::getNOAA.bathy(lon1 = lon_bounds[1], lon2 = lon_bounds[2],
                                   lat1 = lat_bounds[1], lat2 = lat_bounds[2], 
                                   resolution = 4) # resolution default: 4 minutes
## Now I want to convert bathymetry to data frame.
bath_m_fortify = marmap::fortify.bathy(bath_m_raw) 
bath_m = bath_m_fortify %>%
  ## This if else function is keeping only stuff below 0m.
  mutate(depth_m = ifelse(z>0, NA, z)) %>%
  ## Now I'm going to only keep depth information (but negative depth).
  dplyr::select(-z)

head(bath_m)
summary(bath_m)

## Now to plot critical habitats and carcass locations over bathymetry.
rw_map = ggplot()+
  geom_raster(data = bath_m , aes(x = x, y = y, fill = depth_m), alpha=0.75) + 
  scale_fill_gradientn(colors=c("black", "navy", "blue4","lightblue"), 
                       values = scales::rescale(c(-5000, -3000, -300, 0)), 
                       name="Depth (m)") +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "black", color = NA) +
  geom_sf(data=crit_hab, alpha = 0.5, fill='yellow') +
  geom_point(data = carcass, aes(x = Longitude, y = Latitude, color = Carcass_position), size=2) + 
  coord_sf(1.3, xlim = lon_bounds, ylim = lat_bounds) + # Crop map edges
  ylab("Latitude") + xlab("Longitude") + theme_classic()

## Now to view and save the map.
rw_map
ggsave(rw_map, filename='figures/RW_habitats.pdf', device="pdf", height=5, width=7)




