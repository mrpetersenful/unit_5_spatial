### Lesson 5.2: Right whale critical habitats and mortalities
## New packages: sf
## New functions: st_read(), st_as_sf(), st_transform(), st_cast(), geom_sf(), coord


## Spatial vector data: Points, lines, polygons

# In the last class we examined rasters, or gridded spatial data. The other way we typically deal with spatial data is through vector data, which consists of points, lines, and polygons. The most common file format for passing along spatial vector data is the shapefile.

# points: a single x, y coordinate pair (tree, city)
# lines: 2 or more points that are connected (stream, road)
# polygons: 3 or more points that form a closed shape (lake, country)

# Geospatial data in vector format are often stored in a shapefile format. Shapefiles were developed by Esri (the ArcGIS company). Because the structure of points, lines, and polygons are different, each individual shapefile can only contain one vector type (all points, all lines, or all polygons). You will not find a mixture of point, line, and polygon objects in a single shapefile. Shapefiles are actually a collection of files with a common filename prefix and stored in the same directory. In R, the whole set of shapefiles with a common filename prefix can be read in simultaneously with a single load command. 

# Objects stored in a shapefile often have a set of associated attributes that describe the data. For example, a line shapefile that contains the locations of streams, might contain the associated stream name, and other information about each stream line object. 


## sf package

# Let's start y loading the libraries we are going to use. The new sf package gives us tools to work with spatial data in data.frame and tibble formats to play nicely with the tidyverse workflows we have been practicing. An sf (simple features) object is a data.frame (or tibble) with rows of features, columns of attributes, and a special geometry column that contains the spatial aspects of the features. The geometry columns contain a list of sfc (simple features column), which is made up of individual objects of class sfg (simple features geometries) and a CRS (Coordinate Reference System). The sfg objects each represent the geometry of a single feature and contains information about the feature's coordinates, dimension and type of geometry. (Most spatial data, including everything we'll work with in class, is 2-dimensional with X and Y coordinates). The geometry type indicates the shape of a feature. Here are the most common geometry types: 

# Geometry types:

# Point: a single point
# Multipoint: multiple points
# Linestring: sequence of two or more points connected by straight lines
# Multilinestring: multiple lines
# Polygon: a closed ring with zero or more interior holes
# Multipolygon: multiple polygons
# Geometrycollection: any combination of the above types

# All functions and methods in the sf package that operate on spatial data are preceded by the prefix st_ which refers to spatial type. 



## North Atlantic right whales

# We will practice using spatial vector data by recreating a figure I made for a policy-oriented paper on the critically endangered North Atlantic right whale. Right whales are heavily impacted by human activities, and are frequently killed by ship strikes or entanglement in fishing gear. For decades, US and Canadian maritime agencies (especially NOAA and the DFO) implement regulations to reduce ship speeds, modify fishing gear and limit fishing activity in right whale hot spots. In 2017, right whales experienced an Unexpected Mortality Event when 17 carcasses were discovered over the course of a single year. This was the result of an abrupt shift in right whale distribution in response to climate-driven changes to prey availability. By mapping the boundaries of right whale critical habitats in the US and Canada alongside the location of right whale carcasses, we can assess the efficacy of fishery and shipping regulations, and the failure of these regulations in 2017. 

install.packages("sf")
library(sf)
library(broom)
library(tidyverse)
library(mapdata)
library(marmap)


## Reading and organizing disparate spatial datasets































