### Lesson 5.1: Bathymetry and chlorophyll rasters
## New packages: raster, mapdata, marmap
## New functions: raster(), rasterToPoints(), geom_raster(), geom_tile(), map_data(), raster::crop(), raster::extent(), marmap::getNOAA.bathy(), geom_polygon(), coord_fixed(), marmap::fortify.bathy(), marmap::as.raster(), geom_contour(), raster::resample(), raster::stack()

## Raster data

# One of the most common formats for spatial data is raster data. Think of a raster as a flat, 2-D matrix, where each cell is a pixel in your image. Rasters are gridded data products with an extent (the boundaries of your map), an origin (the point closest to (0,0) on your grid), and a resolution (the size of each box/pixel/grid cell in your data). 

# Optionally, your raster can also have a CRS (Coordinate Reference System) associated with it. The CRS ties your raster data to geographic coordinates on the planet Earth. The CRS provides information on how your 2-D data have been flattened from the bumpy 3-D ellipsoid of the Earth's surface. There is a LOT of important information related to understanding CRS, but we won't have time to dive into the world of map projections and transformations in this tutorial. 

# Scientists often have multiple variables attached to each spatial pixel. For example, you could grab satellite chlorophyll data for every day of the month. That means you have 30 rasters, all with the same extent, origin, resolution, and CRS. The easiest way to work with multiple raster variables is to stack them on top of each other so the coordinates line up. This is called a raster stack or a raster brick. 

# Once the rasters are stacked, it's easy to do raster math. For example, you could generate a monthly mean chlorophyll map by taking the mean across all 30 days at each pixel. You can also easily grab all of the variables you need at a specific (x,y). For example, you can plot the time series of chlorophyll over the span of the month at a given point. 


## Satellite remote sensing data

# Since the latter half of the 20th century, satellite remote sensing has become a major data collection pathway in the Earth sciences. Ironically, some of the best and most prolific data available about our planet is collected from space. As a data user, you can access satellite data products at different stages of processing: 

# Data product levels:

# L0: raw data, no geolocation
# L1: remote sensing reflectance - the light in a given frequency band geolocated to the Earth's surface, after accounting for sensor corrections. 
# L2: geophysical products - calculate water leaving-radiances after atmospheric corrections have been applied. Using documented algorithms, these are scientific variables of interest (i.e. chl, sst, cdom, etc.) that have been estimated from the remote sensing reflectance data collected by the satellite. Gridded data are provided at the resolution observed by the satellite. 
# L3: binned/mapped derived products - this is L2 data that has been binned and mapped to create a gridded product with uniform resolution. 

# Researchers will access different levels of products depending on their research interests. Researchers use L0 or L1 products when they are trying to derive new data products, improve existing data products, update atmospheric corrections, or develop methods to calibrate derived products across different satellite missions and instruments. Spatial resolutions viewed from satellite are quite irregular, primarily depending on the angle of the satellite from your location of interest at the time of observation. Researchers who are working on fine spatial scales often work with L2 products so they know the exact resolution and positioning of the pixels they are interested in. L3 products are convenient for researchers working at middle- and large- scales that don't need extremely fine scale spatial information. 

# The images below demonstrate how satellites collect data based on a fixed arclength, or the angle over which data is aggregated at the satellite. Because the arclength is constant, the data collected directly below the satellite (a.k.a. "nadir") has higher spatial resolution than the data collected off-nadir, or not directly below the satellite. 


## Collecting, plotting and manipulating raster data

# In this lesson we will download data from the NASA Ocean Biology Processing Group (OBPG). The OBPG creates popular data products like chlorophyll a and sea surface temperature. Here is a list of products available:

# https://oceancolor.gsfc.nasa.gov/atbd/

# Most satellite data (including all the data from NASA) is free and publicly available. To access the data from the OBPG, you need to create a NASA Earth Data Login profile. You may have already done this in the climate unit: https://urs.earthdata.nasa.gov/

# We can go to the Ocean Data Access Portal: https://oceandata.sci.gsfc.nasa.gov/

# We can collect data from the Modis Aqua mission, which has been running from 2002-present. (If you want to go a bit further back in time, you can get some of the same products from the SeaWiFS mission which ran from 1997-2010). We'll grab Modis Aqua L3 Mapped data that provides a monthly climatology of chlorophyll a. So from the Ocean Data Access Portal, click through the directories:

# MODIS-Aqua -> Mapped -> Monthly_Climatology -> 9_km -> chlor_a

# There you'll see a list of files with names such as A20021822017212.L3m_MC_CHL_chlor_a_9km.nc

# Let's dismantle the file naming convention here: 

# A: Satellite sensor (A for the Modis Aqua mission)
# 2002: start year
# 182: start Julian date (first day of July)
# 2017: end year
# 212: end Julian date (last day of July)
# .L3m: Level 3 monthly
# MC: Monthly Climatology
# CHL: Chlorophyll suite of products
# chlor_a: The specific chlorophyll product (distinct from chlor_ocx which is a legacy algorithm)
# 9km: resolution; the length of the side of one grid cell
# .nc: netCDF file type

# Note: the OBPG is in the process of transitioning naming conventions: https://oceancolor.gsfc.nasa.gov/docs/filenaming-convention/

# For this lesson, let's download and play with this file A20021822017212.L3m_MC_CHL_chlor_a_9km.nc, which is the monthly climatology of chlorophyll a for July (Julian dates 182-212) averaged across the years 2002-2017. This provides chlorophyll a data in mg/m^3. This data, and many other products in government DAACs (Distributed Active Archive Centers), comes in a NetCDF file. NetCDF files contain scientific array-type data and associated metadata. 

library(tidyverse)
install.packages("raster")
library("raster")
install.packages("mapdata")
library(mapdata)
install.packages("marmap")
library(marmap)
library(ncdf4)


# We can use R's raster package to read in the NetCDF file and take a look at the metadata. Let's change the name of the chlorophyll layer because the original name is sooooo long and too hard to type. Then we can convert the raster to a data frame because we are already so comfortable working with data frames. 











