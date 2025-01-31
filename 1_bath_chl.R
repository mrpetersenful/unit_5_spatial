## 22 March 2021
## Lesson 5.3: Bathymetry and chlorophyll data

### Raster data

## One of the most common formats for spatial data is raster data. Think of a raster
## as a flat, 2-D matrix, where each cell is a pixel in your image. Rasters are 
## gridded data products with an extent (the boundaries of your map), an origin (the
## point closest to (0,0) on your grid), and a resuolution (the size of each box/
## pixel/grid cell in your data). 

## Optionally, your raster can also have a CRS (Coordinate Reference System) associated
## with it. The CRS ties your raster data to geographic coordinates on the planet 
## Earth. The CRS provides information on how your 2-D data have been flattened from
## the bumpy 3-D ellipsoid of the Earth's surface. There is a LOT of important 
## information related to understanding CRS, but we won't have time to dive into the
## world of map projections and transformations in this tutorial. 

## Scientists often have multiple variables attached to each spatial pixel. For example, 
## you could grab satellite chlorophyll data for every day of the month. That means 
## you have 30 rasters, all with the same extent, origin, resolution and CRS. The 
## easiest way to work with multiple raster variables is to stack them on top of 
## each other so the coordinates line up. This is called a raster stack or a raster
## brick. 

## Once the rasters are stacked, it's easy to do raster math. For example, you could 
## generate a monthly mean chlorophyll map by taking the mean across all 30 days at
## each pixel. You can also easily grab all of the variables you need at a specific 
## (x,y). For example, you can plot the time series of chlorophyll over the span of
## the month at a given point. 

### Satellite remote sensing data

## Since the latter half of the 20th century, satellite remote sensing has become a
## major data collection pathway in the Earth sciences. Ironically, some of the 
## best and most prolific data available about our planet is collected from space. 
## As a data user, you can access satellite data products at different stages of 
## processing: 

## Data product levels
#### -- L0: raw data, no geolocation
#### -- L1: remote sensing reflectance - the light in a given frequency band geolocated
####      to the Earth's surface, after accounting for sensor corrections. 
#### -- L2: geophysical products - calculate water leaving-radiances after atmospheric
####      corrections have been applied. Using documented algorithms, these are 
####      scientific variables of interest (i.e., chl, sst, cdom, etc.) that have 
####      been estimated from the remote sensing reflectance data collected by the 
####      satellite. Gridded data are provided at the resolution observed by the 
####      satellite. 
#### -- L3: binned/mapped derived products - this is L2 data that has been binned
####      and mapped to create a gridded product with uniform resolution. 

## Researchers will access different levels of products depending on their research
## interests. Researchers use L0 or L1 products when they are trying to derive new 
## data products, improve existing data products, update atmospheric corrections or 
## develop methods to calibrate derived products across different satellite missions
## and instruments. Spatial resolutions vied from satellite are quite irregular, 
## primarily depending on the angle of the satellite from your location of interest 
## at the time of observation. Researchers who are working on fine spatial scales 
## often work with L2 products so they know the exact resolution and positioning of
## the pixels they are interested in. L3 products are convenient for researchers 
## working at middle- and large- scales that don't need extremely fine scale spatial
## information. 

## Satellites collect data based on a fixed arclength, or the angle over which data
## is aggregated at the satellite. Because the arclength is constant, the data 
## collected directly below the satellite (aka "nadir) has higher spatial resolution
## than the data collected off-nadir, or not directly below the satellite. 


### Collecting, plotting and manipulating raster data

## In this lesson, we will download data from the NASA Ocean Biology Processing Group 
## (OBPG). The OBPG creates popular data products like chlorophyll a and sea surface
## temperature.

## Most satellite data (including all the data from NASA) is free and publicly 
## available. To access the data from the OBPG, you need to create a NASA Earth 
## Data Login profile. 

## We can collect data from the Modis Aqua mission, which has been running from 
## 2002-present. (If you want to go a bit farther back in time, you can get some 
## of the same products from the SeaWiFS mission which ran from 1997-2010). We'll 
## grab Modis Aqua L3 Mapped data that provides a monthly climatology of chlorophyll
## a. So from the Ocean Data Access Portal, click through the directories: 

## MODIS-Aqua -> Mapped -> Monthly_Climatology -> 9_km -> chlor_a

## There you'll see a list of files with names such as 
## A20021822017212.L3m_MC_CHL_chlor_a_9km.nc

## Let's dismantle the file naming convention here: 
## A: Satellite Sensor (A for the Modis Aqua mission)
## 2002: Start year
## 182: start Julian date (first day of July)
## 2017: end year
## 212: end Julian date (last day of July)
## .L3m: Level 3 monthly
## MC: Monthly Climatology
## CHL: Chlorophyll suite of products
## chlor_a: The specific chlorophyll product (distinct from chlor_ocx which is a 
##      legacy algorithm)
## 9km: resolution; the length of the side of one grid cell
## .nc: netCDF file type

## Note: the OBPG is in the process of transitioning naming conventions. 

## For this lesson, let's download and play with that file, which is the monthly 
## climatology of chlorophyll a for July (Julian dates 182-212) averaged across the
## years 2002-2017. This provides chlorophyll a data in mg/m^3. This data, and many
## other products in government DAACs (Distributed Active Archive Centers), comes 
## in a NetCDF file. NetCDF files contain scientific array-type data and associated
## metadata. 


library(tidyverse)
library(raster)
library(mapdata)  # map_data worldHires coastline
library(marmap)   # getNOAA.bathy()
# library(ncdf4)  # also good for reading netCDF files

## We can use R's raster package to read in the NetCDF file and take a look at the
## metadata. Let's change the name of the chlorophyll layer because the original 
## name is so long and too hard to type. Then we can convert the raster to a data 
## frame because we are already so comfortable working with data frames. 


## Modis Aqua monthly climatology of chlor_a for July (Julian dates 182-212) 
## averaged over 2002 - 2017
## https://oceandata.sci.gsfc.nasa.gov/
## Data: Modis Aqua -> Mapped -> Monthly Climatology -> 4km -> chlor_a
chl_raster = raster('data/A20021822017212.L3m_MC_CHL_chlor_a_9km.nc') 
# chl_a in mg/m^3 # raster reads in the largest array by default

## First, we want to take a look at the data. 
class(chl_raster)
## Now let's take a look at the raster attributes.
chl_raster

## I'm going to rename my raster layer because it's easier to type:
names(chl_raster) = "chl_a" 

## NOw I'm goin to convert this to a data frame, and then convert to a Spatial 
## Points data frame. 
chl_pts = raster::rasterToPoints(chl_raster, spatial = TRUE) 
chl_df  = data.frame(chl_pts)
head(chl_df)

## I also want to see a range of data to set good limits for the color palette that
## I want to use when I'm creating maps. 
hist(log10(chl_df$chl_a))

## I want to use a rainbow palette. This will reverse the order, and drop the first
## color, which is purple. 
cols = rainbow(7, rev=TRUE)[-1]

## Now, I'm going to plot the global chlorophyll. geom_raster() plots faster, 
## geom_tile() is slightly more precise.

global_chl_map = ggplot() +
  ## Here, I'm using the log10 scale to view the chlorophyll a data, because that's
  ## a typical way to view chl a data. 
  geom_raster(data = chl_df , aes(x = x, y = y, fill = log10(chl_a))) +
  ggtitle("Global chl_a July climatology") +
  theme_classic() +
  scale_fill_gradientn(colors = cols, limits=c(-1.5, 0.75), name="log_10(chl_a)")

## Now view the map, and save. 
global_chl_map
ggsave(global_chl_map, filename='figures/global_chl_July.pdf', device="pdf", 
       height=5, width=9)

## That's pretty cool. There's a lot more chlorophyll in the cold nutrient-rich 
## subpolar waters. And you can see the upwelling zones along the coasts and at
## the equator. 

## Let's zoom in on an area of interest. A lot of Erin's work is in the Gulf of 
## Maine, and there's interesting bathymetry there, too. We can crop a region and
## re-plot the chlorophyll to see more detail. 

## First, let's set the Gulf of Maine map limits.
lon_bounds = c(-72, -62)
lat_bounds = c(39, 47)

## Now, I want to crop the Gulf of Maine data.
chl_GOM_raster = raster::crop(chl_raster, extent(c(lon_bounds, lat_bounds))) 

## Now I want to convert the Gulf of Maine raster data to points (using raster package)
## and then to a data frame. 
chl_GOM_df = data.frame( rasterToPoints(chl_GOM_raster, spatial = TRUE) )
head(chl_GOM_df)
## Now I want to drop the optional column. 
chl_GOM_df = chl_GOM_df %>% dplyr::select(-optional) 

## Now let's grab coastline data from R's WorldHires data in the mapdata package. 
world_map = map_data("worldHires")
## Or use GSHHS Shore line maps for higher resolution data:
## https://eriqande.github.io/2014/12/17/compare-map-resolutions.html

## Now I'm going to map the chlorophyll in the Gulf of Maine.
GOM_chl_map = ggplot() +
  geom_raster(data = chl_GOM_df , aes(x = x, y = y, fill = log10(chl_a))) + 
  ## Now we're adding coastline.
  geom_polygon(aes(x=long, y = lat, group = group), fill = "darkgrey", data=world_map) +
  ## Now we're cropping the map.
  coord_fixed(ratio=1.3, xlim = lon_bounds, ylim = lat_bounds, expand=FALSE) + 
  ggtitle("GOM chl_a July climatology") +
  theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  ## Changing color limits from the global chl a map. 
  scale_fill_gradientn(colors = cols, limits=c(-1, 1.75))  


## Now I want to view and save the image. 
GOM_chl_map
ggsave(GOM_chl_map, filename='figures/GOM_chl_July.pdf', device="pdf", 
       height=5, width=7)

## Great. Now we've got a nice cropped image of the Gulf of Maine with coastlines 
## and the monthly climatology for chlorophyll in mg/m^3 for the month of July from
## 2002-2017. Note that geom_raster makes a nice smooth image, but if you use 
## geom_tile() you'll see a more pixelated (more honest) map of the chlorophyll 
## data, and this may be what you should use in your map, depending on the level 
## of precision you need. We could have downloaded the 4km resolution data instead,
## and it would have just taken a few more seconds to process. If I was putting this 
## in a publication, I'd grab the 4km resolution and use geom_tile(). 

## Using the exact same steps, you could grab daily or monthly values for a specific 
## year corresponding to the time period of some other data that you have collected. 
## You could grab the JUly monthly data each year (rather than the monthly 
## climatology averaged across all years) tod see how chlorophyll has changed over
## time. You could also grab sea surface temperature or a host of other data products
## from the same website as well. It's that easy. 


### Bathymetry

## Now let's download and plot some bathymetry data from NOAA. Some oceanographers 
## in France developed a really nice package to automate reading in NOAA bathymetry
## data from R using the marmap package. This makes bringing bathymetry data into 
## your R environment much easier. With marmap, the default resolution of our 
## raster grid is 4 minutes for the length of the side of a single pixel. You can use
## marmap to download data with a resolution as high as 1 minute. 

## If you happen to need higher resolution bathymetry data for your research, you
## can download netCDF files from GEBCO, which has a global bathymetry and 
## topography resolution of 15 arcseconds: https://download.gebco.net/

## We'll use the marmap package to download and load in our bathymetry data. We'll
## use fortify() to convert the marmap data from its native bathy class to a data
## frame for plotting with ggplot.


## NOAA bathymetry (from marmap package)
## We'll use the same Gulf of Maine map limits. 
lon_bounds = c(-72, -62)
lat_bounds = c(39, 47)

bath_m_raw = marmap::getNOAA.bathy(lon1 = lon_bounds[1], lon2 = lon_bounds[2],
                                   ## Here, I'm saying that the resolution default 
                                   ## should be four minutes. 
                             lat1 = lat_bounds[1], lat2 = 47, resolution = 4) 
class(bath_m_raw) 

## Now I'm going to convert my bathymetry data into a data frame, and use fortify() 
## to convert the marmap data from its native bathy class to a data frame for 
## plotting with ggplot.
bath_m_df = marmap::fortify.bathy(bath_m_raw) 
bath_m = bath_m_df %>%
  ## Here, I'm telling the data that if it's closer to coast than 20 m, then that 
  ## value is NA. 20 m gives us wiggle room for the tides/coastline. 
  mutate(depth_m = ifelse(z>20, NA, z)) %>% 
  dplyr::select(-z)

head(bath_m)
summary(bath_m)

## Now I want to plot the bathymetry map of the Gulf of Maine. 
GOM_bath_map = ggplot()+
  geom_raster(data = bath_m , aes(x = x, y = y, fill = depth_m)) + 
  ## Adding the coastline; group keeps multipolygon coordinates separated into 
  ## distinct groups.
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "darkgrey", color = NA) + 
  ## Here I want to crop the map edges.
  coord_fixed(1.3, xlim = lon_bounds, ylim = lat_bounds, expand=FALSE) + 
  ## This rescales to make 2 different gradients.
  scale_fill_gradientn(colors=c("black", "darkblue", "lightblue"), 
                       values = scales::rescale(c(-6000, -300, 0)), 
                       name="Depth (m)") +
  ylab("Lat") + xlab("Lon") + theme_bw() 

## Print and save.
GOM_bath_map
ggsave(GOM_bath_map, filename='figures/GOM_bath_raster.pdf', device="pdf", 
       height=5, width=7)


## For bathymetry data, it can be nice to add contour lines to the plot. Take a look 
## at these maps with and without the underlying colored bathymetry data.

# plot contours
GOM_bath_map_contours = ggplot()+
  ## Add 100m contour.
  geom_contour(data = bath_m, aes(x=x, y=y, z=depth_m), breaks=c(-100), 
               size=c(0.25), colour="grey") + 
  ## Add 200m contour.
  geom_contour(data = bath_m, aes(x=x, y=y, z=depth_m), breaks=c(-200), 
               size=c(0.5), colour="grey") +
  ## Add 500m contour.
  geom_contour(data = bath_m, aes(x=x, y=y, z=depth_m), breaks=c(-500), 
               size=c(0.75), colour="grey") + 
  ## Add coastline data.
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "black", color = NA) + 
  ## Crop map edges.
  coord_fixed(1.3, xlim = lon_bounds, ylim = lat_bounds, expand=FALSE) + 
  ylab("Latitude") + xlab("Longitude") + theme_classic()

## Print and save. 
GOM_bath_map_contours 
ggsave(GOM_bath_map_contours, filename='figures/GOM_bath_contours.pdf', 
       device="pdf", height=5, width=7)


## Exercise 1.1: 
## Add bathymetric contour lines to the colored bathymetry raster map of the Gulf
## of Maine. Draw the contours at 50m, 250m, and 1000m depths. Are there any basins
## in the GOM that go as deep as 1000m?

## First I want to plot the contours. 

GOM_bath_map_contours = ggplot() +
  geom_raster(data = bath_m, aes(x=x, y=y, fill = depth_m)) +
  scale_fill_gradientn(colors= c("black", "darkblue", "lightblue"), 
                       values=scales::rescale(c(-6000, -300, 0)), 
                       name="Depth (m)") +
  geom_contour(data = bath_m, aes(x=x, y=y, z=depth_m), breaks=c(-50), 
               size=c(0.25), colour="black") + 
  geom_contour(data = bath_m, aes(x=x, y=y, z=depth_m), breaks=c(-250), 
               size=c(0.5), colour="black") + 
  geom_contour(data = bath_m, aes(x=x, y=y, z=depth_m), breaks=c(-1000), 
               size=c(0.75), colour="black") + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "black", color = NA) + 
  coord_fixed(1.3, xlim = lon_bounds, ylim = lat_bounds, expand=FALSE) + 
  ylab("Latitude") + xlab("Longitude") + theme_classic()

GOM_bath_map_contours


### Combine bathymetry and chlorophyll rasters

## We have downloaded and created gorgeous maps with NASA chlorophyll a data and
## NOAA bathymetry data. Now we can layer those two datasets, which will enable 
## us to investigate the relationship between chlorophyll and bathymetry. We can 
## do this by putting both datasets into the raster format, and then resampling 
## them to have the same Coordinate Reference System, extent, origin, and resolution.
## Then the cells (pixels) of our chlorophyll raster grid will line up perfectly 
## with the cells of our bathymetry raster grid. We can stack the raster layers on
## top of each other, and then it's simple to do raster math, run models, or simply
## convert them back into a data frame to use our favorite dplyr functions. 

## The GOM cropped data we made above:
class(chl_GOM_raster)
class(bath_m_raw)  

## So we have two different data types, raster and bathy. We need to make them both
## raster. 

bath_m_raster = marmap::as.raster(bath_m_raw)

## Note that the CRS is the same for both rasters (WGS84) -- they're both American
## data sets. 
chl_GOM_raster
bath_m_raster

## Rename the bathymetry raster layer so its easier to work with.
names(bath_m_raster) = "bath_m"

## Now I want to resample the bath_m to match chl_a. This resamples the raster 1
## extent, origin, and resolution to that of raster 2.
bath_layer_chl_dims = raster::resample(bath_m_raster, chl_GOM_raster)

## If the CRS didn't match up, we'd also project the bath_m raster to match the CRS
## of the chl raster using the following line of code: 

# bath_layer_chl_dims_proj = raster::projectRaster(bath_m_raster, crs = crs(chl_GOM_raster))


## Now that our extent, origin, resolution, and projection match, create the raster
## stack. 
raster_stack = stack(chl_GOM_raster, bath_layer_chl_dims)
raster_stack

## We want to double check that everything is oriented correctly. 
plot(raster_stack) 

## Now that the rasters are all matched up, we can easily do science with them. 
## Converting the raster stack into a data frame (now that all of the spatial chores
## have been taken care of) allows us to use dplyr tools. 


## Convert to a data frame. 
stack_df = data.frame( raster::rasterToPoints(raster_stack))
head(stack_df)
summary(stack_df)
dim(stack_df)

## From O'Reilly et al. 2019, there are chlorophyll a benchmarks for oligo-, meso-
## and eutrophic ocean waters derived from SeaWiFS data. When chlorophyll a is less
## than 0.1 mg/m^3, it's oligotrophic, when it's greater than 1.67 mg/m^3, it's 
## eutrophic. 


oligo_chl_a = 0.1 
eutro_chl_a = 1.67 

## Let's add another column to our data to designate trophic status. 
stack_df = stack_df %>%
  mutate(trophic_index = case_when(chl_a < oligo_chl_a ~ "oligotrophic",
                                   chl_a > oligo_chl_a & chl_a < eutro_chl_a ~ "mesotrophic",
                                   chl_a > eutro_chl_a ~ "eutrophic")) %>%
  mutate(trophic_index = as.factor(trophic_index))

## What portion of our area of interest is classified as oligotrophic, mesotrophic
## and eutrophic?

table(stack_df$trophic_index)  
## This shows that we don't have any oligotrophic regions in the Gulf of Maine. 

n_eutro = dim(stack_df %>% filter(trophic_index == "eutrophic"))[1]
n_meso = dim(stack_df %>% filter(trophic_index == "mesotrophic"))[1]
pct_eutro = n_eutro / (n_eutro + n_meso)
pct_meso = n_meso / (n_eutro + n_meso)

## Now, let's plot a histogram of bathymetric depth (m) for each trophic index.
ggplot() +
  geom_histogram(aes(x=bath_m), data=stack_df %>% filter(!is.na(trophic_index))) +
  facet_wrap(~trophic_index)

## Now, I'm going to plot trophic raster data.
trophic_map = ggplot()+
  geom_raster(data = stack_df, aes(x = x, y = y, fill = trophic_index)) + 
  ## Adding in coastline; group keeps multipolygon coordinates separated into 
  ## distinct groups. 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "darkgrey", color = NA) + 
  ## Cropping the map edges.
  coord_fixed(1.3, xlim = lon_bounds, ylim = lat_bounds, expand=FALSE) +
  ylab("Lat") + xlab("Lon") + theme_bw() 
trophic_map

## So the eutrophic waters in the GOM (for July climatology) are the shallow waters
## near the coast. The rest of the Gulf of Maine is mesotrophic, and there are no
## oligotrophic areas in the GOM. 

## Exercise 1.2: 
## Using the raster stack we created earlier, crop out Cape Cod Bay. Turn this 
## into a data frame. How does the mean chlorophyll a concentration at depths between
## 0 and -50m compare to the mean chlorophyll concentration between -50 and -100m?

CCB_lon = c(-70.8, -69.8)
CCB_lat = c(41.6, 42.2)

CCB_raster = crop(raster_stack, extent(c(CCB_lon, CCB_lat)))
CCB_df = data.frame(raster::rasterToPoints(CCB_raster))

ggplot(aes(x=chl_a), data=CCB_df) +
  geom_histogram()

## Now let's find the mean chl a between depths of 0 and 50m. 
CCB_df %>%
  filter(bath_m <=0, bath_m > -50) %>%
  summarize(chl = mean(chl_a, na.rm=TRUE))
## Our chl a between 0 and 50m is 3.23 mg/m^3.

CCB_df %>%
  filter(bath_m <= -50, bath_m > -100) %>%
  summarize(chl = mean(chl_a, na.rm=TRUE))
## Our chl a between 50 and 100m is 1.51 mg/m^3.

## Now let's plot the raster data. 
CCB_map = ggplot() +
  geom_raster(data = CCB_df, aes(x=x, y=y, fill = bath_m)) +
  ## Adding coastline; group keeps multipolygon coordinates separated into distinct
  ## groups.
  geom_polygon(data= world_map, aes(x=long, y=lat, group=group), fill = "darkgrey",
               color=NA) +
  ## Croup map edges.
  coord_fixed(1.3, xlim = CCB_lon, ylim = CCB_lat, expand=FALSE) +
  ylab("Lat") + xlab("Lon") + theme_bw()
CCB_map
