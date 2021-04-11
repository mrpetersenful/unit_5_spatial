## 31 March 2021
## Lesson 5.3: AIS ship tracks

## Automatic Identification System (AIS) is an automatic tracking system that uses
## tranceivers on vessels to detect their location in space and time. AIS is used 
## secondarily to radar to avoid vessel collisions. Large vessels and passenger 
## ships are required to use AIS. The US Bureau of Ocean Energy Management (BOEM)
## and NOAA provide AIS data collected by the US Coast Guard for free at 
## https://marinecadastre.gov/AIS/

## AIS data is huge! There are a lot of ships in the water and each ping from an 
## AIS unit is a new row of data. I downloaded AIS data from just the month of 
## January 2017 in UTM Zone 17 (i.e., the US East Coast). This is a 4GB spreadsheet
## with 31,884,021 rows. This is too big to include in a GitHub repo or to
## effectively use in class. In the script ais_reduce.R I read in the data, cropped 
## it to a bounding box around Florida and Georgia, removed duplicate rows and rows
## with no vessel CallSign and grabbed only vessels with the status "under way using
## engine". Then I wrote this much smaller file (16 MB) out to a .csv file that has 
## a more manageable size, but is otherwise equivalent to what you could download
## directly from the website. 

## We can use the AIS data to get a look at the amount of vessel traffic that occurs
## in the Southeast US North Atlantic right whale critical breeding habitat. Much of
## this critical habitat has been designated a Seasonal Management Area, where 
## vessels > 65 feet long are required to slow down to <10 knots during the calving
## season (from November 15-April 15 each year). Vessels < 65 feet in length are 
## also encouraged, but not required, to slow down to <10 knots. NOAA issued this
## ruling to reduce the risk of right whale ship strikes. 

## Let's start by loading in the subsetted AIS data, the shapefiles for the right
## whale critical habitat (that we used earlier) and the world coastline data in R's
## mapdata package. We will crop the coastline data to the same bounding box that
## I used when I reduced the size of the original (huge) AIS data.

library(tidyverse)
library(sf)
library(mapdata)
library(marmap)
library(lubridate)

## AIS data; Downloaded for January 2017, UTM Zone 17
## https://marinecadastre.gov/AIS/ (not included in data folder bc ~1GB)
## subsetted to 2017-1-25 data with script ais_reduce.R in this repo
ais_day = read.csv('data/processed_ais/ais_2017-01-25.csv')

## Now add some coastline data
lat_bounds = c(25, 34)
lon_bounds = c( -82, -76)
world_map = map_data("worldHires", ylim = lat_bounds, xlim = lon_bounds)
dim(world_map)

## Now let's read in US critical habitat shapefiles.
# https://www.greateratlantic.fisheries.noaa.gov/educational_resources/gis/data/index.html
USA_crit_hab = st_read('data/North_Atlantic_Right_Whale_Critical_Habitat','North_Atlantic_Right_Whale_Critical_Habitat')
USA_crit_hab_sf = st_transform(USA_crit_hab, crs=4326) #crs="+proj=longlat +datum=WGS84")

## Now let's plot everything that we loaded up. Note that the AIS data is still 
## quite large, with 119,726 rows. Each row represents a lat/lon coordinate where
## an AIS signal was emitted. Putting all of these points on a map takes some 
## processing time, but you can significantly speed it up by never actually printing
## the plot out to the console. Do this by saving the plot to a variable name and
## using ggsave() to write that plot out to a file (like a .pdf, .png or .jpg). Then
## just open up that file to see your map.


head(ais_day)

## Now to plot critical habitats and carcass locations.
ais_map_pts = ggplot()+
  ## Adding coastlines.
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "black", color = NA) + 
  geom_sf(data=USA_crit_hab_sf, alpha = 0.5, color=NA, fill='yellow') +
  geom_point(data = ais_day, aes(x = LON, y = LAT, color = CallSign), size=2) + 
  ## Cropping map edges.
  coord_sf(1.3, xlim = lon_bounds, ylim = lat_bounds) + 
  guides(color=FALSE) +
  ylab("Latitude") + xlab("Longitude") + theme_classic() 

ggsave(ais_map_pts, filename='figures/ais_pts_2017-01-25.pdf', device="pdf", 
       height=5, width=4)


## Now we'll do some analysis. First we can use the ymd_hms() function in the 
## lubridate package to convert the column BaseDateTime into a date/time variable
## (of type POSIXct) that R will actually recognize. Then we can simply use the 
## arrange() function in dplyr to sort our data by date/time to make sure they are
## in chronological order. Once that is done, we can group the AIS points by each
## vessel's unique CallSign and collapse the series of AIS points into ship track
## lines. Note that arranging the data chronologically and grouping by vessel are
## critical steps so that the ship tracks make sense.

dim(ais_day)

ais_day_sf = ais_day %>% 
  mutate(date_time = lubridate::ymd_hms(BaseDateTime)) %>%
  ## This ensures that ship track points are in chronological order.
  arrange(date_time) %>% 
  ## This puts our data into a shapefile with CRS 4326.
  st_as_sf(coords=c("LON", "LAT"), crs=4326) %>% 
  ## This essentially groups by call sign but retains length (in feet) info
  group_by(CallSign, Length) %>%
  ## This collapses the data into multipoint; do_union=FALSE prevents reordering 
  ## the points.
  summarise(do_union=FALSE) %>% 
  st_cast("LINESTRING") 

class(ais_day_sf$BaseDateTime)
class(ais_day_sf$date_time) # POSIXct is the date-time class. Well, actually it's not.

class(ais_day_sf)
glimpse(ais_day_sf)
dim(ais_day_sf)  
## Nice. We went from 119,726 points to 283 lines.

## Now that we've collapsed the AIS points into sensible ship tracks, let's plot the
## results.

ais_day_line_map = ggplot()+
  ## Adding coastlines.
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "black", color = NA) + # add coastline
  geom_sf(data=USA_crit_hab_sf, alpha = 0.5, color=NA, fill='yellow') +
  geom_sf(data=ais_day_sf, aes(color=CallSign)) +
  ## Cropping map edges.
  coord_sf(1.3, xlim = lon_bounds, ylim = lat_bounds) + 
  guides(color=FALSE) +
  ylab("Latitude") + xlab("Longitude") + theme_classic() 

ggsave(ais_day_line_map, filename='figures/ais_lines_2017-01-25.pdf', 
       device="pdf", height=5, width=4)

## Now we can see all of the ships that entered this right whale critical habitat 
## during the heart of calving season (January). We can use sf + dplyr packages to 
## ask really simple questions like: 
#### How many ships entered the critical habitat on 2017-01-25?
#### What are the lengths (in feet) of ships that entered the habitat?
#### What is the total length (in meters) of ship tracks within the critical habitat
####   on 2017-01-25?


#########################################################################
#   Compare 2017-01-25 ship tracks with NARW critical habitats
#########################################################################

library(lwgeom) # where st_make_valid() was originally developed

## Let's find ship tracks that enter RW habitat with spatial join.
ships_RW_join = ais_day_sf %>%
  st_make_valid() %>%
  st_join(USA_crit_hab_sf %>% 
            ## This is a spatial inner join -- only keeps lines that intersect with
            ## polygon.
            dplyr::select(geometry), left=FALSE) 

?st_make_valid

## How many ships went into the RW critical habitat on 2017-1-25?
dim(ships_RW_join)

# What are the lengths of the ships that intersected the RW critical habitat?
ggplot() +
  geom_histogram(data=ships_RW_join, aes(x=Length))

# Plot ship tracks that enter RW habitat
ggplot()+
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "black", color = NA) + # add coastline
  geom_sf(data=USA_crit_hab_sf, alpha = 0.5, color=NA, fill='yellow') +
  geom_sf(data=ships_RW_join, aes(color=CallSign)) +
  coord_sf(1.3, xlim = lon_bounds, ylim = lat_bounds) + # Crop map edges
  guides(color=FALSE) +
  ylab("Latitude") + xlab("Longitude") + theme_classic() 

# Now just grab portions of ship tracks that intersect with RW habitat
ships_RW_intersect = ais_day_sf %>%
  st_make_valid() %>%
  st_intersection(USA_crit_hab_sf %>% select(geometry)) %>%
  mutate(track_length_m = as.numeric(st_length(geometry))) # Calculate length of ship tracks in meters

# Plot ship tracks that intersect with RW habitat
ggplot()+
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "black", color = NA) + # add coastline
  geom_sf(data=USA_crit_hab_sf, alpha = 0.5, color=NA, fill='yellow') +
  geom_sf(data=ships_RW_intersect, aes(color=CallSign)) +
  coord_sf(1.3, xlim = lon_bounds, ylim = lat_bounds) + # Crop map edges
  guides(color=FALSE) +
  ylab("Latitude") + xlab("Longitude") + theme_classic() 

# What are the lengths of the ship tracks that intersected the RW critical habitat?
ships_RW_intersect %>%
  summarize(tot_track_length_m = sum(track_length_m)) %>%
  mutate(tot_track_length_km = tot_track_length_m/1000) %>%
  st_drop_geometry()

