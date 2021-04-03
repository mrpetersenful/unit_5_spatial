## 31 March 2021
## Lesson 5.3: AIS



## ---- message=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(mapdata)
library(marmap)
library(lubridate)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## AIS data; Downloaded for January 2017, UTM Zone 17
# https://marinecadastre.gov/AIS/ (not included in data folder bc ~1GB)
# subsetted to 2017-1-25 data with script ais_reduce.R in this repo
ais_day = read.csv('data/processed_ais/ais_2017-01-25.csv')

# Coastline data
lat_bounds = c(25, 34)
lon_bounds = c( -82, -76)
world_map = map_data("worldHires", ylim = lat_bounds, xlim = lon_bounds)
dim(world_map)

#Read in US critical habitat shapefiles 
# https://www.greateratlantic.fisheries.noaa.gov/educational_resources/gis/data/index.html
USA_crit_hab = st_read('data/North_Atlantic_Right_Whale_Critical_Habitat/','North_Atlantic_Right_Whale_Critical_Habitat') # reads in set of shapefiles
USA_crit_hab_sf = st_transform(USA_crit_hab, crs=4326) #crs="+proj=longlat +datum=WGS84")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(ais_day)

# plot critical habitats and carcass locations
ais_map_pts = ggplot()+
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "black", color = NA) + # add coastline
  geom_sf(data=USA_crit_hab_sf, alpha = 0.5, color=NA, fill='yellow') +
  geom_point(data = ais_day, aes(x = LON, y = LAT, color = CallSign), size=2) + 
  coord_sf(1.3, xlim = lon_bounds, ylim = lat_bounds) + # Crop map edges
  guides(color=FALSE) +
  ylab("Latitude") + xlab("Longitude") + theme_classic() 

ggsave(ais_map_pts, filename='figures/ais_pts_2017-01-25.pdf', device="pdf", height=5, width=4)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dim(ais_day)

ais_day_sf = ais_day %>% 
  mutate(date_time = lubridate::ymd_hms(BaseDateTime)) %>%
  arrange(date_time) %>% # ensure ship tracks points are in chronological order
  st_as_sf(coords=c("LON", "LAT"), crs=4326) %>% # '+proj=longlat +datum=WGS84'
  group_by(CallSign, Length) %>% # essentially groups by call sign but retains length (in feet) info
  summarise(do_union=FALSE) %>% # collapses data into multipoint; do_union=FALSE prevents reordering points
  st_cast("LINESTRING") 

# class(ais_day_sf$BaseDateTime)
# class(ais_day_sf$date_time) # POSIXct is the date-time class

class(ais_day_sf)
glimpse(ais_day_sf)
dim(ais_day_sf)  # went from 119,726 points to 283 lines !!


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plot ship tracks
ais_day_line_map = ggplot()+
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "black", color = NA) + # add coastline
  geom_sf(data=USA_crit_hab_sf, alpha = 0.5, color=NA, fill='yellow') +
  geom_sf(data=ais_day_sf, aes(color=CallSign)) +
  coord_sf(1.3, xlim = lon_bounds, ylim = lat_bounds) + # Crop map edges
  guides(color=FALSE) +
  ylab("Latitude") + xlab("Longitude") + theme_classic() 

ggsave(ais_day_line_map, filename='figures/ais_lines_2017-01-25.pdf', device="pdf", height=5, width=4)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#########################################################################
#   Compare 2017-01-25 ship tracks with NARW critical habitats
#########################################################################

# library(lwgeom) # where st_make_valid() was originally developed
# Find ship tracks that enter RW habitat with spatial join
ships_RW_join = ais_day_sf %>%
  st_make_valid() %>%
  st_join(USA_crit_hab_sf %>% select(geometry), left=FALSE) # spatial inner join (only keep lines that intersect with polygon)

# How many ships went into the RW critical habitat on 2017-1-25?
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

