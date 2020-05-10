
# LOAD, MERGE, AND COMPUTE THE REQUIRED DATA

# LOAD FUNCTIONS + PACKAGES
  # Load functions
  #  source("./src/src1_functions.R") # required
  
  # Install and load all required packages
    install_n_load(c("raster",
                     "sp",
                     "geosphere",
                     "rgeos",
                     "dplyr",
                     "ggplot2",
                     "ggrepel",
                     "maptools",
                     "ggmap")
                   )
  
# LOAD RAW DATA
    # LSOA polygon geospatial data 
      # Contains public sector information licensed under the Open Government Licence v3.
      # Office for National Statistics (2011). 2011 Census: boundary data (England and Wales) [data collection].
      # Retrieved from: https://data.gov.uk/dataset/fa883558-22fb-4a1a-8529-cffdee47d500/lower-layer-super-output-area-lsoa-boundaries
      lsoa_sp = raster::shapefile("./input/England_lsoa_2011_sgen_clipped/england_lsoa_2011_sgen_clipped")
    
    # LSOA population-weighted centroids
      # Contains public sector information licensed under the Open Government Licence v3.
      # Office for National Statistics (2011). 2011 Census: boundary data (England and Wales) [data collection].
      # Retrieved from: https://data.gov.uk/dataset/a40f54f7-b123-4185-952f-da90c56b0564/lower-layer-super-output-areas-december-2011-population-weighted-centroids
      lsoa_cntrds = raster::shapefile("./input/England_lsoa_2011_centroids/england_lsoa_2011_centroids")
      
    # LSOA population + LSOA area sizes
      # Contains public sector information licensed under the Open Government Licence v3.
      # Office for National Statistics (2017). Lower layer Super Output Area population density (National Statistics). Mid-2017: SAPE20DT11: 
      # Retrieved from: https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareapopulationdensity
      lsoa_pop = read.csv("./input/Population_estimates/Mid-2017 Population Density-Table 1.csv",stringsAsFactors = F)
    
    # LSOA IMD scores
      # Contains public sector information licensed under the Open Government Licence v3.
      # Office for National Statistics (2015). English indices of deprivation 2015.  Official Statistics. 2015.
      # Retrieved from: https://www.gov.uk/government/statistics/english-indices-of-deprivation-2015
      lsoa_imd = read.csv("./input/IMD_data/IMD_data.csv",stringsAsFactors = F)
    
    # GREEN SPACES 
      # Ordnance Survey. OS Open Greenspace. 2018. 
      # Licensed under the Open Government Licence v3.
      # Retrieved from: https://www.ordnancesurvey.co.uk/opendatadownload/products.html#OPGRSP
      # due to Github file size restrictions, data had to be trimmed - see below: 
      
    # parkrun event data 
      # Parkrun UK. 2019. 
      # Retrieved from: http://www.parkrun.org.uk/
      event_sp = read.csv("./input/parkrun_data/events_corrected.csv", stringsAsFactors = F)
    
    # LSOA parkrun participation data (raw data from parkrunUK, filtered and processed)
      # Parkrun UK. 2018. 
      runs_per_lsoa = read.csv("./input/parkrun_data/runs_per_lsoa_010118_101218.csv", stringsAsFactors = F)
    

      
# PROCESS DATA
  
  # LSOA DATA 
    # init lsoa_sp: correct projection
      lsoa_sp = spTransform(lsoa_sp,CRS("+proj=longlat"))
    # add population weighted centroids to lsoa_sp
      lsoa_cntrds = spTransform(lsoa_cntrds,CRS("+proj=longlat"))
      lsoa_cntrds_coord = coordinates(lsoa_cntrds)
      lsoa_cntrds_coord_df = data.frame(code = lsoa_cntrds@data$code,
                                        centr_lng = lsoa_cntrds_coord[,1],
                                        centr_lat = lsoa_cntrds_coord[,2],
                                        stringsAsFactors = F
                                        )
      lsoa_sp = merge(lsoa_sp,lsoa_cntrds_coord_df,by=c("code"),sort=F)
      
    # add population + density figures to lsoa_sp
      lsoa_pop = lsoa_pop[,c(1,3,4,5)]
      names(lsoa_pop) = c("code","pop","area_km2","pop_km2")
      lsoa_pop$pop = as.numeric(gsub(",","",lsoa_pop$pop))
      lsoa_pop$pop_km2 = as.numeric(gsub(",","",lsoa_pop$pop_km2))
      lsoa_sp = merge(lsoa_sp,lsoa_pop,by=c("code"),sort=F)
    
    # add time var (observational period = 49 weeks)
      lsoa_sp@data$week = 49
    
    # prepare imd scores 
      lsoa_imd = lsoa_imd %>% 
        dplyr::rename(code = LSOA.code..2011.) %>% 
        dplyr::rename(imd_sc = 'Index.of.Multiple.Deprivation..IMD..Score') %>%
        dplyr::select(c("code","imd_sc"))
      lsoa_sp = merge(lsoa_sp,lsoa_imd,by=c("code"),sort=F)
  
  # PARKRUN EVENTS
      coordinates(event_sp) = ~lng+lat
      projection(event_sp) = "+proj=longlat +ellps=WGS84"
      
  # GREENSPACES 
      # trimming of the original data set (data cannot be provided via this repository due to file size restrictions)
      # greens_sp = raster::shapefile("./input/greenspaces_data/GB_GreenspaceSite")
      # # filter greenspaces
      # greens_sp = spTransform(greens_sp,CRS("+proj=longlat"))
      # # 1. Filter parks on purpose
      # greens_sp = subset(greens_sp,greens_sp@data$function. %in% c("Playing Field","Play Space","Other Sports Facility","Public Park Or Garden"))
      # # 2. Filter parks on min park area
      # MIN_PARK_AREA = 0.1 # km2
      # greens_sp$area_km2 = raster::area(greens_sp)/1000^2 # in square km
      # greens_sp = subset(greens_sp,greens_sp$area_km2 >= MIN_PARK_AREA)
      # # 2. select subset of parks: only English
      # # england_sp = maptools::unionSpatialPolygons(lsoa_sp, IDs = rep("England",times=dim(lsoa_sp)[1]))
      # greens_in_england = rgeos::gIntersects(england_sp,greens_sp, byid = T)
      # greens_sp = subset(greens_sp,greens_in_england[,1])
      # shapefile(greens_sp,filename="./input/greenspaces_data/trimmed_greenspaces.shp") # save trimmed data set
      greens_sp = raster::shapefile("./input/greenspaces_data/trimmed_greenspaces")
  
  # DISTANCE DATA
    # COMPUTE distances matrix (distances between all LSOA and all parkrun events)
      dist_M_full = geosphere::distm(cbind(lsoa_sp$centr_lng,lsoa_sp$centr_lat),coordinates(event_sp))
      rownames(dist_M_full) = lsoa_sp$code
      colnames(dist_M_full) = event_sp$course 
    # distance to the narest event for each LSOA
      lsoa_min_dist = apply(dist_M_full,1,FUN= function(x){round(min(x),0)} )
      lsoa_sp$mn_dstn = lsoa_min_dist / 1000
    # determine name of nearest parkrun event
      lsoa_min_dist_event = apply(dist_M_full,1,FUN= function(x){which(x == min(x))} )
      lsoa_min_dist_event = data.frame(code = row.names(dist_M_full),nrst_evnt =  colnames(dist_M_full)[lsoa_min_dist_event])
      lsoa_sp = merge(lsoa_sp,lsoa_min_dist_event,by="code")
    
    # For the parkrun events: how large is the catchment area?
      # i.e. for how many LSOAs/people is the event the nearest?
      srvd_lsoa = apply(dist_M_full,1,FUN= function(x){which(x == min(x))} )
      srvd_pop = data.frame(code = names(srvd_lsoa),course = colnames(dist_M_full)[as.numeric(srvd_lsoa)])
      srvd_pop$srvd_pop = lsoa_sp@data$pop[match(srvd_pop$code,lsoa_sp@data$code)]
      srvd_pop = aggregate(srvd_pop ~course ,srvd_pop,sum)
      srvd_lsoa = table(srvd_lsoa)
      srvd_lsoa = data.frame(course = colnames(dist_M_full),srvd_lsoa = as.numeric(srvd_lsoa))
      event_sp = merge(event_sp,srvd_lsoa,by="course")
      event_sp = merge(event_sp,srvd_pop,by="course")
    
  
  # PARTICIPATION data
    # prepare runs per week, per 1,000 population 
      lsoa_sp = merge(lsoa_sp,runs_per_lsoa,by="code",sort=F)
      lsoa_sp$run_count[is.na(lsoa_sp$run_count)] = 0
      lsoa_sp$runs_pmil.week = ((lsoa_sp$run_count/49)/lsoa_sp$pop)*1000 # 49 weeks observational period
      
# CLEAN UP
    rm(list=c("lsoa_min_dist_event","srvd_lsoa","srvd_pop","runs_per_lsoa","MIN_PARK_AREA","greens_in_england","lsoa_imd","lsoa_cntrds_coord_df","lsoa_pop","lsoa_cntrds","lsoa_cntrds_coord","events_in_england","prisonruns","inPrisons"))
  
# SAVE the created data sets (saves time for later use)
    save(list=c("dist_M_full","event_sp","greens_sp","lsoa_min_dist","lsoa_sp"),
         file = paste("./output/savegame_", Sys.Date(),".Rdata",sep=""))

    