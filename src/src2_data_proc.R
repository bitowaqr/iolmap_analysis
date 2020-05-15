
# SRC 2: Load and process data

# Load pkgs -------
    install_n_load(c("raster",
                     "sp",
                     "sf",
                     "mapproj",
                     "scales",
                     "geosphere",
                     "rgeos",
                     "reshape2",
                     "dplyr",
                     "ggplot2",
                     "eq5d",
                     "cowplot",
                     "ggrepel",
                     "lwgeom",
                     "maptools",
                     "ggmap",
                     "wCorr")
                   )
  
# Load raw data -------
    
    # LSOA population-weighted centroids
      # Contains public sector information licensed under the Open Government Licence v3.
      # Office for National Statistics (2011). 2011 Census: boundary data (England and Wales) [data collection].
      # Retrieved from: https://data.gov.uk/dataset/a40f54f7-b123-4185-952f-da90c56b0564/lower-layer-super-output-areas-december-2011-population-weighted-centroids
      lsoa_sp = raster::shapefile("./input/England_lsoa_2011_centroids/england_lsoa_2011_centroids")
      
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
    
    # LSOA polygons/shapefiles   
      lsoa_pg = st_read("./input/England_lsoa_2011_sgen_clipped/england_lsoa_2011_sgen_clipped.shp")
      lsoa_pg = st_transform(lsoa_pg,CRS("+proj=longlat"))
      lsoa_pg = st_make_valid(lsoa_pg)
      england_sp = sf::st_union(lsoa_pg,by_feature=F)
      
    # GREEN SPACES 
      # Ordnance Survey. OS Open Greenspace. 2018. 
      # Licensed under the Open Government Licence v3.
      # Retrieved from: https://www.ordnancesurvey.co.uk/opendatadownload/products.html#OPGRSP
      # due to Github file size restrictions, data had to be trimmed - see below: 
      
    # parkrun event data 
      # Parkrun UK. 2019. 
      # Retrieved from: http://www.parkrun.org.uk/
      event_sp = read.csv("./input/parkrun_data/event_info_20181212.csv", stringsAsFactors = F)
      
    
# Process and combine data --------
  
  # LSOA DATA 
      # location
      lsoa_sp = spTransform(lsoa_sp,CRS("+proj=longlat"))
      
      # add population + pop.density
      lsoa_pop = lsoa_pop %>% 
        rename("code" = 1,
               "pop" = 3) %>%
        select(1,3) %>%
        mutate(pop = gsub(",","",pop)) %>%
        mutate(pop = as.numeric(pop))
      lsoa_sp = merge(lsoa_sp,lsoa_pop,by=c("code"),sort=F)
      
      # prepare imd scores 
      lsoa_imd = lsoa_imd %>% 
        dplyr::rename(code = LSOA.code..2011.,
                      imd_sc = 'Index.of.Multiple.Deprivation..IMD..Score',
                      imd_ra = 'Index.of.Multiple.Deprivation..IMD..Rank..where.1.is.most.deprived.'
                      ) %>%
        dplyr::select(c("code","imd_sc","imd_ra")) %>%
        mutate(imd_q5 = cut(imd_sc,quantile(imd_sc,seq(0,1,0.2)),
                            labels=c("Least deprived","Less deprived","Median deprived","More deprived","Most deprived")))
      
      lsoa_sp = merge(lsoa_sp,lsoa_imd,by=c("code"),sort=F)
      
  # PARKRUN EVENTS
      coordinates(event_sp) = ~lng+lat
      projection(event_sp) = "+proj=longlat +ellps=WGS84"
      
  # DISTANCE DATA
      # COMPUTE distances matrix (distances between all LSOA and all parkrun events)
      dist_M_full = geosphere::distm(coordinates(lsoa_sp),
                                     coordinates(event_sp))
      
    # distance to the narest event for each LSOA
      lsoa_min_dist = apply(dist_M_full,1,FUN= function(x){round(min(x),0)} )
      lsoa_sp$mn_dstn = lsoa_min_dist / 1000
      
    # determine name of nearest parkrun event
      lsoa_min_dist_event = apply(dist_M_full,1,FUN= function(x){which(x == min(x))} )
      lsoa_min_dist_event = data.frame(code = lsoa_sp$code,nrst_evnt =  event_sp$course[lsoa_min_dist_event])
      lsoa_sp = merge(lsoa_sp,lsoa_min_dist_event,by="code")
      
      
      # parkrun event stats: how large is the catchment area?
      srvd_lsoa = apply(dist_M_full,1,FUN= function(x){which(x == min(x))} )
      srvd_pop = data.frame(code = lsoa_sp$code,course = event_sp$course[as.numeric(srvd_lsoa)])
      srvd_pop$srvd_pop = lsoa_sp$pop[match(srvd_pop$code,lsoa_sp$code)]
      srvd_pop = aggregate(srvd_pop ~course ,srvd_pop,sum)
      srvd_lsoa = table(srvd_lsoa)
      srvd_lsoa = data.frame(course = event_sp$course,srvd_lsoa = as.numeric(srvd_lsoa))
      event_sp = merge(event_sp,srvd_lsoa,by="course")
      event_sp = merge(event_sp,srvd_pop,by="course")
      
      
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
  
# CLEAN UP
    rm(list=c("lsoa_min_dist_event","srvd_lsoa","srvd_pop","MIN_PARK_AREA","greens_in_england","lsoa_imd","lsoa_pop","lsoa_cntrds"))

# SAVE the created data sets (saves time for later use)
    save(list=c("dist_M_full","event_sp","greens_sp","lsoa_min_dist","lsoa_sp"),
         file = paste("./output/savegame_", Sys.Date(),".Rdata",sep=""))

    