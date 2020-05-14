

# FUNCTIONS
  
# convenient function to install (if neccessary) and load required packages
install_n_load = function(package){
  for(i in 1:length(package)){
    if(eval(parse(text=paste("require(",package[i],")")))==0) {
      install.packages(package[i])
    }
  }
  return (eval(parse(text=paste("require(",package[i],")"))))
}


# compute basic descriptive stats
  summary.stats = function(x,round.dec=2){
    
    res = data.frame(mean.d = mean(x),
                     sd.d = sd(x),
                     median.d = median(x),
                     q25.d = quantile(x,c(.25)),
                     q75.d = quantile(x,c(.75)),
                     min.d = min(x),
                     max.d = max(x)
    )
   round(res,round.dec) 
  }
  
# n % living within x km
  catch_pop = function(x,dist = lsoa_sp$mn_dstn ){
    paste(round(sum(lsoa_sp$pop[dist<=x])/sum(lsoa_sp$pop),4)*100,"/ ", formatC(sum(lsoa_sp$pop[dist<=x]),format="f",big.mark = ","))
    
  }
  
# Create  a static map from results
  create_static_map = function(new_lon,
                               new_lat,
                               rank ,
                               google_api_key = NULL,
                               old_events_sp = event_sp,
                               england_poly = NULL,
                               space.left = NULL,
                               space.right = NULL,
                               space.top = NULL,
                               space.bottom = NULL){
    
    ggmap_bbox <- function(map) {
      if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
      # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
      # and set the names to what sf::st_bbox expects:
      map_bbox <- setNames(unlist(attr(map, "bb")), 
                           c("ymin", "xmin", "ymax", "xmax"))
      
      # Coonvert the bbox to an sf polygon, transform it to 3857, 
      # and convert back to a bbox (convoluted, but it works)
      bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), "WGS84"))
      
      # Overwrite the bbox of the ggmap object with the transformed coordinates 
      attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
      attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
      attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
      attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
      map
    }
    
    
    # set margins for map plot                                                            
    space.left = ifelse(is.null(space.left),0,space.left)
    space.right = ifelse(is.null(space.right),0,space.right)
    space.top = ifelse(is.null(space.top),0,space.top)
    space.bottom = ifelse(is.null(space.bottom),0,space.bottom)
    
    event_loc = as.data.frame(old_events_sp@coords)
    event_loc_new = data.frame(lon = new_lon,
                               lat = new_lat,
                               lab = rank)
    
    lon1 = min(event_loc_new$lon) 
    lon2 = max(event_loc_new$lon)
    lat1 = min(event_loc_new$lat)
    lat2 = max(event_loc_new$lat)
    
    register_google(google_api_key) # you need a google maps api key
    
    # download google map tiles
    map.uk = get_googlemap(center = c(-0.998597,52.701425),
                            zoom = 6,
                            maptype="terrain",
                            source="google",
                            color = "bw",region="uk",
                            style = 'feature:all|element:labels|visibility:off'
    )
    
    england_poly <- st_transform(england_poly, "WGS84")
    
    
    
    map <- ggmap_bbox(map.uk)
    
    static_map = 
      ggmap(map,darken = 0) +
        coord_sf(crs = st_crs("WGS84")) + # force the ggplot2 map to be in 3857
        geom_sf(data = england_poly, aes(), inherit.aes = FALSE,fill="white",col="white",alpha=0.3) +
      # coordinate grid lines
      geom_vline(xintercept = c(-5:1),col="white",size=0.25,alpha=0.5) +
      geom_hline(yintercept = c(50:56),col="white",size=0.25,alpha=0.5) +
      # old event locations
      geom_point(data=event_loc,aes(x=lng,y=lat),      col="darkblue",alpha=1,size=0.9,shape=1) +
      # proposed new event locations
      geom_point(data=event_loc_new,aes(x=lon,y=lat),col="darkred",alpha=1,size=1.2,shape = 6) +
      # label new events by rank
      geom_text_repel(data=event_loc_new,aes(x=lon,y=lat,label=lab),
                      size=2.1,col="darkred",
                      box.padding = unit(0.05, "lines"),
                      segment.color = 'black',
                      segment.alpha = 1,
                      segment.size = 0.2) +
      # set map limits
      coord_sf(xlim=c(lon1-space.left,lon2+space.right),
                ylim=c( lat1-space.bottom,lat2+space.top)) +
      # remove labels
      xlab("") +
      ylab("") 
    
    return(static_map)
  }
  
  
  # Create a legend for the static map
  legend_builder = function(LocAloAlg_res, digits = 3){
      # builds tables for publication with lng and lat coordinates of new event locations
      name.park = substr(LocAloAlg_res$distName1,1,25)
      name.park = ifelse(nchar(LocAloAlg_res$distName1)>20,paste(substr(name.park,1,20),"...",sep=""),name.park)
      name.park[is.na(name.park)] = "unnamed"
      legend_m1 = cbind("Pos"=LocAloAlg_res$pos,
                        "Park"=name.park,
                        "Change (%)"=LocAloAlg_res$change,
                        "Area (km2)"=formatC(LocAloAlg_res$area_km2,digits = 2, format = "f"),
                        "Longitude"=round(LocAloAlg_res$lon,digits),
                        "Lattitude"=round(LocAloAlg_res$lat,digits))
      return(legend_m1)
    }
    
