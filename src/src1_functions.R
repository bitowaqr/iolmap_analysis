

# FUNCTIONS
  
# convenient function to install (if neccessary) and load required packages
  install_n_load = function(package){
      for(i in 1:length(package)){
        if(eval(parse(text=paste("require(",package[i],")")))==0) {
          install.packages(package)
        }
      }
      return (eval(parse(text=paste("require(",package,")"))))
      }
  
# compute basic descriptive stats
  summary.stats = function(x,round.dec=2){
    mean.d = mean(x)
    sd.d = sd(x)
    median.d = median(x)
    q25.d = quantile(x,c(.25))
    q75.d = quantile(x,c(.75))
    min.d = min(x)
    max.d = max(x)
    
    res = data.frame(mean.d,
               sd.d,
               median.d,
               q25.d,
               q75.d,
               min.d,
               max.d
    )
   round(res,round.dec) 
  }
  
# extract some basic glm model results
  extract_glm_values = function(model = m_dist_imd_runs){
    coefs = names(model$coefficients)
    betas = as.numeric(model$coefficients)
    ci.95 = confint.default(model)
    exp.ci.95 = formatC(exp(ci.95),digits=3,format="f")
    ci.95 = formatC(ci.95,digits=3,format="f")
    ci.95 = ifelse(substr(ci.95,1,1)=="-",ci.95,paste(" ",ci.95,sep=""))
    exp.betas = formatC(exp(betas),digits=3,format="f")
    betas = formatC(betas,digits=3,format="f")
    betas = ifelse(substr(betas,1,1)=="-",betas,paste(" ",betas,sep=""))
    ps = ifelse(coef(summary(model))[,'Pr(>|z|)']<0.001,"<0.001",paste(round(coef(summary(model))[,'Pr(>|z|)'],3)))
    names(ps) = NULL
    aic = model$aic
    bic = BIC(model)
    McFadden.pseudo.R = 1-(model$deviance/ model$null.deviance)
    cn = coefs
    c1 = paste(betas," (",ci.95[,1],"; ",ci.95[,2],")",sep="")
    c2 = paste(exp.betas," (",exp.ci.95[,1],"; ",exp.ci.95[,2],")",sep="")
    cn = c(cn,c("aic","bic","R2"))
    gofs = c(round(aic),round(bic),round(McFadden.pseudo.R,2))
    c1 = c(c1,gofs)
    c2= c(c2,c("","",""))
    ps = c(ps,c("","",""))
    tbl = cbind(cn,c1,c2,ps)
    colnames(tbl) = c("Variable","Beta coefficient (95% CI)","Rate ratio (95% CI)","p-value")
    return(tbl)
    
  }
  
######################
## parkrun event Location Allocation
##################
  

# Greedy adding location allocation algorithm
  LocAloAlg = function(candidates = greens_sp, # candidate green spaces for new events
                       event_coordinates.M = event_sp, # established parkrun events
                       lsoa_centr = cbind(lsoa_sp@data$centr_lng, # centroids of LSOAs
                                          lsoa_sp@data$centr_lat),
                       data = imd_runs_df,   # df including predictor variables for glm 
                       model = m_dist_imd_runs, # glm for estimating effect of new events
                       weights = lsoa_sp@data$imd_sc^2, # weights for estimated LSOA effects
                       top_n = 200,    # specifies how many locations are selected
                       objective_function = NULL, # you can specify alternative objective functions
                       maximise = T){ # maximise or minimise objective function?
                       
    
    require(geosphere)
    require(sp)
    
    # extract coordinates sub-function
    extract.coordinates = function(obj){
      if(class(obj)[1] =="SpatialPolygonsDataFrame"| class(obj)[1] =="SpatialPointsDataFrame"){
        res = coordinates(obj)  
      } else {
        if(dim(obj)[2] != 2){
          stop("obj must be either spatial or a matrix with 2 columns")
        }
        res = obj
      }
      return(res)
    }
    
    # objective function sub function default
    if(is.null(objective_function)){
        objective_function = function(m = model,    # takes a glm model 
                                 df = imd_runs_df, # a new df with predictor variables to estimate effect
                                 distances,        # distances (as a seperate predictor)
                                 w = 1){           # and weights (=1 for no weights)
          newdf = data.frame(df,mn_dstn = distances)
          res_temp = predict(m,newdata = newdf,type="response")
          res = res_temp * w
          return(res)
        }
    }
    
    # extract location coordinates
    lsoa.coord = extract.coordinates(lsoa_centr)
    events.coord = extract.coordinates(event_coordinates.M)
    candidates.coord = extract.coordinates(candidates)
    candidate_names = 1:length(candidates.coord[,1]) 
    
    cat("\n Computing distances between lsoa and parkrun events \n")
    # computing the distances between all LSOA and all established events
    # CAVE: high memory demands!
    dist_lsoa_events = distm(lsoa.coord,events.coord)
    # select the shortest distance for each LSOA
    min_dist_lsoa_events = apply(dist_lsoa_events,1,function(x) min(x))
    # evaluate the objective function at baseline
    trans_min_dist_lsoa_events = objective_function(m = model,df = data, distances = min_dist_lsoa_events/1000,w = weights)
    objective_state.append = sum(trans_min_dist_lsoa_events)  
    
    cat("\n Computing distances between lsoa and candidate parks \n")
    # computing the distance between all LSOA and all candidate parks
    # CAVE: high memory demands!
    dist_lsoa_candidates = distm(lsoa.coord,candidates.coord)
    # for each candidate, evaluate the objective function for all LSOAs
    trans_dist_lsoa_candidates = apply(dist_lsoa_candidates,2,function(x) objective_function(m = model,df = data, distances = x/1000,w = weights))
    
    
    cat("\n Running location allocation loop \n")
    # data frame to store results 
    save = data.frame(index=NA,objective=NA,change=rep(NA,times=top_n+1))
    # first row (step 0) shows baseline
    save[1,] = c(0,objective_state.append,0)
    
    
    
    # OUTER LOOP
    for(k in 1:top_n){ # repeat as many times as new locations are to be identified
      timecheck = 0 # progress index
      
      # how many candidate parks are (left) to be searched through
      len = length(candidate_names)
      # INNER LOOP
      for(i in 1:len){ # search through all candidate parks
        
        # show progress as % 
          progress = round((i/len),2)*100
          if(abs(progress - timecheck)>0.1){
            timecheck = progress
            cat("\r Outer loop: ",k, "   /",top_n,"; inner loop: %",progress, sep="")
            #     flush.console()
          }
        
        # evaluate for which LSOA candidate park i would be the (new) nearest event
        reassign_population = dist_lsoa_candidates[,i] < min_dist_lsoa_events
        if(sum(reassign_population)>0){ # if at least one LSOA improves...
          # evaluate the sum of the new objective function
          new_trans_min_dist_pop_combi = rep(NA,times=length(reassign_population))
          new_trans_min_dist_pop_combi[which(reassign_population)] = trans_dist_lsoa_candidates[reassign_population,i]
          new_trans_min_dist_pop_combi[-which(reassign_population)] = trans_min_dist_lsoa_events[!reassign_population]
          prop_objective_state = sum(new_trans_min_dist_pop_combi)
          
          # if the objective improves...
          if(maximise){
            improvement = prop_objective_state > objective_state.append
          } else {
            improvement = prop_objective_state < objective_state.append
          }
          
          
          if(improvement){
            # set candidate i as the optimal candidate for now
            objective_state.append = prop_objective_state
            index = i
            candidate_i = candidate_names[i]
            append.min_dist_pop_event = rep(NA,times=length(reassign_population))
            append.min_dist_pop_event[which(reassign_population)] = dist_lsoa_candidates[reassign_population,i]
            append.min_dist_pop_event[-which(reassign_population)] = min_dist_lsoa_events[!reassign_population]
            # update the objective function for the i+1th candidate park
            append.trans_min_dist_pop_event = new_trans_min_dist_pop_combi
            improvement = F
          }
        }
        
      }
      
      # the park that wins the inner loop is saved in the data frame at step k
      save[k+1,] = c(candidate_i,  # index of winner
                     objective_state.append, # new objective function
                     round((objective_state.append-save$objective[k])/save$objective[k],5)*100) # % change from previous objective function
      
      # add the winning candidate park as an established parkrun event for the next inner loop
      min_dist_lsoa_events = append.min_dist_pop_event
      trans_min_dist_lsoa_events = append.trans_min_dist_pop_event
      # remove winning candidate from the list of candidate parks
      dist_lsoa_candidates = dist_lsoa_candidates[,-index]
      trans_dist_lsoa_candidates = trans_dist_lsoa_candidates[,-index]
      candidate_names = candidate_names[-index]
    }
    
    # prepare results 
      top_consecutive_parks = candidates
      top_consecutive_parks$pos = match(1:length(top_consecutive_parks[,1]),save$index[-1])
      top_consecutive_parks_select = 1:length(top_consecutive_parks[,1]) %in% save$index
      top_consecutive_parks = subset(top_consecutive_parks,top_consecutive_parks_select)
      top_consecutive_parks = top_consecutive_parks[order(top_consecutive_parks$pos),]
      top_consecutive_parks$objective = save$objective[-1]/1e+9 # in millions
      top_consecutive_parks$change = save$change[-1]
      top_consecutive_parks.coord = coordinates(top_consecutive_parks)
      top_consecutive_parks$lon = top_consecutive_parks.coord[,1] 
      top_consecutive_parks$lat = top_consecutive_parks.coord[,2]
    
    return(top_consecutive_parks)
    
  }
  
  
  # Create  a static map from results
  create_static_map = function(LocAloAlg_res,
                               google_api_key = NULL,
                               old_events_sp = event_sp,
                               england_poly = NULL,
                               space.left = NULL,
                               space.right = NULL,
                               space.top = NULL,
                               space.bottom = NULL){
    require(ggplot2)
    require(ggrepel)
    require(ggmap)
    require(maptools)
    
    # set margins for map plot                                                            
    space.left = ifelse(is.null(space.left),0,space.left)
    space.right = ifelse(is.null(space.right),0,space.right)
    space.top = ifelse(is.null(space.top),0,space.top)
    space.bottom = ifelse(is.null(space.bottom),0,space.bottom)
    
    event_loc = as.data.frame(old_events_sp@coords)
    event_loc_new = data.frame(lon = LocAloAlg_res@data$lon,
                               lat = LocAloAlg_res@data$lat, 
                               lab = LocAloAlg_res@data$pos)
    
    lon1 = min(event_loc_new$lon) 
    lon2 = max(event_loc_new$lon)
    lat1 = min(event_loc_new$lat)
    lat2 = max(event_loc_new$lat)
    
    register_google(google_api_key) # you need a google maps api key
    
    # download google map tiles
    map.uk = get_googlemap(center = c(  -0.998597,52.701425),
                            zoom = 6,
                            maptype="terrain",
                            source="google",
                            color = "bw",region="uk",
                            style = 'feature:all|element:labels|visibility:off'
    )
    
    # create a polygon for England
    if(is.null(england_poly)){
      england_poly = maptools::unionSpatialPolygons(lsoa_sp, IDs = rep("England",times=dim(lsoa_sp)[1]))  
    }
    
    # build static map
    static_map = 
      # dark theme
      ggmap(map.uk,darken = 0) +
      # highlight England
      geom_polygon(data=england_poly,aes(x=long,y=lat,group=group),fill="white",col="white",alpha=0.1) +
      # coordinate grid lines
      geom_vline(xintercept = c(-5:1),col="white",size=0.5,alpha=0.5) +
      geom_hline(yintercept = c(50:56),col="white",size=0.5,alpha=0.5) +
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
      coord_map(xlim=c(lon1-space.left,lon2+space.right),
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
    
