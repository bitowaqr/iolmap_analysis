
######################
## parkrun event Location Allocation
## find optimal locations for future parkrun events
##################

    top_n  = 200   # how many new locations are needed?
  
  # setup
    events.coord = coordinates(event_sp)
    
    lsoa.coord = coordinates(lsoa_sp)
    lsoa.n = length(lsoa_sp$code)
    lsoa.pop = lsoa_sp$pop
    min_dist_lsoa_events = lsoa_sp$mn_dstn*1000
    objective_state.append = sum(min_dist_lsoa_events*lsoa.pop)  
    
    candidates.coord = coordinates(greens_sp)
    candidate_index = 1:length(greens_sp$id)
    
  # Compute all distances between LSOA and candidate green spaces 
  # ! ! ! cave: high memory demands ! ! !
    dist_lsoa_candidates = distm(lsoa.coord,candidates.coord)
  
  # data frame to store results 
    res_df = data.frame(index=NA,objective=NA,change=rep(NA,times=top_n+1),n=NA,pop=NA,mean=NA)
    # first row (step 0) shows baseline
    res_df[1,] = c(0,objective_state.append,0,0,0,0)
  
  t1 = Sys.time()
# GREEDY LOC ALO ALG
  # OUTER LOOP
  for(k in 1:top_n){ # repeat as many times as new locations are to be identified
    
    timecheck = 0 # inner loop progress index
    
    len = length(candidate_index) # n candidates remaining
    
    # INNER LOOP
    for(i in 1:len){ # search through all candidate parks
      
      # show progress as % 
      progress = round((i/len),2)*100
      if(abs(progress - timecheck)>0.1){
        timecheck = progress
        cat("\r Outer loop: ",k, "   /",top_n,"; inner loop: %",progress, sep="")
        #     flush.console()
      }
      
      # which lsoa had a lower mn_dstn if candidate i were a parkrun?
      reassign = dist_lsoa_candidates[,i] < min_dist_lsoa_events
      if(sum(reassign)>0){ # if at least one LSOA improves...
        # evaluate the sum of the new objective function
        new_trans_min_dist_pop_combi = rep(NA,times=lsoa.n)
        new_trans_min_dist_pop_combi[reassign] = dist_lsoa_candidates[reassign,i]
        new_trans_min_dist_pop_combi[!reassign] = min_dist_lsoa_events[!reassign]
        prop_objective_state = sum(new_trans_min_dist_pop_combi*lsoa.pop)
        
        if(prop_objective_state < objective_state.append){  # at least 1 lsoa improved?
          # set candidate i as the provisional new optimal candidate
          objective_state.append = prop_objective_state
          index = i
          candidate_i = candidate_index[i]
          append.min_dist_pop_event = rep(NA,times=lsoa.n)
          append.min_dist_pop_event[reassign] = dist_lsoa_candidates[reassign,i]
          append.min_dist_pop_event[!reassign] = min_dist_lsoa_events[!reassign]
          n.affected = sum(reassign)
          pop.affected = sum(lsoa.pop[reassign])
          diff.dist = dist_lsoa_candidates[reassign,i] - min_dist_lsoa_events[reassign]
          mean.dist.change = mean(diff.dist)
        }
      }
      
    } # inner loop end
    
    # the park that wins the inner loop is saved in the data frame at step k
    res_df[k+1,] = c(index = candidate_i,  
                     objective = objective_state.append, 
                     change = round((objective_state.append-res_df$objective[k])/res_df$objective[k],5)*100,
                     n = n.affected,
                     pop = pop.affected,
                     mean = mean.dist.change) # % change from previous objective function
    
    # add the winning candidate park as an established parkrun event for the next inner loop
    min_dist_lsoa_events = append.min_dist_pop_event
    # remove winning candidate from the list of candidate parks
    dist_lsoa_candidates = dist_lsoa_candidates[,-index]
    candidate_index = candidate_index[-index]
  
  } # outer loop end
  
  t2 = Sys.time()
  t2-t1
  # prepare results 
    top_consecutive_parks = greens_sp@data[res_df$index,]
    top_consecutive_parks$objective = res_df$objective[-1]/1e+9 # in millions
    top_consecutive_parks$change = res_df$change[-1]
    top_consecutive_parks$lon = candidates.coord[res_df$index,1]
    top_consecutive_parks$lat = candidates.coord[res_df$index,2]
    top_consecutive_parks$rank = 1:length(top_consecutive_parks$id)
  
    write.csv(top_consecutive_parks,"./output/optimal_locations.csv",row.names = F)


    
    
    
    
# / -------
# 4 CREATE STATIC MAPS FOR PUBLICATION -------
# the functions below use google maps API to retrieve the baseline tile
# to use this function, you need a api key
# For more details, see: https://developers.google.com/maps/documentation/maps-static/intro
    if(file.exists(".key.txt")){
      my_api_key = readChar(".key.txt",nchars = 100)
    } else {
      cat("Google API key required \n For more details, see: https://developers.google.com/maps/documentation/maps-static/intro")
      my_api_key = "ENTER_YOUR_API_KEY"  
    }
    
    pr_map_static = create_static_map(new_lon=top_consecutive_parks$lon,
                                      new_lat=top_consecutive_parks$lat,
                                      rank=top_consecutive_parks$rank,
                                      google_api_key = my_api_key,
                                      old_events_sp = event_sp,
                                      england_poly = england_sp,
                                      space.left = 0.3,
                                      space.right = 0.5,
                                      space.top = 0.3,
                                      space.bottom = 0.3)
    # pr_map_static
    ggsave(plot=pr_map_static,filename = "./output/map_static.jpeg",height = 10,width = 8)
    
    # + legend
    map_legend = legend_builder(top_consecutive_parks)
    write.csv(map_legend,file="./output/map_legend.csv",row.names = F)
    
    
  ## EFFECT OF 200 NEW EVENTS
    # combining new and old events
    after_events = rbind(coordinates(event_sp),
                         cbind(top_consecutive_parks$lon,
                               top_consecutive_parks$lat))
    
    # new distances
    dist_M_new = geosphere::distm(coordinates(lsoa_sp),
                                   coordinates(after_events))
    
  # NEW MIND DISTANCE to the narest new or old event for each LSOA
    lsoa_min_dist_new = apply(dist_M_new,1,FUN= function(x){round(min(x),0)} )
    lsoa_sp$mn_dstn_new = lsoa_min_dist_new / 1000
    lsoa_sp$mn_dstn_diff = lsoa_sp$mn_dstn - lsoa_sp$mn_dstn_new
    
    
    
    
    
    
    