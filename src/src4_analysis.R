
# ANALYSIS src
# Code to produce the results reported in the publication

# 1 SOME DESCRITPIVE STATS ---------------
     
    length(lsoa_sp$code) # n LSOAs
    summary.stats(lsoa_sp$pop) # population
    summary.stats(lsoa_sp$mn_dstn) # distances to the nearest event
    summary.stats(lsoa_sp$run_count) # run counts raw
    summary.stats(lsoa_sp$runs_pmil.week) # run rates
    summary.stats(lsoa_sp$imd_sc) # deprivation indices
    length(event_sp$course) # n parkrun events
    summary.stats(event_sp$srvd_pop) # served population
    summary.stats(event_sp$srvd_lsoa) # served LSOA
    summary.stats( event_sp$Mn_prtc) # mean participants per event
    summary.stats( event_sp$Mn_vlnt) # mean volunteers
    summary.stats(event_sp$Age_yrs) # age of events 
    
    # pop living_within_1km 
    round(sum(lsoa_sp$pop[lsoa_sp$mn_dstn<=1])/sum(lsoa_sp$pop),4)*100
    # pop living_within_5km
    round(sum(lsoa_sp$pop[lsoa_sp$mn_dstn<=5])/sum(lsoa_sp$pop),4)*100
    # pop living_within_10km
    round(sum(lsoa_sp$pop[lsoa_sp$mn_dstn<=10])/sum(lsoa_sp$pop),4)*100
    # pop living_outside_25km
    round(sum(lsoa_sp$pop[lsoa_sp$mn_dstn>25])/sum(lsoa_sp$pop),4)*100

# 2 POISSON REGRESSION MODELS --------
  
# 2.1 RUNS ~ DISTANCE ------------
    m_dist_runs = glm(run_count ~ mn_dstn,
                    data = lsoa_sp@data,
                    family = poisson(link="log"),
                    offset = log(pop*week))
    
    extract_glm_values(m_dist_runs)
    
# 2.2 RUNS ~ IMD  ------------
    m_imd_runs = glm(run_count ~ imd_sc,
                      data = lsoa_sp@data,
                      family = poisson(link="log"),
                      offset = log(pop*week))
    
    extract_glm_values(m_imd_runs)

# 2.3 RUNS ~ IMD + DISTANCE ------------
  m_dist_imd_runs = glm(run_count ~ imd_sc + mn_dstn,
                   data = lsoa_sp@data,
                   family = poisson(link="log"),
                   offset = log(pop*week))
    
    extract_glm_values(m_dist_imd_runs)
  
# 2.4 Illustrating relationships in a plot  ----
  imd_q_10_50_90 = quantile(lsoa_sp@data$imd_sc,c(0.1,0.5,0.9))
  newdf = data.frame(week = 1,
                           pop=1000,
                           imd_sc=rep(imd_q_10_50_90,each=100),
                           mn_dstn = rep(seq(0,35,length.out = 100),3),
                           ## LEGEND FOR PLOT
                           "IMD qunatile" = rep(c("10th percentile","50th percentile (median)","90th percentile"),each=100))
  pred_imd_dist_runs = predict(m_dist_imd_runs,newdata = newdf,type="response")
  
  p_imd_dist_runs = 
    ggplot() +
    geom_point(data=lsoa_sp@data,
               aes(x=mn_dstn,y=1000*(run_count/pop)/week),col="darkgoldenrod2",alpha=.3,size=.5 ) +
    geom_line(data=newdf,aes(x=mn_dstn,
                                   y=pred_imd_dist_runs,
                                   col=IMD.qunatile,
                                   linetype=IMD.qunatile),
              alpha=1,size=1.5,lineend ="round" ) +
    theme_minimal() +
    # 385 (1.2%) observations are not shown (outliers)
     coord_cartesian(xlim=c(0,25),ylim=c(0,7.5))+
    # axes label
    xlab("Distance to the nearest parkrun event in km") +
    ylab("Runs per 1,000 population per week")  +
    # legend title
    labs(col = "LSOA IMD") +
    labs(linetype="LSOA IMD") +
    theme(text=element_text( family="Arial")) +
    # legend position
    theme(legend.position = "bottom") 
  # save plot
  ggsave(plot=p_imd_dist_runs,filename = "./output/p_imd_dist_runs.jpeg",height = 6,width = 10)
  
  
  # / -------
# 3 THE EFFECT OF STARTING 200 NEW EVENTS ------------
  all.events = rbind(coordinates(event_sp),coordinates( opt_loc_runs_imd.prep))
  new.events.names = ifelse(is.na(opt_loc_runs_imd.prep$distName1),"No name",opt_loc_runs_imd.prep$distName1)
  
  # COMPUTE NEW DISTANCE MATRIX 
  new.dist_M = geosphere::distm(cbind(lsoa_sp$centr_lng,lsoa_sp$centr_lat),coordinates(all.events))
  rownames(new.dist_M) = lsoa_sp$code
  colnames(new.dist_M) = c(event_sp$course,new.events.names) 
  
  # AND EXTRACT NEW SHORTEST DISTANCES TO EVENTS
  new.lsoa_min_dist = apply(new.dist_M,1, function(x){round(min(x)/1000,3)} )
  diff_dist =  lsoa_sp$mn_dstn - new.lsoa_min_dist
  
# 3.1  how many affected? By how much?  ------------
  summary.stats(diff_dist) # distance change stats
  sum(diff_dist>0) # changes.n.lsoa
  round(sum(diff_dist!=0) / length(diff_dist),4)*100 # changes.perc.lsoa
  sum(lsoa_sp$pop[diff_dist>0]) # changes n population
  round(sum(lsoa_sp$pop[diff_dist!=0]) / sum(lsoa_sp$pop),4)*100 # changes % pop
  
  
# 3.2 CREATE NEW DATA FRAME   ------------
  # old distances
  df.olddistances = data.frame(mn_dstn = lsoa_sp$mn_dstn,
                               imd_sc = lsoa_sp$imd_sc,
                               pop = lsoa_sp$pop, week = 1)
  
  # new distances
  df.newdistances = data.frame(mn_dstn = new.lsoa_min_dist,
                               imd_sc = lsoa_sp$imd_sc,
                               pop = lsoa_sp$pop, week = 1)
  
  pred_runs_old_dist = predict(m_dist_imd_runs,newdata = df.olddistances,type="response")
  pred_runs_new_dist = predict(m_dist_imd_runs,newdata = df.newdistances,type="response")
  diff_runs  = pred_runs_new_dist- pred_runs_old_dist
  
  # change in runs stats
  summary.stats(diff_runs)
  
  # sum before, after, diff, proportion
  old_runs_sum = sum(pred_runs_old_dist)
  new_runs_sum = sum(pred_runs_new_dist)
  diff_runs_sum = sum(diff_runs)
  new_runs_sum/old_runs_sum
  
# 3.3 change by imd ------------
  ### increase botom 10 % imd
  imd_bottom_10 = quantile(lsoa_sp$imd_sc,probs = c(0.9))
  imd_bottom_10.index = lsoa_sp$imd_sc >= imd_bottom_10
  sum(diff_runs[imd_bottom_10.index])
  ### increase top 10 % imd
  imd_top_10 = quantile(lsoa_sp$imd_sc,probs = c(0.1))
  imd_top_10.index = lsoa_sp$imd_sc <= imd_top_10
  sum(diff_runs[imd_top_10.index])
  
  # proportions and ratio
  sum(diff_runs[imd_bottom_10.index])/(sum(diff_runs))
  sum(diff_runs[imd_top_10.index])/(sum(diff_runs))
  sum(diff_runs[imd_top_10.index])/sum(diff_runs[imd_bottom_10.index])
  
# / -------
# 4 CREATE STATIC MAPS FOR PUBLICATION -------
  # the functions below use google maps API to retrieve the baseline tile
  # to use this function, you need a api key
  # For more details, see: https://developers.google.com/maps/documentation/maps-static/intro
  my_api_key = "ENTER_YOUR_API_KEY"
  
# 4.1 RUN IMD WEIGHTS -----
  pr_map_static = create_static_map(new_lon=opt_loc_runs_imd.prep@data$lon,
                                    new_lat=opt_loc_runs_imd.prep@data$lat,
                                    rank=opt_loc_runs_imd.prep@data$pos,
                                    google_api_key = my_api_key,
                                    old_events_sp = event_sp,
                                    england_poly = england_sp,
                                    space.left = 0.1,
                                    space.right = 0.3,
                                    space.top = 0.2,
                                    space.bottom = 0.2)
  ggsave(plot=pr_map_static,filename = "./output/runs_imd_map_static.jpeg",height = 10,width = 8)
  # + legend
  map1_legend = legend_builder(opt_loc_runs_imd.prep)
  write.csv(map1_legend,file="./output/map1_legend.csv",row.names = F)

# 4.2 RUNS TOTAL-----
  pr2_map_static = create_static_map(new_lon=opt_loc_runs.prep@data$lon,
                                    new_lat=opt_loc_runs.prep@data$lat,
                                    rank=opt_loc_runs.prep@data$pos,
                                    google_api_key = my_api_key,
                                    old_events_sp = event_sp,
                                    england_poly = england_sp,
                                    space.left = 0.1,
                                    space.right = 0.3,
                                    space.top = 0.2,
                                    space.bottom = 0.2)
  ggsave(plot=pr2_map_static,filename = "./output/runs_map_static.jpeg",height = 10,width = 8)
  # + legend
  map2_legend = legend_builder(opt_loc_runs.prep)
  write.csv(map2_legend,file="./output/map2_legend.csv",row.names = F)
  
# 4.3 DISTANCE IMD WEIGHTS -----
  pr3_map_static = create_static_map(new_lon=opt_loc_dist_imd.prep@data$lon,
                                     new_lat=opt_loc_dist_imd.prep@data$lat,
                                     rank=opt_loc_dist_imd.prep@data$pos,
                                     google_api_key = my_api_key,
                                     old_events_sp = event_sp,
                                     england_poly = england_sp,
                                     space.left = 0.1,
                                     space.right = 0.3,
                                     space.top = 0.2,
                                     space.bottom = 0.2)
  ggsave(plot=pr3_map_static,filename = "./output/dist_imd_map_static.jpeg",height = 10,width = 8)
  # + legend
  map3_legend = legend_builder(opt_loc_dist_imd.prep)
  write.csv(map3_legend,file="./output/map3_legend.csv",row.names = F)
  
# 4.4 DISTANCE TOTAL -----
  pr4_map_static = create_static_map(new_lon=opt_loc_dist.prep@data$lon,
                                     new_lat=opt_loc_dist.prep@data$lat,
                                     rank=opt_loc_dist.prep@data$pos,
                                     google_api_key = my_api_key,
                                     old_events_sp = event_sp,
                                     england_poly = england_sp,
                                     space.left = 0.1,
                                     space.right = 0.3,
                                     space.top = 0.2,
                                     space.bottom = 0.2)
  ggsave(plot=pr4_map_static,filename = "./output/dist_map_static.jpeg",height = 10,width = 8)
  # + legend
  map4_legend = legend_builder(opt_loc_dist.prep)
  write.csv(map4_legend,file="./output/map4_legend.csv",row.names = F)

# All done ---  thank you!  
  

  
  