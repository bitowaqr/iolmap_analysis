
### OPTIMAL NEW PARKRUN LOCATIONS
source("./src/src1_functions.R")
source("./src/src2_data_proc.R")
    
    time1 = Sys.time() # keeps track of time
    
    # specify GLM to estimate the effect of new runs on LSOAs
    m_dist_imd_runs = glm(run_count ~ imd_sc + mn_dstn,
                          data = lsoa_sp@data,
                          family = poisson(link="log"),
                          offset = log(pop*week))
    
    # model for distance only
    m_dist_dist = lm(mn_dstn+0 ~ -1 + mn_dstn ,data = lsoa_sp@data)
    
    # new data frame for predictions
    imd_runs_df = lsoa_sp@data[,c("pop","imd_sc","week")]
   
# MAX TOTAL RUNS IMD WEIGHTED
    opt_loc_runs_imd.prep = LocAloAlg(candidates = greens_sp,
                                 event_coordinates.M = event_sp, 
                                 lsoa_centr = cbind(lsoa_sp@data$centr_lng,
                                                    lsoa_sp@data$centr_lat),
                                 data = imd_runs_df,
                                 model = m_dist_imd_runs,
                                 weights = lsoa_sp@data$imd_sc^2, 
                                 top_n = 200)
    shapefile(opt_loc_runs_imd.prep,filename="./output/opt_loc_runs_imd.prep.shp",overwrite = T)
    gc()
    gc()
    
# MAX TOTAL RUNS
    opt_loc_runs.prep = LocAloAlg(candidates = greens_sp,
                             event_coordinates.M = event_sp, 
                             lsoa_centr = cbind(lsoa_sp@data$centr_lng,
                                                lsoa_sp@data$centr_lat),
                             data = imd_runs_df,
                             model = m_dist_imd_runs,
                             weights = 1, 
                             top_n = 200,
                             maximise = T)
    
    shapefile(opt_loc_runs.prep,filename="./output/opt_loc_runs.prep.shp",overwrite = T)
    gc()
    gc()
    
# MIN DISTANCE IMD WEIGHTED
    opt_loc_dist_imd.prep = LocAloAlg(candidates = greens_sp,
                                 event_coordinates.M = event_sp, 
                                 lsoa_centr = cbind(lsoa_sp@data$centr_lng,
                                                    lsoa_sp@data$centr_lat),
                                 data = imd_runs_df,
                                 model = m_dist_dist,
                                 weights = lsoa_sp@data$imd_sc^2, 
                                 top_n = 200,
                                 maximise = F)
    shapefile(opt_loc_dist_imd.prep,filename="./output/opt_loc_dist_imd.prep.shp",overwrite = T)
    
    gc()
    gc()
    
# MIN DISTANCE IMD WEIGHTED
    opt_loc_dist.prep = LocAloAlg(candidates = greens_sp,
                             event_coordinates.M = event_sp, 
                             lsoa_centr = cbind(lsoa_sp@data$centr_lng,
                                                lsoa_sp@data$centr_lat),
                             data = imd_runs_df,
                             model = m_dist_dist,
                             weights = 1, 
                             top_n = 200,
                             maximise = F)
    shapefile(opt_loc_dist.prep,filename="./output/opt_loc_dist.prep.shp",overwrite = T)
    
    gc()
    gc()
    

# DONE
    time2 = Sys.time()
    time2 - time1 # runtime
    