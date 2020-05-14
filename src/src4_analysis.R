
# ANALYSIS src

# 1  Descriptive stats ---------------
     
  # LSOA
    length(lsoa_sp$code) # n LSOAs
    summary.stats(lsoa_sp$pop) # population
    summary.stats(lsoa_sp$imd_sc) # deprivation indices
    summary.stats(lsoa_sp$mn_dstn) # distances to the nearest event
  
  # events
    length(event_sp$course) # n parkrun events
    summary.stats(event_sp$srvd_pop) # served population
    summary.stats(event_sp$srvd_lsoa) # served LSOA
  
  # green spaces    
    summary.stats(greens_sp$area_km2) # green space areas
    
  # pop living_within_xkm catchment 
    catch_pop(1)
    catch_pop(5)
    catch_pop(10)
    catch_pop(20)
    
    # descriptive stats new distances
    sum(lsoa_sp$mn_dstn_diff != 0) # lsoas affected
    sum(lsoa_sp$pop[lsoa_sp$mn_dstn_diff != 0]) # pop affected
    sum(lsoa_sp$pop[lsoa_sp$mn_dstn_diff != 0])/sum(lsoa_sp$pop) # pop %
    summary.stats(lsoa_sp$mn_dstn_diff[lsoa_sp$mn_dstn_diff != 0])
    summary.stats(lsoa_sp$mn_dstn_diff)
    
    
# 2  Association between IMD and Distance ---------------
    
    # as is
    fit1 = lm(mn_dstn ~ imd_sc,lsoa_sp,weights = pop)
    summary(fit1)
    weightedCorr(lsoa_sp$imd_sc, lsoa_sp$mn_dstn,method = "pearson",weights = lsoa_sp$pop)
    weightedCorr(lsoa_sp$imd_sc, lsoa_sp$mn_dstn,method = "spearman",weights = lsoa_sp$pop)
    
    mndst_tbl = by(lsoa_sp$mn_dstn,lsoa_sp$imd_q5,summary.stats)
    mndst_tbl = as.data.frame(do.call("rbind",mndst_tbl))
    # mndst_tbl
    
    # after 200 events are added
    fit2 = lm(mn_dstn_new ~ imd_sc,lsoa_sp,weights = pop)
    summary(fit2)
    weightedCorr(lsoa_sp$imd_sc, lsoa_sp$mn_dstn_new,method = "pearson",weights = lsoa_sp$pop)
    weightedCorr(lsoa_sp$imd_sc, lsoa_sp$mn_dstn_new,method = "spearman",weights = lsoa_sp$pop)
    mndst_new_tbl = by(lsoa_sp$mn_dstn_new,lsoa_sp$imd_q5,summary.stats)
    mndst_new_tbl = as.data.frame(do.call("rbind",mndst_new_tbl))
    # mndst_new_tbl
    
    catch_pop(1,dist = lsoa_sp$mn_dstn_new)
    catch_pop(5,dist = lsoa_sp$mn_dstn_new)
    catch_pop(10,dist = lsoa_sp$mn_dstn_new)
    catch_pop(15,dist = lsoa_sp$mn_dstn_new)
    
# 3  Distance diff table ---------------
    
    diff_tbl = by(lsoa_sp$mn_dstn_diff,lsoa_sp$imd_q5,summary.stats)
    diff_tbl = as.data.frame(do.call("rbind",diff_tbl))
    diff_tbl_lsoa = by(lsoa_sp,lsoa_sp$imd_q5,function(x){length(x$pop[x$mn_dstn_diff != 0])})
    diff_tbl_pop = by(lsoa_sp,lsoa_sp$imd_q5,function(x){sum(x$pop[x$mn_dstn_diff != 0])})
    diff_tbl  = cbind(diff_tbl,afected.pop = as.numeric(diff_tbl_pop),affected.lsoa=as.numeric(diff_tbl_lsoa))
    # diff_tbl
    

# 4  Density plot  ---------------
    # 99 rows removed
    before_after_dist = ggplot(lsoa_sp@data) +
      geom_density(aes(mn_dstn_new),col="black",fill="lightgreen",alpha=0.5) +
      geom_density(aes(mn_dstn),col="black",fill="red",alpha=0.3) +
      scale_x_continuous(label=label_number(suffix = " km"),limits = c(0,25),
                         name = "Distance to the nearest parkrun event") +
      theme_minimal()

      ggsave(plot=before_after_dist,filename = "./output/hist.png",height=6,width=8)
      
      
# 5  Scatter plot  ---------------
      
# relationship dist ~ imd plots
    lsoa_melt = melt(lsoa_sp@data[,c("mn_dstn","mn_dstn_new","imd_ra","imd_sc")],id.vars = c("imd_sc","imd_ra"))
    
    # ggplot(lsoa_melt) +
    #   geom_point(aes(x=(imd_ra),y=(value),col=variable),alpha=0.5,size=0.2) +
    #   geom_smooth(aes(x=(imd_ra),y=(value)),method="lm") +
    #   facet_wrap(~variable) +
    #   theme_minimal() +
    #   theme(legend.position = "none") +
    #   ylim(c(0,25))
    
    levels(lsoa_melt$variable) = c("As is (12th December 2018)","After 200 new events are created")
    
    scatter =  ggplot(lsoa_melt) +
        geom_point(aes(x=(imd_sc),y=(value),col=variable),alpha=0.5,size=0.2) +
        geom_smooth(aes(x=(imd_sc),y=(value)),method="lm") +
        facet_wrap(~variable,) +
        theme_minimal() +
        theme(legend.position = "none")+
      scale_y_continuous(label=label_number(suffix = " km"),limits = c(0,25),
                         name = "Distance to the nearest parkrun event") +
      xlab("Index of multiple deprivation") 
      
    ggsave(plot = scatter,filename = "./output/scatter.png",width = 10,height = 4)
    
  


# fin.
  

  
  