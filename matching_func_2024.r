options(warn=-1)
options(dplyr.summarise.inform = FALSE)

# ppath <- "/gpfs/data1/duncansongp/amberliang/PADDDtracker_DataReleaseV2_May2019/"
# poly <- readOGR(paste(ppath,"PADDDtracker_DataReleaseV2_May2019_Poly.shp",sep=""),verbose = FALSE) %>% spTransform(., CRS("+init=epsg:6933"))
# pts <- readOGR(paste(ppath,"PADDDtracker_DataReleaseV2_May2019_Pts.shp",sep=""), verbose = FALSE) %>% spTransform(., CRS("+init=epsg:6933"))
# poly$Location_K <- toupper(poly$Location_K)
# poly <- poly[poly$Location_K =="Y",]
# poly$EventType <- as.numeric(poly$EventType) 
# poly$EventType[is.na(poly$EventType)] <- 0
# 
# pts$Location_K <- toupper(pts$Location_K)
# pts <- pts[pts$Location_K =="Y",]
# pts$EventType <- as.numeric(pts$EventType) 
# pts$EventType[is.na(pts$EventType)] <- 0
`%notin%` <- Negate(`%in%`)

# Function to allow rbinding dataframes with foreach even when some dataframes 
# may not have any rows
foreach_rbind <- function(d1, d2) {
  if (is.null(d1) & is.null(d2)) {
    return(NULL)
  } else if (!is.null(d1) & is.null(d2)) {
    return(d1)
  } else if (is.null(d1) & !is.null(d2)) {
    return(d2)
  } else  {
    return(rbind(d1, d2))
  }
}

match_wocat <- function(df, pid) {
  
  registerDoParallel(4) #_@_@
  
  options("optmatch_max_problem_size"=Inf)
  
  # Filter out countries without at least one treatment unit or without at
  # least one control unit
  # df <- df %>%
  #   filter(complete.cases(.)) %>%
  
  #the following lines do nothing because only one PA at a time
  #mutate(n_treatment=sum(status),
  #       n_control=sum(!status)) %>%
  #the next line doesn't do 
  #filter(n_treatment >= 1, n_control >= 1)
  
  # Note custom combine to handle iterations that don't return any value
  #test nested foreach loops
  ret <- foreach (this_lc=unique(df$land_cover),
                  .packages=c('optmatch', 'dplyr'),
                  .combine=foreach_rbind, .inorder=FALSE) %dopar% { #_@_@
                    this_d<-df #****Both of these are d_all, make sure these have same info as d_all below****
                    d_wocat <- filter(this_d, status)
                    # Filter out climates and land covers that don't appear in the wocat
                    # sample, and drop these levels from the factors
                    this_d <- filter(this_d, #***** AFTER SCORING EXPORT AS CSV
                                     land_cover %in% unique(d_wocat$land_cover),
                                     wwfbiom %in% unique(d_wocat$wwfbiom),
                                     wwfecoreg %in% unique(d_wocat$wwfecoreg))
                    
                    this_d$land_cover <- droplevels(this_d$land_cover)
                    this_d$wwfbiom <- droplevels(this_d$wwfbiom)
                    this_d$wwfecoreg <- droplevels(this_d$wwfecoreg)
                    # table(this_d$status)
                   dat <- dplyr::select(this_d, lat, lon, UID, status, land_cover, wwfbiom, wwfecoreg, elevation, slope,
                         mean_temp,max_temp,min_temp, prec, d2road, d2city,  popden, tt2city, popcnt) 
                    ps_A <- glm(status ~ mean_temp+max_temp+min_temp + prec + elevation + slope+ d2road + d2city + popden +popcnt+ tt2city,data = dat)
                    dat$propensity_scoreA <- fitted(ps_A)
                    f <- status~ mean_temp + max_temp + min_temp + prec + elevation + slope + d2road + d2city + popden + popcnt + tt2city
                    this_d<- bind_cols(this_d,propensity_scoreA=dat$propensity_scoreA)
                    # Can't stratify by land cover or climate if they only have one level #
                    if (nlevels(this_d$land_cover) >= 2) {
                      f <- update(f, ~ . + strata(land_cover))
                    } else {
                      f <- update(f, ~ . - land_cover)
                    }
                    if (nlevels(this_d$wwfbiom) >= 2) {
                      f <- update(f, ~ . + strata(wwfbiom))
                    } else {
                      f <- update(f, ~ . - wwfbiom)
                    }
                    if (nlevels(this_d$wwfecoreg) >= 2) {
                      f <- update(f, ~ . + strata(wwfecoreg))
                    } else {
                      f <- update(f, ~ . - wwfecoreg)
                    }
                    if (nrow(d_wocat) > 2) {
                      model <- glm(f, data=this_d)
                      dists <- match_on(model, data=this_d) #***NEED PROPENSITY SCORE BEFORE HERE
                    } else {
                      # Use Mahalanobis distance if there aren't enough points to run a glm
                      dists <- match_on(f, data=this_d)
                    }
                    #potentially drop caliper line; will cut down dists matrix but not the speed issue
                    # dists <- caliper(dists, 2)
                    # If the controls are too far from the treatments (due to the caliper) 
                    # then the matching may fail. Can test for this by seeing if subdim 
                    # runs successfully
                    subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                                             error=function(e)return(FALSE))
                    if (subdim_works) {
                      m <- fullmatch(dists, min.controls=1, max.controls=1, data=this_d)
                      prematch_d <- this_d
                      this_d$matched <- m
                      this_d <- this_d[matched(m), ]
                    } else {
                      this_d <- data.frame()
                    }
                    # Need to handle the possibility that there were no matches for this 
                    # treatment, meaning this_d will be an empty data.frame
                    if (nrow(this_d) == 0) {   #if matching w/ ecoreg return no results, match again w/o ecoreg and check matching results
                      this_d<-df
                      d_wocat <- filter(this_d, status)
                      this_d <- filter(this_d,
                                       land_cover %in% unique(d_wocat$land_cover),
                                       wwfbiom %in% unique(d_wocat$wwfbiom))
                      
                      this_d$land_cover <- droplevels(this_d$land_cover)
                      this_d$wwfbiom <- droplevels(this_d$wwfbiom)
                      f <- status ~ mean_temp + max_temp + min_temp + prec + elevation + slope + d2road + d2city + popden + popcnt + tt2city
                      if (nlevels(this_d$land_cover) >= 2) {
                        f <- update(f, ~ . + strata(land_cover))
                      } else {
                        f <- update(f, ~ . - land_cover)
                      }
                      if (nlevels(this_d$wwfbiom) >= 2) {
                        f <- update(f, ~ . + strata(wwfbiom))
                      } else {
                        f <- update(f, ~ . - wwfbiom)
                      }
                      if (nrow(d_wocat) > 2) {
                        model <- glm(f, data=this_d)
                        dists <- match_on(model, data=this_d)
                      } else {
                        dists <- match_on(f, data=this_d)
                      }
                      subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                                               error=function(e)return(FALSE))
                      if (subdim_works) {
                        m <- fullmatch(dists, min.controls=1, max.controls=1, data=this_d)
                        prematch_d <- this_d
                        this_d$matched <- m
                        this_d <- this_d[matched(m), ]
                      } else {
                        this_d <- data.frame()
                      }
                      if(nrow(this_d)==0){
                        log_mes <- paste(pid,"-Matching without wwfecoreg:Failed\n",sep="")
                        filepath<-file.path(f.path3,paste0("/WDPA_matching_log/",iso3))
                        dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
                        cat(log_mes,file=paste(f.path3,"/WDPA_matching_log/",iso3,"/",iso3,"_pa_",pid,"_matching_used_covar_log_wk", gediwk,".txt",sep=""),append=TRUE)
                        return(NULL)
                      } else{
                        log_mes <- paste(pid,"-Matching without wwfecoreg:Succeed\n",sep="")
                        filepath<-file.path(f.path3,paste0("/WDPA_matching_log/",iso3))
                        dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
                        cat(log_mes,file=paste(f.path3,"/WDPA_matching_log/",iso3,"/",iso3,"_pa_",pid,"_matching_used_covar_log_wk", gediwk,".txt",sep=""),append=TRUE)
                        match_results <- list("match_obj" = m, "df" = this_d, "func"=f, "prematch_d"=prematch_d)
                        return(match_results)
                      }
                    } else {
                      log_mes <- paste(pid,"-Matching with wwfecoreg:Succeed\n",sep="")
                      match_results <- list("match_obj" = m, "df" = this_d, "func"=f, "prematch_d"=prematch_d)
                      filepath<-file.path(f.path3,paste0("/WDPA_matching_log/",iso3))
                      dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
                      cat(log_mes,file=paste(f.path3,"/WDPA_matching_log/",iso3,"/",iso3,"_pa_",pid,"_matching_used_covar_log_wk", gediwk,".txt",sep=""),append=TRUE)
                      return(match_results)
                    }
                  }
  
  stopImplicitCluster()  #_@_@
  return(ret)
}

propensity_filter <- function(pa_df, d_control_local){
  pa_df <-pa_df[complete.cases(pa_df), ]  #filter away non-complete cases w/ NA in control set
  d <- dplyr::bind_rows(d_control_local, pa_df)
  ## bring in matching algorithm from STEP5 here to loop through each PA in d_PAs
  #filter controls based on propensity scores 
  d_all <- dplyr::select(d, lat, lon, UID, status, land_cover, wwfbiom, wwfecoreg, elevation, slope,
                         mean_temp,max_temp,min_temp, prec, d2road, d2city,  popden, tt2city, popcnt) 
  #**** HOMEWORK ABOVE Add ps model object and propensity score to above og dataframe, add to wocat function
  d_all$status <- ifelse(d_all$status==TRUE,1,0)
  
  #calculate the propensity scores & filter out controls not overlapping w/ treatment propensity scores
  ps <- glm(status ~ mean_temp+max_temp+min_temp + prec + elevation + slope+ d2road + d2city + popden +popcnt+ tt2city,data = d_all)
  # boxplot(ps)  #check the distribution of propensity scores for treatment and controls
  #filter out the controls with propensity scores outside of the overlapping region
  d_all$propensity_score <- fitted(ps)
  filepath<-file.path(f.path3,paste0(iso3,"_wk",gediwk,"_prefilter_ps/"))
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  # dir.create(paste(f.path3,iso3,"_wk",gediwk,"_prefilter_ps/", sep=""), recursive = TRUE, showWarnings = FALSE)
  write.csv(d_all, paste(filepath,iso3,"_pa_",id_pa,"_pre_filter_ps_wk",gediwk,".csv", sep=""))

  d_sep <- d_all %>% dplyr::group_by(status)
  d_sep_range <- d_all %>% dplyr::group_by(status)%>% 
    dplyr::summarise(propmin= min(propensity_score), promax=max(propensity_score)) 
  # cat(iso3, "Filtering the control sites by overlaping with the treatment PS\n")
  d_filtered <- d_sep %>% 
    filter(status==1 | between(propensity_score,d_sep_range$propmin[2],d_sep_range$promax[2])) %>% 
    ungroup() 
  
  d_filtered$status <- ifelse(d_filtered$status==1,TRUE,FALSE)
  
  return(d_filtered)
}

# isoPadddRas <- function(poly, pts, rtemplate){
#   
#   if((iso3 %in% unique(poly$ISO3166))&&(iso3 %notin% unique(pts$ISO3166))){
#     # print("in poly")
#     polysub <- poly[poly$ISO3166==iso3,]
#     polysubr <- rasterize(polysub,rtemplate,background=NA, field=polysub$EventType)
#     names(polysubr) <- "PADDD"
#     return(polysubr)
#     
#   } else if ((iso3 %notin% unique(poly$ISO3166))&&(iso3 %in% unique(pts$ISO3166))){
#     # print("in pts")
#     ptssub <- pts[pts$ISO3166==iso3,]
#     ptssubr <- rasterize(ptssub@coords[,1:2,drop=FALSE],rtemplate,background=NA, field=ptssub$EventType)
#     names(ptssubr) <- "PADDD"
#     return(ptssubr)
#     
#   } else if ((iso3 %in% unique(poly$ISO3166))&&(iso3 %in% unique(pts$ISO3166))){
#     # print("in both")
#     polysub <- poly[poly$ISO3166==iso3,]
#     polysubr <- rasterize(polysub,rtemplate,background=NA, field=polysub$EventType)
#     ptssub <- pts[pts$ISO3166==iso3,]
#     ptssubr <- rasterize(ptssub@coords[,1:2,drop=FALSE],rtemplate,background=NA, field=ptssub$EventType)
#     m <- merge(polysubr,ptssubr)
#     names(m) <-  "PADDD"
#     return(m)
#   } else {
#     # print("in neither")
#     empr <- rtemplate
#     values(empr) <- NA
#     names(empr) <- "PADDD"
#     return(empr)
#   }
# }
matched2ras <- function(matched_df) {
  
  cat("Converting the matched csv to a raster stack for extraction\n")
  
  # matched_pts <- SpatialPointsDataFrame(coords=matched_df[,c("lon","lat")],
  #                                         proj4string=CRS("+init=epsg:4326"), data=matched_df) %>% 
  #                                         spTransform(., CRS("+init=epsg:6933"))
  matched_pts<- vect(matched_df)
  crs(matched_pts)<-"epsg:4326"
  matched_pts<-project(matched_pts, "epsg:6933")
  
  # Ensure fields are in appropriate formats
  matched_pts$UID <- as.integer(matched_pts$UID)
  matched_pts$pa_id <- as.integer(matched_pts$pa_id)
  matched_pts$status <- as.logical(matched_pts$status)
  matched_pts$wwfbiom <- as.numeric(matched_pts$wwfbiom)
  matched_pts$wwfecoreg <- as.numeric(matched_pts$wwfecoreg)
  
  # Fill NA values with 0
  matched_pts$wwfbiom[is.na(matched_pts$wwfbiom)] <- 0
  matched_pts$wwfecoreg[is.na(matched_pts$wwfecoreg)] <- 0
  
  # Define the raster extent for cropping
  buffer_ext <- ext(buffer(matched_pts, 10000))
  r <- crop(MCD12Q1, buffer_ext)
  continent <- crop(world_region, buffer_ext)
  
  # Assign names to the raster layers
  names(r) <- "pft"
  names(continent) <- "region"
  
  # List of fields to rasterize
  fields <- c("status", "pa_id", "wwfbiom", "wwfecoreg", "UID")
  rasters <- list()
  
  # Rasterize each field
  for (field in fields) {
    r_field <- rasterize(matched_pts, r, field = field)
    rasters[[field]] <- r_field
  }
  
  # Combine rasters into a SpatRaster stack
  rasters <- c(rasters, list(pft = r, region = continent))
  matched_ras <- rast(rasters)
  
  return(matched_ras)
}
                                               

convertFactor <- function(matched0, exgedi){
  exgedi$pft <- as.character(exgedi$pft)
  
  exgedi$pft <- factor(exgedi$pft, levels=sequence(6),
                                 labels = c("ENT",
                                            "EBT",
                                            "ENT",
                                            "DBT",
                                            "GS",
                                            "GS"))
  exgedi$region <- as.character(exgedi$region)
  exgedi$region <- factor(exgedi$region, levels=c(1:7),
                       labels = c("Eu",
                                  "As",
                                  "Au",
                                  "Af",
                                  "As",
                                  "SA",
                                  "US"))
  
  exgedi$stratum <- paste(exgedi$pft, exgedi$region,sep="_")
  
  # exgedi$GOV_TYPE <- exgedi$GOV_TYPE %>% 
  #   factor(levels=seq(length(levels(matched0$GOV_TYPE))),
  #          labels=levels(matched0$GOV_TYPE))
  # 
  # exgedi$OWN_TYPE <- exgedi$OWN_TYPE %>% 
  #   factor(levels=seq(length(levels(matched0$OWN_TYPE))),
  #          labels=levels(matched0$OWN_TYPE))  
  
  # exgedi$DESIG_ENG <- exgedi$DESIG_ENG %>% 
  #   factor(levels=seq(length(levels(matched0$DESIG_ENG))),
  #          labels=levels(matched0$DESIG_ENG))  
  
  exgedi$wwfbiom <- exgedi$wwfbiom %>% 
    factor(levels=seq(length(levels(matched0$wwfbiom))),
           labels=levels(matched0$wwfbiom))  
  
  exgedi$wwfecoreg <- exgedi$wwfecoreg %>% 
    factor(levels=seq(length(levels(matched0$wwfecoreg))),
           labels=levels(matched0$wwfecoreg))  
  
  # tryCatch(exgedi$paddd <- as.character(exgedi$paddd), error=function(e) return(NULL))
  # tryCatch(exgedi$paddd[which(exgedi$paddd=="1")] <- "Downgrade", error=function(e) return(NULL))
  # tryCatch(exgedi$paddd[which(exgedi$paddd=="2")] <- "Degazette", error=function(e) return(NULL))
  # tryCatch(exgedi$paddd[which(exgedi$paddd=="3")] <- "Downsize", error=function(e) return(NULL))
  
  return(exgedi)
}

subdfExport <- function(filtered_df){
  #export invidual pa results
  spt2 <- split(filtered_df, filtered_df$pa_id)
  
  dfl <- lapply(names(spt2), function(x){
  
    if(dim(spt2[[x]])[1]>0){
      control_sub <- spt2[[x]][spt2[[x]]$status==0,]
      treat_sub <-  spt2[[x]][spt2[[x]]$status==1,]
      ncontrol <- nrow(control_sub)
      ntreat <- nrow(treat_sub)
      if (ntreat-ncontrol > 0){
        newtreatid <- sample(ntreat, ncontrol)
        newtreat <- treat_sub[newtreatid,]
        spt2_new <- rbind(newtreat, control_sub)
      } else if (ntreat - ncontrol< 0){
        newcontrolid <- sample(ncontrol, ntreat)
        newcontrol <- control_sub[newcontrolid,]
        spt2_new <- rbind(newcontrol, treat_sub)
      } else if (ntreat-ncontrol==0){
        spt2_new <- spt2[[x]]
      } else if (ntreat==0 || ncontrol==0){
        spt2_new=NA
      }
      biom <- spt2_new$wwfbiom %>% unique() %>% as.character() %>% gsub('\\b(\\pL)\\pL{4,}|.','\\U\\1',.,perl = TRUE)
      if(length(biom)>1){
        biom <- paste(c(biom), collapse="&")
      }
      # print(biom)
      write.csv(spt2_new, file=paste(f.path3,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_PA_",unique(spt2_new$pa_id),"_",biom,".csv", sep=""))
      return(spt2_new)
    }
  })
  
  total_df <- do.call("rbind", dfl) 
  cat("Exported individual PAs results for ", iso3, "\n")
  
  return(total_df)
}

getmode <- function(v,na.rm) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# extract_gedi <- function(matched, mras){
#     lon_bond <- range(matched$lon,na.rm=TRUE)
#     lat_bond <- range(matched$lat,na.rm=TRUE)
#     all_gedil2_f <- list.files(file.path(f.path,"SEN_L2A_Sub"), full.names = FALSE) 
#     all_gedil4_f <- list.files(file.path(f.path,"SEN_L4A_Sub"), full.names = FALSE) 
#     gedil2_f <- all_gedil2_f%>% strsplit( "_") %>% 
#       as.data.frame() %>% 
#       t() %>% as.data.frame(row.names =all_gedil2_f, stringsAsFactors=FALSE,make.names=FALSE) %>% dplyr::select(V3,V4) %>% 
#       mutate(lons=as.numeric(gsub('\\D','', V3)), ew= gsub('\\d','', V3) ) %>% 
#       mutate(lats= as.numeric(gsub('\\D','', V4)), ns= gsub('\\d','', V4) ) %>% 
#       mutate( lons = ifelse(ew!="E", -1*lons, lons)) %>% 
#       mutate( lats = ifelse(ns!="N", -1*lats, lats)) %>% 
#       dplyr::filter( between(lons, floor(lon_bond[1]), ceiling(lon_bond[2]))) %>% 
#       dplyr::filter(between(lats, floor(lat_bond[1]), ceiling(lat_bond[2]))) %>% rownames()
#     gedil4_f <- all_gedil4_f%>% strsplit( "_") %>% 
#       as.data.frame() %>% 
#       t() %>% as.data.frame(row.names =all_gedil4_f, stringsAsFactors=FALSE,make.names=FALSE) %>% dplyr::select(V3,V4) %>% 
#       mutate(lons=as.numeric(gsub('\\D','', V3)), ew= gsub('\\d','', V3) ) %>% 
#       mutate(lats= as.numeric(gsub('\\D','', V4)), ns= gsub('\\d','', V4) ) %>% 
#       mutate( lons = ifelse(ew!="E", -1*lons, lons)) %>% 
#       mutate( lats = ifelse(ns!="N", -1*lats, lats)) %>% 
#       dplyr::filter( between(lons, floor(lon_bond[1]), ceiling(lon_bond[2]))) %>% 
#       dplyr::filter(between(lats, floor(lat_bond[1]), ceiling(lat_bond[2]))) %>% rownames()  #should result in same # of files as the l2 products
    
#     registerDoParallel(cores=round(mproc*0.5))
#     ex_out <- foreach(this_csvid=seq(length(gedil2_f)), .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {
#         ##add the GEDI l4a model prediction for AGB here :
#         cat("Readng in no. ", this_csvid,"csv of ", length(gedil2_f),"csvs for iso3",iso3,"\n")
#         gedi_l2  <- read.csv(paste(f.path,"WDPA_gedi_l2a+l2b_clean2",iso3,gedil2_f[this_csvid], sep="/")) %>%
#           dplyr::select(shot_number,lon_lowestmode, lat_lowestmode,rh_025, rh_050, rh_075, rh_098,cover, pai)
#         l2_latlon <- gedil2_f[this_csvid] %>%  str_split("_") %>% unlist %>% .[3:4] %>% paste(sep="_", collapse ="_") %>% paste("_",.,sep="")
#         l4_pattern <- tryCatch(grep(l2_latlon, gedil4_f, value=TRUE), error=function(cond){return(NA)})
#         gedi_l4  <- tryCatch(read.csv(paste(f.path,"WDPA_gedi_l4a_clean",iso3,l4_pattern, sep="/")), error=function(cond){return(NA)})
                            
#         if (is.na(gedi_l4) || nrow(gedi_l4) < 1){
#           cat("error")
#           gedi_l24 <- gedi_l2
#           gedi_l24$agbd <- NA
#           gedi_l24$agbd_se <- NA
#           gedi_l24$agbd_t <- NA
#           gedi_l24$agbd_t_se <- NA
#         } else {
#           gedi_l4_sub <- gedi_l4 %>%
#             dplyr::select(shot_number, agbd, agbd_se, agbd_t, agbd_t_se)
#           gedi_l24 <- inner_join(gedi_l2, gedi_l4_sub, by="shot_number")
        
#         }
#       gedi_l24[rowSums(is.na(gedi_l24)) > 0, ]   
#       gedi_l24 <- left_join(gedi_l2, gedi_l4, by="shot_number") %>% drop_na()
#       iso_matched_gedi_df <- data.frame()
#       if(nrow(gedi_l24)>0){
#         gedi_l24_sp <- gedi_l24 %>% 
#           SpatialPointsDataFrame(coords=.[,c("lon_lowestmode","lat_lowestmode")],
#                                  proj4string=CRS("+init=epsg:4326"), data=.) %>%spTransform(., CRS("+init=epsg:6933"))
        
#         matched_gedi <- raster::extract(mras,gedi_l24_sp, df=TRUE)
#         matched_gedi_metrics <- cbind(matched_gedi,gedi_l24_sp@data)
#         matched_gedi_metrics_filtered <- matched_gedi_metrics %>% dplyr::filter(!is.na(status)) %>% 
#           convertFactor(matched0 = matched,exgedi = .) 
        
#         iso_matched_gedi_df <- rbind(matched_gedi_metrics_filtered,iso_matched_gedi_df)
        
#       }
#       return(iso_matched_gedi_df)
#   }
#   stopImplicitCluster()
#   cat("Done GEDI for no. ",grep(unique(matched$pa_id), matched_PAs),"pa out of", length(matched_PAs),"\n")
#   return(ex_out)
# }
                                      
#####################
## UPDATED EXTRACT GEDI
#######################
# extract_gedi <- function(matched, mras){
#     gedipath<- file.path("~/shared-buckets/abarenblitt/GEDI_global_PA_v2/")
#     all_gedil2_f <- list.files(file.path(gedipath,"WDPA_gedi_L2A_tiles/"), full.names = FALSE)
#     all_gedil4_f <- list.files(file.path(gedipath,"WDPA_gedi_L4A_tiles/"), full.names = FALSE)  
    
    
#     registerDoParallel(cores=1) #_@_@
#     ex_out <- foreach(this_csvid=seq(length(all_gedil2_f)), 
#                   .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {
#         ##add the GEDI l4a model prediction for AGB here :
#         cat("Readng in no. ", this_csvid,"csv of ", length(all_gedil2_f),"csvs for iso3",iso3,"\n")

#     # for (this_gedi4 in all_gedil4_f[this_csvid]){
#             gedil4_f <- as.data.frame(st_read(paste(gedipath,"WDPA_gedi_L4A_tiles/",all_gedil4_f[this_csvid],sep="")))
#     # }
#     # for (this_gedi2 in all_gedil2_f[this_csvid]){
#             gedil2_f <- as.data.frame(st_read(paste(gedipath,"WDPA_gedi_L2A_tiles/",all_gedil2_f[this_csvid],sep="")))
#     # }
#             if (nrow(gedil4_f) < 1){   #is.na(gedi_l4) || 
#               cat("error")
#               gedi_l24 <- gedil2_f
#               gedi_l24$agbd <- NA
#               gedi_l24$agbd_se <- NA
#               gedi_l24$agbd_t <- NA
#               gedi_l24$agbd_t_se <- NA
#             } else {
#               gedi_l4_sub <- gedil4_f %>%
#                 dplyr::select(shot_number, agbd, agbd_se, agbd_t, agbd_t_se)
#               gedi_l24 <- inner_join(gedil2_f, gedi_l4_sub, by="shot_number")

#             }

#         print(dim(gedi_l24))
#         iso_matched_gedi_df <- data.frame()
#         if(nrow(gedi_l24)>0){
#             gedi_l24_sp <- gedi_l24 %>% 
#                 SpatialPointsDataFrame(coords=.[,c("lon_lowestmode","lat_lowestmode")],
#                                      proj4string=CRS("+init=epsg:4326"), data=.) %>%spTransform(., CRS("+init=epsg:6933"))

#             matched_gedi <- terra::extract(mras,vect(gedi_l24_sp), df=TRUE)
#             matched_gedi_metrics <- cbind(matched_gedi,gedi_l24_sp@data)
#             matched_gedi_metrics_filtered <- matched_gedi_metrics %>% dplyr::filter(!is.na(status)) %>% 
#               convertFactor(matched0 = matched,exgedi = .) 

#             iso_matched_gedi_df <- rbind(matched_gedi_metrics_filtered,iso_matched_gedi_df)
#             print(dim(iso_matched_gedi_df))
#         }

#         return(iso_matched_gedi_df)
# }
    
#     stopImplicitCluster() #_@_@
#     cat("Done GEDI for no. ",grep(unique(matched$pa_id), matched_PAs),"pa out of", length(matched_PAs),"\n")
#     return(ex_out)
    
# } 
                                 
# ############################################################################################   
#  #*********************
#  # NON PARALLEL VERSION
#  # ********************
                                      
extract_gedi <- function(matched, mras){
    results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                            Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/WDPA_gedi_L2A_tiles/",sep=""))
    all_gedil2_f <- sapply(results$Contents, function(x) {x$Key})
    pattern=paste(".gpkg",sep="")
    all_gedil2_f <- grep(pattern, all_gedil2_f, value=TRUE)
    all_gedil2_f <- basename(all_gedil2_f)#[4:6] #Currently specifying working files
    
    results4 <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                                Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/WDPA_gedi_L4A_tiles/",sep=""))
    all_gedil4_f <- sapply(results4$Contents, function(x) {x$Key})
    pattern4=paste(".gpkg",sep="")
    all_gedil4_f <- grep(pattern4, all_gedil4_f, value=TRUE)
    all_gedil4_f <- basename(all_gedil4_f)#[4:6] #Currently specifying working files
  
    # Initialize an empty list to store results
    results_list <- list()
    iso_matched_gedi_df <- NULL # Initialize before loop

            # Iterate over the sequence of indices for your files
    for (this_csvid in seq_along(all_gedil2_f)) {
                cat("Reading in no. ", this_csvid, "csv of ", length(all_gedil2_f), "csvs for iso3", iso3, "\n")
                
                # Read GEDI L4A data
                gedil4_f_path <- paste(gedipath, "WDPA_gedi_L4A_tiles/", all_gedil4_f[this_csvid], sep = "")
                gedil4_f <- as.data.frame(st_read(gedil4_f_path))
                
                # Read GEDI L2A data
                gedil2_f_path <- paste(gedipath, "WDPA_gedi_L2A_tiles/", all_gedil2_f[this_csvid], sep = "")
                gedil2_f <- as.data.frame(st_read(gedil2_f_path))
            
                # Check if GEDI L4A data is empty
                if (nrow(gedil4_f) < 1) {
                    cat("Error: No data for GEDI L4A\n")
                    gedi_l24 <- gedil2_f
                    gedi_l24$agbd <- NA
                    gedi_l24$agbd_se <- NA
                    gedi_l24$agbd_t <- NA
                    gedi_l24$agbd_t_se <- NA
                } else {
                    # Select relevant columns from GEDI L4A
                    gedi_l4_sub <- gedil4_f %>%
                        dplyr::select(shot_number, agbd, agbd_se, agbd_t, agbd_t_se)
                    
                    # Join with GEDI L2A data
                    gedi_l24 <- inner_join(gedil2_f, gedi_l4_sub, by = "shot_number")
                }
            
                print(dim(gedi_l24))
            
                # Initialize empty spatial object for the current iteration
                gedi_l24_sp <- NULL
            
                # Convert to spatial points data frame if there is data
                if (nrow(gedi_l24) > 0) {
                    gedi_l24_sp <- SpatialPointsDataFrame(
                        coords = gedi_l24[, c("lon_lowestmode", "lat_lowestmode")],
                        data = gedi_l24,
                        proj4string = CRS("+init=epsg:4326")
                    ) %>% spTransform(CRS("+init=epsg:6933"))
                matched_gedi <- terra::extract(mras,vect(gedi_l24_sp), df=TRUE)
                matched_gedi_metrics <- cbind(matched_gedi,gedi_l24_sp@data)
                matched_gedi_metrics_filtered <- matched_gedi_metrics %>% dplyr::filter(!is.na(status)) %>% 
                convertFactor(matched0 = matched,exgedi = .) 

            iso_matched_gedi_df <- rbind(matched_gedi_metrics_filtered,iso_matched_gedi_df)
            print(dim(iso_matched_gedi_df))
         }
        
        # Store results in a list
        results_list[[this_csvid]] <- iso_matched_gedi_df
    }
    
    # Combine all results
    if (!is.null(iso_matched_gedi_df)) {
        iso_matched_gedi_df <- do.call(rbind, results_list)
    }
    
    cat("Done GEDI processing\n")
    return(iso_matched_gedi_df)
}

 ############################################################################################   
 #*************************************
 # NON PARALLEL VERSION ADDING IN L2B
 # **************************************

                                               
extract_gedi2b <- function(iso3,tile_id,f.path3,gedipath){
    # Initialize an empty list to store results
  # results_list <- list()
  ### TODO: Are you sure you need the next line?
  iso_matched_gedi_df <- NULL # Initialize before loop
 
   filepath <- file.path(f.path3, paste(iso3,"_gedi_wk_", gediwk, "_Extracted",tile_id,".gpkg", sep = ""))
   if(file.exists(filepath)){
      print("File already exists",sep="")
      } else {
   cat("Reading in no. ", tile, "csv of ", length(all_gedil2_f), "csvs for iso3", iso3, "\n")

    ### Make this it's own function
    # Read GEDI L4A data
    gedil4_f_path <- paste(gedipath, "WDPA_gedi_L4A_tiles/",iso3,"/", all_gedil4_f[tile], sep = "")
    gedil4_f <- st_read(gedil4_f_path, int64_as_string = TRUE,drivers="GPKG")
    
    # Read GEDI L2A data
    gedil2_f_path <- paste(gedipath, "WDPA_gedi_L2A_tiles/",iso3,"/", all_gedil2_f[tile], sep = "")
    gedil2_f <- st_read(gedil2_f_path, int64_as_string = TRUE,drivers="GPKG")
    
    # Check if GEDI L4A data is empty
    if (nrow(gedil4_f) < 1) {
      cat("Error: No data for GEDI L4A\n")
      gedi_l24 <- gedil2_f
      gedi_l24$agbd <- NA
      gedi_l24$agbd_se <- NA
      gedi_l24$agbd_t <- NA
      gedi_l24$agbd_t_se <- NA
    } else {
      # Select relevant columns from GEDI L4A
      ### Drop the geometry, it's redundant
      gedi_l4_sub <- gedil4_f %>% st_drop_geometry() %>%
        dplyr::select(shot_number, agbd, agbd_se, agbd_t, agbd_t_se)
      ### Return here, do the join in the outer function
      # Join with GEDI L2A data
      gedi_l24 <- inner_join(gedil2_f, gedi_l4_sub, by = "shot_number")
    }
    # print(list(gedi_l24$shot_number))
    
    ###TODO: Make this it's own function
    # Read GEDI L2B data
    gedil2b_f_path <- paste(gedipath, "WDPA_gedi_L2B_tiles/",iso3,"/", all_gedil2b_f[tile], sep = "")
    gedil2b_f <- st_read(gedil2b_f_path, int64_as_string = TRUE,drivers="GPKG")
    names(gedil2b_f)[names(gedil2b_f) == "geolocation.lon_lowestmode"] <- "lon_lowestmode"
    names(gedil2b_f)[names(gedil2b_f) == "geolocation.lat_lowestmode"] <- "lat_lowestmode"
    names(gedil2b_f)[names(gedil2b_f) == "land_cover_data.landsat_treecover"] <- "landsat_treecover"

    variables <- c()

    # Loop from 0 to 29
    for (n in 0:29) {
      # Append the desired strings to the vector
      variables <- c(variables, paste("cover_z", n, sep=""))
      variables <- c(variables, paste("pai_z", n, sep=""))
      variables <- c(variables, paste("pavd_z", n, sep=""))
    }
        
    # Check if GEDI L2B data is empty
    if (nrow(gedil2b_f) < 1) {
      cat("Error: No data for GEDI L4A\n")
      gedi_l24b <- gedi_l24
      gedi_l24b <- gedi_l24
      gedi_l24b$landsat_treecover<- NA
      gedi_l24b$pai <- NA
      gedi_l24b$fhd_normal <- NA
    } else {
      # Select relevant columns from GEDI L4A
      ### Drop the geometry, it's redundant
      gedi_l2b_sub <- gedil2b_f %>% st_drop_geometry() %>%
        dplyr::select(shot_number, landsat_treecover, pai, fhd_normal,variables)
      ### Return here, do the join in the outer function
      # Join with GEDI L2A data
      gedi_l24b <- inner_join(gedi_l24, gedi_l2b_sub, by = "shot_number")
    }
    print(dim(gedi_l2b_sub))
    print(head(gedi_l24b))
    # print(colnames(gedi_l24b))
    # names(gedi_l24b)[names(gedi_l24b) == "shot_number"] <- "shotnum"
    # names(gedi_l24b)[names(gedi_l24b) == "lon_lowestmode"] <- "lonlow"
    # names(gedi_l24b)[names(gedi_l24b) == "lat_lowestmode"] <- "latlow"
    # names(gedi_l24b)[names(gedi_l24b) == "landsat_treecover"] <- "landtree"
    # names(gedi_l24b)[names(gedi_l24b) == "geolocation.sensitivity_a2"] <- "geosens"
    # st_write(gedi_l24b, dsn = paste(f.path3, iso3,"_extractStep1/", iso3, 
    #                                        "_gedi_wk_", gediwk, "_Extracted",tile_id,".gpkg", sep = ""))
    st_write(gedi_l24b, dsn = filepath)
    cat(tile_id, "in", iso3, "results are written to directory\n",filepath)
         
    }
    return(filepath)
}       

#Function to access STAC, this is used in the extract_gediPart2 function                                              
stac_to_terra <- function(catalog_url, ...) {
            # fetch STAC items
            stac <- rstac::stac(catalog_url)
            stac_items <- stac |>
              rstac::stac_search(
                ...
              ) |>
              rstac::get_request()
            
            # replace s3:// prefixes with /vsis3
            stac_items$features <- purrr::map(
                stac_items$features,
                ~ {
                    .x$assets <- purrr::map(
                        .x$assets,
                        ~{
                            .x$href <- gsub("s3://", "/vsis3/", .x$href)
                            .x
                        }
                    )
                .x
                }
            )
        
            item_collection_json_file <- tempfile(fileext=".json")
            item_collection_json <- jsonlite::toJSON(stac_items, auto_unbox=TRUE, pretty=TRUE, digits=10)
            print(item_collection_json)
        
            print(paste("writing item collection to", item_collection_json_file))
            write(item_collection_json, item_collection_json_file)
            
            item_collection_dsn <- glue::glue(
                "STACIT:\"{item_collection_json}\":asset=data",
                item_collection_json=item_collection_json_file
            )
                
            terra::rast(item_collection_dsn)
            }
                                               
extract_gediPart2 <- function(matched,mras,extracted,catalog_url){
    iso_matched_gedi_df <- NULL
    results_list <- list()
    # Initialize empty spatial object for the current iteration
    for (this_csvid in seq_along(extracted)) {
        spatial_data <- vect(extracted[this_csvid])

        extent <- sf::st_bbox(spatial_data)

        # load the rasters
            glad_change_rast <- stac_to_terra(
                catalog_url = catalog_url,
                bbox = extent,
                collections = "glad-glclu2020-change-v2",
                datetime = "2020-01-01T00:00:00Z",
            )
            
            glad_rast_2020 <- stac_to_terra(
                catalog_url = catalog_url,
                bbox = extent,
                collections = "glad-glclu2020-v2",
                datetime = "2020-01-01T00:00:00Z",
            )
        
          # extent <- terra::ext(spatial_data)
          # selection <- terra::crop(glad_rast_change, extent)
         spatial_data$glad_change <- terra::extract(glad_change_rast,spatial_data, df=TRUE)[[2]]
         spatial_data$glad_2020 <- terra::extract(glad_rast_2020,spatial_data, df=TRUE)[[2]]

          gedi_l24b_sp <- project(spatial_data, "epsg:6933")
          matched_gedi <- terra::extract(mras,gedi_l24b_sp, df=TRUE)
          matched_gedi_metrics <- cbind(matched_gedi,gedi_l24b_sp)
          # matched_gedi_metrics <- inner_join(matched_gedi_metrics, spatial_data, by = "ID")
          # matched_gedi_metrics <- inner_join(matched_gedi_metrics, matched_glad, by = "ID")
          print(head(matched_gedi_metrics))
          matched_gedi_metrics_filtered <- matched_gedi_metrics %>% dplyr::filter(!is.na(status)) %>% 
          convertFactor(matched0 = matched,exgedi = .) 
          print(head(matched_gedi_metrics_filtered))
          
          iso_matched_gedi_df <- rbind(matched_gedi_metrics_filtered,iso_matched_gedi_df)
          print(dim(iso_matched_gedi_df))
        # }
    
    # Store results in a list
    results_list[[this_csvid]] <- iso_matched_gedi_df
      }
  
      # Combine all results
      if (!is.null(iso_matched_gedi_df)) {
        iso_matched_gedi_df <- do.call(rbind, results_list)
      }
      
      cat("Done GEDI processing\n")
      return(iso_matched_gedi_df)
}                                            

############################################################################################   
                                      

SplitRas <- function(raster,ppside){
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- aggregate(raster,fact=c(h,v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:ncell(agg)){
    e1          <- ext(agg_poly[agg_poly$polis==i,])
    r_list[[i]] <- crop(raster,e1)
  }
  return(r_list)
}

rasExtract2020 <- function(l4_sp){
  # cat(iso3,"converting the matched csv to a raster stack for extraction\n")
  tif2020 <- c("pop_cnt_2020","pop_den_2020","lc2019","tt2cities_2015","wc_prec_2010-2018","wc_tavg_2010-2018","wc_tmax_2010-2018",
               "wc_tmin_2010-2018","wwf_biomes","wwf_ecoreg","dem","slope","d2roads","dcities")
  for (t in 1:length(tif2020)){
    # print(tif2020[t])
    covar2020 <- raster(paste(f.path, "WDPA_input_vars_iso3_v2/",iso3,"/",tif2020[t],".tif", sep=""))
    ras_ex <- raster::extract(covar2020, l4_sp@coords, method="simple", factors=F)
    nm <- names(covar2020)
    l4_sp <- cbind(l4_sp, ras_ex)
    names(l4_sp)[t+6] <- tif2020[t]
  }
  return(l4_sp)
}

# iso_matched_gedi <- foreach(this_csv=gedil2_f, .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {
#   ##add the GEDI l4a model prediction for AGB here :
#   cat("Readng in no. ", match(this_csv, gedil2_f),"csv of ", length(gedil2_f),"csvs for iso3",iso3,"\n")
#   gedi_l2  <- read.csv(paste(f.path,"WDPA_gedi_l2a+l2b_clean",iso3,this_csv, sep="/")) %>%
#     dplyr::select(shot_number,lon_lowestmode, lat_lowestmode, starts_with("rh_"),cover, pai)%>%
#     SpatialPointsDataFrame(coords=.[,c("lon_lowestmode","lat_lowestmode")],
#                            proj4string=CRS("+init=epsg:4326"), data=.) %>%spTransform(., CRS("+init=epsg:6933"))
#   
#   iso_matched_gedi_df <- data.frame()
#   matched_gedi <- raster::extract(mras,gedi_l2, df=TRUE)
#   matched_gedi_metrics <- cbind(matched_gedi,gedi_l2@data)
#   
#   matched_gedi_metrics_filtered <- matched_gedi_metrics %>% dplyr::filter(!is.na(status)) %>% 
#     convertFactor(matched0 = matched,exgedi = .) 
#   
#   matched_gedi_l4a <-matched_gedi_metrics_filtered %>% 
#     dplyr::mutate(
#       LAT=lat_lowestmode,
#       LON=lon_lowestmode,
#       REGION=region,
#       PFT=pft,
#       RH_10=rh_010+100,
#       RH_20=rh_020+100,
#       RH_30=rh_030+100,
#       RH_40=rh_040+100,
#       RH_50=rh_050+100,
#       RH_60=rh_060+100,
#       RH_70=rh_070+100,
#       RH_80=rh_080+100,
#       RH_90=rh_090+100,
#       RH_98=rh_098+100) %>% 
#     modelr::add_predictions(model2, "AGBD")
#   iso_matched_gedi_df <- rbind(matched_gedi_l4a,iso_matched_gedi_df)
#   return(iso_matched_gedi_df)
# }
# stopImplicitCluster()
# 