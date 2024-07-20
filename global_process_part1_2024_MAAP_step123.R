#!/usr/bin/env Rscript

# This global processing script is derived from the global processing notebook 
#the input can be the iso3 code (3-character) for one or multiple countries 

options(warn=-1)
options(dplyr.summarise.inform = FALSE)

packages <- c("sp","rgdal","sf","rgeos","dplyr","plyr","ggplot2","raster","mapview","stringr",
              "maptools","gridExtra","lattice","MASS","foreach","optmatch","doParallel",
              "rlang","tidyr","magrittr","viridis","ggmap","spatialEco","bit64",
              "randomForest","modelr","geojsonio","rgeos") #"hrbrthemes","RItools","Hmisc",
package.check <- lapply(packages, FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
})

#To test, we define the variables manually. For final version, run the commented out section below
#gediwk <- 24
#iso3 <-"ECU"
#mproc <- 1
#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  
  iso3 <- args[1]  #country to process
  gediwk <- args[2]   #the # of weeks GEDI data to use
  mproc <- as.integer(args[3])  #the number of cores to use for matching 
}
#-------------------------------------------------------------------------------

cat("Step 0: Loading global variables for", iso3,"with wk", gediwk, "data \n")

f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"

matching_tifs <- c("wwf_biomes","wwf_ecoreg","lc2000","d2roads", "dcities","dem",
                   "pop_cnt_2000","pop_den_2000","slope", "tt2cities_2000", "wc_prec_1990-1999",
                   "wc_tmax_1990-1999","wc_tavg_1990-1999","wc_tmin_1990-1999" )

ecoreg_key <- read.csv(paste(f.path,"wwf_ecoregions_key.csv",sep=""))
allPAs <- readRDS(paste(f.path,"WDPA_shapefiles/WDPA_polygons/",iso3,"_PA_poly.rds",sep=""))
MCD12Q1 <- raster(paste(f.path,"GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v1_projection_defined_6933.tif",sep=""))
projection(MCD12Q1) <- sp::CRS(paste("+init=epsg:",6933,sep=""))
world_region <- raster(paste(f.path,"GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v1_revised_projection_defined_6933.tif",sep=""))
projection(world_region) <- sp::CRS(paste("+init=epsg:",6933,sep=""))
adm <- readOGR(paste(f.path,"WDPA_countries/shp/",iso3,".shp",sep=""),verbose=F)
adm_prj <- spTransform(adm, "+init=epsg:6933") 
load(paste(f.path,"rf_noclimate.RData",sep=""))
source(paste(f.path,"matching_func.R",sep=""))
source(paste(f.path,"matching_func_2024.R",sep=""))

# STEP1. Create 1km sampling grid with points only where GEDI data is available; first check if grid file exist to avoid reprocessing 
if(!file.exists(paste(f.path,"WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep=""))){
  cat("Step 1: Creating 1km sampling grid filter GEDI data for", iso3,"\n")
  GRID.lats <- raster(file.path(f.path,"EASE2_M01km_lats.tif"))
  GRID.lons <- raster(file.path(f.path,"EASE2_M01km_lons.tif"))
  GRID.lats.adm   <- crop(GRID.lats, adm_prj)
  GRID.lats.adm.m <- raster::mask(GRID.lats.adm, adm_prj)
  GRID.lons.adm   <- crop(GRID.lons, adm_prj)
  GRID.lons.adm.m <- raster::mask(GRID.lons.adm, adm_prj)
  rm(GRID.lats, GRID.lons, GRID.lats.adm, GRID.lons.adm)
  
  #1.3) extract coordinates of raster cells with valid GEDI data in them
  #gedi_folder <- paste(f.path,"WDPA_gedi_l2a+l2b_clean2/",iso3,"/",sep="")
  gedi_folder <- paste(f.path,"WDPA_gedi_L4A_tiles/",sep="")
  #iso3_tiles <- paste("/projects/my-public-bucket/AOIs/vero_1deg_tiles_",iso3,"/",sep="")
  tileindex_df <- read.csv(paste("/projects/my-public-bucket/AOIs/vero_1deg_tileindex/tileindex_",iso3,".csv", sep=""))
  iso3_tiles <- tileindex_df$tileindex
  
  GRID.coords <- data.frame()
  #for(i in 1:length(dir(gedi_folder))){
  #for(i in 1:length(dir(iso3_tiles))){
  for(i in 1:length(iso3_tiles)){
    # print(list.files(gedi_folder)[i])
    
    iso3_tile_in <- paste("tile_num_",iso3_tiles[i],sep="")
    
    if(!file.exists(paste(gedi_folder,iso3_tile_in,"_L4A.gpkg",sep=""))){
      print(paste(iso3_tile_in," does not exist",sep=""))
    } else {
      #print(paste(iso3_tile_in," processing",sep=""))
      #gedi_data <- read.csv(list.files(gedi_folder,full.names=TRUE)[i]) %>%
      gedi_data <- read_sf(paste(gedi_folder,iso3_tile_in,"_L4A.gpkg",sep=""), int64_as_string = TRUE) %>%
        dplyr::select(lon_lowestmode,lat_lowestmode)
      gedi_data <- gedi_data %>% st_drop_geometry()
      gedi_pts  <- SpatialPoints(coords=gedi_data[,c("lon_lowestmode","lat_lowestmode")],
                                 proj4string=CRS("+init=epsg:4326"))
      gedi_pts_prj <- spTransform(gedi_pts, "+init=epsg:6933")
      
      gcount_ras <- rasterize(coordinates(gedi_pts_prj),GRID.lons.adm.m , fun="count",background=NA)
      names(gcount_ras) <- "gshot_counts"
      pxid <- raster::extract(gcount_ras,  gedi_pts_prj)
      gedi_pts_prj %>% 
        SpatialPointsDataFrame(., data=data.frame(pxid)) ->gedi_pts_prj_sp
      gedi_pts_prj_sp$pxid[is.na(gedi_pts_prj_sp$pxid)] <- 0
      gedi_pts_prj_sp[gedi_pts_prj_sp$pxid>5,]->gedi_pts_prj_filtered  #change the numeric threshold to filter with a different min # of GEDI shots in each 1km cell
      
      GRID.lons.overlap <- GRID.lons.adm.m[gedi_pts_prj_filtered]
      GRID.lats.overlap <- GRID.lats.adm.m[gedi_pts_prj_filtered]
      
      x.overlap <- GRID.lons.overlap[!is.na(GRID.lons.overlap)]
      y.overlap <- GRID.lats.overlap[!is.na(GRID.lats.overlap)]
      
      xy.overlap <- cbind(x.overlap,y.overlap)
      xy.overlap.clean <- unique(xy.overlap)
      
      GRID.coords <- rbind(GRID.coords, xy.overlap.clean)
    }
  }
  GRID.for.matching <- SpatialPoints(coords = GRID.coords, proj4string=CRS("+init=epsg:4326"))
  saveRDS(GRID.for.matching, file = paste(f.path,"WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep=""))
} else if (file.exists(paste(f.path,"WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep=""))) {
  cat(paste("STEP 1: Grid file exists, no need to process grids for ",iso3, "\n"))
}

#plot(GRID.for.matching, pch=".", cex=0.25)
  
# STEP2. Clip sampling grid to nonPA areas within country & sample raster layers on nonPA grid
cat("Step 2.0: Reading 1k GRID from RDS for " ,iso3, "\n")
GRID.for.matching <- readRDS(paste(f.path,"WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep="")) 
#print(GRID.for.matching)

if(!file.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_prepped_control_wk",gediwk,".RDS",sep=""))){
  if(!dir.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",sep=""))){
    dir.create(paste(f.path,"WDPA_matching_points/",iso3,"/",sep=""))}    
  cat("Step 2.1: Preparing control dataset for", iso3, "\n")
  GRID.pts.nonPA <- GRID.for.matching %>% spTransform(., "+init=epsg:4326")
  for(i in 1:length(allPAs)){
    PA          <- allPAs[i,]
    PA_prj      <- spTransform(PA, "+init=epsg:6933")
    PA_prj_buff <- gBuffer(PA_prj, width = 10000) #10km buffer
    PA2         <- spTransform(PA_prj_buff, "+init=epsg:4326")
    overlap     <- GRID.pts.nonPA[PA2]
    if(length(overlap)>0){
      GRID.pts.nonPA0 <- st_difference(sf::st_as_sf(GRID.pts.nonPA), sf::st_as_sf(PA2)) ##remove pts inside poly
      GRID.pts.nonPA <- as(GRID.pts.nonPA0$geometry,'Spatial') %>% spTransform(., "+init=epsg:4326")
    } 
    # print(length(GRID.pts.nonPA))
  }
  nonPA_xy  <- coordinates(GRID.pts.nonPA)
  colnames(nonPA_xy)  <- c("x","y")
  nonPA_spdf  <- tryCatch(SpatialPointsDataFrame(nonPA_xy, data=data.frame(nonPA_xy),
                                                 proj4string=CRS("+init=epsg:4326")),
                          error=function(cond){
                            cat("Country too samll, after buffer no grid left, so quit processing country", iso3, dim(nonPA_xy),"\n")
                            writeLines("Country too samll, after buffer no grid left", paste(f.path,"WDPA_log/",iso3,"_log_control.txt", sep=""))
                            return(quit(save="no"))})
  
  for (j in 1:length(matching_tifs)){
    #ras <- raster(paste(f.path, "WDPA_input_vars_iso3/",iso3,"/",matching_tifs[j],".tif", sep=""))
    ras <- raster(paste(f.path, "WDPA_input_vars_GLOBAL/",matching_tifs[j],".tif", sep=""))
    print(matching_tifs[j])
    ras_ex <- raster::extract(ras, nonPA_spdf@coords, method="simple", factors=FALSE)
    nm <- names(ras)
    nonPA_spdf <- cbind(nonPA_spdf, ras_ex)
    names(nonPA_spdf)[j+2] <- matching_tifs[j]
    
  }
  d_control <- nonPA_spdf
  d_control$status <- as.logical("FALSE")
  names(d_control) <- make.names(names(d_control), allow_ = FALSE)
  d_control <- data.frame(d_control) %>%
    dplyr::rename(
      land_cover = lc2000,
      slope = slope,
      elevation = dem,
      popden = pop.den.2000,
      popcnt=pop.cnt.2000,
      min_temp=wc.tmin.1990.1999,
      max_temp=wc.tmax.1990.1999,
      mean_temp = wc.tavg.1990.1999,
      prec = wc.prec.1990.1999,
      tt2city= tt2cities.2000,
      wwfbiom = wwf.biomes,
      wwfecoreg = wwf.ecoreg,
      d2city = dcities,
      d2road = d2roads,
      lon = x,
      lat = y)
  d_control$land_cover <- factor(d_control$land_cover, levels=sequence(7),
                                 labels = c("l1_forest",
                                            "l2_grassland",
                                            "l3_agriculture",
                                            "l4_wetlands",
                                            "l5_artificial",
                                            "l6_other land/bare",
                                            "l7_water"))
  d_control$wwfbiom <- factor(d_control$wwfbiom,
                              levels = as.vector(unique(ecoreg_key[,"BIOME"])),
                              labels = as.vector(unique(ecoreg_key[,"BIOME_NAME"])))
  d_control$wwfecoreg <- factor(d_control$wwfecoreg,
                                levels = as.vector(ecoreg_key[,"ECO_ID"]),
                                labels = as.vector(ecoreg_key[,"ECO_NAME"]))
  
  
  d_control$UID <-  seq.int(nrow(d_control))
  
  saveRDS(d_control, file=paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_prepped_control_wk",gediwk,".RDS",sep="")) 
  
} else if (file.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_prepped_control_wk",gediwk,".RDS",sep=""))){
  cat("Step 2.1: preppred control dataset already exists for", iso3, "no need for reprocessing\n")
}


#STEP3. Loop through all PAs in iso3 country:
# - clip sampling grid to each PA
# - sample raster layers on each PA grid
# - save each PA sample into prepped_pa_##.RDS file

cat("Step 3.0: Reading 1k GRID from RDS for " ,iso3, "\n")
GRID.for.matching <- readRDS(paste(f.path,"WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep="")) 

if(length(dir(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""),pattern = paste(gediwk,".RDS",sep="")))==0){
  if(!dir.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""))){
    dir.create(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""))}
  cat("Step 3.1: Processing prepped PA treatment dataset for ", iso3, "\n")
  for(i in 1:length(allPAs)){
    cat(iso3, i, "out of ", length(allPAs), "\n")
    testPA <- allPAs[i,]
    testPA <- spTransform(testPA, "+init=epsg:4326")
    GRID.pts.testPA <- GRID.for.matching[testPA]
    
    if(length(GRID.pts.testPA)>0){
      testPA_xy <- coordinates(GRID.pts.testPA)
      colnames(testPA_xy) <- c("x","y")
      testPA_spdf <- SpatialPointsDataFrame(testPA_xy, data=data.frame(testPA_xy),
                                            proj4string=CRS("+init=epsg:4326"))
      for (j in 1:length(matching_tifs)){
        #ras <- raster(paste(f.path, "WDPA_input_vars_iso3/",iso3,"/",matching_tifs[j],".tif", sep=""))
        ras <- raster(paste(f.path, "WDPA_input_vars_GLOBAL/",matching_tifs[j],".tif", sep=""))
        ras <- crop(ras, testPA)
        ras_ex <- raster::extract(ras, testPA_spdf@coords, method="simple", factors=F)
        nm <- names(ras)
        testPA_spdf <- cbind(testPA_spdf, ras_ex)
        names(testPA_spdf)[j+2] <- matching_tifs[j]
        
      }
      d_pa <- testPA_spdf
      d_pa$status <- as.logical("TRUE")
      d_pa$DESIG_ENG <- testPA$DESIG_ENG
      d_pa$REP_AREA <- testPA$REP_AREA
      d_pa$PA_STATUS <- testPA$STATUS
      d_pa$PA_STATUSYR <- testPA$STATUS_YR
      d_pa$GOV_TYPE <- testPA$GOV_TYPE
      d_pa$OWN_TYPE <- testPA$OWN_TYPE
      d_pa$MANG_AUTH <- testPA$MANG_AUTH
      names(d_pa) <- make.names(names(d_pa), allow_ = FALSE)
      d_pa <- data.frame(d_pa) %>%
        dplyr::rename(
          land_cover = lc2000,
          slope = slope,
          elevation = dem,
          popden = pop.den.2000,
          popcnt=pop.cnt.2000,
          min_temp=wc.tmin.1990.1999,
          max_temp=wc.tmax.1990.1999,
          mean_temp = wc.tavg.1990.1999,
          prec = wc.prec.1990.1999,
          tt2city= tt2cities.2000,
          wwfbiom = wwf.biomes,
          wwfecoreg = wwf.ecoreg,
          d2city = dcities,
          d2road = d2roads,
          lon = x,
          lat = y)
      d_pa$land_cover <- factor(d_pa$land_cover, levels=sequence(7),
                                labels = c("l1_forest",
                                           "l2_grassland",
                                           "l3_agriculture",
                                           "l4_wetlands",
                                           "l5_artificial",
                                           "l6_other land/bare",
                                           "l7_water"))
      d_pa$wwfbiom <- factor(d_pa$wwfbiom,
                             levels = as.vector(unique(ecoreg_key[,"BIOME"])),
                             labels = as.vector(unique(ecoreg_key[,"BIOME_NAME"])))
      d_pa$wwfecoreg <- factor(d_pa$wwfecoreg,
                               levels = as.vector(ecoreg_key[,"ECO_ID"]),
                               labels = as.vector(ecoreg_key[,"ECO_NAME"]))
      
      d_pa$UID <- seq.int(nrow(d_pa))
      saveRDS(d_pa, file = paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/","prepped_pa_",
                                 testPA$WDPAID,"_wk",gediwk,".RDS", sep="")) 
    }
  }
} else if (length(dir(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""),pattern = paste(gediwk,".RDS",sep="")))>0){
  cat("Step 3.1: prepped PA treatment dataset already exists for ", iso3, "no need for reprocessing\n")
}
