#!/usr/bin/env Rscript

# This global processing script is derived from the global processing notebook 
#the input can be the iso3 code (3-character) for one or multiple countries 

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# List of CRAN packages to be installed
cran_packages <- c(
  "s3"
)

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

library("terra")
library("dplyr")
library("sf")
#install.packages("s3")
library("s3")

#To test, we define the variables manually. For final version, run the commented out section below
#iso3 <-"ECU"
gediwk <- 24

#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  
  iso3 <- args[1]  #country to process
  out <- args[2]
  #gediwk <- args[2]   #the # of weeks GEDI data to use
  #mproc <- as.integer(args[3])  #the number of cores to use for matching
}
#-------------------------------------------------------------------------------

cat("Step 0: Loading global variables for", iso3,"with wk", gediwk, "data \n")

#f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
f.path2 <- "s3://maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/"
gedipath<- "/vsis3/maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/" #Make sure to specify username
veropath<- "/vsis3/maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/" #Make sure to specify username
f.path3<- file.path(out)

adm <- st_read(paste(sub("s3://","/vsis3/", f.path),"WDPA_countries/shp/",iso3,".shp",sep=""))

adm_prj <- project(vect(adm), "epsg:6933")

vect_adm <- vect(adm)

extent <- sf::st_bbox(vect_adm)

source("matching_func_2024.r")

glad <- paste0(iso3,"_GLADCover_reclass2000")
gladEntry <- paste0(iso3,".GLADCover.reclass2000")

gmw <- paste0(iso3,"_gmw1996")
gmwEntry <- paste0(iso3,".gmw1996")

matching_tifs <- c("d2roads", "dcities","dem",
                   "pop_cnt_2000","pop_den_2000","slope", "tt2cities_2000", "wc_prec_1990-1999",
                   "wc_tmax_1990-1999-Copy1","wc_tavg_1990-1999","wc_tmin_1990-1999", glad, gmw)

ecoreg_key <- read.csv(s3_get(paste(f.path,"wwf_ecoregions_key.csv",sep="")))
#unlink(s3_get(paste(f.path,"wwf_ecoregions_key.csv",sep="")))

allPAs <- readRDS(s3_get(paste(f.path,"WDPA_shapefiles/WDPA_polygons/",iso3,"_PA_poly.rds",sep="")))

MCD12Q1 <- rast(paste(veropath,"GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v1_projection_defined_6933.tif",sep=""))
crs(MCD12Q1)  <- "epsg:6933"

world_region <- rast(paste(veropath,"GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v1_revised_projection_defined_6933.tif",sep=""))
crs(world_region)  <- "epsg:6933"

s3_get_files(c(paste(f.path,"WDPA_countries/shp/",iso3,".shp",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".shx",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".prj",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".dbf",sep="")),confirm = FALSE)

load(s3_get(paste(f.path,"rf_noclimate.RData",sep="")))



# STEP1. Create 1km sampling grid with points only where GEDI data is available; first check if grid file exist to avoid reprocessing 
#if(!file.exists(paste(f.path,"WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep=""))){
  cat("Step 1: Creating 1km sampling grid filter GEDI data for", iso3,"\n")
  GRID.lats <- rast(paste(veropath,"EASE2_M01km_lats.tif", sep=""))
  GRID.lons <- rast(paste(veropath,"EASE2_M01km_lons.tif", sep=""))
  GRID.lats.adm   <- crop(GRID.lats, adm_prj)
  GRID.lats.adm.m <- mask(GRID.lats.adm, adm_prj)
  GRID.lons.adm   <- crop(GRID.lons, adm_prj)
  GRID.lons.adm.m <- mask(GRID.lons.adm, adm_prj)
  rm(GRID.lats, GRID.lons, GRID.lats.adm, GRID.lons.adm)

  #1.3) extract coordinates of raster cells with valid GEDI data in them
  gedi_folder <- paste(gedipath,"WDPA_gedi_L4A_tiles/",iso3,"/",sep="")
# gedi_folder <- paste("~/my-public-bucket/GEDI_global_PA_v2/WDPA_gedi_L4A_tiles/",sep="")
  tileindex_df <- read.csv(s3_get(paste(f.path,"vero_1deg_tileindex/tileindex_",iso3,".csv", sep="")))
  iso3_tiles <- tileindex_df$tileindexiso3_tiles <- tileindex_df$tileindex
    
  GRID.coords <- data.frame()
  for(i in 1:length(iso3_tiles)){
    
    iso3_tile_in <- paste("tile_num_",iso3_tiles[i],sep="")
      

    print(paste(iso3_tile_in," processing",sep=""))
    print(paste(gedi_folder,iso3_tile_in,"_L4A.gpkg",sep=""))
    #if(!file.exists(paste(gedi_folder,iso3_tile_in,"_L4A.gpkg",sep=""))){
    #    print(paste(iso3_tile_in," does not exist",sep=""))
    #    } else {
    geopath<- paste0(gedi_folder,iso3_tile_in,"_L4A.gpkg")
    gedi_data <- read_sf(geopath,layer=paste0(iso3_tile_in,"_L4A")) %>% #NO S3 get here, if spatial format, don't use S3 lib
      dplyr::select(lon_lowestmode,lat_lowestmode)
    gedi_data <- gedi_data %>% st_drop_geometry()
    gedi_pts  <- vect(gedi_data, geom=c("lon_lowestmode","lat_lowestmode"), crs="epsg:4326", keepgeom=FALSE)        
    gedi_pts_prj <- project(gedi_pts, "epsg:6933")
        
    gcount_ras <- rasterize(geom(gedi_pts_prj)[,c("x","y")], GRID.lons.adm.m, fun="count", background=NA)
    names(gcount_ras) <- "gshot_counts"
    pxid <- extract(gcount_ras, gedi_pts_prj)
    gedi_pts_prj$pxid <- pxid[,"gshot_counts"]
    gedi_pts_prj_sp <- gedi_pts_prj    
    gedi_pts_prj_sp$pxid[is.na(gedi_pts_prj_sp$pxid)] <- 0
    gedi_pts_prj_filtered <- gedi_pts_prj_sp[gedi_pts_prj_sp$pxid >= 1,]  #change the numeric threshold to filter with a different min # of GEDI shots in each 1km cell
    
    GRID.lons.overlap <- GRID.lons.adm.m[gedi_pts_prj_filtered]
    GRID.lats.overlap <- GRID.lats.adm.m[gedi_pts_prj_filtered]
    
    x.overlap <- GRID.lons.overlap[!is.na(GRID.lons.overlap)]
    y.overlap <- GRID.lats.overlap[!is.na(GRID.lats.overlap)]
    
    xy.overlap <- cbind(x.overlap,y.overlap)
    xy.overlap.clean <- unique(xy.overlap)
    
    GRID.coords <- rbind(GRID.coords, xy.overlap.clean)
    #}
  }
  #GRID.for.matching <- SpatialPoints(coords = GRID.coords, proj4string=CRS("+init=epsg:4326"))
  GRID.for.matching <- vect(GRID.coords, geom=c("x.overlap","y.overlap"), crs = "epsg:4326")
  dir.create(file.path(f.path3, "WDPA_grids"), recursive = TRUE, showWarnings = FALSE)
  filename_out <- paste(f.path3, "/WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep="")
  print(filename_out)
  # s3saveRDS(x=GRID.for.matching, bucket = "s3://maap-ops-workspace/", object = filename_out, region = "us-west-2")
  saveRDS(GRID.for.matching, filename_out)

#} else if (file.exists(paste(f.path,"WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep=""))) {
#  cat(paste("STEP 1: Grid file exists, no need to process grids for ",iso3, "\n"))
#}


# STEP2. Clip sampling grid to nonPA areas within country & sample raster layers on nonPA grid
cat("Step 2.0: Reading 1k GRID from RDS for " ,iso3, "\n")
#GRID.for.matching <- vect(GRID.for.matching)

#if(!file.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_prepped_control_wk",gediwk,".RDS",sep=""))){
#  if(!dir.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",sep=""))){
#      dir.create(paste(f.path,"WDPA_matching_points/",iso3,"/",sep=""))}    
  cat("Step 2.1: Preparing control dataset for", iso3, "\n")
  GRID.pts.nonPA <- project(GRID.for.matching, "epsg:4326") #***Ask Vero why we're doing this 1 at a time, rather than all at once maybe vectorize
  for(i in 1:length(allPAs)){
    PA          <- vect(allPAs[i,])
    PA_prj      <- project(PA, "epsg:6933") #**Why do we have to change the projection twice?
    PA_prj_buff <- buffer(PA_prj, width = 10000) ##10km buffer
    PA2         <- project(PA_prj_buff, "epsg:4326") 
    overlap     <- GRID.pts.nonPA[PA2]
    if(length(overlap)>0){
      GRID.pts.nonPA0 <- st_difference(sf::st_as_sf(GRID.pts.nonPA), sf::st_as_sf(PA2)) ##remove pts inside poly
      GRID.pts.nonPA <- vect(GRID.pts.nonPA0$geometry)
      GRID.pts.nonPA <- project(GRID.pts.nonPA, "epsg:4326")
    } 
    # print(length(GRID.pts.nonPA))
  }
  nonPA_xy  <- geom(GRID.pts.nonPA)[,c("x","y")]
  colnames(nonPA_xy)  <- c("x","y")
  nonPA_spdf  <- tryCatch(vect(nonPA_xy, crs="epsg:4326"),      
                          error=function(cond){
                            cat("Country too small - quit processing ", iso3, dim(nonPA_xy),"\n")
                            return(quit(save="no"))})
    
    
  for (j in 1:length(matching_tifs)){
    ras <- rast(paste(gedipath, "WDPA_input_vars_GLOBAL/",matching_tifs[j],".tif", sep=""))
    print(matching_tifs[j])
    ras_ex <- extract(ras, nonPA_spdf, method="simple", factors=FALSE)
    nm <- names(ras)
    nonPA_spdf$nm <- ras_ex[, matching_tifs[j]]
    names(nonPA_spdf)[j] <- matching_tifs[j]
  }
nonPA_spdf$x <- geom(nonPA_spdf)[,"x"]
  nonPA_spdf$y <- geom(nonPA_spdf)[,"y"]
  
  d_control <- nonPA_spdf
  d_control$status <- as.logical("FALSE")
  names(d_control) <- make.names(names(d_control), allow_ = FALSE)
# STEP2. Clip sampling grid to nonPA areas within country & sample raster layers on nonPA grid
cat("Step 2.0: Reading 1k GRID from RDS for " ,iso3, "\n")
#GRID.for.matching <- vect(GRID.for.matching)

#if(!file.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_prepped_control_wk",gediwk,".RDS",sep=""))){
#  if(!dir.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",sep=""))){
#      dir.create(paste(f.path,"WDPA_matching_points/",iso3,"/",sep=""))}    
  cat("Step 2.1: Preparing control dataset for", iso3, "\n")
  GRID.pts.nonPA <- project(GRID.for.matching, "epsg:4326") #***Ask Vero why we're doing this 1 at a time, rather than all at once maybe vectorize
  for(i in 1:length(allPAs)){
    PA          <- vect(allPAs[i,])
    PA_prj      <- project(PA, "epsg:6933") #**Why do we have to change the projection twice?
    PA_prj_buff <- buffer(PA_prj, width = 10000) ##10km buffer
    PA2         <- project(PA_prj_buff, "epsg:4326") 
    overlap     <- GRID.pts.nonPA[PA2]
    if(length(overlap)>0){
      GRID.pts.nonPA0 <- st_difference(sf::st_as_sf(GRID.pts.nonPA), sf::st_as_sf(PA2)) ##remove pts inside poly
      GRID.pts.nonPA <- vect(GRID.pts.nonPA0$geometry)
      GRID.pts.nonPA <- project(GRID.pts.nonPA, "epsg:4326")
    } 
    # print(length(GRID.pts.nonPA))
  }
  nonPA_xy  <- geom(GRID.pts.nonPA)[,c("x","y")]
  colnames(nonPA_xy)  <- c("x","y")
  nonPA_spdf  <- tryCatch(vect(nonPA_xy, crs="epsg:4326"),      
                          error=function(cond){
                            cat("Country too small - quit processing ", iso3, dim(nonPA_xy),"\n")
                            return(quit(save="no"))})
    
    
  for (j in 1:length(matching_tifs)){
    ras <- rast(paste(gedipath, "WDPA_input_vars_GLOBAL/",matching_tifs[j],".tif", sep=""))
    print(matching_tifs[j])
    ras_ex <- extract(ras, nonPA_spdf, method="simple", factors=FALSE)
    ras_ex[is.na(ras_ex)] <- 0
    nm <- names(ras)
    nonPA_spdf$nm <- ras_ex[, matching_tifs[j]]
    names(nonPA_spdf)[j] <- matching_tifs[j]
  }
  nonPA_spdf$x <- geom(nonPA_spdf)[,"x"]
  nonPA_spdf$y <- geom(nonPA_spdf)[,"y"]
  
  d_control <- nonPA_spdf
  d_control$status <- as.logical("FALSE")
  names(d_control) <- make.names(names(d_control), allow_ = FALSE)
  d_control <- data.frame(d_control) %>%
    dplyr::rename(
      land_cover = gladEntry,
      mangrove = gmwEntry,
      slope = slope,
      elevation = dem,
      popden = pop.den.2000,
      popcnt=pop.cnt.2000,
      min_temp=wc.tmin.1990.1999,
      max_temp=wc.tmax.1990.1999.Copy1,
      mean_temp = wc.tavg.1990.1999,
      prec = wc.prec.1990.1999,
      tt2city= tt2cities.2000,
      d2city = dcities,
      d2road = d2roads,
      lon = x,
      lat = y)
  d_control$land_cover <- factor(d_control$land_cover, levels=sequence(13),
                                 labels = c("desert","semi_arid","dense_short","tree_short","tree_med","tree_tall",
                                           "salt pan", "sparse_veg_wetland","dense_short_wetland", "tree_short_wetland",
                                           "tree_med_wetland","tree_tall_wetland","water"))  
  
  d_control$UID <-  seq.int(nrow(d_control))
  
  filename_out <- paste(f.path3, "/WDPA_grids/",iso3,"_prepped_control_wk",gediwk,".RDS", sep="")
  # s3saveRDS(x = d_control, bucket = "s3://maap-ops-workspace/", object = filename_out, region = "us-west-2")  
  saveRDS(d_control, filename_out)
#} else if (file.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_prepped_control_wk",gediwk,".RDS",sep=""))){
#  cat("Step 2.1: preppred control dataset already exists for", iso3, "no need for reprocessing\n")
#}


#STEP3. Loop through all PAs in iso3 country:
# - clip sampling grid to each PA
# - sample raster layers on each PA grid
# - save each PA sample into prepped_pa_##.RDS file

cat("Step 3.0: Reading 1k GRID from RDS for " ,iso3, "\n")
#GRID.for.matching <- vect(GRID.for.matching)

#if(length(dir(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""),pattern = paste(gediwk,".RDS",sep="")))==0){
#  if(!dir.exists(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""))){
#      dir.create(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""))}
  cat("Step 3.1: Processing prepped PA treatment dataset for ", iso3, "\n")
  for(i in 1:length(allPAs)){
    cat(iso3, i, "out of ", length(allPAs), "\n")
    testPA <- vect(allPAs[i,])
    testPA <- project(testPA, "epsg:4326")
    GRID.pts.testPA <- GRID.for.matching[testPA]
    
    #if(length(GRID.pts.testPA)>0){
    if(length(GRID.pts.testPA)>1){
      testPA_xy <- geom(GRID.pts.testPA)[,c("x","y")]
      colnames(testPA_xy) <- c("x","y")
      testPA_spdf  <- vect(testPA_xy, crs="epsg:4326")
                              
       for (j in 1:length(matching_tifs)){
        ras <- rast(paste(gedipath, "WDPA_input_vars_GLOBAL/",matching_tifs[j],".tif", sep=""))
        # ras[is.na(ras)] <- 0
        ras <- crop(ras, testPA)
        ras_ex <- extract(ras, testPA_spdf, method="simple", factors=F)
        nm <- names(ras)
        testPA_spdf$nm <- ras_ex[, matching_tifs[j]]
        names(testPA_spdf)[j] <- matching_tifs[j]
  }
    testPA_spdf$x <- geom(testPA_spdf)[,"x"]
    testPA_spdf$y <- geom(testPA_spdf)[,"y"]
      
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
              land_cover = gladEntry,
              mangrove = gmwEntry,
              slope = slope,
              elevation = dem,
              popden = pop.den.2000,
              popcnt=pop.cnt.2000,
              min_temp=wc.tmin.1990.1999,
              max_temp=wc.tmax.1990.1999.Copy1,
              mean_temp = wc.tavg.1990.1999,
              prec = wc.prec.1990.1999,
              tt2city= tt2cities.2000,
              d2city = dcities,
              d2road = d2roads,
              lon = x,
              lat = y)

      d_pa$land_cover <- factor(d_pa$land_cover, levels=sequence(13),
                                 labels = c("desert","semi_arid","dense_short","tree_short","tree_med","tree_tall",
                                           "salt pan", "sparse_veg_wetland","dense_short_wetland", "tree_short_wetland",
                                           "tree_med_wetland","tree_tall_wetland","water"))  
      
      d_pa$UID <- seq.int(nrow(d_pa))

      dir.create(file.path(paste(f.path3,"/", iso3,"_testPAs",sep="")), recursive = TRUE, showWarnings = FALSE)
      filename_out <- paste(f.path3,"/",iso3,"_testPAs","/","prepped_pa_", testPA$WDPAID,"_wk",gediwk,".RDS", sep="")
      # s3saveRDS(x=d_pa, bucket = "s3://maap-ops-workspace/", object = filename_out, region = "us-west-2") 
    print(filename_out)
     saveRDS(d_pa, filename_out)
    }
  }
#} else if (length(dir(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""),pattern = paste(gediwk,".RDS",sep="")))==length(allPAs)){
#  cat("Step 3.1: prepped PA treatment dataset already exists for ", iso3, "no need for reprocessing\n")
#}