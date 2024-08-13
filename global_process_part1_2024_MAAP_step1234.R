#!/usr/bin/env Rscript

# This global processing script is derived from the global processing notebook 
#the input can be the iso3 code (3-character) for one or multiple countries 
#This script runs steps 1-4 of the global processing code 
#****Aug 6th To get this to run in your own directory, you need to change the working directory to where your outputs will be stowed

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# List of CRAN packages to be installed
cran_packages <- c(
  "s3","foreach", "aws.s3","stringr"
)

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

library("terra")
library("dplyr")
library("sf")
#install.packages("s3")
library("s3")
library("foreach")
library("stringr")
library("aws.s3")

#To test, we define the variables manually. For final version, run the commented out section below
#iso3 <-"ECU"
gediwk <- 24

#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  
  iso3 <- args[1]  #country to process
  #gediwk <- args[2]   #the # of weeks GEDI data to use
  #mproc <- as.integer(args[3])  #the number of cores to use for matching
}
#-------------------------------------------------------------------------------

cat("Step 0: Loading global variables for", iso3,"with wk", gediwk, "data \n")

#f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
f.path2 <- "s3://maap-ops-workspace/my-public-bucket/GEDI_global_PA_v2/"

matching_tifs <- c("wwf_biomes","wwf_ecoreg","lc2000","d2roads", "dcities","dem",
                   "pop_cnt_2000","pop_den_2000","slope", "tt2cities_2000", "wc_prec_1990-1999",
                   "wc_tmax_1990-1999","wc_tavg_1990-1999","wc_tmin_1990-1999" )

ecoreg_key <- read.csv(s3_get(paste(f.path,"wwf_ecoregions_key.csv",sep="")))
#unlink(s3_get(paste(f.path,"wwf_ecoregions_key.csv",sep="")))

allPAs <- readRDS(s3_get(paste(f.path,"WDPA_shapefiles/WDPA_polygons/",iso3,"_PA_poly.rds",sep="")))

MCD12Q1 <- rast(s3_get(paste(f.path,"GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v1_projection_defined_6933.tif",sep="")))
crs(MCD12Q1)  <- "epsg:6933"

world_region <- rast(s3_get(paste(f.path,"GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v1_revised_projection_defined_6933.tif",sep="")))
crs(world_region)  <- "epsg:6933"

s3_path <- paste("/vsis3/maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/WDPA_countries/shp/",iso3,".shp",sep="")

# s3_get_files(c(paste(f.path,"WDPA_countries/shp/",iso3,".shp",sep=""),
#               paste(f.path,"WDPA_countries/shp/",iso3,".shx",sep=""),
#               paste(f.path,"WDPA_countries/shp/",iso3,".prj",sep=""),
#               paste(f.path,"WDPA_countries/shp/",iso3,".dbf",sep="")),confirm = FALSE)


adm <- st_read(s3_path)


# s3_get_files(c(paste(f.path,"WDPA_countries/shp/",iso3,".shp",sep=""),
#               paste(f.path,"WDPA_countries/shp/",iso3,".shx",sep=""),
#               paste(f.path,"WDPA_countries/shp/",iso3,".prj",sep=""),
#               paste(f.path,"WDPA_countries/shp/",iso3,".dbf",sep="")),confirm = FALSE)

# adm <- read_sf(paste("~/shared-buckets/leitoldv/GEDI_global_PA_v2/WDPA_countries/shp/",iso3,".shp",sep=""))
adm_prj <- project(vect(adm), "epsg:6933")

load(s3_get(paste(f.path,"rf_noclimate.RData",sep="")))
#source(s3_get(paste(f.path,"matching_func.R",sep=""))) 
#source(s3_get(paste(f.path,"matching_func_2024.R",sep="")))


# STEP1. Create 1km sampling grid with points only where GEDI data is available; first check if grid file exist to avoid reprocessing 
#if(!file.exists(paste(f.path,"WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep=""))){
  cat("Step 1: Creating 1km sampling grid filter GEDI data for", iso3,"\n")
  GRID.lats <- rast(s3_get(paste(f.path,"EASE2_M01km_lats.tif", sep="")))
  GRID.lons <- rast(s3_get(paste(f.path,"EASE2_M01km_lons.tif", sep="")))
  GRID.lats.adm   <- crop(GRID.lats, adm_prj)
  GRID.lats.adm.m <- mask(GRID.lats.adm, adm_prj)
  GRID.lons.adm   <- crop(GRID.lons, adm_prj)
  GRID.lons.adm.m <- mask(GRID.lons.adm, adm_prj)
  rm(GRID.lats, GRID.lons, GRID.lats.adm, GRID.lons.adm)
  
  #1.3) extract coordinates of raster cells with valid GEDI data in them
  gedi_folder <- paste("~/my-public-bucket/GEDI_global_PA_v2/WDPA_gedi_L4A_tiles/",sep="")
  tileindex_df <- read.csv(s3_get(paste(f.path2,"vero_1deg_tileindex/tileindex_",iso3,".csv", sep="")))
  iso3_tiles <- tileindex_df$tileindexiso3_tiles <- tileindex_df$tileindex
    
  GRID.coords <- data.frame()
  for(i in 1:length(iso3_tiles)){
    
    iso3_tile_in <- paste("tile_num_",iso3_tiles[i],sep="")

    print(paste(iso3_tile_in," processing",sep=""))
    print(paste(gedi_folder,iso3_tile_in,"_L4A.gpkg",sep=""))
    #if(!file.exists(paste(gedi_folder,iso3_tile_in,"_L4A.gpkg",sep=""))){
    #    print(paste(iso3_tile_in," does not exist",sep=""))
    #    } else {
    print(paste(iso3_tile_in," processing",sep=""))
    gedi_data <- read_sf(paste(gedi_folder,iso3_tile_in,"_L4A.gpkg",sep="")), int64_as_string = TRUE) %>%
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

  filename_out <- paste(f.path2, "WDPA_grids/",iso3,"_grid_wk",gediwk,".RDS", sep="")
  s3saveRDS(x=GRID.for.matching, bucket = "s3://maap-ops-workspace/", object = filename_out, region = "us-west-2")

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
  GRID.pts.nonPA <- project(GRID.for.matching, "epsg:4326")
  for(i in 1:length(allPAs)){
    PA          <- vect(allPAs[i,])
    PA_prj      <- project(PA, "epsg:6933")
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
    ras <- rast(s3_get(paste(f.path, "WDPA_input_vars_GLOBAL/",matching_tifs[j],".tif", sep="")))
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
  
  filename_out <- paste(f.path2, "WDPA_grids/",iso3,"_prepped_control_wk",gediwk,".RDS", sep="")
  s3saveRDS(x = d_control, bucket = "s3://maap-ops-workspace/", object = filename_out, region = "us-west-2")  
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
        ras <- rast(s3_get(paste(f.path, "WDPA_input_vars_GLOBAL/",matching_tifs[j],".tif", sep="")))
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

      filename_out <- paste(f.path2, "WDPA_matching_results/",iso3,"_prepped_pa_",testPA$WDPAID,"_wk",gediwk,".RDS", sep="")
      s3saveRDS(x=d_pa, bucket = "s3://maap-ops-workspace/", object = filename_out, region = "us-west-2")  
    }
  }
#} else if (length(dir(paste(f.path,"WDPA_matching_points/",iso3,"/",iso3,"_testPAs","/",sep=""),pattern = paste(gediwk,".RDS",sep="")))==length(allPAs)){
#  cat("Step 3.1: prepped PA treatment dataset already exists for ", iso3, "no need for reprocessing\n")
#}


##****ADDED to Vero's code 8/6/24, taking out parallels
##*************************************************************************
#STEP4. Set up spatial points data frames (control + each PA) for point matching
# if (file.exists(paste(f.path,"WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_matching_output_wk",gediwk,".RDS", sep=""))){

#foreach_rbin function taken from matching_func.R
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

cat("Step 4: Performing matching for", iso3,"\n")
d_control_local <- readRDS(file=paste(f.path2,"WDPA_grids/",iso3,"_prepped_control_wk",gediwk,".RDS", sep=""))
d_control_local <-d_control_local[complete.cases(d_control_local), ]  #filter away non-complete cases w/ NA in control set

if(!dir.exists(paste("WDPA_matching_results/",iso3,"_wk",gediwk,"/",sep=""))){
  # cat("Matching result dir does not EXISTS\n")
  dir.create(file.path(paste("WDPA_matching_results/",iso3,"_wk",gediwk,"/",sep="")))
  d_PAs <- list.files(paste("WDPA_matching_results/",iso3,"_testPAs/", sep=""), pattern=paste("wk",gediwk,sep=""), full.names=FALSE)
} else if (dir.exists(paste("WDPA_matching_results/",iso3,"_wk",gediwk,"/",sep=""))){   #if matching result folder exists, check for any PAs w/o matched results
  pattern1 = c(paste("wk",gediwk,sep=""),"RDS")
  matched_PAid <- list.files(paste("WDPA_matching_results/",iso3,"_wk",gediwk,"/",sep=""), full.names = FALSE, pattern=paste0(pattern1, collapse="|"))%>%
    readr::parse_number() %>% unique()
  d_PAs<- list.files(paste("WDPA_matching_results/",iso3,"_testPAs/", sep=""), pattern=paste("wk",gediwk,sep=""), full.names=FALSE)
  d_PA_id <- d_PAs %>% readr::parse_number()
  runPA_id1 <- d_PA_id[!(d_PA_id %in% matched_PAid)]
  
  matched_all <- list.files(paste("WDPA_matching_results/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = FALSE)
  # registerDoParallel(3)
  matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %do% {   #non-NA matched results
    matched_PAs=c()
    # print(this_rds)
    if(nchar(iso3)>3){
      id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[4]  
    } else {
      id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[3]
    }
    matched <- readRDS(paste("WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""))
    if(!is.null(matched)){
      if(nrow(matched)!=0){
        matched_PAs=c(matched_PAs,this_rds) 
      }
    }else {
      # print(this_rds)
      matched_PAs=matched_PAs
    }
    return(matched_PAs)
  }
  # stopImplicitCluster()
  
  if(!is.null(matched_PAs)){
    fullmatch_ids <- matched_PAs %>% readr::parse_number()
    runPA_id2 <- d_PA_id[!(d_PA_id %in% fullmatch_ids)]
    runPA_id <- c(runPA_id1,runPA_id2)
    
  } else{
    fullmatch_ids <- d_PAs %>% readr::parse_number()
    runPA_id2 <- fullmatch_ids#d_PA_id[!(d_PA_id %in% fullmatch_ids)]
    runPA_id <- c(runPA_id1,runPA_id2)
    
  }
  
  if (length(runPA_id)>0){
    # Pattern2 <-  paste(runPA_id, collapse="|")
    t <- d_PA_id %in% runPA_id
    runPA <-  d_PAs[t]
    d_PAs <- runPA
  } else {
    d_PAs <- NULL
  }
  write.csv(d_PAs, paste(f.path"WDPA_matching_results/", iso3, "_wk_", gediwk, "_null_matches_rerun.csv",sep=""))
  cat("Step 4: need to rerun ", length(d_PAs),"PAs\n")
}

# registerDoParallel(mproc)
# cat("Parallel processing",getDoParWorkers(),"PAs \n")
# startTime <- Sys.time()
foreach(this_pa=d_PAs,.combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr')) %do% {
  pa <- this_pa
  id_pa <-pa %>%str_split("_") %>% unlist %>% .[3]
  # cat(id_pa, "in",iso3,"\n")
  cat("No.", match(pa,d_PAs),"of total",length(d_PAs),"PAs in ", iso3, "\n" )
  d_pa <- readRDS(paste("WDPA_matching_results/",iso3,"_testPAs/",pa, sep=""))
  d_filtered_prop <- tryCatch(propensity_filter(d_pa, d_control_local), error=function(e) return(NA))  #return a df of control and treatment after complete cases and propensity filters are applied
  # cat("Propensity score filtered DF dimension is",dim(d_filtered_prop),"\n")
  d_wocat_all <- tryCatch(filter(d_filtered_prop, status),error=function(e) return(NA))
  d_control_all <- tryCatch(filter(d_filtered_prop, !status),error=function(e) return(NA))
  
  n_control <- dim(d_control_all)[1]
  # ids_all <- d_control_all$UID   #seq(1,n_control)
  ids_all0 <- tryCatch(d_control_all["UID"], error=function(e) return(NA))
  ids_all <- d_control_all["UID"]
  set.seed(125)
  # cat("Using number of cores:",getDoParWorkers(),"\n")
  N <- ceiling(nrow(d_wocat_all)/300)
  l <- tryCatch(split(d_wocat_all, sample(1:N, nrow(d_wocat_all), replace=TRUE)),error=function(e) return(NULL))
  # l <- tryCatch(split(d_wocat_all, (as.numeric(rownames(d_wocat_all))-1) %/% 300),error=function(e) return(0))
  
  if (length(l)<900 && length(l)>0 ){
    pa_match <- data.frame()
    for (pa_c in 1:length(l)){
      ids_all <- d_control_all["UID"]
      cat("chunk",pa_c,"out of ",length(l), "chunks of PA", id_pa,"\n")
      
      d_wocat_chunk <- l[[pa_c]]
      # #sample the control dataset to the size of the sample dataset, keep unsampled ids to iterate until full number of matches found
      n_treatment <- dim(d_wocat_chunk)[1]
      
      t <- ifelse(floor(n_control/n_treatment)<=7, ifelse(floor(n_control/n_treatment)<1, 1,floor(n_control/n_treatment)),7)   #floor(n_control/n_treatment))
      n_sample <- round(n_treatment*t)    #now the n_control is 1.4 times the number of n_treatment, 7 will set the if ststament below to flase
      m_all2_out <- data.frame()
      Bscore <- data.frame()
      n_matches <- 0
      tryCatch(
        while(n_matches < n_treatment){
          n_ids <- length(ids_all)
          # cat("n ids",n_ids,"\n")
          if(n_ids > n_sample){
            set.seed(125)
            sample_ids_bar <- sample(ids_all, n_sample)
            sample_ids <- sample(ids_all0, n_sample)
            d_control_sample <- d_control_all[d_control_all["UID"] %in% sample_ids,]
            ids_all <-setdiff(ids_all, sample_ids_bar)    #ids_all[-sample_ids]
            # cat("protected uid", head(d_wocat_chunk$UID),"\n")
            # All approaches
            new_d <- tryCatch(rbind(d_wocat_chunk,d_control_sample),error=function(e) return(NULL))
            # new_d <- tryCatch(rbind(d_wocat_chunk,d_control_all),error=function(e) return(NULL))
            #create a smaller distance matrix
            m_all <- tryCatch(match_wocat(new_d, pid=id_pa),error=function(e) return(NULL))
            # m_all <- match_wocat(new_d)
            m_all2 <- tryCatch(m_all[1,],error=function(e) return(NULL))
            # m_all2 <- m_all[1,]
            n_matches_temp <- tryCatch(nrow(m_all2$df),error=function(e) return(NULL))
            # n_matches_temp <- nrow(m_all2$df)
            if(!is.null(n_matches_temp)){
              # n_matches <- n_matches + nrow(m_all2$df)
              m_all2$df$pa_id <- rep(id_pa,n_matches_temp)
              m_all2_out <- rbind(m_all2_out, m_all2$df)
              matched_protected <- m_all2$df %>% dplyr::filter(status==TRUE)
              matched_control <- m_all2$df %>% dplyr::filter(status==FALSE)
              cat("matched_protected", nrow(matched_protected),"\n")
              n_matches <- n_matches + nrow(matched_protected)
              d_wocat_chunk <- d_wocat_chunk[-(match(matched_protected["UID"],d_wocat_chunk["UID"])),]
              # d_control_all <- d_control_all[-(match(matched_control$UID,d_control$UID)),]
            } 
            # ids_all <-setdiff(ids_all, sample_ids)
            ids_all0 <-setdiff(ids_all0, matched_control["UID"])
            # else {
            #   n_treatment <- 0  #if not macthes are found in this sampling
            # }
          } else {n_treatment <- n_matches}
        }, error=function(e) return(NULL))
      # ids_all0 <-setdiff(ids_all0, matched_control$UID)
      match_score <- m_all2_out
      cat(table(match_score$status),"\n")
      pa_match <- rbind(pa_match,match_score)
    }
  } else if (length(l)>=900){
    registerDoParallel(4)
    pa_match <- foreach(pa_c=1:length(l), .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr'))%dopar%{
      # cat("Matching treatment chunk", pa_c, "out of", length(l), "for PA", id_pa,"\n")
      cat("chunk",pa_c,"out of ",length(l), "chunks of PA", id_pa,"\n")
      # cat("head control",head(ids_all0),"\n")
      d_wocat_chunk <- l[[pa_c]]
      # #sample the control dataset to the size of the sample dataset, keep unsampled ids to iterate until full number of matches found
      n_treatment <- dim(d_wocat_chunk)[1]
      # cat( "n control", length(ids_all0),"\n")
      t <- ifelse(floor(n_control/n_treatment)<=7, ifelse(floor(n_control/n_treatment)<1, 1,floor(n_control/n_treatment)),7)   #floor(n_control/n_treatment))
      n_sample <- round(n_treatment*t)    #now the n_control is 1.4 times the number of n_treatment, 7 will set the if ststament below to flase
      m_all2_out <- data.frame()
      Bscore <- data.frame()
      n_matches <- 0
      
      tryCatch(
        while(n_matches < n_treatment){
          
          n_ids <- length(ids_all0)
          # cat("n ids",n_ids,"\n")
          if(n_ids > n_sample){
            set.seed(125)
            sample_ids_bar <- sample(ids_all, n_sample)
            sample_ids <- sample(ids_all0, n_sample)
            d_control_sample <- d_control_all[d_control_all["UID"] %in% sample_ids,]
            ids_all <-setdiff(ids_all, sample_ids)    #ids_all[-sample_ids]
            # cat("protected uid", head(d_wocat_chunk$UID),"\n")
            # All approaches
            new_d <- tryCatch(rbind(d_wocat_chunk,d_control_sample),error=function(e) return(NULL))
            # new_d <- tryCatch(rbind(d_wocat_chunk,d_control_all),error=function(e) return(NULL))
            
            #create a smaller distance matrix
            m_all <- tryCatch(match_wocat(new_d, pid=id_pa),error=function(e) return(NULL))
            # m_all <- match_wocat(new_d)
            m_all2 <- tryCatch(m_all[1,],error=function(e) return(NULL))
            # m_all2 <- m_all[1,]
            n_matches_temp <- tryCatch(nrow(m_all2$df),error=function(e) return(NULL))
            # n_matches_temp <- nrow(m_all2$df)
            if(!is.null(n_matches_temp)){
              # n_matches <- n_matches + nrow(m_all2$df)
              m_all2$df$pa_id <- rep(id_pa,n_matches_temp)
              m_all2_out <- rbind(m_all2_out, m_all2$df)
              matched_protected <- m_all2$df %>% dplyr::filter(status==TRUE)
              matched_control <- m_all2$df %>% dplyr::filter(status==FALSE)
              cat("matched_protected", nrow(matched_protected),"\n")
              n_matches <- n_matches + nrow(matched_protected)
              d_wocat_chunk <- d_wocat_chunk[-(match(matched_protected["UID"],d_wocat_chunk["UID"])),]
              # d_control_all <- d_control_all[-(match(matched_control$UID,d_control$UID)),]
              # 
            } 
            ids_all0 <-setdiff(ids_all0, matched_control["UID"])
            # cat( "n control", length(ids_all0),"\n")
            
            # else {
            #   n_treatment <- 0  #if not macthes are found in this sampling
            # }
          } else {n_treatment <- n_matches}
        }, error=function(e) return(NULL))
      # ids_all0 <-setdiff(ids_all0, matched_control$UID)
      match_score <- m_all2_out
      # cat(table(match_score$status),"\n")
      return(match_score)
    }
    # stopImplicitCluster()
  } else{
    pa_match <- NULL
  }
    s3saveRDS(x = pa_match, bucket = "s3://maap-ops-workspace/", object=paste(f.path2, "WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""), region = "us-west-2")
  # cat("Results exported for PA", id_pa,"\n")
  rm(pa_match)
  return(NULL)
}

# tElapsed <- Sys.time()-startTime
# # cat(tElapsed, "for matching all PAs in", iso3,"\n")
# stopImplicitCluster()
cat("Done matching for",iso3,". Finishing...\n")

