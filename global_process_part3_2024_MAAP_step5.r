#When running Rscript, install these packages first
#conda install -c conda-forge r-paws r-terra r-optmatch r-sf r-dplyr r-plyr r-ggplot2 r-mapview r-stringr r-maptools r-gridExtra r-lattice r-MASS r-foreach r-doParallel r-rlang r-tidyr r-magrittr r-aws.s3 r-rgeos r-rlemon r-svd r-sparsem r-survival r-RItools

gediwk <- 24
mproc <- 2

#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  
  iso3 <- args[1]  #country to process
  out <- args[2]
  # range <- args[3]
#  flag <- args[2]  #"run all" PAs or "run remaining" only
  #gediwk <- args[2]   #the # of weeks GEDI data to use
  #mproc <- as.integer(args[3])  #the number of cores to use for matching
}
#-------------------------------------------------------------------------------


# options(repos = c(CRAN = "https://cloud.r-project.org"))

# # # List of CRAN packages to be installed
# cran_packages <- c(
#   "s3","optmatch", "RItools"
# )

# # Install CRAN packages
# install.packages(cran_packages, dependencies = TRUE)

# options(warn=-1)
# options(dplyr.summarise.inform = FALSE)

packages <- c("sf","rstac","dplyr","plyr","ggplot2","mapview","stringr","terra",
              "foreach","optmatch","doParallel","RItools",
              "rlang","tidyr","magrittr","aws.s3","s3")

package.check <- lapply(packages, FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
})

s3 <- paws::s3()

# iso3<- "GNB"
gediwk<-24

f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
f.path2 <- "s3://maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/Matching_Results/"
gedipath<- "/vsis3/maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/" #Make sure to specify username
f.path3<- file.path(out)
# f.path3<- "~/output/WDPA_matching_results/"
# f.path4<- "~/output/WDPA_matching_results/" #Rename folder to "output" since DPS looks for this, move up in the code, set as default but allow an argument to change output file

# f.path3<- file.path(out) #Rename folder to "output" since DPS looks for this, move up in the as default but allow an argument to change output file


source("matching_func_2024.r")

cat("Step 0: Loading global variables to process country", iso3,"with GEDI data until week", gediwk, "\n")

allPAs <- readRDS(s3_get(paste(f.path,"WDPA_shapefiles/WDPA_polygons/",iso3,"_PA_poly.rds",sep="")))

MCD12Q1 <- rast(s3_get(paste(f.path,"GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v1_projection_defined_6933.tif",sep="")))
crs(MCD12Q1)  <- "epsg:6933"

world_region <- rast(s3_get(paste(f.path,"GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v1_revised_projection_defined_6933.tif",sep="")))
crs(world_region)  <- "epsg:6933"

s3_path <- paste("/vsis3/maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/WDPA_countries/shp/",iso3,".shp",sep="") #Redo this 

adm <- st_read(s3_path)

adm_prj <- project(vect(adm), "epsg:6933")

load(s3_get(paste(f.path,"rf_noclimate.RData",sep="")))

# #GLAD Landcover Raster 
# terra::setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
# terra::setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", "TIF")

# MAAP STAC API URL
catalog_url <- "https://stac.maap-project.org"

vect_adm <- vect(adm)

extent <- sf::st_bbox(vect_adm)

#Call GLAD rasters from STAC once per folder
  glad_change <- stac_to_terra(
  catalog_url = catalog_url,
  bbox = extent,
  collections = "glad-glclu2020-change-v2",
  datetime = "2020-01-01T00:00:00Z",
         )
         
  glad_2020 <- stac_to_terra(
   catalog_url = catalog_url,
   bbox = extent,
   collections = "glad-glclu2020-v2",
   datetime = "2020-01-01T00:00:00Z",
            )
reclass_matrix <- matrix(c(
   0,  1,  1, 
   2, 18, 2,
   19, 24, 3,
   25, 32, 4,
   33, 42, 5,
   43, 48, 6,
   49,99, 99,
   100, 101, 7,
   102, 118, 8,
   119, 124, 9,
   125, 132, 10,
   133, 142, 11,
   143, 148, 12,
   149, 199, 99,
   200, 207, 13,
   208, 255, 99
), ncol = 3, byrow = TRUE)

glad_rast_2020 <- classify(glad_2020, reclass_matrix)

reclass_matrix2 <- matrix(c(
   0,  1,  1, 
   2, 18, 2,
   19, 24, 3,
   25, 48, 4,
   49, 72, 5,
   73, 96, 6,
   97,99, 99,
   100, 101, 7,
   102, 118, 8,
   119, 124, 9,
   125, 148, 10,
   149, 172, 11,
   173, 196, 12,
   197, 207, 99,
   212, 239, 99
), ncol = 3, byrow = TRUE)

glad_change_rast <- classify(glad_change, reclass_matrix2)

#End GLAD Codes

flag <- "run all"

#---------Pull appropriate GEDI Tiles-------#
# TO DO: SEPARATE TILES BY COUNTRY
results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                            Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/WDPA_gedi_L2A_tiles/",iso3,"/",sep=""))
all_gedil2_f <- sapply(results$Contents, function(x) {x$Key})
pattern=paste(".gpkg",sep="")
all_gedil2_f <- grep(pattern, all_gedil2_f, value=TRUE)
all_gedil2_f <- basename(all_gedil2_f)

results4 <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                            Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/WDPA_gedi_L4A_tiles/",iso3,"/",sep=""))
all_gedil4_f <- sapply(results4$Contents, function(x) {x$Key})
pattern4=paste(".gpkg",sep="")
all_gedil4_f <- grep(pattern4, all_gedil4_f, value=TRUE)
all_gedil4_f <- basename(all_gedil4_f)
results2b <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                        Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/WDPA_gedi_L2B_tiles/",iso3,"/",sep=""))
all_gedil2b_f <- sapply(results2b$Contents, function(x) {x$Key})
pattern=paste(".gpkg",sep="")
all_gedil2b_f <- grep(pattern, all_gedil2b_f, value=TRUE)
all_gedil2b_f <- basename(all_gedil2b_f)



#---------------STEP5. GEDI PROCESSING - using GEDI shots to extract the treatment/control status, also extract the MODIS PFT for AGB prediction---------------- 
# if (file.exists(paste(f.path,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_gedi_extracted_matching_wk",gediwk,".RDS", sep=""))){
cat(paste("Step 5: Performing WK ",gediwk,"GEDI extraction for", iso3,"\n"))
#matched_all <-read.csv(paste(f.path,"WDPA_extract4_residual_PAs/", iso3, "_wk_", gediwk, "_null_matches_rerun.csv",sep="")) 
# matched_all<-list.files(paste(f.path2,iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
# matched_all

results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/Matching_Results/",iso3,"/",iso3,"_wk24/",sep=""))
  matched_all <- sapply(results$Contents, function(x) {x$Key})
  pattern=paste(".RDS")
  matched_all <-grep(pattern, matched_all, value=TRUE)

#Adding limits for runs to only run specified number of PAs
# new<- unlist(regmatches(range, gregexpr("[[:digit:]]+", range)))
# start<-new[1]
# stop<- new[2]
# matched_all<-matched_all[start:stop]

matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','terra')) %do% {   #non-NA matched results
  matched_PAs=c()
  print(this_rds)
  if(nchar(iso3)>3){
    id_pa <- basename(this_rds)%>%readr::parse_number() %>% unique()  
  } else {
    id_pa <- basename(this_rds)%>%readr::parse_number() %>% unique()
  }
  matched <- readRDS(s3_get(paste(f.path2,iso3,"/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep="")))
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
length(matched_PAs)

if(flag=="run all"){  #determine how many PAs to run the extraction process
  matched_PAs <- matched_PAs
  cat("Step 5: runing extraction on all", length(matched_PAs),"of non-NA matched results in", iso3,"\n")
} else if (flag=="run remaining"){
  # pattern1 = c(paste("wk",gediwk,sep=""),"RDS")
  # extracted_PAid <- list.files(paste(f.path2,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",sep=""), full.names = F, pattern=paste0(pattern1, collapse="|"))
    results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/WDPA_GEDI_extract/",sep=""))
  extracted_PAid <- sapply(results$Contents, function(x) {x$Key})
  pattern=paste("wk_",gediwk,sep="")
  extracted_PAid <- basename(grep(pattern, extracted_PAid, value=TRUE))%>%
    readr::parse_number() %>% unique()
  matched_PA_id <- matched_PAs %>% readr::parse_number()
  runPA_id <- matched_PA_id[!(matched_PA_id %in% extracted_PAid)]
  if (length(runPA_id)>0){
    Pattern2 <-  paste(runPA_id, collapse="|")
    runPA <-  matched_PAs[grepl(Pattern2,matched_PAs)]
    # runPA_ind <- str_detect(matched_PAs, paste(runPA_id, collapse = "|"))
    matched_PAs <-runPA
  } else {
    matched_PAs <- NULL
    cat("Step 5 already done for", iso3, "\n")
  }
}

## Changed error catching and loop now works ##
#Oct 15 Updated to split loop into multiple extract_gedi functions

for (tile in seq_along(all_gedil2_f)){
    tile_id <- basename(all_gedil2_f[tile]) %>% readr::parse_number()
    iso_test<-tryCatch({
        extract_gedi2b(iso3 = iso3,tile_id = tile_id,f.path3 = f.path3,gedipath = gedipath)
        }, error = function(e) {
            cat("Error extracting GEDI data for tile:", e$message, "\n")
            return(NULL)
        })
    }




for (this_rds in matched_PAs) {
    
    # Extract PA ID from the filename
    id_pa <- basename(this_rds) %>% readr::parse_number() %>% unique()
    
    # Construct the path to read the RDS file
    rds_path <- paste(f.path2,iso3,"/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep="")
     
    # Read the RDS file with error handling
    matched <- tryCatch({
        readRDS(s3_get(rds_path))
    }, error = function(e) {
        cat("Error reading RDS file for PA", id_pa, ":", e$message, "\n")
        return(NULL)
    })
    
    # Skip iteration if matched is NULL or has no rows
    if (is.null(matched) || nrow(matched) == 0) {
        cat("Matched result is null or empty for PA", id_pa, "quitting...\n")
        next  # Skip to the next iteration
    }
    
    # Convert matched data frame to raster stack with error handling
    mras <- tryCatch({
        matched2ras(matched)
    }, error = function(e) {
        cat("Error converting matched data to raster stack for PA", id_pa, ":", e$message, "\n")
        return(NULL)
    })
    
    # Check if raster stack is valid
    if (is.null(mras) || any(table(mras$status[]) == 0)) {
        cat("Rasterized results unbalanced or null for PA", id_pa, "quitting...\n")
        next  # Skip to the next iteration
    }
    
    # Start timing for extraction
    startTime <- Sys.time()

    # iso_test<-extract_gedi2b(iso3 = iso3)
    # Extract GEDI data with error handling
    extracted<-list.files(f.path3, pattern=".gpkg", full.names = TRUE)
    iso_matched_gedi <- tryCatch({
        extract_gediPart2(matched = matched, mras = mras,extracted = extracted, glad_change_rast=glad_change_rast,glad_rast_2020=glad_rast_2020)
    }, error = function(e) {
        cat("Error extracting GEDI data for PA", id_pa, ":", e$message, "\n")
        return(NULL)
    })
    
    # Check if extraction results are valid
    if (is.null(iso_matched_gedi)) {
        cat("GEDI extraction result is null for PA", id_pa, "quitting...\n")
        next  # Skip to the next iteration
    }
    
    # Calculate elapsed time
    # tElapsed <- Sys.time() - startTime
    # cat(tElapsed, "for extracting all PAs in", iso3, "\n")
    cat("Done GEDI for PA", match(this_rds, matched_PAs), "out of", length(matched_PAs), "\n")

    variables <- c()

    # Loop from 0 to 29
    for (n in 0:29) {
      # Append the desired strings to the vector
      variables <- c(variables, paste("cover_z", n, sep=""))
      variables <- c(variables, paste("pai_z", n, sep=""))
      variables <- c(variables, paste("pavd_z", n, sep=""))
      }
        
        
    selected_columns <- c("pa_id", "status", "shot_number", "glad_change","glad_2020",
                      "UID","fhd_normal","pai","landsat_treecover","rh20","rh70", "rh10", "rh60","rh100",  
                          "rh90", "rh50", "rh40", "rh98","rh80","rh30","rh25","rh75")
    # Process and select columns
    iso_matched_gedi <- iso_matched_gedi %>%
        dplyr::select(c(selected_columns, variables))
        
    # Determine continent mode
    continent <- unique(iso_matched_gedi$region) %>% getmode()
    
    # Print dimensions of output dataframe
    cat('Output df dimensions:', dim(iso_matched_gedi), "\n")
    
    # Create output directory if it does not exist
    dir.create(file.path(f.path3, "WDPA_GEDI_extract"), recursive = TRUE, showWarnings = FALSE)
    
    iso_matched_gedi_sf <- st_as_sf(iso_matched_gedi, wkt = "geometry", crs = 4326)  # Assuming "geometry" column contains WKT format
    
    
    # # Save results to RDS and CSV files
    # saveRDS(iso_matched_gedi, file = paste(f.path3, "/WDPA_GEDI_extract/", iso3, "_pa_", id_pa, 
    #                                        "_gedi_wk_", gediwk, "_conti_", "biome_", pabiome, ".RDS", sep = ""))
    st_write(iso_matched_gedi_sf, paste(f.path3, "/WDPA_GEDI_extract/", iso3, "_pa_", id_pa, 
                                              "_iso_matched_gedi_sub_wk_", gediwk, ".gpkg", sep = ""))
    
    cat(id_pa, "in", iso3, "results are written to directory\n")
}