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
  range <- args[3]
#  flag <- args[2]  #"run all" PAs or "run remaining" only
  #gediwk <- args[2]   #the # of weeks GEDI data to use
  #mproc <- as.integer(args[3])  #the number of cores to use for matching
}
#-------------------------------------------------------------------------------

packages <- c("sf","rstac","dplyr","plyr","ggplot2","mapview","stringr","terra",
              "foreach","optmatch","doParallel","RItools","httr","jsonlite",
              "rlang","tidyr","magrittr","aws.s3","s3")

package.check <- lapply(packages, FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
})

s3 <- paws::s3()

gediwk<-24

f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
f.path2 <- "s3://maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/Matching_Results/"
gedipath<- "/vsis3/maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/"
gedipath2 <- paste("/vsis3/maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/WDPA_GEDI_L2AL2BL4AL4C_Filtered2A/",iso3,sep="")
# gedipath2 <- paste("s3://maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/WDPA_GEDI_L2AL2BL4AL4C/",iso3,sep="")
f.path3<- file.path(out)

source("matching_func_2024.r")

# Add this function near the top of your script (after the source statement):
get_matching_tile_files <- function(iso3, tile_ids) {
  # Get files via S3 API
  s3_results <- s3$list_objects_v2(
    Bucket = "maap-ops-workspace", 
    Prefix = paste0("shared/abarenblitt/GEDI_global_PA_v2/WDPA_GEDI_L2AL2BL4AL4C/", iso3, "/")
  )
  
  matching_files <- c()
  
  if (length(s3_results$Contents) > 0) {
    # Extract key names and filter for GPKG files
    all_keys <- sapply(s3_results$Contents, function(x) x$Key)
    gpkg_files <- grep("\\.gpkg$", all_keys, value = TRUE)
    
    # Convert to full vsis3 paths
    full_paths <- paste0("/vsis3/maap-ops-workspace/", gpkg_files)
    
    # Find files containing the tile IDs
    for (tile_id in tile_ids) {
      pattern <- paste0("combined_tile_", tile_id, "\\.gpkg")
      matches <- grep(pattern, full_paths, value = TRUE)
      matching_files <- c(matching_files, matches)
    }
  }
  
  return(matching_files)
}

cat("Step 0: Loading global variables to process country", iso3,"with GEDI data until week", gediwk, "\n")

allPAs <- readRDS(s3_get(paste(f.path,"WDPA_shapefiles/WDPA_polygons/",iso3,"_PA_poly.rds",sep="")))

MCD12Q1 <- rast(s3_get(paste(f.path,"GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v1_projection_defined_6933.tif",sep="")))
crs(MCD12Q1)  <- "epsg:6933"

world_region <- rast(s3_get(paste(f.path,"GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v1_revised_projection_defined_6933.tif",sep="")))
crs(world_region)  <- "epsg:6933"

s3_path <- paste("/vsis3/maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/WDPA_countries/shp/",iso3,".shp",sep="")

adm <- st_read(paste(sub("s3://","/vsis3/", f.path),"WDPA_countries/shp/",iso3,".shp",sep=""))

adm_prj <- project(vect(adm), "epsg:6933")

load(s3_get(paste(f.path,"rf_noclimate.RData",sep="")))

# MAAP STAC API URL
catalog_url <- "https://stac.maap-project.org"

vect_adm <- vect(adm)

extent <- sf::st_bbox(vect_adm)

#Adding code for AOI tiles to speed up code
file_path <- paste0(f.path,"vero_1deg_tileindex/tileindex_", iso3, ".csv")

# Read the CSV file into a data frame
tileindex_df <- read.csv(s3_get(file_path))

# Extract the s3path column
json_files <- tileindex_df$s3path

# Assign to AOIs variable
AOIs <- json_files

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


                       
# glad_rast_2020 <- rast(paste(gedipath, "WDPA_input_vars_GLOBAL/",iso3,"_glad_2020.tif", sep=""))
# glad_change_rast <- rast(paste(gedipath, "WDPA_input_vars_GLOBAL/",iso3,"_glad_change.tif", sep=""))

flag <- "run all"

#---------Pull appropriate GEDI Tiles-------#
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

results4c <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                        Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/WDPA_gedi_L4C_tiles/",iso3,"/",sep=""))
all_gedil4c_f <- sapply(results4c$Contents, function(x) {x$Key})
pattern=paste(".gpkg",sep="")
all_gedil4c_f <- grep(pattern, all_gedil4c_f, value=TRUE)
all_gedil4c_f <- basename(all_gedil4c_f)

#---------------STEP5. GEDI PROCESSING---------------- 
cat(paste("Step 5: Performing WK ",gediwk,"GEDI extraction for", iso3,"\n"))

results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/Matching_Results/",iso3,"/",iso3,"_wk24/",sep=""))
  matched_all <- sapply(results$Contents, function(x) {x$Key})
  pattern=paste(".RDS")
  matched_all <-grep(pattern, matched_all, value=TRUE)

# Adding limits for runs to only run specified number of PAs
new<- unlist(regmatches(range, gregexpr("[[:digit:]]+", range)))
start<-new[1]
stop<- new[2]
matched_all<-matched_all[start:stop]

matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','terra')) %do% {   
  if(nchar(iso3)>3){
    id_pa <- basename(this_rds)%>%readr::parse_number() %>% unique()  
  } else {
    id_pa <- basename(this_rds)%>%readr::parse_number() %>% unique()
  }
  
  matched <- tryCatch({
    readRDS(s3_get(paste(f.path2,iso3,"/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep="")))
  }, error = function(e) {
    return(NULL)
  })
  
  if(!is.null(matched)){
    if(nrow(matched)!=0){
      return(this_rds)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

# Remove NULL values from the result
if (length(matched_PAs) > 0) {
  matched_PAs <- matched_PAs[!is.null(matched_PAs)]
}

if(flag=="run all"){  
  matched_PAs <- matched_PAs
  cat("Step 5: running extraction on all", length(matched_PAs),"of non-NA matched results in", iso3,"\n")
} else if (flag=="run remaining"){
  results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
    Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/WDPA_GEDI_extract/",iso3,"/WDPA_GEDI_extract/",sep=""))
  extracted_PAid <- sapply(results$Contents, function(x) {x$Key})
  
  pattern=paste("wk_",gediwk,sep="")
  extracted_PAid <- basename(grep(pattern, extracted_PAid, value=TRUE))%>%
    readr::parse_number() %>% unique()
  
  matched_PA_id <- basename(matched_PAs) %>% readr::parse_number()
  runPA_id <- matched_PA_id[!(matched_PA_id %in% extracted_PAid)]
  
  if (length(runPA_id)>0){
    Pattern2 <-  paste(runPA_id, collapse="|")
    runPA <-  matched_PAs[grepl(Pattern2,matched_PAs)]
    matched_PAs <-runPA
  } else {
    matched_PAs <- NULL
    cat("Step 5 already done for", iso3, "\n")
  }
}

## Process GEDI tiles
for (tile in seq_along(all_gedil2_f)){
    tile_id <- basename(all_gedil2_f[tile]) %>% readr::parse_number()
    iso_test<-tryCatch({
        extract_gedi2b(iso3 = iso3,tile_id = tile_id,f.path3 = f.path3,gedipath = gedipath,gedipath2 = gedipath2)
        }, error = function(e) {
            cat("Error extracting GEDI data for tile:", e$message, "\n")
            return(NULL)
        })
    }

# Use the function to get all bboxes
bboxes <- extract_all_bboxes(AOIs)

# Print summary
cat("Successfully extracted", length(bboxes), "bounding boxes out of", length(AOIs), "AOIs\n")

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

    # Check if the "matched" data frame intersects with any of the bboxes
    intersecting_aois <- check_dataframe_intersection(matched, bboxes)

    # Print results
    if (length(intersecting_aois) > 0) {
      cat("\nThe 'matched' data frame intersects with", length(intersecting_aois), "AOI(s):\n")
      
      # Extract and display the URLs of intersecting AOIs
      for (i in seq_along(intersecting_aois)) {
        cat(i, ": ", basename(intersecting_aois[[i]]$url), "\n", sep="")
      }
      
      # Extract just the URLs if needed
      intersecting_urls <- sapply(intersecting_aois, function(x) x$url)
      
    } else {
      cat("\nNo intersections found between the 'matched' data frame and any AOIs.\n")
    }
        
    id_extraction <- extract_tile_ids_from_aois(intersecting_aois)

    # Display the results
    cat("Found", sum(!is.na(id_extraction$full_results$tile_id)), 
        "tile IDs out of", nrow(id_extraction$full_results), "URLs\n")
    cat("Found", length(id_extraction$unique_ids), "unique tile IDs\n\n")

    # Show the unique IDs
    cat("Unique tile IDs:", paste(id_extraction$unique_ids, collapse=", "), "\n")

    # Get just the unique IDs as a vector for further use
    unique_tile_ids <- id_extraction$unique_ids

    # # Get all GeoPackage files
    # gpkg_files <- list.files(gedipath2, pattern="\\.gpkg$", full.names = TRUE)

    # # Find files containing any of our ID numbers
    # extracted <- match_files_with_numbers(gpkg_files, unique_tile_ids)

    # Get matching tile files using S3 API (list.files doesn't work with vsis3)
    extracted <- get_matching_tile_files(iso3, unique_tile_ids)
                                  
                                  
    # Extract GEDI data with error handling
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
    
    cat("Done GEDI for PA", match(this_rds, matched_PAs), "out of", length(matched_PAs), "\n")

    variables <- c()

    # Loop from 0 to 29
    for (n in 0:29) {
      # Append the desired strings to the vector
      variables <- c(variables, paste("cover_z", n, sep=""))
      variables <- c(variables, paste("pai_z", n, sep=""))
      variables <- c(variables, paste("pavd_z", n, sep=""))
      }
        
    selected_columns <- c("pa_id", "status","land_cover","mangrove", "shot_number", "glad_change","glad_2020",
                      "UID","fhd_normal","pai","landsat_treecover","rh20","rh70", "rh10", "rh60","rh100",  
                          "rh90", "rh50", "rh40", "rh98","rh80","rh30","rh25","rh75","wsci")
    # Process and select columns
    iso_matched_gedi <- iso_matched_gedi %>%
        dplyr::select(c(selected_columns, variables))
        
    # Determine continent mode
    continent <- unique(iso_matched_gedi$region) %>% getmode()
    
    # Print dimensions of output dataframe
    cat('Output df dimensions:', dim(iso_matched_gedi), "\n")
    
    # Create output directory if it does not exist
    dir.create(file.path(f.path3, "WDPA_GEDI_extract"), recursive = TRUE, showWarnings = FALSE)
    
    iso_matched_gedi_sf <- st_as_sf(iso_matched_gedi, wkt = "geometry", crs = 4326)
    
    st_write(iso_matched_gedi_sf, paste(f.path3, "/WDPA_GEDI_extract/", iso3, "_pa_", id_pa, 
                                              "_iso_matched_gedi_sub_wk_", gediwk, ".gpkg", sep = ""),append=FALSE)
    
    cat(id_pa, "in", iso3, "results are written to directory\n")
}