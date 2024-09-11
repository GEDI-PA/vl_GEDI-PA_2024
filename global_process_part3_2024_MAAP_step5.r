#When running Rscript, install these packages first
#conda install -c conda-forge r-paws r-terra r-optmatch r-sp r-sf r-rgeos r-dplyr r-plyr r-ggplot2 r-mapview r-stringr r-maptools r-gridExtra r-lattice r-MASS r-foreach r-doParallel r-rlang r-tidyr r-magrittr r-aws.s3 r-rgeos r-rlemon r-svd r-sparsem r-survival

gediwk <- 24
mproc <- 2

#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  
  iso3 <- args[1]  #country to process
  out <- args[2]
  #gediwk <- args[2]   #the # of weeks GEDI data to use
  # mproc <- as.integer(args[2])  #the number of cores to use for matching
}
#-------------------------------------------------------------------------------


# options(repos = c(CRAN = "https://cloud.r-project.org"))

# # List of CRAN packages to be installed
# cran_packages <- c(
#   "s3","optmatch", "RItools"
# )

# # Install CRAN packages
# install.packages(cran_packages, dependencies = TRUE)

# options(warn=-1)
# options(dplyr.summarise.inform = FALSE)

packages <- c("sp","sf","rgeos","dplyr","plyr","ggplot2","mapview","stringr","terra",
              "foreach","optmatch","doParallel","RItools",
              "rlang","tidyr","magrittr","aws.s3","s3")

package.check <- lapply(packages, FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
})

s3 <- paws::s3()

# iso3<- "GNB"
gediwk<-24

f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
gedipath<- "/vsis3/maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/" #Make sure to specify username
f.path3<- file.path(out,"WDPA_matching_results/") #Rename folder to "output" since DPS looks for this, move up in the as default but allow an argument to change output file


source("matching_func_2024.r")

cat("Step 0: Loading global variables to process country", iso3,"with GEDI data until week", gediwk, "\n")

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

s3_path <- paste("/vsis3/maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/WDPA_countries/shp/",iso3,".shp",sep="") #Redo this 

adm <- st_read(s3_path)

adm_prj <- project(vect(adm), "epsg:6933")

load(s3_get(paste(f.path,"rf_noclimate.RData",sep="")))

flag <- "run all"

#---------------STEP5. GEDI PROCESSING - using GEDI shots to extract the treatment/control status, also extract the MODIS PFT for AGB prediction---------------- 
# if (file.exists(paste(f.path,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_gedi_extracted_matching_wk",gediwk,".RDS", sep=""))){
cat(paste("Step 5: Performing WK ",gediwk,"GEDI extraction for", iso3,"\n"))
#matched_all <-read.csv(paste(f.path,"WDPA_extract4_residual_PAs/", iso3, "_wk_", gediwk, "_null_matches_rerun.csv",sep="")) 
matched_all<-list.files(paste(f.path3,iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
matched_all

matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','terra')) %do% {   #non-NA matched results
  matched_PAs=c()
  print(this_rds)
  if(nchar(iso3)>3){
    id_pa <- basename(this_rds)%>%readr::parse_number() %>% unique()  
  } else {
    id_pa <- basename(this_rds)%>%readr::parse_number() %>% unique()
  }
  matched <- readRDS(paste(f.path3,iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""))
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
  pattern1 = c(paste("wk",gediwk,sep=""),"RDS")
  extracted_PAid <- list.files(paste(f.path3,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",sep=""), full.names = F, pattern=paste0(pattern1, collapse="|"))%>%
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
#Sep 9 2024, code works but needs to be updated for running as DPS

for (this_rds in matched_PAs) {
    
    # Extract PA ID from the filename
    id_pa <- basename(this_rds) %>% readr::parse_number() %>% unique()
    
    # Construct the path to read the RDS file
    rds_path <- paste(f.path3, iso3, "_wk", gediwk, "/", iso3, "_pa_", id_pa, "_matching_results_wk24.RDS", sep = "")
    
    # Read the RDS file with error handling
    matched <- tryCatch({
        readRDS(rds_path)
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
    
    # Extract GEDI data with error handling
    iso_matched_gedi <- tryCatch({
        extract_gedi2b(matched = matched, mras = mras)
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
    tElapsed <- Sys.time() - startTime
    cat(tElapsed, "for extracting all PAs in", iso3, "\n")
    cat("Done GEDI for PA", match(this_rds, matched_PAs), "out of", length(matched_PAs), "\n")
    
    # Process and select columns
    iso_matched_gedi <- iso_matched_gedi %>%
        dplyr::select("pa_id", "status", "wwfbiom", "wwfecoreg", "shot_number", "lon_lowestmode",
                      "lat_lowestmode", "rh25", "rh50", "rh75", "rh90", "rh98","UID")
    
    # Determine biome name
    if (length(unique(iso_matched_gedi$wwfbiom)) > 1) {
        pabiome <- iso_matched_gedi$wwfbiom %>%
            unique() %>%
            gsub('\\b(\\pL)\\pL{2,}|.', '\\U\\1', ., perl = TRUE) %>%
            str_c(collapse = "+")
    } else if (length(unique(iso_matched_gedi$wwfbiom)) == 1) {
        pabiome <- iso_matched_gedi$wwfbiom %>%
            unique() %>%
            gsub('\\b(\\pL)\\pL{2,}|.', '\\U\\1', ., perl = TRUE)
    } else {
        pabiome <- iso_matched_gedi$wwfbiom %>% unique()
    }
    
    # Determine continent mode
    continent <- unique(iso_matched_gedi$region) %>% getmode()
    
    # Print dimensions of output dataframe
    cat('Output df dimensions:', dim(iso_matched_gedi), "\n")
    
    # Create output directory if it does not exist
    dir.create(file.path(f.path3, "WDPA_GEDI_extract"), recursive = TRUE, showWarnings = FALSE)
    
    # Save results to RDS and CSV files
    saveRDS(iso_matched_gedi, file = paste(f.path3, "WDPA_GEDI_extract/", iso3, "_pa_", id_pa, 
                                           "_gedi_wk_", gediwk, "_conti_", "biome_", pabiome, ".RDS", sep = ""))
    write.csv(iso_matched_gedi, file = paste(f.path3, "WDPA_GEDI_extract/", iso3, "_pa_", id_pa, 
                                              "_iso_matched_gedi_sub_wk_", gediwk, ".csv", sep = ""))
    
    cat(id_pa, "in", iso3, "results are written to directory\n")
}
