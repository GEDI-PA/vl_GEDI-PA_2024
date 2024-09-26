#!/usr/bin/env Rscript

# This global processing script is derived from the global processing notebook 
#the input can be the iso3 code (3-character) for one or multiple countries 

library("terra")
library("dplyr")
library("sf")
#install.packages("s3")
library("s3")

library("sp")
library("foreach")
library("stringr")
library("aws.s3")
#conda install -c conda-forge r-optmatch #r-ggmap r-hrbrthemes r-Hmisc
#library("optmatch")
#install.packages("doParallel")
library("doParallel")
#install.packages("RItools")    
#library("RItools")

#f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
source(s3_get(paste(f.path,"vl_GEDI-PA_2024/matching_func_2024.R",sep="")))

#To test, we define the variables manually. For final version, run the commented out section below
#iso3 <-"ECU"
gediwk <- 24
mproc <- 1

#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  
  iso3 <- args[1]  #country to process
  flag <- args[2]  #"run all" PAs or "run remaining" only
  #gediwk <- args[2]   #the # of weeks GEDI data to use
  #mproc <- as.integer(args[3])  #the number of cores to use for matching
}
#-------------------------------------------------------------------------------

cat("Step 0: Loading global variables for", iso3,"with wk", gediwk, "data \n")

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

s3_get_files(c(paste(f.path,"WDPA_countries/shp/",iso3,".shp",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".shx",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".prj",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".dbf",sep="")),confirm = FALSE)
adm <- st_read(s3_get(paste(f.path,"WDPA_countries/shp/",iso3,".shp",sep="")))

#s3_path <- paste("/vsis3/maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/WDPA_countries/shp/",iso3,".shp",sep="") #Redo this for the gpkg
#adm <- st_read(s3_path)
adm_prj <- project(vect(adm), "epsg:6933")

load(s3_get(paste(f.path,"rf_noclimate.RData",sep="")))
#source(s3_get(paste(f.path,"matching_func.R",sep="")))
#source(s3_get(paste(f.path,"vl_GEDI-PA_2024/matching_func_2024.R",sep="")))
#source("/projects/my-public-bucket/GEDI_global_PA_v2/vl_GEDI-PA_2024/matching_func_2024.R")

#flag <- "run all"
#flag <- "run remaining"

#---------------STEP5. GEDI PROCESSING ---------------- 
#using GEDI shots to extract the treatment/control status, also extract the MODIS PFT for AGB prediction

# if (file.exists(paste(f.path,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_gedi_extracted_matching_wk",gediwk,".RDS", sep=""))){
cat(paste("Step 5: Performing WK ", gediwk, "GEDI extraction for", iso3,"\n"))
#matched_all <-read.csv(paste(f.path,"WDPA_extract4_residual_PAs/", iso3, "_wk_", gediwk, "_null_matches_rerun.csv",sep="")) 
matched_all <- read.csv(s3_get(paste(f.path,"WDPA_matching_results/",iso3,"_wk24_filelist.csv",sep="")))[,2]
print(length(matched_all))

matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','terra')) %do% {  #non-NA matched results
  matched_PAs=c()
  print(this_rds)
  if(nchar(iso3)>3){
    id_pa <- basename(this_rds)%>%readr::parse_number() %>% unique()  
  } else {
    id_pa <- basename(this_rds)%>%readr::parse_number() %>% unique()
  }
  matched <- readRDS(s3_get(paste(f.path,"WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_pa_",id_pa,"_matching_results_wk",gediwk,".RDS",sep="")))
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
print(length(matched_PAs))

#f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"

#flag <- "run all"
#flag <- "run remaining"

if(flag=="run all"){  #determine how many PAs to run the extraction process
  matched_PAs <- matched_PAs
  cat("Step 5: running extraction on all", length(matched_PAs),"of non-NA matched results in", iso3,"\n")
} else if (flag=="run remaining"){
    #pattern1 = c(paste("wk",gediwk,sep=""),"RDS")
    #extracted_PAid <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,"/",sep=""), full.names = F, pattern=paste0(pattern1, collapse="|"))%>% readr::parse_number() %>% unique()
  
    s3 <- paws::s3()
    bucket <- "maap-ops-workspace"
    prefix <- paste("shared/leitoldv/GEDI_global_PA_v2/WDPA_extract/",iso3,"_wk",gediwk,"/",sep="")
    
    extracted_PAs <- s3$list_objects_v2(Bucket = bucket, Prefix = prefix)
    extracted_fnames <- sapply(extracted_PAs$Contents, function(x) {x$Key})
    pattern <- paste(".RDS",sep="")
    extracted_fnames <- grep(pattern, extracted_fnames, value=TRUE)
    extracted_fnames <- basename(extracted_fnames)

    extracted_PAid <- extracted_fnames %>% readr::parse_number() %>% unique()

  matched_PA_id <- matched_PAs %>% readr::parse_number()
  matched_PAdf <- cbind(matched_PAs, matched_PA_id)
    for(i in 1:nrow(matched_PAdf)){
    matched_PAdf[i,"matched_PA_id"] <- readr::parse_number(matched_PAs[i])
    }
  #runPA_id <- matched_PA_id[!(matched_PA_id %in% extracted_PAid)]
  runPA_id <- matched_PAdf[matched_PAdf[,"matched_PA_id"] != extracted_PAid, "matched_PA_id"]
  if (length(runPA_id)>0){
    #Pattern2 <-  paste(runPA_id, collapse="|")
    #runPA <-  matched_PAs[grepl(Pattern2,matched_PAs)]
    # runPA_ind <- str_detect(matched_PAs, paste(runPA_id, collapse = "|"))
    #matched_PAs <- runPA
    matched_PAs <- matched_PAdf[matched_PAdf[,"matched_PA_id"] != extracted_PAid, "matched_PAs"]  
  } else {
    matched_PAs <- NULL
    cat("Step 5 already done for", iso3, "\n")
  }
}

print(paste("remaining PAs to be extracted =",length(matched_PAs),sep=" "))  ##remaining PAs to be extracted

registerDoParallel(cores=round(mproc))
getDoParWorkers()
startTime <- Sys.time()

foreach(this_rds=matched_PAs, .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr')) %dopar% {
  cat("Extracting for no. ", match(this_rds,matched_PAs),"pa out of", length(matched_PAs),"\n")
  id_pa <- basename(this_rds) %>% readr::parse_number() %>% unique()
  matched <- readRDS(s3_get(paste(f.path,"WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_pa_",id_pa,"_matching_results_wk",gediwk,".RDS",sep="")))
#  matched$pa_id <- rep(id_pa, nrow(matched))
  print(paste("PA id",unique(matched$pa_id),sep=" "))

  if (is.null(matched)==TRUE  | nrow(matched)==0) {
    cat("Matched result is null for PA", id_pa, "quitting...\n")
  } else if (!is.null(matched)==TRUE){
    mras  <- tryCatch(matched2ras(matched),
                      error=function(cond){
                        message(cond)
                        cat("Matched result is likely null for country", iso3,"pa", id_pa, "dimension of the match is", dim(matched),"\n")
                        return(NULL)}) #convert the macthed df to a raster stack 
    print(table(mras$status[]))
    cat("Rasterized results are balanced for PA", id_pa, "\n")
    
    if(table(mras$status[])[2]==0 | table(mras$status[])[1]==0 | is.null(mras)){
      cat("Rasterized results unbalanced for PA", id_pa, "quitting...\n")
    } else {
      startTime <- Sys.time()
      #iso_matched_gedi<- extract_gedi(matched=matched, mras = mras)  #run filtered csvs on mras for extarction
      iso_matched_gedi <- extract_gedi(matched=matched, mras=mras, iso3=iso3)
      tElapsed <- Sys.time()-startTime
      cat(tElapsed, "for extracting all PAs in", iso3,"\n")
      iso_matched_gedi <-  iso_matched_gedi %>%
            dplyr::select("pa_id","status","wwfbiom","wwfecoreg","pft","region",
                          "shot_number","lat_lowestmode","lon_lowestmode",
                          "rh98","agbd","agbd_se") #"filename","rh25","rh50","rh75",
    if (length(unique(iso_matched_gedi$wwfbiom)) >1){
        pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub("\\b(\\p{L})\\p{L}{2,}|.","\\U\\1",.,perl = TRUE)%>% str_c( collapse = "+")
    } else if (length(unique(iso_matched_gedi$wwfbiom))==1){
        pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\p{L})\\p{L}{2,}|.','\\U\\1',.,perl = TRUE)
    } else {
        pabiome <- iso_matched_gedi$wwfbiom %>% unique()
    }
    # papaddd <- unique(iso_matched_gedi$PADDD) %>% getmode()
    continent <- unique(iso_matched_gedi$region) %>% getmode()
    print(paste('output df',dim(iso_matched_gedi)))

    dir.create(file.path(paste("output/WDPA_extract/",iso3,"_wk",gediwk,"/",sep="")),recursive=TRUE)
    saveRDS(iso_matched_gedi, file=paste("output/WDPA_extract/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_gedi_wk_",gediwk,"_conti_","biome_",pabiome,".RDS", sep=""))
#    write.csv(iso_matched_gedi, file=paste("output/WDPA_extract/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_iso_matched_gedi_sub_wk_",gediwk,".csv", sep=""))

    cat("PA#",id_pa,"in",iso3,"result is written to dir\n")
    rm(iso_matched_gedi)
    }
  }
  return(NULL)
}

stopImplicitCluster()
tElapsed <- Sys.time()-startTime
cat(tElapsed, "for extracting all PAs in", iso3,"\n")
cat("Done GEDI extraction for pa in ",iso3,"\n")    

