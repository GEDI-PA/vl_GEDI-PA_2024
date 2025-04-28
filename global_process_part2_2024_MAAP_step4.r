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
library("foreach")
library("stringr")
library("aws.s3")
library("optmatch")
library("doParallel")

s3 <- paws::s3()

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
  # mproc <- as.integer(args[3])  #the number of cores to use for matching
}
#-------------------------------------------------------------------------------

cat("Step 0: Loading global variables for", iso3,"with wk", gediwk, "data \n")

#f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
f.path2 <- "s3://maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/Matching_Results"
gedipath<- "/vsis3/maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/" #Make sure to specify username
f.path3<- file.path(out)

allPAs <- readRDS(s3_get(paste(f.path,"WDPA_shapefiles/WDPA_polygons/",iso3,"_PA_poly.rds",sep="")))

MCD12Q1 <- rast(s3_get(paste(f.path,"GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v1_projection_defined_6933.tif",sep="")))
crs(MCD12Q1)  <- "epsg:6933"

world_region <- rast(s3_get(paste(f.path,"GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v1_revised_projection_defined_6933.tif",sep="")))
crs(world_region)  <- "epsg:6933"

s3_get_files(c(paste(f.path,"WDPA_countries/shp/",iso3,".shp",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".shx",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".prj",sep=""),
              paste(f.path,"WDPA_countries/shp/",iso3,".dbf",sep="")),confirm = FALSE)

adm <- st_read(paste(sub("s3://","/vsis3/", f.path),"WDPA_countries/shp/",iso3,".shp",sep=""))

adm_prj <- project(vect(adm), "epsg:6933")

load(s3_get(paste(f.path,"rf_noclimate.RData",sep="")))

source("matching_func_2024.r")

#STEP4. Set up spatial points data frames (control + each PA) for point matching
# if (file.exists(paste(f.path,"WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_matching_output_wk",gediwk,".RDS", sep=""))){

cat("Step 4: Performing matching for", iso3,"\n")
d_control_local <- readRDS(s3_get(paste(f.path2,"/",iso3,"/WDPA_grids/",iso3,"_prepped_control_wk",gediwk,".RDS", sep="")))
d_control_local <-d_control_local[complete.cases(d_control_local), ]  #filter away non-complete cases w/ NA in control set

#All this file creation should be a function

if(!dir.exists(paste(f.path3,"/",iso3,"_wk",gediwk,"/",sep=""))){
  # cat("Matching result dir does not EXISTS\n")
  dir.create(file.path(f.path3,paste0(iso3,"_wk",gediwk)),recursive=TRUE)

  results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                              Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/Matching_Results","/",iso3,"/",iso3,"_testPAs/",sep=""))
  d_PAs <- sapply(results$Contents, function(x) {x$Key})
    
} else if (dir.exists(paste(f.path3,"/",iso3,"_wk",gediwk,"/",sep=""))){   #if matching result folder exists, check for any PAs w/o matched results

  pattern1 = c(paste("wk",gediwk,sep=""),"RDS")
  matched_PAid <- list.files(paste(f.path3,"/",iso3,"_wk",gediwk,"/",sep=""), full.names =    FALSE, pattern=paste0(pattern1, collapse="|"))%>%
  readr::parse_number() %>% unique()
  results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                              Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/Matching_Results","/",iso3,"/",iso3,"_testPAs/",sep=""))
  d_PAs <- sapply(results$Contents, function(x) {x$Key})
  pattern=paste(".RDS",sep="")
  d_PAs <- grep(pattern, d_PAs, value=TRUE)
  d_PAsBase <-basename(d_PAs)
  d_PA_id <- d_PAsBase %>% readr::parse_number()
  runPA_id1 <- d_PA_id[!(d_PA_id %in% matched_PAid)]
    
  matched_all <- list.files(paste(f.path3,"/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = FALSE)
  
  # registerDoParallel(3) #_@_@
  matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','terra')) %do% { #_@_@ #non-NA matched results
    matched_PAs=c()
    # print(this_rds)
    if(nchar(iso3)>3){
      id_pa <- this_rds  %>% readr::parse_number() %>% unique()
    } else {
      id_pa <- this_rds  %>% readr::parse_number() %>% unique()
    }
    matched <- readRDS(paste(f.path3,"/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""))
    print(matched)
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
  # stopImplicitCluster() #_@_
  
  if(!is.null(matched_PAs)){
    fullmatch_ids <- matched_PAs %>%readr::parse_number() %>% unique()
    runPA_id <- d_PA_id[!(d_PA_id %in% fullmatch_ids)]
    # runPA_id <- c(runPA_id1,runPA_id2)
    print(runPA_id)
    
  } else{
    fullmatch_ids <- d_PAsBase %>%readr::parse_number() %>% unique()
    runPA_id <- fullmatch_ids#d_PA_id[!(d_PA_id %in% fullmatch_ids)]
    # runPA_id <- c(runPA_id1,runPA_id2)
    print(runPA_id)
  }
  
  if (length(runPA_id)>0){
    # Pattern2 <-  paste(runPA_id, collapse="|")
    t <- d_PA_id %in% runPA_id
    runPA <-  d_PAs[t]
    d_PAs <- runPA
  } else {
    d_PAs <- NULL
  }
  d_PAs
  write.csv(d_PAs, paste(f.path3,"/", iso3, "_wk_", gediwk, "_null_matches_rerun.csv",sep=""))
  cat("Step 4: need to rerun ", length(d_PAs),"PAs\n")
}


# registerDoParallel(1) #_@_@
# cat("Parallel processing",getDoParWorkers(),"PAs \n")
# startTime <- Sys.time()

results <- s3$list_objects_v2(Bucket = "maap-ops-workspace", 
                              Prefix=paste("shared/abarenblitt/GEDI_global_PA_v2/Matching_Results","/",iso3,"/",iso3,"_testPAs/",sep=""))
d_PAs <- sapply(results$Contents, function(x) {x$Key})
d_PAs <- sapply(results$Contents, function(x) {x$Key})
pattern=paste(".RDS",sep="")
d_PAs <- grep(pattern, d_PAs, value=TRUE)
d_PAs<- basename(d_PAs)

foreach(this_pa=d_PAs,.combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','optmatch')) %do% { #_@_@
  pa <- basename(this_pa)
    print(pa)
  id_pa <-basename(this_pa)%>%readr::parse_number() %>% unique() #%>%str_split("_") %>% unlist %>% .[4] #With new files, check where PA ID is in string
  # cat(id_pa, "in",iso3,"\n")
    path <- paste(f.path2,"/",iso3,"/",iso3,"_testPAs/",pa, sep="")
  cat("No.", match(pa,d_PAs),"of total",length(d_PAs),"PAs in ", iso3, "\n" )
  d_pa <- readRDS(s3_get(path))
  d_pa$mangrove[is.na(d_pa$mangrove)] <- 0  
  d_filtered_prop <- tryCatch(propensity_filter(d_pa, d_control_local), error=function(e) return(NA))  #return a df of control and treatment after complete cases and propensity filters are applied
  # cat("Propensity score filtered DF dimension is",dim(d_filtered_prop),"\n")
  d_wocat_all <- tryCatch(filter(d_filtered_prop, status),error=function(e) return(NA))
  d_control_all <- tryCatch(filter(d_filtered_prop, !status),error=function(e) return(NA))
  
  n_control <- dim(d_control_all)[1]
  # ids_all <- d_control_all$UID   #seq(1,n_control)
  ids_all0 <- tryCatch(d_control_all[['UID']], error=function(e) return(NA))
  ids_all <- d_control_all[['UID']]
  set.seed(125)
  # cat("Using number of cores:",getDoParWorkers(),"\n")
  N <- ceiling(nrow(d_wocat_all)/300)
  l <- tryCatch(split(d_wocat_all, sample(1:N, nrow(d_wocat_all), replace=TRUE)),error=function(e) return(NULL))
  # l <- tryCatch(split(d_wocat_all, (as.numeric(rownames(d_wocat_all))-1) %/% 300),error=function(e) return(0))
 
  if (length(l)<900 && length(l)>0 ){
    pa_match <- data.frame()
    for (pa_c in 1:length(l)){
      ids_all <- d_control_all[['UID']]
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
            d_control_sample <- d_control_all[d_control_all[['UID']] %in% sample_ids,]
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
              d_wocat_chunk <- d_wocat_chunk[-(match(matched_protected$UID,d_wocat_chunk[['UID']])),]
              # d_control_all <- d_control_all[-(match(matched_control$UID,d_control$UID)),]
            } 
            # ids_all <-setdiff(ids_all, sample_ids)
            ids_all0 <-setdiff(ids_all0, matched_control$UID)
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
    # registerDoParallel(4) #_@_@
    pa_match <- foreach(pa_c=1:length(l), .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','optmatch'))%do%{  #_@_@
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
            d_control_sample <- d_control_all[d_control_all[['UID']] %in% sample_ids,]
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
              d_wocat_chunk <- d_wocat_chunk[-(match(matched_protected$UID,d_wocat_chunk[['UID']])),]
              # d_control_all <- d_control_all[-(match(matched_control$UID,d_control$UID)),]
              # 
            } 
            ids_all0 <-setdiff(ids_all0, matched_control$UID)
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
    # stopImplicitCluster() #_@_@
  } else{
    pa_match <- NULL
  }
  # s3saveRDS(x = pa_match, bucket = "maap-ops-workspace", object=paste(f.path3,iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""), region = "us-west-2")
  path<- paste(f.path3,"/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep="")
  saveRDS(pa_match, file=path)
  print(pa_match)
  print(path)
  cat("Results exported for PA", id_pa,"\n")
  rm(pa_match)                                    
  return(NULL)
}

# tElapsed <- Sys.time()-startTime
# cat(tElapsed, "for matching all PAs in", iso3,"\n")
# stopImplicitCluster() #_@_@
cat("Done matching for",iso3,". Finishing...\n")