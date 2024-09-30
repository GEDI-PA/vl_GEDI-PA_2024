#!/usr/bin/env Rscript

# This global processing script is derived from the global processing notebook 
#the input can be the iso3 code (3-character) for one or multiple countries 

#install.packages("s3")
#install.packages("doParallel")
#install.packages("RItools")    

library("terra")
library("dplyr")
library("sf")
library("s3")
library("sp")
library("foreach")
library("stringr")
library("aws.s3")
#library("doParallel")
#library("RItools")

#f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
#iso3 <- "BaE"
gediwk <- 24

#source(s3_get(paste(f.path,"matching_func.R",sep="")))
source(paste(f.path,"vl_GEDI-PA_2024/matching_func_2024.R",sep=""))

#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  
  iso3 <- args[1]  #country to process
#  flag <- args[2]  #"run all" PAs or "run remaining" only
  #gediwk <- args[2]   #the # of weeks GEDI data to use
  #mproc <- as.integer(args[3])  #the number of cores to use for matching
}
#-------------------------------------------------------------------------------

## STEP6: [FIGURE 4B] Calculating per pa summary stats, 1 pa per row, contain shot#/PA---------------------------- 
cat(paste("Step 6: calculating per pa summary stats for", iso3,"\n"))

gedi_paf <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)

for (this_paf in gedi_paf){
  pa_metrics <- readRDS(this_paf) %>% unique()
  if (length(table(pa_metrics$status))<2) {
    cat(iso3, this_paf, "has 0 protected or treatment \n")
  } else if (table(pa_metrics$status)[1]!=0 && table(pa_metrics$status)[2]!=0) {
    #filter the datafrme by number of shots in each cell, which is equivalent to the number of occurence of each unique UID code
##    tt <- table(pa_metrics$UID)
##    qcellid <- table(pa_metrics$UID)[tt>5] %>% names()
##    pa_metrics_filtered <- pa_metrics %>% dplyr::filter(UID %in% qcellid)
    #calc summary stats for each country 
    #pa_stats_summary <- pa_metrics_filtered %>%
    pa_stats_summary <- pa_metrics %>%
      group_by(status) %>% 
      dplyr::mutate(pa_id=as.character(pa_id)) %>%
      dplyr::summarise(pa_id=na.omit(unique(pa_id)),
                       count=length(rh98),
                       meanrh98=mean(rh98, na.rm=TRUE), sdrh98=sd(rh98, na.rm=TRUE),
                       meanagbd=mean(agbd, na.rm=TRUE), sdagbd=sd(agbd, na.rm=TRUE),
                       wwfbiom=getmode(wwfbiom), wwfecoreg=getmode(wwfecoreg),
                       REGION=getmode(region), PFT=getmode(pft)) %>% 
      tidyr::pivot_wider(names_from=status, values_from= setdiff(names(.),c("pa_id", "status")))
    #writeLine to a large txt file where world pas stats are
    pa_stats_summary$iso3 <- iso3

    if(!file.exists(paste(f.path,"WDPA_extract/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""))){
      print("pa_stats_summary file does not exist")
      write.csv(pa_stats_summary, file=paste(f.path,"WDPA_extract/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""), row.names = FALSE)
    } else if (file.exists(paste(f.path,"WDPA_extract/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""))){
      print("pa_stats_summary file exists so appending to existing file")
      write.table(pa_stats_summary, file=paste(f.path,"WDPA_extract/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""),
                  sep=",", append=TRUE , row.names=FALSE, col.names=FALSE)
    #will not overwrite but append to existing files
    } 
  }
}

cat("Done summarizing pa-level stats for region", iso3, "\n")  
