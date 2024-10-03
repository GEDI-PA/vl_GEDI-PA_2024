#install.packages("s3")
#install.packages("doParallel")
#install.packages("RItools")    

library("terra")
library("dplyr")
library("sf")
#library("s3")
library("sp")
library("foreach")
library("stringr")
library("aws.s3")
#library("doParallel")
#library("RItools")

#f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
gediwk <- 24
#iso3 <- "Bpt"

iso3 <- "BaE"
gedi_paf_BaE <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
length(gedi_paf_BaE)

iso3 <- "BaW"
gedi_paf_BaW <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
length(gedi_paf_BaW)

iso3 <- "Bca"
gedi_paf_Bca <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
length(gedi_paf_Bca)

iso3 <- "Bce"
gedi_paf_Bce <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
length(gedi_paf_Bce)

iso3 <- "Bma"
gedi_paf_Bma <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
length(gedi_paf_Bma)

iso3 <- "Bpp"
gedi_paf_Bpp <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
length(gedi_paf_Bpp)

iso3 <- "Bpt"
gedi_paf_Bpt <- list.files(paste(f.path,"WDPA_extract/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
length(gedi_paf_Bpt)

#------------
gedi_paf_ALL <- c(gedi_paf_BaE, gedi_paf_BaW, gedi_paf_Bca,
                  gedi_paf_Bce, gedi_paf_Bma, gedi_paf_Bpp, gedi_paf_Bpt)
print(length(gedi_paf_ALL))

#------------
gedi_paid_ALL <- c()

for(i in 1:length(gedi_paf_ALL)){

    this_paf <- gedi_paf_ALL[i]
    this_paid <- basename(this_paf) %>% readr::parse_number() %>% unique()
    gedi_paid_ALL <- c(gedi_paid_ALL, this_paid)
    #print(this_paid)
}

print(length(gedi_paid_ALL))
print(length(unique(gedi_paid_ALL)))

#------------
gedi_paid_DUP <- gedi_paid_ALL[which(duplicated(gedi_paid_ALL))]
length(gedi_paid_DUP)

gedi_paid_UNQ <- setdiff(gedi_paid_ALL, gedi_paid_DUP)
length(gedi_paid_UNQ)

#------------
#for(i in 1:length(gedi_paid_DUP)){
#
#DUP_id <- gedi_paid_DUP[i]
#print(DUP_id)
#DUP_paf <- gedi_paf_ALL[which(gedi_paid_ALL == DUP_id)]
#print(length(DUP_paf))
#
#    if(length(DUP_paf) == 2){
#   
#        pa_metrics1 <- readRDS(DUP_paf[1]) %>% unique()
#        pa_metrics2 <- readRDS(DUP_paf[2]) %>% unique()    
#        pa_metrics_ALL <- rbind(pa_metrics1, pa_metrics2)    
#        wwfbiom <- strsplit(strsplit(DUP_paf[1], split=c("_conti_biome_"))[[1]][2], split=".RDS")[[1]]
#    write.csv(pa_metrics_ALL,
#          #file=paste(f.path,"/WDPA_extract/pa_stats_ALL/BRA_pa_",DUP_id,"_gedi_wk24_",wwfbiom,".csv",sep=""))
#    }
#
#    else if(length(DUP_paf) == 3){
#        pa_metrics1 <- readRDS(DUP_paf[1]) %>% unique()
#        pa_metrics2 <- readRDS(DUP_paf[2]) %>% unique()
#        pa_metrics3 <- readRDS(DUP_paf[3]) %>% unique()
#        pa_metrics_ALL <- rbind(pa_metrics1, pa_metrics2, pa_metrics3)
#        wwfbiom <- strsplit(strsplit(DUP_paf[1], split=c("_conti_biome_"))[[1]][2], split=".RDS")[[1]]
#    write.csv(pa_metrics_ALL,
#     file=paste(f.path,"/WDPA_extract/pa_stats_ALL/BRA_pa_",DUP_id,"_gedi_wk24_",wwfbiom,".csv",sep=""))
#    }
#
#    else if(length(DUP_paf) > 3){
#        print(paste("there are more than 3 files with pa_id = ", DUP_id, sep=""))
#    }
#}


#------------
for(i in 1:length(gedi_paid_UNQ)){

    UNQ_id <- gedi_paid_UNQ[i]
    UNQ_paf <- gedi_paf_ALL[which(gedi_paid_ALL == UNQ_id)]
    wwfbiom <- strsplit(strsplit(UNQ_paf[1], split=c("_conti_biome_"))[[1]][2], split=".RDS")[[1]]
    #print(UNQ_paf)

    if(file.exists(paste(f.path,"/WDPA_extract/pa_stats_ALL/BRA_pa_",UNQ_id,"_gedi_wk24_",wwfbiom,".csv",sep=""))){
        print(paste("csv file already exists for",UNQ_paf,sep=""))
    } else {

    print(UNQ_id)
    pa_metrics <- readRDS(UNQ_paf) %>% unique()
    wwfbiom <- strsplit(strsplit(UNQ_paf, split=c("_conti_biome_"))[[1]][2], split=".RDS")[[1]]
    write.csv(pa_metrics,
        file=paste(f.path,"/WDPA_extract/pa_stats_ALL/BRA_pa_",UNQ_id,"_gedi_wk24_",wwfbiom,".csv",sep=""))
    }
}

print(length(list.files(paste(f.path,"/WDPA_extract/pa_stats_ALL/",sep=""))))






