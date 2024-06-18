#!/usr/bin/Rscript
library(vroom)

#chr_list <-   "/tank/Filo_utilities/filotools/utility/chr_list.txt"
#stats_path <- "/tank/USB3/LB_ANALYSIS_2022/BAM/STATS"
#out_path <-   "/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/READLENGTH_COUNTS/"
#sample_name <- "PLB03-T00-1C00-B01-L02-FAT85244-0"
#keep_S <- F
#keep_H <- F
#mapq = 20
#id_path <-   "/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/BOTH_BARCODES/"

chr_list <- commandArgs(trailingOnly=TRUE)[1]
stats_path <-    commandArgs(trailingOnly=TRUE)[2]
out_path <-      commandArgs(trailingOnly=TRUE)[3]
sample_name <-   commandArgs(trailingOnly=TRUE)[4]
keep_S <-   as.logical(commandArgs(trailingOnly=TRUE)[5])
keep_H <-   as.logical(commandArgs(trailingOnly=TRUE)[6])
mapq <-   commandArgs(trailingOnly=TRUE)[7]
id_path <-   commandArgs(trailingOnly=TRUE)[8]

stats_path <- file.path(stats_path,paste(sample_name,".stats", sep=""))

#### loading chromosome list
chr <- read.table(chr_list, stringsAsFactors = F)

#### loading stat files
stats <- as.data.frame(vroom(stats_path, show_col_types = FALSE, col_names = F))
colnames(stats) <- paste("V", seq(1,length(colnames(stats))), sep="")

if (id_path != "NO"){
  both <- as.data.frame(vroom(file.path(id_path, paste(sample_name,".idf", sep="")), delim=" ", show_col_types = FALSE, col_names = F))
  
  if (length(colnames(both)) == 1){
    colnames(both) <- "V4"
    stats <- merge(stats, both, by="V4", all=F)
  } else {
    colnames(both) <- c("id_filter", "filter")
    stats <- cbind(stats, both)
    if (all(stats$id_filter == stats$V4)){
      stats <- subset(stats, stats$filter == 1)
    } else {
      print(paste("FATAL ERROR: read IDs for sample",sample_name, "not corresponding to fildering IDs"))
    }
  }
}

#### keeping reads mapped on selected chromosomes, with MAPQ > 20, readlength < 700 and without 5'/3' hard/soft clipping
stats <- subset(stats, stats$V8 > mapq & stats$V1 %in% chr[,1])

if (!(keep_S)){
  stats <- subset(stats, stats$V10 == 0  & stats$V12 == 0) 
}

if (!(keep_H)){
  stats <- subset(stats, stats$V9 == 0  & stats$V13 == 0) 
}


if (!(file.exists(out_path))){
  dir.create(out_path)
}

#### creating count table for each readlength value
length_counts <- table(stats$V11)
saveRDS(length_counts, file.path(out_path,paste(sample_name,".readlength_counts.R", sep="")))
