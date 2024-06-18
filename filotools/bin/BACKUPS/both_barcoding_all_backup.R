#!/usr/bin/Rscript
usage <- 
"calculates read counts for each readlength bin (1bp) from aligments (.stat files) or .fastq files
usage: filotools readl [options] <input_file(s)> 
  can also glob, example:
  filotool readl *.bam 
       
	options:
  -o: specify output folder [default = <stats_input_folder>/MOTIF_COUNTS]
-t: [INT] number of threads for multiprocessing (default: 4)
-s: [CHAR] subset the dataset based on read ids 

filtering options (only for bam/stats)
-l: [CHAR] specify chr_list files (default: config)
-q: [INT] specify minimum mapping quality (default: 20)
-k: [INT] specify max readlength (default: 700)
-H: don't filter out hard clipped reads
		-S: don't filter out soft clipped reads

Config file is used by default for:
  
  chr_list

-O: force output folder as readl_path (for bam/stats) and raw_readl_path (for fastq/fastq.gz)
-I: force input  folder as stats_path (for bam/stats) and fastq_path (for fastq/fastq.gz)
-C: combination of -O and -I
"


library(getopt)
library(foreach)
library(doParallel)
library(vroom)
library(data.table)

source(file.path("/tank/Filo_utilities/filotools","/bin/functions.R"))

#source(file.path(Sys.getenv("basedir"),"/bin/functions.R"))


#config_path <- Sys.getenv("config_path")
#config_path <- "/tank/Filo_utilities/filotools/config_files/default.config"
config_path <- "/tank/USB3/LB_ANALYSIS_2022/BAM_FINAL/final_config.config"
confl <- get_config(config_path)
confl <- set_defaults(confl)


#setwd("/tank/USB3/LB_ANALYSIS_2022/BAM_FINAL")
#opt <- c()
#opt$pattern <- "*.bam"
#opt$filter_tresholds <- "60:60"
#opt$samplesheet <- "/tank/USB3/LB_ANALYSIS_2022/BAM_FINAL/samplesheet.tsv"
#opt$fastq_path <- 
#opt$barcoding_path <-


spec = matrix(c(
  'output_dir'         , 'o', 2, "character",
  'samplesheet'        , 's', 2, "character",
  'filter_tresholds'   , 'f', 2, "character",
  'stats_path'         , 'i', 2, "character", #### facoltativo altrimenti si sfrutta il config
  
  #### questi sono obbligatori solo se samplesheet non è fornito
  'pattern'            , 'p', 2, "character",
  'fastq_path'         , 'F', 2, "character",
  'merged_fastq_path'  , 'M', 2, "character",
  'barcoding_path'     , 'B', 2, "character",
  
  'force_config'       , 'C', 0, "double",
  
  'threads'            , 't', 2, "integer",
  'help'               , 'h', 0, "double"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

opt <- inherit_config(opt,confl,sets="all",checknull=T)


print(opt)

if ( !is.null(opt$help) ) {
  cat()
  q(status=1)
}


sample_list <- list()

if (!is.null(opt$samplesheet)){
  
  check_def(c("stats_path"))
  
  samplesheet <- read.table(opt$samplesheet, stringsAsFactors = F)
  colnames(samplesheet) <- c("sname", "bc", "summ")
  basenames <- unique(samplesheet$sname)
  
  for (i in basenames){
    sample_list[[i]] <- list()
  }
  
  for (i in 1:nrow(samplesheet)){
    temp_list <- list(samplesheet[i,"bc"])
    names(temp_list) <- samplesheet[i,"summ"]
    sample_list[[samplesheet[i,"sname"]]] <- append(sample_list[[samplesheet[i,"sname"]]], temp_list)
  }
  
  temp_list <- NULL
  libraries <- unique(samplesheet$summ)

} else {
  check_def(c("pattern","fastq_path","barcoding_path"))
  #### se non esiste il campioni nel fastq path allora controlla anche ,"merged_fastq_path"
}

print("passato")

#list(samplesheet[i,"summ"] = samplesheet[i,"bc"])
#list("drago" = samplesheet[i,"bc"])
#
#sample_list[[1]]
#sample_list[[1]]
#
#as.character(samplesheet[i,"sname"])
#
#
#
#
#
#
#
#
#opt[[c("fastq_path","barcoding_path")]]
#opt[[c("barcoding_path")]]
#
#
#
#
#
#
#
#
#
#
#if ( !is.null(opt$help) ) {
#  #cat(getopt(spec, usage=TRUE))
#  q(status=1)
#}
#
#if (grepl(".bam$", opt$pattern, perl=T)){
#  setwd(file.path(getwd(), "STATS"))
#  opt$pattern <- gsub(".bam$", ".stats", opt$pattern, perl=T)
#}
#
#opt <- inherit_config(opt,confl,sets="barcoding_path",checknull=T)
#
#if ( is.null(opt$filter_tresholds    ) ) {thr <- c(60,60)} else { thr <- as.numeric(strsplit(opt$filter_tresholds, ":")[[1]]) }
#if ( is.null(opt$output_dir          ) ) { opt$output_dir <- file.path(getwd(), "BOTH_BARCODES") } 
#if ( is.null(opt$threads             ) ) { opt$threads  <- "24"} 
#
#is.null(opt[["summary"]]             )
#slist <- list()
#
#fili <- list.files(getwd(), pattern=opt$pattern)
#basenames <- gsub(".stats$", "", fili, perl=T)
#
##namesplit <- strsplit(basenames, "-")
##
##field <- function(n){
##  sapply( namesplit, `[`, n)
##}
#
#if ( is.null(opt$summary) ){
#  for (ba in basenames){
#    
#  }
#  libraries <- paste(field(5), field(6), paste("count",field(7), sep=""), sep="_")
#  
#  lib_list <- list()
#  for (l in unique(libraries)){
#    lib_list[[l]] <- vroom(file.path(lib_path,l, "GUPPY_DEM_EITHER/barcoding_summary.txt"))
#  }
#  
#} else {
#  libraries <- rep("user_defined", length(fili))
#  lib_list[["user_defined"]] <- vroom(opt$summary)
#}
#
##if ( is.null(opt$summary) ){
##  libraries <- paste(field(5), field(6), paste("count",field(7), sep=""), sep="_")
##  
##  lib_list <- list()
##  for (l in unique(libraries)){
##    lib_list[[l]] <- vroom(file.path(lib_path,l, "GUPPY_DEM_EITHER/barcoding_summary.txt"))
##  }
##  
##} else {
##  libraries <- rep("user_defined", length(fili))
##  lib_list[["user_defined"]] <- vroom(opt$summary)
##}
#
#
#
#
#barcodes <- substring(field(4), 2)
#
#
#
##cl <-makeCluster(opt$threads, type ="PSOCK")
#registerDoParallel(opt$threads)  # use multicore, set to the number of our cores
##foreach (i=1:nrow(samples_info)) %dopar% {
##foreach (i=1:length(fili)) %dopar% {
#for (i in seq(1,length(fili))) {
#  l <- libraries[i]
#  summa <- lib_list[[l]][which(lib_list[[l]]$barcode_arrangement == paste("barcode",barcodes[i],sep="") &
#                                 lib_list[[l]]$barcode_front_score > thr[1] & 
#                                 lib_list[[l]]$barcode_rear_score > thr[2]    ),]
#  
#  stats <- vroom(file.path(getwd(), fili[i]), col_names = F)
#  dir.create(opt$output_dir, showWarnings = F)
#  write.table(cbind(stats$X4, as.numeric(stats$X4 %chin% summa$read_id)), file=file.path(opt$output_dir, paste(basenames[i],".idf",sep="")), quote=F, row.names = F, col.names = F)
#}
#
#
##aio <- read.table(file=file.path(opt$output_dir, paste(basenames[i],".idf",sep="")))
##?write.table()