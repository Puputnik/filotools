#!/usr/bin/Rscript
library(optparse)
library(foreach)
library(doParallel)
library(vroom)
library(data.table)
library(tools)

#source(file.path("/tank/Filo_utilities/filotools","/bin/functions.R"))
source(file.path(Sys.getenv("basedir"),"/bin/functions.R"))

print(commandArgs(trailingOnly=TRUE))


config_path <- Sys.getenv("config_path")
#config_path <- "/tank/Filo_utilities/filotools/config_files/default.config"
#config_path <- "/tank/USB3/LB_ANALYSIS_2022/BAM_FINAL/final_config.config"
confl <- get_config(config_path)
confl <- set_defaults(confl)

option_list = list(
  make_option(c("-o", "--output_dir") ,        action="store",        default=NA,      type='character', help="path to output folder (mandatory)"),
  make_option(c("-f", "--filter_tresholds"),   action="store",        default="60:60", type="character", help="front and rear threshold for barcode detection, fomatted as <front_thr>:<rear_thr> [default: 60:60]"),
  
  make_option(c("-s", "--samplesheet"),        action="store",        default=NA,      type='character', help="path to tab separated file containing <sample_id_without_extension>\t<BXX>\t<path_to_barcoding_summary.txt> where XX is the two digits number of the barcode. if a sample is a result of a merging of samples from different runs, provide each barcode/barcoding summary pair in a separate line referencing to the same sample_id" ),
  make_option(c("-i", "--stats_path"),         action="store",        default=NA,      type='character', help="path to folder containing .stats files (mandatory if using samplesheet, can be inherited from config_file using -C)"),
  
  
  make_option(c("-p", '--pattern'           ), action="store",        default="*.stats",      type="character", help= "pattern to retrieve samples from working directory. IMPORTANT: launch this tool from within the folder containing .bam or .stats files (used only if samplesheet is not provided). [default= '*.stats']"),
  make_option(c("-F", '--fastq_path'        ), action="store",        default=NA,      type="character", help= "path in which fastq files are stored (used only if samplesheet is not provided)"),
  make_option(c("-M", '--merged_fastq_path' ), action="store",        default=NA,      type="character", help= "path in which merged fastq files are stored (used only if samplesheet is not provided)"),
  make_option(c("-B", '--barcoding_path'    ), action="store",        default=NA,      type="character", help= "path in which barcoded results are stored. this path should include one subfolder for each run with filo's naming layout LXX_FATXXXXX_countX. (used only if samplesheet is not provided)"),
  make_option(c("-S", '--barcoding_suffix'    ), action="store",      default="GUPPY_DEM_EITHER/barcoding_summary.txt", type="character", help="path from barcoding_path subfolders to barcoding_summary.txt file [default: GUPPY_DEM_EITHER/barcoding_summary.txt]. (used only if samplesheet is not provided)"),
  
  make_option(c("-C", '--force_config'    ),   action="store_true",   default=F,       type="character", help="forces the use of fastq_path/barcoding_path/merged_fastq_path (only if samplesheet is not provided) and stats_path (always) from config file " )
  #make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show feets")
  
)

usage_string <- "creates .ids files with read ids (first column) and info about barcoding status (second column): 0=single barcoded, 1=double barcoded.

if using samplesheet

filotools bothbar -i <stats_path> -s <samplesheet_path> -o <output_dir>


if not using samplesheet
IMPORTANT: all file/library names must follow Filo's layout. Launch this tool from within the folder containing .bam or .stats files

filotools bothbar -p '*.stats' -F <fastq_path> -M <merged_fastq_path> -B <barcoding_path> -o <output_dir>
filotools bothbar -p '*.stats' -C -o <output_dir>


"

opt = parse_args(OptionParser(option_list=option_list, 
                              usage=usage_string,
                              add_help_option = T))

thr <- strsplit(opt$filter_tresholds,":")[[1]]

#opt$samplesheet <- "/tank/USB3/LB_ANALYSIS_2022/BAM_FINAL/samplesheet.tsv"
#opt$samplesheet <- NA
#opt$pattern <- "*.stats"

sample_list <- list()

if (!is.na(opt$samplesheet)){
  
  if (opt$force_config) {opt <- inherit_args(opt, confl,sets=c("stats_path"), checknull = T)}
  check_def(c("stats_path"))
  #set_ignored(c("pattern","fastq_path","barcoding_path","merged_fastq_path","barcoding_suffix"))
  
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
  print("confl")
  print(confl)
  if (opt$force_config) {opt <- inherit_args(opt, confl,sets=c("fastq_path","barcoding_path"), checknull = T)}
  print("opt")
  print(opt)
  check_def(c("pattern","fastq_path","barcoding_path"))
  sname_tot <- c()
  for (f in list.files(getwd(),opt$pattern)){
    sname_tot <- c(sname_tot, file_path_as_absolute(f))
  }
  if (is.null(sname_tot)){
    #print()
    stop("no input files found. check your pattern")
  }
  filename <- basename(sname_tot)
  exts <- file_ext(filename)
  dirs <- dirname(sname_tot)
  
  if (length(unique(exts)) > 1 | length(unique(dirs)) > 1) {
    cat("too many extensions or input directories!!!") 
    cat("extensions:",unique(exts))
    cat("directories:",unique(dirs))
    stop()
  }
  
  
  if (opt$force_config) {opt <- inherit_args(opt, confl,sets=c("stats_path"), checknull = T)}
  if (is.na(opt$stats_path)){
    #### infer!
    if (grepl(".bam$", unique(exts), perl=T)){
      opt$stats_path <- file.path(unique(dirs),"STATS")
    } else if (grepl(".stats$", unique(exts), perl=T)){
      opt$stats_path <- file.path(unique(dirs))
    } else {
      cat("unknown ext")
      stop()
    }
  } 
  check_def(c("stats_path"))
  
  basenames <- file_path_sans_ext(filename)
  
  libraries <- c()
  for (i in basenames){
    sample_list[[i]] <- list()
    
    if (any(grepl(i,list.files(opt$fastq_path)))){
      ##### infer dal singolo usando field
      temp_list <- list(field(i,4))
      names(temp_list) <- file.path(opt$barcoding_path, paste(field(i,5),"_",field(i,6),"_count",field(i,7),sep=""), opt$barcoding_suffix)
      sample_list[[i]] <- append(sample_list[[i]], temp_list)
      
      libraries <- c(libraries,names(temp_list))
    
    } else {
      if (opt$force_config) {opt <- inherit_args(opt, confl,sets=c("merged_fastq_path"), checknull = T)}
      check_def("merged_fastq_path")
      if (any(grepl(i,list.files(opt$merged_fastq_path)))){
        logg <- read.table(file.path(opt$merged_fastq_path,"logs_merge",paste(i,".log", sep="")), stringsAsFactors = F)
        if (logg$V5[1] != "Complete"){stop(cat("merge sample",i,"incomplete!! exiting"))}
        
        sub_files <- strsplit(logg$V4[1],",")[[1]]
        for (s in sub_files){
          temp_list <- list(field(s,4))
          names(temp_list) <- file.path(opt$barcoding_path, paste(field(s,5),"_",field(s,6),"_count",field(s,7),sep=""), opt$barcoding_suffix)
          sample_list[[i]] <- append(sample_list[[i]], temp_list)
          libraries <- c(libraries,names(temp_list))
        }
      }
    }
  }
}

cat("loading summaries:\n")
##### load summaries
lib_list <- list()
for (l in unique(libraries)){
  cat("loading",l,"\n")
  temp_summ <- vroom(l,show_col_types = FALSE)
  temp_summ <- temp_summ[which(temp_summ$barcode_front_score > thr[1] & temp_summ$barcode_rear_score >  thr[2]),]
  #lib_list[[l]] <- vroom(l)
  lib_list[[l]] <- temp_summ
  temp_summ <- NULL
}

if (is.na(opt$output_dir)){opt$output_dir <- file.path(opt$stats_path,"BOTH_BARCODES")}

cat("filtering samples:\n")
##### perform filtering
for (i in basenames){
  cat("filtering", i,"\n")
  summa <- c()
  for(l in names(sample_list[[i]])){
    summa <- rbind(summa, lib_list[[l]][which(lib_list[[l]]$barcode_arrangement %chin% paste("barcode", substring(sample_list[[i]][[l]], 2),sep="")),])
  }
  stats <- vroom(file.path(opt$stats_path, paste(i,".stats",sep="")), col_names = F, show_col_types = FALSE)
  dir.create(opt$output_dir, showWarnings = F)
  write.table(cbind(stats$X4, as.numeric(stats$X4 %chin% summa$read_id)), file=file.path(opt$output_dir, paste(i,".idf",sep="")), quote=F, row.names = F, col.names = F)
}


#all(vroom("/tank/USB3/LB_ANALYSIS_2022/BAM_FINAL/STATS/BOTH_BARCODES/ELB04-T00-1121-B04-L01-FAS34889-A.idf") == vroom("/tank/USB3/LB_ANALYSIS_2022/BAM_FINAL/STATS/BOTH_BARCODES/sasone"))
