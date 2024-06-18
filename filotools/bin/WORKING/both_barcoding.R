#!/usr/bin/Rscript
library(getopt)
library(foreach)
library(doParallel)
library(vroom)
library(data.table)

setwd("/tank/USB3/LB_ANALYSIS_2022/BAM/")
opt <- c()
opt$pattern <- "*1.bam"
opt$filter_tresholds <- "60:70"

spec = matrix(c(
  'pattern'            , 'p', 1, "character",
  'output_dir'         , 'o', 2, "character",
  'summary'            , 's', 2, "character",
  'barcode'            , 'b', 2, "character",
  'read_ids'           , 'i', 2, "character",
  'filter_tresholds'   , 'f', 2, "character",
  'merged'             , 'm', 0, "double",
  
  'threads'            , 't', 2, "integer",
  'help'               , 'h', 0, "double"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

patt <- opt$pattern




if (grepl(".bam$", patt, perl=T)){
  setwd(file.path(getwd(), "STATS"))
  patt <- gsub(".bam$", ".stats", patt, perl=T)
}

fili <- list.files(getwd(), patter=patt)

if ( is.null(opt$output_dir          ) ) { } 
if ( is.null(opt$summary             ) ) {lib_path <- "/tank/LIQUID_BIOPSY_GRANT/MEGALODON/" } 
if ( is.null(opt$filter_tresholds    ) ) {thr <- c(60,60)} else { thr <- as.numeric(strsplit(opt$filter_tresholds, ":")[[1]]) }
if ( is.null(opt$output_dir          ) ) { opt$output_dir <- file.path(getwd(), "BOTH_BARCODES") } 
if ( is.null(opt$threads             ) ) { opt$threads  <- "24"} 

#'output_dir'         , 'o', 2, "character",
#'summary'            , 's', 2, "character",
#'barcode'            , 'b', 2, "character",
#'read_ids'           , 'i', 2, "character",
#'filter_tresholds'   , 'f', 2, "character",
#'merged'             , 'm', 0, "double",
#'threads'            , 't', 2, "integer",


slist <- list()
basenames <- gsub(".stats$", "", fili, perl=T)
for (ba in basenames){
  
}

namesplit <- strsplit(basenames, "-")
field <- function(n){
  sapply( namesplit, `[`, n)
}
#samples_info <- cbind(fili,basenames,paste(field(5), field(6), paste("count",field(7), sep=""), sep="_"), field(4))
#colnames(samples_info) <- c("file", "basename", "library", "barcode")
libraries <- paste(field(5), field(6), paste("count",field(7), sep=""), sep="_")
barcodes <- substring(field(4), 2)

lib_list <- list()

#for (l in unique(samples_info[,"library"])){
for (l in unique(libraries)){
  lib_list[[l]] <- vroom(file.path(lib_path,l, "GUPPY_DEM_EITHER/barcoding_summary.txt"))
}

registerDoParallel(numCores)  # use multicore, set to the number of our cores
foreach (i=1:nrow(samples_info)) %dopar% {
  l <- libraries[i]
  if (!(is.null(opt$barcode          ))) {thr <- c(60,60)}
  
  summa <- lib_list[[l]][which(lib_list[[l]]$barcode_arrangement == paste("barcode",barcodes[i],sep="") &
                                 lib_list[[l]]$barcode_front_score > thr[1] & 
                                 lib_list[[l]]$barcode_rear_score > thr[2]    ),]
  
  stats <- vroom(file.path(getwd(), fili[i]), col_names = F)
  write.table(cbind(stats$X4, as.numeric(stats$X4 %chin% summa$read_id)), file=file.path(outdir, paste(basenames[i],".idf",sep="")))
  ?write.table()
}








sandro$ciano[which(sandro$ciano$X1 > 1),]




registerDoParallel(numCores)  # use multicore, set to the number of our cores
foreach (i=1:3) %dopar% {
  sqrt(i)
}

#print(opt)


ciano <- data.frame(matrix(seq(1,6), nrow=3))
cieno <- data.frame(matrix(seq(1,8), nrow=4))

sandro <- list(ciano, cieno)

names(sandro) <- c("ciano", "cieno")





#if(x < 500) { stop("flag -s should be used in combination with -b")}

## set some reasonable defaults for the options that are needed,
## but were not specified.
#if ( is.null(opt$mean    ) ) { opt$mean    = 0     }
#if ( is.null(opt$sd      ) ) { opt$sd      = 1     }
#if ( is.null(opt$count   ) ) { opt$count   = 10    }
#if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }
#
## print some progress messages to stderr, if requested.
#if ( opt$verbose ) { write("writing...",stderr()) }