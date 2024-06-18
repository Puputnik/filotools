set_ignored <- function(opt,sets){
  for (i in sets){
    if (!(is.na(opt[[i]]))){
      cat("ignoring option", i)
      opt[[i]] <- NA
    }
  }
  
}

field <- function(namen,n){
  namesplit <- strsplit(namen,"-")
  return(sapply( namesplit, `[`, n))
}


get_config <- function(config_path){
  conft <- read.table(config_path, stringsAsFactors = F)
  confl <- list()
  for (i in 1:nrow(conft)){
    confl[[conft[i,1]]] <- conft[i,2]
  }
  return(confl)
}

#inherit_args <- function(opt,confl,sets="all",exclude=NULL,checknull=T){
#  #### if checknull F --> overwrite anyway, if T --> overwrite only if null
#  if (sets == "all"){ sets <- names(confl) }
#  if (!(is.null(exclude))){ sets <- sets[which(!(sets %in% exclude))]}
#  
#    for (i in sets){
#      if(!(checknull)){
#        opt[[i]] <- confl[[i]]
#      } else {
#        if (is.null(opt[[i]])){
#          opt[[i]] <- confl[[i]]
#        }
#      }
#    } 
#  return(opt)
#}

inherit_args <- function(target,replacement,sets="all",exclude=NULL,checknull=T){
  #### if checknull F --> overwrite anyway, if T --> overwrite only if null
  if (sets == "all"){ sets <- names(replacement) }
  if (!(is.null(exclude))){ sets <- sets[which(!(sets %in% exclude))]}
  
  for (i in sets){
    if(!(checknull)){
      target[[i]] <- replacement[[i]]
    } else {
      if (is.null(target[[i]]) | is.na(target[[i]])){
        target[[i]] <- replacement[[i]]
      }
    }
  } 
  return(target)
}


check_def <- function(sets){
  nulli <- c()
  for (i in sets){
    if (is.null(opt[[i]]) | is.na(opt[[i]])){
      nulli <- c(nulli, i)
    }
  }
  if(length(nulli) > 0){
    cat("missing mandatory info:")
    print(nulli)
    stop()
  }
}


set_defaults <- function(opt, sets="all"){
  listd  <- list("stats_path"          =c("bam_path","STATS"), 
                 "motif_path"          =c("stats_path","MOTIF"),
                 "motif_count_path"    =c("stats_path","MOTIF_COUNTS"),
                 "readl_path"          =c("stats_path","READLENGTH_COUNTS"),
                 "raw_readl_path"      =c("fastq_path","READL_RAW"),
                 "merged_fastq_path"   =c("fastq_path","MERGED")
  ) 
  
  if (sets == "all"){ sets <- names(listd)}
  
  for (i in sets){
    if (is.null(opt[[i]])){
      opt[[i]] <- file.path(opt[[listd[[i]][1]]],listd[[i]][2])
    } else if (opt[[i]] == "default"){
      opt[[i]] <- file.path(opt[[listd[[i]][1]]],listd[[i]][2])
    }
  }
  return(opt)
}

set_as_default <- function(opt, sets){
  for (i in sets){
    if (is.null(opt[[i]])){
      opt[[i]] <- "default"
    }
  }
}


##### merge and filterings

define_mainl <- function(confl){
  mainl <- list()
  mainl$stats <- list("path" = confl$stats_path, "suffix"=".stats","cols"=c("chr", "start", "end", "stats_readID", "zero", "strand", "flag", "mapQ", "LHclip", "LSclip", "readl", "RSclip", "RHclip", "mate","TLEN"))
  mainl$motif <- list(fun=function(x){x <- x[,-4]; x[,3] <- toupper(x[,3]); return(x)}, "path" = confl$motif_path, "suffix"=".motif","cols"=c("motif_readID", "motif_coord", "motif", "ignore"))
  mainl$time <- list("path" = file.path(confl$stats_path, "SEQ_TIME"), "suffix"=".time","cols"=c("time_readID", "time", "run"))
  mainl$barcodes <- list("path" = file.path(confl$stats_path, "BOTH_BARCODES"), "suffix"=".idf","cols"=c("barcodes_readID", "status"))
  mainl$mod <- list("path" = file.path(confl$stats_path, "MOD"), "suffix"=".5mC","cols"=c("mod_readID", "mod_chr", "mod_pos", "CGn", "CGmod", "Cn", "mod_pos_Cmod", "mod_pos_allCG", "modvalues_CG" ))
  return(mainl)
}

merge_subfiles <- function(sets, sample){
  tot <- c()
  for (s in sets){
    load <- as.data.frame(vroom(paste(mainl[[s]]$path,"/", sample, mainl[[s]]$suffix, sep=""), col_names = F, show_col_types = F, comment="#"))
    colnames(load) <- mainl[[s]]$cols
    
    if (!(is.null(mainl[[s]]$fun))){
      load <- mainl[[s]]$fun(load)
    }
    
    if (which(sets == s) == 1){
      tot <- load
    } else {
      tot <- cbind(tot, load)
    }
    if (!( all(tot[,paste(sets[1],"readID",sep="_")] == tot[,paste(s,"readID",sep="_")]) )){print("IDS NOT MATCHING!")} # else {print("samba")}
  }
  return(tot)
}




noclip5 <- expression((tot$strand == "+" & tot$LHclip == 0 & tot$LSclip == 0) | (tot$strand == "-" & tot$RHclip == 0 & tot$RSclip == 0))
noclip <- expression(tot$LHclip == 0 & tot$LSclip == 0 & tot$strand == "-" & tot$RHclip == 0 & tot$RSclip == 0)


##### misc

show_merged <- function(files, confl=confl){
  #files <- tools::file_path_sans_ext(files)
  files <- gsub("\\..*$", "", files, perl=T)
  for (f in files){
    logo <- file.path(confl$merged_fastq_path,"logs_merge",paste(f,".log",sep=""))
    if (file.exists(logo)){
      print(read.table(logo, stringsAsFactors = F, header=F))
    } else {
      print(paste(f,"not merged"))
    }
  }
}

list.filo <- function(path = ".", pattern = NULL, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE, ext=FALSE){
  filinames <- list.files(path = path, pattern = pattern, all.files = all.files, full.names = FALSE, recursive = recursive, ignore.case = ignore.case, include.dirs = include.dirs, no.. = no..)
  fili <- list.files(path = path, pattern = pattern, all.files = all.files, full.names = TRUE, recursive = recursive, ignore.case = ignore.case, include.dirs = include.dirs, no.. = no..)
  if (!(ext)){
    filinames <- gsub("\\..*$", "", filinames, perl=T)
  }
  names(fili) <- filinames
  return(fili)
}



#set_defaults <- function(opt){
#  
#  if (is.null(opt$stats_path) | (opt$stats_path == "default")){
#    opt$stats_path <- file.path(opt$bam_path, "STATS")
#  }
#  
#  if (is.null(opt$motif_path) | (opt$motif_path == "default")){
#    opt$motif_path <- file.path(opt$stats_path,"MOTIF")
#    
#  }
#  
#  if (is.null(opt$motif_count_path) | (opt$motif_count_path == "default")){
#    opt$motif_count_path <- file.path(opt$stats_path,"MOTIF_COUNTS")
#    
#  }
#  
#  if (is.null(opt$readl_path) | (opt$readl_path == "default")){
#    opt$readl_path <- file.path(opt$stats_path,"READLENGTH_COUNTS")
#  }
#  
#  if (is.null(opt$raw_readl_path) | (opt$raw_readl_path == "default")){
#    opt$raw_readl_path <- file.path(opt$fastq_path, "READL_RAW")
#  }
#  
#  if (is.null(opt$merged_fastq_path) | (opt$merged_fastq_path == "default")){
#    opt$merged_fastq_path <- file.path(opt$fastq_path,"MERGED")
#  }
#  return(opt)
#}