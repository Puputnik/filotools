#!/usr/bin/Rscript

source(file.path(Sys.getenv("basedir"),"/bin/functions.R"))
#source(file.path("/tank/Filo_utilities/filotools","/bin/functions.R"))

exts <- read.table(file.path(Sys.getenv("basedir"),"/utility/extensions.tsv"), stringsAsFactors = F)
#exts <- read.table(file.path("/tank/Filo_utilities/filotools/","/utility/extensions.tsv"), stringsAsFactors = F)

config_path <- Sys.getenv("config_path")
#config_path <- "/tank/Filo_utilities/filotools/config_files/default.config"
#config_path <- "/tank/USB3/LB_ANALYSIS_2022/BAM_FINAL/final_config.config"

confl <- get_config(config_path)
confl <- set_defaults(confl)

#print(confl)

summ <- c()
for (n in 1:nrow(exts)){
  #print(n)
  ref <- exts[n,1]
  
  if (ref %in% names(confl)){
  ext <- exts[n,2]
  ext <- strsplit(ext,",")[[1]]
  pat <- confl[[ref]]
  sam <- c()
  fil <- c()
  for (e in ext){
    fil_temp <- list.files(pat, paste("*",e,sep=""))
    whi <- which(fil_temp %in% fil)
    if (length(whi) > 0){
      fil_temp <- fil_temp[-whi]
    }
    
    if (length(fil_temp) > 0){
      fil <- c(fil, fil_temp)
      sam <- c(sam, gsub(paste(e,"$",sep=""),"", fil_temp, perl=T))
    }
  }
  if (length(sam) > 0){
    summ <- rbind(summ, cbind(file.path(pat,fil), sam, pat))
  }
  }
}

summ <- data.frame(summ, stringsAsFactors = F)
colnames(summ) <- c("file", "sample",  "path")

samples <- unique(summ$sample)
paths <- unique(summ$path)

c=0
while (c < length(paths)){
  c=c+1
  new_paths <- list.dirs(paths[c], recursive = F)
  subfiles <- setdiff(list.files(paths[c],include.dirs = FALSE, full.names = T),new_paths)
  
  whi <- which(new_paths %in% paths)
  if (length(whi > 0)){
    new_paths <- new_paths[-whi]
  }
  paths <- c(paths, new_paths)
  
  
  for (s in samples){
    if(grepl(s, paths[c])){
      summ <- rbind(summ, c(paths[c], s,"DIR"))
    }
    
    matching <- subfiles[which(grepl(s, subfiles))]
    #matching <- subfiles[which(grepl("asdfa", subfiles))]
    
    for (m in matching){
      #print("ciao")
      summ <- rbind(summ, c(m, s,paths[c]))
    }
  }
}

write.table(summ, quote=F, row.names=F, col.names = F, sep="\t")

#c=120
#c=0
#while (c < length(paths)){
#  c=c+1
#  for (s in samples){
#    if(grepl(s, paths[c])){
#      print(s)
#      sumo <- rbind(summ, c(paths[c], s,"DIR"))
#    }
#  }
#}
#
#paths
