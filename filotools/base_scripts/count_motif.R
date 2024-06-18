#!/usr/bin/Rscript
library(vroom)

chr_list <- commandArgs(trailingOnly=TRUE)[1]
stats_path <-    commandArgs(trailingOnly=TRUE)[2]
motif_path <- commandArgs(trailingOnly=TRUE)[3]
out_path <-      commandArgs(trailingOnly=TRUE)[4]
sample_name <-   commandArgs(trailingOnly=TRUE)[5]
keep_S <-   as.logical(commandArgs(trailingOnly=TRUE)[6])
keep_H <-   as.logical(commandArgs(trailingOnly=TRUE)[7])
mapq <-   commandArgs(trailingOnly=TRUE)[8]
id_path <-   commandArgs(trailingOnly=TRUE)[9]
max_length <-   commandArgs(trailingOnly=TRUE)[10]

debug <- c()

stats_path <- file.path(stats_path,paste(sample_name,".stats", sep=""))
motif_path <- file.path(motif_path,paste(sample_name,".motif", sep=""))

#### loading chromosome list
chr <- read.table(chr_list, stringsAsFactors = F)

#### loading stat files
stats <- as.data.frame(vroom(stats_path, show_col_types = FALSE, col_names = F))
colnames(stats) <- paste("V", seq(1,length(colnames(stats))), sep="")
debug <- c(debug, paste("Stats total reads:", length(rownames(stats))))

#### loading motif files
motif <- as.data.frame(vroom(motif_path, show_col_types = FALSE, col_names = F, delim=" "))
motif$X3 <- toupper(motif$X3)
stats <- cbind(stats, motif)
debug <- c(debug, paste("Read names corresponding:", all(stats$V4 == stats$X1)))

#### filtering by ids
if (id_path != "NO"){
  both <- as.data.frame(vroom(file.path(id_path, paste(sample_name,".idf", sep="")), delim="\t", show_col_types = FALSE, col_names = F))
  
  if (length(colnames(both)) == 1){
    colnames(both) <- "V4"
    stats <- merge(stats, both, by="V4", all=F)
  } else {
    colnames(both) <- c("id_filter", "filter")
    stats <- cbind(stats, both)
    stats <- subset(stats, stats$filter == 1)
  }
}


stats <- subset(stats, stats$V8 > mapq & stats$V1 %in% chr[,1] & stats$V11 < max_length)

if (!(keep_S)){
  stats <- subset(stats, (stats$V6 == "+" & stats$V10 == 0) | (stats$V6 == "-" & stats$V12 == 0)) 
}


if (!(keep_H)){
  stats <- subset(stats, (stats$V6 == "+" & stats$V9 == 0) | (stats$V6 == "-" & stats$V13 == 0)) 
}


##### keeping reads mapped on selected chromosomes, with MAPQ > 20, readlength < 700 and without 5' hard/soft clipping
#stats <- subset(stats, stats$V8 > 20 & stats$V1 %in% chr[,1] & stats$V11 < 700 & ((stats$V6 == "+" & stats$V9 == 0 & stats$V10 == 0) | (stats$V6 == "-" & stats$V12 == 0 & stats$V13 == 0)))
debug <- c(debug, paste("Stats filtered reads:", length(rownames(stats))))

####obtaining 4bp motif counts
counts <- table(stats$X3)

if (!(file.exists(out_path))){
  dir.create(out_path)
}

saveRDS(counts, file.path(out_path,paste(sample_name,".motif.R", sep="")))
write.table(debug, file.path(out_path,paste(sample_name,".motif.log", sep="")), row.names = F, col.names = F, quote = F)

