#!/usr/bin/Rscript
##### argument 1: path of the pore occupancy csv, argument 2: path for output pdf (if ignored, same as arg1.pdf)
library(plyr)
library(ggplot2)
library(getopt)

spec = matrix(c(
  'input_path'                    , 'i', 1, "character",
  'output_path'                   , 'o', 1, "character",
  'bins'                          , 'b', 2, "character",
  'fixed_time'                    , 'f', 2, "character",
  'min_time'                      , 'm', 2, "integer",
  'max_time'                      , 'M', 2, "integer",
  'time_format'                   , 't', 2, "character",
  'nolegend'                      , 'L', 2, "character",
  'help'                          , 'h', 0, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

if (!(is.null(opt$help))){
  cat("
        recreate MinKNOW-like duty time plots from  pore_activity_*.csv files
        
        usage:
        filotools dutyplot -i <input_pore_activity.csv> [options]
        
        options:
        
        -o: [CHAR] specify path for output file [default: input_pore_activity.pdf]
        -b: [INT]  specify number of runtime bins [default: 30]
        -f: [INT]  specify a fixed lenght (in minutes) for each runtime bin (overrides -b) [default: null]
        -m: [INT]  specify the start time of the plot (in minutes) [default: 0]
        -M: [INT]  specify the end time of the plot (in minutes) [default: total runtime]
        -L:        suppress legend
        
        ")
  stop()
}

#pored <- "/tank/LIQUID_BIOPSY_MAIN/RUNS/Mol_Cancer/M2/20191023_1139_GA10000_FAL01577_6ba7621a/reports/duty_time.csv"
#pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L02_FAT85244_count1/L02_FAT85244_count1/20220728_1857_MN33722_FAT85244_7c13fdcb/pore_activity_FAT85244_0752585b.csv"
#pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L02_FAT85244_count1/L02_FAT85244_count1/20220728_1857_MN33722_FAT85244_7c13fdcb/pore_activity_FAT85244_0752585b.csv"
#pored <- commandArgs(trailingOnly=TRUE)[1]
#out <- commandArgs(trailingOnly=TRUE)[2]
#bins <- commandArgs(trailingOnly=TRUE)[3]

pored <- opt$input_path
out <- opt$output_path

if (is.null(out)){
  out <- paste(tools::file_path_sans_ext(basename(pored)), ".pdf", sep="")
}
titlename <- tools::file_path_sans_ext(basename(out))

if (is.null(opt$bins)){opt$bins = 30}

pored <- read.table(pored, sep=",", header=T, fill=T)
pore <- pored

orig_labels <- rev(c("strand"                      ,
                 "adapter"                     ,
                 "pore"                        ,
                 "unavailable"                 ,
                 "unblocking"                  ,
                 "no_pore"                     ,
                 "unknown_positive"            ,
                 "multiple"                    ,
                 "saturated"                   ,
                 "unknown_negative"            ,
                 "zero"                        ,
                 "disabled"                    ,
                 "unclassified"                ,
                 "unclassified_following_reset",
                 "pending_manual_reset"        ,
                 "pending_mux_change"          ,
                 "locked"                      ))

target_labels <- rev(
                    c("Strand",
                     "Adapter",
                     "Single pore",
                     "Unavailable",
                     "Active feedback",
                     "No Pore From Scan",
                     "Out of range 2",
                     "Possible multiple",
                     "Saturated",
                     "Out of range 1",
                     "Zero",
                     "Channel Disabled",
                     "unclassified",
                     "unclassified_following_reset",
                     "pending_manual_reset",               
                     "pending_mux_change",                 
                     "Pending Reselection")
                  )

lab_colors <- rev(c("#02ff00", 
              "#ede797",   
              "#07cb01",   
              "#54b8b1",   
              "#a53e97",   
              "#8fc6e7",   
              "#0084a9",   
              "#f47e20",   
              "#333333",   
              "#b5aea7",   
              "#4ca9c3",
              "#ccd5d8",
              "#455456",   
              "#455456",	
              "#455456",   
              "#455456",   
              "#ccd5d8"))  
pore$Channel.State <- mapvalues(pore$Channel.State, orig_labels, target_labels)
pore$Channel.State <- factor(pore$Channel.State, levels=target_labels)

#levels(pore$Channel.State)

missing <- which(!(target_labels %in% pore$Channel.State))
if (length(missing > 0)){
  print(paste("missing labels:", paste(target_labels[missing], collapse=", ")))
  target_labels <- target_labels[-missing]
  lab_colors <- lab_colors[-missing]
}
#missing <- c(4,5,6)

if (!(is.null(opt$max_time))){
  pore <- subset(pore, pore$Experiment.Time..minutes. <= opt$max_time)
}
if (!(is.null(opt$min_time))){
  pore <- subset(pore, pore$Experiment.Time..minutes. >= opt$min_time)
}

if (!(is.null(opt$fixed_time))){opt$bins = round(max(pore$Experiment.Time..minutes)/as.numeric(opt$fixed_time),digits = 0) }

p <- ggplot(pore, aes(fill=Channel.State, y=State.Time..samples., x=Experiment.Time..minutes.)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_discrete(name = "Channel State") + 
      scale_x_binned(
         name = waiver(),
         n.breaks = as.numeric(opt$bins),
         nice.breaks = TRUE)  +
      scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0)) + 
      scale_fill_manual(labels= target_labels, values=lab_colors) + 
      ggtitle(titlename) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

if (!(is.null(opt$nolegend))){p <- p + theme(legend.position="none") + theme(axis.title.y=element_blank(),axis.title.x=element_blank())}
#if (!(is.null(opt$time_format))){suppressMessages(library(scales))
#                             p <- p + scale_x_datetime(labels = date_format("%Y-%m-%d %H")) }

pdf(out, 10,5)  
print(p)
dev.off()
