#!/usr/bin/Rscript
##### argument 1: path of the pore occupancy csv, argument 2: path for output pdf (if ignored, same as arg1.pdf)
library(plyr)
library(ggplot2)

#pored <- "/tank/LIQUID_BIOPSY_MAIN/RUNS/Mol_Cancer/M2/20191023_1139_GA10000_FAL01577_6ba7621a/reports/duty_time.csv"
#pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L02_FAT85244_count1/L02_FAT85244_count1/20220728_1857_MN33722_FAT85244_7c13fdcb/pore_activity_FAT85244_0752585b.csv"
pored <- commandArgs(trailingOnly=TRUE)[1]
out <- commandArgs(trailingOnly=TRUE)[2]
#bins <- commandArgs(trailingOnly=TRUE)[3]

if (is.na(out)){
  out <- paste(tools::file_path_sans_ext(basename(pored)), ".pdf", sep="")
}
titlename <- tools::file_path_sans_ext(basename(out))

#if (is.na(bins)){
#  bins = 30
#}

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

levels(pore$Channel.State)

missing <- which(!(target_labels %in% pore$Channel.State))
if (length(missing > 0)){
  print(paste("missing labels:", paste(target_labels[missing], collapse=", ")))
  target_labels <- target_labels[-missing]
  lab_colors <- lab_colors[-missing]
}
#missing <- c(4,5,6)
p <- ggplot(pore, aes(fill=Channel.State, y=State.Time..samples., x=Experiment.Time..minutes.)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_discrete(name = "Channel State") + 
      scale_x_binned(
         name = waiver(),
         n.breaks = 30,
         nice.breaks = TRUE)  +
      scale_fill_manual(labels= target_labels, values=lab_colors) + 
      ggtitle(titlename) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(out)  
print(p)
dev.off()
