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
  'nolegend'                      , 'L', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)



#pored <- "/tank/LIQUID_BIOPSY_MAIN/RUNS/Mol_Cancer/M2/20191023_1139_GA10000_FAL01577_6ba7621a/reports/duty_time.csv"
#pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L02_FAT85244_count1/L02_FAT85244_count1/20220728_1857_MN33722_FAT85244_7c13fdcb/pore_activity_FAT85244_0752585b.csv"
pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L02_FAT85244_count1/L02_FAT85244_count1/20220728_1857_MN33722_FAT85244_7c13fdcb/pore_activity_FAT85244_0752585b.csv"
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


pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L08_FAS49889_count0/L08_FAS49889_count0/20221013_1130_MN33722_FAS49889_800ee41c/pore_activity_FAS49889_c1f8fe6b.csv"
pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L09_FAS49797_count0/L09_FAS49797_count0/20221013_1611_MN33722_FAS49797_49e888e0/pore_activity_FAS49797_ca62b1ec.csv"

pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L04_FAT85247_count0/L03_FAT85247_count0/20220825_1543_MN33722_FAT85247_1e5dd683/pore_activity_FAT85247_b99b7964.csv"
pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L05_FAT85247_count0/L05_FAT85247_count0/20220826_1817_MN33722_FAT85247_1c64a32d/pore_activity_FAT85247_0177f118.csv"
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

unique(pore$Channel.State)

canali <- c("Strand",
            "Adapter",
            #"Single pore",
            "Unavailable",
            "Active feedback")
#canali <- c("Strand",
#            "Adapter")


tot <- c()
co=0
for (c in canali){
  co=co+1
  pare <- subset(pore, pore$Channel.State == c)
  if (co==1){
    tot <- pare[,-1]
    colnames(tot)[2] <- c
  } else {
    if (all(pare$Experiment.Time..minutes. == tot$Experiment.Time..minutes.)){
      colnames(pare)[3] <- c
      pare <- pare[,3]
      tot <- cbind(tot, pare)
    } else {
      print("diobello")
    }
  }
}
colnames(tot) <- c("time",canali)
for (i in rownames(tot)){
  tot[i,"ratio_all"] <- tot[i,canali[1]]/sum(tot[i,canali])
}
tot <- subset(tot, tot$time <= 180)
L08 <- tot
totals <- sum(L08$Strand)/sum(sum(L08$Strand),sum(L08$Adapter), sum(L08$Unavailable), sum(L08$`Active feedback`) )
totals
L09 <- tot

L09 <- subset(L09, L09$time <= 180)

LTOT <- cbind(L09[,c("time","ratio_all")],run="0.5X+1.2X")
LTOT <- rbind(LTOT, cbind(L08[,c("time","ratio_all")],run="0.5X"))

#LTOT <- cbind(L09[,c(1,4)],run="0.5X+1.2X")
#LTOT <- cbind(L09[,c(1,6)],run="0.5X+1.2X")
#LTOT <- rbind(LTOT, cbind(L08[,c(1,4)],run="0.5X"))
#LTOT <- rbind(LTOT, cbind(L08[,c(1,6)],run="0.5X"))
LTOT$ratio_norm <- LTOT$ratio_all/totals


LTOT <- cbind(L09[,c("time","ratio_all")],run="1.6X+1.2X")
LTOT <- rbind(LTOT, cbind(L08[,c("time","ratio_all")],run="1.6X"))

p <- ggplot(LTOT, aes(x=time, y=ratio_all, col=run)) + geom_line()
p <- ggplot(LTOT, aes(x=time, y=ratio_norm, col=run)) + geom_line() + ylim(0,1.6)
pdf("~/LOW_INPUT.pdf",7,4.5)
png("~/05_norm_ylim.png",700,450)
png("~/16_norm.png",700,450)
png("~/HIG_INPUT.png",700,450)
print(p)
dev.off()

i=1
opt$max_time = 180
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
