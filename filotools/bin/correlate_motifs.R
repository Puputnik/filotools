#!/usr/bin/Rscript
suppressMessages(library(getopt))
library(ggplot2)
library(tools)
suppressMessages(library(tidyverse))
#opt <- c()
#setwd("/tank/USB3/LB_ANALYSIS_2022/BAM")
#opt$Xaxis <- "ULB10-T00-F100-B02-L03-FAT85247-0.bam"
#opt$Yaxis <- "ULB10-T00-F100-B02-L04-FAT85247-0.bam"
#print("boh")
#stop()

#print("ciao")

spec = matrix(c(
  'Xaxis'                         , 'x', 1, "character",
  'Yaxis'                         , 'y', 1, "character",
  'output_path'                   , 'o', 2, "character",
  'output_dir'                    , 'd', 2, "character",
  #'input_path'                    , 'i', 2, "character",
  'config'                        , 'c', 2, "character",
  'force_config'                  , 'C', 0, "character",
  'motif_list'                    , 'l', 2, "character",
  'split_corr_by_starting_base'   , 's', 0, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

Xname <- file_path_sans_ext(basename(opt$Xaxis))
Yname <- file_path_sans_ext(basename(opt$Yaxis))

#if (is.null(opt$input_path)){
#  if (!(is.null(opt$force_config))){
#      print("still making this")
#      stop()
#  } else {
#    exts <- unique(file_ext(c(opt$Xaxis, opt$Yaxis)))
#    dirs <- unique(dirname(file.path(getwd(),c(opt$Xaxis, opt$Yaxis))))
#    
#    if (length(exts) > 1 | length(dirs) > 1){
#      print("error: too many extensions and/or input dirs, cannot infer input path for .motif.R files")
#      stop() 
#    } else {
#      if (exts == "bam"){
#        #opt$input_path=file.path(dirs,"STATS/MOTIF_COUNTS")
#        opt$Xaxis=file.path(dirs,"STATS/MOTIF_COUNTS",paste(Xname, ".motif.R",sep=""))
#        opt$Yaxis=file.path(dirs,"STATS/MOTIF_COUNTS",paste(Yname, ".motif.R",sep=""))
#      } else if (exts == "stats"){
#        #opt$input_path=file.path(dirs,"MOTIF_COUNTS")
#        opt$Xaxis=file.path(dirs,"MOTIF_COUNTS",paste(Xname, ".motif.R",sep=""))
#        opt$Yaxis=file.path(dirs,"MOTIF_COUNTS",paste(Yname, ".motif.R",sep=""))
#      } 
#    }
#  }
#} else {
#  opt$Xaxis=file.path(opt$input_path,paste(Xname, ".motif.R",sep=""))
#  opt$Yaxis=file.path(opt$input_path,paste(Yname, ".motif.R",sep=""))
#}


if (is.null(opt$split_corr_by_starting_base)){separate=""} else {separate="_separate"}
#### ORIGNAL PATH STUFF
#if (is.null(opt$output_path)){
#  opt$output_path <- file.path(getwd(), paste(Xname,"_",Yname,separate,"_corr_motif.pdf", sep=""))
#  } else {
#  if (file_ext(opt$output_path) == ""){
#    opt$output_path <- file.path(opt$output_path, paste(Xname,"_",Yname,separate,"_corr_motif.pdf", sep=""))
#  }
#}


if (is.null(opt$output_path)){
  if (is.null(opt$output_dir)){
    pat_out <- getwd()
  } else {
    pat_out <- opt$output_dir
  }
  opt$output_path <- file.path(pat_out, paste(Xname,"_",Yname,separate,"_corr_motif.pdf", sep=""))
} else {
  if (file_ext(opt$output_path) == ""){
    opt$output_path <- file.path(opt$output_path, paste(Xname,"_",Yname,separate,"_corr_motif.pdf", sep=""))
  }
}



  df1 <- data.frame(readRDS(opt$Xaxis))
  df2 <- data.frame(readRDS(opt$Yaxis))
  
  df1c <- sum(df1$Freq)
  df2c <- sum(df2$Freq)
  
  #Mergin#
  mrgd <- merge (df1, df2, by = "Var1")
  init <- c(substr(mrgd$Var1, 1,1))
  rownames(mrgd) <- mrgd$Var1
  
  mrgd$Freq.x <- mrgd$Freq.x/sum(mrgd$Freq.x)
  mrgd$Freq.y <- mrgd$Freq.y/sum(mrgd$Freq.y)
  
  max <- max(mrgd$Freq.x,mrgd$Freq.y, na.rm = TRUE)
  
  reg <- lm(Freq.y ~ Freq.x, data=mrgd)
  regsum <- summary(reg)
  resd <- regsum$residuals
  #negresd <- resd[which(resd < 0)]
  #posresd <- resd[which(resd > 0)]
  r2<- round(regsum$adj.r.squared,digits=3)

  
  if (all(names(resd) == rownames(mrgd))){
      mrgd <- cbind(mrgd, init, resd)
      
      #### PLOT
      
      
      
      if (is.null(opt$split_corr_by_starting_base)){
        
        
        plot <- ggplot(data = mrgd, aes(x= Freq.x, y = Freq.y))+ 
                       xlim(0, max)+
                       ylim(0, max)+
                       geom_point(size = 0.6)+
                       geom_abline(linetype = "dotdash", colour = "red3", size = 0.3  )+
                       geom_smooth(method = lm, se = FALSE, data=mrgd, size = 0.4)+
                       ggtitle(paste("R2=", r2, sep="" ))+
                       labs(x= paste(basename(Xname), df1c), y= paste(basename(Yname), df2c))
      } else {
      
        plot <- ggplot(data = mrgd, 
                       mapping = aes(x= Freq.x, y = Freq.y, colour = init ))+ 
                       xlim(0, max)+
                       ylim(0, max)+
                       geom_point(size = 0.6)+
                       geom_abline(linetype = "dotdash", colour = "red3", size = 0.3  )+
                       geom_smooth(method = lm, se = FALSE, data=mrgd, size = 0.4)+
                       #geom_smooth(method = lm, se = FALSE, data = subset(mrgd,init == "A"), size = 0.4)+
                       #geom_smooth(method = lm, se = FALSE, data = subset(mrgd,init == "C"), size = 0.4)+
                       #geom_smooth(method = lm, se = FALSE, data = subset(mrgd,init == "G"), size = 0.4)+
                       #geom_smooth(method = lm, se = FALSE, data = subset(mrgd,init == "T"), size = 0.4)+
                       ggtitle(paste("R2=", r2, sep="" ))+
                       labs(x= paste(basename(Xname), df1c), y= paste(basename(Yname), df2c))
        plot <- plot + scale_color_brewer(palette = "Set1")
        

      }
    
      pdf(opt$output_path,5.7,5)
      print(plot) 
      dev.off()
  
  } else {
      
    print("error in motif order")
    
}


  
  