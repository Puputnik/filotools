#!/usr/bin/Rscript

#path1 <- "/tank/USB3/LB_ANALYSIS_2022/FIRST_1M_READS/STATS/MOTIF/MOTIF_COUNTS/BC08.HACF1M.motif.R"
#path2 <- "/tank/USB3/LB_ANALYSIS_2022/FIRST_1M_READS/STATS/MOTIF/MOTIF_COUNTS/BC09.HACF1M.motif.R"
#outname <- commandArgs(trailingOnly=TRUE)[3]

path1 <- commandArgs(trailingOnly=TRUE)[1]
path2 <- commandArgs(trailingOnly=TRUE)[2]
outname <- commandArgs(trailingOnly=TRUE)[3]

if (is.na(path1) ) {
  stop("INSERIRISCI IL PERCORSO DEI CAMPIONI", call.=FALSE)
} else if (is.na(outname)) {
  outname = "output.pdf"
}

pdf(outname)

  df1 <- data.frame(readRDS(path1))
  df2 <- data.frame(readRDS(path2))
  
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
  
    #GGPLOTTIN
    library(tidyverse)
    
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
      labs(x= basename(path1), y= basename(path2))
    
    plot <- plot + scale_color_brewer(palette = "Set1")
    print(plot) 
    
    dev.off()
  
  } else {
      
    print("error in motif order")
    
  }


  
  