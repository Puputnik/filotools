#!/usr/bin/Rscript
suppressMessages(library(ShortRead))
#overwrite = F

pat <- commandArgs(trailingOnly=TRUE)[1]
patsave <- commandArgs(trailingOnly=TRUE)[2]
f <- commandArgs(trailingOnly=TRUE)[3]

#pat <- "/tank/LIQUID_BIOPSY_GRANT/FASTQ/GUPPY_DEM_EITHER"
#pat <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L06_FAS60961_count0_partial/no_sample/20220907_1751_MN33722_FAS60961_4d879b95/fastq_pass"
#patsave <- file.path(pat, "READL_RAW")
dir.create(patsave, showWarnings = FALSE)
#tot <- data.frame()
#fili <- list.files(pat, ".*FAT85247.*")
#fili <- list.files(pat, "*fastq")
#for (f in fili){
  #f=fili[14]
  print(f)
  outf <- file.path(patsave,paste(tools::file_path_sans_ext(tools::file_path_sans_ext(basename(f))), "_RAW.readlength_counts.R", sep= ""))
  print(outf)
  #if (file.exists(outf) & !(overwrite)){
  #  tab <- readRDS(outf)
  #} else {
    Var1 <- width(readFastq(file.path(pat, f)))
    tab <- table(Var1)
    saveRDS(tab, outf)
  #}
  
#  tab <- cbind(data.frame(tab),f)
#  tab$Freq <- tab$Freq/sum(tab$Freq)
#  tot <- rbind(tot , tab)
#}

#tot$Var1 <- as.numeric(as.character(tot$Var1))
#totb <- tot
#tot <- subset(totb, grepl("B03", totb$f))
#
#p<- ggplot(tot,aes(x=Var1,weight=Freq,col=f)) + 
#  geom_density(size=0.3, bw=2) + xlim(0,425) + xlab("fragment length")  + 
#  theme(text = element_text(size = 20)) +
#  guides(color = guide_legend(override.aes = list(size = 3) ) )+
#  theme(legend.position = "bottom",legend.direction = "vertical")
#print(p)
