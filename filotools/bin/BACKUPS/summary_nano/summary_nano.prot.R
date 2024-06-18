

patte <- c("/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L01_FAS34889_count0/L01_FAS34889_count0/20220718_1639_MN33722_FAS34889_35f5914c/pore_activity_FAS34889_427854b8.csv",
"/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L01_FAS34889_count1/L01_FAS34889_count1/20220719_1142_MN33722_FAS34889_7f769c2c/pore_activity_FAS34889_8b36481d.csv",
"/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L02_FAT85244_count0/L02_FAT85244_count0/20220728_1653_MN33722_FAT85244_560cd7c7/pore_activity_FAT85244_24ec271c.csv",
"/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L02_FAT85244_count1/L02_FAT85244_count1/20220728_1857_MN33722_FAT85244_7c13fdcb/pore_activity_FAT85244_0752585b.csv")

for (pored in patte){
tar <- strsplit(pored, "/")[[1]][6]
tar <- "CAREGGI_02R9"
pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK112/CAREGGI_02R9/careggi_with9_4/no_sample/20220413_1147_MN33722_FAQ75633_5affa651/duty_time_FAQ75633_1c8082d3.csv"
tar <- "CAREGGI_02R904X"
pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK112/CAREGGI_02R9/careggi_0_4Xbeads/0_4Xbeads_9_4/20220419_1223_MN33722_FAP06081_932beeae/duty_time_FAP06081_e96f0418.csv"
tar <- "CAREGGI_02R90_fewpores"
pored <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK112/CAREGGI_02R9/careggi9_4biswithfewpores/careggi94withfewpores/20220415_1000_MN33722_FAQ75633_7fd15cb6/duty_time_FAQ75633_ed1fa246.csv"

pored <- "/tank/LIQUID_BIOPSY_MAIN/RUNS/Mol_Cancer/M2/20191023_1139_GA10000_FAL01577_6ba7621a/reports/duty_time.csv"

#pored <- read.table(pored, sep=",", header=T)
pored <- read.table(pored, sep=",", header=T, fill=T)

#unique(pore$Channel.State)
colnames(pored) 

pore <- subset(pored, pored$Channel.State %in% c("adapter", "pore", "strand"))

pdf(paste("/tank/USB3/LB_ANALYSIS_2022/",tar,".pdf", sep=""))
p <- ggplot(pore, aes(fill=Channel.State, y=State.Time..samples., x=Experiment.Time..minutes.)) + 
  geom_bar(position="fill", stat="identity")
print(p)
dev.off()

pore <- subset(pored, pored$Channel.State %in% c("adapter", "pore", "strand", "unblocking"))
pore <- pored
pdf(paste("/tank/USB3/LB_ANALYSIS_2022/",tar,".una.pdf", sep=""))


pored <- "/tank/LIQUID_BIOPSY_MAIN/RUNS/Mol_Cancer/M2/20191023_1139_GA10000_FAL01577_6ba7621a/reports/duty_time.csv"
pored <- read.table(pored, sep=",", header=T, fill=T)
pore <- pored
unique(pore$Channel.State)

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
                 "unclassified"                ,
                 "unclassified_following_reset",
                 "pending_manual_reset"        ,
                 "pending_mux_change"          ,
                 "locked"                      ))

target_labels <- rev(c("Strand",
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
                   "unclassified",
                   "unclassified_following_reset",
                   "pending_manual_reset",               
                   "pending_mux_change",                 
                   "Pending Reselection"))

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
              "#455456",   
              "#455456",	
              "#455456",   
              "#455456",   
              "#ccd5d8"))  

pore$Channel.State <- factor(pore$Channel.State, levels=orig_labels)

#scale_color_manual(labels = c("T999", "T888"), values = c("blue", "red")) +
  

p <- ggplot(pore, aes(fill=Channel.State, y=State.Time..samples., x=Experiment.Time..minutes.)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_discrete(name = "Channel State", labels = target_labels) +
  scale_x_binned(
    name = waiver(),
    n.breaks = 30,
    nice.breaks = TRUE)  +
  scale_fill_manual(labels= target_labels, values=lab_colors) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
print(p)
dev.off()

}





library(lubridate)
seconds_to_period(86400)
minutes_to_period(86400)














pare <- subset(pore, pore$Experiment.Time..minutes. == "30")
pare <- subset(pore, pore$Experiment.Time..minutes. == "500")
pare <- subset(pore, pore$Experiment.Time..minutes. == "1100")

pare <- cbind(pare, pare$State.Time..samples./sum(pare$State.Time..samples.))
colnames(pare)[4] <- "rapp"

plot(pare$rapp,xaxt='n', xlab="", ylim=c(0,0.35))
points(pare$rapp,xaxt='n', xlab="", col="red")
points(pare$rapp,xaxt='n', xlab="", col="blue")
axis(side=1,at=seq(1,17),labels=pare$Channel.State, las=2)


report <- "/tank/LIQUID_BIOPSY_GRANT/RUNS/LSK109/L01_FAS34889_count0/L01_FAS34889_count0/20220718_1639_MN33722_FAS34889_35f5914c/throughput_FAS34889_427854b8.csv"
report <- read.table(report, sep=",", header=T)
upper <- length(rownames(report))
lower <- upper - ((60*4)+15)

upper=1000
lower=940


tru <- report[upper,2] - report[lower,2]
rpm <- tru/(upper-lower)
rph <- rpm*60


c = length(rownames(report))
while (c >= 1){
  if (c > 1){
  report[c,"permin"] <- report[c,2]-report[c-1,2]
  } else {
  report[c, "permin"] <- report[c,2]
  }
  c=c-1
}
plot(report[,"permin"]*60)
