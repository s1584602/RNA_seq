#install and load packages ----
install.packages("ggpubr")

library(tidyverse)
library("ggpubr")

#install functions ----
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


#load raw data ----


# load mean read depth data, 10 nt gap around pre-miRNA, 95 nt long flanks
##move column names one column to right
##replace hsa-mir with hsa-miR
##remove " 100margin"
##replace "=-ifn" with NA
##meanraw <- read.table(file = "clipboard", sep = "\t", header=TRUE)
##write.csv(meanraw[,c(1,8:56)],"meanraw.csv")
meanraw <- read.csv("meanraw.csv", header=TRUE)[,2:51]

# load miRBase data about chrompsomal location and direction miRNAs
miRBase_overview<- read.table(file = "miRBase_overview.txt", 
                              sep = "\t", header=TRUE)
miRBase_overview$pri<-miRBase_overview$ID%>%gsub("hsa-mir-", "hsa-miR-",.)
miRBase_overview<-miRBase_overview[,c(9,4,7:8)]

#load mature miRNA data
mature_KO<-read.csv("PA1__tRNA_WT1__results_KO.csv")
mature_H<-read.csv("PA1__tRNA_WT1__results_H.csv")
rownames(mature_KO)<-mature_KO[,1]
rownames(mature_H)<-mature_H[,1]
colnames(mature_KO)<-c("miRNA","baseMean","mature_FC_KO","lfcSE", "pvalue","mature_KO_padj")
colnames(mature_H)<-c("miRNA","baseMean","mature_FC_H","lfcSE", "pvalue","mature_H_padj")
mature_KO$pri<-mature_KO$miRNA%>%gsub("-3p", "",.)%>%gsub("-5p", "",.)
mature_H$pri<-mature_H$miRNA%>%gsub("-3p", "",.)%>%gsub("-5p", "",.)



#load data conrad et al. (2014)
##conradraw <- read.table(file = "clipboard", sep = "\t", header=TRUE)
##conradraw$pri<-conradraw$miRNA%>%gsub("hsa-mir-", "hsa-miR-",.)
##write.csv(conradraw,"conradraw.csv")
conradraw<-read.csv("conradraw.csv", header=TRUE)[,2:40]
rownames(conradraw)<-conradraw$pri

#load data Rice at al. (2020)
##riceraw <- read.table(file = "clipboard", sep = "\t", header=TRUE)
##riceraw$pri<-riceraw$miRNA%>%gsub("hsa-mir-", "hsa-miR-",.)
##write.csv(riceraw,"riceraw.csv")
riceraw<-read.csv("riceraw.csv", header=TRUE)[,2:20]
rownames(riceraw)<-riceraw$pri


#load TAM2.0 database

TAM2.0_angiogenesis<-strsplit("hsa-let-7b	hsa-mir-19b-1	hsa-mir-34a	hsa-mir-486-2	hsa-mir-133a-1
                              hsa-mir-891a	hsa-mir-221	hsa-mir-1246	hsa-mir-320a	hsa-mir-20a
                              hsa-mir-133b	hsa-mir-126	hsa-let-7f-1	hsa-let-7f-2	hsa-mir-17
                              hsa-mir-190b	hsa-mir-149	hsa-mir-382	hsa-mir-92a-1	hsa-mir-486-1
                              hsa-mir-129-1	hsa-mir-134	hsa-mir-22	hsa-mir-16-1	hsa-mir-19b-2
                              hsa-mir-190a	hsa-mir-31	hsa-mir-30a	hsa-mir-378e	hsa-mir-378g
                              hsa-mir-718	hsa-mir-16-2	hsa-mir-640	hsa-mir-1275	hsa-mir-146a
                              hsa-mir-222	hsa-mir-130a	hsa-mir-184	hsa-mir-145	hsa-mir-370
                              hsa-mir-378a	hsa-mir-10b	hsa-mir-363	hsa-mir-18a	hsa-mir-21
                              hsa-mir-378b	hsa-mir-122	hsa-mir-320e	hsa-mir-133a-2	hsa-mir-492
                              hsa-mir-320b-2	hsa-mir-210	hsa-mir-15a	hsa-mir-205	hsa-mir-320b-1
                              hsa-mir-378f	hsa-mir-93	hsa-mir-19a	hsa-mir-27b	hsa-mir-497
                              hsa-mir-92a-2	hsa-mir-296	hsa-mir-150	hsa-mir-378c	hsa-mir-143","[[:space:]]")[[1]]%>%
  gsub("mir", "miR",.)


# load genomic origin miRNAs (manually from UCSC GRCh38/hg38,Dec. 2013)
##origin<- read.table(file = "clipboard", sep = "\t", header=TRUE)
##write.csv(origin,"origin.csv")
origin<- read.csv(file = "origin.csv", header=TRUE)[,2:5]

# load mirtron data (Wen et al., 2015)
##mirtron<- read.table(file = "clipboard", sep = "\t", header=TRUE)
###write.csv(mirtron,"mirtron.csv")
mirtron<- read.csv(file = "mirtron.csv", header=TRUE)[,2:3]

#MPI_overview: orientation of primary RNA; MPI reorder from genome direction to transcript direction ----

rownames(meanraw)<-meanraw$miRNA
meanraw$pri<-meanraw$miRNA
MPI_raw_miRBase<-left_join(meanraw,miRBase_overview)

MPI_overview_raw<-data.frame(MPI_raw_miRBase$pri)

MPI_overview_raw$S13_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S13_LEFT_MEAN,MPI_raw_miRBase$S13_RIGHT_MEAN)
colnames(MPI_overview_raw)<-c("pri","S13_upstream")
MPI_overview_raw$S13_pre<-MPI_raw_miRBase$S13_MIDDLE_MEAN
MPI_overview_raw$S13_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S13_LEFT_MEAN,MPI_raw_miRBase$S13_RIGHT_MEAN)

MPI_overview_raw$S14_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S14_LEFT_MEAN,MPI_raw_miRBase$S14_RIGHT_MEAN)
MPI_overview_raw$S14_pre<-MPI_raw_miRBase$S14_MIDDLE_MEAN
MPI_overview_raw$S14_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S14_LEFT_MEAN,MPI_raw_miRBase$S14_RIGHT_MEAN)

MPI_overview_raw$S15_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S15_LEFT_MEAN,MPI_raw_miRBase$S15_RIGHT_MEAN)
MPI_overview_raw$S15_pre<-MPI_raw_miRBase$S15_MIDDLE_MEAN
MPI_overview_raw$S15_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S15_LEFT_MEAN,MPI_raw_miRBase$S15_RIGHT_MEAN)

MPI_overview_raw$S16_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S16_LEFT_MEAN,MPI_raw_miRBase$S16_RIGHT_MEAN)
MPI_overview_raw$S16_pre<-MPI_raw_miRBase$S16_MIDDLE_MEAN
MPI_overview_raw$S16_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S16_LEFT_MEAN,MPI_raw_miRBase$S16_RIGHT_MEAN)

MPI_overview_raw$S17_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S17_LEFT_MEAN,MPI_raw_miRBase$S17_RIGHT_MEAN)
MPI_overview_raw$S17_pre<-MPI_raw_miRBase$S17_MIDDLE_MEAN
MPI_overview_raw$S17_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S17_LEFT_MEAN,MPI_raw_miRBase$S17_RIGHT_MEAN)

MPI_overview_raw$S18_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S18_LEFT_MEAN,MPI_raw_miRBase$S18_RIGHT_MEAN)
MPI_overview_raw$S18_pre<-MPI_raw_miRBase$S18_MIDDLE_MEAN
MPI_overview_raw$S18_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S18_LEFT_MEAN,MPI_raw_miRBase$S18_RIGHT_MEAN)

MPI_overview_raw$S19_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S19_LEFT_MEAN,MPI_raw_miRBase$S19_RIGHT_MEAN)
MPI_overview_raw$S19_pre<-MPI_raw_miRBase$S19_MIDDLE_MEAN
MPI_overview_raw$S19_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S19_LEFT_MEAN,MPI_raw_miRBase$S19_RIGHT_MEAN)

MPI_overview_raw$S20_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S20_LEFT_MEAN,MPI_raw_miRBase$S20_RIGHT_MEAN)
MPI_overview_raw$S20_pre<-MPI_raw_miRBase$S20_MIDDLE_MEAN
MPI_overview_raw$S20_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S20_LEFT_MEAN,MPI_raw_miRBase$S20_RIGHT_MEAN)

MPI_overview_raw$S21_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S21_LEFT_MEAN,MPI_raw_miRBase$S21_RIGHT_MEAN)
MPI_overview_raw$S21_pre<-MPI_raw_miRBase$S21_MIDDLE_MEAN
MPI_overview_raw$S21_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S21_LEFT_MEAN,MPI_raw_miRBase$S21_RIGHT_MEAN)

MPI_overview_raw$S22_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S22_LEFT_MEAN,MPI_raw_miRBase$S22_RIGHT_MEAN)
MPI_overview_raw$S22_pre<-MPI_raw_miRBase$S22_MIDDLE_MEAN
MPI_overview_raw$S22_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S22_LEFT_MEAN,MPI_raw_miRBase$S22_RIGHT_MEAN)

MPI_overview_raw$S23_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S23_LEFT_MEAN,MPI_raw_miRBase$S23_RIGHT_MEAN)
MPI_overview_raw$S23_pre<-MPI_raw_miRBase$S23_MIDDLE_MEAN
MPI_overview_raw$S23_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S23_LEFT_MEAN,MPI_raw_miRBase$S23_RIGHT_MEAN)

MPI_overview_raw$S24_upstream<-ifelse(MPI_raw_miRBase$Strand=="+",MPI_raw_miRBase$S24_LEFT_MEAN,MPI_raw_miRBase$S24_RIGHT_MEAN)
MPI_overview_raw$S24_pre<-MPI_raw_miRBase$S24_MIDDLE_MEAN
MPI_overview_raw$S24_downstream<-ifelse(MPI_raw_miRBase$Strand=="-",MPI_raw_miRBase$S24_LEFT_MEAN,MPI_raw_miRBase$S24_RIGHT_MEAN)

MPI_overview_raw$total<-rowSums(MPI_overview_raw[,2:37])
MPI_overview_raw<-subset(MPI_overview_raw,total>0)

#Normalise to both flanks ----
flanksums<-colSums(MPI_overview_raw[,c(2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37)])

flanknormfactors<-c(rep(mean(flanksums[1:2]),3),
                    rep(mean(flanksums[3:4]),3),
                    rep(mean(flanksums[5:6]),3),
                    rep(mean(flanksums[7:8]),3),
                    rep(mean(flanksums[9:10]),3),
                    rep(mean(flanksums[11:12]),3),
                    rep(mean(flanksums[13:14]),3),
                    rep(mean(flanksums[15:16]),3),
                    rep(mean(flanksums[17:18]),3),
                    rep(mean(flanksums[19:20]),3),
                    rep(mean(flanksums[21:22]),3),
                    rep(mean(flanksums[23:24]),3))
MPI_overview<-as.data.frame(mapply('/', MPI_overview_raw[2:37], flanknormfactors)*mean(flanksums[1:8]))
MPI_overview$pri<-MPI_overview_raw$pri
MPI_overview$flank_norm<-rowSums(MPI_overview[,c(1,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28,30,31,33,34,36)])/24
MPI_overview$upstream<-rowSums(MPI_overview[,c(1,4,7,10,13,16,19,22,25,28,31,34)])/12
MPI_overview$downstream<-rowSums(MPI_overview[,c(3,6,9,12,15,18,21,24,27,30,33,36)])/12

MPI_overview$WT_flank<-rowSums(MPI_overview[,c(1,3,4,6,7,9,10,12)])/8
MPI_overview$KO_flank<-rowSums(MPI_overview[,c(13,15,16,18,19,21,22,24)])/8
MPI_overview$H_flank<-rowSums(MPI_overview[,c(25,27,28,30,31,33,34,36)])/8 
MPI_overview$flank_FC_KO<-MPI_overview$KO_flank/MPI_overview$WT_flank
MPI_overview$flank_FC_H<-MPI_overview$H_flank/MPI_overview$WT_flank    
rownames(MPI_overview)<-MPI_overview$pri

#calculate the MPI for each sample ----


#calculate MPI
MPI_overview$S13_MPI<-log((MPI_overview$S13_pre+1)/(((MPI_overview$S13_upstream+1)+(MPI_overview$S13_downstream+1))/2),2)
MPI_overview$S14_MPI<-log((MPI_overview$S14_pre+1)/(((MPI_overview$S14_upstream+1)+(MPI_overview$S14_downstream+1))/2),2)
MPI_overview$S15_MPI<-log((MPI_overview$S15_pre+1)/(((MPI_overview$S15_upstream+1)+(MPI_overview$S15_downstream+1))/2),2)
MPI_overview$S16_MPI<-log((MPI_overview$S16_pre+1)/(((MPI_overview$S16_upstream+1)+(MPI_overview$S16_downstream+1))/2),2)
MPI_overview$S17_MPI<-log((MPI_overview$S17_pre+1)/(((MPI_overview$S17_upstream+1)+(MPI_overview$S17_downstream+1))/2),2)
MPI_overview$S18_MPI<-log((MPI_overview$S18_pre+1)/(((MPI_overview$S18_upstream+1)+(MPI_overview$S18_downstream+1))/2),2)
MPI_overview$S19_MPI<-log((MPI_overview$S19_pre+1)/(((MPI_overview$S19_upstream+1)+(MPI_overview$S19_downstream+1))/2),2)
MPI_overview$S20_MPI<-log((MPI_overview$S20_pre+1)/(((MPI_overview$S20_upstream+1)+(MPI_overview$S20_downstream+1))/2),2)
MPI_overview$S21_MPI<-log((MPI_overview$S21_pre+1)/(((MPI_overview$S21_upstream+1)+(MPI_overview$S21_downstream+1))/2),2)
MPI_overview$S22_MPI<-log((MPI_overview$S22_pre+1)/(((MPI_overview$S22_upstream+1)+(MPI_overview$S22_downstream+1))/2),2)
MPI_overview$S23_MPI<-log((MPI_overview$S23_pre+1)/(((MPI_overview$S23_upstream+1)+(MPI_overview$S23_downstream+1))/2),2)
MPI_overview$S24_MPI<-log((MPI_overview$S24_pre+1)/(((MPI_overview$S24_upstream+1)+(MPI_overview$S24_downstream+1))/2),2)

#mean MPI
MPI_overview$WT_MPI<-rowMeans(MPI_overview[,46:49])
MPI_overview$KO_MPI<-rowMeans(MPI_overview[,50:53])
MPI_overview$H_MPI<-rowMeans(MPI_overview[,54:57])


#fold change MPI
MPI_overview$MPI_FC_KO<-MPI_overview$KO_MPI-MPI_overview$WT_MPI
MPI_overview$MPI_FC_H<-MPI_overview$H_MPI-MPI_overview$WT_MPI


# allow mature miRNA join of pri duplicates
MPI_overview$duplicates<-ifelse(grepl("^.+(-1|-2|-3|-4|-5|-9)$",MPI_overview$pri)==TRUE,
                                MPI_overview$pri,
                                NA)%>%
  gsub('.{2}$', '',.)

MPI_overview[,c(37,63)]

##save complete table
write.csv(MPI_overview,"MPI_overview.csv")


#merge summarised data MPI with mature miRNA data ----

mature_overview<-left_join(mature_KO[,c(1:3,6)],mature_H[,c(1:3,6:7)])
mature_overview$duplicates<-mature_overview$pri

MPI_merge_mature_full<-full_join(subset(MPI_overview[,c(37:45,58:63)],flank_norm>2&upstream/downstream<4&downstream/upstream<4),
                                 left_join(mature_KO[,c(1:3,6)],mature_H[,c(1:3,6:7)]))



MPI_merge_no_duplicates<-left_join(subset(MPI_overview[,c(37:45,58:63)],
                                          flank_norm>2&upstream/downstream<4&downstream/upstream<4&
                                            (is.na(MPI_overview[,63]))),
                                   mature_overview[,1:7],
                                   by="pri")
MPI_merge_no_duplicates<-MPI_merge_no_duplicates[!is.na(MPI_merge_no_duplicates$miRNA),]


MPI_merge_duplicates<-na.omit(
  bind_rows(
    left_join(subset(MPI_overview[,c(37:45,58:63)],
                     flank_norm>2&upstream/downstream<4&downstream/upstream<4&
                       (!is.na(MPI_overview[,63]))),
              mature_overview[,c(1:6,8)],
              by="duplicates"),
    left_join(subset(MPI_overview[,c(37:45,58:63)],
                     flank_norm>2&upstream/downstream<4&downstream/upstream<4&
                       (!is.na(MPI_overview[,63]))),
              mature_overview[,1:7],
              by="pri")))

###MPI_merge_mature<-bind_rows(MPI_merge_no_duplicates,MPI_merge_duplicates)

MPI_merge_mature_full<-full_join(bind_rows(MPI_merge_no_duplicates,MPI_merge_duplicates),
                                 subset(mature_overview, miRNA %in% 
                                          setdiff(mature_overview$miRNA,
                                                  c(MPI_merge_duplicates$miRNA,MPI_merge_no_duplicates$miRNA))))[,c(1:2,5:21)]




##save complete tables
write.csv(mature_overview,"mature_overview.csv")
write.csv(MPI_merge_mature_full,"MPI_merge_mature_full.csv")




#merge summarised data with external sources ----
MPI_merge<-left_join(left_join(MPI_merge_mature_full,
                               full_join(miRBase_overview,
                                         full_join(conradraw[,c(6,39)],
                                                   riceraw[,c(10:12,16:19)]))),
                     origin)

MPI_merge$mirtron_origin<-ifelse(MPI_merge$origin2=="3p tailed mirtron","3p tailed mirtron",
                                 ifelse(MPI_merge$origin2=="5p tailed mirtron","5p tailed mirtron",
                                        ifelse(MPI_merge$origin2=="mirtron","mirtron",
                                               NA)))
MPI_merge$mirtron_origin<-replace_na(MPI_merge$mirtron_origin, "no_mirtron")


write.csv(MPI_merge,"MPI_merge.csv")


pri_distinct<-distinct(MPI_merge,pri,.keep_all = TRUE)%>% drop_na(flank_norm)
write.csv(pri_distinct,"pri_distinct.csv")


pri_distinct_nm<-pri_distinct%>%subset(mirtron_origin=="no_mirtron")

pri_distinct_nm$subset1<-ifelse(pri_distinct_nm$origin=="3UTR"|
                                  pri_distinct_nm$origin=="5UTR"|
                                  pri_distinct_nm$origin=="exon",
                                "exon",
                                ifelse(pri_distinct_nm$origin=="intergenic"|
                                         pri_distinct_nm$origin=="own primary transcript"|
                                         pri_distinct_nm$origin=="primary cluster",
                                       "intergenic/cluster",pri_distinct_nm$origin))




write.csv(pri_distinct_nm,"pri_distinct_nm.csv")

mature_distinct<-distinct(MPI_merge,miRNA,.keep_all = TRUE)
write.csv(mature_distinct,"mature_distinct.csv")