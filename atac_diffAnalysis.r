#!/bin/r

library(DiffBind)
bamdir="/Users/ehresms/computational/ATAC/clean_bam"
beddir="/Users/ehresms/computational/ATAC/peaks"
files=list.files(bamdir)
dataset=data.frame(matrix(nrow=length(files), ncol = 7))
colnames(dataset)<-c("SampleID", "Condition", "Treatment", "Replicate", "bamReads", "Peaks", "PeakCaller")

n=1
for (file in files){
  
  name=strsplit(file, split = "_")[[1]][1]
  bamfile=paste0(bamdir, "/", file)
  
  bedfile=paste0(beddir, "/", name, "/", name , "_peaks.narrowPeak")
  
  if (grepl(name, pattern = "NT")){
    treatment="NT"
  } else {
    treatment="TGFB"
  }
  
  replicate=reverse(strsplit(reverse(name), split = "-")[[1]][1])
  
  if( grepl(name, pattern = "D0") ){
    condition="D0"
  } else {
    condition=strsplit(name, split = "-")[[1]][2]
  }
  
  dataset[n,]<-c(name, condition, treatment, replicate, bamfile, bedfile, "narrowPeak")
  n=n+1
}

test<-dba(sampleSheet = dataset, DBA = test)

test<-dba.count(test)
plot(test)
test<-dba.normalize(test, method = "DBA_DESEQ2")

test<-dba.contrast(test, design="~Condition + Treatment")


## Start from dba object made on compute canada
test<-readRDS("/Users/ehresms/computational/ATAC/testdbaobject.rds")

n=1
for (file in files){ ##change location of bam and bed files
  
  name=strsplit(file, split = "_")[[1]][1]
  bamfile=paste0(bamdir, "/", file)
  bedfile=paste0(beddir, "/", name, "/", name , "_peaks.narrowPeak")
  
  test$samples$bamReads[n]<-bamfile
  test$samples$Peaks[n]<-bedfile
  test$class[grep(test$class, pattern = paste0(name, "_"))]<-bamfile
  
  
  n=n+1
}  
  
## Get consensus peaks with DiffBind

test<-dba.analyze(test, design= "~Treatment + Condition")
test<-dba.contrast(test, design= "~Condition + Treatment",  contrast = c("Treatment", "TGFB", "NT"))
dba.plotVolcano(test, contrast = 11) ##contrast tgfb v nt

repObj<-dba.report(test, contrast=11, bDB = TRUE)

dba.plotProfile(test)
profiles<-dba.plotProfile(test, sites = repObj)
dba.plotProfile(profiles, labels = c("DBA_CONDITION", "DBA_TREATMENT"))


## Get Diffenretial opening with DESeq2 
library(DESeq2)
library(DEGreport)

counts_raw<-dba.peakset(test, bRetrieve = TRUE, DataType = DBA_DATA_FRAME ) ## get raw count values
countmat<-as.matrix(counts_raw[,c(4:ncol(counts_raw))])
countmat<-countmat[,c(5,6, 1:4, 7:22)]
countmat_timeseries<-matrix(nrow=nrow(countmat), ncol=ncol(countmat)+2)

countmat_timeseries[,3:24]<-countmat
countmat_timeseries[,1:2]<-countmat[,1:2]
countmat_timeseries[,5]<-countmat_timeseries[,6]

countdata<-countmat_timeseries
peaknames<-paste0(counts_raw$CHR, "_", counts_raw$START, "_", counts_raw$END)
row.names(countdata)<-peaknames
colnames(countdata)<-c("MCF10A.D0.TGFB.1", "MCF10A.D0.TGFB.2", colnames(counts_raw)[8:9], colnames(counts_raw)[4:7], colnames(counts_raw)[10:25])

coldata<-as.data.frame(matrix(nrow=ncol(countdata), ncol = 2))
row.names(coldata)<-colnames(countdata)
colnames(coldata)<-c("Condition", "Time")
coldata$Condition<-c(rep(c("NT", "NT", "TGFB", "TGFB"), 6))
coldata$Time<-c(rep("D0", 4), rep("D0.5", 4), rep("D1", 4), rep("D2", 4), rep("D4", 4), rep("D6", 4))

dds<-DESeqDataSetFromMatrix(countData = round(countdata) , colData = coldata, design = ~ Time + Condition + Time:Condition)

dds <- DESeq(dds, test="LRT", reduced = ~ Time + Condition) ## Get significant results across timecourse
res<-results(dds)
res$padj[is.na(res$padj)]<-1
sig_res<-res[res$padj<0.05,]


##Get clusters with DEGreport
norm_counts<-log(countdata[row.names(sig_res),], 2)

clusters<-degPatterns(norm_counts, metadata = coldata, time = "Time", col = "Condition") ##get clusters

cluster_groups<-clusters$df

##Make ggplots for clusters
library(ggplot2)
ggplot(res_cluster_15, aes(Time, value)) + geom_boxplot(aes(fill = Condition), notch = TRUE, lwd = 1) + 
geom_point(position = position_jitterdodge(dodge.width = 0.8), aes(fill = Condition), alpha = 2/10) + scale_fill_manual( values = c("White", "Red"))


res_data<-clusters$raw

res_data$Time[res_data$Time == "D0"]<-1
res_data$Time[res_data$Time == "D0.5"]<-2
res_data$Time[res_data$Time == "D1"]<-3
res_data$Time[res_data$Time == "D2"]<-4
res_data$Time[res_data$Time == "D4"]<-5
res_data$Time[res_data$Time == "D6"]<-6

res_data$Time<-as.numeric(res_data$Time)

##Make graphs for all data

ggplot(data = res_data, aes(x = Time, y = value, group = genes, color = Condition)) +
   facet_wrap(~ cluster) + geom_smooth(aes(group=Condition), method = "loess") +
  scale_color_manual(values = c("blue", "red")) +theme_classic() +
  scale_x_continuous(breaks=c(1:6), labels=c("D0", "D0.5", "D1", "D2", "D4", "D6")) + ylab(label = "Peak count")



ggplot(data = res_data, aes(x = Time, y = value, group = genes, color = Condition)) +
  facet_wrap(~ cluster) + geom_path() +
  scale_color_manual(values = c("blue", "red")) +theme_classic() +
  scale_x_continuous(breaks=c(1:6), labels=c("D0", "D0.5", "D1", "D2", "D4", "D6")) + ylab(label = "Peak count")


ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/test_grif.tiff", device = "tiff", width = 20, height = 10, units = "in")


##Make graphs for selected clusters
number_inclusters<-data.frame(matrix(nrow = length(unique(res_data$cluster)), ncol = 1))
for (cluster in unique(res_data$cluster)){
  number=nrow(res_data[res_data$cluster == cluster,])
  number_inclusters[cluster, ]<-number
}

row.names(number_inclusters)<-paste0("Cluster_", row.names(number_inclusters))
colnames(number_inclusters)<-"number"
number_inclusters$number[is.na(number_inclusters)]<-1

big_clusters<-number_inclusters[number_inclusters > 500*6]
big_clusters<-row.names(number_inclusters)[number_inclusters$number %in% big_clusters]
big_cluster_names<-as.numeric(sapply(strsplit(big_clusters, split = "_"), "[", 2))

res_selected<-res_data[res_data$cluster %in% big_cluster_names,]


ggplot(data = res_selected, aes(x = Time, y = value, group = genes, color = Condition)) +
  facet_wrap(~ cluster) + geom_smooth(aes(group=Condition), method = "loess") +
  scale_color_manual(values = c("black", "red")) + theme_classic() +
  scale_x_continuous(breaks=c(1:6), labels=c("D0", "D0.5", "D1", "D2", "D4", "D6")) + ylab(label = "Relative Peak Density")



## Get peak locations per cluster


peak_locations<-data.frame(matrix(nrow = nrow(norm_counts), ncol = 2))
peak_locations[,1]<-sapply(strsplit(row.names(countdata[row.names(sig_res),]), split = "_"), "[", 1)
peak_locations[,2]<-sapply(strsplit(row.names(countdata[row.names(sig_res),]), split = "_"), "[", 2)
peak_locations[,3]<-sapply(strsplit(row.names(countdata[row.names(sig_res),]), split = "_"), "[", 3)
peak_locations[,4]<-paste0("peak_", 1:nrow(peak_locations))
peak_locations[,5]<-"."
peak_locations[,6]<-"."
peak_locations[,7]<-row.names(countdata[row.names(sig_res),])


selected_starts<-unique(res_selected$genes[res_selected$cluster == 15]) ## change this to get new set of peaks

selected_locations<-peak_locations[peak_locations$V7 %in% selected_starts,]

write.table(selected_locations[,1:6], col.names = FALSE, row.names = FALSE,
            file = "/Users/ehresms/computational/ATAC/diffbind/cluster_15.bed", sep = "\t",
            quote = FALSE)
