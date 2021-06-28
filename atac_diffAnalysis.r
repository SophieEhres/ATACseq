#!/bin/r

library(DiffBind)
library(stringr)
library(ggplot2)
library(rtracklayer)
library(DESeq2)
library(DEGreport)
library(plyr)
library(AnnotationHub)
library(basejump)
library(scales)
library(RColorBrewer)
library(GenomicFeatures)
library(ensembldb)
library(ggforce)
library(scico)

`%notin%` <- Negate(`%in%`)

bamdir="/Users/ehresms/computational/ATAC/clean_bam"
beddir="/Users/ehresms/computational/ATAC/bed"
files=list.files(bamdir)
files<-grep(files, pattern = "bai", invert = TRUE, value = TRUE)
dataset=data.frame(matrix(nrow=length(files), ncol = 7))
colnames(dataset)<-c("SampleID", "Condition", "Treatment", "Replicate", "bamReads", "Peaks", "PeakCaller")

n=1
for (file in files){
  
  name=strsplit(file, split = "_")[[1]][1]
  bamfile=paste0(bamdir, "/", file)
  day=str_split(name, pattern = "-")[[1]][2:3]
  if (day[1] == "D0"){
    day="D0"
  } else {
    day<-paste0(day[1], "-", day [2])
  }
  
  bedfile=grep(grep(list.files(beddir, pattern = ".bed"), pattern = day, value = TRUE), pattern = "sorted", invert = TRUE, value = TRUE)
  bedfile=paste0(beddir, "/", bedfile)
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


test<-dba(sampleSheet = dataset)

test<-dba.count(test)
plot(test)
test<-dba.normalize(test, method = "DBA_ALL_METHODS")

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
    
## Get consensus peaks with DiffBind [NOT USED]
  test<-dba.analyze(test, design= "~Treatment + Condition")
  test<-dba.contrast(test, design= "~Condition + Treatment",  contrast = c("Treatment", "TGFB", "NT"))
  dba.plotVolcano(test, contrast = 11) ##contrast tgfb v nt
  
  repObj<-dba.report(test, contrast=11, bDB = TRUE)
  
  dba.plotProfile(test)
  profiles<-dba.plotProfile(test, sites = repObj)
  dba.plotProfile(profiles, labels = c("DBA_CONDITION", "DBA_TREATMENT"))

##Export consensus bed file
  peakset<-dba.peakset(test, bRetrieve = TRUE)
  peakset<-as.data.frame(peakset)
  peakset$seqnames<-sapply(str_split(peakset$seqnames, pattern = "chr"), "[", 2 )
  export.bed(peakset, "/Users/ehresms/computational/ATAC/diffbind/consensus_peakset.bed")

## Get Differential opening with DESeq2 
  counts_raw<-dba.peakset(test, bRetrieve = TRUE, DataType = DBA_DATA_FRAME) ## get raw count values
  
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

##Get results with different FDR value
  res<-results(dds)
  res$padj[is.na(res$padj)]<-1
  sig_res<-res[res$padj<0.01,]

##Get clusters with DEGreport
  norm_counts<-log(countdata[row.names(sig_res),], 2)
  
  clusters<-degPatterns(norm_counts, metadata = coldata, time = "Time", col = "Condition") ##get clusters
  
  cluster_groups<-clusters$df

##Make dataframe for plots
  res_data<-clusters$raw
  
  res_data$Time[res_data$Time == "D0"]<-1
  res_data$Time[res_data$Time == "D0.5"]<-2
  res_data$Time[res_data$Time == "D1"]<-3
  res_data$Time[res_data$Time == "D2"]<-4
  res_data$Time[res_data$Time == "D4"]<-5
  res_data$Time[res_data$Time == "D6"]<-6
  
  res_data$Time<-as.numeric(res_data$Time)

##Make line graphs for all data

  ggplot(data = res_data, aes(x = Time, y = value, group = genes, color = Condition)) +
     facet_wrap(~ cluster) + geom_smooth(aes(group=Condition), method = "loess") +
    scale_color_manual(values = c("black", "red")) +theme_classic() +
    scale_x_continuous(breaks=c(1:6), labels=c("D0", "D0.5", "D1", "D2", "D4", "D6")) + ylab(label = "Peak count") +
    theme(text = element_text(family = "Arial", size = 20))
  
  ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/ATACclusters_allsig_smooth.tiff", device = "tiff")
  
  
  ggplot(data = res_data, aes(x = Time, y = value, group = genes, color = Condition)) +
    facet_wrap(~ cluster) + geom_path() +
    scale_color_manual(values = c("black", "red")) +theme_classic() +
    scale_x_continuous(breaks=c(1:6), labels=c("D0", "D0.5", "D1", "D2", "D4", "D6")) + ylab(label = "Peak count") +
    theme(text = element_text(family = "Arial", size = 20))
  
  
  ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/ATACclusters_allsig_lines.tiff", device = "tiff")


##Make line graphs for selected clusters
##Get "Big" Clusters
  number_inclusters<-data.frame(matrix(nrow = length(unique(res_data$cluster)), ncol = 1))
  for (cluster in unique(res_data$cluster)){
    number=nrow(res_data[res_data$cluster == cluster,])
    number_inclusters[cluster, ]<-number
  }
  
  row.names(number_inclusters)<-paste0("Cluster_", row.names(number_inclusters))
  colnames(number_inclusters)<-"number"
  number_inclusters$number[is.na(number_inclusters)]<-1
  
  big_clusters<-number_inclusters[number_inclusters > 500*6] ## "500" genes per cluster. Change this number to get different clusters
  big_clusters<-row.names(number_inclusters)[number_inclusters$number %in% big_clusters]
  big_cluster_names<-as.numeric(sapply(strsplit(big_clusters, split = "_"), "[", 2))
  
  res_selected<-res_data[res_data$cluster %in% big_cluster_names,]
  
  number_big_clusters<-as.data.frame(matrix(nrow = length(number_inclusters[big_clusters,]), ncol = 2))
  colnames(number_big_clusters)<-c("number", "cluster")
  
  number_big_clusters$number<-paste0("n=", number_inclusters[big_clusters,"number"])
  number_big_clusters$cluster<-big_cluster_names
  labels_bigclusters<-as.list(number_big_clusters$number)
  names(labels_bigclusters)<-number_big_clusters$cluster
  
  cluster_labeller <- function(variable,value){
    return(labels_bigclusters[value])
  }
  
  
  ggplot(data = res_selected, aes(x = Time, y = value, group = genes, color = Condition)) +
    geom_smooth(aes(group=Condition), method = "loess") +
    facet_wrap(~ cluster, labeller = as_labeller(labels_bigclusters)) +
    scale_color_manual(values = c("black", "red")) + theme_classic() +
    scale_x_continuous(breaks=c(1:6), labels=c("D0", "D0.5", "D1", "D2", "D4", "D6")) + ylab(label = "Relative Peak Density") +
    theme(text = element_text(family = "Arial", size = 40)) 
  
  
  ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/ATACclusters_bigclusters_smooth.tiff", device = "tiff")


## Get peak locations per cluster and export BED file
  
  peak_locations<-data.frame(matrix(nrow = nrow(norm_counts), ncol = 2))
  peak_locations[,1]<-sapply(strsplit(row.names(countdata[row.names(sig_res),]), split = "_"), "[", 1)
  peak_locations[,2]<-sapply(strsplit(row.names(countdata[row.names(sig_res),]), split = "_"), "[", 2)
  peak_locations[,3]<-sapply(strsplit(row.names(countdata[row.names(sig_res),]), split = "_"), "[", 3)
  peak_locations[,4]<-paste0("peak_", 1:nrow(peak_locations))
  peak_locations[,5]<-"."
  peak_locations[,6]<-"."
  peak_locations[,7]<-row.names(countdata[row.names(sig_res),])
  colnames(peak_locations)<-c("chr", "start", "end", "peak_ID", "score", "strand", "peakName")
  
  
  selected_starts<-unique(res_selected$genes[res_selected$cluster == 15]) ## change this to get new set of peaks
  
  selected_locations<-peak_locations[peak_locations$peakName %in% selected_starts,]
  selected_locations$chr<-sapply(str_split(selected_locations$chr, pattern = "chr"), "[", 2)
  
  
  write.table(selected_locations[,1:6], col.names = FALSE, row.names = FALSE,
              file = "/Users/ehresms/computational/ATAC/diffbind/cluster_15.bed", sep = "\t",
              quote = FALSE)

##Make GRanges for gene and transcript
  gene_gr<-makeGRangesFromEnsembl(organism = "Homo sapiens", level = "gene", genomeBuild = "hg38")
  transcript_gr<-makeGRangesFromGFF("/Users/ehresms/computational/genomes/human/hg38/Homo_sapiens.GRCh38.104.chr.gtf", level = "transcript")


##Make GRanges for peakset, significant peakset, and selected cluster
  peakset_gr<-as.data.frame(peakset)
  peakset_gr<-makeGRangesFromDataFrame(peakset_gr)
  
  sigPeakset_gr<-as.data.frame(peak_locations)
  sigPeakset_gr<-makeGRangesFromDataFrame(sigPeakset_gr)
  
  selected_gr<-selected_locations
  selected_gr<-makeGRangesFromDataFrame(selected_gr)


##Get nearest transcript and gene
test_gr<-selected_gr ##modify here to get different gr results

##Make color dataframe 
#Color Palette
  color_palette<-"Spectral"
  colors<-colorRampPalette(colors = brewer.pal(8, color_palette))
  
  #Gene
  gene_df<-as.data.frame(gene_gr)
  type_gene<-unique(gene_df$geneBiotype)
  threshold<-nrow(gene_df)*3/100 #3% of all genes
  
  for (type in type_gene) {
    if (nrow(gene_df[gene_df$geneBiotype == type,]) < threshold) {
      gene_df$geneBiotype[gene_df$geneBiotype == type]<-"other"
      
    }
    
  }
  
  type_gene_corrected<-unique(gene_df$geneBiotype)
  type_gene_count<-length(type_gene_corrected)
  color_gene_df<-as.data.frame(matrix(nrow = type_gene_count, ncol = 2))
  colnames(color_gene_df)<-c("type", "color")
  color_gene_df$color<-colors(type_gene_count)
  color_gene_df$type<-type_gene_corrected
  
  
  #Transcript
  transcript_df<-as.data.frame(transcript_gr)
  type_transcript<-unique(transcript_df$txBiotype)
  threshold<-nrow(transcript_df)*0.5/100 #3% of all transcripts
  
  for (type in type_transcript) {
    if (nrow(transcript_df[transcript_df$txBiotype == type,]) < threshold) {
      transcript_df$txBiotype[transcript_df$txBiotype == type]<-"other"
      
    }
    
  }
  
  type_transcript_corrected<-unique(transcript_df$txBiotype)
  type_transcript_count<-length(type_transcript_corrected)
  color_transcript_df<-as.data.frame(matrix(nrow = type_transcript_count, ncol = 2))
  colnames(color_transcript_df)<-c("type", "color")
  color_transcript_df$color<-colors(type_transcript_count)
  color_transcript_df$type<-type_transcript_corrected
  
##Get nearest feature
  #Gene
  nearest_gene<-distanceToNearest(test_gr, gene_gr)
  test_df<-as.data.frame(nearest_gene)
  test_df$distance=test_df$distance +1
  
  #Transcript
  nearest_transcript<-distanceToNearest(test_gr, transcript_gr)
  nearest_transcript_df<-as.data.frame(nearest_transcript)
  nearest_transcript_df$distance=nearest_transcript_df$distance +1
  
  ##Make Plots
  ggplot(test_df, aes(x = distance)) + geom_density() + 
    scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + theme_classic() + xlab(label = "Distance to nearest gene (kb)")
  
  ggplot(nearest_transcript_df, aes(x = distance)) + geom_density() + 
    scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + theme_classic() + xlab(label = "Distance to nearest transcript (kb)")

  ##Get type of nearest feature [transcript and gene]
  #Choose gene or transcript
  test<-gene
  
  if(test == "gene"){
    test_df<-nearest_gene_df
    test_gr<-gene_gr
    color_df<-color_gene_df
    
    }
    elif(test == "transcript"){
    test_df<-nearest_transcript_df
    test_gr<-transcript_gr
    color_df<-color_transcript_df
    }
  
  
  nearest_type_gr<-test_gr[test_df$subjectHits,]
  nearest_type<-as.data.frame(nearest_type_gr$geneBiotype)
  nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
  colnames(nearest_type_count)<-c("type", "count")
  nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
  row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))
  
  for (n in 1:nrow(nearest_type_count)){
    type=nearest_type_count$type[n]
    number<-length(nearest_type[nearest_type$value == type,])
    nearest_type_count$count[n]<-number
  }

  other_types<-nearest_type_count$type[nearest_type_count$type %notin% color_df$type]
  other_types_count<-sum(nearest_type_count$count[nearest_type_count$type %in% other_types])
  type_df<-nearest_type_count[nearest_type_count$type %in% color_df$type,]
  others<-c("other", other_types_count)
  names(others)<-c("type", "count")
  
  nearest_type_count<-rbind(nearest_type_count, others)
  nearest_type_count_df<-nearest_type_count[nearest_type_count$type %in% c(color_df$type, "other"),]
  nearest_type_count_df<-nearest_type_count_df[order(nearest_type_count_df$count),]
  nearest_type_plot_df<-merge.data.frame(nearest_type_count_df, color_df, by.all = "type")
  nearest_type_plot_df$count<-as.numeric(nearest_type_plot_df$count)
  nearest_type_plot_df$type<-factor(nearest_type_plot_df$type, levels = nearest_type_plot_df$type[order(nearest_type_plot_df$count)])
  

  ggplot(nearest_type_plot_df[order(nearest_type_plot_df$count, decreasing = TRUE),], aes(x = "", y = count, fill = type)) + 
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", direction = 1) + theme_bw() + 
    scale_fill_manual(values = nearest_type_plot_df$color[order(nearest_type_plot_df$count, decreasing = TRUE)])
  
  
## close to features

close_nearest_df<-nearest_df[nearest_df$distance<1000,]
nearest_selected_gr<-gene_gr[close_nearest_df$subjectHits,]
nearest_type<-as.data.frame(nearest_selected_gr$geneBiotype)
nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
colnames(nearest_type_count)<-c("type", "count")
nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))

for (n in 1:nrow(nearest_type_count)){
  type=nearest_type_count$type[n]
  number<-length(nearest_type[nearest_type$value == type,])
  nearest_type_count$count[n]<-number
}
colors<-colorRampPalette(brewer.pal(8, "Spectral"))

ggplot(nearest_type_count, aes(x = "type", y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  coord_polar("y") + theme_bw() + scale_fill_manual(values = colors(nrow(nearest_type_count)))


##far from features

far_nearest_df<-nearest_df[nearest_df$distance>1000,]
nearest_selected_gr<-gene_gr[far_nearest_df$subjectHits,]
nearest_type<-as.data.frame(nearest_selected_gr$geneBiotype)
nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
colnames(nearest_type_count)<-c("type", "count")
nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))

for (n in 1:nrow(nearest_type_count)){
  type=nearest_type_count$type[n]
  number<-length(nearest_type[nearest_type$value == type,])
  nearest_type_count$count[n]<-number
}
colors<-colorRampPalette(brewer.pal(8, "Spectral"))

ggplot(nearest_type_count, aes(x = "type", y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  coord_polar("y") + theme_bw() + scale_fill_manual(values = colors(nrow(nearest_type_count)))

##Same thing for all locations (significant)

##Density nearest
colnames(peak_locations)<-c("chr", "start", "end", "peak_ID", "score", "strand", "peakName")
all_gr<-peak_locations
all_gr$chr<-sapply(str_split(all_gr$chr, pattern = "chr"), "[", 2)
all_gr<-makeGRangesFromDataFrame(all_gr)
nearest_all<-distanceToNearest(all_gr, gene_gr)

nearest_all_df<-as.data.frame(nearest_all)
nearest_all_df$distance<-nearest_all_df$distance + 1
ggplot(nearest_all_df, aes(x = distance)) + geom_density() + 
  scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + theme_classic() + xlab(label = "Distance to nearest gene")  +
  theme(text = element_text(family = "Arial", size = 40)) + ylim(0, 2)

ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/Peak_distanceToFeature_sigPeaks.tiff", device = "tiff", dpi = 300)

## feature type
nearest_all_gr<-gene_gr[nearest_all_df$subjectHits,]
nearest_type<-as.data.frame(nearest_all_gr$geneBiotype)
nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
colnames(nearest_type_count)<-c("type", "count")
nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))

for (n in 1:nrow(nearest_type_count)){
  type=nearest_type_count$type[n]
  number<-length(nearest_type[nearest_type$value == type,])
  nearest_type_count$count[n]<-number
}
colors<-colorRampPalette(brewer.pal(8, "Spectral"))

ggplot(nearest_type_count, aes(x = "type", y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  coord_polar("y") + theme_void() + scale_fill_manual(values = colors(nrow(nearest_type_count))) + theme(text = element_text(family = "Arial", size = 40))


ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/NearestFeatureType_sigPeaks.tiff", device = "tiff",
        dpi = 300)
##far from features - all

far_nearest_df<-nearest_all_df[nearest_all_df$distance>1000,]
nearest_selected_gr<-gene_gr[far_nearest_df$subjectHits,]
nearest_type<-as.data.frame(nearest_selected_gr$geneBiotype)
nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
colnames(nearest_type_count)<-c("type", "count")
nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))

for (n in 1:nrow(nearest_type_count)){
  type=nearest_type_count$type[n]
  number<-length(nearest_type[nearest_type$value == type,])
  nearest_type_count$count[n]<-number
}
colors<-colorRampPalette(brewer.pal(8, "Spectral"))

ggplot(nearest_type_count, aes(x = "type", y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  coord_polar("y") + theme_void() + scale_fill_manual(values = colors(nrow(nearest_type_count))) + theme(text = element_text(family = "Arial", size = 40))


ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/NearestFeatureTypeFar_sigPeaks.tiff", device = "tiff",
     dpi = 300)

##near to features - all

near_nearest_df<-nearest_all_df[nearest_all_df$distance<1000,]
nearest_selected_gr<-gene_gr[near_nearest_df$subjectHits,]
nearest_type<-as.data.frame(nearest_selected_gr$geneBiotype)
nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
colnames(nearest_type_count)<-c("type", "count")
nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))

for (n in 1:nrow(nearest_type_count)){
  type=nearest_type_count$type[n]
  number<-length(nearest_type[nearest_type$value == type,])
  nearest_type_count$count[n]<-number
}
colors<-colorRampPalette(brewer.pal(8, "Spectral"))

ggplot(nearest_type_count, aes(x = "type", y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  coord_polar("y") + theme_void() + scale_fill_manual(values = colors(nrow(nearest_type_count))) + theme(text = element_text(family = "Arial", size = 40))


ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/NearestFeatureTypeClose_sigPeaks.tiff", device = "tiff",
       dpi = 300)

##Get density of nearest feature for whole peakset [TRANSCRIPT GTF]

nearest_all<-distanceToNearest(peakset_gr, transcript_gr)

nearest_all_df<-as.data.frame(nearest_all)
nearest_all_df$distance<-nearest_all_df$distance + 1
ggplot(nearest_all_df, aes(x = distance)) + geom_density() + 
  scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + theme_classic() + xlab(label = "Distance to nearest transcript (kb)")  +
  theme(text = element_text(family = "Arial", size = 40))  + ylim(0, 2)

ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/Peak_distanceToFeature_peakset.tiff", device = "tiff",
        dpi = 300)

##Get density of nearest feature for whole peakset [GENE GTF]
peakset_gr<-as.data.frame(peakset)
peakset_gr<-makeGRangesFromDataFrame(peakset_gr)
nearest_all<-distanceToNearest(peakset_gr, transcript_gr)

nearest_all_df<-as.data.frame(nearest_all)
nearest_all_df$distance<-nearest_all_df$distance + 1
ggplot(nearest_all_df, aes(x = distance)) + geom_density() + 
  scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + theme_classic() + xlab(label = "Distance to nearest transcript")  +
  theme(text = element_text(family = "Arial", size = 40))  + ylim(0, 2)


## feature type
nearest_all_gr<-gene_gr[nearest_all_df$subjectHits,]
nearest_type<-as.data.frame(nearest_all_gr$geneBiotype)
nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
colnames(nearest_type_count)<-c("type", "count")
nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))

for (n in 1:nrow(nearest_type_count)){
  type=nearest_type_count$type[n]
  number<-length(nearest_type[nearest_type$value == type,])
  nearest_type_count$count[n]<-number
}
colors<-colorRampPalette(brewer.pal(8, "Spectral"))

ggplot(nearest_type_count, aes(x = "type", y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  coord_polar("y") + theme_void() + scale_fill_manual(values = colors(nrow(nearest_type_count))) + theme(text = element_text(family = "Arial", size = 40))


ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/NearestFeatureType_peakset.tiff", device = "tiff",
       dpi = 300)
##far from features - all

far_nearest_df<-nearest_all_df[nearest_all_df$distance>1000,]
nearest_selected_gr<-gene_gr[far_nearest_df$subjectHits,]
nearest_type<-as.data.frame(nearest_selected_gr$geneBiotype)
nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
colnames(nearest_type_count)<-c("type", "count")
nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))

for (n in 1:nrow(nearest_type_count)){
  type=nearest_type_count$type[n]
  number<-length(nearest_type[nearest_type$value == type,])
  nearest_type_count$count[n]<-number
}
colors<-colorRampPalette(brewer.pal(8, "Spectral"))

ggplot(nearest_type_count, aes(x = "type", y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  coord_polar("y") + theme_void() + scale_fill_manual(values = colors(nrow(nearest_type_count))) + theme(text = element_text(family = "Arial", size = 40))


ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/NearestFeatureTypeFar_peakset.tiff", device = "tiff",
       dpi = 300)

##near to features - all

near_nearest_df<-nearest_all_df[nearest_all_df$distance<1000,]
nearest_selected_gr<-gene_gr[near_nearest_df$subjectHits,]
nearest_type<-as.data.frame(nearest_selected_gr$geneBiotype)
nearest_type_count<-as.data.frame(matrix(nrow=length(unique(nearest_type$value)), ncol = 2))
colnames(nearest_type_count)<-c("type", "count")
nearest_type_count[,1]<-as.vector(unique(nearest_type$value))
row.names(nearest_type_count)<-as.vector(unique(nearest_type$value))

for (n in 1:nrow(nearest_type_count)){
  type=nearest_type_count$type[n]
  number<-length(nearest_type[nearest_type$value == type,])
  nearest_type_count$count[n]<-number
}
colors<-colorRampPalette(brewer.pal(8, "Spectral"))

ggplot(nearest_type_count, aes(x = "type", y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  coord_polar("y") + theme_void() + scale_fill_manual(values = colors(nrow(nearest_type_count))) +
  theme(text = element_text(family = "Arial", size = 40))


ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/NearestFeatureTypeClose_peakset.tiff", device = "tiff", dpi = 300)


### Get closest tss [BY TRANSCRIPT]

tss_transcript<-transcript_gr
end(tss_transcript[strand(tss_transcript)=="+",])  =start(tss_transcript[strand(tss_transcript)=="+",])
start(tss_transcript[strand(tss_transcript)=="-",])=end(tss_transcript[strand(tss_transcript)=="-",])
tss_transcript=tss_transcript[!duplicated(tss_transcript),]

test_gr<-peakset_gr ##Change this for other tested peakset
nearest_tss<-distanceToNearest(test_gr, tss_transcript) 
nearest_tss_df<-as.data.frame(nearest_tss)
nearest_tss_df$distance<-nearest_tss_df$distance + 1

ggplot(nearest_tss_df, aes(x = distance)) + geom_density() + 
  scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + theme_classic() + xlab(label = "Distance to nearest tss")  +
  theme(text = element_text(family = "Arial", size = 40))  + ylim(0, 2)
ggsave(plot = last_plot(), filename = "/Users/ehresms/computational/ATAC/diffbind/DistanceToTSS_peakset.tiff", device = "tiff", dpi = 300)

nearest_tss_df$distance_transcript<-nearest_all_df$distance

ggplot(nearest_tss_df, aes(x = distance, y = distance_transcript)) + stat_density_2d(geom = "raster", aes(fill = after_stat(density)), contour = FALSE) +
  xlab(label = "Distance to nearest tss") + ylab (label = "Distance to nearest transcript") +
  scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + theme_classic() + scale_y_log10(labels=trans_format("log10", math_format(10^.x)))


## Get annotation from homer [BY CLUSTER OR PEAKSET]

anno<-read.table("/Users/ehresms/computational/ATAC/diffbind/peakset_anno.txt", sep = "\t", header = TRUE) ##change here for different annotation file
colnames(anno)[1]<-"PeakID"
anno_df<-anno[,c("PeakID", "Annotation", "Gene.Type", "Distance.to.TSS")] ## make smaller dataframe
anno_df$Annotation<-sapply(str_split(anno_df$Annotation, pattern =" "), "[", 1) 
anno_df$Distance.to.TSS<-anno_df$Distance.to.TSS+1
anno_feature<-data.frame(matrix(nrow = length(unique(anno_df$Annotation)), ncol = 2)) ##make annotation count dataframe
colnames(anno_feature)<-c("Annotation", "Count")

distance_df<-anno_df[,c("PeakID", "Distance.to.TSS")]
distance_min<-min(distance_df$Distance.to.TSS)
distance_max<-max(distance_df$Distance.to.TSS)

distance_df$scaled_distance<-((distance_df$Distance.to.TSS - distance_min)/(distance_max - distance_min))

for (n in 1:nrow(anno_feature)) { ##fill annotation count dataframe
  feature = unique(anno_df$Annotation)[n]
  feature_count<-length(anno_df$Annotation[anno_df$Annotation == feature])
  anno_feature$Annotation[n]<-feature
  anno_feature$Count[n]<-feature_count
}

anno_type<-data.frame(matrix(nrow = length(unique(anno_df$Gene.Type)), ncol = 2)) ##make annotation count dataframe
colnames(anno_type)<-c("Gene.Type", "Count")

for (n in 1:nrow(anno_type)) { ##fill annotation count dataframe
  feature = unique(anno_df$Gene.Type)[n]
  feature_count<-length(anno_df$Gene.Type[anno_df$Gene.Type == feature])
  anno_type$Gene.Type[n]<-feature
  anno_type$Count[n]<-feature_count
}


ggplot(data = distance_df, aes(x = scaled_distance)) + geom_density()

  
scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + theme_classic() + scale_y_log10(labels=trans_format("log10", math_format(10^.x)))

ggplot(data = anno_feature, aes(x = "Annotation", y = Count, fill = Annotation)) + geom_bar(stat = "identity") +
  coord_polar("y")

ggplot(data = anno_type, aes(x = "Gene.Type", y = Count, fill = Gene.Type)) + geom_bar(stat = "identity") +
  coord_polar("y") + theme_classic()

ggsave(plot=last_plot(), "/Users/ehresms/computational/ATAC/diffbind/GeneType_Peakset.tiff", device= "tiff")
