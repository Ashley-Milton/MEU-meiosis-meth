#### Setting working directory and loading packages ####

setwd("/path/to/working/dir")

library(tidyverse)
library(methylKit)
library(factoextra)

#### Reading in CpG report files ####

file.list <- list("HY7HJDSX2_F11-NEB-UDI_merged.deduplicated.CpG_report.txt",
                  "HY7HJDSX2_G11-NEB-UDI_merged.deduplicated.CpG_report.txt",
                  "HY7HJDSX2_H11-NEB-UDI_merged.deduplicated.CpG_report.txt",
                  "HY7HJDSX2_A12-NEB-UDI_merged.deduplicated.CpG_report.txt",
                  "HCTJJDSX3_1-NEB-UDI_merged.deduplicated.CpG_report.txt",
                  "HCTJJDSX3_2-NEB-UDI_merged.deduplicated.CpG_report.txt",
                  "HCTJJDSX3_3-NEB-UDI_merged.deduplicated.CpG_report.txt",
                  "HCTJJDSX3_4-NEB-UDI_merged.deduplicated.CpG_report.txt")


all.samples=methRead(file.list,
                     sample.id = list("spg1","pd1","rs1","sperm1","spg2","pd2","rs2","sperm2"), 
                     assembly = "MEU_final_v7", 
                     treatment = c(0,1,2,3,0,1,2,3), 
                     context = "CpG",
                     pipeline = 'bismarkCytosineReport',
                     mincov = 1)

#### Filtering, normalisation and unite ####

filtered.text=filterByCoverage(all.samples, lo.perc = NULL, hi.count = NULL, hi.perc = NULL, lo.count = NULL) 
normalised.reads=normalizeCoverage(filtered.text)

meth=methylKit::unite(normalised.reads, destrand = T) #unites data
head(meth)

unique(meth$chr)

#### Putting split up chromosomes back together ####

#Need to change chromosome names and also adjust position values for 1.B and 2.B

#First: adjust positions
meth$start <- ifelse(meth$chr=="chr1.B", meth$start+400000000, meth$start)
meth$end <- ifelse(meth$chr=="chr1.B", meth$end+400000000, meth$end)
meth$start <- ifelse(meth$chr=="chr2.B", meth$start+300000000, meth$start)
meth$end <- ifelse(meth$chr=="chr2.B", meth$end+300000000, meth$end)

# Filter for a specific chromosome pattern - checking to see that additions worked, and everything is as it should be
filtered_df <- meth[grepl("^chr2$", meth$chr), ]

# Get the range of values for the 'start' column
start_range <- range(filtered_df$start)

# Print the range
print(start_range)

#Second: adjust chr names
meth$chr <- gsub('chr1.B', 'chr1', meth$chr)
meth$chr <- gsub('chr2.B', 'chr2', meth$chr)
unique(meth$chr)

#### Creating the percent methylation dataframe ####

##Getting coverage values

covspg1 <- meth$coverage1
covpd1 <- meth$coverage2
covrs1 <- meth$coverage3
covsperm1 <- meth$coverage4
covspg2 <- meth$coverage5
covpd2 <- meth$coverage6
covrs2 <- meth$coverage7
covsperm2 <- meth$coverage8

#Using percent methylation function then converting output to dataframe

percMeth <- percMethylation(meth, rowids=T, save.txt = T)

percMeth1 <- as.data.frame(percMeth)

percMeth1 <- tibble::rownames_to_column(percMeth1, var = "a")

percMeth1 <- percMeth1 %>% separate(a, c("chr", "start", "end"))

#Adding coverage values back to percent methylation dataframe

percMeth1$covpd1 <- covpd1
percMeth1$covpd2 <- covpd2
percMeth1$covrs1 <- covrs1
percMeth1$covrs2 <- covrs2
percMeth1$covsperm1 <- covsperm1
percMeth1$covsperm2 <- covsperm2
percMeth1$covspg1 <- covspg1
percMeth1$covspg2 <- covspg2

#Reordering the columns

percMeth2 <- percMeth1 %>% relocate(covpd1, .after = pd1)
percMeth2 <- percMeth2 %>% relocate(covpd2, .after = pd2)
percMeth2 <- percMeth2 %>% relocate(covspg1, .after = spg1)
percMeth2 <- percMeth2 %>% relocate(covspg2, .after = spg2)
percMeth2 <- percMeth2 %>% relocate(covsperm1, .after = sperm1)
percMeth2 <- percMeth2 %>% relocate(covsperm2, .after = sperm2)
percMeth2 <- percMeth2 %>% relocate(covrs1, .after = rs1)
percMeth2 <- percMeth2 %>% relocate(covrs2, .after = rs2)
head(percMeth2)

#Saving the dataframe
write.table(percMeth2, file = "/path/to/out/dir/Tammar_final_MEUv7_PercMeth_mincov1_150623.txt", quote = F, row.names = F, sep = "\t", col.names = T)


#### Sequencing stats ####

##Coverage histograms

#Make smaller dataframe
chr <- meth$chr
start <- meth$start
end <- meth$end
spg1 <- meth$coverage1
pd1 <- meth$coverage2
rs1 <- meth$coverage3
sperm1 <- meth$coverage4
spg2 <- meth$coverage5
pd2 <- meth$coverage6
rs2 <- meth$coverage7
sperm2 <- meth$coverage8
coverage <- data.frame(chr, start, end, spg1, pd1, rs1, sperm1, spg2, pd2, rs2, sperm2)

#Make it long format
covlong <- pivot_longer(coverage, 4:11)

#Assign colours to samples
colours <- data.frame(
  name = c("spg1", "pd1", "rs1", "sperm1", "spg2", "pd2", "rs2", "sperm2"),
  col = c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#6E016B")
)


covlong <- covlong %>% 
  left_join(colours, by = "name")
covlong$col <- paste0("'", covlong$col, "'")


covlong %>%
  ggplot(aes(x=value, fill=factor(name, levels=c("spg1", "spg2", "pd1", "pd2", "rs1", "rs2", "sperm1", "sperm2")))) +
  geom_histogram(colour="black", binwidth = 2) +
  xlim(c(0,80)) +
  theme_bw() +
  ggtitle("Coverage histogram") +
  scale_fill_manual(values = c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#6E016B")) +
  facet_wrap(~factor(name, levels=c("spg1", "pd1", "rs1", "sperm1", "spg2", "pd2", "rs2", "sperm2")), nrow = 2)+
  theme(legend.position = "none")

ggsave("/path/to/out/dir/coverage_histogram_lo.count1_final_MEUv7.png", dpi = 300, height = 8, width = 20, units = "cm")


##% CpG methylation histograms 
png("/path/to/out/dir/Hist_perc_CpG_meth_final_MEUv7.png", res = 300, height = 6, width = 10, units = "in")
par(mfrow = c(2, 4))
getMethylationStats(all.samples[[1]],plot=T,both.strands = FALSE) #Histogram of % CpG methylation for sample 1 - Change [[1]] to other numbers to see other samples
getMethylationStats(all.samples[[2]],plot=T,both.strands = FALSE)
getMethylationStats(all.samples[[3]],plot=T,both.strands = FALSE)
getMethylationStats(all.samples[[4]],plot=T,both.strands = FALSE)
getMethylationStats(all.samples[[5]],plot=T,both.strands = FALSE)
getMethylationStats(all.samples[[6]],plot=T,both.strands = FALSE)
getMethylationStats(all.samples[[7]],plot=T,both.strands = FALSE)
getMethylationStats(all.samples[[8]],plot=T,both.strands = FALSE)
dev.off()

##CpG coverage histograms
png("/path/to/out/dir/Hist_CpG_cov_final_MEUv7.png", res = 300, height = 6, width = 10, units = "in")
par(mfrow = c(2, 4))
getCoverageStats(all.samples[[1]],plot=TRUE,both.strands = FALSE) #Histogram of CpG Coverage for sample 1
getCoverageStats(all.samples[[2]],plot=TRUE,both.strands = FALSE)
getCoverageStats(all.samples[[3]],plot=TRUE,both.strands = FALSE)
getCoverageStats(all.samples[[4]],plot=TRUE,both.strands = FALSE)
getCoverageStats(all.samples[[5]],plot=TRUE,both.strands = FALSE)
getCoverageStats(all.samples[[6]],plot=TRUE,both.strands = FALSE)
getCoverageStats(all.samples[[7]],plot=TRUE,both.strands = FALSE)
getCoverageStats(all.samples[[8]],plot=TRUE,both.strands = FALSE)
dev.off()

##Sample correlation results
getCorrelation(meth,plot=F) #scatterplot of correlation between samples

####Clustering samples#
clusterSamples(meth, dist="correlation", method="ward.D", plot=FALSE)
clusterSamples(meth, dist="correlation", method="ward", plot=T) #plots a dendrogram

png("/path/to/out/dir/Dendrogram_final_MEUv7.png", res = 300, height = 6, width = 8, units = "in")
clusterSamples(meth, dist="correlation", method="ward", plot=T)
dev.off()

##PCA
par(mfrow = c(1, 1))

#Saving PCA
png("/path/to/out/dir/PCA_plot_final_MEUv7.png", res = 300, height = 6, width = 8, units = "in")
PCASamples(meth)
dev.off()

PCASamples(meth, screeplot=T) #plots screeplot
pca <- PCASamples(meth, obj.return = T) #PCA Plot, saves as object

#Better screeplot with % variance labels
fviz_eig(pca, addlabels = T)

##Saving screeplot
png("/path/to/out/dir/Screeplot_pca_variances_final_MEUv7.png", res = 300, height = 6, width = 8, units = "in")
fviz_eig(pca, addlabels = T) #better screeplot with % variance labels
dev.off()
