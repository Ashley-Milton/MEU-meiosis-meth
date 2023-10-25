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

#### Putting split up chromosomes back together ####

#Need to change chromosome names and also adjust position values for 1.B and 2.B

#First: adjust positions
meth$start <- ifelse(meth$chr=="chr1.B", meth$start+400000000, meth$start)
meth$end <- ifelse(meth$chr=="chr1.B", meth$end+400000000, meth$end)
meth$start <- ifelse(meth$chr=="chr2.B", meth$start+300000000, meth$start)
meth$end <- ifelse(meth$chr=="chr2.B", meth$end+300000000, meth$end)

#Second: adjust chr names
meth$chr <- gsub('chr1.B', 'chr1', meth$chr)
meth$chr <- gsub('chr2.B', 'chr2', meth$chr)
unique(meth$chr)

#### Creating the summary Cs and Ts dataframe ####

methdf <- getData(meth)
methdf$c <- ifelse(methdf$chr == "chrX", "X", "A")

#Summing total numCs and numTs for each cell type
summary <- methdf %>% group_by(c) %>% summarise(numCs_spg1 = sum(numCs1),
                   numTs_spg1 = sum(numTs1),
                   numCs_pd1 = sum(numCs2),
                   numTs_pd1 = sum(numTs2),
                   numCs_rs1 = sum(numCs3),
                   numTs_rs1 = sum(numTs3),
                   numCs_sperm1 = sum(numCs4),
                   numTs_sperm1 = sum(numTs4),
                   numCs_spg2 = sum(numCs5),
                   numTs_spg2 = sum(numTs5),
                   numCs_pd2 = sum(numCs6),
                   numTs_pd2 = sum(numTs6),
                   numCs_rs2 = sum(numCs7),
                   numTs_rs2 = sum(numTs7),
                   numCs_sperm2 = sum(numCs8),
                   numTs_sperm2 = sum(numTs8)) %>% 
  pivot_longer(2:17)
  
summary <- summary %>% separate(name, into = c("type", "name"), sep = "_") %>% pivot_wider(names_from = "type", values_from = "value") %>% mutate(percmeth = numCs/(numCs+numTs))

summary$name <- gsub('spg1', 'Spg 1', summary$name)
summary$name <- gsub('pd1', 'P/D 1', summary$name)
summary$name <- gsub('rs1', 'RS 1', summary$name)
summary$name <- gsub('sperm1', 'Sperm 1', summary$name)
summary$name <- gsub('spg2', 'Spg 2', summary$name)
summary$name <- gsub('pd2', 'P/D 2', summary$name)
summary$name <- gsub('rs2', 'RS 2', summary$name)
summary$name <- gsub('sperm2', 'Sperm 2', summary$name)

write.table(summary, file="/path/to/out/dir/methylKit_summary_meth.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#Plot

summary %>% ggplot(aes(x = c, y = percmeth, fill=c))+
  geom_col()+
  facet_wrap(~factor(name, levels=c("Spg 1", "P/D 1", "RS 1", "Sperm 1", "Spg 2", "P/D 2", "RS 2", "Sperm 2")), nrow = 1)+
  theme_classic()+
  scale_fill_manual(values = c("#cdccff", "#3e68db"), name = "", labels = c("Autosomes", "X chromosome"))+
  theme(axis.title.x = element_blank())

ggsave("/path/to/out/dir/methylKit_summary_meth.pdf", dpi=800, width = 10, height = 2.5, units = "in")

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