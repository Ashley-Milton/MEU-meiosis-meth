#Milton et al.


# Setup, loading packages and methylation data ----------------------------

setwd("/path/to/working/dir")

options(scipen = 999)

#Loading packages
library(tidyverse)
library(GenomicRanges)
library(scales)
library(rcompanion)
library(geomtextpath)

#Reading in the data
data <- read.table(file = "data/Tammar_final_MEUv7_PercMeth_mincov1_150623.txt", header = T)

#Adding an empty strand column
data$strand <- "*"

#Averaging percent methylation for replicates of same cell type
data$spg <- (data$spg1+data$spg2)/2
data$pd <- (data$pd1+data$pd2)/2
data$rs <- (data$rs1+data$rs2)/2
data$sperm <- (data$sperm1+data$sperm2)/2

#Adding coverage values for replicates of same cell type
data$covspg <- data$covspg1+data$covspg2
data$covpd <- data$covpd1+data$covpd2
data$covrs <- data$covrs1+data$covrs2
data$covsperm <- data$covsperm1+data$covsperm2

# Loading and manipulating GFF --------------------------------------------

#Reading in GFF
GeneGFF <- read.delim("data/Tammar_Male_v7_Genes_Reformat.sorted.gff3", header=FALSE, comment.char="#")

#Filtering for genes
GeneGFF <- subset(GeneGFF, grepl("gene", V3))
GeneGFF <- GeneGFF[c(1,4,5,9,7)]
TSS <- GeneGFF

#Calculating 20 kb windows around gene start sites
TSS$V4 <- ifelse(GeneGFF$V7=="+", TSS$V4 - 10000, TSS$V4)
TSS$V5 <- ifelse(GeneGFF$V7=="+", TSS$V4 + 20000, TSS$V5)
TSS$V4 <- ifelse(GeneGFF$V7=="-", TSS$V5 - 10000, TSS$V4)
TSS$V5 <- ifelse(GeneGFF$V7=="-", TSS$V5 + 10000, TSS$V5)
TSS$V4 <- ifelse(TSS$V4<0, 0,TSS$V4)

colnames(TSS) <- c("chr", "start", "end", "gene","strand")


# GRanges intersect between methylation data and 20 kb TSS windows --------

df <- data
df1 <- TSS
gr <- makeGRangesFromDataFrame(df, TRUE)
gr1 <- makeGRangesFromDataFrame(df1, TRUE)
olaps <- findOverlaps(gr, gr1, type="any")
x <- left_join(tibble(queryHits = olaps@from, subjectHits = olaps@to), df %>% rowid_to_column(var = "queryHits"))
x <- left_join(x, df1 %>% rowid_to_column(var="subjectHits"), by = "subjectHits")

#Calculating position relative to TSS
x$posB <- ifelse(grepl("0",x$start.y,fixed = TRUE), 
                 ifelse(grepl("+",x$strand.y,fixed = TRUE), x$start.x - (x$end.y - 10000) , (x$end.y - 10000) - x$start.x),
                 ifelse(grepl("+",x$strand.y,fixed = TRUE), x$start.x - (x$start.y + 10000) , (x$start.y + 10000) - x$start.x)
)

#Adding a chromosome type column - X or autosome
x$c <- ifelse(x$chr.x=="chrX", "X", "A")

# Filtering for coverage >19 ----------------------------------------------

#Combined replicates
methylation <- x
methylation$spg <- ifelse(methylation$covspg1 > 19 & methylation$covspg2 > 19, methylation$spg, NA)
methylation$pd <- ifelse(methylation$covpd1 > 19 & methylation$covpd2 > 19, methylation$pd, NA)
methylation$rs <- ifelse(methylation$covrs1 > 19 & methylation$covrs2 > 19, methylation$rs, NA)
methylation$sperm <- ifelse(methylation$covsperm1 > 19 & methylation$covsperm2 > 19, methylation$sperm, NA)
methylation <- methylation %>% select(chr.x, start.x, end.x, spg, pd, rs, sperm, c)
methlong <- pivot_longer(methylation, 4:7)
methlong$name <- gsub('spg', 'Spg', methlong$name)
methlong$name <- gsub('pd', 'P/D', methlong$name)
methlong$name <- gsub('rs', 'RS', methlong$name)
methlong$name <- gsub('sperm', 'Sperm', methlong$name)

#Separate replicates
methylation2 <- x
methylation2$spg <- ifelse(methylation2$covspg1 > 19 & methylation2$covspg2 > 19, methylation2$spg, NA)
methylation2$pd <- ifelse(methylation2$covpd1 > 19 & methylation2$covpd2 > 19, methylation2$pd, NA)
methylation2$rs <- ifelse(methylation2$covrs1 > 19 & methylation2$covrs2 > 19, methylation2$rs, NA)
methylation2$sperm <- ifelse(methylation2$covsperm1 > 19 & methylation2$covsperm2 > 19, methylation2$sperm, NA)
methylation2 <- methylation2 %>% select(chr.x, start.x, end.x, spg1, pd1, rs1, sperm1, spg2, pd2, rs2, sperm2, c)
methlong2 <- pivot_longer(methylation2, 4:11)
methlong2$name <- gsub('spg1', 'Spg 1', methlong2$name)
methlong2$name <- gsub('pd1', 'P/D 1', methlong2$name)
methlong2$name <- gsub('rs1', 'RS 1', methlong2$name)
methlong2$name <- gsub('sperm1', 'Sperm 1', methlong2$name)
methlong2$name <- gsub('spg2', 'Spg 2', methlong2$name)
methlong2$name <- gsub('pd2', 'P/D 2', methlong2$name)
methlong2$name <- gsub('rs2', 'RS 2', methlong2$name)
methlong2$name <- gsub('sperm2', 'Sperm 2', methlong2$name)


# DNA methylation histogram -----------------------------------------------

#DNA methylation perc. histogram separate replicates
#FIG. S1(B)
methlong2 %>% na.omit() %>% 
  ggplot(aes(x=value)) +
  geom_histogram(fill="grey20") +
  theme_classic() +
  xlab("Methylation (%)") +
  ylab("Count (x10^5)")+
  facet_wrap(~factor(name, levels=c("Spg 1", "P/D 1", "RS 1", "Sperm 1", "Spg 2", "P/D 2", "RS 2", "Sperm 2")), nrow = 2)+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-5))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7)

ggsave(file = "figures/DNA_meth_perc_hist_FILTERED_v2.pdf", dpi=800, height = 8, width = 13)

#DNA methylation perc. histogram combined replicates
#FIG. 1(C)
methlong %>% na.omit() %>% 
  ggplot(aes(x=value)) +
  geom_histogram(fill="grey20") +
  theme_classic() +
  xlab("Methylation (%)") +
  ylab("Count (x10^5)")+
  facet_wrap(~factor(name, levels=c("Spg", "P/D", "RS", "Sperm")), nrow = 1)+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-5), n.breaks = 4)+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7)

ggsave(file = "figures/DNA_meth_perc_hist_combined_FILTERED.pdf", dpi=800, height = 3.5, width = 13)

# Coverage histogram ------------------------------------------------------

coverage2 <- x %>% select(chr.x, start.x, end.x, covspg1, covpd1, covrs1, covsperm1, covspg2, covpd2, covrs2, covsperm2, c)
covlong2 <- pivot_longer(coverage2, 4:11)
covlong2$name <- substring(covlong2$name, 4)
covlong2$name <- gsub('spg1', 'Spg 1', covlong2$name)
covlong2$name <- gsub('pd1', 'P/D 1', covlong2$name)
covlong2$name <- gsub('rs1', 'RS 1', covlong2$name)
covlong2$name <- gsub('sperm1', 'Sperm 1', covlong2$name)
covlong2$name <- gsub('spg2', 'Spg 2', covlong2$name)
covlong2$name <- gsub('pd2', 'P/D 2', covlong2$name)
covlong2$name <- gsub('rs2', 'RS 2', covlong2$name)
covlong2$name <- gsub('sperm2', 'Sperm 2', covlong2$name)

#FIG. S1(A)
covlong2 %>% 
  ggplot(aes(x=value)) +
  geom_histogram(fill="grey20", binwidth = 2) +
  xlim(c(0,70)) +
  theme_classic() +
  xlab("Coverage")+
  ylab("Count (x10^5)")+
  facet_wrap(~factor(name, levels=c("Spg 1", "P/D 1", "RS 1", "Sperm 1", "Spg 2", "P/D 2", "RS 2", "Sperm 2")), nrow = 2)+
  geom_vline(xintercept = 20, col="#eb1b24", linetype = "dashed")+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-5))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7)

ggsave(file = "figures/Cov_hist_v2.pdf", dpi=600, height = 8, width = 13)

# CpG sampling relative to TSS --------------------------------------------

plotdata <- x
plotdata$spg <- ifelse(plotdata$covspg1 > 19 & plotdata$covspg2 > 19, plotdata$spg, NA)
plotdata$pd <- ifelse(plotdata$covpd1 > 19 & plotdata$covpd2 > 19, plotdata$pd, NA)
plotdata$rs <- ifelse(plotdata$covrs1 > 19 & plotdata$covrs2 > 19, plotdata$rs, NA)
plotdata$sperm <- ifelse(plotdata$covsperm1 > 19 & plotdata$covsperm2 > 19, plotdata$sperm, NA)
plotdata <- plotdata %>% select(chr.x, start.x, end.x, spg, pd, rs, sperm, posB, c)
plotdata <- pivot_longer(plotdata, 4:7)
plotdata <- plotdata %>% na.omit()

#FIG. S1(E)
ggplot(plotdata, aes(posB, fill=c))+ geom_histogram()+
  scale_fill_manual(values = c("#cdccff", "#3e68db"), name = "", labels = c("Autosomes", "X chromosome"))+
  theme_classic()+
  xlab("Position relative to TSS (kb)")+
  ylab("Count (x10^3)")+
  facet_wrap(~c, scales = "free_y", labeller = labeller("A" = "Autosomes", "X" = "X chromosomes"))+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-3))+
  scale_x_continuous(labels = unit_format(unit="", scale = 1e-3))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7,
        legend.position = "bottom")

ggsave(file = "figures/CpG_relative_TSS_FILTERED_v2.pdf", dpi=800, height = 3.5, width = 6.5)

# Taking median value at each position relative to TSS --------------------

z <- x

meddata <- z %>% filter(covspg1 > 19 & covspg2 > 19) %>% select(posB, c, spg)
meddata <- meddata %>% drop_na() %>% group_by(posB, c) %>% summarise(medianpercmeth=median(spg))
meddata$sample <- "spg"
medianspg <- meddata

meddata <- z %>% filter(covpd1 > 19 & covpd2 > 19) %>% select(posB, c, pd)
meddata <- meddata %>% drop_na() %>% group_by(posB, c) %>% summarise(medianpercmeth=median(pd))
meddata$sample <- "pd"
medianpd <- meddata

meddata <- z %>% filter(covrs1 > 19 & covrs2 > 19) %>% select(posB, c, rs)
meddata <- meddata %>% drop_na() %>% group_by(posB, c) %>% summarise(medianpercmeth=median(rs))
meddata$sample <- "rs"
medianrs <- meddata

meddata <- z %>% filter(covsperm1 > 19 & covsperm2 > 19) %>% select(posB, c, sperm)
meddata <- meddata %>% drop_na() %>% group_by(posB, c) %>% summarise(medianpercmeth=median(sperm))
meddata$sample <- "sperm"
mediansperm <- meddata


#Putting them together
allmedians <- rbind(medianspg, medianpd, medianrs, mediansperm)
allmedians <- allmedians %>% ungroup()
allmedians$sample <- gsub('spg', 'Spg', allmedians$sample)
allmedians$sample <- gsub('pd', 'P/D', allmedians$sample)
allmedians$sample <- gsub('rs', 'RS', allmedians$sample)
allmedians$sample <- gsub('sperm', 'Sperm', allmedians$sample)

allmedians$sample_c <- paste(allmedians$sample, allmedians$c, sep = "_")

# TSS median plots --------------------------------------------------------

#FIG. 1(E)
allmedians %>% ggplot(aes(x=posB, y=medianpercmeth, color=c)) +
  geom_smooth(aes(fill=c), se=T, method="loess", span = 0.35) +
  # geom_point(size=0.10, alpha=0.4)+
  labs(x="Position relative to TSS (kb)", y="Median methylation (%)") + 
  ylim(0,100)+
  coord_cartesian(x=c(-10000,10000))+
  theme_classic()+
  scale_color_manual(values = c("#cdccff", "#3e68db"), name = "", labels = c("Autosomes", "X chromosome"))+
  scale_fill_manual(values = c("#cdccff", "#3e68db"), name = "", labels = c("Autosomes", "X chromosome"))+
  theme(legend.position="bottom")+
  facet_wrap(~factor(sample, levels = c("Spg", "P/D", "RS", "Sperm")), nrow=1)+
  scale_x_continuous(labels = unit_format(unit="", scale = 1e-3))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7)

ggsave(file = "figures/TSS_median_plots_v5.pdf", dpi=800, height = 6, width = 13)

#FIG. S1(C)
allmedians %>% ggplot(aes(x=posB, y=medianpercmeth, color=c)) +
  geom_point(size=0.10)+
  labs(x="Position relative to TSS (kb)", y="Median methylation (%)") + 
  ylim(0,100)+
  coord_cartesian(x=c(-10000,10000))+
  theme_classic()+
  scale_color_manual(values = c("#333333", "#3E68DB"), name = "", labels = c("Autosomes", "X chromosome"))+
  scale_fill_manual(values = c("#333333", "#3E68DB"), name = "", labels = c("Autosomes", "X chromosome"))+
  theme(legend.position="bottom")+
  facet_wrap(~factor(sample, levels = c("Spg", "P/D", "RS", "Sperm")), nrow=1)+
  scale_x_continuous(labels = unit_format(unit="", scale = 1e-3))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7)

ggsave(file = "figures/TSS_median_plots_points.pdf", dpi=800, height = 6, width = 13)


# Creating 500 row windows ------------------------------------------------

#Make separate dfs for each cell type
spg <- methylation %>% select(chr.x, start.x, end.x, spg, c)
pd <- methylation %>% select(chr.x, start.x, end.x, pd, c)
rs <- methylation %>% select(chr.x, start.x, end.x, rs, c)
sperm <- methylation %>% select(chr.x, start.x, end.x, sperm, c)

chunk_size <- 500

spg_summary <- spg %>% na.omit() %>% mutate(group_index = (row_number() - 1) %/% chunk_size) %>% 
  group_by(chr.x, group_index) %>% 
  summarize(count = n(), average_value = mean(spg, na.rm = T)) %>% select(-count) %>% mutate(name = "Spg")

pd_summary <- pd %>% na.omit() %>% mutate(group_index = (row_number() - 1) %/% chunk_size) %>% 
  group_by(chr.x, group_index) %>% 
  summarize(count = n(), average_value = mean(pd, na.rm = T)) %>% select(-count) %>% mutate(name = "P/D")

rs_summary <- rs %>% na.omit() %>% mutate(group_index = (row_number() - 1) %/% chunk_size) %>% 
  group_by(chr.x, group_index) %>% 
  summarize(count = n(), average_value = mean(rs, na.rm = T)) %>% select(-count) %>% mutate(name = "RS")

sperm_summary <- sperm %>% na.omit() %>% mutate(group_index = (row_number() - 1) %/% chunk_size) %>% 
  group_by(chr.x, group_index) %>% 
  summarize(count = n(), average_value = mean(sperm, na.rm = T)) %>% select(-count) %>% mutate(name = "Sperm")

#Add together, add levels to name
windows_boxplot_data <- rbind(spg_summary, pd_summary, rs_summary, sperm_summary)
windows_boxplot_data$name <- factor(windows_boxplot_data$name, levels=c("Spg", "P/D", "RS", "Sperm"))

#Adding c (chromosome type) column
windows_boxplot_data$c <- ifelse(windows_boxplot_data$chr.x=="chrX", "X", "A")

windows_boxplot_data$c <- gsub('A', 'Autosomes', windows_boxplot_data$c)
windows_boxplot_data$c <- gsub('X', 'X chromosome', windows_boxplot_data$c)

# 500 row window boxplots -------------------------------------------------

#FIG. 1(D)
windows_boxplot_data %>% ggplot(aes(y=average_value, x=c, fill=c))+
  geom_boxplot(notch = T)+
  xlab("")+
  ylab("Mean methylation (%)")+
  theme_classic()+
  theme(legend.position = "bottom")+
  facet_wrap(~name, nrow = 1)+
  scale_fill_manual(values = c("#cdccff", "#3e68db"))+
  stat_summary(fun.y=mean, geom = "point", show.legend = F, size=0.5)+
  theme(legend.text = element_text(size = 14),
        legend.title = element_blank(),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 0.7)

# ggsave(file = "figures/500_CpG_boxplots_newmethod_FILTERED.pdf", dpi=800, height = 5, width = 13)

#Mood's median test on this data

windows_boxplot_data$name_c <- paste(windows_boxplot_data$name, windows_boxplot_data$c, sep = "_")

moods <- pairwiseMedianTest(average_value~name_c, data=windows_boxplot_data, method = "holm")

moods_sig <- moods %>% filter(p.adjust < 0.05)

#write.table(moods, file="data/Moods_median_500_row_boxplots.txt", quote=F, sep="\t", row.names = F, col.names = T)


# Filtering for coverage >9 -----------------------------------------------

filtered_methylation <- data

#Filtering for coverage of at least 10 in individual replicates - leaves an NA if cov not sufficient
filtered_methylation$spg <- ifelse(filtered_methylation$covspg1 > 9 & filtered_methylation$covspg2 > 9, filtered_methylation$spg, NA)
filtered_methylation$pd <- ifelse(filtered_methylation$covpd1 > 9 & filtered_methylation$covpd2 > 9, filtered_methylation$pd, NA)
filtered_methylation$rs <- ifelse(filtered_methylation$covrs1 > 9 & filtered_methylation$covrs2 > 9, filtered_methylation$rs, NA)
filtered_methylation$sperm <- ifelse(filtered_methylation$covsperm1 > 9 & filtered_methylation$covsperm2 > 9, filtered_methylation$sperm, NA)

#Simplifying and making long format
filtered_methylation <- filtered_methylation %>% select(chr, start, end, spg, pd, rs, sperm)
filtered_methlong <- pivot_longer(filtered_methylation, 4:7)

#Renaming samples for use in plots
filtered_methlong$name <- gsub('spg', 'Spg', filtered_methlong$name)
filtered_methlong$name <- gsub('pd', 'P/D', filtered_methlong$name)
filtered_methlong$name <- gsub('rs', 'RS', filtered_methlong$name)
filtered_methlong$name <- gsub('sperm', 'Sperm', filtered_methlong$name)

# Rsx locus plots ---------------------------------------------------------

#Rsx start and end positions obtained by BLASTing Monodelphis domestica Rsx to this Tammar wallaby genome

plotdata <- filtered_methlong %>% na.omit() %>% filter(chr == "chrX" & start > 54720626-15000 & end < 54759097+15000)

#FIG. 1(F)
ggplot(data=plotdata, aes(x=start, y=value))+
  geom_smooth(se=F, method = "loess", span = 0.25, color = "#3e68db", fill = "#3e68db")+
  facet_wrap(~factor(name, levels = c("Spg", "P/D", "RS", "Sperm")), nrow=1)+
  theme_classic()+
  ylab("Methylation (%)")+
  xlab("Position on X chromosome (kb)")+
  scale_x_continuous(labels = unit_format(unit="", scale = 1e-3), breaks = seq(54720626-6000+5000, 54759097+9500, by = 15000))+
  scale_y_continuous(limits = c(-8, 100))+
  geom_textsegment(aes(x=54720626, y=-5, xend=54759097, yend=-5, label = "Rsx"), size=4.5, arrow = arrow(end = "first", type = "closed", length = unit(0.1, units = "inches")))+
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 16),
        aspect.ratio = 0.7)+
  coord_cartesian(xlim = c(54720626-6000, 54759097+9500))


ggsave(file = "figures/Rsx_meth_plots_points_v2.pdf", dpi=800, height = 6, width = 13)

#FIG. S1(D)
ggplot(data=plotdata, aes(x=start, y=value))+
  geom_point(color = "#3e68db", fill = "#3e68db")+
  facet_wrap(~factor(name, levels = c("Spg", "P/D", "RS", "Sperm")), nrow=1)+
  theme_classic()+
  ylab("Methylation (%)")+
  xlab("Position on X chromosome (kb)")+
  scale_x_continuous(labels = unit_format(unit="", scale = 1e-3), breaks = seq(54720626-6000+5000, 54759097+9500, by = 15000))+
  scale_y_continuous(limits = c(-8, 100))+
  geom_textsegment(aes(x=54720626, y=-5, xend=54759097, yend=-5, label = "Rsx"), size=4.5, arrow = arrow(end = "first", type = "closed", length = unit(0.1, units = "inches")))+
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 16),
        aspect.ratio = 0.7)+
  coord_cartesian(xlim = c(54720626-6000, 54759097+9500))


ggsave(file = "figures/Rsx_meth_plots_points_v2.pdf", dpi=800, height = 6, width = 13)
