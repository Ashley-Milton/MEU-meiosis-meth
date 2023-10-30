#Milton et al.


# Setup, loading packages and methylation data ----------------------------

setwd("/path/to/working/dir")

options(scipen = 999)

#Loading packages
library(tidyverse)
library(GenomicRanges)
library(scales)
library(geomtextpath)
library(gridExtra)
library(rcompanion)

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

#Adding an autosome or X chromosome column
data$c <- ifelse(data$chr=="chrX", "X", "A")

data_original <- data

#Re-naming cell types
colnames(data) <- gsub("spg", "Spg ", colnames(data))
colnames(data) <- gsub("pd", "P/D ", colnames(data))
colnames(data) <- gsub("rs", "RS ", colnames(data))
colnames(data) <- gsub("sperm", "Sperm ", colnames(data))

#Removing unnecessary space from some column names
data <- data %>% dplyr::rename('Spg' = 'Spg ', 'covSpg' = 'covSpg ',
                'P/D' = 'P/D ', 'covP/D' = 'covP/D ',
                'RS' = 'RS ', 'covRS' = 'covRS ',
                'Sperm' = 'Sperm ', 'covSperm' = 'covSperm ')

#Checking
colnames(data)

# Filtering original dataframe for coverage -------------------------------

#Combined replicates unfiltered
ufdata <- data
ufdata <- ufdata %>% dplyr::select(chr, start, end, Spg, 'P/D', RS, Sperm, c)
uflong <- pivot_longer(ufdata, 4:7)

#Separate replicates unfiltered
ufdata2 <- data
ufdata2 <- ufdata2 %>% dplyr::select(chr, start, end, 'Spg 1', 'P/D 1', 'RS 1', 'Sperm 1', 'Spg 2', 'P/D 2', 'RS 2', 'Sperm 2', c)
uflong2 <- pivot_longer(ufdata2, 4:11)

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

# Filtering intersect dataframe for coverage ------------------------------

#Combined replicates
methylation <- x
methylation$Spg <- ifelse(methylation$'covSpg 1' > 14 & methylation$'covSpg 2' > 14, methylation$Spg, NA)
methylation$'P/D' <- ifelse(methylation$'covP/D 1' > 14 & methylation$'covP/D 2' > 14, methylation$'P/D', NA)
methylation$RS <- ifelse(methylation$'covRS 1' > 14 & methylation$'covRS 2' > 14, methylation$RS, NA)
methylation$Sperm <- ifelse(methylation$'covSperm 1' > 14 & methylation$'covSperm 2' > 14, methylation$Sperm, NA)
methylation <- methylation %>% select(chr.x, start.x, end.x, Spg, 'P/D', RS, Sperm, posB, c)
methlong <- pivot_longer(methylation, 4:7)

# DNA methylation histograms ----------------------------------------------

#DNA methylation perc. histogram separate replicates
#FIG. S1(B) ##THIS MAY NOT RUN LOCALLY, COMPUTATIONALLY INTENSIVE
uflong3 <- uflong2 %>% mutate(type = str_extract(name, "^[^ ]+"))

plot_S1B_A <- uflong3 %>% filter(c == "A") %>% na.omit() %>% 
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position = "identity", alpha=0.8, colour="black", linewidth=0.3) +
  scale_fill_manual(values = c("#cdccff","#cdccff","#cdccff","#cdccff","#cdccff","#cdccff","#cdccff","#cdccff"))+
  theme_classic() +
  xlab("Methylation (%)") +
  ylab("Count (x10^5)")+
  facet_wrap(~factor(type, levels=c("Spg", "P/D", "RS", "Sperm")), nrow = 1)+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-5))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7,
        legend.position = "none")

plot_S1B_X <- uflong3 %>% filter(c == "X") %>% na.omit() %>% 
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position = "identity", alpha=0.8, colour="black", linewidth=0.3) +
  scale_fill_manual(values = c("#3e68db","#3e68db","#3e68db","#3e68db","#3e68db","#3e68db","#3e68db","#3e68db"))+
  theme_classic() +
  xlab("Methylation (%)") +
  ylab("Count (x10^5)")+
  facet_wrap(~factor(type, levels=c("Spg", "P/D", "RS", "Sperm")), nrow = 1)+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-5, accuracy = 0.1))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7,
        legend.position = "none")

grid <- grid.arrange(plot_S1B_A, plot_S1B_X, nrow=2)
grid

ggsave(files = "figures/Fig.S1(B)_DNA_meth_perc_hist_wg_A_vs_X_unfiltered_separate_transparency.pdf", plot = grid, dpi=800, height=16, width=13)

#DNA methylation perc. histogram combined replicates
#FIG. 1(C) ##THIS MAY NOT RUN LOCALLY, COMPUTATIONALLY INTENSIVE
plot_1C_A <- uflong %>% filter(c == "A") %>% na.omit() %>% 
  ggplot(aes(x=value)) +
  geom_histogram(fill="#cdccff") +
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

plot_1C_X <- uflong %>% filter(c == "X") %>% na.omit() %>% 
  ggplot(aes(x=value)) +
  geom_histogram(fill="#3e68db") +
  theme_classic() +
  xlab("Methylation (%)") +
  ylab("Count (x10^5)")+
  facet_wrap(~factor(name, levels=c("Spg", "P/D", "RS", "Sperm")), nrow = 1)+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-5, accuracy = 0.1), n.breaks = 4)+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7)

grid <- grid.arrange(plot_1C_A, plot_1C_X, nrow=2)

grid

ggsave(files = "figures/Fig.1(C)_DNA_meth_perc_hist_wg_A_vs_X_unfiltered.pdf", plot = grid, dpi=800, height=8, width=13)

# Coverage histograms -----------------------------------------------------

#Creating long format dataframe for coverage plots
#Separate replicates
wgcov2 <- data
wgcov2 <- wgcov2 %>% dplyr::select(chr, start, end, 'covSpg 1', 'covP/D 1', 'covRS 1', 'covSperm 1', 'covSpg 2', 'covP/D 2', 'covRS 2', 'covSperm 2', c)
wgcovlong2 <- pivot_longer(wgcov2, 4:11)
wgcovlong2$name <- substring(wgcovlong2$name, 4)
wgcovlong3 <- wgcovlong2 %>% mutate(type = str_extract(name, "^[^ ]+"))

#FIG. S1(A) ##THIS MAY NOT RUN LOCALLY, COMPUTATIONALLY INTENSIVE
plot_S1A_A <- wgcovlong3 %>% filter(c == "A") %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(binwidth = 2, position = "identity", alpha=0.8, colour="black", linewidth=0.3) +
  scale_fill_manual(values = c("#cdccff","#cdccff","#cdccff","#cdccff","#cdccff","#cdccff","#cdccff","#cdccff"))+
  xlim(c(0,70)) +
  theme_classic() +
  xlab("Coverage")+
  ylab("Count (x10^5)")+
  facet_wrap(~factor(type, levels=c("Spg", "P/D", "RS", "Sperm")), nrow = 1)+
  geom_vline(xintercept = 15, col="#eb1b24", linetype = "dashed")+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-5, accuracy = 0.1))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7,
        legend.position = "none")

plot_S1A_X <- wgcovlong3 %>% filter(c == "X") %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(binwidth = 2, position = "identity", alpha=0.8, colour="black", linewidth=0.3) +
  scale_fill_manual(values = c("#3e68db","#3e68db","#3e68db","#3e68db","#3e68db","#3e68db","#3e68db","#3e68db"))+
  xlim(c(0,70)) +
  theme_classic() +
  xlab("Coverage")+
  ylab("Count (x10^5)")+
  facet_wrap(~factor(type, levels=c("Spg", "P/D", "RS", "Sperm")), nrow = 1)+
  geom_vline(xintercept = 15, col="#eb1b24", linetype = "dashed")+
  scale_y_continuous(labels = unit_format(unit="", scale = 1e-5))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7,
        legend.position = "none")

grid <- grid.arrange(plot_S1A_A, plot_S1A_X, nrow=2)

grid

ggsave(files = "figures/Fig.S1(A)_cov_hist_wg_A_vs_X_separate_unfiltered_transparency.pdf", plot = grid, dpi=800, height=16, width=13)

# CpG sampling relative to TSS --------------------------------------------

#No coverage filtering of intersected dataframe
plotdata <- x
plotdata <- plotdata %>% select(chr.x, start.x, end.x, Spg, 'P/D', RS, Sperm, posB, c)
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

ggsave(file = "figures/Fig.S1(E)_CpG_relative_TSS_unfiltered.pdf", dpi=800, height = 3.5, width = 6.5)

#FIG. S1(F)
methlong %>% na.omit() %>% ggplot(aes(posB, fill=c))+ geom_histogram()+
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

ggsave(file = "figures/Fig.S1(F)_CpG_relative_TSS_cov14.pdf", dpi=800, height = 3.5, width = 6.5)

# TSS plots ---------------------------------------------------------------

#Filtering for coverage above 14 in both biological replicates,
#filtering for at least two datapoints per relative position, and
#taking median value at each position relative to TSS
z <- x

meddata <- z %>% filter(z$'covSpg 1' > 14 & z$'covSpg 2' > 14) %>% select(posB, c, Spg)
meddata <- meddata %>% drop_na() %>% group_by(posB, c) %>% mutate(count=n()) %>% filter(count > 2) %>% summarise(medianpercmeth=median(Spg))
meddata$sample <- "Spg"
medianspg <- meddata

meddata <- z %>% filter(z$'covP/D 1' > 14 & z$'covP/D 2' > 14) %>% select(posB, c, 'P/D') %>% dplyr::rename('PD' = 'P/D')
meddata <- meddata %>% drop_na() %>% group_by(posB, c) %>% mutate(count=n()) %>% filter(count > 2) %>% summarise(medianpercmeth=median(PD))
meddata$sample <- "P/D"
medianpd <- meddata

meddata <- z %>% filter(z$'covRS 1' > 14 & z$'covRS 2' > 14) %>% select(posB, c, RS)
meddata <- meddata %>% drop_na() %>% group_by(posB, c) %>% mutate(count=n()) %>% filter(count > 2) %>% summarise(medianpercmeth=median(RS))
meddata$sample <- "RS"
medianrs <- meddata

meddata <- z %>% filter(z$'covSperm 1' > 14 & z$'covSperm 2' > 14) %>% select(posB, c, Sperm)
meddata <- meddata %>% drop_na() %>% group_by(posB, c) %>% mutate(count=n()) %>% filter(count > 2) %>% summarise(medianpercmeth=median(Sperm))
meddata$sample <- "Sperm"
mediansperm <- meddata

#Putting cell type dataframes together
allmedians <- rbind(medianspg, medianpd, medianrs, mediansperm)
allmedians <- allmedians %>% ungroup()
allmedians$sample_c <- paste(allmedians$sample, allmedians$c, sep = "_")

# TSS median plots

#FIG. 1(E)
allmedians %>% ggplot(aes(x=posB, y=medianpercmeth, color=c)) +
  geom_smooth(aes(fill=c), se=F, method="loess", span = 0.35) +
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

ggsave(file = "figures/TSS_median_plots_cov14_count2.pdf", dpi=800, height = 6, width = 13)

#FIG. S1(C)
allmedians %>% ggplot(aes(x=posB, y=medianpercmeth, color=c)) +
  geom_point(size=0.10)+
  labs(x="Position relative to TSS (kb)", y="Median methylation (%)") + 
  ylim(0,100)+
  coord_cartesian(x=c(-10000,10000))+
  theme_classic()+
  scale_color_manual(values = c("#333333", "#3E68DB"), name = "", labels = c("Autosomes", "X chromosome"))+
  scale_fill_manual(values = c("#333333", "#3E68DB"), name = "", labels = c("Autosomes", "X chromosome"))+
  # scale_color_manual(values = c("#3E68DB", "#333333"), name = "", labels = c("X chromosome", "Autosomes"))+
  # scale_fill_manual(values = c("#3E68DB", "#333333"), name = "", labels = c("X chromosome", "Autosomes"))+
  theme(legend.position="bottom")+
  facet_wrap(~factor(sample, levels = c("Spg", "P/D", "RS", "Sperm")), nrow=1)+
  scale_x_continuous(labels = unit_format(unit="", scale = 1e-3))+
  theme(legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        aspect.ratio = 0.7)

ggsave(file = "figures/TSS_median_plots_points_cov14_count2.pdf", dpi=800, height = 6, width = 13)

# Organising data into 500 CpG windows ------------------------------------

#Make separate dataframes for each cell type - unfiltered by coverage
#Renaming is to avoid issues with slashes
spg <- data %>% rename('Spg' = 'spg') %>% select(chr, start, end, spg, c)
pd <- data %>% rename('P/D' = 'pd') %>% select(chr, start, end, pd, c)
rs <- data %>% rename('RS' = 'rs') %>% select(chr, start, end, rs, c)
sperm <- data %>% rename('Sperm' = 'sperm') %>% select(chr, start, end, sperm, c)

#Define window side
chunk_size <- 500

spg_summary <- spg %>% na.omit() %>% mutate(group_index = (row_number() - 1) %/% chunk_size) %>% 
  group_by(chr, group_index) %>% 
  summarize(count = n(), average_value = mean(spg, na.rm = T)) %>% select(-count) %>% mutate(name = "Spg")

pd_summary <- pd %>% na.omit() %>% mutate(group_index = (row_number() - 1) %/% chunk_size) %>% 
  group_by(chr, group_index) %>% 
  summarize(count = n(), average_value = mean(pd, na.rm = T)) %>% select(-count) %>% mutate(name = "P/D")

rs_summary <- rs %>% na.omit() %>% mutate(group_index = (row_number() - 1) %/% chunk_size) %>% 
  group_by(chr, group_index) %>% 
  summarize(count = n(), average_value = mean(rs, na.rm = T)) %>% select(-count) %>% mutate(name = "RS")

sperm_summary <- sperm %>% na.omit() %>% mutate(group_index = (row_number() - 1) %/% chunk_size) %>% 
  group_by(chr, group_index) %>% 
  summarize(count = n(), average_value = mean(sperm, na.rm = T)) %>% select(-count) %>% mutate(name = "Sperm")

#Add together, add levels to name
windows_boxplot_data <- rbind(spg_summary, pd_summary, rs_summary, sperm_summary)
windows_boxplot_data$name <- factor(windows_boxplot_data$name, levels=c("Spg", "P/D", "RS", "Sperm"))

#Adding c (chromosome type) column
windows_boxplot_data$c <- ifelse(windows_boxplot_data$chr=="chrX", "X", "A")

windows_boxplot_data$c <- gsub('A', 'Autosomes', windows_boxplot_data$c)
windows_boxplot_data$c <- gsub('X', 'X chromosome', windows_boxplot_data$c)

#Calculating summary stats
windows_boxplot_stats <- windows_boxplot_data %>% group_by(c, name) %>% summarise(mean = mean(average_value), median = median(average_value), count = n())
windows_boxplot_stats

write.table(windows_boxplot_stats, file="data/Fig.1(D)_500_CpG_boxplots_wg_unfiltered_stats.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# 500 CpG window boxplots -------------------------------------------------

#FIG. 1(D)
windows_boxplot_data %>% ggplot(aes(y=average_value, x=c, fill=c))+
  geom_boxplot(notch = T, outlier.shape = NA)+
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

ggsave(file = "figures/Fig.1(D)_500_CpG_boxplots_wg_unfiltered.pdf", dpi=800, height = 5, width = 13)

#Mood's median test on this data

windows_boxplot_data$name_c <- paste(windows_boxplot_data$name, windows_boxplot_data$c, sep = "_")

moods <- pairwiseMedianTest(average_value~name_c, data=windows_boxplot_data, method = "holm")

moods_sig <- moods %>% filter(p.adjust < 0.05)

#write.table(moods, file="data/Moods_median_500_row_boxplots.txt", quote=F, sep="\t", row.names = F, col.names = T)

# Rsx locus plots ---------------------------------------------------------

#Filtering data (combined replicates) for coverage >9 (>14 is too stringent for single-locus plots)
singlelocus <- data

#Filtering for coverage of at least 10 in individual replicates - leaves an NA if cov not sufficient
singlelocus$Spg <- ifelse(singlelocus$'covSpg 1' > 9 & singlelocus$'covSpg 2' > 9, singlelocus$Spg, NA)
singlelocus$'P/D' <- ifelse(singlelocus$'covP/D 1' > 9 & singlelocus$'covP/D 2' > 9, singlelocus$'P/D', NA)
singlelocus$RS <- ifelse(singlelocus$'covRS 1' > 9 & singlelocus$'covRS 2' > 9, singlelocus$RS, NA)
singlelocus$Sperm <- ifelse(singlelocus$'covSperm 1' > 9 & singlelocus$'covSperm 2' > 9, singlelocus$Sperm, NA)

#Simplifying and making long format
singlelocus <- singlelocus %>% select(chr, start, end, Spg, 'P/D', RS, Sperm)
singlelocuslong <- pivot_longer(singlelocus, 4:7)

#Rsx start and end positions obtained by BLASTing Monodelphis domestica Rsx to this Tammar wallaby genome

Rsxplotdata <- singlelocuslong %>% na.omit() %>% filter(chr == "chrX" & start > 54720626-15000 & end < 54759097+15000)

#FIG. 1(F)
ggplot(data=Rsxplotdata, aes(x=start, y=value))+
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

ggsave(file = "figures/Rsx_meth_plots.pdf", dpi=800, height = 6, width = 13)

#FIG. S1(D)
ggplot(data=Rsxplotdata, aes(x=start, y=value))+
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

ggsave(file = "figures/Rsx_meth_plots_points.pdf", dpi=800, height = 6, width = 13)
