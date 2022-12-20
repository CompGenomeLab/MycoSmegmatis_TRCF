install.packages(
  c("ggplot2",
    "maditr",
    "dplyr",
    "reshape2",
    "tidyverse", 
    "ggpubr", 
    "matrixStats", 
    "heatmaply",
    "ggridges",
    "tidyr",
    "cowplot"))

library(ggplot2)
library(maditr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(matrixStats)
library(heatmaply)
library(rstatix)
library(ggridges)
library(dplyr)
library(tidyr)
library(cowplot)

setwd("./")

######### Input Files #########
project <- "myco"
lengthDistributionFile <- paste("results", project, "length.tsv", sep="/")
nucleotideContentFile <- paste("results", project, "nucleotide_content.tsv", sep="/")
rnaFile <- paste("results", project, "rna_readCounts.tsv", sep="/")
mappableTSNTS_file <- paste("results", project, "mappable_TSNTS.tsv", sep="/")
readCountsFile <- paste("results", project, "readCountsTSNTS.tsv", sep="/")
randomReadsFile <- paste("results", project, "random_readCountsTSNTS.tsv", sep="/")
###############################

######### Output Files #########
out_p_lenNucCor <- paste("results", "figures", "lenNucCor.pdf", sep="/") # Fig 3
out_p_TCRfig <- paste("results", "figures", "TCR.png", sep="/") # Fig 5
########################################


######### Length Distribuion #########
d <- read.table(lengthDistributionFile, sep='\t', header = F)
names(d) <- c('len', 'count', 'ind', 'product', 'rep', 'strain', 'time', 'title')
summary <- d %>% group_by(title) %>% summarize("total"=sum(count))
d <- merge(d,summary,by=c("title"))
d$freq <- 100*d$count/d$total
d$strain <- factor(d$strain, 
                   levels = c("WT","mfd","uvrD","uvrD_mfd"),
                   labels = c("WT","Delta~mfd","Delta~uvrD","Delta~uvrD~Delta~mfd"))
d$rep <- factor(d$rep, 
                levels = c(1,2),
                labels = c("Rep.~1","Rep.~2"))
p_lenDist <- ggplot(d, aes(x=len, y=freq)) +
  geom_col(aes(fill=factor(len))) +
  scale_fill_manual(values = c(rep("#344D67",12), "#ED7D3A", rep("#344D67",13))) +
  facet_grid(rows = vars(product), cols=vars(strain,rep), labeller = label_parsed) +
  scale_x_continuous(limits = c(-1, 25), breaks = seq(0, 25, by = 10)) +
  scale_y_continuous(limits = c(-1, 35), breaks = seq(0, 50, by = 10)) + 
  xlab("Read Length") +
  ylab("Frequency (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none",
        strip.text.y = element_text(color="white"))
########################################

######### Nucleotide Frequency #########
d <- read.table(nucleotideContentFile, sep='\t', header = F)
names(d) <- c('pos', 'nucleotide', 'count', 'percent', 'length', 'ind', 'product', 'rep', 'strain', 'time', 'title')
d$length <- factor(d$length, levels = rev(c(8,9,10,11,12,13)))
d$strain <- factor(d$strain, levels = c("WT","mfd","uvrD","uvrD_mfd"))
d$nucleotide <- factor(d$nucleotide, levels = c("A", "G", "C", "T"))
p_nucFreq <- ggplot(d, aes(x=pos, y=percent, fill=nucleotide)) +
  geom_col() +
  scale_fill_manual(values = rev(c("#344D67", "#6ECCAF", "#ADE792", "#F3ECB0")), name="") +
  facet_grid(rows = vars(length), cols=vars(strain,rep)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 25), labels = scales::percent) +
  xlab("Position") +
  ylab("Frequency (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position="top",
        legend.key.size = unit(0.2, 'cm'))
########################################

######### RNA Data Preparation #########
rna_dat <- read.table(rnaFile, sep='\t', header = F)
rna_titles <- c("chr", "start", "end", "name", "score", "strand", "count", "replicate", "time", "title")
names(rna_dat) <- rna_titles
rna_summary_d <- rna_dat %>% 
  group_by(replicate, time) %>% 
  summarise(totalRead = sum(count))

rna_dm <- merge(rna_dat, rna_summary_d, by= c("replicate", "time"), all.x = TRUE)
rna_dm$RPKM <- (rna_dm$count / (rna_dm$end - rna_dm$start + 1)) *1000000/rna_dm$totalRead
rna_md <- rna_dm %>% reshape2::dcast(name ~ title, value.var = "RPKM")
rna_md <- rna_md %>% mutate(quartile_A = ntile(`16h_1`, 4))
rna_md <- rna_md %>% mutate(quartile_B = ntile(`16h_2`, 4))
rna_md$consensusQ <- rna_md$quartile_A + rna_md$quartile_B
########################################


######### Mappable Sites  #########
map <- read.table(mappableTSNTS_file, sep='\t', header = F)
titles <- c("chr", "start", "end", "name", "score", "strand", "TS", "NTS")
names(map) <- titles
#ggplot(map) + geom_histogram(aes(x=TS+NTS)) + xlim(c(0,50))
minSiteCutoff = 10
highDamageSite = filter(map, TS>=minSiteCutoff & NTS>=minSiteCutoff)
genesWithMeaningfulMappableSites = highDamageSite$name
###################################


######### Main TS/NTS Analysis  #########

dat <- read.table(readCountsFile, sep='\t', header = F)
r <- read.table(randomReadsFile, sep='\t', header = F)

titles <- c("chr", "start", "end", "name", "score", "strand", "TS", "NTS", "induction", "product", "replicate", "strain", "time", "title")
names(dat) <- titles
titles_r <- c("chr", "start", "end", "name", "score", "strand", "TSr", "NTSr", "induction", "product", "replicate", "strain", "time", "title")
names(r) <- titles_r
dr <- subset(r, select=c("chr", "start", "end", "title", "TSr", "NTSr"))
d <- merge(dat,dr,by=c("chr", "start", "end", "title"))
d <- filter(d, induction=="True")

d$TSp1 <- d$TS + 1
d$NTSp1 <- d$NTS + 1
d$TSNTSa <- log2(d$TSp1/d$NTSp1)

d$TSrp1 <- d$TSr + 1
d$NTSrp1 <- d$NTSr + 1
d$TSNTSr <- log2(d$TSrp1/d$NTSrp1)

d$TSNTS <- d$TSNTSa - d$TSNTSr
d$total <- log2((d$TSp1 + d$NTSp1)/(d$TSrp1 + d$NTSrp1))
d$strain <- factor(d$strain, levels=c("WT", "mfd", "uvrD", "uvrD_mfd"))

summary_d <- d %>% 
  group_by(strain, replicate, induction) %>% 
  summarise(mean_val = mean(TSNTS), med_val = median(TSNTS), totalRepair = sum(TS) + sum(NTS))

dm <- merge(d, summary_d, by= c("strain", "replicate", "induction"), all.x = TRUE)
dm$RPKM <- (dm$TS + dm$NTS / (dm$end - dm$start + 1)) *1000000/dm$totalRepair
# dmf <- filter(dm, RPKM>8)
dmf <- filter(dm, TRUE)


# Merge with RNA data
reptr <- merge(d, rna_md, by= c("name"), all.x = TRUE)
reptr <- transform(reptr, minSite = pmin(TSr, NTSr, TS, NTS))
genebased <- reptr %>% reshape2::dcast(name ~ title, value.var = "minSite")
genebased$row_minimum = rowMins(as.matrix(genebased[,c(-1)]))

# ggplot(genebased) + geom_histogram(aes(x=row_minimum), binwidth = 1) + xlim(c(-1,100))
rowMin_whitelist = filter(genebased, row_minimum>=10)
whitelist <- filter(rowMin_whitelist, name %in% genesWithMeaningfulMappableSites)$name
q1reptrt <- filter(reptr, consensusQ==8, name %in% whitelist)

q1_summary <- q1reptrt %>% 
  group_by(strain, replicate, induction) %>% 
  summarise(mean_val = mean(TSNTS), med_val = median(TSNTS), totalRepair = sum(TS) + sum(NTS))

q1_fil <- filter(q1reptrt, TRUE)
q1_fil_d <- q1_fil %>% reshape2::dcast(name ~ title, value.var = "TSNTS")


#########################################


############ Correlation ############
results <- q1_fil %>% filter(induction=="True") %>% cohens_d(TSNTS ~ title, var.equal = TRUE)

reprt_d <- filter(reptr, induction=="True") %>% reshape2::dcast(name ~ strain + replicate, sep=" ", value.var = "total")
reprt_d <- reprt_d[,c(-1)]
my_cor <- cor(reprt_d)
my_cor_m <- reshape::melt(my_cor, value.name="Correlation")
my_cor_m$Correlation <- my_cor_m$value

heatmap_plot <- ggplot(data = my_cor_m, aes(x = X1, y = X2)) +
  geom_tile(aes(fill = Correlation)) +
  scale_fill_gradientn(colors = c("#FDE725FF", "#2B748EFF", "#440154FF")) +
  theme(
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#####################################

######### Length Distribuion, Nucleotide Frequency and Correlation #########
p_lenNucCor <- ggdraw() +
  draw_plot(p_lenDist, x = 0, y = 0.75, width = 1, height = .25) +
  draw_plot(p_nucFreq, x = 0, y = 0.25, width = 1, height = .50) +
  draw_plot(heatmap_plot, x = 0.25, y = 0, width = 0.5, height = .25) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0.20), y = c(1, 0.75, 0.25))
ggsave2("p_lenNucCor.pdf", p_lenNucCor, units="cm", width = 24, height = 28)
########################################



d_merged_rep.TS <- d %>% reshape2::dcast(name ~ strain, value.var = "TS", fun.aggregate = sum)
d_merged_rep.TS_m <- d_merged_rep.TS %>% reshape2::melt(value.name="TSmerged", variable.name="strain")
d_merged_rep.NTS <- d %>% reshape2::dcast(name ~ strain, value.var = "NTS", fun.aggregate = sum)
d_merged_rep.NTS_m <- d_merged_rep.NTS %>% reshape2::melt(value.name="NTSmerged", variable.name="strain")
d_merged_rep.TSr <- d %>% reshape2::dcast(name ~ strain, value.var = "TSr", fun.aggregate = sum)
d_merged_rep.TSr_m <- d_merged_rep.TSr %>% reshape2::melt(value.name="TSrmerged", variable.name="strain")
d_merged_rep.NTSr <- d %>% reshape2::dcast(name ~ strain, value.var = "NTSr", fun.aggregate = sum)
d_merged_rep.NTSr_m <- d_merged_rep.NTSr %>% reshape2::melt(value.name="NTSrmerged", variable.name="strain")

my_merge <- function(df1, df2){                       
  return (merge(df1, df2, by = c("name", "strain")))
}
a <- my_merge(d_merged_rep.TS_m, d_merged_rep.NTS_m)
b <- my_merge(d_merged_rep.TSr_m, d_merged_rep.NTSr_m)
d_merged_rep <- my_merge(a,b)

rep1 <- filter(d, replicate==1)

masterd <- subset(merge(rep1, d_merged_rep, by = c("name", "strain")), select = -c(replicate, induction, title, TS, NTS, TSp1, NTSp1, TSr, NTSr, TSNTS, total, TSrp1, NTSrp1, time, product, score, TSNTSa, TSNTSr))

colnames(masterd)[colnames(masterd) == "TSmerged"] = "TS"
colnames(masterd)[colnames(masterd) == "NTSmerged"] = "NTS"
colnames(masterd)[colnames(masterd) == "TSrmerged"] = "TSr"
colnames(masterd)[colnames(masterd) == "NTSrmerged"] = "NTSr"

masterd$TSp1 <- masterd$TS + 1
masterd$NTSp1 <- masterd$NTS + 1
masterd$TSNTSa <- log2(masterd$TSp1/masterd$NTSp1)

masterd$TSrp1 <- masterd$TSr + 1
masterd$NTSrp1 <- masterd$NTSr + 1
masterd$TSNTSr <- log2(masterd$TSrp1/masterd$NTSrp1)

masterd$TSNTS <- masterd$TSNTSa - masterd$TSNTSr
masterd$total <- log2((masterd$TSp1 + masterd$NTSp1)/(masterd$TSrp1 + masterd$NTSrp1))

master_reptr <- merge(masterd, rna_md, by= c("name"), all.x = TRUE)

master_reptr <- transform(master_reptr, minSite = pmin(TSr, NTSr, TS, NTS))
master_genebased <- master_reptr %>% reshape2::dcast(name ~ strain, value.var = "minSite")
master_genebased$row_minimum = rowMins(as.matrix(master_genebased[,c(-1)]))

# ggplot(master_genebased) + geom_histogram(aes(x=row_minimum), binwidth = 1) + xlim(c(-1,100))

master_rowMin_whitelist = filter(master_genebased, row_minimum>=10)
master_whitelist <- filter(master_rowMin_whitelist, name %in% genesWithMeaningfulMappableSites)$name

master_q1reptrt <- filter(master_reptr, consensusQ==8, name %in% master_whitelist)
master_q1_summary <- master_q1reptrt %>% 
  group_by(strain) %>% 
  summarise(mean_val = mean(TSNTS), med_val = median(TSNTS), totalRepair = sum(TS) + sum(NTS))

master_q1_fil <- filter(master_q1reptrt, TRUE)
master_q1_fil_d <- master_q1_fil %>% reshape2::dcast(name + start + end ~ strain, value.var = "TSNTS")
master_q1_fil_d$len <- master_q1_fil_d$end - master_q1_fil_d$start


#### DENSITY PLOTS #####

f <- filter(q1_fil, strain %in% c("WT", "mfd"))

rna_dm <- merge(rna_dat, rna_summary_d, by= c("replicate", "time"), all.x = TRUE)
rna_dm$RPKM <- (rna_dm$count / (rna_dm$end - rna_dm$start + 1)) *1000000/rna_dm$totalRead
rna_md <- rna_dm %>% reshape2::dcast(name ~ title, value.var = "RPKM")
rna_md <- rna_md %>% mutate(quantile_A = ntile(`16h_1`, 10))
rna_md <- rna_md %>% mutate(quantile_B = ntile(`16h_2`, 10))
rna_md_f <- filter(rna_md, quantile_A==quantile_B)
rna_md_f$quantile <- rna_md_f$quantile_A

reptr <- merge(d, rna_md_f, by= c("name"), all.x = TRUE)

fd <- reptr %>% reshape2::dcast(name + quantile ~ strain + replicate, value.var = "TSNTS")

fd$base_var_1 <- fd$WT_1
fd$test_var_1 <- fd$mfd_1
fd$base_var_2 <- fd$WT_2
fd$test_var_2 <- fd$mfd_2

fd$rep1diff_1 <- fd$test_var_1 - fd$base_var_1
fd$rep2diff_1 <- fd$test_var_2 - fd$base_var_2
fd$cons_1 <- fd$rep1diff_1/abs(fd$rep1diff_1) * fd$rep2diff_1/abs(fd$rep2diff_1)

fd$base_var_1 <- fd$WT_1
fd$test_var_1 <- fd$uvrD_1
fd$base_var_2 <- fd$WT_2
fd$test_var_2 <- fd$uvrD_2

fd$rep1diff_2 <- fd$test_var_1 - fd$base_var_1
fd$rep2diff_2 <- fd$test_var_2 - fd$base_var_2
fd$cons_2 <- fd$rep1diff_2/abs(fd$rep1diff_2) * fd$rep2diff_2/abs(fd$rep2diff_2)

fd$base_var_1 <- fd$WT_1
fd$test_var_1 <- fd$uvrD_mfd_1
fd$base_var_2 <- fd$WT_2
fd$test_var_2 <- fd$uvrD_mfd_2

fd$rep1diff_3 <- fd$test_var_1 - fd$base_var_1
fd$rep2diff_3 <- fd$test_var_2 - fd$base_var_2
fd$cons_3 <- fd$rep1diff_3/abs(fd$rep1diff_3) * fd$rep2diff_3/abs(fd$rep2diff_3)

fd$base_var_1 <- fd$uvrD_1
fd$test_var_1 <- fd$mfd_1
fd$base_var_2 <- fd$uvrD_2
fd$test_var_2 <- fd$mfd_2

fd$rep1diff_4 <- fd$test_var_1 - fd$base_var_1
fd$rep2diff_4 <- fd$test_var_2 - fd$base_var_2
fd$cons_4 <- fd$rep1diff_4/abs(fd$rep1diff_4) * fd$rep2diff_4/abs(fd$rep2diff_4)

fd$base_var_1 <- fd$mfd_1
fd$test_var_1 <- fd$uvrD_mfd_1
fd$base_var_2 <- fd$mfd_2
fd$test_var_2 <- fd$uvrD_mfd_2

fd$rep1diff_5 <- fd$test_var_1 - fd$base_var_1
fd$rep2diff_5 <- fd$test_var_2 - fd$base_var_2
fd$cons_5 <- fd$rep1diff_5/abs(fd$rep1diff_5) * fd$rep2diff_5/abs(fd$rep2diff_5)

fd$base_var_1 <- fd$uvrD_1
fd$test_var_1 <- fd$uvrD_mfd_1
fd$base_var_2 <- fd$uvrD_2
fd$test_var_2 <- fd$uvrD_mfd_2

fd$rep1diff_6 <- fd$test_var_1 - fd$base_var_1
fd$rep2diff_6 <- fd$test_var_2 - fd$base_var_2
fd$cons_6 <- fd$rep1diff_6/abs(fd$rep1diff_6) * fd$rep2diff_6/abs(fd$rep2diff_6)

# fdf <- filter(fd, fd$cons_4>0, fd$quantile >0, fd$quantile <= 10)
fdf <- filter(fd, fd$quantile >0, fd$quantile <= 10)

theme <-   theme_classic() +
  theme(legend.position="none", 
        axis.text.y=element_blank(), 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

p1 <- ggplot(
  fdf, 
  aes(x = rep1diff_1, y = reorder(quantile, desc(quantile)), fill = quantile)
) +
  geom_density_ridges_gradient(scale = 2, size = 0.2, rel_min_height = 0.01) +
  xlim(-4,4) +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  scale_fill_viridis_c(name = "TCR Difference") +
  ylab("Transcription Quantile") +
  xlab("Δmfd TCR relative to WT\n[log2(Δmfd/WT)]") +
  theme_classic() +
  theme(legend.position="none", axis.line.y = element_blank())

p2 <- ggplot(
  fdf, 
  aes(x = rep1diff_2, y = reorder(quantile, desc(quantile)), fill = quantile)
) +
  geom_density_ridges_gradient(scale = 2, size = 0.2, rel_min_height = 0.01) +
  xlim(-4,4) +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  scale_fill_viridis_c(name = "TCR Difference") +
  ylab("") +
  xlab("ΔuvrD TCR relative to WT\n[log2(ΔuvrD/WT)]") +
  theme


p3 <- ggplot(
  fdf, 
  aes(x = rep1diff_3, y = reorder(quantile, desc(quantile)), fill = quantile)
) +
  geom_density_ridges_gradient(scale = 2, size = 0.2, rel_min_height = 0.01) +
  xlim(-4,4) +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  scale_fill_viridis_c(name = "TCR Difference") +
  ylab("") +
  xlab("ΔuvrD_Δmfd TCR relative to WT\n[log2(ΔuvrD_Δmfd/WT)]") +
  theme


p4 <- ggplot(
  fdf, 
  aes(x = rep1diff_4, y = reorder(quantile, desc(quantile)), fill = quantile)
) +
  geom_density_ridges_gradient(scale = 2, size = 0.2, rel_min_height = 0.01) +
  xlim(-4,4) +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  scale_fill_viridis_c(name = "TCR Difference") +
  ylab("") +
  xlab("Δmfd TCR relative to ΔuvrD\n[log2(Δmfd/ΔuvrD]") +
  theme


p5 <- ggplot(
  fdf, 
  aes(x = rep1diff_5, y = reorder(quantile, desc(quantile)), fill = quantile)
) +
  geom_density_ridges_gradient(scale = 2, size = 0.2, rel_min_height = 0.01) +
  xlim(-4,4) +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  scale_fill_viridis_c(name = "TCR Difference") +
  ylab("") +
  xlab("ΔuvrD_Δmfd TCR relative to Δmfd\n[log2(ΔuvrD_Δmfd/Δmfd]") +
  theme


p6 <- ggplot(
  fdf, 
  aes(x = rep1diff_6, y = reorder(quantile, desc(quantile)), fill = quantile)
) +
  geom_density_ridges_gradient(scale = 2, size = 0.2, rel_min_height = 0.01) +
  xlim(-4,4) +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  scale_fill_viridis_c(name = "TCR Difference") +
  ylab("") +
  xlab("ΔuvrD_Δmfd TCR relative to ΔuvrD\n[log2(ΔuvrD_Δmfd/ΔuvrD]") +
  theme


theme0 <-   theme_bw() +
  theme(legend.position="none")

master_ggpaired_ <- function(cond1_, cond2_, color_pair, p, star){
  
  test_method = "t.test"
  # test_method = "wilcox.test"
  st <- stat_compare_means(paired=TRUE, method=test_method, label = "p.format", label.x.npc = "left", label.y.npc = "bottom")
  point_size = 0
  width_ = 0.25
  line_size = 0.05
  pal <- get_palette(c("#00AFBB", "#FC4E07"), 4)
  return (ggpaired(master_q1_fil_d, cond1 = cond1_ , cond2 = cond2_,
                   fill = "condition", palette = color_pair, 
                   line.size = line_size, 
                   point.size = point_size,
                   width = width_,
                   ylab="", xlab=""
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    ylim(-4,6) + theme0 +
    annotate(
      "text", label = paste('p <',format(p$p.value, scientific = TRUE, trim = TRUE, digits=3)),
      x = 1.5, y = -4, size = 4, colour = "black"
    ) +
    annotate("text", label = star, x = 1.5, y = -3.5, size = 4, colour = "black")
    )
}
pal <- get_palette(c("#00AFBB", "#FC4E07"), 4)

starfun <- function(p) {
  if (p<1e-6) {
    star = "***"
  } else if (p<1e-3) {
    star = "**"
  } else if (p<0.05) {
    star = "*"
  } else {
    star = "ns"
  }
  return(star)
}

cond1_ = "WT"
cond2_ = "mfd"
t <- t.test(master_q1_fil_d$WT, master_q1_fil_d$mfd, paired = TRUE, alternative = "two.sided")
p <- t[3]
color_pair = get_palette(c(pal[1], pal[2]), 2)
master_twoStrain1 <- master_ggpaired_(cond1_, cond2_, color_pair, p, starfun(p)) + 
  ylab("Relative TCR (TS/NTS)") +
  scale_x_discrete(labels=c('WT', 'Δmfd'))
# ggsave("master_WT-mfd.pdf", plot=master_twoStrain1)


theme <-   theme_bw() +
  theme(legend.position="none", 
        axis.text.y=element_blank(), 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

cond1_ = "WT"
cond2_ = "uvrD"
t <- t.test(master_q1_fil_d$WT, master_q1_fil_d$uvrD, paired = TRUE, alternative = "two.sided")
p <- t[3]
color_pair = get_palette(c(pal[1], pal[3]), 2)
master_twoStrain2 <- master_ggpaired_(cond1_, cond2_, color_pair, p, starfun(p)) + theme + 
  theme(axis.text.y=element_text(color="white")) +
  scale_x_discrete(labels=c('WT', 'ΔuvrD'))

cond1_ = "WT"
cond2_ = "uvrD_mfd"
t <- t.test(master_q1_fil_d$WT, master_q1_fil_d$uvrD_mfd, paired = TRUE, alternative = "two.sided")
p <- t[3]
color_pair = get_palette(c(pal[1], pal[4]), 2)
master_twoStrain3 <- master_ggpaired_(cond1_, cond2_, color_pair, p, starfun(p)) + theme +
  scale_x_discrete(labels=c('WT', 'ΔuvrD Δmfd'))

cond1_ = "uvrD"
cond2_ = "mfd"
t <- t.test(master_q1_fil_d$uvrD, master_q1_fil_d$mfd, paired = TRUE, alternative = "two.sided")
p <- t[3]
color_pair = get_palette(c(pal[3], pal[2]), 2)
master_twoStrain4 <- master_ggpaired_(cond1_, cond2_, color_pair, p, starfun(p)) + theme +
  scale_x_discrete(labels=c('ΔuvrD', 'Δmfd'))

cond1_ = "mfd"
cond2_ = "uvrD_mfd"
t <- t.test(master_q1_fil_d$mfd, master_q1_fil_d$uvrD_mfd, paired = TRUE, alternative = "two.sided")
p <- t[3]
color_pair = get_palette(c(pal[2], pal[4]), 2)
master_twoStrain5 <- master_ggpaired_(cond1_, cond2_, color_pair, p, starfun(p)) + theme +
  scale_x_discrete(labels=c('Δmfd', 'ΔuvrD Δmfd'))

cond1_ = "uvrD"
cond2_ = "uvrD_mfd"
t <- t.test(master_q1_fil_d$uvrD, master_q1_fil_d$uvrD_mfd, paired = TRUE, alternative = "two.sided")
p <- t[3]
color_pair = get_palette(c(pal[3], pal[4]), 2)
master_twoStrain6 <- master_ggpaired_(cond1_, cond2_, color_pair, p, starfun(p)) + theme +
  scale_x_discrete(labels=c('ΔuvrD', 'ΔuvrD Δmfd'))
# master_twoStrain6

TCRfig <- ggarrange(master_twoStrain1, 
                  master_twoStrain2, 
                  master_twoStrain3, 
                  master_twoStrain4,
                  master_twoStrain5,
                  master_twoStrain6,
                  p1, p2, p3, p4, p5, p6, 
                  labels = c("A", "", "", "", "", "",
                             "B", "", "", "", "", ""),
                  ncol = 6, nrow = 2)

ggexport(TCRfig, filename=out_p_TCRfig, width=5000, height=3000, res=300)
