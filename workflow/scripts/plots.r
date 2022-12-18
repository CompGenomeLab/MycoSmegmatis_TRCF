install.packages(c("ggplot2","maditr","reshape2","tidyverse", "ggpubr"))

library(ggplot2)
library(maditr)
library(tidyverse)
library(reshape2)
library(ggpubr)

setwd("/home/oadebali/smeg_umi")

dat <- read.table("results/myco/readCountsTSNTS.tsv", sep='\t', header = F)
r <- read.table("results/myco/random/readCountsTSNTS.tsv", sep='\t', header = F)

titles <- c("chr", "start", "end", "name", "score", "strand", "TS", "NTS", "induction", "product", "replicate", "strain", "time", "title")
names(dat) <- titles
titles_r <- c("chr", "start", "end", "name", "score", "strand", "TSr", "NTSr", "induction", "product", "replicate", "strain", "time", "title")
names(r) <- titles_r

dr <- subset(r, select=c("chr", "start", "end", "title", "TSr", "NTSr"))

d <- merge(dat,dr,by=c("chr", "start", "end", "title"))



d$TSp1 <- d$TS + 1
d$NTSp1 <- d$NTS + 1
d$TSNTSa <- log2(d$TSp1/d$NTSp1)

d$TSrp1 <- d$TSr + 1
d$NTSrp1 <- d$NTSr + 1
d$TSNTSr <- log2(d$TSrp1/d$NTSrp1)

d$TSNTS <- d$TSNTSa - d$TSNTSr
d$total <- log2((d$TSp1 + d$NTSp1)/(d$TSrp1 + d$NTSrp1))


# d$strain <- recode(d$strain, mfd = 'Mfd', 
#                         UvrD_Mfd = 'UvrD Mfd')

# d$strain <- factor(ded$strain, levels=c("WT", "UvrD", "Mfd", "UvrD Mfd" ))
d$strain <- factor(d$strain, levels=c("WT", "mfd", "uvrD", "uvrD_mfd"))

summary_d <- d %>% 
  group_by(strain, replicate, induction) %>% 
  summarise(mean_val = mean(TSNTS), med_val = median(TSNTS), totalRepair = sum(TS) + sum(NTS))


dm <- merge(d, summary_d, by= c("strain", "replicate", "induction"), all.x = TRUE)
dm$RPKM <- (dm$TS + dm$NTS / (dm$end - dm$start + 1)) *1000000/dm$totalRepair
# dmf <- filter(dm, RPKM>8)
dmf <- filter(dm, TRUE)


# RNA part
rna_dat <- read.table("results/rna/readCounts.tsv", sep='\t', header = F)
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
# rna_topgenelist = filter(rna_md, consensusQ==8)
# topGenes = rna_topgenelist$name
####

reptr <- merge(d, rna_md, by= c("name"), all.x = TRUE)
q1reptrt <- filter(reptr, consensusQ==8)
q1_summary <- q1reptrt %>% 
  group_by(strain, replicate, induction) %>% 
  summarise(mean_val = mean(TSNTS), med_val = median(TSNTS), totalRepair = sum(TS) + sum(NTS))


q1_fil <- filter(q1reptrt, title=="WTind_5m_R1" | title=="MfdInd_5m_R1")
q1_fil_d <- q1_fil %>% reshape2::dcast(name ~ title, value.var = "TSNTS")
ggpaired(q1_fil_d, cond1 = "WTind_5m_R1" , cond2 = "MfdInd_5m_R1",
    fill = "condition", palette = "jco")



p1 <- ggplot(d, aes(TSNTS)) + 
  geom_histogram() + 
  # geom_vline(xintercept=0, color="red") +
  geom_vline(data = summary_d, aes(xintercept=mean_val), color="blue") +
  geom_vline(data = summary_d, aes(xintercept=med_val), color="green") +
  facet_grid(~strain~replicate~induction) +
  ylab("Gene count") +
  xlab("log2-transformed TS/NTS")
ggsave("results/plots/p1_histogram.pdf", plot=p1)


p1q <- ggplot(q1reptrt, aes(TSNTS)) + 
  geom_histogram() + 
  # geom_vline(xintercept=0, color="red") +
  geom_vline(data = q1_summary, aes(xintercept=mean_val), color="blue") +
  geom_vline(data = q1_summary, aes(xintercept=med_val), color="green") +
  facet_grid(~strain~replicate~induction) +
  ylab("Gene count") +
  xlab("log2-transformed TS/NTS")
ggsave("results/plots/p1q_histogram.pdf", plot=p1q)

md <- d %>% reshape2::dcast(name ~ strain + replicate, value.var = "TSNTS")
md2 <- d %>% reshape2::dcast(name ~ strain + replicate, value.var = "total")


p2 <- ggplot(md, aes(x=WT_1, y=mfd_1, color=total)) + geom_point() + 
  geom_smooth(method="lm", color="purple") +
  geom_abline(slope=1, intercept=0, color="blue") +
  ggtitle("log2-transformed TS/NTS values")
ggsave("results/plots/p2_scatter.pdf", plot=p2)

p3 <- ggplot(dm, aes(log2(RPKM))) + 
  geom_histogram() + 
  # geom_vline(xintercept=0, color="red") +
  facet_grid(~strain~replicate) +
  ylab("Gene count") +
  xlab("log2-transformed RPKM")
ggsave("results/plots/p3_RPKMhist.pdf", plot=p3)


p4 <- ggplot(dmf) + geom_point(aes(x=log2(RPKM), y=TSNTS, fill=total, text=name)) +
  geom_smooth(aes(x=log2(RPKM), y=TSNTS)) +
  facet_grid(~strain~replicate) +
  ylab("log2-transformed TS/NTS") +
  xlab("log2-transformed RPKM")
ggsave("results/plots/p4_scatterRPKMvsTSNTS.pdf", plot=p4)


library(plotly)
pp <- ggplotly(p4)
htmlwidgets::saveWidget(pp, "index.html")





gd <- merge(d, summary_d, by= c("strain", "replicate"), all.x = TRUE)

gd$TS_RPKM <- (gd$TS/ (gd$end - gd$start + 1)) *1000000/gd$totalRepair
gd$NTS_RPKM <- (gd$NTS/ (gd$end - gd$start + 1)) *1000000/gd$totalRepair

genes = list("MSMEI_5474", "MSMEI_5475", "MSMEI_5476", "MSMEI_5477", "MSMEI_5478",
             "MSMEI_1328", #rpoB
             "MSMEI_0963" #pnp
)

g <- filter(gd, name %in% genes)
mg <- reshape2::melt(g, id.vars = c("name", "strain", "replicate"), measure.vars = c("TS_RPKM", "NTS_RPKM"))

mgn <- mg %>%
  mutate_all(~case_when(name == "MSMEI_5474" ~ "mfd", 
                        name == "MSMEI_5475" ~ "MSMEI_5475", 
                        name == "MSMEI_5476" ~ "MSMEI_5476", 
                        name == "MSMEI_5477" ~ "glmU", 
                        name == "MSMEI_5478" ~ "prs", 
                        name == "MSMEI_1328" ~ "rpoB", 
                        name == "MSMEI_0963" ~ "pnp"))

mg$gene <- mgn$name 
  
p5 <- ggplot(mg) + geom_col(aes(x=strain, y=value*1000, fill=variable),position="dodge") + 
  scale_fill_grey() +
  scale_y_continuous(trans='log2') +
  facet_grid(~gene~replicate) +
  ylab("RPKB") + 
  xlab("Strain") +
  theme_bw()
p5
ggsave("results/plots/p5_bars.pdf", plot=p5)

dmf<- filter(dm, name %in% genes)


dmfn <- dmf %>%
  mutate_all(~case_when(name == "MSMEI_5474" ~ "mfd", 
                        name == "MSMEI_5475" ~ "MSMEI_5475", 
                        name == "MSMEI_5476" ~ "MSMEI_5476", 
                        name == "MSMEI_5477" ~ "glmU", 
                        name == "MSMEI_5478" ~ "prs", 
                        name == "MSMEI_1328" ~ "rpoB", 
                        name == "MSMEI_0963" ~ "pnp"))

dmf$gene <- dmfn$name
  
p6 <- ggplot(dmf) + geom_col(aes(x=strain, y=RPKM),position="dodge") + 
  scale_fill_grey() +
  # scale_y_continuous(trans='log2') +
  facet_grid(~gene~replicate) +
  ylab("RPKM") + 
  xlab("Strain") +
  theme_bw()

p6


# 
# ggplot(md, aes(x=WT_2, Mfd_2)) + geom_point() + 
#   geom_smooth(method="lm", color="purple") +
#   geom_abline(slope=1, intercept=0, color="blue")
# 
# ggplot(md, aes(x=WT_1, UvrD_1)) + geom_point() + 
#   geom_smooth(method="lm", color="purple") +
#   geom_abline(slope=1, intercept=0, color="blue")
# 
# ggplot(md, aes(x=WT_2, UvrD_2)) + geom_point() + 
#   geom_smooth(method="lm", color="purple") +
#   geom_abline(slope=1, intercept=0, color="blue")
# 
# ggplot(md, aes(x=WT_1, UvrD_Mfd_1)) + geom_point() + 
#   geom_smooth(method="lm", color="purple") +
#   geom_abline(slope=1, intercept=0, color="blue")
# 
# ggplot(md, aes(x=WT_2, UvrD_Mfd_2)) + geom_point() + 
#   geom_smooth(method="lm", color="purple") +
#   geom_abline(slope=1, intercept=0, color="blue")
# 
# ggplot(md, aes(x=mfd_1, UvrD_1)) + geom_point() + 
#   geom_smooth(method="lm", color="purple") +
#   geom_abline(slope=1, intercept=0, color="blue")
# 
# ggplot(md, aes(x=mfd_2, UvrD_2)) + geom_point() + 
#   geom_smooth(method="lm", color="purple") +
#   geom_abline(slope=1, intercept=0, color="blue")
# 
# ggplot(md, aes(x=mfd_1-UvrD_Mfd_1)) + geom_histogram()
