
options(scipen=999)

library(data.table)
library(tidyverse)
library(ggpubr)
library(ggforce)
library(scales)
library(patchwork)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

load("data/PCAWG_donorinfo/cohorts_palette.RData")
col_ref <- setNames(cohorts[1:37, "colcode"], cohorts[1:37, "cohort"])
names(col_ref) <- gsub(pattern = "AdenoCa", replacement = "AdenoCA", names(col_ref))
rm(cohorts)


snvcount <- fread("results/count_snv_per_donor_noncodingsnvs.tsv")
g <- snvcount[,.(Donor_ID, Tumor_Type, Color_code, total_snv_obs=as.numeric(snv_count))]

#sigact <- fread("data/PCAWG_signatures/SA_COMPOSITE_SNV.activity_relative.WHITELIST_SET.20200226.txt")
#g <- sigact[,.(Donor_ID=icgc_donor_id, Tumor_Type, dcc_project_code, total_snv_obs=as.numeric(total_snv_obs))]
#g[, Tumor_Type := gsub("_", "-", Tumor_Type)]

g[, `:=`(median_snv = median(total_snv_obs),
         mean_snv = mean(total_snv_obs),
         n_genomes_per_type = .N), by = Tumor_Type]


g[, sum(total_snv_obs)]
reflength <- sum(k11c_dt$tot_hg19)

g[, sum(total_snv_obs)/2583/reflength*10^6] ==
  mean(g$total_snv_obs)/reflength*10^6

# sina
mean_snv_rate <- g[, mean(total_snv_obs)/reflength*10^6]
p_tmb_sina <- g %>% 
  ggplot(aes(x = reorder(Tumor_Type, rank(mean_snv)), 
             y = total_snv_obs/reflength*10^6, fill = Color_code)) +
  geom_hline(yintercept=mean_snv_rate, linetype = "dashed", alpha = 0.5) +
  annotate(geom="text", x=-Inf, y=10, vjust=-1, hjust=0, label=paste("mean rate =\n",format(round(mean_snv_rate,2),nsmall=2),"SNV/patient/Mb")) +
  geom_sina(scale = FALSE, shape = 21, col = "black") +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  scale_y_log10("SNV/Mb",labels = function(x) prettyNum(x, big.mark = ","), breaks = c(0.1, 10, 10^3)) + #round(mean_snv_rate, 0)
  scale_fill_identity(guide = "none") +
  #scale_fill_manual(values = col_ref, guide = "none") +
  coord_flip(clip="off") +
  labs(x = "") +
  theme_pubr() +
  theme(strip.background = element_blank(), strip.text = element_blank())

gsum <- unique(g[,.(Tumor_Type, Color_code, mean_snv, median_snv, n_genomes_per_type)], 
               by = "Tumor_Type")[order(rank(mean_snv))]

p_ngenomes_bar <- gsum %>%
  ggplot(aes(x=reorder(Tumor_Type, rank(mean_snv)), y=n_genomes_per_type, fill = Color_code)) +
  geom_col(col = "black") +
  geom_text(aes(label = n_genomes_per_type), hjust = -.2) +
  scale_x_discrete(labels = gsum[, format(round(mean_snv/reflength*10^6, 2), nsmall = 2)]) +
  scale_y_continuous(breaks = c(0, 150, 300)) +
  scale_fill_identity(guide = "none") +
  #scale_fill_manual(values = col_ref, guide = "none") +
  labs(y="Genomes", x ="mean SNV/patient/Mb") +
  coord_flip(clip="off") +
  theme_pubr() +
  theme(axis.ticks.y = element_blank(), panel.grid.minor = element_blank())

p_align <- p_tmb_sina + p_ngenomes_bar + plot_spacer() +
  plot_layout(nrow = 1, widths = c(30, 15,1))

p_align
ggsave(plot=p_align, device="pdf", width=6, height=12, useDingbats = FALSE,
       filename="plots/dataset_tumortype_intro_v4.pdf")
