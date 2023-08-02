# author:     Gustav Poulsgaard
# date created: 30-06-2021
# purpose:    analyse VAF distribution of SNVs inside and outside
#               select alpha satellite repeats.
# details:    compute the within sample mean difference in VAF
# include:    samples with >20 SNVs in alpha satellite repeats (<=311 genomes)

# 
options(scipen=999)
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

# 
library(data.table)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(scales)

# load data
dt <- fread("data/signature_snv_filtered_annotated4.tsv")
setnames(dt, c("Chromosome", "End_position"), c("chr", "end"))
dt[, start:=end] ; setkey(dt, chr, start, end)

rmsk <- fread("../reference_work/encodetracks/rmsk_RepeatMasker_hg19.bed")
rmsk <- rmsk[genoName%in%paste0("chr", 1:22)]
setnames(rmsk, old=c("genoName","genoStart","genoEnd"),
         new=c("chr","start","end"))
setkey(rmsk, chr, start, end)

activity_dt <- fread("data/PCAWG_signatures/SA_COMPOSITE_SNV.activity_relative.WHITELIST_SET.20200226.txt")

# overlap rmsk annotations
t1 <- Sys.time()
dt1 <- foverlaps(dt, rmsk) # 2 min
t2 <- Sys.time() ; difftime(t2,t1)

# define region of interest
dt1[, region := ifelse(repName=="ALR/Alpha","alphaSatellite", "genomewide")]
summary(factor(dt1$region)) # 
dt1[is.na(region), region:="genomewide"]
# filter (include)
count_alr_donor <- dt1[,.(n=sum(region=="alphaSatellite")), by=Donor_ID]

# dt2 <- dt1[Donor_ID %in% count_alr_donor[n>=20, Donor_ID], 
#     .(VAF=mean(VAF)), by=c("Donor_ID","region")][order(Donor_ID)]
# 
# dt2[, VAF_diff:=VAF[region=="genomewide"]-VAF[region=="alphaSatellite"], by = Donor_ID]
# 
# dt2_l <- dt2 %>% spread(key="region",value="VAF") %>% 
#   gather(key="region", value="vaf", 2:4) %>% as.data.table()
# 
# dt2_l[, mean(vaf), by = region]
# p <- ggplot(dt2_l, aes(x = vaf, fill=region)) + 
#   geom_histogram(bins=30, alpha=0.5) +
#   geom_vline(data=dt2_l[, .(vaf=mean(vaf)), by = region], aes(xintercept=vaf))+
#   facet_wrap(region=="VAF_diff"~., scales = "free", ncol = 1) +
#   labs(title = "mean VAF per genome") +
#   theme_pubr(legend = "right") +
#   theme(strip.text = element_blank())
# 
# ggsave(plot=p, device = "pdf", filename="plots/vaf_distribution/vaf_alphaSatellite_vs_genomewide.pdf")
# 




# p2 <- dt1[Donor_ID%in%activity_dt[BI_COMPOSITE_SNV_SBS17b_P>=.05, icgc_donor_id] &
#       Donor_ID%in%count_alr_donor[n>=20, Donor_ID], 
#     .(VAF=mean(VAF), .N), by=c("Donor_ID","region")][order(Donor_ID)] %>%
#   ggplot(aes(x=factor(region,c("genomewide","alphaSatellite")), y=VAF, group=Donor_ID))+
#   geom_boxplot(aes(group=region)) +
#   geom_point(alpha=0.1)+
#   geom_line(alpha=0.1)+
#   coord_cartesian(ylim = c(0.3, .45)) +
#   labs(title = "Genomes w. S17b and >=20 alphaSatellite SNVs")
# ggsave(plot=p2, device = "pdf", filename="plots/vaf_distribution/test2.pdf")


# dt3 <- dt1[Donor_ID%in%activity_dt[BI_COMPOSITE_SNV_SBS17b_P>=.05, icgc_donor_id] &
#       Donor_ID%in%count_alr_donor[n>=20, Donor_ID] & 
#       BI_COMPOSITE_SNV_SBS17b_P>=.5, .(VAF=mean(VAF), .N), by=c("Donor_ID","region")][order(Donor_ID)]

dt3 <- dt1[Donor_ID%in%activity_dt[BI_COMPOSITE_SNV_SBS17b_P>=.05, icgc_donor_id] &
      BI_COMPOSITE_SNV_SBS17b_P>=.5, .(VAF=mean(VAF), .N), 
    by=c("Donor_ID","region")][N>=20][order(Donor_ID)]

plot_vaf <- function(signature, save_plot=TRUE){
  x <- dt1[Donor_ID%in%activity_dt[get(signature)>=.05, icgc_donor_id] &
        get(signature)>=.5, .(VAF=mean(VAF), .N),
        by=c("Donor_ID","region")][N>=20][order(Donor_ID)]
  y <- x[Donor_ID %in% x[, .(n_region=.N), by=Donor_ID][n_region==2, Donor_ID]]
  
  t.test.var <- t.test(x=y[region=="genomewide", VAF], y=y[region=="alphaSatellite", VAF], paired=TRUE)
  t.test.text <- paste0(t.test.var$method, "\nnumber of samples = ",y[,uniqueN(Donor_ID)],
                        "\nmean of the differences = ", round(t.test.var$estimate,4),
                        "\np-value < ", scientific(t.test.var$p.value, 2))
  y_stat <- y[, .(VAF=median(VAF), n_snv=sum(N), n_genome=.N), by=region]
  p_boxes <- y %>%
    ggplot(aes(x=factor(region,c("genomewide","alphaSatellite")), y=VAF, group=Donor_ID))+
    geom_boxplot(aes(group=region)) +  geom_point(alpha=0.1)+  geom_line(alpha=0.1)+
    annotate(geom="segment",x=1,xend=2,
             y=max(y$VAF)*1.02, yend=max(y$VAF)*1.02, color="black")+
    annotate(geom="text", x=1.5, y=max(y$VAF)*1.02, label=t.test.text, vjust=-.1)+
    scale_y_continuous(breaks=c(seq(0.3,0.45,length.out=4), y_stat$VAF),
                       labels=function(x) round(x,2)) +
    coord_cartesian(ylim = c(0.3, .45)) +
    labs(title="Sample mean VAF per genomic region", x="region",
         subtitle = paste0("Genomes w. ",str_sub(signature,20,-3),
                           " and >=20 alphaSatellite SNVs and pr[", 
                           str_sub(signature,20,-3),"|SNVs]>=50%")) +
    theme_pubr()
  y_diff <- y[,.(VAF_diff=VAF[region=="genomewide"]-VAF[region=="alphaSatellite"]), by = Donor_ID]
  y_diffstat <- y_diff[, summary(VAF_diff)]
  
  p_histo <- y_diff %>%
    ggplot(aes(x=VAF_diff)) + geom_histogram(col="black",fill="grey", bins=30) +
    geom_vline(xintercept=0, color="black", linetype="solid")+
    geom_vline(xintercept = y_diffstat["Mean"], linetype="dashed") +
    annotate(geom="text", x=y_diffstat["Mean"], y=Inf, label=round(y_diffstat["Mean"],4), hjust=-.1) +
    coord_cartesian(xlim = c(-round(y_diffstat["Max."],2), round(y_diffstat["Max."],2))) +
    labs(title="Difference in sample mean VAF between regions", x = "VAF[genomewide] - VAF[alphaSatellite]")+
    theme_pubr()
  
  p <- p_boxes + p_histo + plot_layout(ncol=1, heights = c(3,1))
  if(!save_plot) return(p)
  filename_var <- paste0("s",str_sub(signature, 21,-3),"_vaf_genomewide_vs_alphaSatellite.pdf")
  ggsave(plot=p, device = "pdf", filename=paste0("plots/vaf_distribution/",filename_var), useDingbats=FALSE)
}

plot_vaf("BI_COMPOSITE_SNV_SBS17b_P")
plot_vaf("BI_COMPOSITE_SNV_SBS61_S")
plot_vaf("BI_COMPOSITE_SNV_SBS62_S")
plot_vaf("BI_COMPOSITE_SNV_SBS63_S")
plot_vaf("BI_COMPOSITE_SNV_SBS78_S")


dt3 <- dt3[Donor_ID %in% dt3[, .(n_region=.N), by=Donor_ID][n_region==2, Donor_ID]]


t.test.var <- t.test(x=dt3[region=="genomewide", VAF], y=dt3[region=="alphaSatellite", VAF], paired=TRUE)
t.test.text <- paste0(t.test.var$method, "\nnumber of samples = ",dt3[,uniqueN(Donor_ID)],
                      "\nmean of the differences = ", round(t.test.var$estimate,4),
                      "\np-value < ", scientific(t.test.var$p.value, 2))

dt3_stat <- dt3[, .(VAF=median(VAF), n_snv=sum(N), n_genome=.N), by=region]
p3_boxes <- dt3 %>%
  ggplot(aes(x=factor(region,c("genomewide","alphaSatellite")), y=VAF, group=Donor_ID))+
  geom_boxplot(aes(group=region)) +  geom_point(alpha=0.1)+  geom_line(alpha=0.1)+
  annotate(geom="segment",x=1,xend=2,
           y=max(dt3$VAF)*1.02, yend=max(dt3$VAF)*1.02, color="black")+
  annotate(geom="text", x=1.5, y=max(dt3$VAF)*1.02, label=t.test.text, vjust=-.1)+
  scale_y_continuous(breaks=c(seq(0.3,0.45,length.out=4), dt3_stat$VAF),
                     labels=function(x) round(x,2)) +
  coord_cartesian(ylim = c(0.3, .45)) +
  labs(title="Sample mean VAF per genomic region", x="region",
       subtitle = "Genomes w. S17b and >=20 alphaSatellite SNVs and pr[S17b|SNVs]>=50%") +
  theme_pubr()

dt3_diff <- dt3[,.(VAF_diff=VAF[region=="genomewide"]-VAF[region=="alphaSatellite"]), by = Donor_ID]
dt3_diffstat <- dt3_diff[, summary(VAF_diff)]
p3_histo <- dt3_diff %>%
  ggplot(aes(x=VAF_diff)) + geom_histogram(col="black",fill="grey", bins=30) +
  geom_vline(xintercept=0, color="black", linetype="solid")+
  geom_vline(xintercept = dt3_diffstat["Mean"], linetype="dashed") +
  annotate(geom="text", x=dt3_diffstat["Mean"], y=30, label=round(dt3_diffstat["Mean"],4), hjust=-.1) +
  coord_cartesian(xlim = c(-round(dt3_diffstat["Max."],2), round(dt3_diffstat["Max."],2))) +
  labs(title="Difference in sample mean VAF between regions", x = "VAF[genomewide] - VAF[alphaSatellite]")+
  theme_pubr()

p3 <- p3_boxes + p3_histo + plot_layout(ncol=1, heights = c(3,1))
ggsave(plot=p3, device = "pdf", filename="plots/vaf_distribution/vaf_genomewide_vs_alphaSatellite_s17b.pdf", useDingbats=FALSE)


