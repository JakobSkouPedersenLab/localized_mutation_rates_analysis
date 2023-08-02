# reduce ascertainment bias in mutrates of 11-mers with hotspots
# 20-05-2021


# 1. include non-mutated SNVs
# 2. when subsetting with hotspots:
#       for the mutation rate:
#           include less recurrent SNVs
#           exclude SNVs as or more recurrent than selected threshold


options(scipen = 999)
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr) ; library(patchwork)
  library(scales)
})

#d <- fread("mutprocesses_project/data/snvdata_w_tfbs_hmm_overlap.tsv")
d <- fread("data/signature_snv_filtered_annotated4.tsv", select=c("kmer_11","Donor_ID","eid", "n_pc"))
k11c_dt_all <- fread(file = "results/11mer_count_hg19_PCAWG_v3.tsv")
k11c_dt <- k11c_dt_all[tot_hg19!=0]
#k11c_dt <- k11c_dt_all[n_snv!=0] # gdk
#k11c_dt[, `:=`(triN = NULL, freq_hg19=NULL, center=NULL, mutrate=NULL)]

k11c_dt[, `:=`(mrate = n_snv/tot_hg19/2583*10^6,
               mrate_singletons = rec1/tot_hg19/2583*10^6,
               mrate_rec5p = (n_snv-rec5-rec6-rec7p)/tot_hg19/2583*10^6)]
k11c_dt[, .(n_snv, tot_hg19, mrate, kmer_11, group = "all (1+)", correction = "none")]
k11c_dt[, .(mrate=sum(n_snv)/sum(tot_hg19)/2583*10^6,group = "all (1+)", correction = "none")]
k11c_dt[grep("(2|3|4|5|6|7)", rec), .(mrate=sum(n_snv)/sum(tot_hg19)/2583*10^6,group = "all (1+)", correction = "none")]
k11c_dt[grep("(5|6|7)", rec), .(mrate=sum(n_snv)/sum(tot_hg19)/2583*10^6,group = "all (1+)", correction = "none")]

dlist <- list()
# all snvs / complete span
dlist[[1]] <- k11c_dt[, .(mrate, n_snv, tot_hg19, kmer_11, group = "all (1+)", correction = "none")]
dlist[[2]] <- k11c_dt[grep("(2|3|4|5|6|7)", rec), .(mrate, n_snv, tot_hg19, kmer_11, group = "hotspots (2+)", correction = "none")]
dlist[[3]] <- k11c_dt[grep("(5|6|7)", rec), .(mrate, n_snv, tot_hg19, kmer_11, group = "hotspots (5+)", correction = "none")]
# select snvs / complete span
dlist[[4]] <- k11c_dt[, .(mrate=mrate_singletons, n_snv=rec1, tot_hg19, kmer_11, group = "all (1+)", correction = "snv")]
dlist[[5]] <- k11c_dt[grep("(2|3|4|5|6|7)", rec), .(mrate=mrate_singletons, n_snv=rec1, tot_hg19, kmer_11, group = "hotspots (2+)", correction = "snv")]
dlist[[6]] <- k11c_dt[grep("(5|6|7)", rec), .(mrate=mrate_rec5p, n_snv=rec1+rec2+rec3+rec4, tot_hg19, kmer_11, group = "hotspots (5+)", correction = "snv")]
# select snvs / select span
hs2_npos <- d[n_pc>=2, .(hotspot_span=uniqueN(eid)), by=kmer_11]
dlist[[7]] <- merge(k11c_dt, hs2_npos, by = "kmer_11", all.x=TRUE)
dlist[[7]][, span_cor := ifelse(tot_hg19==hotspot_span | is.na(hotspot_span), tot_hg19, tot_hg19-hotspot_span)]
dlist[[7]] <- dlist[[7]][, .(mrate=rec1/span_cor/2583*10^6, n_snv=rec1, tot_hg19=span_cor, kmer_11, group = "all (1+)", correction = "snv_span")]

dlist[[8]] <- merge(k11c_dt[grep("(2|3|4|5|6|7)", rec)], hs2_npos, by = "kmer_11", all.x=TRUE)
dlist[[8]][, span_cor := ifelse(tot_hg19==hotspot_span | is.na(hotspot_span), tot_hg19, tot_hg19-hotspot_span)]
dlist[[8]] <- dlist[[8]][, .(mrate=rec1/span_cor/2583*10^6, n_snv=rec1, tot_hg19=span_cor, kmer_11, group = "hotspots (2+)", correction = "snv_span")]

hs5_npos <- d[n_pc>=5, .(hotspot_span=uniqueN(eid)), by=kmer_11]
dlist[[9]] <- merge(k11c_dt[grep("(5|6|7)", rec)], hs5_npos, by = "kmer_11", all.x=TRUE)
dlist[[9]][, span_cor := ifelse(tot_hg19==hotspot_span | is.na(hotspot_span), tot_hg19, tot_hg19-hotspot_span)]
dlist[[9]] <- dlist[[9]][, .(mrate=(n_snv-rec5-rec6-rec7p)/span_cor/2583*10^6, n_snv=(n_snv-rec5-rec6-rec7p), tot_hg19=span_cor, kmer_11, group = "hotspots (5+)", correction = "snv_span")]


#dcor <- rbindlist(dlist)
dt <- rbindlist(dlist)
rm(dlist)

#dcor[, mrate_log10:=round(log10(mrate),1)]
dt[, mrate_log10:=round(log10(mrate),1)]

dt[, .(mrate_uw=mean(mrate), mrate_w=sum(n_snv)/sum(tot_hg19)/2583*10^6,
       n_snv=sum(n_snv), tot_hg19=as.numeric(sum(tot_hg19))), by=c("group","correction")]

dt[, sum(mrate==0)/.N*100, by = c("group","correction")]

p_dcor <- dcor %>% #dplyr::sample_n(size=10^5) %>%
  ggplot(aes(x=mrate+0.1)) +
  geom_histogram(bins = 30) +
  #geom_vline(data=d_kmer_long[, mean(rate), by = measure], mapping = aes(xintercept=V1)) +
  #geom_text(data=d_kmer_long[, mean(rate), by = measure], 
  #          mapping = aes(x=V1,y=Inf, label=round(V1,0)), vjust = 2, hjust = -.1) +
  scale_x_log10() + #labels = comma
  scale_y_continuous(name = "count")+#, breaks = 0:3*10^5, limits = c(0,3*10^5)) +
  facet_wrap(factor(correction, levels=c("none","snv","snv_span"))~group, scale="free_y") +
  theme_pubr()

ggsave(plot=p_dcor, device="pdf", width=10, height=10, filename="plots/hist_11mer_rec_split_mutrates_corrected_for_hotspots_v2.pdf")



# necessary statistics before plotting
dt_stat <- dt[, .(mrate_uw=mean(mrate), mrate_w=sum(n_snv)/sum(tot_hg19)/2583*10^6,
                        n_snv=sum(n_snv), tot_hg19=as.numeric(sum(tot_hg19))), by=c("group","correction")]
dt_stat[, fc:=mrate_w/mrate_w[group=="all (1+)"], by = correction]
d_aggr <- dt[, .(tot_span=sum(tot_hg19)), by=c("mrate_log10","group","correction")]
d_aggr[, `:=`(mrate=10^mrate_log10, correction=factor(correction, levels=rev(c("none","snv","snv_span"))))]

# version 1
p_aggr <- d_aggr %>%
  ggplot(aes(x=mrate+0.1, y=tot_span, fill=correction)) +
  geom_col(width=0.09, col="black") +
  geom_vline(data=dt_stat, aes(xintercept=mrate_w)) +
  geom_text(data=dt_stat,aes(x=mrate_w, y=Inf, label=paste0(format(round(mrate_w,2), nsmall = 2),
                                                                    "\n(",format(round(fc,1), nsmall = 1),"x)")), vjust = 2, hjust = -.1) +
  scale_x_log10(name = "Mutation rate (SNV/patient/Mb+0.1)",breaks=10^(-1:4), labels=c("0.1","1","10","100", "1.000", "10.000")) +
  #scale_y_log10(name = "Genomic span (Mb)", labels=function(x) x/10^6) +
  scale_y_continuous(name = "Genomic span (Mb)", labels=function(x) x/10^6) +
  scale_fill_manual(name="Correction",values=c("none"=NA,"snv"="grey30","snv_span"="grey70"),
                    guide = guide_legend(reverse = TRUE)) +
  #facet_grid(factor(correction, levels=c("none","snv","snv_span"))~group, scale="free_y") +
  facet_wrap(factor(correction, levels=c("none","snv","snv_span"))~group, scale="free_y") +
  labs(title = "Mutation rate correction of all (+1) or hotspot-annotated (2+ or 5+) 11-mers",
       subtitle = "Quantification of ascertainment bias in mutation rate of hotspot-selected 11-mers",
       caption = "Mutation rate correction method:
       none (mut.rate = all SNVs / all 11-mer instances)
       conservative (mut.rate = less recurrent SNVs / all 11-mer instances)
       fair (mut.rate = less recurrent SNVs / less recurrent 11-mer instances)") +
  theme_pubr(legend = "right")

ggsave(plot=p_aggr, device="pdf", width=10, height=10, filename="plots/hist_11mer_rec_split_wmutrates_corrected_for_hotspots.pdf")

# version 2
p_aggr <- d_aggr %>%
  ggplot(aes(x=mrate+0.1, y=tot_span, fill=correction)) +
  geom_col(width=0.09, col="black") +
  geom_vline(data=dt_stat, aes(xintercept=mrate_w)) +
  geom_text(data=dt_stat,aes(x=100, y=0, label=paste0(format(round(mrate_w,2), nsmall = 2),
                                                            "\n(",format(round(fc,1), nsmall = 1),"x)")),
            vjust = 0, hjust = -1) +
  scale_x_log10(name = "Mutation rate (SNV/patient/Mb+0.1)",breaks=10^(-1:4), labels=c("0.1","1","10","100", "1.000", "10.000")) +

  scale_y_continuous(name = "Genomic span (Mb)", labels=function(x) x/10^6) +
  scale_fill_manual(name="Correction",values=c("none"=NA,"snv"="grey30","snv_span"="grey70"),
                    guide = guide_legend(reverse = TRUE)) +
  facet_grid(factor(correction, levels=c("none","snv","snv_span"))~group, scale="free_x") +
  coord_flip() +
  labs(title = "Mutation rate correction of all (+1) or hotspot-annotated (2+ or 5+) 11-mers",
       subtitle = "Quantification of ascertainment bias in mutation rate of hotspot-selected 11-mers",
       caption = "Mutation rate correction method:
       none (mut.rate = all SNVs / all 11-mer instances)
       conservative (mut.rate = less recurrent SNVs / all 11-mer instances)
       fair (mut.rate = less recurrent SNVs / less recurrent 11-mer instances)") +
  theme_pubr(legend = "right")

ggsave(plot=p_aggr, device="pdf", width=10, height=10, filename="plots/hist_11mer_rec_split_wmutrates_corrected_for_hotspots_v2.pdf")


