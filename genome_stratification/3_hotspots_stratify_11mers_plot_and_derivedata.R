# info ====
# author:   Gustav Poulsgaard
# date:     16-04-2021
# purpose:  stratify signature-assigned 11-mer data based on recurrence levels
#             plot and save the derived data


# define environment ====
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

options(scipen = 999)
library(tidyverse, quietly = TRUE) # data manipulation
library(data.table, quietly = TRUE) # data manipulation
library(Biostrings, quietly = TRUE) # DNA seq manipulation (consensusMatrix)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE) # assess hg19 seq
library(ggseqlogo, quietly = TRUE) # plot logo plots
library(ggpubr, quietly = TRUE) # pretty themes
library(scales, quietly = TRUE) # pretty numbers
library(patchwork, quietly = TRUE) # plot arrangement

# load data ====
cat(paste0("\nload data\t", Sys.time(), "\n"))
# color references
load("data/PCAWG_donorinfo/cohorts_palette.RData")
col_ref <- setNames(cohorts[1:37, "colcode"], cohorts[1:37, "cohort"]) ; rm(cohorts)
names(col_ref) <- gsub("AdenoCa","AdenoCA", names(col_ref))
# suggested signature etiologies
sigetio <- fread("data/reference_signatures/SignatureAnalyzer_SBS_etiology.tsv")
sbs_columns <- sigetio$BI_COMPOSITE
# signature loads
sigactivity <- fread("data/PCAWG_signatures/SA_COMPOSITE_SNV.activity_relative.WHITELIST_SET.20200226.txt")
# 11-mer stats
k11c_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")

cat(paste0("\nload heavy data\t", Sys.time(),"\n"))
# 11-mer signature stats
kmerstats <- fread(file = "results/signatures_stratify_genomes_and_11mers/kmer_stats_stratified_by_signature.tsv")
ndonors_dt <- kmerstats[, str_split(str_c(donors, collapse = ","), ","), by=baselinesig][, .(n_donor=uniqueN(V1)), by=baselinesig]
kmerstats2 <- merge(kmerstats, ndonors_dt, by = "baselinesig")
rm(kmerstats)
# snv data
snv <- fread("data/signature_snv_filtered_annotated4.tsv", 
             select=c("Donor_ID","kmer_11"))

cat(paste0("\ncompute statistics\t", Sys.time(),"\n"))
kmer_allsnvcounts <- snv[, .(n_snv=.N), by=c("kmer_11","Donor_ID")]

# cancer types
cancertypes_dt <- sigactivity[,1:2]
cancertypes_dt[, Tumor_Type:=gsub("_","-",Tumor_Type)]
cancertypes_dt[, colcode:=col_ref[Tumor_Type]]

# compute expected pfm from hg19 autosomes
b <- Reduce("+",lapply(1:22, function(chr) oligonucleotideFrequency(Hsapiens[[chr]], width = 1)))
pfm_exp <- matrix( (b/sum(b)), nrow = 4, ncol = 11, dimnames = list(names(b)))
rm(b)

# define plotting functions ====

plot_logo_bits_hg19auto <- function(dt, limy = c(0, 1)) {
  # 1. loop over recurrence to access consensusmatrix
  # 2. loop over positions (1-11)
  # 3. loop over bases (1-4)
  n <- length(dt$consensus$censusmatrix_full)
  
  pfm_obs <- lapply(1:n, function(r) dt$consensus$censusmatrix_full[[r]])
  D_KL <- lapply(1:n, function(r) colSums(pfm_obs[[r]]*log2(pfm_obs[[r]]/pfm_exp), na.rm = TRUE))
  
  compute_height <- function(obs, ic) sapply(1:11, function(i) obs[,i]*ic[i])
  heights <- lapply(1:n, function(r) compute_height(pfm_obs[[r]], D_KL[[r]]))
  
  ggseqlogo(heights, method = "custom", seq_type = "dna", ncol = 1)+
    # ggtitle(paste0("rec: ", n_pc[i]),
    #         paste0("n_mut=",comma(nmut[i]),
    #                "\nn_pos=",comma(npos[i]),
    #                "\nspan=",comma(span[i])))+
    scale_x_continuous(expand = c(0,0), breaks = 1:11, labels = 1:11-6, position = "bottom") +
    scale_y_continuous(name = "bits", expand = c(0, 0), breaks = c(0, .5, 1, 1.5, 2),
                       labels = c(0, "", 1, "", 2)) +
    geom_hline(yintercept = c(0)) +
    coord_cartesian(ylim = c(0, 2),  clip = "off") +
    facet_wrap(~seq_group, ncol = 1) +
    theme(#axis.text.x = element_blank(),
      #axis.title.y = element_text(size = 10),
      axis.line.y = element_line(color = "black"), axis.ticks.y = element_line(color = "black"),
      axis.title.x = element_blank(), strip.text = element_blank())
}

plot_mutationrate_test3 <- function(pd_list, plot = TRUE){
  signature_name <- pd_list$mutrate$maxsig_kmer[1]
  recurrence_min <- min(pd_list$mutrate$rec)
  recurrence_opt <- length(pd_list$mutrate$rec)
  
  snv <- readRDS(file = paste0("mutprocesses_project/temp/snv_activitycohort_", str_sub(signature_name, 18, -3),".rds"))
  npatients <- sigactivity[get(signature_name)>=0.05, .N]
  snv[, mrate := (n_mut_kmer/tot_hg19)/npatients*10^6]
  
  snv <- merge(snv,kmer_rec[,.(kmer, rec1, rec2_4, rec5plus)], 
               by.x = "kmer_11", by.y = "kmer", all.x = TRUE)
  snv[, rec := ifelse(rec5plus>0, "5+", ifelse(rec2_4>0, "2-4", ifelse(rec1>1, "1", NA)))]
  
  
  mrate.cohort_mean_var <- snv[, mean(mrate)]
  mrate.sig11mer_mean_var <- snv[maxsig_kmer==signature_name, mean(mrate)]
  
  dt_sig11mer <- snv[maxsig_kmer==signature_name]
  
  dt_hotspots <- rbind(dt_sig11mer[, .(kmer_11, mrate, n_mut_kmer, tot_hg19, rec = "1+")],
                       dt_sig11mer[rec!="1", .(kmer_11, mrate, n_mut_kmer, tot_hg19, rec = "2+")],
                       dt_sig11mer[rec=="5+", .(kmer_11, mrate, n_mut_kmer, tot_hg19, rec = "5+")])
  dt_hotspots <- dt_hotspots[rec >= recurrence_min]
  
  dt_hotspots[, `:=`(mrate.cohort_mean = mrate.cohort_mean_var,
                     mrate.sig11mer_mean = mrate.sig11mer_mean_var)]
  dt_hotspots[, mrate.hotspot_mean := mean(mrate), by = rec]
  
  dt_hotspots_stat <- dt_hotspots[, .(mrate.rec_mean = mean(mrate), 
                                      kmer.rec_n = .N,
                                      mrate.cohort_mean = mrate.cohort_mean_var,
                                      mrate.sig11mer_mean = mrate.sig11mer_mean_var), by = rec]
  #dt3 <- dt_hotspots[rec >= recurrence_min, .(mrate = mean(mrate), .N), by = rec]
  
  if(!plot) return(dt_hotspots)
  ggplot(data = dt_hotspots, aes(mrate)) +
    geom_histogram(aes(y = stat(count), fill = stat(x)), bins = 30, col = "black") +
    # mutation rate mean for previous data subset
    geom_vline(data = dt_hotspots_stat, aes(xintercept = mrate.sig11mer_mean), col = "grey", linetype = "dashed", size = 1) +
    geom_text(data = dt_hotspots_stat, aes(x = min(dt_hotspots$mrate), y = Inf, label = paste0("n = ",comma(kmer.rec_n))), 
              vjust = 1, hjust = 0) +
    
    # mutation rate mean for each recurrence group
    geom_text(data = dt_hotspots_stat, aes(x = mrate.rec_mean, y = Inf, label = paste("mean =",comma(mrate.rec_mean, 1))), vjust = 1) +
    geom_vline(data = dt_hotspots_stat, aes(xintercept = mrate.rec_mean), col = "black", linetype = "solid", size = 1) +
    # axis and color scaling
    scale_x_log10(name = "mutation rate (SNV / Mb / patient)", labels = function(x) comma(x, 1), expand = c(0,0)) + #, breaks = c(1, 100, 10^4)
    scale_y_continuous(name = "count", labels = function(x) comma(x, 1), expand = c(0,0)) +
    scale_fill_gradient2(low = "blue", high = "brown", midpoint = log10(dt_hotspots_stat$mrate.sig11mer_mean[1])) +
    # facet and theme
    coord_cartesian(clip = "off") +
    facet_grid(rec~., scales = "free_y") +
    theme_pubr(legend = "none") +
    theme(strip.background = element_blank(), strip.text = element_blank())
}

plot_cancercohorts_test <- function(dt) {
  recurrence_opt <- length(dt$mutrate$rec)
  ggplot(data = dt$cohort[[1]], aes(x = factor(rec, levels = rev(c("1", "2-4", "5+"))), y = N, fill = Tumor_Type)) + 
    geom_col(width = .3, position = position_fill()) +
    scale_y_continuous(breaks = c(0,.25,.5,.75,1), labels = c(0,"",.5,"",1),
                       expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = col_ref) +
    labs(x = "Recurrence", y = "SNVs / tumor type")+
    geom_text(inherit.aes = FALSE, data = dt$cohort[[2]], 
              aes(x = rec, y = 0.5, label = paste0(comma(dt$mutrate$n_mut_kmer), " SNVs\n",
                                                   Donor_ID," patients  |  ", Tumor_Type," tumor types")),
              size = 2, nudge_x = ifelse(recurrence_opt==1, 0.2, 0.3), 
              hjust = 0.5, vjust = 1) +
    coord_flip(clip = "off") +
    theme_pubr(legend = "none") + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8),
          strip.background = element_blank(), panel.background = element_blank(),
          plot.background = element_blank())
}

full_plot_test <- function(sig_name, saveggplot = FALSE){
  plotdata_list <- readRDS(file = paste0("mutprocesses_project/results/mutablemotif_data/summary_kmer_signature_mutrate/plotdata_list_", 
                                         sig_name, "_rec1_2-4_5.rds"))
  
  
  filename_string <- paste0(str_sub(sig_name, 18, -3),
                            "_mutablemotif_rec_",
                            gsub(pattern = "\\+", replacement = "", x = paste0(unique(plotdata_list$mutrate$rec), collapse = "_")),
                            "_newlayout.pdf")
  #if(!saveggplot) return(plot_mutationrate_test3(plotdata_list, T))
  
  # p_newlayout <- cowplot::plot_grid(
  #   plot_logo_bits_hg19auto(plotdata_list) + 
  #     ggtitle(sig_name),
  #   plot_logo_bits(plotdata_list),
  #   plot_logo_freq(plotdata_list),
  #   plot_mutationrate_test3(plotdata_list), # test vs test2 vs test3
  #   plot_cancercohorts_test(plotdata_list), # test vs ""
  #   align = "h", axis = "tb", nrow = 1, rel_widths = c(8, 8, 10, 6, 6))
  #require(patchwork)
  p_newlayout <- cowplot::plot_grid(
    plot_logo_bits_hg19auto(plotdata_list) + 
      ggtitle(sig_name),
    # plot_logo_bits(plotdata_list),
    # plot_logo_freq(plotdata_list),
    plot_mutationrate_test3(plotdata_list), # test vs test2 vs test3
    plot_cancercohorts_test(plotdata_list), # test vs ""
    align = "h", axis = "tb", nrow = 1, rel_widths = c(8, 10, 6))
  
  if(saveggplot){
    path <- "mutprocesses_project/plots/hotspots_stratify_11mers/"
    ggsave(plot = p_newlayout, device = "pdf",
           filename = paste0(path, filename_string),
           width = 12, height = 4, units = "in")
    cat(paste0("\nsaved ggplot\n\t", filename_string, "\nto\n\t", path, "\n"))
  }
  return(p_newlayout)
}

logohighlight_plot <- function(sig_name, recurrence = "5+", saveggplot = FALSE){
  plotdata_list <- readRDS(file = paste0("mutprocesses_project/results/mutablemotif_data/summary_kmer_signature_mutrate/plotdata_list_", 
                                         sig_name, "_rec1_2-4_5.rds"))
  npatients <- plotdata_list$cohort[[2]][rec==1, Donor_ID]
  pd_l_rec <- list(mutrate = plotdata_list$mutrate[rec==recurrence],
                   cohort = list(plotdata_list$cohort[[1]][rec==recurrence],
                                 plotdata_list$cohort[[2]][rec==recurrence]),
                   consensus = plotdata_list$consensus[rec==recurrence])

    p_newlayout <- plot_logo_bits_hg19auto(pd_l_rec) +
    plot_mutationrate_test3(pd_l_rec) +
    plot_cancercohorts_test(pd_l_rec) + 
    patchwork::plot_layout(nrow = 1, widths = c(8, 10, 6))
  

  
  if(saveggplot){
    filename_string <- paste0("mutprocesses_project/plots/hotspots_stratify_11mers/",
                              str_sub(sig_name,18,-3),"_highlight_rec5_logo.pdf")
    ggsave(plot = p_newlayout, device = "pdf", width = 13, height = 2,
           filename = filename_string)
    cat(paste0("\nsaved ggplot\n\t", filename_string, "\n\n"))
  }
  return(p_newlayout)
}


###### mutation rate distribution stratified by signature and recurrence


# define functions ====
subset_data <- function(signature, recurrence_regex){
  # subset data
  s3_dt <- kmerstats2[baselinesig==signature & baselinesig==maxsig_kmer & grepl(recurrence_regex, rec)]
  s3_dt[, mrate:=n_snv/tot_hg19/n_donor*10^6]
  s3_dt[, mrate_log10:=round(log10(mrate),1)]
  return(s3_dt)
}
pfm <- function(seq, as.prob=F) consensusMatrix(DNAStringSet(seq), as.prob = as.prob)[1:4,]
plot_logo_bits_hg19auto_v2 <- function(subset_dt) {
  # 1. loop over recurrence to access consensusmatrix
  # 2. loop over positions (1-11)
  # 3. loop over bases (1-4)
  
  pfm_obs <- subset_dt[, pfm(rep(kmer_11, n_snv), as.prob=T)]
  D_KL <- colSums(pfm_obs*log2(pfm_obs/pfm_exp), na.rm = TRUE)
  
  compute_height <- function(obs, ic) sapply(1:11, function(i) obs[,i]*ic[i])
  heights <- compute_height(pfm_obs, D_KL)
  
  ggseqlogo(heights, method = "custom", seq_type = "dna", ncol = 1)+
    # ggtitle(paste0("rec: ", n_pc[i]),
    #         paste0("n_mut=",comma(nmut[i]),
    #                "\nn_pos=",comma(npos[i]),
    #                "\nspan=",comma(span[i])))+
    scale_x_continuous(expand = c(0,0), breaks = 1:11, labels = 1:11-6, position = "bottom") +
    scale_y_continuous(name = "bits", expand = c(0, 0), breaks = c(0, .5, 1, 1.5, 2),
                       labels = c(0, "", 1, "", 2)) +
    geom_hline(yintercept = c(0)) +
    coord_cartesian(ylim = c(0, 2),  clip = "off") +
    facet_wrap(~seq_group, ncol = 1) +
    theme(#axis.text.x = element_blank(),
      #axis.title.y = element_text(size = 10),
      axis.line.y = element_line(color = "black"), axis.ticks.y = element_line(color = "black"),
      axis.title.x = element_blank(), strip.text = element_blank())
}
plot_mutrate_aggregated_by_span <- function(subset_dt){
  # aggregate mutationrates by genome span
  kmer_aggr_dt <- as.data.table(aggregate(tot_hg19 ~ mrate_log10, data = subset_dt, FUN = sum))
  kmer_aggr_dt[, mrate:=10^mrate_log10]
  
  # compute weighted mean and axis breaks
  mrate_weightedmean <- subset_dt[, sum(n_snv)/sum(tot_hg19)/n_donor[1]*10^6]
  x_breaks <- 10^(-1:5)
  y_breaks <- pretty(seq(from = 0, to = max(kmer_aggr_dt$tot_hg19), length.out = 4), n = 4)
  y_breaks_sec <- y_breaks/sum(k11c_dt$tot_hg19)
  gspanlabel <- subset_dt[, paste0(format(round(sum(tot_hg19)/10^6,2),nsmall=2, big.mark=",")," Mb\n(",format(round(sum(tot_hg19)/sum(k11c_dt$tot_hg19)*100,1),nsmall=1),"%)")]
  # plotting
  p_histo <- ggplot(kmer_aggr_dt, aes(x=mrate, y=tot_hg19, fill=log10(mrate))) +
    geom_col(col="black", width = 0.09)+
    geom_vline(xintercept = 5.96, linetype = "dashed", col = "black", size = 1) +
    annotate("text", x=5.96, y=Inf, vjust=2, hjust=-.1, label="1x") +
    geom_vline(xintercept = mrate_weightedmean, linetype = "dashed", col = "black", size = 1) +
    annotate("text", x=mrate_weightedmean, y=Inf, vjust=1, hjust=-.1,
             label=paste0(format(round(mrate_weightedmean,2),nsmall=2),"\n(",
                          format(round(mrate_weightedmean/5.96,1),nsmall=1),"x)")) +
    annotate("text", x=10^4, y=Inf, vjust=1, hjust=-.1, label=gspanlabel) +
    scale_x_log10(name="mutation rate (SNV/patient/Mb)",
                  breaks=x_breaks,
                  labels=function(x) prettyNum(x, big.mark=",")) +
    scale_y_continuous(name="Genomic span (Mb)", 
                       labels = function(y) round(y/10^6, 1), expand=c(0,0),
                       breaks = y_breaks,
                       sec.axis = sec_axis(name = "Genomic span (% of WG)", 
                                           trans = function(y) y/sum(k11c_dt$tot_hg19),
                                           breaks = y_breaks_sec,
                                           labels = function(y) format(round(y*100, 2), nsmall=2)))+
    scale_fill_gradient2(low = "blue", high = "brown", midpoint = log10(5.96)) +
    #labs(title=paste("Signature",str_sub(dt_aggr$sig[1], 20, -3))) +
    coord_cartesian(xlim=range(x_breaks), ylim=c(0,max(y_breaks))) +
    theme_pubr(legend = "none") +
    theme(axis.title.y.right = element_text(angle=90))
  return(p_histo)
}
plot_cancertype_gspan <- function(subset_dt){
  donors_var <- subset_dt$donors %>% str_c(collapse=",") %>% 
    str_split(",") %>% unlist() %>% unique()
  kmercounts_per_donor3 <- kmer_allsnvcounts[kmer_11%in%subset_dt$kmer_11 & Donor_ID%in%donors_var] %>%
    merge(cancertypes_dt, by.x="Donor_ID", by.y="icgc_donor_id") %>%
    merge(k11c_dt[,.(tot_hg19, kmer_11)], by="kmer_11")
  
  kmercounts_per_donor4 <- kmercounts_per_donor3[, .(tot_hg19=sum(tot_hg19)),by=c("Tumor_Type")]
  kmercounts_per_donor6 <- kmercounts_per_donor3[, .(n_snv=sum(n_snv)),by=c("Donor_ID","Tumor_Type","colcode")] %>%
    merge(kmercounts_per_donor4, by = "Tumor_Type")
  kmercounts_per_donor6 <- kmercounts_per_donor6[order(Tumor_Type)][,xval:=1:.N]
 
  p_types <- kmercounts_per_donor6 %>%
    ggplot(aes(factor(xval), tot_hg19))+
    #geom_col(width=1, col="black")+ # add black edges
    geom_col(aes(fill=colcode), width=1)+
    scale_x_discrete(name="Genomes",expand=c(0,0))+
    scale_y_continuous(name="Genomic span (Mb)", labels=function(x) comma(x/10^6), expand=c(0,0),
                       sec.axis=sec_axis(name="Genomic span (% of WG)",
                                         trans=~./sum(k11c_dt$tot_hg19)*100,
                                         labels=function(x) format(round(x,2), nsmall=2)))+
    scale_fill_identity()+
    facet_grid(.~Tumor_Type, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_pubr() +
    theme(axis.title.y.right = element_text(angle=90), axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"), 
          strip.text = element_text(angle=90), strip.background = element_blank())
  return(p_types)
} #ooo
plot_cancertype_nsnv <- function(subset_dt){
  donors_var <- subset_dt$donors %>% str_c(collapse=",") %>% 
    str_split(",") %>% unlist() %>% unique()
  kmercounts_per_donor3 <- kmer_allsnvcounts[kmer_11%in%subset_dt$kmer_11 & Donor_ID%in%donors_var] %>%
    merge(cancertypes_dt, by.x="Donor_ID", by.y="icgc_donor_id") %>%
    merge(k11c_dt[,.(tot_hg19, kmer_11)], by="kmer_11")
  
  kmercounts_per_donor6 <- kmercounts_per_donor3[, .(n_snv=sum(n_snv), tot_hg19=sum(tot_hg19)),
                                                 by=c("Donor_ID","Tumor_Type","colcode")]
  kmercounts_per_donor6 <- kmercounts_per_donor6[order(Tumor_Type, n_snv)][,xval:=1:.N]
  #return(kmercounts_per_donor6)
  
  p_types <- kmercounts_per_donor6 %>%
    ggplot(aes(factor(xval), n_snv))+
    #geom_col(width=1, col="black")+ # add black edges
    geom_col(aes(fill=colcode), width=1)+
    scale_x_discrete(name=paste0("Genomes (n=",length(donors_var),")"),expand=c(0,0))+
    scale_y_log10(name="SNV count",labels=comma,expand=c(0,0))+
    #scale_y_continuous(name="Genomic span (Mb)", labels=function(x) comma(x/10^6), expand=c(0,0),
    #                   sec.axis=sec_axis(name="Genomic span (% of WG)",
    #                                     trans=~./sum(k11c_dt$tot_hg19)*100,
    #                                     labels=function(x) format(round(x,2), nsmall=2)))+
    scale_fill_identity()+
    facet_grid(.~Tumor_Type, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_pubr() +
    theme(axis.title.y.right = element_text(angle=90), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"), 
          strip.text = element_text(angle=0), strip.background = element_blank())
  return(p_types)
} #ooo
plot_cancercohorts_v2 <- function(subset_dt) {
  donors_var <- subset_dt$donors %>% str_c(collapse=",") %>% 
    str_split(",") %>% unlist() %>% unique()
  donor_stat <- cancertypes_dt[icgc_donor_id%in%donors_var, 
                               .N, by=c("Tumor_Type","colcode")]

  p_types <- donor_stat %>%
    ggplot(aes("", N))+
    #geom_col(width=1, col="black")+ # add black edges
    geom_col(aes(fill=colcode), width=1)+
    scale_x_discrete(name="",expand=c(0,0))+ # name=paste0("Genomes (n=",length(donors_var),")"),
    scale_y_continuous(name="",expand=c(0,0), breaks=ceiling(seq(0,length(donors_var),length.out=4)))+
    scale_fill_identity()+
    coord_flip(clip = "off") +
    theme_pubr()+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
    # theme(axis.title.y.right = element_text(angle=90), axis.ticks.x = element_blank(),
    #       axis.line.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"), 
    #       strip.text = element_text(angle=0), strip.background = element_blank())
  return(p_types)

}

plot_pipeline_rec1.2.5 <- function(signature){
  recurrence_regexes=c("(1)","(2|3|4|5|6|7+)","(5|6|7+)")
  
  p_list <- list()
  
  for(i in 1:3){
    rec <- recurrence_regexes[i]
    x <- subset_data(signature, rec)
    p_list[[paste0("type_",rec)]] <- plot_cancertype_nsnv(x) #p_types 
    p_list[[paste0("logo_",rec)]] <- plot_logo_bits_hg19auto_v2(x) # p_logo
    p_list[[paste0("mrate_",rec)]] <- plot_mutrate_aggregated_by_span(x) + # p_histo
      ggtitle(label = paste0("Signature ",str_sub(signature,20,-3)," rec=",rec))
  }
  # # singletons
  # x1 <- subset_data(signature, recurrence_regexes[1])
  # p_histo1 <- plot_mutrate_aggregated_by_span(x1)
  # p_logo1 <- plot_logo_bits_hg19auto_v2(x1)
  # p_types1 <- plot_cancertype_nsnv(x1) + ggtitle(label = paste0(signature,))
  # # all hotspots
  # x2 <- subset_data(signature, recurrence_regexes[2])
  # p_histo2 <- plot_mutrate_aggregated_by_span(x2)
  # p_logo2 <- plot_logo_bits_hg19auto_v2(x2)
  # p_types2 <- plot_cancertype_nsnv(x2)
  # # 5+ hotspots
  # x5p <- subset_data(signature, recurrence_regexes[3])
  # p_histo5p <- plot_mutrate_aggregated_by_span(x5p)
  # p_logo5p <- plot_logo_bits_hg19auto_v2(x5p)
  # p_types5p <- plot_cancertype_nsnv(x5p)
  # plot alignment
  design = "
  133\n233\n233
  466\n566\n566
  799\n899\n899"
  
  #p_align <- p_list[[1]] + p_list[[2]] + p_list[[3]] + 
  # plot_layout(design="
  #             133
  #             233
  #             233")
  
  p_align <- Reduce("+", p_list) + plot_layout(design=design)
  # p_align <-
  #   p_types1 + p_logo1 + p_histo1 +
  #   p_types2 + p_logo2 + p_histo2 +
  #   p_types5p + p_logo5p + p_histo5p +
  #   plot_layout(design=design)
  return(p_align)
} #ooo
plot_pipeline_rec1.2.5_flat <- function(signature){
  recurrence_regexes=c("(1|2|3|4|5|6|7+)","(2|3|4|5|6|7+)","(5|6|7+)")
  
  p_list <- list()
  
  for(i in 1:3){
    rec <- recurrence_regexes[i]
    x <- subset_data(signature, rec)
    p_list[[paste0("logo_",rec)]] <- plot_logo_bits_hg19auto_v2(x) # p_logo
    p_list[[paste0("mrate_",rec)]] <- plot_mutrate_aggregated_by_span(x) + # p_histo
      ggtitle(label = paste0("Signature ",str_sub(signature,20,-3)," rec=",rec))
    #p_list[[paste0("type_",rec)]] <- plot_cancertype_nsnv(x) #p_types
    p_list[[paste0("type_",rec)]] <- plot_cancercohorts_v2(x)
  }
  # plot alignment
  p_align <- Reduce("+", p_list) + plot_layout(nrow=3, ncol=3, widths = c(2, 4, 1))
  # p_align <-
  #   p_types1 + p_logo1 + p_histo1 +
  #   p_types2 + p_logo2 + p_histo2 +
  #   p_types5p + p_logo5p + p_histo5p +
  #   plot_layout(design=design)
  return(p_align)
}

# p_align <- plot_pipeline_rec1.2.5_flat("BI_COMPOSITE_SNV_SBS17b_P")
# 
# 
# ggsave(plot=p_align, device="pdf", width=14, height=7,
#        filename = "plots/hotspots_stratify_11mers/test_logohisto.pdf")


cat(paste0("\nplot\t", Sys.time(),"\n"))
# sbs_ord <- unique(kmerstats2$baselinesig)[order(factor(unique(kmerstats2$baselinesig), levels=sbs_columns))]
# pdf(file = "plots/hotspots_stratify_11mers/allsigs_overview.pdf", width = 14, height = 7)
# lapply(sbs_ord, plot_pipeline_rec1.2.5_flat)
# dev.off()


sbs_select_ord <- c("BI_COMPOSITE_SNV_SBS1_P", "BI_COMPOSITE_SNV_SBS17a_P", "BI_COMPOSITE_SNV_SBS17b_P", 
"BI_COMPOSITE_SNV_SBS9_P", "BI_COMPOSITE_SNV_SBS19_P", "BI_COMPOSITE_SNV_SBS72_P", # unknown
"BI_COMPOSITE_SNV_SBS36_P","BI_COMPOSITE_SNV_SBS70_P", # Unknown

"BI_COMPOSITE_SNV_SBS2_P","BI_COMPOSITE_SNV_SBS13_P","BI_COMPOSITE_SNV_SBS69_P", # APOBEC

"BI_COMPOSITE_SNV_SBS40_P", "BI_COMPOSITE_SNV_SBS68_P",# dinucleotide OR homopolymer runs
"BI_COMPOSITE_SNV_SBS77_P","BI_COMPOSITE_SNV_SBS71_P", "BI_COMPOSITE_SNV_SBS64_P", "BI_COMPOSITE_SNV_SBS22_P", # borders of hp runs
"BI_COMPOSITE_SNV_SBS16_P", "BI_COMPOSITE_SNV_SBS18_P","BI_COMPOSITE_SNV_SBS33_P",# dinucleotide repeats

"BI_COMPOSITE_SNV_SBS7a_S","BI_COMPOSITE_SNV_SBS7b_S","BI_COMPOSITE_SNV_SBS7c_S","BI_COMPOSITE_SNV_SBS38_S", # UV
"BI_COMPOSITE_SNV_SBS55_S","BI_COMPOSITE_SNV_SBS65_S","BI_COMPOSITE_SNV_SBS67_S", "BI_COMPOSITE_SNV_SBS75_S",# UV
"BI_COMPOSITE_SNV_SBS6_S", "BI_COMPOSITE_SNV_SBS14_S", "BI_COMPOSITE_SNV_SBS15_S", "BI_COMPOSITE_SNV_SBS14_S", # MMR
"BI_COMPOSITE_SNV_SBS21_S", "BI_COMPOSITE_SNV_SBS26_S", "BI_COMPOSITE_SNV_SBS44_S", "BI_COMPOSITE_SNV_SBS73_S",# MMR
"BI_COMPOSITE_SNV_SBS74_S", "BI_COMPOSITE_SNV_SBS76_S", "BI_COMPOSITE_SNV_SBS79_S", # MMR
"BI_COMPOSITE_SNV_SBS10a_S", "BI_COMPOSITE_SNV_SBS61_S", "BI_COMPOSITE_SNV_SBS62_S",# POLE
"BI_COMPOSITE_SNV_SBS63_S", "BI_COMPOSITE_SNV_SBS66_S", "BI_COMPOSITE_SNV_SBS78_S", # POLE
"BI_COMPOSITE_SNV_SBS35_P" # platinum
)
#unique(kmerstats2$baselinesig)[order(factor(unique(kmerstats2$baselinesig), levels=sbs_columns))]
pdf(file = "plots/hotspots_stratify_11mers/selectsigs_overview_v2.pdf", width = 14, height = 7)
lapply(sbs_select_ord, plot_pipeline_rec1.2.5_flat)
dev.off()

cat(paste0("\nDONE\t", Sys.time(),"\n"))

p_test <- plot_pipeline_rec1.2.5_flat("BI_COMPOSITE_SNV_SBS17b_P")
ggsave(p_test, device="pdf", width=10, height=5,
       filename="plots/hotspots_stratify_11mers/test.pdf")
# s3_dt <- subset_data("BI_COMPOSITE_SNV_SBS17b_P", "(5|6|7+)")
# p_histo <- plot_mutrate_aggregated_by_span(s3_dt)
# p_logo <- plot_logo_bits_hg19auto_v2(s3_dt)
# p_types <- plot_cancertype_nsnv(s3_dt)
# 
# design <- "
# 22#
# 331
# 331
# "
# #p_align <- p_types + p_logo + p_histo + plot_layout(design=design)
# p_align <- p_types +  p_logo + p_histo +plot_layout(design="
#                                                     1333
#                                                     2333
#                                                     2333")

#####

# # apply plotting functions and save plots ====
# full_plot_test("BI_COMPOSITE_SNV_SBS17b_P", saveggplot = TRUE)
# full_plot_test("BI_COMPOSITE_SNV_SBS7a_S", saveggplot = TRUE)
# 
# 
# logohighlight_plot("BI_COMPOSITE_SNV_SBS62_S", saveggplot = TRUE) # Pol epsilon
# 
# logohighlight_plot("BI_COMPOSITE_SNV_SBS9_P", saveggplot = TRUE) # Pol eta
# logohighlight_plot("BI_COMPOSITE_SNV_SBS19_P", saveggplot = TRUE) # unknown
# logohighlight_plot("BI_COMPOSITE_SNV_SBS72_P", saveggplot = TRUE) # lymph
# logohighlight_plot("BI_COMPOSITE_SNV_SBS18_P", saveggplot = TRUE) # ROS
# logohighlight_plot("BI_COMPOSITE_SNV_SBS30_P", saveggplot = TRUE) # unknown
# 
# logohighlight_plot("BI_COMPOSITE_SNV_SBS37_P", saveggplot = TRUE)
# 
# logohighlight_plot("BI_COMPOSITE_SNV_SBS17a_P", saveggplot = TRUE)
# #logohighlight_plot("BI_COMPOSITE_SNV_SBS17b_P", saveggplot = TRUE)
# 
# logohighlight_plot("BI_COMPOSITE_SNV_SBS68_P", saveggplot = TRUE) # homopolymer boundary T->A
# logohighlight_plot("BI_COMPOSITE_SNV_SBS28_P", saveggplot = TRUE)
# 
# logohighlight_plot("BI_COMPOSITE_SNV_SBS71_P", saveggplot = TRUE)
# logohighlight_plot("BI_COMPOSITE_SNV_SBS64_P", saveggplot = TRUE)
# 
# #logohighlight_plot("BI_COMPOSITE_SNV_SBS7a_S", saveggplot = TRUE)
# 
# # save derived data ====
# derivedData_hotspot <- function(sig_name, summariseData = TRUE, saveData = FALSE){
#   path_input <- "mutprocesses_project/results/mutablemotif_data/summary_kmer_signature_mutrate/"
#   input <- readRDS(file = paste0(path_input, "plotdata_list_",sig_name, "_rec1_2-4_5.rds"))
#   path_output <- "mutprocesses_project/results/hotspots_stratify_11mers/"
#   derivedData <- plot_mutationrate_test3(input, plot = FALSE)
#   derivedData[, `:=`(signature = sig_name, analysis = "hotspot")]
#   output <- derivedData
#   filename_output <- paste0(str_sub(sig_name,18,-3), "_hotspot")
#   if(summariseData) {
#     output <- derivedData[, .(pfm_obs = list(consensusMatrix(DNAStringSet(kmer_11))[1:4,]),
#                               mrate.cohort_mean = mrate.cohort_mean[1], 
#                               mrate.sig11mer_mean = mrate.sig11mer_mean[1],
#                               mrate.hotspot_mean = mrate.hotspot_mean[1]),
#                           by = c("rec", "analysis","signature")]
#     filename_output <- paste0(filename_output, "_summarydata")
#   }
#   if(saveData) saveRDS(output, file = paste0(path_output, filename_output, ".rds"))
#   return(output)
# }
# 
# 
# derivedSummaryData_allsigs <- lapply(sbs_columns[-c(35,58,60)],
#     function(x) derivedData_hotspot(sig_name = x, summariseData = TRUE, saveData = FALSE))
# 
# derivedSum_dt <- rbindlist(derivedSummaryData_allsigs)
# saveRDS(derivedSum_dt, file = "mutprocesses_project/results/hotspots_stratify_11mers/allsignatures_hotspot_summarydata.rds")
# 
# 
# derivedData_dt <- rbindlist(
#   lapply(sbs_columns[-c(35,58,60)],
#          function(x) derivedData_hotspot(sig_name = x,
#                                          summariseData = FALSE,
#                                          saveData = FALSE)
#          )
#   )
# 
# fwrite(derivedData_dt, sep = "\t", file = "mutprocesses_project/results/hotspots_stratify_11mers/allsignatures_hotspot_deriveddata.tsv")
