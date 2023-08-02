# Author          Gustav Poulsgaard
# Date created    16-04-2021
# Purpose
#   calc summary of kmer-instances in- and outside of genomic elements

# Setup environment ====
cat(paste0("Set environment\t",Sys.time(),"\n"))
options(scipen = 999)
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")
suppressPackageStartupMessages({
library(data.table) ; library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19) ; library(Biostrings)
library(ggplot2) ; library(ggseqlogo) ; library(ggpubr) ; library(patchwork)
library(scales)
})
  
# Load data ====
cat(paste0("Load data\t",Sys.time(),"\n"))
donor_kmer_feature_dt <- fread("results/genomicregions_stratify_11mers/allgenomes_allkmers_overlap_enchmmfeatures.tsv")
kmerstats <- fread(file = "results/signatures_stratify_genomes_and_11mers/kmer_stats_stratified_by_signature.tsv")

# patient data relative signature loads
#sigactivity <- fread(file = "data/PCAWG_signatures/SA_COMPOSITE_SNV.activity_relative.WHITELIST_SET.20200226.txt")

#derivedData_dt <- fread(file = "mutprocesses_project/results/hotspots_stratify_11mers/allsignatures_hotspot_deriveddata.tsv")

k11c_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")
kmer_dict <- fread(file = "results/kmer_dict/kmer_background_countfeature_revcomplcollapse.tsv")
setnames(kmer_dict, old = c("count_hg19"), new = c("tot_hg19"))

cat(paste0("Compute statistics\t",Sys.time(),"\n"))
features_of_interest <- colnames(kmer_dict)[3:ncol(kmer_dict)]
kmer_dict[, count_hg19 := sum(.SD), by = seq, .SDcols = features_of_interest]
kmer_dict_long <- melt.data.table(data = kmer_dict, id.vars = c("seq", "tot_hg19"),  #"count_hg19", 
                                  measure.vars = colnames(kmer_dict)[3:17], variable.name = "feature", 
                                  value.name = "refcount_feature")[order(seq, feature)]

match_enchmm_text <- "(?<=[:digit:]{1,2}_)[:alpha:]*[:punct:]*[:alpha:]*"
kmer_dict_long[, `:=`(kmer_11=seq, feature_collapse=str_extract(feature, match_enchmm_text))]
kmer_dict_long_collapse <- kmer_dict_long[, .(refcount_feature=sum(refcount_feature)), by=c("kmer_11","feature_collapse", "tot_hg19")]
rm(kmer_dict_long)

# compute hg19 base distribution to use as background for logos
b <- Reduce("+",lapply(1:22, function(chr) oligonucleotideFrequency(Hsapiens[[chr]], width = 1)))
pfm_exp <- matrix( (b/sum(b)), nrow = 4, ncol = 11, dimnames = list(names(b)))
rm(b)


# filtering
subset_enchmmcount <- function(signature, recurrence_regex){
  # donor criteria
  donor_var <- kmerstats[grepl(recurrence_regex, rec) & 
                           baselinesig==maxsig_kmer &
                           baselinesig==signature, donors] %>% 
    str_c(collapse = ",") %>% str_split(",") %>% unlist() %>% unique()
  # kmer criteria
  kmer_var <- kmerstats[grepl(recurrence_regex, rec) & 
                          baselinesig==maxsig_kmer &
                          baselinesig==signature, kmer_11]
  # filter on the two criterias
  out <- donor_kmer_feature_dt[Donor_ID%in%donor_var & kmer_11%in%kmer_var]
  out[,signature:=signature]
  return(out)
}

subset_enchmmcount2 <- function(signature, recurrence_regex){
  # donor criteria
  donor_var <- kmerstats[grepl(recurrence_regex, rec) & 
                           baselinesig==maxsig_kmer &
                           baselinesig==signature, donors] %>% 
    str_c(collapse = ",") %>% str_split(",") %>% unlist() %>% unique()
  # kmer criteria
  kmer_var <- kmerstats[grepl(recurrence_regex, rec) & 
                          baselinesig==maxsig_kmer &
                          baselinesig==signature, kmer_11]
  # filter on the two criterias
  x <- donor_kmer_feature_dt[Donor_ID%in%donor_var & kmer_11%in%kmer_var]
  x[, signature:=signature]
  
  match_enchmm_text <- "(?<=[:digit:]{1,2}_)[:alpha:]*[:punct:]*[:alpha:]*"
  x[, feature_collapse:=str_extract(feature, match_enchmm_text)]
  x2 <- x[, .(n_snv=sum(N), donors=paste(Donor_ID,collapse=",")), by=c("kmer_11","feature_collapse", "signature")]
  
  y <- merge(x2, kmer_dict_long_collapse, by=c("kmer_11","feature_collapse"))
  y[, n_donor:=uniqueN(x$Donor_ID)]
  y[, mrate:=n_snv/refcount_feature/n_donor*10^6]
  y[, mrate_log10:=round(log10(mrate),1)]
  
  return(y)
}


# regular expression explanation:
# (?<=[:digit:]{1,2}_)
# look ahead of target for 1 or 2 digits followed by 1 underscore 
# [:alpha:]*[:punct:]*[:alpha:]*
# match >=0 letters, followed >=0 special characters, followed by >=0 letters

# feature_collapse_var_ord <- unique(str_extract(feature_var_ord, match_enchmm_text))
feature_collapse_var_ord<- c("Active_Promoter","Weak_Promoter","Poised_Promoter","Strong_Enhancer","Weak_Enhancer",
     "Insulator","Txn_Transition","Txn_Elongation","Weak_Txn","Repressed","Heterochrom/lo","Repetitive/CNV")
#y[, feature_collapse:=factor(str_extract(feature, match_enchmm_text), levels=feature_collapse_var_ord)]


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
plot_mutrate_aggregated_by_span <- function(subset_dt, region_title=NULL){
  # aggregate mutationrates by genome span
  kmer_aggr_dt <- as.data.table(aggregate(refcount_feature ~ mrate_log10, data = subset_dt, FUN = sum))
  kmer_aggr_dt[, mrate:=10^mrate_log10]
  
  #
  signature_var <- str_sub(unique(subset_dt$signature),20,-3)
  #if(is.null(region_title)) region_title <- unique(subset_dt$feature_collapse)
  feature_var <- unique(subset_dt$feature_collapse)
  # compute weighted mean and axis breaks
  mrate_weightedmean <- subset_dt[, sum(n_snv)/sum(refcount_feature)/n_donor[1]*10^6]
  x_breaks <- 10^(-1:5)
  #x_breaks <- 10^(0:4)
  y_breaks <- pretty(seq(from = 0, to = max(kmer_aggr_dt$refcount_feature), length.out = 4), n = 4)
  y_breaks_sec <- y_breaks/sum(k11c_dt$tot_hg19) # y_breaks/sum(k11c_dt$tot_hg19)

  # plotting
  p_histo <- ggplot(kmer_aggr_dt, aes(x=mrate, y=refcount_feature, fill=log10(mrate))) +
    geom_col(col="black", width = 0.09)+
    geom_vline(xintercept = 5.96, linetype = "dashed", col = "black", size = 1) +
    annotate("text", x=5.96, y=Inf, vjust=2, hjust=-.1, label="1x") +
    geom_vline(xintercept = mrate_weightedmean, linetype = "dashed", col = "black", size = 1) +
    annotate("text", x=mrate_weightedmean, y=Inf, vjust=2, hjust=-.1, label=paste0(round(mrate_weightedmean/5.96,1),"x")) +
    scale_x_log10(name="mutation rate (SNV/patient/Mb)",
                  breaks=x_breaks,
                  labels=function(x) ifelse(x>1, prettyNum(x, big.mark=","), signif(x, 2))) +
    scale_y_continuous(name="Genomic span (Mb)", 
                       breaks = y_breaks,
                       labels = y_breaks/10^6, expand=c(0,0),
                       sec.axis = sec_axis(name = "Genomic span (%)", 
                                           trans = ~./sum(k11c_dt$tot_hg19),
                                           breaks = y_breaks_sec,
                                           labels = signif(y_breaks_sec*100,2)))+ #,labels = function(y) format(round(y*100, 2), nsmall=2))
    scale_fill_gradient2(low = "blue", high = "brown", midpoint = log10(5.96)) +
    #labs(title=paste("Signature",str_sub(dt_aggr$sig[1], 20, -3))) +
    coord_cartesian(xlim=range(x_breaks), ylim=c(0,max(y_breaks))) +
    labs(subtitle=paste0(signature_var," ",feature_var, "  (",round(mrate_weightedmean, 1)," SNV/patient/Mb)")) +
    theme_pubr(legend = "none") +
    theme(axis.title.y.right = element_text(angle=90))
  return(p_histo)
}

plot_pipeline <- function(signature, recurrence_regex="(5|6|7+)"){
  x <- subset_enchmmcount(signature=signature, recurrence_regex=recurrence_regex)

  match_enchmm_text <- "(?<=[:digit:]{1,2}_)[:alpha:]*[:punct:]*[:alpha:]*"
  x[, feature_collapse:=str_extract(feature, match_enchmm_text)]
  x2 <- x[, .(n_snv=sum(N), donors=paste(Donor_ID,collapse=",")), by=c("kmer_11","feature_collapse", "signature")]
  
  
  y <- merge(x2, kmer_dict_long_collapse, by=c("kmer_11","feature_collapse"))
  y[, n_donor:=uniqueN(x$Donor_ID)]
  y[, mrate:=n_snv/refcount_feature/n_donor*10^6]
  y[, mrate_log10:=round(log10(mrate),1)]
  
  range_features <- unique(y$feature_collapse)[order(factor(unique(y$feature_collapse), 
                                                            levels=feature_collapse_var_ord))]
  p_list <- list()
  y_gw <- y[, .(n_snv=sum(n_snv), feature_collapse="Genome-wide"), by=c("kmer_11","tot_hg19","n_donor")]
  y_gw[, `:=`(refcount_feature=tot_hg19, mrate=n_snv/tot_hg19/n_donor*10^6)]
  y_gw[, mrate_log10:=round(log10(mrate),1)]
  suppressWarnings({
    p_list[["logo_genomewide"]] <- plot_logo_bits_hg19auto_v2(y_gw)
    p_list[["mrate_genomewide"]] <- plot_mutrate_aggregated_by_span(y_gw)
    })
  
  for(ifeat in range_features){
    mrate_w <- y[feature_collapse==ifeat, sum(n_snv)/sum(refcount_feature)/n_donor*10^6]
    
    suppressWarnings({
      p_list[[paste0("logo_",ifeat)]] <- plot_logo_bits_hg19auto_v2(y[feature_collapse==ifeat])
      p_list[[paste0("mrate_",ifeat)]] <- plot_mutrate_aggregated_by_span(y[feature_collapse==ifeat])
    })
    }
  
  p_align <- Reduce("+",p_list) + plot_layout(ncol=2, byrow = TRUE, widths = c(1,3))
  return(p_align)
}

plot_pipeline2 <- function(signature, recurrence_regex="(5|6|7+)"){
  y <- subset_enchmmcount2(signature=signature, recurrence_regex=recurrence_regex)
  
  range_features <- unique(y$feature_collapse)[order(factor(unique(y$feature_collapse), 
                                                            levels=feature_collapse_var_ord))]
  p_list <- list()
  
  y_gw <- unique(y[, .(n_snv=sum(n_snv), tot_hg19=tot_hg19[!is.na(tot_hg19)], 
                       feature_collapse="Genome-wide"), by = c("kmer_11","n_donor")])
  y_gw[, `:=`(refcount_feature=tot_hg19, mrate=n_snv/tot_hg19/n_donor*10^6)]
  y_gw[, mrate_log10:=round(log10(mrate),1)]
  suppressWarnings({
    p_list[["logo_genomewide"]] <- plot_logo_bits_hg19auto_v2(y_gw)
    p_list[["mrate_genomewide"]] <- plot_mutrate_aggregated_by_span(y_gw)
  })
  
  for(ifeat in range_features){
    mrate_w <- y[feature_collapse==ifeat, sum(n_snv)/sum(refcount_feature)/n_donor*10^6]
    
    suppressWarnings({
      p_list[[paste0("logo_",ifeat)]] <- plot_logo_bits_hg19auto_v2(y[feature_collapse==ifeat])
      p_list[[paste0("mrate_",ifeat)]] <- plot_mutrate_aggregated_by_span(y[feature_collapse==ifeat])
    })
  }
  
  p_align <- Reduce("+",p_list) + plot_layout(ncol=2, byrow = TRUE, widths = c(1,3))
  return(p_align)
}

select_sig <- kmerstats[baselinesig==maxsig_kmer, unique(baselinesig)]

#x <- subset_enchmmcount2(select_sig[1], recurrence_regex="(1|2|3|4|5|6|7+)")

cat(paste0("Apply plotting pipeline\t",Sys.time(),"\n"))
for(i_signature in select_sig){
  cat(paste0("\t",i_signature,"\t", Sys.time(),"\n"))
  rec <- "(5|6|7+)"
  p_align <- plot_pipeline2(signature=i_signature, recurrence_regex=rec)
  filename <- paste0("plots/genomicregions_stratify_11mers/enchmm_stratified_11mers_rec5plus/logo_mratehisto_",
                     str_sub(i_signature,20,-3), "_rec=",rec,".pdf")
  ggsave(plot=p_align, device="pdf",width=10,height=20,file=filename)
  cat(paste0("Saved\t",filename,"\n",Sys.time(),"\n"))
}


# Specialized plotting for main figure:

# signatures  regions
# 17b         gw, strongenh, ins, rep, hetero
# 7a          gw, activeprom, strongenh, (ins, weaktxn), rep, hetero 
# 62          
# 72

plot_pipeline_mainfigure <- function(signature, recurrence_regex="(5|6|7+)", features_select=NULL){
  x <- subset_enchmmcount(signature=signature, recurrence_regex=recurrence_regex)
  
  match_enchmm_text <- "(?<=[:digit:]{1,2}_)[:alpha:]*[:punct:]*[:alpha:]*"
  x[, feature_collapse:=str_extract(feature, match_enchmm_text)]
  x2 <- x[, .(n_snv=sum(N), donors=paste(Donor_ID,collapse=",")), by=c("kmer_11","feature_collapse", "signature")]
  
  
  y <- merge(x2, kmer_dict_long_collapse, by=c("kmer_11","feature_collapse"))
  y[, n_donor:=uniqueN(x$Donor_ID)]
  y[, mrate:=n_snv/refcount_feature/n_donor*10^6]
  y[, mrate_log10:=round(log10(mrate),1)]
  
  p_list <- list()
  
  theme_rm_axis <- theme(axis.text.x = element_blank(), axis.title = element_blank())
  # genome-wide on top
  y_gw <- y[, .(n_snv=sum(n_snv), feature_collapse="Genome-wide"), by=c("kmer_11","tot_hg19","n_donor")]
  y_gw[, `:=`(refcount_feature=tot_hg19, mrate=n_snv/tot_hg19/n_donor*10^6)]
  y_gw[, mrate_log10:=round(log10(mrate),1)]
  suppressWarnings({
    p_list[["logo_genomewide"]] <- plot_logo_bits_hg19auto_v2(y_gw) + theme_rm_axis
    p_list[["mrate_genomewide"]] <- plot_mutrate_aggregated_by_span(y_gw) + theme_rm_axis
  })
  
  if(is.null(features_select)) {
    features_select <- unique(y$feature_collapse)[order(factor(unique(y$feature_collapse),levels=feature_collapse_var_ord))]
    }
  # middle features
  for(ifeat in head(features_select,-1)){
    
    mrate_w <- y[feature_collapse==ifeat, sum(n_snv)/sum(refcount_feature)/n_donor*10^6]
    
    suppressWarnings({
      p_list[[paste0("logo_",ifeat)]] <- plot_logo_bits_hg19auto_v2(y[feature_collapse==ifeat]) + theme_rm_axis
      p_list[[paste0("mrate_",ifeat)]] <- plot_mutrate_aggregated_by_span(y[feature_collapse==ifeat]) + theme_rm_axis
    })
  }
  # bottom feature with x-axis text
  bottom_feature <- tail(features_select, 1)
  suppressWarnings({
    p_list[[paste0("logo_",bottom_feature)]] <- plot_logo_bits_hg19auto_v2(y[feature_collapse==bottom_feature])
    p_list[[paste0("mrate_",bottom_feature)]] <- plot_mutrate_aggregated_by_span(y[feature_collapse==bottom_feature])
  })
  
  p_align <- Reduce("+",p_list) + plot_layout(ncol=2, byrow = TRUE, widths = c(1,3))
  return(p_align)
}

plot_pipeline2_mainfigure <- function(signature, recurrence_regex="(5|6|7+)", features_select=NULL){
  y <- subset_enchmmcount2(signature=signature, recurrence_regex=recurrence_regex)
  
  range_features <- unique(y$feature_collapse)[order(factor(unique(y$feature_collapse), 
                                                            levels=feature_collapse_var_ord))]
  p_list <- list()
  
  y_gw <- unique(y[, .(n_snv=sum(n_snv), tot_hg19=tot_hg19[!is.na(tot_hg19)], 
                       feature_collapse="Genome-wide"), by = c("kmer_11","n_donor")])
  y_gw[, `:=`(refcount_feature=tot_hg19, mrate=n_snv/tot_hg19/n_donor*10^6)]
  y_gw[, mrate_log10:=round(log10(mrate),1)]
  
  p_list <- list()
  
  theme_rm_axis <- theme(axis.text.x = element_blank(), axis.title = element_blank())
  # genome-wide on top
  y_gw <- y[, .(n_snv=sum(n_snv), feature_collapse="Genome-wide"), by=c("kmer_11","tot_hg19","n_donor")]
  y_gw[, `:=`(refcount_feature=tot_hg19, mrate=n_snv/tot_hg19/n_donor*10^6)]
  y_gw[, mrate_log10:=round(log10(mrate),1)]
  suppressWarnings({
    p_list[["logo_genomewide"]] <- plot_logo_bits_hg19auto_v2(y_gw) + theme_rm_axis
    p_list[["mrate_genomewide"]] <- plot_mutrate_aggregated_by_span(y_gw) + theme_rm_axis
  })
  
  if(is.null(features_select)) {
    features_select <- unique(y$feature_collapse)[order(factor(unique(y$feature_collapse),levels=feature_collapse_var_ord))]
  }
  # middle features
  for(ifeat in head(features_select,-1)){
    
    mrate_w <- y[feature_collapse==ifeat, sum(n_snv)/sum(refcount_feature)/n_donor*10^6]
    
    suppressWarnings({
      p_list[[paste0("logo_",ifeat)]] <- plot_logo_bits_hg19auto_v2(y[feature_collapse==ifeat]) + theme_rm_axis
      p_list[[paste0("mrate_",ifeat)]] <- plot_mutrate_aggregated_by_span(y[feature_collapse==ifeat]) + theme_rm_axis
    })
  }
  # bottom feature with x-axis text
  bottom_feature <- tail(features_select, 1)
  suppressWarnings({
    p_list[[paste0("logo_",bottom_feature)]] <- plot_logo_bits_hg19auto_v2(y[feature_collapse==bottom_feature])
    p_list[[paste0("mrate_",bottom_feature)]] <- plot_mutrate_aggregated_by_span(y[feature_collapse==bottom_feature])
  })
  
  p_align <- Reduce("+",p_list) + plot_layout(ncol=2, byrow = TRUE, widths = c(1,3))
  return(p_align)
}

signature_feature_list <- list(
  "17b_P"=list("features"=c("Strong_Enhancer", "Insulator", "Repetitive/CNV", "Heterochrom/lo"), "rec"="(5|6|7+)"),
  "7a_S"=list("features"=c("Active_Promoter", "Strong_Enhancer", "Insulator", "Repetitive/CNV", "Heterochrom/lo"), "rec"="(5|6|7+)"),
  "62_S"=list("features"=c("Poised_Promoter","Strong_Enhancer","Insulator","Txn_Transition","Repetitive/CNV", "Heterochrom/lo"), "rec"="(5|6|7+)"),
  "72_P"=list("features"=c("Active_Promoter","Weak_Promoter","Poised_Promoter", "Strong_Enhancer","Repressed","Heterochrom/lo"), "rec"="(5|6|7+)"))

save_to_folder <- "plots/genomicregions_stratify_11mers/enchmm_stratified_11mers_rec5plus/mainfig/"

for(i_sig in names(signature_feature_list)){
  sig_var <- i_sig
  rec_var <- signature_feature_list[[i_sig]]$rec
  fea_var <- signature_feature_list[[i_sig]]$features
  #cat(paste("\n",i_sig))
  #cat(paste("\n",rec_var))
  #cat(paste("\n",fea_var))
  
  p_align <- plot_pipeline2_mainfigure(
    signature = paste0("BI_COMPOSITE_SNV_SBS", i_sig),
    recurrence_regex = rec_var,
    features_select = fea_var)
  
  height_var <- length(fea_var)+1
  filename_var <- paste0("mrate_logo_S", str_sub(sig_var, 1,-3), "_rec",rec_var,"_selectfeatures.pdf")
  ggsave(plot = p_align, device = "pdf", width = 12, height = height_var,
         filename = paste0(save_to_folder, filename_var))
}
