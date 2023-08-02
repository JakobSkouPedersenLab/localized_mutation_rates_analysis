# created:    14-06-2021
# author:     Gustav Poulsgaard
# purpose:    generate overview insights from all previous analyses
#             and plot these insights.

# define environment ====
options(scipen = 999)
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

suppressPackageStartupMessages({
library(data.table) ; library(tidyverse) # data manipulation
library(ggpubr) ; library(patchwork) ; library(ggseqlogo) # plotting
library(Biostrings) ; library(BSgenome.Hsapiens.UCSC.hg19) # sequence manipulation
library(scales)}) # easy number formatting


# load data ====
cat(paste0("Load data\t",Sys.time(),"\n"))
donor_kmer_feature_dt <- fread("results/genomicregions_stratify_11mers/allgenomes_allkmers_overlap_enchmmfeatures.tsv")
kmerstats <- fread(file = "results/signatures_stratify_genomes_and_11mers/kmer_stats_stratified_by_signature.tsv")

k11c_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")

kmer_dict <- fread(file = "results/kmer_dict/kmer_background_countfeature_revcomplcollapse.tsv")
setnames(kmer_dict, old = c("count_hg19"), new = c("tot_hg19"))

cat(paste0("Compute statistics\t",Sys.time(),"\n"))
features_of_interest <- colnames(kmer_dict)[3:ncol(kmer_dict)]
kmer_dict[, count_hg19 := sum(.SD), by = seq, .SDcols = features_of_interest]
kmer_dict_long <- melt.data.table(data = kmer_dict, id.vars = c("seq", "tot_hg19"),
                                  measure.vars = colnames(kmer_dict)[3:17], variable.name = "feature", 
                                  value.name = "refcount_feature")[order(seq, feature)]
match_enchmm_text <- "(?<=[:digit:]{1,2}_)[:alpha:]*[:punct:]*[:alpha:]*"
kmer_dict_long[, `:=`(kmer_11=seq, feature_collapse=str_extract(feature, match_enchmm_text))]
kmer_dict_long_collapse <- kmer_dict_long[, .(refcount_feature=sum(refcount_feature)), by=c("kmer_11","feature_collapse", "tot_hg19")]
rm(kmer_dict_long)


# compute hg19 base distribution to use as background for logos
# quick and dirty:
b <- Reduce("+",lapply(1:22, function(chr) oligonucleotideFrequency(Hsapiens[[chr]], width = 1)))
pfm_exp <- matrix( (b/sum(b)), nrow = 4, ncol = 11, dimnames = list(names(b)))
rm(b)
# timely (>5 min):
# repseqs_all <- k11c_dt[, rep(kmer_11, tot_hg19)] # 37 sec
# #rid <- unique(round(runif(n = , min = 1, max = k11c_dt[, sum(tot_hg19)]),0))
# #repseqs <- repseqs_all[rid]
# id_intervals <- round(seq.int(1, length(repseqs_all), length.out = 51),0)
# # cat(paste0("|", paste0(rep("-",50),collapse=""), "|\n"))
# pfm_list <- lapply(1:5, function(i){
#          cat("=") ; i_start <- id_intervals[i] ; i_end <- id_intervals[i+1]
#          pfm <- consensusMatrix(x = repseqs_all[i_start:i_end]) #, width=11
#          return(pfm)})
# pfm_exp <- Reduce("+", pfm_list)

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


wpfm <- function(dt) consensusMatrix(DNAStringSet(dt[, rep(kmer_11, n_snv)]))[1:4,]
logo_list <- list()

signature_var <- "BI_COMPOSITE_SNV_SBS3_P"
rec_var <- "(5|6|7+)"
region_var <- "Repetitive/CNV"
dt <- data.table(step = 1:4, n_snv = 0, tot_hg19 = 0, n_kmer = 0, wmrate = 0, 
                 overall_wmrate = k11c_dt[, sum(n_snv)/sum(tot_hg19)/2583*10^6],
                 subset = c("active signature in genomes", "signature-assigned 11-mers",
                            "hotspot-annotated 11-mers", "regional instances of 11-mers"),
                 activity_sig = signature_var, assigned_sig = "none", rec = "1+", region = "Genome-wide")
# step 1
dt1 <- kmerstats[baselinesig==signature_var]
donor_var <- dt1[, donors] %>% str_c(collapse = ",") %>% str_split(",") %>% unlist() %>% unique()
values_var <- as.numeric(dt1[, .(sum(n_snv),sum(tot_hg19),uniqueN(kmer_11),sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6)])
dt[1, `:=`(n_snv=values_var[1],tot_hg19=values_var[2],n_kmer=values_var[3],wmrate=values_var[4])]
logo_list[[dt[1,subset]]] <- wpfm(dt1) # logo

# step 2
dt2 <- kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer]
values_var <- as.numeric(dt2[,.(sum(n_snv),sum(tot_hg19),uniqueN(kmer_11),sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6)])
dt[2, `:=`(n_snv=values_var[1],tot_hg19=values_var[2],n_kmer=values_var[3],wmrate=values_var[4],
           assigned_sig = signature_var)]
logo_list[[dt[2,subset]]] <- wpfm(dt2)
# step 3
dt3 <- kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer & grepl(rec_var, rec)]
values_var <- as.numeric(dt3[,.(sum(n_snv),sum(tot_hg19),uniqueN(kmer_11),sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6)])
dt[3, `:=`(n_snv=values_var[1],tot_hg19=values_var[2],n_kmer=values_var[3],wmrate=values_var[4],
           assigned_sig = signature_var, rec = rec_var)]
logo_list[[dt[3,subset]]] <- wpfm(dt3)
# step 4
subset_enchmmcount2 <- function(signature, recurrence_regex){
  # donor criteria
  donor_var <- kmerstats[#grepl(recurrence_regex, rec) & 
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
  
  y <- merge(x2, kmer_dict_long_collapse, by=c("kmer_11","feature_collapse"), all.x=TRUE)
  y[, n_donor:=uniqueN(x$Donor_ID)]
  y[, mrate:=n_snv/refcount_feature/n_donor*10^6]
  y[, mrate_log10:=round(log10(mrate),1)]
  
  return(y)
}
x <- subset_enchmmcount2(signature_var, recurrence_regex=rec_var)
#unique(x[,.(n_snv=sum(n_snv), tot_hg19=tot_hg19[!is.na(tot_hg19)]), by = kmer_11])[,.(sum(n_snv),sum(tot_hg19),sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6)]

# kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer & grepl("(5|6|7+)", rec),
#           .(sum(n_snv),sum(tot_hg19),sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6)]
region_var <- as.character(x[, sum(n_snv)/sum(refcount_feature), by = feature_collapse][which.max(V1), 1])
dt4 <- x[feature_collapse==region_var]
values_var <- as.numeric(dt4[,.(sum(n_snv), sum(refcount_feature),uniqueN(kmer_11), sum(n_snv)/length(donor_var)/sum(refcount_feature)*10^6)])
dt[4, `:=`(n_snv=values_var[1],tot_hg19=values_var[2],n_kmer=values_var[3],wmrate=values_var[4],
           assigned_sig = signature_var, rec = rec_var, region = region_var)]
logo_list[[dt[4,subset]]] <- wpfm(dt4)




dt[, fc:=wmrate/overall_wmrate]
dt[, fraction_gw:=tot_hg19/sum(k11c_dt$tot_hg19)]
fixround <- function(x, digits) gsub(" ","",format(round(x, digits), nsmall = digits))
dt[, mrate_label:=paste0(fixround(wmrate,2),"\n(",fixround(fc,2),"x)")]
dt[, gspan_label:=paste0(prettyNum(signif(tot_hg19/10^6,3), big.mark = ",")," Mb\n(",signif(fraction_gw*100,4),"%)")]

p_wmrate <- dt %>% ggplot(aes(x = factor(step,rev(step)), y = wmrate, fill=log10(wmrate))) +
  geom_col(col="black", width=0.75) + geom_text(aes(y=mean(range(wmrate)),label=mrate_label)) +
  geom_hline(yintercept=unique(dt[,overall_wmrate]), linetype="dashed", col="grey") +
  scale_x_discrete(name="",expand=c(0,0))+
  scale_y_continuous(name="Mutation rate (SNV/patient/Mb)",expand=c(0,0), position = "right")+
  scale_fill_gradient2(low="blue", mid="white",high="brown", midpoint=log10(unique(dt[,overall_wmrate]))) +
  coord_flip() + theme_pubr(legend="none") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y=element_blank())

p_logo <- ggplot() + geom_logo(data = logo_list) +
  scale_x_continuous(name="Position (relative to mutation)",
                     breaks=1:11,labels = function(x) x-6, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(title = paste0("Signature ", str_sub(signature_var, 21, -3), " rec", rec_var, " ",region_var)) +
  facet_wrap(~seq_group, ncol = 1) + 
  coord_cartesian(clip = "off")+
  theme_logo() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.y = element_line(), axis.ticks.y = element_line(),
        axis.title = element_text(size=12), axis.text = element_text(size=12), 
        strip.text = element_text(size=14))

p_gspan <- dt %>% ggplot(aes(x = factor(step,rev(step)), y = tot_hg19, fill=log10(wmrate))) +
  geom_col(col="black", width=0.75) + geom_text(aes(y=15*10^8,label=gspan_label)) +
  geom_hline(yintercept=sum(k11c_dt$tot_hg19), linetype="dashed", col="grey") +
  scale_x_discrete(name="",expand=c(0,0))+
  scale_y_continuous(name="Genomic span (Mb)", breaks = 0:2*10^9, position = "right",
                     expand=c(0,0), labels = function(x) comma(x/10^6)) +
  scale_fill_gradient2(low="blue", mid="white",high="brown", midpoint=log10(unique(dt[,overall_wmrate]))) +
  coord_flip() + theme_pubr(legend="none") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y=element_blank())


p_align <- p_wmrate + p_logo + p_gspan + plot_layout(nrow = 1, widths = c(1, 2, 1))

ggsave(plot = p_align, device = "pdf", width = 10, height = 5,
       filename = "plots/factorization_steps/v3/test.pdf")


# generalization now: ====
fixround <- function(x, digits) gsub(" ","",format(round(x, digits), nsmall = digits, big.mark = ","))
wpfm <- function(dt) consensusMatrix(DNAStringSet(dt[, rep(kmer_11, n_snv)]))[1:4,]
subset_enchmmcount2 <- function(signature, recurrence_regex){
  # donor criteria
  donor_var <- kmerstats[#grepl(recurrence_regex, rec) & 
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
  
  y <- merge(x2, kmer_dict_long_collapse, by=c("kmer_11","feature_collapse"), all.x=TRUE)
  y[, n_donor:=uniqueN(x$Donor_ID)]
  y[, mrate:=n_snv/refcount_feature/n_donor*10^6]
  y[, mrate_log10:=round(log10(mrate),1)]
  
  return(y)
}

##> results ====
results <- function(signature_var, rec_var, region_var=NULL){
  # template data.table
  dt <- data.table(step = 1:4, n_snv = 0, tot_hg19 = 0, n_kmer = 0, wmrate = 0, 
                   overall_wmrate = k11c_dt[, sum(n_snv)/sum(tot_hg19)/2583*10^6],
                   subset = c("active signature in genomes", "signature-assigned 11-mers",
                              "hotspot-annotated 11-mers", "regional instances of 11-mers"),
                   activity_sig = signature_var, assigned_sig = "none", rec = "1+", region = "Genome-wide")
  
  logo_list <- list()
  # step 1
  dt1 <- kmerstats[baselinesig==signature_var]
  donor_var <- dt1[, donors] %>% str_c(collapse = ",") %>% str_split(",") %>% unlist() %>% unique()
  values_var <- as.numeric(dt1[, .(sum(n_snv),sum(tot_hg19),uniqueN(kmer_11),sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6)])
  dt[1, `:=`(n_snv=values_var[1],tot_hg19=values_var[2],n_kmer=values_var[3],wmrate=values_var[4])]
  logo_list[[dt[1,subset]]] <- wpfm(dt1)
  
  # step 2
  dt2 <- kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer]
  values_var <- as.numeric(dt2[,.(sum(n_snv),sum(tot_hg19),uniqueN(kmer_11),sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6)])
  dt[2, `:=`(n_snv=values_var[1],tot_hg19=values_var[2],n_kmer=values_var[3],wmrate=values_var[4],
             assigned_sig = signature_var)]
  logo_list[[dt[2,subset]]] <- wpfm(dt2)
  
  # step 3
  dt3 <- kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer & grepl(rec_var, rec)]
  values_var <- as.numeric(dt3[,.(sum(n_snv),sum(tot_hg19),uniqueN(kmer_11),sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6)])
  dt[3, `:=`(n_snv=values_var[1],tot_hg19=values_var[2],n_kmer=values_var[3],wmrate=values_var[4],
             assigned_sig = signature_var, rec = rec_var)]
  logo_list[[dt[3,subset]]] <- wpfm(dt3)
  
  # step 4
  x <- subset_enchmmcount2(signature_var, recurrence_regex=rec_var)
  if(is.null(region_var)) region_var <- as.character(x[, sum(n_snv)/sum(refcount_feature), by = feature_collapse][which.max(V1), 1])
  dt4 <- x[feature_collapse==region_var]
  values_var <- as.numeric(dt4[,.(sum(n_snv), sum(refcount_feature),uniqueN(kmer_11), sum(n_snv)/length(donor_var)/sum(refcount_feature)*10^6)])
  dt[4, `:=`(n_snv=values_var[1],tot_hg19=values_var[2],n_kmer=values_var[3],wmrate=values_var[4],
             assigned_sig = signature_var, rec = rec_var, region = region_var)]
  logo_list[[dt[4,subset]]] <- wpfm(dt4)
  
  # annotations to plots
  dt[, fc:=wmrate/overall_wmrate]
  dt[, fraction_gw:=tot_hg19/sum(k11c_dt$tot_hg19)]
  dt[, mrate_label:=paste0(fixround(wmrate,2),"\n(",fixround(fc,2),"x)")]
  dt[, gspan_label:=paste0(prettyNum(signif(tot_hg19/10^6,3), big.mark = ",")," Mb\n(",signif(fraction_gw*100,4),"%)")]
  
  # combine results (data.table and logo list) in a list
  data_list <- list("dt"=dt, "pfm"=logo_list)
  return(data_list)
}

r_s3 <- results(signature_var = "BI_COMPOSITE_SNV_SBS3_P", "(5|6|7+)")

r_s7a <- results(signature_var = "BI_COMPOSITE_SNV_SBS7a_S", "(5|6|7+)")

p <- r_s7a[[1]] %>%
  ggplot(aes(x = tot_hg19, y = wmrate)) +
  geom_point() +
  geom_line(size=1) +
  scale_x_reverse()
ggsave(p, file="plots/factorization_steps/test/test1.pdf", device="pdf")


##> plots ====

plot_overview <- function(results_out){
  # make the list's entries easily accessible
  dt <- results_out$dt ; logo_list <- results_out$pfm
  # define signature, recurrence, and region
  signature_var <- dt$activity_sig[1] ; rec_var <- dt$rec[3] ; region_var <- dt$region[4]
  
  p_wmrate <- dt %>% ggplot(aes(x = factor(step,rev(step)), y = wmrate, fill=log10(wmrate))) +
    geom_col(col="black", width=0.75) + 
    geom_text(aes(y=mean(range(wmrate)),label=mrate_label, col=as.character(step))) +
    geom_hline(yintercept=unique(dt[,overall_wmrate]), linetype="dashed", col="grey") +
    scale_x_discrete(name="",expand=c(0,0))+
    scale_y_continuous(name="Mutation rate (SNV/patient/Mb)",expand=c(0,0), 
                       position = "right", labels = comma)+
    scale_color_manual(values = setNames(c(rep("black",3),"white"), 1:4)) +
    scale_fill_gradient2(low="blue", mid="white",high="brown", midpoint=log10(unique(dt[,overall_wmrate]))) +
    coord_flip() + theme_pubr(legend="none") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.text.y=element_blank())
  
  p_logo <- ggplot() + geom_logo(data = logo_list) +
    scale_x_continuous(name="Position (relative to mutation)",
                       breaks=1:11,labels = function(x) x-6, expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(title = paste0("Signature ", str_sub(signature_var, 21, -3), " rec", rec_var, " ",region_var)) +
    facet_wrap(~seq_group, ncol = 1) + 
    coord_cartesian(clip = "off", ylim = c(0,2))+
    theme_logo() +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.y = element_line(), axis.ticks.y = element_line(),
          axis.title = element_text(size=12), axis.text = element_text(size=12), 
          strip.text = element_text(size=14))
  
  p_gspan <- dt %>% ggplot(aes(x = factor(step,rev(step)), y = tot_hg19, fill=log10(wmrate))) +
    geom_col(col="black", width=0.75) + geom_text(aes(y=15*10^8,label=gspan_label)) +
    geom_hline(yintercept=sum(k11c_dt$tot_hg19), linetype="dashed", col="grey") +
    scale_x_discrete(name="",expand=c(0,0))+
    scale_y_continuous(name="Genomic span (Mb)", breaks = 0:2*10^9, position = "right",
                       expand=c(0,0), labels = function(x) comma(x/10^6)) +
    scale_fill_gradient2(low="blue", mid="white",high="brown", midpoint=log10(unique(dt[,overall_wmrate]))) +
    coord_flip() + theme_pubr(legend="none") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.text.y=element_blank())
  
  
  p_align <- p_wmrate + p_logo + p_gspan + plot_layout(nrow = 1, widths = c(1, 2, 1))
  
  
}

p_s7a <- plot_overview(r_s7a)

ggsave(plot = p_s7a, device = "pdf", width = 12, height = 5,
       filename = "plots/factorization_steps/v3/test.pdf")


signatures <- kmerstats[,unique(baselinesig)]

lapply(signatures[4:57], function(i_sig){
  r <- results(i_sig,rec_var = "(5|6|7+)") ; p <- plot_overview(r) 
  filename_var <- paste0("sig", str_sub(r$dt$activity_sig[1],20,-3),"_rec", r$dt$rec[3], 
                         "_", gsub("[[:punct:]]", "",r$dt$region[4]), ".pdf")
  ggsave(plot=p, device="pdf", width=12, height=5,
         filename=paste0("plots/factorization_steps/v3/", filename_var))
  })

