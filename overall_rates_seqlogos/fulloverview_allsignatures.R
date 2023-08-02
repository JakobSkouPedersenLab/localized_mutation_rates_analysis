# created:    08-07-2021
# author:     Gustav Poulsgaard
# purpose:    generate full overview from all previous analyses
#             and plot in mega figure.

# define environment ====
options(scipen = 999)
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

suppressPackageStartupMessages({
  library(data.table) ; library(tidyverse) # data manipulation
  library(ggpubr) ; library(patchwork) ; library(ggseqlogo) # plotting
  library(Biostrings) ; library(BSgenome.Hsapiens.UCSC.hg19) # sequence manipulation
  library(scales)# easy number formatting
  }) 

# load data
k11c_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")

donor_kmer_feature_dt <- fread("results/genomicregions_stratify_11mers/allgenomes_allkmers_overlap_enchmmfeatures.tsv")
kmerstats <- fread("results/signatures_stratify_genomes_and_11mers/kmer_stats_stratified_by_signature.tsv")

kmer_dict <- fread(file = "results/kmer_dict/kmer_background_countfeature_revcomplcollapse.tsv")
setnames(kmer_dict, old = c("count_hg19"), new = c("tot_hg19"))



features_of_interest <- colnames(kmer_dict)[3:ncol(kmer_dict)]
kmer_dict[, count_hg19 := sum(.SD), by = seq, .SDcols = features_of_interest] # TODO: why so long runtime?
kmer_dict_long <- melt.data.table(data = kmer_dict, id.vars = c("seq", "tot_hg19"),
                                  measure.vars = colnames(kmer_dict)[3:17], variable.name = "feature", 
                                  value.name = "refcount_feature")[order(seq, feature)]
match_enchmm_text <- "(?<=[:digit:]{1,2}_)[:alpha:]*[:punct:]*[:alpha:]*"
kmer_dict_long[, `:=`(kmer_11=seq, feature_collapse=str_extract(feature, match_enchmm_text))]
kmer_dict_long_collapse <- kmer_dict_long[, .(refcount_feature=sum(refcount_feature)), by=c("kmer_11","feature_collapse", "tot_hg19")]
rm(kmer_dict_long, features_of_interest, kmer_dict)

# b <- Reduce("+",lapply(1:22, function(chr) oligonucleotideFrequency(Hsapiens[[chr]], width = 1)))
# pfm_exp <- matrix( (b/sum(b)), nrow = 4, ncol = 11, dimnames = list(names(b)))
# rm(b)


# convenient rounding and comma
fixround <- function(x, digits=1) gsub(" ","",format(round(x, digits), nsmall = digits, big.mark = ","))
# position frequency matrix
PFM <- function(seq) consensusMatrix(DNAStringSet(seq))[1:4,]
compressPFM <- function(pfm) paste(sapply(1:4, function(i) paste(c(rownames(pfm)[i],pfm[i,]), collapse="\t")), collapse = "\n")
#
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
# define results
racerbil <- function(signature_var){
  cat("#1")
  donor_var <- kmerstats[baselinesig==signature_var, donors] %>% str_c(collapse = ",") %>% str_split(",") %>% unlist() %>% unique()
  
  s1 <- kmerstats[baselinesig==signature_var,.(step=1,
                                               n_snv = sum(n_snv), 
                                         tot_hg19 = sum(tot_hg19),
                                         n_kmer = uniqueN(kmer_11),
                                         n_donor = length(donor_var),
                                         wmrate = sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6,
                                         wpfm = compressPFM(PFM(rep(kmer_11, n_snv))))]
  
  cat("#2")
  s2 <- kmerstats[baselinesig==signature_var & 
              baselinesig==maxsig_kmer,.(step=2,
                                         n_snv = sum(n_snv), 
                                         tot_hg19 = sum(tot_hg19),
                                         n_kmer = uniqueN(kmer_11),
                                         n_donor = length(donor_var),
                                         wmrate = sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6,
                                         wpfm = compressPFM(PFM(rep(kmer_11, n_snv))))]
  
  
  #m2 <- kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer & grepl("(2|3|4|5|6|7+)", rec), .(mrate=sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6, rec="(2|3|4|5|6|7+)")]
  #m5 <- kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer & grepl("(5|6|7+)", rec), .(mrate=sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6, rec="(5|6|7+)")]
  cat("#3a")
  s3_2p <- kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer & grepl("(2|3|4|5|6|7+)", rec), 
            .(step=3,
              n_snv = sum(n_snv), 
             tot_hg19 = sum(tot_hg19),
             n_kmer = uniqueN(kmer_11),
             n_donor = length(donor_var),
             wmrate = sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6,
             wpfm = compressPFM(PFM(rep(kmer_11, n_snv))),
             rec = "(2|3|4|5|6|7+)")]
  cat("b")
  s3_5p <- kmerstats[baselinesig==signature_var & baselinesig==maxsig_kmer & grepl("(5|6|7+)", rec), 
                     .(step=3,
                       n_snv = sum(n_snv), 
                       tot_hg19 = sum(tot_hg19),
                       n_kmer = uniqueN(kmer_11),
                       n_donor = length(donor_var),
                       wmrate = sum(n_snv)/length(donor_var)/sum(tot_hg19)*10^6,
                       wpfm = compressPFM(PFM(rep(kmer_11, n_snv))),
                       rec = "(5|6|7+)")]
  s3 <- rbind(s3_2p, s3_5p)
  rec_var <- s3[, rec[wmrate==max(wmrate)]]
  mrate_step3_var <- s3[, max(wmrate)]
  
  cat("#4")
  x <- subset_enchmmcount2(signature = signature_var, recurrence_regex = rec_var)
  s4 <- x[,.(step=4,
             n_snv = sum(n_snv), 
       tot_hg19 = sum(refcount_feature),
       n_kmer = uniqueN(kmer_11),
       n_donor = length(donor_var),
       wmrate = sum(n_snv)/length(donor_var)/sum(refcount_feature)*10^6,
       wpfm = compressPFM(PFM(rep(kmer_11, n_snv))),
       rec = rec_var), by = feature_collapse]#[wmrate>=mrate_step3_var,]
  s <- rbind(s1,s2,s3,s4, fill=TRUE)
  s[, signature:=signature_var]
  cat("\n")
  return(s)
}

#d17b <- racerbil("BI_COMPOSITE_SNV_SBS17b_P")
#d17b#[step!=4 | wmrate>=wmrate[step==3]]


library(parallel)
t0 <- Sys.time()
dt <- rbindlist( mclapply(unique(kmerstats$baselinesig), racerbil, mc.cores = 10) )
t1 <- Sys.time() ; difftime(t1, t0)

fwrite(dt, file="results/fulloverview_allsignatures.tsv", sep="\t")

# annotations to plots
dt[, fc:=wmrate/overall_wmrate]
dt[, fraction_gw:=tot_hg19/sum(k11c_dt$tot_hg19)]
dt[, mrate_label:=paste0(fixround(wmrate,2),"\n(",fixround(fc,2),"x)")]
dt[, gspan_label:=paste0(prettyNum(signif(tot_hg19/10^6,3), big.mark = ",")," Mb\n(",signif(fraction_gw*100,4),"%)")]

# combine results (data.table and logo list) in a list
data_list <- list("dt"=dt, "pfm"=logo_list)


dt <- fread("results/fulloverview_allsignatures.tsv")
read.table(text = compressPFM(logo_list[[1]]), row.names = 1)
