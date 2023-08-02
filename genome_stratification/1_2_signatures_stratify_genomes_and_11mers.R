# Date created: 28-05-2021
# Authors:      Gustav Poulsgaard
# Purpose:      Derive data for plotting of signature-cohorts and signature-11mers


timefun <- function(ts,te,unit="sec",d=1) round(difftime(te,ts,units=unit),digits = d)

time_setenvironment <- Sys.time()
cat(paste0("\nSet environment\t", time_setenvironment, "\t"))

options(scipen=999)
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")
library(data.table, quietly = TRUE)
library(stringr, quietly = TRUE)

time_loaddata <- Sys.time()
cat(paste0(timefun(time_setenvironment, time_loaddata), " sec\n"))
cat(paste0("\nLoad data\t", time_loaddata, "\t"))

sigact <- fread("data/PCAWG_signatures/SA_COMPOSITE_SNV.activity_relative.WHITELIST_SET.20200226.txt")
sbs_columns <- grep("BI_COMPOSITE", colnames(sigact), value=TRUE)

kmer_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")
# exclude columns to load faster
column_names <- colnames(fread("data/signature_snv_filtered_annotated4.tsv", nrow = 0))
column_select <- column_names[!column_names%in%c("snpOverlap","byFrequency","sample","normid")]
snv <- fread("data/signature_snv_filtered_annotated4.tsv", select = column_select)
setnames(snv, old = "BI_COMPOSITENV_SBS72_P", new = "BI_COMPOSITE_SNV_SBS72_P")


donor_sigact <- rbindlist(sapply(X = sbs_columns,
       FUN = function(x) sigact[get(x)>=0.05, .(icgc_donor_id, Tumor_Type, active_sig=x)],
       simplify = FALSE))
#[, .(active_sigs = paste0(active_sig, collapse=","), n_sigact = .N), by = icgc_donor_id]

time_definepipeline <- Sys.time()
cat(paste0(timefun(time_loaddata, time_definepipeline, unit = "min"), " min\n"))
# define pipeline
cat(paste0("\nDefine pipeline\t", time_definepipeline, "\t"))


# step 1
filter_by_signatureactivity <- function(signature){
  donor_target <- donor_sigact[active_sig==signature, icgc_donor_id]
  snv_activityfilter <- snv[Donor_ID%in%donor_target]
  snv_activityfilter[, active_signature := signature]
  return(snv_activityfilter)
}
# step 2
compute_and_filterby_mean_kmer_signatureload <- function(dt){
  pos_mean <- dt[, lapply(.SD, mean), .SDcols = sbs_columns, by = c("eid", "kmer_11")]
  kmer_mean <- pos_mean[, lapply(.SD, mean), .SDcols = sbs_columns, by = kmer_11] # 7 sec
  kmer_long <- as.data.table(tidyr::gather(kmer_mean, key = "signature", value = "prob", sbs_columns))
  kmer_long[, maxprob:=max(prob), by=kmer_11]
  kmer_long_max <- kmer_long[prob==maxprob, .(kmer_11, maxsig_kmer=signature, maxprob_kmer=maxprob)]
  return(kmer_long_max)
}
# step3
annotate_kmers <- function(dt_snv, dt_kmersigsum){
  snv_kmer <- dt_snv[, .(n_snv=.N, donors=paste0(unique(Donor_ID),collapse=","), 
                         baselinesig=unique(active_signature)), by=kmer_11]
  x <- merge(snv_kmer, kmer_dt[,.(kmer_11, tot_hg19, rec)], by="kmer_11")
  y <- merge(x, dt_kmersigsum, by="kmer_11")
  return(y)
}
# step 1-3
pipeline <- function(signature, save=FALSE, verbose=FALSE){
  if(verbose) {t1<-Sys.time() ; cat(paste0("\tsignature:\t", signature, "\n\tstart:\t",t1,"\n\tstep1:\t"))}
  snvs_pts_with_active_signature <- filter_by_signatureactivity(signature)
  if(verbose) {t2<-Sys.time() ; cat(paste0(timefun(t1,t2,d=0), " sec\n\tstep2:\t"))}
  kmer_counts <- compute_and_filterby_mean_kmer_signatureload(snvs_pts_with_active_signature)
  if(verbose) {t3<-Sys.time() ; cat(paste0(timefun(t2,t3,d=0), " sec\n\tstep3:\t"))}
  kmers_annotated <- annotate_kmers(snvs_pts_with_active_signature, kmer_counts)
  if(verbose) {t4<-Sys.time() ; cat(paste0(timefun(t3,t4,d=0), " sec\n\ttotal:\t",timefun(t1,t4,unit="min")," min\n"))}
  if(save) {
    path <- "results/signatures_stratify_genomes_and_11mers/temp/"
    filename <- paste0("annotated_kmer_subset_",str_sub(signature,18,-3), ".tsv")
    fwrite(x=kmers_annotated, file=filename, sep="\t")
    cat(paste0("saved ", filename, " to ", path, "\n"))
    } else return(kmers_annotated)
}

# run pipeline
time_runpipeline <- Sys.time()
cat(paste0(timefun(time_definepipeline, time_runpipeline, unit = "sec"), " sec\n"))
cat(paste0("\nRun pipeline\t", time_runpipeline, "\t"))

lapply(sbs_columns, pipeline, verbose = TRUE)

time_end <- Sys.time()
cat(paste0(timefun(time_runpipeline, time_end)))
cat(paste0("\nResults saved\t", time_end, "\nTotal time elapsed:\t", 
           timefun(time_setenvironment, time_end, unit = "min"), " min\n\n"))


#kmer_sigs <- rbindlist(lapply(sbs_columns, pipeline, verbose = TRUE))

time_saveresults <- Sys.time()
cat(paste0(timefun(time_runpipeline, time_saveresults, unit = "min"), " min\n"))


file_list <- system("ls results/signatures_stratify_genomes_and_11mers/temp", intern = TRUE)

all_subset <- rbindlist(lapply(file_list, fread))
fwrite(x=all_subset, sep="\t", file="results/signatures_stratify_genomes_and_11mers/annotated_kmer_subsets.tsv")

#cat(paste0("\nSave results\t", time_saveresults, "\t"))
#cat(paste0("\n\tresults/signatures_stratify_genomes_and_11mers/kmer_stats_stratified_by_signature.tsv"))
#fwrite(kmer_sigs, sep="\t",
#       file="results/signatures_stratify_genomes_and_11mers/kmer_stats_stratified_by_signature.tsv")
#time_end <- Sys.time()
#cat(paste0("\nResults saved\t", time_end, "\nTotal time elapsed:\t", 
#           timefun(time_setenvironment, time_end, unit = "min"), " min\n\n"))
