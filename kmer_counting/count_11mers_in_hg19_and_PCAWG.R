# Date created:   25-05-2021
# Author:         Gustav Poulsgaard
# Purpose:        Generate file of 11-mer annotations
#                 11-mer, count in hg19, mutation count in PCAWG, recurrence counts


options(scipen=999)
setwd("PCAWG/faststorage/Gustav/mutprocesses_project")
library(data.table) # data storage and manipulation
library(Biostrings) # oligonucleotideFrequency(), reverseComplement(), DNAStringSet()
library(BSgenome) # 
library(BSgenome.Hsapiens.UCSC.hg19) # hg19 sequence (Hsapiens)
library(stringr) # str_sub()
library(tidyr) # spread()


# PCAWG 11-mer counting ====
# load data
#snv <- fread("data/snvdata_w_tfbs_hmm_overlap.tsv") # contain coding mutations
snv <- fread("data/signature_snv_filtered_annotated4.tsv", select=c("eid","n_pc","kmer_11"))
snv[, rec:=ifelse(n_pc>=7, "7+", n_pc)]
# snv <- fread("../data/signatures/signature_snv_autosome_n_project_header_filt_annot.tsv", 
#              select = c("eid", "rec","flankSeq_ag"))
# set(snv, NULL, "kmer_11", str_sub(snv$flankSeq_ag, 11-5, 11+5))
# set(snv, NULL, "flankSeq_ag", NULL)
# summarise SNV count by kmer and recurrence
kmer_rec <- snv[, .N, by = c("kmer_11","rec")]
# summarise SNV count and recurrence by kmer (paste0 with collapse is heavy; takes a couple of min to run)
kmer <- kmer_rec[, .(n_snv=sum(N), rec_overview=paste0(sort(rec), collapse=",")), by = kmer_11] 
setnames(kmer, old = "rec_overview", new = "rec")
# spread the counts by recurrence into each column from 1 through 7+
kmer_wide <- tidyr::spread(kmer_rec, key = "rec", value="N", fill=0, sep="")
setnames(kmer_wide, old = "rec7+", new = "rec7p")
# merge summaries by kmer
kmer2 <- merge(kmer, kmer_wide, by = "kmer_11")
setkeyv(kmer2, "kmer_11")

# hg19 11-mer counting ====
# count 11-mer occurrences in hg19 (chromosome 1-22)
ref_kmer <- Reduce("+", lapply(1:22, function(chr) oligonucleotideFrequency(Hsapiens[[chr]], width = 11)))
ref_kmer_dt <- data.table(kmer_11 = names(ref_kmer), tot_hg19 = as.numeric(ref_kmer), key = "kmer_11")
ref_kmer_dt[, center_nuc := str_sub(kmer_11,6,6)]
rc <- function(x) as.character(reverseComplement(DNAStringSet(x)))
ref_kmer_dt[!(center_nuc%in%c("C","T")), kmer_11 := rc(kmer_11)]
# summarise occurrence by center pyrimidine
ref_kmer_dt2 <- ref_kmer_dt[, .(tot_hg19 = sum(tot_hg19)), by = kmer_11]
setkeyv(ref_kmer_dt2, "kmer_11")

# merge PCAWG and hg19 11-mer counts ====
# merge reference occurrence with snv counts
kmer_dt <- merge(ref_kmer_dt2, kmer2, all = TRUE)
# replace NAs with zeros
for(colname in names(kmer_dt)) set(x=kmer_dt, i=which(is.na(kmer_dt[[colname]])), j=colname, value=0)
#set(kmer_dt, NULL, "rec2p", rowSums(kmer_dt[, c(paste0("rec",2:6),"rec7p")]))
#set(kmer_dt, NULL, "rec5p", rowSums(kmer_dt[, c("rec5","rec6","rec7p")]))
fwrite(kmer_dt, file = "results/11mer_count_hg19_PCAWG_v3.tsv", sep = "\t")
