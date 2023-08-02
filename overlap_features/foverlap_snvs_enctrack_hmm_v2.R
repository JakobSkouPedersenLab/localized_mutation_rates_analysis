# Author          Gustav Poulsgaard
# Date created    03-06-2021
# Purpose         overlap enchmm features with SNVs
#                 and summarise results by donor id and sequence context

# set environment ====
options(scipen = 999)
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

# load packages ====
library(tidyverse)
library(data.table)
library(scales)

# load data ====
kmer_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")

## encode track file
hmmtrack <- fread("../reference_work/encodetracks/trimmed_tracks/wgEncodeBroadHmmH1hescHMM_trim.bed",
                  drop =  5:9, col.names = c("chr", "start", "end", "feature"))
setkey(hmmtrack, chr, start, end)

## snv data
### exclude columns to load faster
column_select <- c("Chromosome","End_position","Donor_ID","kmer_11")
snv <- fread("data/signature_snv_filtered_annotated4.tsv", select = column_select)
setnames(snv, c("Chromosome", "End_position"), c("chr", "end"))
snv[, start:=end] ; setkey(snv, chr, start, end)

# overlap data (20 sec) ====
fov <- foverlaps(x = snv, y = hmmtrack, type = "any", by.x = c("chr", "start", "end"), mult="first")

# summarise ====
donor_kmer_feature_dt <- fov[, .N, by=c("Donor_ID","kmer_11","feature")]

# save results ====
fwrite(donor_kmer_feature_dt, sep="\t", 
       file="results/genomicregions_stratify_11mers/allgenomes_allkmers_overlap_enchmmfeatures.tsv")

# compute fraction of genome taken up by each feature
# hmmtrack[, length:=end-start]
# hmmsum <- hmmtrack[, .(total_span=sum(length)), by=feature]
# hmmsum[, fraction_span:=total_span/sum(total_span)]
