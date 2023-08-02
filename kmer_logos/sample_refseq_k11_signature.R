# author: Gustav Alexander Poulsgaard
# date: 2022-02-23
# purpose: create background set of 11-mers to be used with kplogo

# define environment ====
library(data.table) # read, write, and manipulate data
# set working directory
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project/")

# load reference 11-mers ====
# read data
kmer11_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")
# select columns
k11 <- kmer11_dt[, .(kmer = kmer_11,
                     freq = tot_hg19/sum(tot_hg19),
                     triN = substr(kmer_11, 5,7))]
# set key
setkeyv(k11, "triN")

# load reference signatures ====
refsig <- fread("data/reference_signatures/SignatureAnalyzer_COMPOSITE_SBS_W96.signature.031918.txt")
# define trinucleotide
refsig[, triN := substr(feature, 8, 10)]
# store signature names
sbs_names <- grep("SBS", colnames(refsig), value = TRUE)
# summarize signature probabilities across trinucleotides
sigprob <- refsig[, lapply(.SD, sum), by = triN, .SDcols = sbs_names]
# set key
setkeyv(sigprob, "triN")

# join data ====
k11mer_freq_prob <- k11[sigprob]
k11mer_weights <- k11mer_freq_prob[, lapply(.SD, function(x) freq*x*10^6),
                                   .SDcols = sbs_names, by = kmer]

# sample background 11-mers ====
# set seed to the ISO date of script creation
set.seed(20220223)
# sample 11-mers
background_set <- k11mer_weights[, sample(x = kmer, size = 10^6, replace = TRUE,
                                          prob = BI_COMPOSITE_SNV_SBS17b_P)]
# check length of background set
length(background_set)


# save background set ====
fwrite(x = list(background_set),
       file = "results/kplogo/s17b_logo/input/backgroundset_s17b.tsv")
