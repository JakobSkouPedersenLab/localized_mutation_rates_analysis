
# number of available cores
n_cores <- 1 #as.numeric(commandArgs(trailingOnly=TRUE))


require(data.table, quietly = T, warn.conflicts = F)
require(stringr, quietly = T, warn.conflicts = F)
require(BSgenome.Hsapiens.UCSC.hg19, quietly = T, warn.conflicts = F)
require(Biostrings, quietly = T, warn.conflicts = F)
require(parallel, quietly = T, warn.conflicts = F)


# strategy for countning 11-mers inside genomic regions.

# (1) create list of subjects to search in
#     split human genome by the different regions-of-interest
#     - what about boundaries?
#       + include 5 bp extra at each boundary to catch an 11-mer
#           with center at the end of a region
#     exclude sex chromosomes
# (2) create list of patterns to search for (count)
# (3) count the presence of each pattern in each subject

setwd("~/PCAWG/faststorage/Gustav")

# load data
hmmtrack <- fread("reference_work/encodetracks/trimmed_tracks/wgEncodeBroadHmmH1hescHMM_trim.bed", 
        drop =  c(5, 6, 7, 8), col.names = c("chr", "start", "end", "feature", "rgb"))[!(chr %in% c("chrX","chrY"))]
setkeyv(hmmtrack, c("chr","start", "end"))

# order by number
features_of_interest <- unique(hmmtrack$feature)[order(as.numeric(str_extract(unique(hmmtrack$feature), "[:digit:]+")))]
# order by size
features_of_interest <- features_of_interest[c(15,14,4,3,1,9,5,12,10,2,8,6,11,7,13)]


# (1) ====
# split genome into regions-of-interest
hmmtrack[feature==features_of_interest[15], min(start)]

fetchRegionSeqs <- function(region){
  region_intervals <- hmmtrack[feature==region]
  n <- nrow(region_intervals)
  seqs <- sapply(1:n, function(i) {
    i_range <- region_intervals[i, c("chr","start", "end")]
    i_range[, `:=`(start_m5 = start-5, end_p5 = end+5)]
    i_range[, start_m5_min0 := ifelse(start_m5<0, 0, start_m5)]
    seq <- Hsapiens[[i_range$chr]][i_range$start_m5_min0:i_range$end_p5]
    return(seq)
    })
  return(seqs)
  }


# (2) ====
cat("\nCount ~4 million 11-mers\n")
t1 <- Sys.time()
gw_kmer_count <- Reduce("+", mclapply(1:22, function(chr) oligonucleotideFrequency(Hsapiens[[chr]], width = 11), mc.cores = n_cores))
dt <- data.table(seq = names(gw_kmer_count), total_count = gw_kmer_count)
rm(gw_kmer_count)
data.table::set(dt, j = features_of_interest, value = NA)
t2 <- Sys.time()
cat(paste0("\n(", as.character(round(difftime(t2, t1, units = "min"),0)), " min)\n" ))

# (3) ====

# iterate (1)-(3)

mclapply(1:14, function(i) {
  
  cat(paste("done:",paste0(features_of_interest[ as.logical(dt[, lapply(.SD, function(x) all(!is.na(x))), .SDcols = features_of_interest]) ] , 
  collapse = ", ")))
  
  t1 <- Sys.time()
  subject_set <- DNAStringSet( fetchRegionSeqs(region = features_of_interest[i]) )
  t2 <- Sys.time()
  set(dt, j=features_of_interest[i], value=vcountPDict(PDict(dt$seq), subject_set, collapse = 1))
  t3 <- Sys.time()
  elap1.2 <- as.character(round(difftime(t2, t1, units = "min"),0))
  elap2.3 <- as.character(round(difftime(t3, t2, units = "min"),0))
  if(i==1) cat("\n\titeration\t\t(step2 + step3)")
  cat(paste0("\n\t",i,"\n\t",features_of_interest[i],"\t(", elap1.2," + ",elap2.3, " min)\n\n"))
  },
  mc.cores = n_cores)

fwrite(dt, sep = "\t", file = "mutprocesses_project/results/kmer_dict/kmer_background_countfeature.tsv")
cat("\n\tfile:\t\t\tkmer_background_countfeature.tsv\n\tsuccessfully saved in:\tmutprocesses_project/results/kmer_dict\n")

dt <- fread("mutprocesses_project/results/kmer_dict/kmer_background_countfeature.tsv")


cat(paste0("\nHere goes the Heterochromatin job\nstart:\t\t", Sys.time()))

features_of_interest

t1 <- Sys.time()
subject_set_heterochrom <- DNAStringSet( fetchRegionSeqs(region = "13_Heterochrom/lo") )
t2 <- Sys.time()
elap1.2 <- as.character(round(difftime(t2, t1, units = "min"),0))
cat(paste0("\n fetchRegionSeqs from heterochromatin (nsubjects = ", 
           length(subject_set_heterochrom),")\n(", elap1.2," min)\n"))


irange <- round(seq(0, nrow(dt), length.out = 33), 0)

out_list <- mclapply(1:32, function(i) {
  t3 <- Sys.time()
  out <- data.table(seq = dt$seq[(irange[i]+1):irange[i+1]], 
    `13_Heterochrom/lo` = vcountPDict(PDict(dt$seq[(irange[i]+1):irange[i+1]]), subject_set_heterochrom, collapse = 1))
  t4 <- Sys.time()
  elap2.3 <- as.character(round(difftime(t4, t3, units = "min"),0))
  if(i==1) cat("\n\titeration\t\t(step2 + step3)")
  cat(paste0("\n\t", irange[i]+1, ":", irange[i+1],"\n\t(", elap2.3, " min)\n\n"))
  return(out)
},
mc.cores = 32)

length(out_list)
str(out_list)
dt_13 <- rbindlist(out_list)
setnames(dt_13, "count", "13_Heterochrom/lo")

dt1 <- cbind(dt[, -c("13_Heterochrom/lo")], dt_13[,-c("seq")])

fwrite(dt1, sep = "\t", file = "mutprocesses_project/results/kmer_dict/kmer_background_countfeature.tsv")
cat("\n\tfile successfully updated\t\t\tkmer_background_countfeature.tsv\n")
