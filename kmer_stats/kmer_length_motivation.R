# motivate the length of k-mers

# nullomers
# expected number of mutations in each k-mer
# potentially computational efficiency

# k = 1, 3, ..., 13

#
#

options(scipen=999)

library(data.table)
library(stringr)
library(scales)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)


dt <- data.table(k = seq(1,13, by = 2))
n_snv <- 41318050

kmerFreq <- function(klength){
  out <- Reduce("+",lapply(1:22, function(chr) oligonucleotideFrequency(x = Hsapiens[[chr]], width = klength)))
  return(out)
}
revcom_kmers <- function(x){
  fwd_seqs <- names(x)
  rev_seqs <- as.character(reverseComplement(DNAStringSet(fwd_seqs)))
  dt <- data.table(fwd = fwd_seqs, rev = rev_seqs, fwd_count = as.numeric(x))
  k <- str_length(dt$fwd[1])
  p <- (k+1)/2
  # dt[, Y_core := ifelse(str_sub(fwd, p,p)%in%c("C","T")]
  dt[, seq_Ycore := ifelse(str_sub(fwd, p, p)%in%c("C","T"), fwd, rev)]
  dt[, tot_count := sum(fwd_count), by = seq_Ycore]
  out <- dt[rev==seq_Ycore, .(k = k, seq = seq_Ycore, count = tot_count)]
  return(out)
}

dt <- rbindlist(lapply(dt$k, function(k) {
  cat(k)
  x <- kmerFreq(klength = k)
  cat("#")
  y <- revcom_kmers(x = x)
  return(y)
  }))

dim(dt)
# number of families
dt_sum <- dt[, .(n_fam_exp = .N,
             # number of hg19 11-mer families
             n_fam_obs = sum(count>0),
             # number of absent families (nullomers)
             n_fam_null = sum(count==0),
             # expected number of snvs per family
             n_snv_unif = n_snv/sum(count>0),
             # summarised instances
             n_inst_min = min(count),
             n_inst_q1 = quantile(count, .25),
             n_inst_median = median(count),
             n_inst_mean = mean(count),
             n_inst_q3 = quantile(count, .75),
             n_inst_max = max(count)),
         by = k]

fwrite(dt_sum, sep = "\t",
       file = "~/PCAWG/faststorage/Gustav/mutprocesses_project/results/kmer_length_stats.tsv")

dt[k==13&count==0, .(.N, consensusString(seq, threshold = 0.5)),
   by = str_sub(seq, 6,8)][order(-rank(N))]


dt[, .(max(count), seq[count==max(count)]), by = k]

