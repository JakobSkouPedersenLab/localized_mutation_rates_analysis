#!/usr/bin/env Rscript
# Author: Gustav Poulsgaard
# Date: 15-08-2022
# Purpose: produce 11-mer foreground sets, background set, and run kplogo

options(scipen = 999)
cat("##############################################\n")
cat(paste0("#### Start of script: ",Sys.time()," ####\n"))
cat("##############################################\n")
# define environment ====
library(data.table)
library(stringr)
library(bit64)
project_dir <- "/faststorage/project/PCAWG/Gustav/mutprocesses_project"
setwd(project_dir)

# load data required to generate background ====
if(!exists("k11")) {
  kmer11_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")
  # select columns
  k11 <- kmer11_dt[, .(kmer = kmer_11,
                       #count_11mer = tot_hg19,
                       freq_11mer = tot_hg19/sum(tot_hg19),
                       triN = substr(kmer_11, 5,7))][,
                                                     `:=`(#count_triN = sum(count_11mer),
                                                       freq_triN = sum(freq_11mer)), by = triN]
  # set key
  setkeyv(k11, "triN")
  rm(kmer11_dt)
}
if(!exists("sigprob") | !exists("sbs_names")) {
  refsig <- fread("data/reference_signatures/SignatureAnalyzer_COMPOSITE_SBS_W96.signature.031918.txt")
  # store signature names
  sbs_names <- grep("SBS", colnames(refsig), value = TRUE)
  # define trinucleotide
  refsig[, triN := substr(feature, 8, 10)]
  # summarize signature probabilities across trinucleotides
  sigprob <- refsig[, lapply(.SD, sum), by = triN, .SDcols = sbs_names]
  # set key
  setkeyv(sigprob, "triN")
  rm(refsig)
}
if(!exists("k11mer_weights")){
  # join data
  k11mer_freq_prob <- k11[sigprob]
  # compute weights
  # P(mut sig & 11mer) = P(mut sig & trinuc) / P(trinuc) * P(11mer)
  k11mer_weights <- k11mer_freq_prob[,lapply(.SD, function(x) (x/freq_triN)*freq_11mer),
                                     .SDcols = sbs_names,
                                     by = kmer]
  rm(k11mer_freq_prob)
}
# create background ====
sample_11mers <- function(SIGNATURE_NAME, seed=20220406) {
  # v2 (06-04-2022)
  # set seed to the ISO date of script creation
  if(!is.null(seed)) set.seed(seed = seed)
  # sample 11-mers
  background_set <- k11mer_weights[, .(kmer, rmultinom(n = 1, size = 10^6, prob = get(SIGNATURE_NAME)))]
  # check length of background set
  #length(background_set)
  return(background_set[, rep(kmer, V1)])
}

# load data required to generate foregrounds ====
cohort_dt <- fread("~/PCAWG/faststorage/Gustav/mutprocesses_project/results/signatures_stratify_genomes_and_11mers/kmer_stats_stratified_by_signature.tsv")
if(!exists("hsdt")) hsdt <- fread("results/hotspots_stratify_11mers/allsignatures_hotspot_deriveddata.tsv")
gregion_dt <- fread("~/PCAWG/faststorage/Gustav/mutprocesses_project/results/genomicregions_stratify_11mers/allsignatures_genomicregion_deriveddata.tsv")
gregion_dt[,feature:=str_extract(gsub("(_|/)", "",feature_hmm), "[A-Za-z_]+")]

cat("####################################\n")
cat("##### data loaded successfully #####\n")
cat("####################################\n")
## define workflow with functions ====
# 1.
create_environment <- function(sigsearch, verbose=FALSE){
  setwd(project_dir)
  # define signature-of-interest ====
  signamefull <- grep(sigsearch, sbs_names, value = TRUE, ignore.case = TRUE)
  signame <- tolower(substr(signamefull, 18, nchar(signamefull)-2))
  # create directory with signature name and subdirectories for in- and output
  dir <- paste0(signame, "_logo")
  if(!dir %in% system("ls results/kplogo/", intern = T)) {
    # feedback directory structure
    if(verbose) {cat(paste0("directory not found.\ncreating new directories:
  results/kplogo/", dir, "/
   - input
   - output
     1. cohort
     2. signature
     3. hotspot
     4. genomicregion\n"))}
    
    # create directories
    setwd("results/kplogo") ; system(command = paste("mkdir", dir))
    setwd(dir) ; system(command = "mkdir input output")
    setwd("output") ; system(command = "mkdir cohort signature hotspot genomicregion")
    setwd("../../../..")
  } else cat(paste("directories for", signame,"already exist\n"))
  
}
# 2.
create_signature_background <- function(sigsearch){
  # working directory should be sbsXX_logo/
  # define signature-of-interest ====
  signamefull <- grep(sigsearch, sbs_names, value = TRUE)
  signame <- tolower(substr(signamefull, 18, nchar(signamefull)-2))
  # create background ====
  # define name of background file
  background.file <- paste0("results/kplogo/",signame,"_logo/","background_", signame, ".tsv")
  # background.file <- paste0("input/backgroundset_",signame,".tsv")
  # sample 11-mer background
  background_set <- sample_11mers(SIGNATURE_NAME = signamefull)
  # save background
  fwrite(x = list(background_set), file = background.file)
  # feedback
  cat("saved background set in:\n\t", background.file, "\n")
  return(background.file)
}
# 3.
create_foreground <- function(sigsearch, level, recurrence="5+", region=""){
  # level = 1, 2, 3, or 4
  # level >= 3: recurrence = 1+, 2+, or 5+
  # level == 4: genomicregion = ...
  
  # define signature-of-interest
  signamefull <- grep(sigsearch, sbs_names, value = TRUE, ignore.case = TRUE)
  signame <- tolower(substr(signamefull, 18, nchar(signamefull)-2))
  cat(paste(signame,"\n", signamefull,"\n"))
  # set new working directory
  dir <- paste0(signame, "_logo")
  # setwd(paste0(project_dir,"/results/kplogo/", dir))
  
  # create foreground
  # foreground.file <- paste0("input/foreground_",signame,"_hotspot5.tsv")
  foreground.file <- paste0("results/kplogo/",dir,"/input/foreground_", signame, "_",
                            switch(level, "1_cohort", "2_11mer", paste0("3_hotspot",recurrence),
                                   paste0("4_region_", gsub("(/|_)", "", region))), ".tsv")
  foreground_set <- switch(level,
                           # 1 (signature cohort)
                           cohort_dt[baselinesig==signamefull, rep(kmer_11, n_snv)],
                           # 2 (signature 11-mers)
                           hsdt[signature==signamefull & rec=="1+"][, rep(kmer_11, n_mut_kmer)],
                           # 3 (signature hotspots)
                           hsdt[signature==signamefull & rec==recurrence][, rep(kmer_11, n_mut_kmer)],
                           # 4 (genomic region)
                           gregion_dt[signature==signamefull & feature==region][, rep(kmer_11, refcount_feature)])
  
  #foreground_set <- hsdt[signature==signamefull & rec==recurrence][, rep(kmer_11, n_mut_kmer)]
  fwrite(x = list(foreground_set), file = foreground.file)
  cat("saved foreground set in:\n\t", foreground.file, "\n")
  return(foreground.file)
}
# 4.
run_kplogo <- function(sigsearch, foreground.file, background.file,
                       output.directory, output.prefix, statistics="s"){
  # if(!(statistics %in% c("p", "b", "f", "s"))){
    # stop("Invalid statistics. Should be p: raw p-value (default), b: Bonferroni corrected p, f: FDR, s: statisitcs")}
  
  time_start <- Sys.time()
  # define signature-of-interest ====
  signamefull <- grep(sigsearch, sbs_names, value = TRUE)
  signame <- tolower(substr(signamefull, 18, nchar(signamefull)-2))
  
  # name prefix (and location) of output files
  # outprefix <- paste0("results/kplogo/", signame,
  #     "_logo/output/kplogo.out.", signame, ".hotspot5")
  # outprefix <- paste0(output.directory, "/", output.prefix)
  outprefix <- output.prefix
  
  setwd(paste0(project_dir, "/", output.directory))
  
  # define command to run kpLogo ====
  cmd <- paste(
    "kpLogo", paste0("../../input/", foreground.file),
    "-bgfile", paste0("../../", background.file), # background sequence file
    "-o", outprefix, # prefix for all output files, default=kpLogo
    "-alphabet ACGT",
    "-max_k 5", # consider all kmers of length 1,2,...,INT. default=4 
    "-shift 0", # max shift (to right) allowed for kmer positions
    "-minCount 0.01", # minimum number of sequences to have this kmer to include in output
    # if smaller than 1, treat as fraction of sequences (default=5)
    "-stack_order 1", # stack residuals by frequency (1) or reverse (-1) or alphabet (0)
    "-pseudo 1", # pseudocount added to background counts. default=1.0. Ignored by -markov
    "-fix 0.75", # fix a position with a specific residual if it occurs in more than maxFreq of the sequences
    paste("-plot", statistics) # which statistics to plot: p: raw p-value (default), b: Bonferroni corrected p, f: FDR, s: statisitcs
    #"-pc 0.01", # Bonferroni corrected p-value cut-off, default=0.05
    #"-FDR" # adjust p value by FDR method ( default is Bonferroni correction)
    # fixed residuals will be plotted as 1.1x of hight of the position with the highest total height
  )
  
  # send command to system ====
  system(command = cmd)
  time_end <- Sys.time()
  cat(paste0("\n",signame," DONE in ", round(difftime(time_end, time_start, units = "secs"),0), " secs\n",Sys.time(),"\n"))
  setwd(project_dir)
}

# create infrastructure ====
# creating directories, foregrounds, and background
# for (sig in sbs_names) {
#   # 1. create environment
#   create_environment(sig)
#   # 2. create background
#   background.file = create_signature_background(sig)
#   # 3. create foreground
#   for( lvl in 1:4 ){
#     # 3a. one foreground for each analysis level
#     if(lvl %in% 1:3) foreground.file = create_foreground(sigsearch = sig, level = lvl)
#     if(lvl==4){
#       for( reg in unique(gregion_dt[signature==sig,feature]) ){
#         tryCatch({
#           if( length(reg) == 0 ) stop("empty data")
#           # 3b. one foreground for each genome region in level 4 analysis
#           if(reg != "Genome-wide") create_foreground(sigsearch = sig, level = 4, region = reg)
#           }, error=function(e){}
#         )
#       }
#     # if(lvl > 4) stop("Level is out of range (1-4)")
#     }
#   }
# }


# save infrastructure in a data table ====
fileoverview_dt <- data.table(sbs = system("ls results/kplogo/", intern = T))
fileoverview_dt <- fileoverview_dt[, .(input = system(paste0("ls results/kplogo/", sbs, "/input"), intern = T)), by = sbs]
fileoverview_dt[, statistic := "s"]
#fileoverview_dt <- fileoverview_dt[, .(statistic = c("s", "p")), by=.(sbs, input)]
fileoverview_dt[, signame := unlist(strsplit(sbs, "_"))[1], by = .(input)]
fileoverview_dt[, level := as.integer(unlist(strsplit(input, "_"))[3]), by=.(input)]
fileoverview_dt[, level_name := switch(as.integer(level), "cohort",
                                       "signature", "hotspot", "genomicregion"), by = .(input,level)]
fileoverview_dt[, output.prefix := paste0(gsub("(foreground_|.tsv|\\+)", "",input),"_",statistic)]
fileoverview_dt[, output.dir := paste0("results/kplogo/",sbs,"/output/",level_name)]
cat("######################################\n")
cat("# file overview created successfully #\n")
cat("######################################\n")


donorsigs <- fread("data/PCAWG_signatures/SA_COMPOSITE_SNV.activity_relative.WHITELIST_SET.20200226.txt")
donorsigs_count <- lapply(colnames(donorsigs)[8:67], function(x) {
  donors <- donorsigs[get(x)>0.05, "icgc_donor_id", with=FALSE]$icgc_donor_id
  ndonors <- length(donors)
  data.table(baselinesig=x, ndonors=ndonors, donors=list(donors))
}) %>% rbindlist()
setkey(donorsigs_count, baselinesig)
dt_sumstat <- fread(file = "results/stats4final_viz/summarystatistics_allsigs.tsv")
# test against mutation rate in previous analysis step ====
for(i in 1:nrow(dt_sumstat)) {
  signature <- dt_sumstat[i, assigned_sig]
  switch(EXPR = dt_sumstat[i]$step,
         dt_sumstat[i, wmrate_prev:=wmrate], #1
         dt_sumstat[i, wmrate_prev:=dt_sumstat[step==1 & assigned_sig==signature, wmrate]], #2
         dt_sumstat[i, wmrate_prev:=dt_sumstat[step==2 & assigned_sig==signature, wmrate]], #3
         dt_sumstat[i, wmrate_prev:=dt_sumstat[step==3 & assigned_sig==signature, wmrate]]) #4
}

dt_sumstat[!is.na(tot_hg19), pval_prev := binom.test(x = as.double(n_snv), n = as.double(tot_hg19),
                                                     p = wmrate_prev/10^6*donorsigs_count[assigned_sig, ndonors], alternative = "greater")$p.value,
           by=.(step, assigned_sig, region)]
dt_sumstat[, pval_bonf_prev := pval_prev*nrow(dt_sumstat)]
dt_sumstat[pval_bonf_prev>1, pval_bonf_prev := 1]
dt_sumstat[, ppval_bonf_prev := -log10(pval_bonf_prev)]
dt_sumstat[is.infinite(ppval_bonf_prev), ppval_bonf_prev:=dt_sumstat[is.finite(ppval_bonf_prev), max(ppval_bonf_prev)]]


# dt_sumstat[( step!=4 | (step==4 & ppval_bonf_prev>2) )]

#selected_sigs <- tolower(str_extract(unique(d$assigned_sig), "SBS[0-9a-d]+"))
# dt_sumstat[step==4 & ppval_bonf_prev>=2, .(n=.N, n_sig=uniqueN(assigned_sig), n_reg=uniqueN(sub("(/|-|_)","",region)))]
# with corrected.p <= 0.01 ( -log10(corrected.p)>=2 ):
#     203 region-signature combos (49 signatures)
selected_sigs_regions <- dt_sumstat[step==4 & ppval_bonf_prev>=2,
  paste0(step, "_", tolower(str_extract(assigned_sig, "SBS[0-9a-d]+")),
  "_", gsub("(/|-|_)","",region))]

fileoverview_dt[level!=4,sbs_step_genomicregion:=paste0(level,"_",signame, "_", "Genomewide")]
fileoverview_dt[level==4,sbs_step_genomicregion:=paste0(level,"_",signame, "_",
                                           str_extract(input, "(?<=region_)[A-z]+"))]
#fileoverview_dt[sbs_step_genomicregion %in% selected_sigs_regions]

# run kplogo on infrastructure ====

# DONE
# fileoverview_dt[signame == "sbs10a" & level<=3,
#                 run_kplogo(sigsearch = signame,
#                            foreground.file = input, 
#                            background.file = paste0("background_",signame,".tsv"),
#                            output.directory = output.dir,
#                            output.prefix = output.prefix,
#                            statistics = statistic),
#                 by = .(input, statistic)]
# DONE
# fileoverview_dt[signame == "sbs10a" & level==4,
#                 run_kplogo(sigsearch = signame,
#                            foreground.file = input, 
#                            background.file = paste0("background_",signame,".tsv"),
#                            output.directory = output.dir,
#                            output.prefix = output.prefix,
#                            statistics = statistic),
#                 by = .(input, statistic)]

# DONE
# lapply(X = unique(fileoverview_dt$signame)[1:4],
#        FUN = function(x) {
#          fileoverview_dt[signame == x & level<=3,
#                          run_kplogo(sigsearch = signame,
#                                     foreground.file = input, 
#                                     background.file = paste0("background_",signame,".tsv"),
#                                     output.directory = output.dir,
#                                     output.prefix = output.prefix,
#                                     statistics = statistic),
#                          by = .(input, statistic)]
#        })


for (x in selected_sigs_regions) {
  s <- Sys.time()
  tryCatch({
    fileoverview_dt[sbs_step_genomicregion == x,
                    run_kplogo(sigsearch = signame,
                               foreground.file = input,
                               background.file = paste0("background_",signame,".tsv"),
                               output.directory = output.dir,
                               output.prefix = output.prefix,
                               statistics = "s"),
                    by = .(input)]
  }, error=function(e){})
  
  e <- Sys.time()
  d <- difftime(e, s)
  cat("####################################################\n")
  cat(paste0("Time spent:\t",d,"\n\tSignature\t", x, " is DONE\n"))
  cat("####################################################\n\n")
}


cat("############################################\n")
cat(paste0("#### End of script\n",Sys.time()," ####\n"))
cat("############################################\n")