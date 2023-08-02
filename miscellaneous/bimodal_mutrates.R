# look into the bimodality in mutation rates of 11-mers with 5+ hotspots


# define environment ====
setwd("~/PCAWG/faststorage/Gustav/")
#setwd("~/GenomeDK/PCAWG/faststorage/Gustav/")
options(scipen = 999)
library(tidyverse) # data manipulation
library(data.table) # data manipulation
library(Biostrings) # DNA seq manipulation (consensusMatrix)
library(BSgenome.Hsapiens.UCSC.hg19) # assess hg19 seq
library(ggseqlogo) # plot logo plots
library(ggpubr) # pretty ggplot themes
library(scales) # pretty numbers
library(patchwork) # plot arrangement

# load data ====
derivedData_dt <- fread("mutprocesses_project/results/hotspots_stratify_11mers/allsignatures_hotspot_deriveddata.tsv")

# compute expected pfm from hg19 autosomes
b <- Reduce("+",lapply(1:22, function(chr) oligonucleotideFrequency(Hsapiens[[chr]], width = 1)))
pfm_exp <- matrix( (b/sum(b)), nrow = 4, ncol = 11, dimnames = list(names(b)))
rm(b)

# filter and plot ====
p_logo_s62_high <- derivedData_dt[analysis=="hotspot" & grepl("SBS62",signature) & rec=="5+" & mrate>=mrate.sig11mer_mean,
               consensusMatrix(kmer_11, as.prob = T)] %>% ggseqlogo() + ggtitle(label = "SBS62 | rec5+ | high mutrates")
p_logo_s62_low <- derivedData_dt[analysis=="hotspot" & grepl("SBS62",signature) & rec=="5+" & mrate<mrate.sig11mer_mean,
               consensusMatrix(kmer_11, as.prob = T)] %>% ggseqlogo() + ggtitle(label = "SBS62 | rec5+ | low mutrates")

# save plot ====
ggsave(plot=p_logo_s62_high+p_logo_s62_low+plot_layout(ncol=1), 
       device="pdf", width = 7, height = 7, 
       filename = "mutprocesses_project/plots/hotspots_stratify_11mers/bimodal_mutrates/s62_rec5+_highlow.pdf")

# filter and plot ====


pfm <- function(seq, as.prob=F) consensusMatrix(DNAStringSet(seq), as.prob = as.prob)[1:4,]
plot_logo_bits_hg19auto_v2 <- function(pfm_obs) {
  # 1. loop over recurrence to access consensusmatrix
  # 2. loop over positions (1-11)
  # 3. loop over bases (1-4)
  
  # pfm_obs <- subset_dt[, pfm(rep(kmer_11, n_snv), as.prob=T)]
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


derivedData_dt_s7a <- derivedData_dt[analysis=="hotspot" &
                                       grepl("SBS7a",signature) & rec=="5+"]

p_logo_s7a_high <- derivedData_dt_s7a[mrate>=mrate.sig11mer_mean,
              consensusMatrix(kmer_11, as.prob = T)] %>%
  plot_logo_bits_hg19auto_v2() + ggtitle(label = "SBS7a | rec5+ | high mutrates")


p_logo_s7a_low <- derivedData_dt[analysis=="hotspot" & 
                                   grepl("SBS7a",signature) & rec=="5+" & 
                                   mrate<mrate.sig11mer_mean,
              consensusMatrix(kmer_11, as.prob = T)] %>%
  plot_logo_bits_hg19auto_v2() + ggtitle(label = "SBS7a | rec5+ | high mutrates")
  # ggseqlogo() + ggtitle(label = "SBS7a | rec5+ | low mutrates")
# save plot ====
ggsave(plot=p_logo_s7a_low+p_logo_s7a_high+plot_layout(nrow=1), 
       device="pdf", width = 7, height = 2, 
       filename = "mutprocesses_project/plots/hotspots_stratify_11mers/bimodal_mutrates/s7a_rec5+_highlow_v2.pdf")

