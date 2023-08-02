## set options
options(scipen = 999)
## load packages
library(data.table)
library(tidyverse)
library(ggpubr)
library(ggforce)
library(ggseqlogo)
library(scales)
library(gridExtra)
library(grid)
library(gtable)
library(Biostrings)

# load data ====
load("/faststorage/project/PCAWG/Gustav/reference_work/cohorts_palette.RData") # cohorts
col_ref <- setNames(cohorts[1:37, "colcode"], cohorts[1:37, "cohort"]) ; rm(cohorts)

load("/faststorage/project/PCAWG/Gustav/reference_work/signature_etiologies.RData") # proposed_etiologies
sbs_columns <- proposed_etiologies$BI_COMPOSITE #; rm(proposed_etiologies)

load("/faststorage/project/PCAWG/Gustav/reference_work/motifs/kmers_k11_wConsensusMatrices.RData") # k11c_dt
k11c_dt[, freq_hg19 := tot_hg19/sum(tot_hg19)]

load("/faststorage/project/PCAWG/Gustav/results/motif_analysis/sbs_k11_pfm_weighted_by_signatures.Rdata") # sbs_k11


## read data
s <- Sys.time() ; cat("\tStart:\t",as.character(s), "\n\tExpected end:\t", as.character(s+60*5), "\n")
snvsigs_maf <- fread(file = "/faststorage/project/PCAWG/Gustav/data/signatures/signature_snv_autosome_n_project_header.maf", key = "Variant_Classification")
e <- Sys.time() ; difftime(e,s) # 4.2 - 6.9 min

#2. Filter data ====

#cds <- c("Silent", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_SNP")
noncds <- c("IGR", "5'Flank", "5'UTR", "Intron", "3'UTR", "RNA", "lincRNA")
blacklistGenes <- data.table(Gene="TERT", Chromosome="chr5", Start_position=1253287, End_position=1295162+2000, Strand="-", 
                             Source="UCSCGenomeBrowser_hg19", Comment="Known drivers in promoter region. Oncogene.")


s <- Sys.time()
cs_dt <- snvsigs_maf[noncds] #; if(exists("cs_dt")) rm(snvsigs_maf)
cs_dt <- cs_dt[ !(Chromosome==blacklistGenes$Chromosome & Start_position>=blacklistGenes$Start_position & End_position<=blacklistGenes$End_position)]
e <- Sys.time() ; difftime(e,s) # 2.9 min
rm(noncds, blacklistGenes)

# manipulate columns ====
cs_dt <- cs_dt[order(as.numeric(str_sub(Chromosome, 4, -1)), End_position)]
cs_dt[, dist_all := diff(End_position)]
cs_dt[n_pc >= 2, dist_rec2 := diff(End_position)]
cs_dt[n_pc >= 3, dist_rec3 := diff(End_position)]
cs_dt[n_pc >= 4, dist_rec4 := diff(End_position)]
cs_dt[n_pc >= 5, dist_rec5 := diff(End_position)]
cs_dt[n_pc >= 6, dist_rec6 := diff(End_position)]
cs_dt[n_pc >= 7, dist_rec7 := diff(End_position)]

cs_dt[, bin_pos := as.numeric(str_sub(Chromosome, 4, -1)) + round(End_position/100, 0)/10^7]

cs_dt[, `:=`(bin_n = .N),
      by = bin_pos]
s <- Sys.time()
cs_dt[, bin_n_uniq := uniqueN(eid), by = bin_pos]
e <- Sys.time()
e-s
cs_dt[, `:=`(bin_dist_mean = mean(dist_all)), 
      by = bin_pos]

cs_dt[bin_n_uniq>15 & rec >= 2 & 
        str_sub(flankSeq_ag, 11-3, 11+2) == "AACTTA" & 
        grepl("Eso", Tumor_Type),
      c("eid", "bin_pos", "n_pc", 
        "Donor_ID", "Tumor_Type", 
        "flankSeq_ag", "bin_n_uniq", "dist_all")]

cs_dt[, is_AACTT := str_sub(flankSeq_ag, 11-3, 11+1) == "AACTT"]
#cs_dt[is_AACTT==1, dist_AACTT := diff(End_position)]
cs_dt[is_AACTT==1 & rec >= 2, 
      dist_AACTT_hs := pmin(c(NA, diff(End_position)), c(diff(End_position), NA)) ]

cs_dt[is.finite(dist_AACTT_hs) & 
        dist_AACTT_hs <= 10 & dist_AACTT_hs >= 0][order(dist_AACTT_hs)][,c("eid", "bin_pos", "n_pc",
       "Donor_ID", "Tumor_Type",
       "flankSeq_ag", "dist_AACTT_hs")]


unique(cs_dt[dist_rec5 < 20 & bin_n_uniq > 2 & 
        str_sub(flankSeq_ag, 11-2, 11+2) == "ACTTA",
      c("eid", "bin_pos", "n_pc", 
        "Donor_ID", "Tumor_Type", 
        "flankSeq_ag", "bin_n_uniq", 
        "dist_all", "dist_rec5")], by = "eid")[, dist_bin := c(NA, diff(bin_pos))][order(dist_bin)]

cs_dt[Chromosome == "chr1" &
        End_position >= 217449090 & End_position <= 217449150, 
      c("eid", "bin_pos", "n_pc",
        "Donor_ID", "Tumor_Type",
        "flankSeq_ag", "bin_n_uniq",
        "dist_all", "dist_rec5")]
# bin

## 1.248382
cs_dt[Chromosome == "chr1" &
        End_position >= 248381705 & End_position <= 248381770] %>% 
  ggplot(aes(Donor_ID, End_position, col = Tumor_Type, shape = rec)) + 
  geom_point() + 
  geom_text_repel(aes(label = str_sub(flankSeq_ag, 11-1, 11+1))) +
  coord_flip() 

## 3.111499  liver, eso
cs_dt[bin_pos >= 3.11149925 & bin_pos <= 3.111499350] %>% 
  ggplot(aes(Donor_ID, End_position, col = Tumor_Type, shape = rec)) + 
  geom_point() + 
  geom_text_repel(aes(label = str_sub(flankSeq_ag, 11-1, 11+1))) +
  coord_flip() 

## 4.030240 eso
cs_dt[bin_pos >= 4.0302397 & bin_pos <= 4.0302400] %>% 
  ggplot(aes(Donor_ID, End_position, col = Tumor_Type, shape = rec)) + 
  geom_point() + 
  geom_text_repel(aes(label = str_sub(flankSeq_ag, 11-1, 11+1))) +
  coord_flip() 

# 14.024895 liver
cs_dt[bin_pos >= 14.024895 & bin_pos <= 14.0248955]


cs_dt[bin_dist_mean>0 & bin_dist_mean < 10 & bin_n>10] %>%
  ggplot(aes(bin_pos, 1/bin_dist_mean,
             col = str_sub(flankSeq_ag, 11-1, 11+1)=="CTT")) + 
  geom_point() +
  facet_wrap(factor(round(bin_pos, 0), levels = c(1:22))~., 
             scales = "free_x") +
  theme(strip.background = element_blank())


cs_dt[Chromosome == "chr2" &
        End_position>=4*10^5 & 
        End_position<=5*10^5 &
        dist_all <= 20] %>%
  ggplot(aes(End_position, dist_all, 
             col = str_sub(flankSeq_ag, 11-1, 11+1)=="CTT",
             shape = str_sub(flankSeq_ag, 11-1, 11+1)=="CTT")) + 
  geom_point(alpha = 0.5)+
  theme_pubr(legend = "none")

cs_dt[Chromosome == "chr2" &
        End_position>=4.1*10^5 & 
        End_position<=4.25*10^5]




rc <- function(x) as.character(reverseComplement(DNAStringSet(x)))

s <- Sys.time()
cs_dt[, `:=`(normid = NULL,
             Tumor_Type = Project_Code,
             eid = paste(str_sub(Chromosome, 4, -1), End_position, sep = "_"),
             flankSeq_ag = ifelse(Reference_Allele %in% c("C","T"), toupper(ref_context), rc(ref_context)),
             rec = ifelse(n_pc>=7, "7+", n_pc))]
e <- Sys.time() ; difftime(e,s) # 5.1 - 5.4 min


## count PCAWG 11-mers
cs_dt[, kmer := str_sub(flankSeq_ag, 11-5, 11+5)]
pcawg_all11mers <- cs_dt[, .(tot_mut = .N), by = c("kmer", "rec")]
pcawg_all11mers[, rec_bin := ifelse(rec==1, "rec1",
                                ifelse(rec>=2 & rec<=4, "rec2_4",
                                       "rec5plus"))]

pcawg_all11mers_rec <- pcawg_all11mers[, .(tot_mut = sum(tot_mut)), by = c("rec_bin", "kmer")] %>% spread(key = "rec_bin", value = "tot_mut", fill = 0)
pcawg_all11mers_rec[, tot_mut := rec1 + rec2_4 + rec5plus]

## set keys
setkey(pcawg_all11mers_rec, kmer) ; setkey(k11c_dt, kmer)

## merge
k11c_dt2 <- pcawg_all11mers_rec[k11c_dt[tot_hg19>0]]
k11c_dt2[is.na(tot_mut), tot_mut := 0]
k11c_dt2[, `:=`(cm = NULL,
                center = str_sub(kmer, 6, 6),
                mutrate = ifelse(tot_mut>0, tot_mut/tot_hg19, 0))]
fwrite(k11c_dt2, file = "PCAWG/faststorage/Gustav/mutprocesses_project/results/11mer_count_hg19_PCAWG.tsv")
## bin data
#k11c_dt2 %>% 

# visualize data ====
## histogram
### count mutations
p_hist_mut <- gghistogram(k11c_dt2, x = "tot_mut") + 
  scale_x_continuous(limits = c(1, 3.7*10^4), trans = "log10", labels = comma) +
  scale_y_continuous(labels = comma) # limits = c(0, 0.35*10^6), 

### count occurrence in reference
p_hist_ref <- gghistogram(k11c_dt2, x = "tot_hg19") +
  scale_x_continuous(limits = c(1, 4.7*10^6), trans = "log10", labels = comma) +
  scale_y_continuous(labels = comma) # limits = c(0, 0.35*10^6), 

## scatter plot
p_scatter <- ggplot(k11c_dt2,
                    aes(tot_hg19, tot_mut + 0.2, col = center)) + 
  geom_point(aes(alpha = mutrate)) +
  geom_abline(intercept = c(0, -1, -2, -3), slope = 1, linetype = "solid", col = "grey") +
  theme_pubr(legend = "right")+
  scale_x_continuous(limits = c(1, 4674610), trans = "log10",
                     labels = comma) +
  scale_y_continuous(limits = c(.1, 36897), trans = "log10", labels = comma)
  

## arrange
ggarrange(p_hist_ref, NULL,
          p_scatter, p_hist_mut + rotate(),
          ncol = 2, nrow = 2, align = "hv",
          widths = c(2, 1), heights = c(1, 2))


## heatmap density
smoothScatter(k11c_dt2[1:100000, lapply(.SD, log10), .SDcols = c("tot_hg19", "tot_mut")])
  


ggplot(k11c_dt2[1:100000], aes(tot_hg19, tot_mut)) +
  stat_bin2d() +
  #stat_density_2d(geom = "raster", aes(fill = (..density..^(1/4))), contour = FALSE) + 
  geom_abline(intercept = c(0, -1, -2, -3), slope = 1, linetype = "solid", col = "grey") +
  #scale_fill_gradient2(midpoint = 0.25, low = "grey", mid = "yellow", high = "red") +
  scale_x_log10() + scale_y_log10() +
  theme_pubr()
