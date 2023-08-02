# report APOBEC findings
# characterize the two substypes/motifs according to span,
#   genomic distribution, mutation rate, patient mutability, etc.. (1 day)
#   The statistics in plain text will do.

# APOBEC motifs:
# APOBEC3A: ytCaT
# APOBEC3B: rtCaT


# define environment ====
options(scipen=999)
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project/")
library(data.table)
library(stringr)
library(tidyverse)
library(Biostrings)
library(ggseqlogo)
library(ggpubr)
library(ggforce)



# load data ====
sbs_etio <- fread("data/reference_signatures/SignatureAnalyzer_SBS_etiology.tsv")
sbs_etio[grepl("APOBEC", Etiology_full)]
k11c_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")

cohort_dt <- fread("~/PCAWG/faststorage/Gustav/mutprocesses_project/results/signatures_stratify_genomes_and_11mers/kmer_stats_stratified_by_signature.tsv")
if(!exists("hsdt")) hsdt <- fread("results/hotspots_stratify_11mers/allsignatures_hotspot_deriveddata.tsv")
#gregion_dt <- fread("~/PCAWG/faststorage/Gustav/mutprocesses_project/results/genomicregions_stratify_11mers/allsignatures_genomicregion_deriveddata.tsv")

colnames(cohort_dt)

cohort_apobec <- cohort_dt[grepl("(SBS2_|SBS13|SBS69)", baselinesig)]
cohort_apobec[, APOBEC := "FALSE"]
cohort_apobec[str_sub(kmer_11, 6-2, 6+1) %in% c("CTCA","TTCA"), APOBEC:="A3A"]
cohort_apobec[str_sub(kmer_11, 6-2, 6+1) %in% c("ATCA","GTCA"), APOBEC:="A3B"]

cohort_apobec[, uniqueN(unlist(list(str_split(donors, pattern = ",")))), by=baselinesig]
cohort_apobec[grepl("SBS2_", baselinesig), n_donors := 303]
cohort_apobec[grepl("SBS13", baselinesig), n_donors := 330]
cohort_apobec[grepl("SBS69", baselinesig), n_donors := 468]
cohort_apobec[, wmrate:=n_snv/tot_hg19/n_donors*10^6]
# cohort_apobec[,colnames(cohort_apobec)[-3], with=F]


hsdt_apobec <- hsdt[grepl("(SBS2_|SBS13|SBS69)", signature)]
# Y = Pyrimidine = C or T
# R = Purine = A or G
hsdt_apobec[, APOBEC := "FALSE"]
hsdt_apobec[str_sub(kmer_11, 6-2, 6+1) %in% c("CTCA","TTCA"), APOBEC:="A3A"]
hsdt_apobec[str_sub(kmer_11, 6-2, 6+1) %in% c("ATCA","GTCA"), APOBEC:="A3B"]
hsdt_apobec[grepl("SBS2_", signature), n_donors := 303]
hsdt_apobec[grepl("SBS13", signature), n_donors := 330]
hsdt_apobec[grepl("SBS69", signature), n_donors := 468]
hsdt_apobec[, wmrate:=n_mut_kmer/tot_hg19/n_donors*10^6]


# relative contribution =======
cohort_apobec[,.(wmrate=sum(n_snv)/sum(tot_hg19)/unique(n_donors)*10^6,
                 span=paste(round(sum(tot_hg19)/10^6,2), "Mb"),
                 n_donors=unique(n_donors)),
              by=.(baselinesig, APOBEC)][,
  r:=wmrate/wmrate[APOBEC=="FALSE"], by=baselinesig][,]

cohort_apobec[baselinesig==maxsig_kmer,
              .(wmrate=sum(n_snv)/sum(tot_hg19)/unique(n_donors)*10^6,
                 span=paste(round(sum(tot_hg19)/10^6,2), "Mb"),
                 n_donors=unique(n_donors)),
              by=.(baselinesig, APOBEC)][,
  r:=wmrate[APOBEC=="A3A"]/wmrate[APOBEC=="A3B"], by=baselinesig][,]

hsdt_apobec[rec=="5+",.(wmrate=sum(n_mut_kmer)/sum(tot_hg19)/unique(n_donors)*10^6,
                 span=paste(round(sum(tot_hg19)/10^6,2), "Mb"),
                 n_donors=unique(n_donors)),
              by=.(signature, APOBEC)][,
  r:=wmrate[APOBEC=="A3A"]/wmrate[APOBEC=="A3B"], by=signature][,]

# test for normality
p_qq <- cohort_apobec[, .(mrate_log=log10(wmrate)), by=.(APOBEC, baselinesig)] %>%
  dplyr::sample_n(size=10^4) %>%
  ggplot(aes(sample = mrate_log)) +
  stat_qq() + stat_qq_line() +
  facet_grid(baselinesig~APOBEC, scales="free")
ggsave(plot = p_qq, filename = "plots/apobec/mrate_qqplots.png", device = "png")

p_hist <- cohort_apobec[, .(mrate_log=log10(wmrate)), by=.(APOBEC, baselinesig)] %>%
  ggplot(aes(x = mrate_log)) +
  geom_histogram() +
  facet_grid(APOBEC~baselinesig, scales="free")
ggsave(plot = p_hist, filename = "plots/apobec/mrate_hist.png", device = "png")

p_violin <- cohort_apobec[, .(mrate_log=log10(wmrate)), by=.(APOBEC, baselinesig)] %>%
  ggplot(aes(x = APOBEC, y = mrate_log, fill=APOBEC)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = .5, fill="transparent") +
  scale_fill_brewer(palette = "Pastel1") +
  facet_grid(baselinesig~., scales="free")
ggsave(plot = p_violin, filename = "plots/apobec/mrate_violin.png", device = "png")

# p_qq <- hsdt_apobec[rec=="1+", .(mrate_log=log10(mrate)), by=.(rec, APOBEC, signature)] %>%
#   ggplot(aes(sample = mrate_log)) +
#   stat_qq() + stat_qq_line() +
#   facet_grid(signature~APOBEC, scales="free")

# qqplots looks fine for each APOBEC group in rec 1+

# which signature is most interesting to test?

stat_cohort <- cohort_apobec[, .(step="cohort",
                  n=.N,
                  wmrate=sum(n_snv)/sum(tot_hg19)/unique(n_donors)*10^6,
                  span=sum(as.numeric(tot_hg19)),
                  consensus_motif=consensusString(rep(kmer_11,n_snv))),
              by=.(APOBEC, signature=baselinesig, n_donors)]

stat_11mer <- cohort_apobec[baselinesig==maxsig_kmer,
              .(step="11mer",
                n=.N,
                wmrate=sum(n_snv)/sum(tot_hg19)/unique(n_donors)*10^6,
                span=sum(as.numeric(tot_hg19)),
                consensus_motif=consensusString(rep(kmer_11,n_snv))),
              by=.(APOBEC, signature=baselinesig, n_donors)]

stat_hotspot <- hsdt_apobec[rec=="5+",
            .(step="hotspot",
              n=.N,
              wmrate=sum(n_mut_kmer)/sum(tot_hg19)/unique(n_donors)*10^6,
              span=sum(tot_hg19),
              consensus_motif=consensusString(rep(kmer_11,n_mut_kmer))),
            by=.(APOBEC, signature, n_donors)]

stat_apobec <- rbind(stat_cohort, stat_11mer, stat_hotspot)
stat_apobec[, step:=factor(step, levels=c("cohort","11mer","hotspot"))]
stat_apobec[, span_mb:=paste(round(span/10^6,2), "Mb")]
stat_apobec[, consensus_motif_trim:=paste0(
  str_sub(consensus_motif, 1,5),
  "[",str_sub(consensus_motif, 6,6),">N]",
  str_sub(consensus_motif, 7,11))]
stat_apobec[, consensus_motif_trim:=gsub("(^[?]+|[?]+$)", "", consensus_motif_trim)]
stat_apobec[, sig:=str_sub(signature, 18,-3)]

stat_apobec[, .(APOBEC, sig, n_donors, step,
                n=format(n, big.mark=","), span_mb, wmrate
                # motif=consensus_motif_trim
                )][,
  fc := wmrate/wmrate[step=="cohort"], by=.(sig, APOBEC)][order(sig, APOBEC, step)]


hsdt_apobec[rec=="5+" & APOBEC=="A3A",
            consensusMatrix(rep(kmer_11, n_mut_kmer), as.prob=T)]

# testing mutation rate changes =====
#### from cohort to hotspots ####

cohort_apobec[grepl("SBS2", baselinesig) & APOBEC=="A3A",
              .(step="cohort", kmer_11, wmrate, span = as.numeric(tot_hg19))]
hsdt_apobec[grepl("SBS2", signature) & rec=="5+" & APOBEC=="A3A",
            .(step="hotspot", kmer_11, wmrate, span=as.numeric(tot_hg19))]

stat_apobec[sig=="SBS69" & APOBEC!="FALSE"][order(APOBEC, step)]

test_apobec_mrate <- function(sig){
  p_a3a <- t.test(
    x = cohort_apobec[grepl(sig, baselinesig) & APOBEC=="A3A", wmrate],
    y = hsdt_apobec[grepl(sig, signature) & rec=="5+" & APOBEC=="A3A", wmrate]
    )$p.value
  
  p_a3b <- t.test(
    x = cohort_apobec[grepl(sig, baselinesig) & APOBEC=="A3B", wmrate],
    y = hsdt_apobec[grepl(sig, signature) & rec=="5+" & APOBEC=="A3B", wmrate]
  )$p.value
  
  c(paste("A3A: p =", round(p_a3a,5)), paste("A3B: p =", round(p_a3b,5)))
}

test_apobec_mrate("SBS2")
test_apobec_mrate("SBS13")
test_apobec_mrate("SBS69")


t.test(
  x = cohort_apobec[grepl("SBS2", baselinesig) & APOBEC=="A3A", wmrate],
  y = cohort_apobec[grepl("SBS2", baselinesig) & APOBEC=="A3B", wmrate]
) %>% {log10(.$p.value)}
t.test(
  x = cohort_apobec[grepl("SBS13", baselinesig) & APOBEC=="A3A", wmrate],
  y = cohort_apobec[grepl("SBS13", baselinesig) & APOBEC=="A3B", wmrate]
) %>% {log10(.$p.value)}
t.test(
  x = cohort_apobec[grepl("SBS69", baselinesig) & APOBEC=="A3A", wmrate],
  y = cohort_apobec[grepl("SBS69", baselinesig) & APOBEC=="A3B", wmrate]
) %>% {log10(.$p.value)}

# hotspots
t.test(
  x = hsdt_apobec[grepl("SBS2", signature) & rec=="5+" & APOBEC=="A3A", wmrate],
  y = hsdt_apobec[grepl("SBS2", signature) & rec=="5+" & APOBEC=="A3B", wmrate]
)
t.test(
  x = hsdt_apobec[grepl("SBS13", signature) & rec=="5+" & APOBEC=="A3A", wmrate],
  y = hsdt_apobec[grepl("SBS13", signature) & rec=="5+" & APOBEC=="A3B", wmrate]
)
t.test(
  x = hsdt_apobec[grepl("SBS69", signature) & rec=="5+" & APOBEC=="A3A", wmrate],
  y = hsdt_apobec[grepl("SBS69", signature) & rec=="5+" & APOBEC=="A3B", wmrate]
)


# plotting linegraph =====

scientific_10 <- function(x) parse(text=paste0(gsub("e", " %*% 10^", scales::scientific_format()(x)),"*'%'"))
axis_label <- function(x) ifelse(x*100<0.01 & x!=0, scientific_10(x*100), paste0(signif(x*100,2),"%"))
library(scales)
library(ggrepel)
p_linegraph <- stat_apobec %>%
  ggplot(aes(x = span, y = wmrate, fill=APOBEC, label=consensus_motif_trim)) +
  geom_line(color="grey") +
  geom_point(shape=21, col="black", size=3) +
  geom_text_repel(vjust=0) +
  annotation_logticks(sides = "l") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_continuous(name="Genomic span", trans=trans_reverser("log10"),
                     breaks = 10^seq(0, 9, by=3),
                     labels = c("1 bp", "1 kb", "1 Mb", "1 Gb"),
                     sec.axis = sec_axis(name = "", #"Genomic span (%)",
                                         trans = ~./sum(k11c_dt$tot_hg19),
                                         breaks = 10^seq(9, 0, by=-3)/sum(k11c_dt$tot_hg19),
                                         labels = axis_label)) +
  scale_y_log10(name="Mutation rate (SNV/patient/Mb)",
                breaks=10^(0:6),
                labels=trans_format("log10", math_format(10^.x)),
                sec.axis = dup_axis()) +
  # labs(title = paste("Signature", str_extract(sig, "(?<=SBS)[0-9a-z]+"))) +
  facet_grid(.~signature, scales="free") +
  # facet_wrap(.~factor(str_extract(assigned_sig, "SBS[0-9a-z]+"), 
                      # levels = str_extract(refsig_etio$BI_COMPOSITE, "SBS[0-9a-z]+"))) +
  coord_cartesian(xlim = c(2.6*10^9, 100), ylim=c(1, max(stat_apobec$wmrate)*1.5), clip = "on") +
  theme_pubr()

ggsave(plot = p_linegraph, device="png", width=14, height=7,
       filename="plots/apobec/linegraph_sbs2_apobec.png")

plot_linegraph_v2 <- function(sig="SBS", ppval_threshold=10, single_region=TRUE, keep.strip=FALSE){
  
  theme_linegraph <- theme(strip.background = element_blank(),
                           strip.text = element_blank())
  if(keep.strip) theme_linegraph <- theme(strip.placement = "outside")
  
  scientific_10 <- function(x) parse(text=paste0(gsub("e", " %*% 10^", scales::scientific_format()(x)),"*'%'"))
  axis_label <- function(x) ifelse(x*100<0.01 & x!=0, scientific_10(x*100), paste0(signif(x*100,2),"%"))
  
  
  d <- dt_sumstat[grepl(sig, assigned_sig)]
  if(single_region) d <- d[step<4 | (step==4 & wmrate==max(wmrate[step==4], na.rm=T))]
  
  d %>%
    ggplot(aes(tot_hg19, wmrate)) +
    # plotting data
    # draw lines
    geom_line(color="grey") + # data=d, 
    # geom_line(data=d[step<4], color="grey") +
    # geom_segment(data=d[step>=3 & !(ppval_bonf_prev<=ppval_threshold & step==4), 
    #                      .(x=tot_hg19[step==3], y=wmrate[step==3],
    #   xend=tot_hg19[step==4], yend=wmrate[step==4]), by=assigned_sig],
    #   aes(x=x, y=y, xend=xend, yend=yend), color="grey") +
    # add points
    geom_point(aes(fill=-log10(pval_bonf_prev)), shape=21, col="black", size=3) +
    # add text annotations
    ## label with analysis step / region
    geom_text(data=d[step==1], aes(label=label), vjust=1, hjust=0, size=3) +
    geom_text(data=d[step%in%2:4], aes(label=label), vjust=-1, size=3) +
    # geom_text(data=d[step==4,
    #                  .(region=region[which.max(wmrate)],wmrate=wmrate[which.max(wmrate)], tot_hg19=tot_hg19[which.max(wmrate)]),
    #                  by=.(assigned_sig)], aes(label=region), vjust=-1, size=3) +
    ## label with mutation rate increase
    geom_text(data=d[step!=4], aes(label=label_increase), vjust=1, hjust=0, size=4) +
    geom_text(data=d[step==4], aes(label=label_increase), vjust=1, hjust=1, size=4) +
    # controling scales
    annotation_logticks(sides = "l") +
    scale_x_continuous(name="Genomic span", trans=trans_reverser("log10"),
                       breaks = 10^seq(0, 9, by=3),
                       labels = c("1 bp", "1 kb", "1 Mb", "1 Gb"),
                       #labels = function(x) format(x/10^6, drop0trailing = TRUE, big.mark = ","),
                       sec.axis = sec_axis(name = "", #"Genomic span (%)",
                                           trans = ~./sum(k11c_dt$tot_hg19),
                                           breaks = 10^seq(9, 0, by=-3)/sum(k11c_dt$tot_hg19),
                                           labels = axis_label)) + #function(x) paste0(signif(x*100, 2),"%"))) +
    scale_y_log10(name="Mutation rate (SNV/patient/Mb)",
                  breaks=10^(0:6),
                  labels=trans_format("log10", math_format(10^.x)),
                  sec.axis = dup_axis()) +
    scale_fill_gradient(name=expression(-log[10]~`p-value`), #"-log10(p-value x 817)", 
                        low = "white", high="brown", na.value = "brown",
                        limits=c(0,320), breaks=c(0,300), #midpoint = 320/2, 
                        guide = guide_colorbar(direction = "horizontal",
                                               barheight = unit(2,"mm"),
                                               barwidth = unit(1.5, "cm"),
                                               ticks = FALSE,
                                               title.position = "top",
                                               title.theme = element_text(size=10, hjust=.5),
                                               frame.colour = "black")) +
    labs(title = paste("Signature", str_extract(sig, "(?<=SBS)[0-9a-z]+"))) +
    facet_wrap(.~factor(str_extract(assigned_sig, "SBS[0-9a-z]+"), 
                        levels = str_extract(refsig_etio$BI_COMPOSITE, "SBS[0-9a-z]+"))) +
    coord_cartesian(xlim = c(2.6*10^9, 100), ylim=c(1, max(d$wmrate)*1.5), clip = "off") +
    theme_pubr() +
    theme_linegraph +
    theme(legend.position = c(.8,.15), legend.background = element_blank(),
          axis.text.y.right = element_blank(), axis.title.y.right = element_blank(), axis.ticks.y.right = element_blank())
}




# donor overlap ====
donor_labels <- fread("data/PCAWG_donorinfo/donor_labels.tsv")
donors_apobec <- cohort_apobec[,.(donors=unique(unlist(list(str_split(donors, pattern = ","))))),
              by=baselinesig]
donors_apobec <- merge(donors_apobec, donor_labels, by.x="donors", by.y="Donor_ID")
donors_apobec[, .N, by=.(baselinesig, Tumor_Type)][order(baselinesig, -rank(N))]


l_donors <- lapply(unique(donors_apobec$baselinesig),
       function(signature) {
         donors_apobec[baselinesig==signature, donors]
       })
names(l_donors) <- unique(donors_apobec$baselinesig)


donors_apobec[, .N, by=baselinesig]
# S2    303
# S13   330
# S69   468

donor_labels[Donor_ID %in%
               intersect(l_donors[["BI_COMPOSITE_SNV_SBS2_P"]],
                         l_donors[["BI_COMPOSITE_SNV_SBS13_P"]]),
             .N, by=Tumor_Type][order(-rank(N))]
# 162
uniqueN(union(l_donors[["BI_COMPOSITE_SNV_SBS2_P"]], l_donors[["BI_COMPOSITE_SNV_SBS13_P"]]))
330+303-162

donor_labels[Donor_ID %in%
               intersect(l_donors[["BI_COMPOSITE_SNV_SBS13_P"]],
                         l_donors[["BI_COMPOSITE_SNV_SBS69_P"]]),
             .N, by=Tumor_Type][order(-rank(N))]
# 71

donor_labels[Donor_ID %in%
               intersect(l_donors[["BI_COMPOSITE_SNV_SBS2_P"]],
                         l_donors[["BI_COMPOSITE_SNV_SBS69_P"]]),
             .N, by=Tumor_Type][order(-rank(N))]
# 71

donor_labels[Donor_ID %in% Reduce(intersect, l_donors), .N, by=Tumor_Type]
# 15 patients are in all three cohorts
# 12 Breast-AdenoCA, 3 Panc-AdenoCA


# plot logo ====

stat_apobec_hotspot <- hsdt_apobec[rec=="5+",
            .(step="hotspot",
              n=.N,
              wmrate=sum(n_mut_kmer)/sum(tot_hg19)/unique(n_donors)*10^6,
              span=sum(tot_hg19),
              consensus_motif=consensusString(rep(kmer_11, n_mut_kmer)),
              pfm_motif=list(consensusMatrix(rep(kmer_11, n_mut_kmer), as.prob=T))),
            by=.(APOBEC, signature, n_donors)]
stat_apobec_hotspot

sequences <- list(
"A3A"=hsdt_apobec[rec=="5+" & APOBEC=="A3A", rep(kmer_11, n_mut_kmer)],
"A3B"=hsdt_apobec[rec=="5+" & APOBEC=="A3B", rep(kmer_11, n_mut_kmer)],
"none"=hsdt_apobec[rec=="5+" & APOBEC==FALSE, rep(kmer_11, n_mut_kmer)])

str(sequences)

library(Biostrings)

b <- Reduce("+", lapply(1:22, function(chr)
  letterFrequency(BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr]],
                  letters=c("A","C","G","T"), as.prob=FALSE)))

consensusMatrix(sequences[['A3A']], as.prob=TRUE) %>% round(1)
consensusMatrix(sequences[['A3B']], as.prob=TRUE) %>% round(1)
consensusMatrix(sequences[['none']], as.prob=TRUE) %>% round(1)

#### Information content ####
bits <- function(C, background="uniform"){
  PPM <- function(C) consensusMatrix(C, as.prob=TRUE)+0.001
  ICtotal <- function(C) log2(length(C))
  U <- function(C) -sum(PPM(C) * log2(PPM(C)), na.rm = TRUE)
  ICfinal <- function(C) ICtotal(C) - U(C)
  IC <- function(C) PPM(C) * ICfinal(C)
  if(background=="uniform") return(IC(C))
  B <- b/sum(b)
  IC <- function(C) pmax(PPM(C) * log2(PPM(C) / B), 0)
  return(IC(C))
}

bits(sequences[['A3A']]) %>% round(1)
bits(sequences[['A3B']]) %>% round(1)
bits(sequences[['none']]) %>% round(1)

bits(sequences[['A3A']], background="non-uniform") %>% round(1)
bits(sequences[['A3B']], background="non-uniform") %>% round(1)
bits(sequences[['none']], background="non-uniform") %>% round(1)

paste(c("A","C","G","T")[apply(bits(sequences[['A3A']]), 2, which.max)], collapse = "")
paste(c("A","C","G","T")[apply(bits(sequences[['A3B']]), 2, which.max)], collapse = "")

# sequences <- stat_apobec_hotspot[, pfm_motif]
# names(sequences) <- stat_apobec_hotspot[, paste(APOBEC,"-",str_sub(signature, 18,-3))]

library(ggseqlogo)
p <- ggseqlogo(data=sequences, facet='wrap', method = "bits")

ggsave(plot = p, filename = "plots/apobec/ggseqlogo_apobec_signature.png",
       device = "png")


# sig.13 have the highest weighted mutation rate (18.5) and genomic span (291 Mb)

cohort_apobec[wmrate>=mean(wmrate)+2*sd(wmrate), consensusString(kmer_11), by=.(APOBEC, baselinesig)]
cohort_apobec[wmrate>=mean(wmrate)+2*sd(wmrate), summary(wmrate)]
cohort_apobec[, sd(wmrate)]

cohort_apobec_stat <-
  cohort_apobec[wmrate>=mean(wmrate)+2*sd(wmrate), .(n=.N, wmrate=sum(n_snv)/sum(tot_hg19)/unique(n_donors)*10^6,
    span=sum(as.numeric(tot_hg19)),
    consensus_motif=consensusString(kmer_11)),
    #pfm_motif=list(consensusMatrix(kmer_11, as.prob=TRUE))),
  by=.(APOBEC, baselinesig)]

cohort_apobec_stat[APOBEC=="A3A", baselinesig]
sequences <- cohort_apobec_stat[, pfm_motif]
names(sequences) <- cohort_apobec_stat[, paste(APOBEC,"-",baselinesig)]

library(ggseqlogo)
p <- ggseqlogo(data=sequences, facet='wrap', method = "bits")

ggsave(plot = p, filename = "plots/apobec/ggseqlogo_apobec_signature.png",
       device = "png")

hsdt_apobec_stat <- hsdt_apobec[rec=="5+",
            .(n=.N, wmrate=sum(n_mut_kmer)/sum(tot_hg19)/unique(n_donors)*10^6,
              span=sum(tot_hg19),
              consensus_motif=consensusString(rep(kmer_11, n_mut_kmer))),
            by=.(APOBEC,signature)]

hsdt_apobec_stat %>% dcast(APOBEC ~ signature, value.var="wmrate")
hsdt_apobec_stat %>% dcast(APOBEC ~ signature, value.var="span")
hsdt_apobec_stat %>% dcast(APOBEC ~ signature, value.var="consensus_motif")


# test difference in mutation rate between APOBEC-motifs

test2 <- hsdt_apobec[rec=="1+" & grepl("SBS2", signature) & grepl("A3",APOBEC),
            t.test(log10(mrate)~APOBEC)]
test13 <- hsdt_apobec[rec=="1+" & grepl("SBS13", signature) & grepl("A3",APOBEC),
            t.test(log10(mrate)~APOBEC)]

hsdt_apobec[rec=="1+" & grepl("SBS13", signature) & grepl("A3",APOBEC),
            .(10^(mean(log10(mrate)))), by=APOBEC]

test69 <- hsdt_apobec[rec=="1+" & grepl("SBS69", signature) & grepl("A3",APOBEC),
            t.test(log10(mrate)~APOBEC)]

Reduce("/",10^test2$estimate) # 1.26-times higher A3A mrate
Reduce("/",10^test13$estimate) # 1.27-times higher A3A mrate
Reduce("/",10^rev(test69$estimate)) # 1.33-times higher A3B mrate



# test difference in means between APOBEC groups
hsdt_apobec[rec=="1+", summary(aov(log10(mrate) ~ APOBEC))]
# there's a difference between the group
hsdt_apobec[rec=="1+", TukeyHSD(aov(log10(mrate) ~ APOBEC))]


###########
hsdt_apobec[rec=="1+" & grepl("SBS2", signature) & grepl("A3",APOBEC),
            t.test(log10(mrate) ~ APOBEC)]

# t.test of difference in log-mrate
test_s2 <- hsdt_apobec[rec=="1+" & grepl("SBS2", signature),
  t.test(log10(mrate[APOBEC=="A3A"]), log10(mrate[APOBEC=="A3B"]))]

test_s13 <- hsdt_apobec[rec=="1+" & grepl("SBS13", signature),
  t.test(log10(mrate[APOBEC=="A3A"]), log10(mrate[APOBEC=="A3B"]))]

test_s69 <- hsdt_apobec[rec=="1+" & grepl("SBS69", signature),
  t.test(log10(mrate[APOBEC=="A3A"]), log10(mrate[APOBEC=="A3B"]))]

(test_s69$conf.int)






x <- unique(hsdt_apobec[, .(n=.N, 
                            mrate = mean(mrate),
                            #mrate=sum(n_mut_kmer)/sum(tot_hg19)/n_donors*10^6, 
                            sd=sd(sum(n_mut_kmer)/sum(tot_hg19)/n_donors, na.rm = T),
                            tot_hg19=sum(tot_hg19)), by = .(APOBEC, rec, signature)])


y <- hsdt_apobec[, ks.test(mrate[APOBEC=="A3A"], mrate[APOBEC=="A3B"]), by = .(rec, signature)][, .(rec, signature, statistic, p.value)]
hsdt_apobec_stat <- merge(x, y, by = c("signature", "rec"))[, signature:=str_sub(signature, 20, -3)]

hsdt_apobec_stat[, .(signature=factor(signature, levels=paste0("S", c(2,13,69))),
      rec, APOBEC, mrate=round(mrate,1), span = tot_hg19/10^6, p.value)][order(signature)]


# table view of mutation rates
hsdt_apobec_stat[rec=="1+", .(signature=factor(signature, levels=paste0("S", c(2,13,69))),
  APOBEC, mrate=round(mrate,1))] %>% spread(key = "APOBEC", value = "mrate")

# table view of genomic span
hsdt_apobec_stat[rec=="1+", .(signature=factor(signature, levels=paste0("S", c(2,13,69))),
  APOBEC, span=paste("(",round(tot_hg19/10^6,1),"Mb )"))] %>% spread(key = "APOBEC", value = "span")

hsdt_apobec_stat[rec=="1+", .(signature=factor(signature, levels=paste0("S", c(2,13,69))),
  APOBEC, label=paste0(round(mrate,1),ifelse(p.value<0.01,"*","")," (",round(tot_hg19/10^6,1)," Mb)"))] %>%
  spread(key = "APOBEC", value = "label")

# span
# mutation rate
# patient mutability

