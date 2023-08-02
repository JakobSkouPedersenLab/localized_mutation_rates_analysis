# set environment ====
options(scipen = 9)
  
suppressPackageStartupMessages({
  library(data.table) ; library(tidyverse)
  library(Hmisc)
  library(scales)
  library(ggforce) ; library(ggseqlogo) ; library(ggrepel) ; library(ggpubr) ; library(patchwork)
})

setwd("~/GenomeDK/PCAWG/faststorage/Gustav/mutprocesses_project")

# load data ====
dt_sumstat <- fread(file = "results/stats4final_viz/summarystatistics_allsigs.tsv")
dt_distrib <- fread(file = "results/stats4final_viz/rate2span_distrib_allsigs.tsv")
dt_bitlogo <- fread(file = "results/stats4final_viz/motifdata_bitlogo_allsigs.tsv")
dt_freqlogo<- fread(file = "results/stats4final_viz/motifdata_freqlogo_allsigs.tsv")
k11c_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")
refsig_etio <- fread("data/reference_signatures/SignatureAnalyzer_SBS_etiology.tsv")
donorsigs <- fread("data/PCAWG_signatures/SA_COMPOSITE_SNV.activity_relative.WHITELIST_SET.20200226.txt")
donorsigs_count <- lapply(colnames(donorsigs)[8:67], function(x) {
  donors <- donorsigs[get(x)>0.05, "icgc_donor_id", with=FALSE]$icgc_donor_id
  ndonors <- length(donors)
  data.table(baselinesig=x, ndonors=ndonors, donors=list(donors))
}) %>% rbindlist()
setkey(donorsigs_count, baselinesig)


# test against baseline mutation rate (5.96) ====
# dt_sumstat[!is.na(tot_hg19), pval_baseline := binom.test(x = as.double(n_snv), n = as.double(tot_hg19), 
#                          p = 5.96/10^6*donorsigs_count[assigned_sig, ndonors],
#                          alternative = "greater")$p.value,
#            by=.(step, assigned_sig, region)]

# test against cohort's own mutation rate ====
# for(i in 1:nrow(dt_sumstat)) {
#   if(!is.na(dt_sumstat[i, tot_hg19])){
#     signature <- dt_sumstat[i, assigned_sig]
#     dt_sumstat[i, pval_cohort := binom.test(x = as.double(n_snv), n = as.double(tot_hg19),
#                       p = dt_sumstat[step==1 & assigned_sig==signature,
#                                      wmrate]/10^6*donorsigs_count[assigned_sig, ndonors],
#                       alternative = "greater")$p.value,
#                by=.(step, assigned_sig, region)]
#   } else {dt_sumstat[i, pval_cohort := 1]}
# }
# 

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
# compute confidence intervals for mutation rates ====
dt_sumstat[!is.na(tot_hg19),
  c("cil_prev", "ciu_prev") := {bci=Hmisc::binconf(x = as.double(n_snv),
    n = as.double(tot_hg19), alpha = 1-10^(-5), method = "wilson", return.df = T)
  list(bci$Lower/donorsigs_count[assigned_sig, ndonors]*10^6,
       bci$Upper/donorsigs_count[assigned_sig, ndonors]*10^6)
  }, by=.(step, assigned_sig, region)]


# dt_sumstat[!is.na(tot_hg19), c("pval_prev","cil_prev", "cih_prev"):={bt=binom.test(x = as.double(n_snv), n = as.double(tot_hg19),
#   p = wmrate_prev/10^6*donorsigs_count[assigned_sig, ndonors], alternative = "greater")
#   list(pval_prev = bt$p.value, cil_prev = bt$conf.int[1], cih_prev = bt$conf.int[2])},
#            by=.(step, assigned_sig, region)]


# add labels ====
dt_distrib[, label := switch(EXPR = step, "cohort","signature","hotspot",
                             region), by=.(assigned_sig, step)]
dt_sumstat[, label := switch(EXPR = step, "cohort","signature","hotspot",
                             region), by=.(assigned_sig, step)]
dt_bitlogo[, label := switch(EXPR = step, "cohort","signature","hotspot",
                             region), by=.(assigned_sig, step)]
dt_freqlogo[, label := switch(EXPR = step, "cohort","signature","hotspot",
                             region), by=.(assigned_sig, step)]
dt_sumstat[, label_rate := paste0(round(wmrate,2))]#,"\n(",round(cil_prev,2),"-",round(ciu_prev,2),")")]

# compute mutation rate increase ====
dt_sumstat[, wmrate_increase:=wmrate/wmrate_prev]
dt_sumstat[, label_increase:=paste0(format(wmrate_increase,digits = 1),"x")]
dt_sumstat[step==1, label_increase:=""]

# p-value adjustment ====
dt_sumstat[, pval_bonf_prev := pval_prev*nrow(dt_sumstat)]
dt_sumstat[pval_bonf_prev>1, pval_bonf_prev := 1]
dt_sumstat[, ppval_bonf_prev := -log10(pval_bonf_prev)]
dt_sumstat[is.infinite(ppval_bonf_prev), ppval_bonf_prev:=dt_sumstat[is.finite(ppval_bonf_prev), max(ppval_bonf_prev)]]


# define plotting ====
plot_linegraph <- function(sig="SBS", ppval_threshold=10, single_region=TRUE, keep.strip=FALSE){
  
  theme_linegraph <- theme(strip.background = element_blank(),
                           strip.text = element_blank())
  if(keep.strip) theme_linegraph <- theme(strip.placement = "outside")
  
  scientific_10 <- function(x) parse(text=paste0(gsub("e", " %*% 10^", scales::scientific_format()(x)),"*'%'"))
  axis_label <- function(x) ifelse(x*100<0.01 & x!=0, scientific_10(x*100), paste0(signif(x*100,2),"%"))
  
  
  d <- dt_sumstat[grepl(sig, assigned_sig)]
  if(single_region) d <- d[step<4 | (step==4 & wmrate==max(wmrate[step==4], na.rm=T))]
  
  d[!(ppval_bonf_prev<=ppval_threshold & step==4)] %>%
    ggplot(aes(tot_hg19, wmrate)) +
    # plotting data
    # draw lines
    geom_line(data=d, color="grey") +
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
    geom_text(data=d[step%in%2:3], aes(label=label), vjust=-1, size=3) +
    geom_text(data=d[step==4 & !ppval_bonf_prev<=ppval_threshold,
      .(region=region[which.max(wmrate)],wmrate=wmrate[which.max(wmrate)], tot_hg19=tot_hg19[which.max(wmrate)]),
      by=.(assigned_sig)], aes(label=region), vjust=-1, size=3) +
    ## label with mutation rate increase
    geom_text(data=d[step!=4], aes(label=label_increase), vjust=1, hjust=0, size=4) +
    geom_text(data=d[step==4 & !ppval_bonf_prev<=ppval_threshold,
                     .(label_increase=label_increase[which.max(wmrate)], wmrate=wmrate[which.max(wmrate)], tot_hg19=tot_hg19[which.max(wmrate)]),
                     by=.(assigned_sig)], aes(label=label_increase), vjust=1, hjust=0, size=4) +
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


# IMPORTANT. CHOOSE TO SET Y-AXIS LIMITS FOR MAIN PLOT (y_range_fun) ====
plot_barchart <- function(sig="SBS", astep=NA, keep.strip=FALSE){
  d <- dt_distrib[grepl(sig, assigned_sig)]
  x_breaks <- 10^(-1:5)
  y_breaks <- d[, c(0, max(tot_hg19)/2, max(tot_hg19)), by=label]$V1
  y_breaks_sec <- y_breaks/sum(k11c_dt$tot_hg19)
  
  theme_barchart <- theme(strip.background.x = element_blank(),
                           strip.text.x = element_blank())
  if(keep.strip) theme_barchart <- theme(strip.placement = "outside")

  
  mrate_cohort <- dt_sumstat[grepl(sig, assigned_sig) & step==1, wmrate]
  max_mrate_region <- dt_sumstat[grepl(sig, assigned_sig) &
               step==4, region[which.max(wmrate)]]
  
  if(!is.na(astep)) d <- d[step==astep]
  
  label_lvls <- rev(c("cohort","signature","hotspot", max_mrate_region))
  
  scientific_10 <- function(x) parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  axis_label <- function(x) ifelse(x<0.1 & x!=0, scientific_10(x), signif(x,2))
  y_breaks_fun <- function(x) {
    xmax <- max(x, na.rm = TRUE)
    xmax_log10 <- floor(log10(xmax))
    xmax_round <- round(xmax, -xmax_log10)
    if(xmax_log10==0) {breaks <- c(0, xmax/2, xmax) 
    } else{ breaks <- c(0, xmax_round/2, xmax_round) }
    return(breaks)
    }
  y_range_fun <- function(x) {
    xmax <- max(x, na.rm = TRUE)
    xmax_log10 <- floor(log10(xmax))
    xmax_round <- round(xmax, -xmax_log10)
    if(xmax_round==0) xmax_round <- round(xmax, -(xmax_log10-1))
    limits <- c(0, xmax_round*1.25)
    return(limits)
  }
  
  genome_scale <- function(x){
    ifelse(x==0, x, 
      #ifelse(x/10^9>=.1, paste(x/10^9, "Gb"),
        ifelse(x/10^6>=.1, paste(x/10^6,"Mb"),
          ifelse(x/10^3>=.1,paste(x/10^3, "kb"),
                 paste(x, "bp"))))#)
    }
  
  # return(d[region=="Genome-wide" | region==max_mrate_region,
  #          y_breaks_fun(tot_hg19), by=label])
  
  d[region=="Genome-wide" | region==max_mrate_region] %>%
    ggplot(aes(x=mrate, y=tot_hg19, fill=log10(mrate))) +
    geom_col(col="black", width = 0.09) +
    geom_text(data = d[region=="Genome-wide" | region==max_mrate_region,
                .(mrate=.1, tot_hg19=Inf), by=label], mapping = aes(label=label),
              vjust=1, hjust=0, size=3) +
    scale_x_log10(name="Mutation rate (SNV/patient/Mb)",
                  breaks=x_breaks,
                  labels=trans_format("log10", math_format(10^.x))) +
                    #function(x) ifelse(x>1, prettyNum(x, big.mark=","), signif(x, 2))) +
    scale_y_continuous(
      name="Genomic span", 
      breaks = y_breaks_fun, #limits = y_range_fun,
      labels = genome_scale, expand=c(0,0),
      sec.axis = sec_axis(name = "Fraction (%)", trans = ~./sum(k11c_dt$tot_hg19),
                          breaks = pretty_breaks(n=2), labels = function(x) axis_label(x*100)))+
    scale_fill_gradient2(name=expression(log[10]~rate),
      #name="log10(mutation rate)",
      low = "blue", high = "brown",
                         midpoint = log10(mrate_cohort), limits=c(-1.1,5.1),
                         breaks=c(-1,round(log10(mrate_cohort),1),5),
                         labels=function(x) format(x, drop0trailing = TRUE),
                         guide=guide_colorbar(direction = "horizontal",
                                              barheight = unit(2,"mm"),
                                              barwidth = unit(1.5, "cm"),
                                              ticks.colour = "black",
                                              title.position = "top",
                                              title.theme = element_text(size=10, hjust = .5),
                                              frame.colour = "black")) +
    # facet_grid(paste(step, region)~str_extract(assigned_sig, "SBS[0-9a-z]+"), scales="free_y") +
    facet_wrap(factor(label,label_lvls)~., #str_extract(assigned_sig, "SBS[0-9a-z]+")
               scales="free_y", strip.position = "left", ncol=1) + # grid: switch = "y"
    coord_cartesian(xlim=range(x_breaks)) +
    theme_pubr() +
    #theme_barchart +
    theme(text = element_text(size=12),
          axis.title.y.right = element_text(angle=90), strip.placement = "inside",
          strip.background = element_blank(), strip.text = element_blank(), 
          panel.spacing = unit(5, "mm"),
          legend.position = c(.8,.15), legend.background = element_blank())
}

plot_bitlogo <- function(heights){
  ggseqlogo(heights, method = "custom", seq_type = "dna", ncol = 1)+
    scale_x_continuous(name="Position (relative to mutation)", expand = c(0,0),
                       breaks = 1:11, labels = 1:11-6, position = "bottom") +
    scale_y_continuous(
      name="Information content (bits)", expand = c(0, 0),
      breaks = c(0, .5, 1, 1.5, 2), labels = c(0, "", 1, "", 2)) +
    geom_hline(yintercept = c(0)) +
    geom_text(inherit.aes = FALSE, data=data.frame(x=1.1,y=1.9,seq_group=names(heights)),
              aes(x=x, y=y, label = seq_group), hjust=0, size=3) +
    coord_cartesian(ylim = c(0, 2),  clip = "off") +
    facet_wrap(~seq_group, ncol = 1) +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12),
          axis.line.y = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      strip.text = element_blank(), panel.spacing = unit(5, "mm"))
      #strip.text = element_text(hjust=0, vjust=-1))
}
plot_bitlogos_facets <- function(sig="SBS", genomicregion=TRUE){
  max_mrate_region <- dt_sumstat[grepl(sig, assigned_sig) & step==4, region[which.max(wmrate)]]
  subsetnames <- rev(c("cohort","signature","hotspot",max_mrate_region))
  if(!genomicregion) subsetnames <- rev(c("cohort","signature","hotspot"))
  heights <- lapply(subsetnames, function(x){
    d <- dt_bitlogo[grepl(sig, assigned_sig) & label==x &
                      (region=="Genome-wide" | region==max_mrate_region),
                    c("base", paste0("V", 1:11)), with=FALSE]
    m <- as.matrix(d, rownames = "base")
    colnames(m) <- NULL
    return(m)
  })
  names(heights) <- subsetnames
  # return(heights)
  p_bitlogo <- plot_bitlogo(heights)
  return(p_bitlogo)
  }

# MAIN FIGURE
# brute force merge plots with patchwork ====
unique(dt_sumstat$assigned_sig)[c(20,7,38,48)]
p_patch <- plot_linegraph("BI_COMPOSITE_SNV_SBS17b_P", single_region = T) +
  plot_barchart("BI_COMPOSITE_SNV_SBS17b_P") +
  plot_bitlogos_facets("BI_COMPOSITE_SNV_SBS17b_P") +
  plot_linegraph("BI_COMPOSITE_SNV_SBS7a_S", single_region = T) +
  plot_barchart("BI_COMPOSITE_SNV_SBS7a_S") +
  plot_bitlogos_facets("BI_COMPOSITE_SNV_SBS7a_S") +
  plot_linegraph("BI_COMPOSITE_SNV_SBS62_S", single_region = T) +
  plot_barchart("BI_COMPOSITE_SNV_SBS62_S") +
  plot_bitlogos_facets("BI_COMPOSITE_SNV_SBS62_S") +
  plot_linegraph("BI_COMPOSITE_SNV_SBS72_P", single_region = T) +
  plot_barchart("BI_COMPOSITE_SNV_SBS72_P") +
  plot_bitlogos_facets("BI_COMPOSITE_SNV_SBS72_P") +
  plot_annotation(tag_levels = "a") +
  plot_layout(ncol = 3, nrow = 4, byrow = TRUE)

ggsave(plot=p_patch, device="pdf", width = 30, height = 40, units = "cm", useDingbats=FALSE,
       filename="plots/factorization_steps/v4_new/fig5_patch.pdf")


# EXTENDED DATA
# individual plots ====
for (SBS in unique(dt_sumstat$assigned_sig)) {
  p <- plot_linegraph_v2(SBS, single_region = T) +
    plot_barchart(SBS) +
    plot_bitlogos_facets(SBS) +
    plot_annotation(title = paste("Signature", str_extract(SBS, "(?<=SBS)[0-9a-z]+"))) +
    plot_layout(nrow=1, widths=c(1,1,1))
  # plot(p)
  ggsave(plot=p, device="pdf", width = 30, height = 10, units = "cm",
         filename=paste0("plots/factorization_steps/v4_new/linegraph_barchart_logo/sbs",
                         str_extract(SBS, "(?<=SBS)[0-9a-z]+"), "_fig5.pdf"))
}



# extended data figures ====
# re-run 05-01-2023.
for (i in 1:9) {
  start <- seq(1, 60, by=6)[i]
  end <- seq(1, 60, by=6)[i]+5
  signatures <- unique(dt_sumstat$assigned_sig)[start:end]
  sigs <- str_extract(signatures, "(?<=SBS)[0-9a-z]+")
  sigs_collapse <- paste0(sigs, collapse = "_")
  filename_var <- paste0("plots/factorization_steps/v4_new/supfig_v2/extended_data_sbs_",
                         sigs_collapse,".pdf")
  print(sigs_collapse)
  p <- ggarrange(plotlist =
              lapply(signatures, function(SBS){
                p <- plot_linegraph_v2(SBS, single_region = T) + plot_barchart(SBS) +
                  plot_bitlogos_facets(SBS) + plot_layout(ncol=3, byrow=TRUE)
                return(p)}),
            ncol = 1)
  ggsave(plot = p, device = "pdf", width=30, height=60, units="cm",
         filename=filename_var)
}

signatures <- unique(dt_sumstat$assigned_sig)[55:57]
sigs <- str_extract(signatures, "(?<=SBS)[0-9a-z]+")
sigs_collapse <- paste0(sigs, collapse = "_")
filename_var <- paste0("plots/factorization_steps/v4_new/extended_data_sbs_",
                       sigs_collapse,".pdf")

p <- ggpubr::ggarrange(plotlist =
                 lapply(signatures, function(SBS){
                   p <- plot_linegraph(SBS, single_region = T) + plot_barchart(SBS) +
                     plot_bitlogos_facets(SBS) + plot_layout(ncol=3, byrow=TRUE)
                   return(p)}),
               ncol = 1)
ggsave(plot = p, device = "pdf", width=30, height=30, units="cm",
       filename=filename_var)


# 
plot_bitlogos_facets("BI_COMPOSITE_SNV_SBS17b_P")
dt_bitlogo[step%in%1:3 & grepl("SBS17b", assigned_sig)]




k11_sig17b_hs <- fread("results/kplogo/sbs17b_logo/input/foreground_sbs17b_3_hotspot5+.tsv", header=FALSE, col.names = "kmer")
k11_sig17b_sig <- fread("results/kplogo/sbs17b_logo/input/foreground_sbs17b_2_11mer.tsv", header=FALSE, col.names = "kmer")
k11_sig17b_coh <- fread("results/kplogo/sbs17b_logo/input/foreground_sbs17b_1_cohort.tsv", header=FALSE, col.names = "kmer")
freqmatrix <- lapply(list(k11_sig17b_hs, k11_sig17b_sig, k11_sig17b_coh),
       function(x) Biostrings::consensusMatrix(x = x$kmer, as.prob = TRUE))
ggseqlogo(freqmatrix, method = "p")

plot_bitlogos_facets("SBS17b")

dt_freqlogo
dt_bitlogo


# extended logos: freqlogo + bitslogo ====

plot_freqlogo <- function(sig, genomicregion=TRUE){
  max_mrate_region <- dt_sumstat[grepl(sig, assigned_sig, ignore.case = TRUE) &
                                   step==4, region[which.max(wmrate)]]
  subsetnames <- rev(c("cohort","signature","hotspot",max_mrate_region))
  if(!genomicregion) subsetnames <- rev(c("cohort","signature","hotspot"))
  heights <- lapply(subsetnames, function(x){
    d <- dt_freqlogo[grepl(sig, assigned_sig) & label==x &
          (region=="Genome-wide" | region==max_mrate_region),
                    c("base", paste0("V", 1:11)), with=FALSE]
    m <- as.matrix(d, rownames = "base")
    colnames(m) <- NULL
    return(m)
  })
  names(heights) <- subsetnames
  #return(heights)
  ggseqlogo(heights, method = "custom", seq_type = "dna", ncol = 1)+
    scale_x_continuous(name="Position (relative to mutation)", expand = c(0,0),
                       breaks = 1:11, labels = 1:11-6, position = "bottom") +
    scale_y_continuous(
      name="Letter frequency (%)", expand = c(0, 0),
      breaks = c(0, .25, .5, .75, 1), labels = c(0, "", 50, "", 100) ) +
    geom_hline(yintercept = c(0)) +
    geom_text(inherit.aes = FALSE, data=data.frame(x=1.1,y=1.1,seq_group=names(heights)),
              aes(x=x, y=y, label = seq_group), hjust=0, size=3) +
    coord_cartesian(ylim = c(0, 1.1),  clip = "off") +
    facet_wrap(~seq_group, ncol = 1) +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12),
          axis.line.y = element_line(color = "black"),
          axis.ticks.y = element_line(color = "black"),
          strip.text = element_blank(), panel.spacing = unit(5, "mm"))
  #strip.text = element_text(hjust=0, vjust=-1))
}

plot_freqlogo("SBS17b", genomicregion = F) +
  plot_bitlogos_facets("SBS17b", genomicregion = F) + plot_layout(ncol=2)


for(SBS in unique(dt_sumstat$assigned_sig)) {
  p <- plot_freqlogo(SBS, genomicregion = F) + ggtitle(label = "Frequency") + 
    plot_bitlogos_facets(SBS, genomicregion = F) + ggtitle(label = "Information") +
    plot_annotation(title = paste("Signature", str_extract(SBS, "(?<=SBS)[0-9a-z]+"))) +
    # tag_levels = list(c("Frequency", "Information"))
    plot_layout(nrow=1, widths=c(1,1))
  # plot(p)
  ggsave(plot=p, device="pdf", width = 16, height = 13, units = "cm", # height=7.5
         filename=paste0("plots/factorization_steps/v4_new/extended_logo/sbs",
                         str_extract(SBS, "(?<=SBS)[0-9a-z]+"), "_freq_bits.pdf"))
}




# genomic regions
genome_scale <- function(x){
  ifelse(x==0, x,
         ifelse(x/10^6>=.1, paste(format(round(x/10^6,1), big.mark=",", digits = 1, nsmall=0), "Mb"),
                ifelse(x/10^3>=.1, paste(round(x/10^3,1), "kb"), paste(x, "bp"))))
}
# create labels for supplementary logo figures
dt_sumstat[, `:=`(span = genome_scale(as.numeric(tot_hg19)),
                  span_frac = tot_hg19/sum(k11c_dt$tot_hg19))]

dt_sumstat[assigned_sig %in% dt_sumstat[(step>=3 & ppval_bonf_prev>=2), unique(assigned_sig)]][,
                                                                                               .(assigned_sig, step, label,
                                                                                                 span=paste0(span," (", ifelse(span_frac*100>=0.01, round(span_frac*100,2), "<0.01"), "%)"),
                                                                                                 mut.rate=paste0(round(wmrate,2), ifelse(label_increase=="","",paste0(" (",gsub(" ","",label_increase),")"))),
                                                                                                 p.adjust.log10=round(ppval_bonf_prev))][grepl("SBS17b", assigned_sig) & !(step==4 & p.adjust.log10<2)]
dt_sumstat[step!=4 | (step==4 & ppval_bonf_prev>2) & label!="", .N]

dt_sumstat[grepl("SBS17b", assigned_sig) & step==4, label]

dt_freqlogo[grepl("SBS17b", assigned_sig) & step==4,
            c("base", paste0("V", 1:11)), with=FALSE]

#lvls <- str_sub(dt_sumstat[label!="", .(unique(paste(step, label)))][order(V1)]$V1, 3, -1)

lab_dt <- dt_sumstat[assigned_sig %in% dt_sumstat[(step>=3 & ppval_bonf_prev>=2), unique(assigned_sig)] &
                       step!=4 | (step==4 & ppval_bonf_prev>2)][,
    .(assigned_sig, step, label,
      span=paste0(" span ",span," (", ifelse(span_frac*100>=0.01, round(span_frac*100,2), "<0.01"), "%)"),          
      mut.rate=paste0(" rate ",round(wmrate,2)),
      mut.rate.increase=ifelse(label_increase=="","",paste0(" factor increase ",gsub(" ","",label_increase))),
      p.adjust.log10=ifelse(label_increase=="","", paste0(" -log10(p) = ", round(ppval_bonf_prev))))][order(assigned_sig, step, label)]

lab_dt[, info_label:= paste0(label,"\n", span, "\n", mut.rate, " SNV/Mb/patient \n",mut.rate.increase,
                             "\n", p.adjust.log10)]

ggplot(lab_dt[grepl("SBS17b", assigned_sig)][, i:=.N:1][,], aes(x=0, y=i, label=info_label))+
  geom_text(hjust=0, vjust=1) + coord_cartesian(ylim = c(-.5,6)) + theme_void()


plot_regions_freq_bits <- function(sig){
  region_names <- c(dt_sumstat[grepl(sig, assigned_sig) & 
                                 step!=4 & label!="", label],
                    dt_sumstat[grepl(sig, assigned_sig) & 
                               step==4 & ppval_bonf_prev>2
                               & label!="", sort(label)])
  
  #lab_dt[grepl("SBS17b", assigned_sig) & !(step==4 & p.adjust.log10>2)]
  
  #region_names <- lab_dt[grepl(sig, assigned_sig), info_label]

  sbs <- str_extract(sig, "SBS[0-9a-z]+")
  # extract freqlogo heights
  freq_heights <- lapply(region_names, function(x){
    d <- dt_freqlogo[grepl(sig, assigned_sig) & label==x,
                     c("base", paste0("V", 1:11)), with=FALSE]
    m <- as.matrix(d, rownames = "base")
    colnames(m) <- NULL
    return(m)
  })
  names(freq_heights) <- region_names
  
  # extract bitlogo heights
  bit_heights <- lapply(region_names, function(x){
    d <- dt_bitlogo[grepl(sig, assigned_sig) & label==x,
                    c("base", paste0("V", 1:11)), with=FALSE]
    m <- as.matrix(d, rownames = "base")
    colnames(m) <- NULL
    return(m)
  })
  names(bit_heights) <- region_names
  
  # plot freqlogo
  p_freq <- ggseqlogo(freq_heights, method = "custom", seq_type = "dna", ncol = 1)+
    scale_x_continuous(name="Position (relative to mutation)", expand = c(0,0),
                       breaks = 1:11, labels = 1:11-6, position = "bottom") +
    scale_y_continuous(
      name="percent (%)", expand = c(0, 0),
      breaks = c(0, .25, .5, .75, 1), labels = c(0, "", 50, "", 100) ) +
    geom_hline(yintercept = c(0)) +
    # geom_text(inherit.aes = FALSE, data=data.frame(x=1.1,y=1.1,seq_group=names(freq_heights)),
    #           aes(x=x, y=y, label = seq_group), hjust=0, size=3) +
    ggtitle(label = "Letter frequency") +
    coord_cartesian(ylim = c(0, 1),  clip = "off") +
    facet_wrap(~seq_group, ncol = 1) +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12),
          axis.line.y = element_line(color = "black"),
          axis.ticks.y = element_line(color = "black"),
          strip.text = element_blank(), #text(angle = 180),
          panel.spacing = unit(10, "mm"))
  # plot bitslogo
  p_bits <- ggseqlogo(bit_heights, method = "custom", seq_type = "dna", ncol = 1)+
    scale_x_continuous(name="Position (relative to mutation)", expand = c(0,0),
                       breaks = 1:11, labels = 1:11-6, position = "bottom") +
    scale_y_continuous(
      name="bits", expand = c(0, 0),
      breaks = c(0, .5, 1, 1.5, 2), labels = c(0, "", 1, "", 2)) +
    geom_hline(yintercept = c(0)) +
    # geom_text(inherit.aes = FALSE, data=data.frame(x=1.1,y=1.9,seq_group=names(bit_heights)),
              # aes(x=x, y=y, label = seq_group), hjust=0, size=3) +
    ggtitle(label = "Information content") +
    coord_cartesian(ylim = c(0, 2),  clip = "off") +
    facet_wrap(~seq_group, ncol = 1) +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12),
          axis.line.y = element_line(color = "black"),
          axis.ticks.y = element_line(color = "black"),
          strip.text = element_blank(), panel.spacing = unit(10, "mm"))
  
  labsig_dt <- lab_dt[grepl(sig, assigned_sig)][, i:=(.N:1)]
  p_labels <- ggplot(labsig_dt, aes(x=0, y=i, label=info_label))+
    geom_text(hjust=0, vjust=.5, size=3) + coord_cartesian(xlim=c(-.1,.3),ylim = c(1,max(labsig_dt$i)), clip = "off") + theme_void()
  #size=ifelse(max(labsig_dt$i)>6, 3, 4)
  return(p_labels + p_freq + p_bits + plot_annotation(title = sbs) + plot_layout(ncol=3))
}


plot_regions_freq_bits(unique(dt_sumstat$assigned_sig)[48])


dt_sumstat[(step>=3 & ppval_bonf_prev>=5 & wmrate_increase>=2), .(regions=paste0(region, " (", gsub(" ", "", label_increase),")", collapse=", ")), by=assigned_sig]








# this was run 29-08-2022 ====

setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project/")
for (signature in dt_sumstat[(step>=3 & ppval_bonf_prev>=2), unique(assigned_sig)]){
  n <- dt_sumstat[grepl(signature, assigned_sig) & 
               step==4 & ppval_bonf_prev>2
             & label!="", .N] + 3
  p_freq_bits <- plot_regions_freq_bits(signature)
  sbs <- tolower(str_extract(signature, "SBS[0-9a-z]+"))
  ggsave(plot=p_freq_bits, device = "png", width = 7, height = 4/3*n, dpi = 600,
         filename=paste0("plots/factorization_steps/v4_new/extended_logo/genomicregion/",sbs,"_genomicregion_bits_freq.png"))
}








# pLogo ====
#m <- function(dt) matrix(as.numeric(as.matrix(dt[,c(1:11)])), nrow = 4, ncol=11, byrow=F, dimnames = list(c("A","C","G","T")))
m <- function(dt) matrix(as.numeric(as.matrix(dt)), nrow = 4, ncol=11, byrow=F, dimnames = list(c("A","C","G","T")))
dt_sumstat[, sbs:=tolower(str_extract(assigned_sig, "(?<=SNV_)SBS[0-9a-z]+"))]
dt_sumstat[, dir_label1 := switch(EXPR = step, "cohort","signature","hotspot","genomicregion"), by=.(assigned_sig, step)]
#dt_sumstat[, dir_label2 := switch(EXPR = step, "cohort","11mer","hotspot5", region), by=.(assigned_sig, step)]

dt_pwmfiles <- dt_sumstat[, .(filename=unique(dir(path=paste0("results/kplogo/",sbs,"_logo/output/",dir_label1),
                 pattern=paste0(sbs,"_",step,"_[0-9A-Za-z//_]*s.pwm.txt"), full.names=TRUE))),
           by=.(assigned_sig, step)]

l_pfm <- lapply(1:nrow(dt_pwmfiles), function(i) {
  filename_var <- dt_pwmfiles[i, filename]
  d <- fread(filename_var, drop = "V12")[,
    `:=`(sig=dt_pwmfiles[i]$assigned_sig, step=dt_pwmfiles[i]$step,
         filename=dt_pwmfiles[i]$filename)]
  return(d)
  # o <- m(d)
})

dt_pwm <- rbindlist(l_pfm)

# pwm_sbs17b_cohort <- fread("results/kplogo/sbs17b_logo/output/cohort/sbs17b_1_cohort_s.pwm.txt")
dt_sbs17b_cohort <- fread("results/kplogo/sbs17b_logo/output/cohort/sbs17b_1_cohort_p.most.significant.each.position.txt")

matrix(pwm_sbs17b_cohort[,c(1:11)], nrow = 4, ncol=11, byrow = )

lapply(1:3, function(i) m(dt_pwm[grepl("SBS17b", sig) & step==i, c(1:11)]))

ggplot() +
  geom_logo(data=lapply(1:3, function(i) m(dt_pwm[grepl("SBS17b", sig) & step==i, c(1:11)])),
            method = "custom") +
  geom_hline(yintercept = 0) +
  scale_x_continuous(name="Position (relative to mutation)", expand = c(0,0),
                     breaks = 1:11, labels = 1:11-6, position = "bottom") +
  scale_y_continuous(name="Z-score statistic", 
                     breaks = pretty_breaks(n=7),
                     labels = function(x) format(x, big.mark = ",")) +
  facet_grid(seq_group~., scales="free_y") +
  theme_logo() +
  theme(axis.line = element_line(), axis.ticks.y = element_line())

