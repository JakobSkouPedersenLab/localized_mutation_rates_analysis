# author Gustav Poulsgaard
# created mar 24, 2020
# modified apr 29, 2021
# purpose 
# (plots for figure 1)

# define environment ====
## set options
options(scipen = 999)
## set working directory
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")
## load packages
library(data.table)
library(tidyverse)
library(ggpubr)
library(scales)

# load data ====
k11c_dt_all <- fread(file = "results/11mer_count_hg19_PCAWG_v3.tsv")
k11c_dt <- k11c_dt_all[tot_hg19!=0]
rm(k11c_dt_all)
k11c_dt[, `:=`(tot_mut = n_snv, kmer = kmer_11)]
#k11c_dt <- k11c_dt_all#[tot_mut!=0] # gdk
#k11c_dt[, `:=`(triN = NULL, freq_hg19=NULL, center=NULL, mutrate=NULL)]

# modify data ====
## define constant (used instead of zero)
cpseudo = 0.5
## 
k11c_dt[, mrate := tot_mut/tot_hg19/2583*10^6]
k11c_dt[, tot_mut_nonzero := ifelse(tot_mut==0, cpseudo, tot_mut)]
k11c_dt[, mrate_nonzero := ifelse(mrate==0, 0.1, mrate)]

# visualize data ====
## define cases
cases <- c("AAAACTTACGG", "GTATTCCGTAA") # new (18-05-2021)
#cases <- c("AAAACTTATCC", "AACTTCCGGGT")

## set axis limits and breaks
x_limits_ref = c(1, 4674610)
x_breaks_ref = c(1, 10^2, 10^4, 10^6)

y_limits_mut =c(cpseudo*0.7, 36897) # include non-mutated
#y_limits_mut = c(0.7, 36897) # exclude non-mutated
y_breaks_mut = c(cpseudo, 1, 10, 100, 10^3, 10^4)

x_breaks_rate = 10^c(0, -1, -2, -3)
breaks_rate_patient = c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7)

## define constant (density to the power of cdens)
cdens = 1/8

# optimized for exclusion of non-mutated
smoothscat_test2 <- function(dt, legend_shape = guide_legend(), axis_text = element_text(), unadj_mutrate = FALSE){
  if(!is.data.table(dt)) dt <- as.data.table(dt)
  mutrate_labels <- paste0(c(comma(10^3), 10^2, 10^1, 10^0, 10^(-1)), " SNV/patient/Mb")
  
  ggplot(dt, aes(tot_hg19, tot_mut)) +
    ## background grid
    #geom_hline(yintercept = c(2.583, 25.83, 258.3, 2583, 25830), alpha = .05) + 
    #geom_vline(xintercept = c(1, 10, 100, 10^3, 10^4, 10^5, 10^6), alpha = .05) +
    ## modify x axis
    scale_x_continuous(name = "k-mer genomic occurrence", 
                       trans = "log10",
                       expand = c(0, 0),
                       limits = x_limits_ref,
                       breaks = x_breaks_ref,
                       labels = comma(x_breaks_ref, accuracy = 1)) +
    ## modify y axis
    scale_y_continuous(name = "k-mer mutation count",
                       trans = "log10", 
                       expand = c(0, 0),
                       limits = y_limits_mut, # max(k11c_dt$tot_mut)
                       breaks = y_breaks_mut,
                       labels = function(x) comma(x, accuracy = 1)) +
    ## modify filling
    scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256), guide = FALSE) +
    ## density by tiles
    stat_density2d(geom = "tile", aes(fill= ..density..^cdens, alpha = 1), contour = FALSE) + 
    ## points at cases
    geom_point(inherit.aes = FALSE,
               data = dt[kmer %in% cases][, kmer := factor(kmer, levels = cases)], 
               mapping = aes(x = tot_hg19, y = tot_mut, shape = kmer), 
               col = "black", alpha = 1, size = 3) +
    ## modify shapes
    scale_shape_manual(name = "", values = setNames(c(0, 1), cases), 
                       guide = legend_shape) +
    ## diagonal dashed lines
    geom_abline(intercept = log10(breaks_rate_patient*2583), slope = 1, alpha = 0.25, linetype = "dashed", col = "navy") +
    ## text per patient mutation rate
    annotate("text", 
             x = 10^c(3, 4, 5, 6, 6),
             y = 2*(2583*10^c(0, 0, 0, 0, -1)),
             label =  mutrate_labels,
             col = "black", angle = 45) +
    # show the total number of 11-mers included
    #annotate("text", x=1.5, y=10^4, hjust = 0, vjust = 0.1,
    #         label = dt[, paste("n =", comma(.N))]) +
    theme_pubr() +
    guides(alpha = FALSE) +
    theme(legend.position = c(.15, .85),
          legend.background = element_blank(),
          legend.direction = "vertical",
          legend.margin = margin(6, 6, 6, 6),
          axis.text = axis_text, axis.title = axis_text)
}

# optimized for inclusion of non-mutated
smoothscat_test3 <- function(dt, legend_shape = guide_legend(), axis_text = element_text(), unadj_mutrate = FALSE){
  if(!is.data.table(dt)) dt <- as.data.table(dt)
  mutrate_labels <- paste0(c(comma(10^3), 10^2, 10^1, 10^0, 10^(-1)), " SNV/patient/Mb")
  
  ggplot(dt, aes(tot_hg19, tot_mut_nonzero)) +
    ## background grid
    #geom_hline(yintercept = c(2.583, 25.83, 258.3, 2583, 25830), alpha = .05) + 
    #geom_vline(xintercept = c(1, 10, 100, 10^3, 10^4, 10^5, 10^6), alpha = .05) +
    ## modify x axis
    scale_x_continuous(name = "k-mer genomic occurrence", 
                       trans = "log10",
                       expand = c(0, 0),
                       limits = x_limits_ref,
                       breaks = x_breaks_ref,
                       labels = comma(x_breaks_ref, accuracy = 1)) +
    ## modify y axis
    scale_y_continuous(name = "k-mer mutation count",
                       trans = "log10", 
                       expand = c(0, 0),
                       limits = y_limits_mut, # max(k11c_dt$tot_mut)
                       breaks = y_breaks_mut,
                       labels = function(x) comma(x, accuracy = 1)) +
    ## modify filling
    scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256), guide = FALSE) +
    ## density by tiles
    stat_density2d(geom = "tile", aes(fill= ..density..^cdens, alpha = 1), contour = FALSE) + 
    ## points at cases
    geom_point(inherit.aes = FALSE,
               data = dt[kmer %in% cases][, kmer := factor(kmer, levels = cases)], 
               mapping = aes(x = tot_hg19, y = tot_mut, shape = kmer), 
               col = "black", alpha = 1, size = 3) +
    ## modify shapes
    scale_shape_manual(name = "", values = setNames(c(0, 1), cases), 
                       guide = legend_shape) +
    ## diagonal dashed lines
    geom_abline(intercept = log10(breaks_rate_patient*2583), slope = 1, alpha = 0.25, linetype = "dashed", col = "navy") +
    ## text per patient mutation rate
    annotate("text", 
             x = 10^c(3, 4, 5, 6, 6),
             y = 2*(2583*10^c(0, 0, 0, 0, -1)),
             label =  mutrate_labels,
             col = "black", angle = 45) +
    # show the total number of 11-mers included
    #annotate("text", x=1.5, y=10^4, hjust = 0, vjust = 0.1,
    #         label = dt[, paste("n =", comma(.N))]) +
    theme_pubr() +
    guides(alpha = FALSE) +
    theme(legend.position = c(.15, .85),
          legend.background = element_blank(),
          legend.direction = "vertical",
          legend.margin = margin(6, 6, 6, 6),
          axis.text = axis_text, axis.title = axis_text)
}

## smooth scatter
#p_smoothscat <- smoothscat_test2(k11c_dt)
p_smoothscat <- smoothscat_test3(k11c_dt)
#p_smoothscat

## histograms
### count mutations
p_hist_mut <- ggplot(k11c_dt, aes(tot_mut_nonzero)) + # tot_mut_nonzero
  ## histogram with 30 bins
  geom_histogram(bins = 30, col = "black", fill = "grey90") +
  ## modify x axis (vertical when rotated)
  scale_x_continuous(trans = "log10",
                     expand = c(0, 0),
                     limits = y_limits_mut,
                     breaks = y_breaks_mut, #c(cpseudo*.85, y_breaks_mut),
                     labels = y_breaks_mut) + #c(0, y_breaks_mut)) +
  ## modify y axis (horizontal when rotated)
  scale_y_continuous(name = "mutation count\n(count x1,000)",
                     expand = c(0, 0),
                     limits = c(0, 3.1*10^5),
                     breaks = c(1, 2, 3)*10^5,
                     labels = c("100", "200", "300")) +
  ## modify theme
  theme_pubr() + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  rotate()

### count occurrence in reference
p_hist_ref <- ggplot(k11c_dt, aes(tot_hg19)) +
  ## histogram with 30 bins
  geom_histogram(bins = 30, col = "black", fill = "grey90") +
  ## modify x axis 
  scale_x_continuous(trans = "log10", 
                     expand = c(0, 0),
                     #limits = x_limits_ref, 
                     breaks = x_breaks_ref) +
  ## modify y axis 
  scale_y_continuous(name = "genomic occurrence\n(count x1,000)",
                     breaks = c(1, 2, 3)*10^5,
                     labels = c("100", "200", "300"),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = x_limits_ref, ylim = c(0, 4*10^5)) +
  ## modify theme
  theme_pubr() + theme(axis.text.x = element_blank(), axis.title.x = element_blank())



## arrange
library(patchwork)
p_arranged <- p_hist_ref + plot_spacer() + p_smoothscat + p_hist_mut + 
  patchwork::plot_layout(ncol = 2, nrow = 2, byrow = T, widths = c(2,1), heights = c(1,2))

ggsave(plot = p_arranged, device = "pdf",
       filename = "plots/kmer_general_mutrates/smoothscatter_histograms_11mer_refcount_pcawgmuts_v8.pdf",
       width = 10.1, height = 7.5, #1.34 ratio
       units = "in", useDingbats = FALSE)
rm(p_hist_ref,p_hist_mut,p_smoothscat, p_arranged)

###
### mutation rate
# include non-mutated (nonzero)
hist_rate <- function(dt, y_breaks=NULL, ylimits=NULL){
  mean_mutrate <- dt[, mean(ifelse(is.na(mrate), 0, mrate))]

  mutrate_labels <- c(comma(10^3), 10^2, 10^1, 10^0, 10^(-1))
  
  ggplot(data = dt, mapping = aes(mrate_nonzero)) + # mrate+0.09
    ## histogram with 30 bins
    geom_histogram(bins = 30, mapping = aes(fill = ..count../sum(..count..)), col = "black") +
    scale_fill_gradientn(colours = colorRampPalette(blues9)(256), guide = FALSE) +
    ## dashed line indicating median mutrate
    geom_vline(xintercept = mean_mutrate, linetype = "dotted", col = "black", size = 1) +

    annotate(geom="text", x=Inf, y=Inf, vjust = 1, hjust = 1.1,
             #x=mean_mutrate, y=Inf, vjust = 1, hjust = -.1,
             label=paste0("Mean 11-mer mutation rate\n", round(mean_mutrate,2), " SNV/patient/Mb"))+
    ## add cases
    geom_segment(inherit.aes = FALSE,
                 data = dt[kmer %in% cases][, kmer := factor(kmer, levels = cases)],
                 mapping = aes(x = mrate, xend = mrate, y = c(9*10^4, 15*10^4), yend = c(3*10^4, 8*10^4)),
                 arrow = arrow(type="closed",angle = 10)) +
    geom_point(inherit.aes = FALSE,
               data = dt[kmer %in% cases][, kmer := factor(kmer, levels = cases)], 
               mapping = aes(x = mrate, y = c(10*10^4, 16*10^4), shape = kmer), 
               col = "black", fill = "white", alpha = 1, size = 3) +
    geom_text(inherit.aes = FALSE, data = dt[kmer %in% cases],
              mapping = aes(x = mrate, y = c(10*10^4, 16*10^4), label = kmer), 
              angle = 0, vjust = 0, hjust = 0, nudge_y = 10^4) +
    scale_shape_manual(name = "", values = setNames(c(0, 1), cases), 
                       guide = "none") +
    ## modify x axis (vertical when rotated)
    scale_x_continuous(name = "k-mer mutation rate\n(SNV/patient/Mb)",
                       trans = "log10",
                       #breaks = breaks_rate_patient*10^6,
                       breaks = c(0.12, 1, 10, 100, 1000),
                       #labels = function(x) ifelse(x<1, round(x, 1), comma(x, 1)),
                       labels = c("0", "1", "10", "100", "1,000"),
                       expand = c(0, 0)) +
    ## modify y axis (horizontal when rotated)
    scale_y_continuous(name = "count (x1,000)",
                       expand = c(0, 0),
                       breaks = y_breaks, #c(100, 200, 300)*1000,
                       labels = function(x) ifelse(x<1000, round(x/1000,3), comma(x/1000, accuracy = 1))) +
    ## modify theme
    coord_cartesian(xlim=range(breaks_rate_patient)*10^6, ylim=ylimits) +
    theme_pubr()
}
p_hist_rate <- hist_rate(k11c_dt, y_breaks = c(100, 200, 300)*1000, ylimits = c(0, 3*10^5))
#p_hist_rate
#ggsave(p_hist_rate, device = "pdf",
#       filename = "plots/kmer_general_mutrates/histogram_mutrate_11mer_all_v3.pdf", 
#       width = 7.3, height = 3.7, units = "in", useDingbats = FALSE)


# hist_genome_distribution <- function(x){p <- plot(x); return(p)}


k11c_aggr <- aggregate(tot_hg19 ~ mrate_log10, data = k11c_dt[, .(tot_hg19, mrate_log10=round(log10(mrate_nonzero),1))], FUN = sum)
setDT(k11c_aggr)
k11c_aggr[, `:=`(mrate = 10^mrate_log10)] #, fraction_of_genome = tot_hg19/sum(tot_hg19))]

mean_mutrate_weighted <- k11c_dt[, sum(tot_mut)/sum(tot_hg19)/2583*10^6]
  
p_aggre_genome_mrate_distrib <- k11c_aggr %>%
  ggplot(aes(mrate+0.02, tot_hg19, fill=tot_hg19)) +
  geom_col(col="black", width = 0.09) +
  geom_vline(xintercept = mean_mutrate_weighted, linetype = "dashed", col = "black", size = 1) +
  annotate(geom="text", x=Inf, y=Inf, vjust = 1, hjust = 1.1,
           label=paste0("Mean genome mutation rate\n", round(mean_mutrate_weighted,2), " SNV/patient/Mb"))+
  scale_x_continuous(name = "k-mer mutation rate\n(SNV/patient/Mb)",
                     trans = "log10",
                     breaks = c(0.12, 1, 10, 100, 1000),
                     labels = c("0", "1", "10", "100", "1,000"),
                     expand = c(0,0)) +
  
  scale_y_continuous(name="Genomic span (Mb)", labels = function(y) round(y/10^6, 1), expand=c(0,0), 
                     sec.axis = sec_axis(name = "Genomic span (%)", trans = function(y) y/sum(k11c_aggr$tot_hg19),
                                         breaks = (1:4*10^8)/sum(k11c_aggr$tot_hg19), 
                                         labels = function(y) format(round(y*100, 1), nsmall=1)))+
  scale_fill_gradientn(colours = colorRampPalette(blues9)(256), guide = FALSE) +
  coord_cartesian(xlim=range(breaks_rate_patient)*10^6, ylim=c(0,400*10^6)) +
  theme_pubr() +
  theme(axis.title.y.right = element_text(angle = 90))

p_mratehist_genomespanmrate <- p_hist_rate+theme(axis.title.x = element_blank()) + 
  p_aggre_genome_mrate_distrib + plot_layout(ncol=1)
ggsave(p_mratehist_genomespanmrate, device = "pdf",
       filename = "plots/kmer_general_mutrates/histx2_mrate_count_genomesum_v3.pdf", 
       width = 7, height = 7, units = "in", useDingbats = FALSE)


## recurrence subsets ====
### smoothscatter
# p_smoothscat_all <- smoothscat_test2(k11c_dt, legend_shape = FALSE) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# p_smoothscat_rec2plus <- smoothscat_test2(k11c_dt[rec2_4>0], legend_shape = FALSE) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# p_smoothscat_rec5plus <- smoothscat_test2(k11c_dt[rec5plus>0], legend_shape = FALSE)
# 
# ### histogram
# p_hist_rate_all <- hist_rate_test(k11c_dt, y_breaks = c(100, 200, 300)*1000, ylimits = c(0, 3*10^5)) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# p_hist_rate_rec2plus <- hist_rate_test(k11c_dt[rec2_4>0], y_breaks = c(10, 20, 30)*1000, ylimits = c(0, 3.5*10^4)) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# p_hist_rate_rec5plus <- hist_rate_test(k11c_dt[rec5plus>0], y_breaks = c(0.1, 0.2, 0.3)*1000, ylimits = c(0, 3*10^2))
# 
# p_smoothscat_histrate_subset <-
#   p_smoothscat_all + p_smoothscat_rec2plus + p_smoothscat_rec5plus + 
#   p_hist_rate_all + p_hist_rate_rec2plus + p_hist_rate_rec5plus +
#   patchwork::plot_layout(ncol = 2, nrow = 3, byrow = FALSE)
# 
# ggsave(plot = p_smoothscat_histrate_subset, device = "pdf",
#        filename = "plots/kmer_general_mutrates/smoothscat_histrate_11mer_refcount_subsets_pcawgmuts_v4.pdf",
#        width = 11, height = 12, units = "in", useDingbats = FALSE)
