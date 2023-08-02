# binomial test for all 11-mers with Bonferroni correction

library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)
library(patchwork)

setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

# READ DATA ====================
kmer11_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")
kmer11_dt[, mrate := n_snv/tot_hg19/2583*10^6]

# BINOMIAL TEST ===================
# make a one-sided binomial test
# we are not interested in the mutation rates less than expected
kmer11_dt[tot_hg19!=0 & n_snv<=tot_hg19, pval:=binom.test(n_snv, tot_hg19, 5.96/10^6*2583, alternative = "greater")$p.value, by=kmer_11]
kmer11_dt[n_snv>tot_hg19, pval:=binom.test(tot_hg19, tot_hg19, 5.96/10^6*2583, alternative = "greater")$p.value, by=kmer_11]

# BONFERRONI CORRECTION ===================
# kmer11_dt[, p.bonf:=p.adjust(p = pval, method = "bonferroni")]

# TRANSFORM SCALE OF P-VALUES =================
# kmer11_dt[, ppval:=-log10(p.bonf)]
kmer11_dt[, ppval:=-log10(pval)]

# SIMULATE MUTATION RATES ================================
p_interval <- c(2,5,9)

sim <- purrr::map_df(.x = p_interval, ~{
  famsize <- unique(round(10^seq(0, log10(max(kmer11_dt$tot_hg19)), length.out = 100),0))
  mrate <- sapply(famsize, function(size){
    n_snv = qbinom(p = 1-(10^(-.x)),
                   size = size,
                   prob = 5.96/10^6*2583,
                   lower.tail = TRUE)
    rate = n_snv/size/2583*10^6
    return(rate)
  })
  
  cbind(mrate, famsize) %>% 
    as.data.table() %>%
    dplyr::mutate(p = .x)
}) %>% dplyr::mutate(p = factor(p))



mrate_increase <- c(2, 5)
dt_labels <- data.table(x=2.5,
                        y=5.96*mrate_increase,
                        factor=paste0(mrate_increase,"x"))


# coarse-grained genomic bins
span_breaks <- 10^(0:7)[c(1,3,5,8)]

# fine-grained genomic bins
#famsize = unique(round(10^seq(0, log10(max(kmer11_dt$tot_hg19)), length.out = 50),0))
famsize = c(unique(round(10^seq(0, log10(5*10^4), length.out = 30),0)), 5*10^6)
# famsize = c(1, 2, 3, 6, 10, 20, 30, 60, 100, 170, 300, 530,
#             10^3, 1.7*10^3, 3*10^3, 5*10^3, 9*10^3, 30*10^3, 50*10^3, 50*10^5)
kmer11_dt[, span_bin := cut(x = tot_hg19,
                            breaks = famsize,
                            right = F,
                            labels = paste0(famsize[1:(length(famsize)-1)],"-",famsize[2:length(famsize)]))]


rate_signif_span <- function(factor=1, ppvalue=2){
  d <- copy(kmer11_dt[tot_hg19>0, .(n=sum(mrate>=5.96*factor, na.rm=T),
                                    x=sum(mrate>=5.96*factor & ppval>=ppvalue, na.rm=T),
                                    # pval<=10^(-ppvalue), na.rm=T),
                                    mrate=mean(mrate, na.rm = T),
                                    span=sum(as.double(tot_hg19)), span_min=min(tot_hg19), span_max=max(tot_hg19),
                                    n_kmer=.N,
                                    factor_increase=factor, pval_threshold=ppvalue), by=span_bin])
  d[, frac:=x/n]
  d[!is.finite(frac), frac:=0]
  return(d[order(as.numeric(span_min))])
}

rate_signif_span(factor = 5, ppvalue = 2)

d <- rbindlist(lapply(p_interval,function(p)
  rbindlist(lapply(mrate_increase, function(f)
    rate_signif_span(factor = f, ppvalue = p)))))

d[factor_increase%in%c(2) &
    pval_threshold%in%c(2), .(span_min, span_max)]

# fine-grained genomic bins ====
p_margin_v2 <- d %>%
  ggplot(aes(x = span_max, #factor(span_bin, levels=paste0(famsize[1:(length(famsize)-1)],"-",famsize[2:length(famsize)])),
             y = frac, col = factor(pval_threshold), group = factor(pval_threshold))) +
  geom_line(size = 1) +
  scale_color_brewer(name=expression(-log[10]~`p-value`), palette="Reds",
                     guide = "none") +
  scale_x_log10(name="Genomic span", breaks = 10^seq(0, 6, by=3),
                expand=c(0,0), labels = c("1 bp", "1 kb", "1 Mb")) +
  scale_y_continuous(name="Fraction significant (%)", breaks=c(0,.5,1),
                     labels=function(x) round(x*100),
                     sec.axis=sec_axis(name="Factor increase", trans="identity")) +
  annotation_logticks(sides="b") +
  facet_grid(factor(paste0("≥",factor_increase,"x"),levels=paste0("≥",c(10,5,2,1),"x"))~.) +
  theme_pubr() +
  theme(strip.placement = "outside", strip.text.y = element_text(angle=0),
        panel.spacing = unit(5, "mm"),
        axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank())
p_margin_v2


# DENSITY CLOUD WITH SIGNIFICANCE LINES ========================

kmer11_dt[, mrate_pseudo:=ifelse(mrate==0, 0.1, mrate)]
p_dens <- kmer11_dt %>%
  ggplot(aes(x=tot_hg19, y=mrate_pseudo)) +
  stat_density2d(geom = "tile", aes(fill=..density..^(1/6), alpha=..density..^(1/6)),
                 contour = FALSE, show.legend = FALSE) + 
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
  scale_alpha_continuous(range=c(0,1), guide="none") +
  geom_hline(yintercept = 5.96) +
  geom_hline(yintercept = 5.96*mrate_increase, col="grey", linetype="dashed") +
  geom_text(inherit.aes = FALSE, data = dt_labels, aes(x=x, y=y, label=factor), vjust=0, hjust=1) +
  annotation_logticks(sides = "lb") +
  geom_line(data=sim, aes(famsize, mrate, color=factor(p)), size=1) +
  scale_x_log10(name="Genomic span", expand=c(0,0),
                breaks = 10^seq(0, 6, by=3), minor_breaks=span_breaks,
                labels = c("1 bp", "1 kb", "1 Mb"), limits=c(1,4.7*10^6)) +
  scale_y_log10(name="Mutation rate (SNV/Mb/patient)", expand=c(0,0),
                breaks=10^(-1:4), limits=c(0.1, 10^3),
                labels=c(0, 1, 10, expression(10^2), expression(10^3), expression(10^4))) + #trans_format("log10", math_format(10^.x))
  scale_color_brewer(name="p-value", palette="Reds", #expression(-log[10]~`p-value`)
                     labels=math_format(10^-.x),
                     guide = guide_legend(reverse = T, title.position = "top", label.position = "left")) +
  theme_pubr() +
  theme(legend.position = c(.85, .25), legend.background = element_blank(),
        legend.key = element_blank(), legend.box.background = element_blank(),
        legend.direction = "vertical", legend.margin = margin(6, 6, 6, 6))

p_dens


# merge plots with patchwork ====
p_merged <- p_dens + p_margin_v2 +
  plot_layout(ncol=1, heights = c(2,1)) +
  plot_annotation(tag_levels = 'a')

# save with ggsave ====
# ggsave(plot = p_merged, device = "pdf", width = 17, height = 16, units = "cm", useDingbats=FALSE,
#        filename = "plots/binomtest/binomtest_11mers_fulldataset_v2.pdf")

# save with cairo_pdf ====
grDevices::cairo_pdf(filename = "plots/binomtest/binomtest_11mers_fulldataset_cairo_v4.pdf",
                     width = 17/2.54, height = 16/2.54)
p_merged
dev.off()
