

library(data.table)

setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")
sigs <- fread("data/reference_signatures/SignatureAnalyzer_COMPOSITE_SBS_W96.signature.031918.txt")
sbs_etio <- fread("data/reference_signatures/SignatureAnalyzer_SBS_etiology.tsv")

cossim <- function(a, b)
  sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2)))

sbs_names <- grep("BI_COMPOSITE", colnames(sigs), value = TRUE)
simmatrix <- sapply(sbs_names, function(i)
  sapply(sbs_names, function(ii)
    cossim(sigs[, i, with=F], sigs[, ii, with=F])))
simmatrix_full <- simmatrix
simmatrix[upper.tri(simmatrix)] <- NA


hc <- hclust(as.dist(simmatrix))


hc$labels[hc$order]

sim_dt <- simmatrix %>%
  #as.dist(diag=TRUE) %>%
  #as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "signature_i") %>%
  gather(key = "signature_ii", value = "cossim", c(2:61)) %>%
  as.data.table()

sim2_dt <- merge(
  x = merge(sim_dt, sbs_etio[, .(BI_COMPOSITE, Etiology_i=Etiology)], by.x="signature_i", by.y="BI_COMPOSITE"),
  y = sbs_etio[, .(BI_COMPOSITE, Etiology_ii=Etiology)],
  by.x="signature_ii", by.y="BI_COMPOSITE")

#sim2_dt[, signature_i := factor(signature_i, levels = hc$labels[hc$order])]
#sim2_dt[, signature_ii := factor(signature_ii, levels = rev(hc$labels[hc$order]))]
sim2_dt[, `:=`(signature_i = factor(signature_i, levels = rev(sbs_names)),
               signature_ii = factor(signature_ii, levels = sbs_names))]

sigs_of_interest <- paste("(S7a","67","75","S7b","62","10a",
            "S9_","72","17","19","68","28","30)", sep = "|")
#sigs_of_interest <- paste("(S62","S10a","S9","S72","S37","17",
#  "19","68","28", "30","18","S7a","67","75","7b)", sep = "|")

etio_lvls <- c("5mCdeamin","Artifact Mode")
sim2_dt[, sig_comparison:=(paste(pmin(as.character(signature_ii), as.character(signature_i)), 
                  pmax(as.character(signature_ii), as.character(signature_i)))),
        by = .(signature_ii)]
sim2_dt <- sim2_dt[order(signature_ii, Etiology_ii)]
sim2_dt[, x:=cossim]
sim2_dt[duplicated(sig_comparison), x:=0]

sim2_dt[, lab:=""]
sim2_dt[x>=0.5 & signature_i!=signature_ii,
        lab:=round(cossim*100,0)]
sim2_dt[, textcolor:=ifelse(x>=0.8, "white", "black")]

p_cossim <- sim2_dt %>% #[grepl(sigs_of_interest, signature_i)] %>% 
  ggplot(aes(signature_ii, signature_i, fill=x)) +
  geom_tile() +
  geom_text(aes(label=lab, col=textcolor), size=2)+
  scale_x_discrete(labels=function(x) str_sub(x, 21, -3)) +
  scale_y_discrete(labels=function(x) str_sub(x, 21, -3)) +
  scale_fill_gradient2(low="white", high="brown", midpoint = 0.4)+
  scale_color_identity() +
  #facet_grid(Etiology_ii~Etiology_i, switch = "both",
  #           scales = "free", space = "free") +
  coord_cartesian(clip="off") +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "mm"),
        axis.text.x = element_text(angle=90, vjust=0.5),
        strip.background = element_blank(),
        strip.text.x = element_text(angle=90, vjust=0.5, size=5),
        strip.text.y = element_text(angle=0, vjust=0.5, size=5),
        legend.position = "none", strip.placement = "outside")
p_cossim

ggsave(p_cossim, device = "pdf", width = 7, height = 7, units = "in",
       filename = "plots/signature_similarities/cossim_allsignatures.pdf")



# focused plot
sim2_dt <- simmatrix_full %>%
  #as.dist(diag=TRUE) %>%
  #as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "signature_i") %>%
  gather(key = "signature_ii", value = "cossim", c(2:61)) %>%
  as.data.table() %>%
  merge(sbs_etio[, .(BI_COMPOSITE, Etiology_i=Etiology)], by.x="signature_i", by.y="BI_COMPOSITE") %>%
  merge(sbs_etio[, .(BI_COMPOSITE, Etiology_ii=Etiology)],by.x="signature_ii", by.y="BI_COMPOSITE")

sim2_dt[, `:=`(signature_i = factor(signature_i, levels = rev(sbs_names)),
               signature_ii = factor(signature_ii, levels = sbs_names))]
sim2_dt <- sim2_dt[order(signature_ii, Etiology_ii)]
sim2_dt[, x:=cossim]
#sim2_dt[duplicated(sig_comparison), x:=0]

sim2_dt[, lab:=""]
sim2_dt[x>=0.5 & signature_i!=signature_ii,
        lab:=round(cossim*100,0)]
sim2_dt[, textcolor:=ifelse(x>=0.795, "white", "black")]

p_cossim_focus <- sim2_dt[grepl(sigs_of_interest, signature_i)] %>% 
  ggplot(aes(signature_ii, signature_i, fill=x)) +
  geom_tile() +
  geom_text(aes(label=lab, col=textcolor), size=2)+
  scale_x_discrete(labels=function(x) str_sub(x, 21, -3)) +
  scale_y_discrete(labels=function(x) str_sub(x, 21, -3)) +
  scale_fill_gradient2(low="white", high="brown", midpoint = 0.4)+
  scale_color_identity() +
  #facet_grid(Etiology_ii~Etiology_i, switch = "both",
  #           scales = "free", space = "free") +
  coord_cartesian(clip="off") +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "mm"),
        axis.text.x = element_text(angle=90, vjust=0.5),
        strip.background = element_blank(),
        strip.text.x = element_text(angle=90, vjust=0.5, size=5),
        strip.text.y = element_text(angle=0, vjust=0.5, size=5),
        legend.position = "none", strip.placement = "outside")

p_cossim_focus

ggsave(p_cossim_focus, device = "pdf", width = 7, height = 2.38, units = "in",
       filename = "plots/signature_similarities/cossim_focussignatures.pdf")



sim2_dt[grepl(sigs_of_interest, signature_ii) &
          grepl(sigs_of_interest, signature_i) &
         signature_i!=signature_ii & #(Etiology_i!=Etiology_ii) &
          cossim>=.5][order(Etiology_i, lab)][,
      .(str_sub(signature_ii,21,-3), str_sub(signature_i, 21, -3),
        Etiology_ii, Etiology_i, lab,
        isEtioSame=Etiology_i==Etiology_ii)][order(-rank(as.integer(lab)))]

sim2_dt[grepl("APOBEC", Etiology_i) & grepl("APOBEC", Etiology_ii) &
          signature_i!=signature_ii]

sim2_dt[grepl(sigs_of_interest, signature_ii) &
          grepl(sigs_of_interest, signature_i) &
          signature_i!=signature_ii][cossim>=0.75,
  lapply(.SD, head, n=1), by=signature_i][order(-rank(cossim))]


# do patients share the similar signatures?

e <- fread("data/PCAWG_signatures/SA_COMPOSITE_SNV.activity_relative.WHITELIST_SET.20200226.txt")
sum(e[BI_COMPOSITE_SNV_SBS36_P>=0.05, icgc_donor_id] %in%
e[BI_COMPOSITE_SNV_SBS10a_S>=0.05, icgc_donor_id])

sum(e[BI_COMPOSITE_SNV_SBS11_S>=0.05, icgc_donor_id] %in%
      e[BI_COMPOSITE_SNV_SBS30_P>=0.05, icgc_donor_id])

sum(e[BI_COMPOSITE_SNV_SBS78_S>=0.05, icgc_donor_id] %in%
      e[BI_COMPOSITE_SNV_SBS28_P>=0.05, icgc_donor_id])

sum(e[BI_COMPOSITE_SNV_SBS11_S>=0.05, icgc_donor_id] %in%
      e[BI_COMPOSITE_SNV_SBS7b_S>=0.05, icgc_donor_id])

s4_co <- apply(e[BI_COMPOSITE_SNV_SBS4_P>=0.05, c(8:67)], 2, function(x) sum(x>0.05))
s4_co_sim <- simmatrix_full["BI_COMPOSITE_SNV_SBS4_P", names(s4_co[s4_co>0])]
s4_co_sim[s4_co_sim>=0.5]

s61_co <- apply(e[BI_COMPOSITE_SNV_SBS61_S>=0.05, c(8:67)], 2, function(x) sum(x>0.05))
s61_co_sim <- simmatrix_full["BI_COMPOSITE_SNV_SBS61_S", names(s61_co[s61_co>0])]
s61_co_sim[s61_co_sim>=0.5]

