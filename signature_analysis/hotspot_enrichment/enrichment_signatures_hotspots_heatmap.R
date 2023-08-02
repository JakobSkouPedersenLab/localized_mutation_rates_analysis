# define environment ====

## set working directory
setwd("~/PCAWG/faststorage/Gustav/")
## load packages
library(data.table)
library(tidyverse)

# load data ====
## read summary data
recmean <- fread("mutprocesses_project/results/signaturemean_by_recurrence.tsv")
load("reference_work/signature_etiologies.RData")

# wrangle data ====
## modify etiology names
proposed_etiologies <- proposed_etiologies %>% 
  mutate(Etiology_short = ifelse(Etiology == "ACID", "AA",
                          ifelse(Etiology == "5mCdeamin", "Age",
                          ifelse(Etiology == "Artifact Mode", "Artifact",
                          ifelse(Etiology == "BER", "BERd",
                          ifelse(Etiology == "HR", "HRd",
                          ifelse(Etiology == "LYMPH", "Lymph",
                          ifelse(Etiology == "MMR", "MMRd",
                          ifelse(Etiology == "MMR + POLE", "MMRd+POLEd", 
                          ifelse(Etiology == "POLE", "POLEd",
                          ifelse(Etiology == "POLI", "POLId",
                          ifelse(Etiology == "TEMOZOLOMIDE", "TMZ",
                          
                          ifelse(Etiology == "TOBACCO", "Tobacco",
                          ifelse(Etiology == "PLATINUM", "PT",
                          ifelse(Etiology == "PROSTATE", "Prostate",
                          
                                 Etiology)))))))))))))))
sort_etio <- c("Age", "Artifact",  # 
               "BERd","HRd","MMRd","MMRd+POLEd","POLEd","POLId", # DNA repair deficiencies
               "Lymph", "Prostate", # associated tissue
               "Unknown",
               "APOBEC", "INDEL", "ROS", # intrinsic exposures
               "AA", "PT", "TMZ","Tobacco", "UV") # extrinsic exposures

## create reference variable for signature names
sbs_columns <- colnames(recmean)[grepl("BI_COMPOSITE_SNV_SBS", colnames(recmean))]
## gather data in long format rather than wide
recmean_long <- recmean %>% gather("signature", "prob", sbs_columns) %>% as.data.table()
## calculate fold change
recmean_long[, prob_rec1 := prob[rec==1], by = signature]
recmean_long[, prob_fc := log2((prob+0.01)/(prob_rec1+0.01))]
recmean_long <- merge(x = recmean_long, 
                      y = proposed_etiologies[, c("BI_COMPOSITE", "Etiology_short")],
                      by.x = "signature", by.y = "BI_COMPOSITE")

# visualize ====
p_heatmap <- recmean_long %>% 
  # set the order of signature names
  mutate(signature = factor(signature, levels = sbs_columns)) %>%
  # plot
  ggplot(aes(reorder(signature, prob_fc, FUN = function(x) -mean(x)), rec, fill = prob_fc)) +
  geom_tile() +
  
  scale_x_discrete(name = "Mutational signature (single base substitution)", labels = function(x) str_sub(x, 20, -3), ) +
  #ylab(label = "    Recurrence      Association") +
  ylab(label = expression(~~~~~~Recurrence~~~~~~italic(Association))) +
  scale_fill_gradient2(name = expression(italic(fold~change)==log[2]~frac(italic(obs), italic(exp))), low = muted("blue"), high = "brown",
                       guide = guide_colorbar(title.position = "left", 
                                              title.hjust = 0.5,
                                              title.theme = element_text(size = 10),
                                              label.position = "top", 
                                              barwidth = 30, 
                                              barheight = 0.5)) +
  #guides(fill = guide_legend(title.position = "top", label.position = "bottom")) +

  theme_pubr(legend = "bottom") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(hjust = 0),
        panel.spacing = unit(0.1, "lines"), 
        panel.border = element_rect(color = "black", fill = NA),
        strip.text = element_text(angle = 90, hjust = 0, face = "italic"), 
        strip.background = element_blank(),
        legend.box.spacing = unit(0, "mm")) +
  facet_grid(.~factor(Etiology_short, levels = sort_etio), scales = "free_x", space = "free_x")

# save plot ====
setwd("~/PCAWG/faststorage/Gustav/")
gdk_path_to_plot <- "mutprocesses_project/plots/signature_contribution_in_hotspots/"
ggsave(plot = p_heatmap, device = "pdf", width = 10, height = 3.6, units = "in",
       filename = paste0(gdk_path_to_plot, "enrichment_signatures_hotspots_rec1to7_heatmap.pdf"))
local_path_to_plot <- "~/OneDrive - Aarhus Universitet/Projects/Localized_mutational_processes/plots/signature_contribution_in_hotspots/"
ggsave(plot = p_heatmap, device = "pdf", width = 10, height = 3.6, units = "in",
       filename = paste0(local_path_to_plot,"enrichment_signatures_hotspots_rec1to7_heatmap.pdf"))

       