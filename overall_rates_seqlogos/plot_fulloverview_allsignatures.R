
setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

library(data.table)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggpubr)
library(ggseqlogo)
library(patchwork)

# load data ====
dt <- fread("results/fulloverview_allsignatures.tsv")

feature_ord <- c("","Active_Promoter", "Poised_Promoter", "Weak_Promoter", "Strong_Enhancer",
  "Weak_Enhancer", "Insulator", "Txn_Transition", "Txn_Elongation", "Weak_Txn",
  "Repressed", "Repetitive/CNV","Heterochrom/lo")

rec_ord <- c("(1|2|3|4|5|6|7+)", "(2|3|4|5|6|7+)", "(5|6|7+)")

rec_feature_ord <- expand.grid(rec_ord, feature_ord, stringsAsFactors = FALSE) %>% 
  mutate(v3=paste(Var1, Var2),
         Var1=factor(Var1, rec_ord),
         Var2=factor(Var2, feature_ord)) %>% 
  arrange(Var1, Var2) %>% .$v3

load("../results/motif_analysis/sbs_k11_pfm_weighted_by_signatures.Rdata")
b1 <- Reduce("+", lapply(1:22, function(chr) oligonucleotideFrequency(Hsapiens[[chr]], width = 1)))
pfm_k1_hg19 <- matrix(data = ( b1/sum(b1) ), nrow = 4, ncol = 11, dimnames = list(names(b1)))


# add dummy data
dt_dummy <- dt[step==4, .(step=4, rec="(5|6|7+)", wmrate=0,
                          wpfm=paste0(c("A","C","G","T"),paste0(rep("\t1", 11), collapse=""), collapse = "\n"),
                          feature_collapse=feature_ord[-1][ !(feature_ord[-1] %in% feature_collapse) ]), by=signature]
dt_full <- rbind(dt, dt_dummy[!is.na(feature_collapse)], fill=TRUE)

# define plotting functions ====
plot_mrate <- function(d) {
  p <- d %>% ggplot(aes(x = as.character(step),
                        y = factor(paste(rec, feature_collapse), rec_feature_ord),
                        fill = log10(wmrate), 
                        label = ifelse(wmrate==0,"",
                                       format(wmrate, big.mark = ",", nsmall = 1, digits = 2)))) +
    geom_tile(col="black") +
    geom_text(aes(fontface=ifelse(wmrate==max(wmrate), "bold","plain"))) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0), position="right",
                     labels=function(x) {
                       ifelse(d$step<=2,
                              paste("Signature",str_sub(d$signature, 21,-3)),
                              ifelse(d$step==3, x,
                              gsub(pattern = "[$\\(].*[\\)]", replacement = "", x)))}) +
    scale_fill_gradient2(low = "white", mid = "white", high = "brown", 
                         midpoint = log10(5.96), guide = "none", limits = log10(range(dt$wmrate))) +
    labs(x="", y="") +
    coord_flip(clip = "off") +
    theme_pubr()+
    theme(axis.title = element_blank(), axis.text.x = element_text(angle=90, hjust=0),
          axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin=unit(c(0,-.5,-.5,0), "cm"))
  return(p)
}
plot_logo <- function(d) {
  pfm_list <- lapply(1:d[, .N], 
         function(i) {
           pfm <- d[i, wpfm] %>%
             fread(colClasses = c("character",rep("integer",11))) %>% 
             as.matrix() %>% .[1:4, 2:12] %>%
             apply(MARGIN = 2, as.numeric)
           
           pfm_obs <- pfm/sum(pfm[,1])
           pfm_exp <- pfm_k1_hg19
           D_KL <- colSums(pfm_obs*log2(pfm_obs/pfm_exp), na.rm = TRUE)
           compute_height <- function(obs, ic) sapply(1:11, function(i) obs[,i]*ic[i])
           heights <- compute_height(pfm_obs, D_KL)

           rownames(heights) <- c("A","C","G","T")
           return(heights)
         })
  
  if(length(pfm_list)>1){
    names(pfm_list) <- d[, paste(rec, feature_collapse)]
    pfm_list <- pfm_list[as.character(sort(factor(names(pfm_list), rec_feature_ord)))]
  }
  
  p <- ggseqlogo(pfm_list, nrow = 1, method="custom", seq_type="dna") +
    scale_x_continuous(name = "", breaks=1:11, labels = rep("",11)) +
    scale_y_continuous(breaks = 0:2, labels = function(x) format(x,nsmall = 0)) +
    coord_cartesian(ylim=c(0,2.2)) +
    theme(strip.text = element_blank(), 
          plot.margin = unit(c(-.5,-.5,-.5,-.5), "cm"))
    #scale_x_continuous(breaks=1:11, labels=1:11-6)
  return(p)
}
plot_logo_v2 <- function(d, asseenbysig=FALSE) {
  # create list of pfms (position frequency matrices)
  pfm_list <- lapply(1:d[, .N],
                     function(i) {
                       pfm <- d[i, wpfm] %>%
                         fread(colClasses = c("character",rep("integer",11))) %>% 
                         as.matrix() %>% .[1:4, 2:12] %>%
                         apply(MARGIN = 2, as.numeric)
                       
                       pfm_obs <- pfm/sum(pfm[,1])
                       if(asseenbysig) {
                         pfm_exp <- sbs_k11[[d[i,signature]]]
                         } else pfm_exp <- pfm_k1_hg19
                       
                       D_KL <- colSums(pfm_obs*log2(pfm_obs/pfm_exp), na.rm = TRUE)
                       compute_height <- function(obs, ic) sapply(1:11, function(i) obs[,i]*ic[i])
                       heights <- compute_height(pfm_obs, D_KL)
                       
                       rownames(heights) <- c("A","C","G","T")
                       return(heights)
                     })
  
  # name each element of the list
  if(length(pfm_list)>1) names(pfm_list) <- d[, paste(rec, feature_collapse)]
  
  # # <test (add dummy data to missing elements)>
  # missing_features <- feature_ord[-1][ !(feature_ord[-1] %in% d$feature_collapse) ]
  # #return(missing_features)
  # if(length(missing_features)>0){
  #   dummy_matrix <- matrix(data=rep(0.001,44),nrow=4, ncol=11, dimnames = list(c("A","C","G","T")))
  #   for(i in 1:length(missing_features)) {
  #     name_missing_features <- paste(unique(d$rec), missing_features[i])
  #     pfm_list[[name_missing_features]] <- dummy_matrix
  #   }
  # }
  # return(pfm_list)
  # # </test>
  
  if(length(pfm_list)>1) pfm_list <- pfm_list[as.character(sort(factor(names(pfm_list), rec_feature_ord)))]
  
  p <- ggseqlogo(pfm_list, nrow = 1, method="custom", seq_type="dna") +
    scale_x_continuous(name = "", breaks=1:11, labels = rep("",11)) +
    scale_y_continuous(name="bits",breaks = 0:2, labels = function(x) format(x,nsmall = 0)) +
    coord_cartesian(ylim=c(0,2.2)) +
    theme(strip.text = element_blank(), axis.title.y = element_text(size = 6),
          plot.margin = unit(c(-.5,-.5,-.5,-.5), "cm"))
  return(p)
}
plot_logo_v2(dt[step==4 & signature == "BI_COMPOSITE_SNV_SBS16_P"], asseenbysig = T)

plot_span <- function(d) {
  span_kb <- function(x) paste(format(x/10^3, small.mark = ".", nsmall = 1, digits = 2),"Kb")
  span_mb <- function(x) paste(format(x/10^6, big.mark = ",", nsmall = 0, digits = 0),"Mb")
  
  p <- d %>%
    ggplot(aes(x = as.character(step),
               y = factor(paste(rec, feature_collapse), rec_feature_ord), fill = log10(tot_hg19),
               label = ifelse(is.na(tot_hg19), "",
                              paste0(ifelse(tot_hg19/10^6>=1,
                                     span_mb(tot_hg19),span_kb(tot_hg19)),
                              "\n(", round(tot_hg19/sum(b1)*100,3),"%)")))) +
    geom_tile(col="black") +
    geom_text() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "white", 
                         midpoint = log10(sum(b1)), guide = "none", limits = log10(range(dt$tot_hg19))) +
    labs(x="", y="") + 
    coord_flip(clip = "off") +
    theme_pubr() +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin=unit(c(-.5,-.5,0,-.5), "cm"))
  return(p)
}

plot_full <- function(signature_name){
  p_list <- lapply(1:4, function(i){
    p_mrate <- dt[step==i & signature == signature_name] %>% plot_mrate()
    p_logo <- dt[step==i & signature == signature_name] %>% plot_logo()
    p_span <- dt[step==i & signature == signature_name] %>% plot_span()
    return(list("mrate"=p_mrate, "logo"=p_logo, "span"=p_span))
    })
  #return(p_list)
  #p_arranged <- Reduce("+", unlist(p_list, recursive = F)) + plot_layout(nrow=1)
  # p_arranged <-
  #   #1
  #   (p_list[[1]][["mrate"]] + p_list[[1]][["logo"]] + p_list[[1]][["span"]] +
  #   #2
  #   p_list[[2]][["mrate"]] + p_list[[2]][["logo"]] + p_list[[2]][["span"]] +
  #     plot_layout(ncol=2, nrow=3, byrow=FALSE)) +
  #   #3
  #   (p_list[[3]][["mrate"]] + p_list[[3]][["logo"]] + p_list[[3]][["span"]] +
  #   plot_layout(ncol=))
    
    #4
    #(p_list[[4]][["mrate"]] + p_list[[4]][["logo"]] + p_list[[4]][["span"]] + plot_layout(nrow=1, widths = c(1,1,1))) + 
    
    #layout
    #plot_layout(nrow=1)
  
  p_unlist <- unlist(p_list, recursive = F)
  design_var <- "
  ADGGJJJJJJJJ
  BEHHKKKKKKKK
  CFIILLLLLLLL"
  p_arranged <- Reduce("+", p_unlist) + plot_layout(design = design_var)
  return(p_arranged)
}
plot_full_v2 <- function(signature_name){
  p_list <- lapply(1:4, function(i){
    p_mrate <- dt_full[step==i & signature == signature_name] %>% plot_mrate()
    p_logo1 <- dt_full[step==i & signature == signature_name] %>% plot_logo_v2(asseenbysig = FALSE)
    p_logo2 <- dt_full[step==i & signature == signature_name] %>% plot_logo_v2(asseenbysig = TRUE)
    p_span <- dt_full[step==i & signature == signature_name] %>% plot_span()
    return(list("mrate"=p_mrate, "logo_bits"=p_logo1, "logo_sigbias"=p_logo2, "span"=p_span))
  })
  
  p_unlist <- unlist(p_list, recursive = F)
  # design_var <- "
  # ADGGJJJJJJJJ
  # BEHHKKKKKKKK
  # CFIILLLLLLLL"
  design_var2 <- "
  AEIIMMMMMMMM
  BFJJNNNNNNNN
  CGKKOOOOOOOO
  DHLLPPPPPPPP"
  p_arranged <- Reduce("+", p_unlist) + plot_layout(design = design_var2)
  return(p_arranged)
}
plot_full_v3 <- function(signature_name, stackingfriendly=TRUE){
  
  theme_hstrip <- theme()
  theme_0strip <- theme()
  
  if(stackingfriendly) {
    theme_hstrip <- theme(axis.text.x = element_text(angle=0))
    theme_0strip <- theme(axis.text.x = element_blank())
    }
  
  p_list <- lapply(1:4, function(i){
    if(i <= 2){
      p_mrate <- dt_full[step==i & signature == signature_name] %>% plot_mrate() + theme_hstrip
      p_logo1 <- dt_full[step==i & signature == signature_name] %>% plot_logo_v2(asseenbysig = FALSE)
      p_logo2 <- dt_full[step==i & signature == signature_name] %>% plot_logo_v2(asseenbysig = TRUE)
      p_span <- dt_full[step==i & signature == signature_name] %>% plot_span()
    } else {
      p_mrate <- dt_full[step==i & signature == signature_name] %>% plot_mrate() + theme_0strip
      p_logo1 <- dt_full[step==i & signature == signature_name] %>% plot_logo_v2(asseenbysig = FALSE)
      p_logo2 <- dt_full[step==i & signature == signature_name] %>% plot_logo_v2(asseenbysig = TRUE)
      p_span <- dt_full[step==i & signature == signature_name] %>% plot_span()
      }
    return(list("mrate"=p_mrate, "logo_bits"=p_logo1, "logo_sigbias"=p_logo2, "span"=p_span))
  })
  
  p_unlist <- unlist(p_list, recursive = F)
  design_var2 <- "
  AEIIMMMMMMMM
  BFJJNNNNNNNN
  CGKKOOOOOOOO
  DHLLPPPPPPPP"
  p_arranged <- Reduce("+", p_unlist) + plot_layout(design = design_var2)
  return(p_arranged)
}
plot_full_v2("BI_COMPOSITE_SNV_SBS3_P")
plot_full("BI_COMPOSITE_SNV_SBS16_P")

# plot and save ====
pdf(width=16, height=4, file="plots/factorization_steps/supplemental/fulloverview_allsignatures.pdf")
lapply(unique(dt$signature), plot_full)
dev.off()

# newest edition right here (v2)
pdf(width=16, height=4, file="plots/factorization_steps/supplemental/fulloverview_allsignatures_v2.pdf")
lapply(unique(dt$signature), plot_full_v2)
dev.off()

# remove text (v3)
p_list_stackfriendly <- lapply(unique(dt$signature), plot_full_v3)
pdf(width=16, height=2.667, #width=6.299, height=1.059,
    file="plots/factorization_steps/supplemental/fulloverview_allsignatures_stackfriendly.pdf")
p_list_stackfriendly
dev.off()



# testing v2 ====
par(mfrow=c(1,3))

p_list[[2]]

gridExtra::grid.arrange(p_list[[1]], p_list[[2]], p_list[[3]], ncol=1)

Reduce("+", lapply(1:2, function(i) as_ggplot(p_list[[i]])))

as_ggplot(p_list[[1]])

gridExtra::gtable_rbind(p_list[[1]], p_list[[2]])




Reduce("+", p_list[1:2])# + plot_layout(ncol=1, nrow=3)
ggarrange(plotlist = p_list[1:3])



# testing ====
Reduce("+", lapply(1:3,
  function(i){
    plot_logo_v2(dt[step==i &
                      signature == "BI_COMPOSITE_SNV_SBS17b_P"],
                 asseenbysig = FALSE)
    })) + plot_layout(ncol=1)


dt[step==1 & signature == "BI_COMPOSITE_SNV_SBS17b_P", wpfm]%>%
  fread(colClasses = c("character",rep("integer",11))) %>% 
  as.matrix() %>% .[1:4, 2:12] %>%
  apply(MARGIN = 2, as.numeric)

consensusString_fromSignature <- function(signature_name, asseenbysig=TRUE){
  pfm_list_out <- lapply(1:4, function(i)
    {
    d <- dt[step==i & signature == signature_name]
    
    pfm_list <- lapply(1:d[, .N],
                       function(i)
                         {
                         pfm <- d[i, wpfm] %>%
                           fread(colClasses = c("character",rep("integer",11))) %>% 
                           as.matrix() %>% .[1:4, 2:12] %>%
                           apply(MARGIN = 2, as.numeric)
                         pfm_obs <- pfm/sum(pfm[,1])
                         rownames(pfm_obs) <- c("A","C","G","T"); return(consensusString(pfm_obs))
                         if(asseenbysig)
                           {
                           pfm_exp <- sbs_k11[[d[i,signature]]]
                           } else pfm_exp <- pfm_k1_hg19
                         
                         D_KL <- colSums(pfm_obs*log2(pfm_obs/pfm_exp), na.rm = TRUE)
                         compute_height <- function(obs, ic) sapply(1:11, function(i) obs[,i]*ic[i])
                         heights <- compute_height(pfm_obs, D_KL)
                         
                         rownames(heights) <- c("A","C","G","T")
                         return(heights)
                         })
    
    if(length(pfm_list)>1){
      names(pfm_list) <- d[, paste(rec, feature_collapse)]
      pfm_list <- pfm_list[as.character(sort(factor(names(pfm_list), rec_feature_ord)))]
      }
    })
  return(pfm_list_out)
  }

consensusString_fromSignature("BI_COMPOSITE_SNV_SBS17b_P")


#####

if(length(pfm_list)>1){
  names(pfm_list) <- d[, paste(rec, feature_collapse)]
  pfm_list <- pfm_list[as.character(sort(factor(names(pfm_list), rec_feature_ord)))]
}
#####

pfm <- dt[step==3 &signature == "BI_COMPOSITE_SNV_SBS17b_P"][2, wpfm] %>% 
  fread(colClasses = c("character",rep("integer",11))) %>%
  as.matrix() %>% .[1:4, 2:12] %>%
  apply(MARGIN = 2, as.numeric)

pfm_obs <- pfm/sum(pfm[,1])
rownames(pfm_obs) <- c("A","C","G","T")
consensusString(pfm_obs, threshold=0.5)

pfm_exp <- sbs_k11[["BI_COMPOSITE_SNV_SBS17b_P"]]

D_KL <- colSums(pfm_obs*log2(pfm_obs/pfm_exp), na.rm = TRUE)
compute_height <- function(obs, ic) sapply(1:11, function(i) obs[,i]*ic[i])
heights <- compute_height(pfm_obs, D_KL)
ggseqlogo(heights, method="custom") + coord_cartesian(ylim=c(0,2))
ggseqlogo(pfm_obs)

consensusString(round(pfm_exp,4))
