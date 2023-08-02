

options(scipen=999)

library(data.table)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("~/PCAWG/faststorage/Gustav/mutprocesses_project")

k11c_dt <- fread("results/11mer_count_hg19_PCAWG_v3.tsv")


rmsk <- fread("../reference_work/encodetracks/rmsk_RepeatMasker_hg19.bed", sep = "\t")
rmsk <- rmsk[genoName%in%paste0("chr", 1:22)]
setnames(rmsk, old=c("genoName","genoStart","genoEnd"),
         new=c("chr","start","end"))
setkey(rmsk, chr, start, end)
####
rmsk_length <- rmsk[, .(length=sum(end-start)), by=c("repClass")][order(-rank(length))]
rmsk_length <- rbind(rmsk_length, data.table(repClass=NA, length=sum(k11c_dt$tot_hg19)-sum(rmsk_length$length)))
rmsk_length[, frac_length:=length/sum(length)]
rmsk_length[, sum(frac_length), by=repClass=="<NA>"] # 49.92% repeat genome


enchmm <- fread("../reference_work/encodetracks/wgEncodeBroadHmmH1hescHMM.bed")
match_enchmm_text <- "(?<=[:digit:]{1,2}_)[:alpha:]*[:punct:]*[:alpha:]*"
enchmm[, name_collapse:=str_extract(name, match_enchmm_text)]
setnames(enchmm, old = c("chrom","chromStart", "chromEnd"), new = c("chr","start","end"))
setkeyv(enchmm, c("chr","start","end"))
enchmm[, length:=end-start]
#enchmm_stat <- enchmm[, .(total_length=sum(length)), by = name_collapse]
#enchmm_stat[, rel_length_feature:=total_length/sum(total_length)]

## snv data
### exclude columns to load faster
column_select <- c("Chromosome","End_position","Donor_ID","kmer_11")
snv <- fread("data/signature_snv_filtered_annotated4.tsv", select = column_select)
setnames(snv, c("Chromosome", "End_position"), c("chr", "end"))
snv[, start:=end] ; setkey(snv, chr, start, end)

# overlap data (40 sec) ====
t1 <- Sys.time()
fov <- foverlaps(x=snv, y=rmsk, type="any", by.x = c("chr","start","end"), mult="first")
t2 <- Sys.time() ; difftime(t2,t1)

# chromHMM
t1 <- Sys.time()
fov2 <- foverlaps(x=fov[,.(chr,start=i.start,end=i.end, repName,repClass,repFamily,
                           Donor_ID, kmer_11)], 
                  y=enchmm[,c(2:5,11:12)], type="any", by.x = c("chr","start","end"), mult="first")
t2 <- Sys.time() ; difftime(t2,t1)

fov2[grepl("Repetitive",name_collapse)&kmer_11=="GAAACTTCTTT", .N, 
     by=c("repClass","repFamily", "repName")][,frac:=N/sum(N)][order(-rank(frac))]
fov2[grepl("ALR/Alpha",repName)&kmer_11=="GAAACTTCTTT", .N, 
     by=name_collapse][,frac:=N/sum(N)][order(-rank(frac))]
fov2[kmer_11=="GAAACTTCTTT", .N, 
     by=c("repClass","repFamily", "repName")][,frac:=N/sum(N)][order(-rank(frac))]

fov[,.N]
fov_rep <- fov[, .(n_snv=.N, n_donor=uniqueN(Donor_ID)), by = repClass]
fov_rep2 <- merge(fov_rep, rmsk_length, by="repClass", all.x=TRUE)
fov_rep2[, mrate:=n_snv/n_donor/length*10^6]
fov_rep2[order(-rank(mrate))]

# balanced distribution
# of SNVs inside and outside repeats
fov[, is.rep:=!is.na(repClass)]
gw_stat <- rbind(fov[, .(n_snv=.N), by = is.rep][, .(is.rep, var="n_snv", value=n_snv, frac=n_snv/sum(n_snv))],
data.table(is.rep=c(TRUE,FALSE), var="length", value=c(sum(rmsk_length$length),
  sum(k11c_dt$tot_hg19)-sum(rmsk_length$length)))[, frac:=value/sum(value)][,])
gw_stat
tidyr::spread(gw_stat[, 1:3], key="var", value="value")[, n_snv/length/2583*10^6,by=is.rep]
#


fov_stat <- fov[, .(n_snv=.N, n_donor=uniqueN(Donor_ID)), by=c("kmer_11","repClass")][order(-rank(n_snv))]

# overall mutation rate per Repeat Class ====
fov_stat2 <- merge(fov_stat, rmsk_length, by="repClass", all.x=TRUE)
fov_stat2[is.na(repClass), length:=sum(k11c_dt$tot_hg19)-sum(rmsk_length$length)]
fov_stat2[, mrate_11mer_rmsk:=n_snv/length/2583*10^6]
fov_stat2[, mrate_rmsk:=sum(n_snv)/length/2583*10^6, by=repClass]
unique(fov_stat2[,.(repClass, mrate_rmsk, frac_length)])[order(-rank(mrate_rmsk))]

# relative occurrence of GAAACTTCTTT in Repeat Class ====
matches <- rbindlist(lapply(1:22, function(chr) {
  pattern_var <- "GAAACTTCTTT"
  pattern_rc_var <- as.character(reverseComplement(DNAStringSet(pattern_var)))
  
  match_var <- matchPattern(pattern = pattern_var, subject = Hsapiens[[chr]])
  match_rc_var <- matchPattern(pattern = pattern_rc_var, subject = Hsapiens[[chr]])
  
  dt <- as.data.table(match_var@ranges)[, `:=`(chr=paste0("chr",chr), pattern=pattern_var)]
  dt_rc <- as.data.table(match_rc_var@ranges)[, `:=`(chr=paste0("chr",chr), pattern=pattern_rc_var)]
  dt_c <- rbind(dt, dt_rc)
}))
setkeyv(matches, c("chr","start","end"))
fov_ref <- foverlaps(x = matches, y = rmsk[,.(chr, start, end, repClass, repFamily, repName)], mult = "first")
# GAAACTTCTTT mutation rate in Repeat Class
y <- merge(fov_stat[kmer_11=="GAAACTTCTTT"], fov_ref[, .(n_ref=.N),by=c("repClass")], by="repClass", all.x=TRUE)
y[, mrate_rmsk:=n_snv/n_ref/n_donor*10^6]
y[order(-rank(mrate_rmsk))]

###
# overall mutation rate per Repeat Family ====
fam_stat <- fov[, .(n_snv=.N, n_donor=uniqueN(Donor_ID), length=as.numeric(sum(end-start))), by=c("repClass","repFamily", "repName")]
fam_stat[is.na(repClass), length:=sum(k11c_dt$tot_hg19)-rmsk_fam_length[is.na(repClass), sum(length)]]
fam_stat[,mrate:=n_snv/length/n_donor*10^6]
fam_stat[n_donor>=500][order(-rank(mrate))][1:20]


fov_fam <- fov[, .(n_snv=.N, n_donor=uniqueN(Donor_ID)), by=c("kmer_11","repClass", "repFamily", "repName")][order(-rank(n_snv))]
fov_fam2 <- merge(fov_fam, fov_ref[, .(n_ref=.N),by=c("repClass","repFamily", "repName")],
      by=c("repClass","repFamily", "repName"), all.x=TRUE)
fov_fam2[, mrate:=n_snv/n_ref/n_donor*10^6]


fov_fam2[kmer_11=="GAAACTTCTTT"][order(-rank(mrate))]
rmsk_fam_length

library(stringr)
fov[kmer_11=="GAAACTTCTTT", .N, by=chr][order(as.numeric(str_sub(chr, 4,-1)))]
#####

# Genome-wide mutation rate (5.96 SNV/patient/Mb)
fov_stat2[, sum(n_snv)]/sum(k11c_dt$tot_hg19)/2583*10^6
# RepeatMasker mutation rate (6.06 SNV/patient/Mb)
fov_stat2[!is.na(repClass), sum(n_snv)]/sum(rmsk_length$length)/2583*10^6

fov_stat2[, .(n_snv=sum(n_snv)), by=c("repClass", "length")][, mrate:=n_snv/length*10^6][order(-rank(mrate))]
fov_stat2[kmer_11=="GAAACTTCTTT"]
snv[kmer_11=="GAAACTTCTTT", .N]

library(ggplot2)
library(ggpubr)
library(scales)

### mutation rate ====
fov_fam1 <- merge(fov[kmer_11=="GAAACTTCTTT", .(n_snv=.N, n_donor=uniqueN(Donor_ID)), by=c("kmer_11","repClass")],
                  fov_ref[, .(n_ref=.N), by=repClass], by="repClass")
fov_fam1[, mrate_class:=n_snv/n_ref/n_donor*10^6]
fov_fam1[is.na(repClass),repClass:="Genome-wide"]

# mutation count ====
p <- fov_fam1 %>%
  ggplot(aes(reorder(repClass, n_snv), n_snv,
             label=ifelse(n_snv>=40,paste0(comma(mrate_class,accuracy = 1),"\nSNV/patient/Mb"),""))) +
  geom_col() +
  geom_text(nudge_y = 20) +
  labs(title = "Mutability of GAAAC[T>N]TCTTT", 
       subtitle = "within repetitive elements or outside (Genome-wide)",
       x="region") +
  scale_y_continuous(name="mutation count", breaks = seq(0,1200, length.out = 5),
                     sec.axis=sec_axis(trans=~./1450, name="fraction")) +
  coord_cartesian(clip="off")+
  theme_pubr() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))

ggsave(plot = p, device="pdf", width = 200, height = 100, units = "mm",
       filename = "plots/supplfig3_mut_GAAACTTCTTT_alphasatellite_enrichment.pdf")

# referrence occurrence ====
p_ref <- fov_ref[, .N, by = c("repClass","repFamily","repName")][,
          N_class:=sum(N), by=repClass][is.na(repClass),repClass:="Genome-wide"] %>%
  ggplot(aes(reorder(repClass,N_class), N)) +
  geom_col() +
  labs(title = "Reference occurrence of GAAACTTCTTT", 
       subtitle = "within repetitive elements or outside (Genome-wide)",
       x = "region") +
  scale_y_continuous(name="reference count", 
                     sec.axis=sec_axis(trans=~./10745, name="fraction")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))

ggsave(plot=p_ref, device="pdf", width=200, height=100, units="mm",
       filename="plots/supplfig3_ref_GAAACTTCTTT_alphasatellite_enrichment.pdf")


