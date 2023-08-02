

setwd("/faststorage/project/PCAWG/Gustav/")
library(data.table)

cs_dt <- fread("data/signatures/signature_snv_autosome_n_project_header_filt_annot.tsv", key = "Donor_ID")
sbs_columns <- colnames(cs_dt)[grepl("BI_COMPOSITE_SNV_SBS", colnames(cs_dt))]
  
posmean <- cs_dt[, lapply(.SD, mean), .SDcols = sbs_columns, by = c("eid","rec")]
recmean <- posmean[, lapply(.SD, mean), .SDcols = sbs_columns, by = c("rec")]
fwrite(recmean, file = "mutprocesses_project/results/signaturemean_by_recurrence.tsv")
