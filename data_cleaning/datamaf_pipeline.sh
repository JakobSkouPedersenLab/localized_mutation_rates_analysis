# Date created:   28-05-2021
# Author:         Gustav Poulsgaard
#
# Purpose:        Modify mutation data with signature annotations (.maf)
#
# Pipeline:       COPY -> UNZIP -> FILTER -> ANNOTATE
#   FILTER:         exclude chrX (chr23) and chrY (chr24)
#                   exclude protein-coding regions (n[autosomal_snv]=343,923)
#   ANNOTATE:       eid (chr_position)
#                   strand-agnostic kmer_11 (substr(ref_context, 6, 11))


# make a COPY of the signature file from the PCAWG/signatures-folder to own PCAWG/Gustav/data/signatures-folder
cp ~/PCAWG/faststorage/signatures/signature.sort.maf.bed.gz ~/PCAWG/faststorage/Gustav/mutprocesses_project/data

# set working directory to own signatures-folder
cd ~/PCAWG/faststorage/Gustav/mutprocesses_project/data

# UNZIP the copy
gunzip signature.sort.maf.bed.gz

# FILTER
# create a subset of the file containing only autosomal SNVs (SNP) and remove protein-coding-regions

#cut -f5 signature.sort.maf.bed | sort | uniq -c
# Variant_Classification (field $6)
# non-protein-coding annotations (n[snv]=41,318,716): 3'UTR, 5'Flank, 5'UTR, IGR, Intron, lincRNA, RNA
#       n[snv]      Variant_Classification
#       251992      3'UTR
#       529668      5'Flank
#       21992526    IGR
#       13686967    Intron
#       3148362     lincRNA
#       55795       5'UTR
#       1653406     RNA
# protein-coding annotations (n[snv]=343,923)
#       n[snv]      Variant_Classification
#       686         De_novo_Start_InFrame
#       1378        De_novo_Start_OutOfFrame
#       221474      Missense_Mutation
#       16747       Nonsense_Mutation
#       300         Nonstop_Mutation
#       91323       Silent
#       11584       Splice_Site
#       431         Start_Codon_SNP

# wc -l signature.sort.maf.bed
# 46,607,263

awk 'BEGIN{OFS=FS="\t"}
$1!="chrX" &&
$1!="chrY" &&
$6=="SNP" && 
$5 ~ /3*UTR|5*Flank|5*UTR|IGR|Intron|lincRNA|RNA/' signature.sort.maf.bed > signature_snv_autosome_noncds.maf

# wc -l signature_snv_autosome_noncds.maf
# 41,318,716

# ANNOTATE

# count occurrences of Chromosome ($1), Start_position ($2) and End_position ($3). 
# For SNVs only one position is necessary
# Not sure why FNR == 0 works...
awk 'BEGIN{FS=OFS="\t"}
NR==FNR {T[$1,$2,$3]++; next} 
FNR==0 {print $0; next}
{print $0, T[$1,$2,$3]}' signature_snv_autosome_noncds.maf signature_snv_autosome_noncds.maf > signature_snv_filtered_anno1.maf

# wc -l signature_snv_filtered_anno1.maf
# 41,318,716


# first, create function revcomp()
# second (BEGIN), create array c that holds each base's complementary partner
# thrid (BEGIN), create header
#   "Chromosome\tStart_position\tEnd_position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tmaja\tTumor_Seq_Allele2\tsnpOverlap\tbyFrequency\tsample\tnormid\tposAndType\tref_context\tVAF\tBI_COMPOSITE_SNV_SBS1_P\tBI_COMPOSITE_SNV_SBS2_P\tBI_COMPOSITE_SNV_SBS3_P\tBI_COMPOSITE_SNV_SBS4_P\tBI_COMPOSITE_SNV_SBS5_P\tBI_COMPOSITE_SNV_SBS6_S\tBI_COMPOSITE_SNV_SBS7a_S\tBI_COMPOSITE_SNV_SBS7b_S\tBI_COMPOSITE_SNV_SBS7c_S\tBI_COMPOSITE_SNV_SBS8_P\tBI_COMPOSITE_SNV_SBS9_P\tBI_COMPOSITE_SNV_SBS10a_S\tBI_COMPOSITE_SNV_SBS11_S\tBI_COMPOSITE_SNV_SBS12_P\tBI_COMPOSITE_SNV_SBS13_P\tBI_COMPOSITE_SNV_SBS14_S\tBI_COMPOSITE_SNV_SBS15_S\tBI_COMPOSITE_SNV_SBS16_P\tBI_COMPOSITE_SNV_SBS17a_P\tBI_COMPOSITE_SNV_SBS17b_P\tBI_COMPOSITE_SNV_SBS18_P\tBI_COMPOSITE_SNV_SBS19_P\tBI_COMPOSITE_SNV_SBS21_S\tBI_COMPOSITE_SNV_SBS22_P\tBI_COMPOSITE_SNV_SBS26_S\tBI_COMPOSITE_SNV_SBS28_P\tBI_COMPOSITE_SNV_SBS30_P\tBI_COMPOSITE_SNV_SBS33_P\tBI_COMPOSITE_SNV_SBS35_P\tBI_COMPOSITE_SNV_SBS36_P\tBI_COMPOSITE_SNV_SBS37_P\tBI_COMPOSITE_SNV_SBS38_S\tBI_COMPOSITE_SNV_SBS39_P\tBI_COMPOSITE_SNV_SBS40_P\tBI_COMPOSITE_SNV_SBS60_P\tBI_COMPOSITE_SNV_SBS55_S\tBI_COMPOSITE_SNV_SBS44_S\tBI_COMPOSITE_SNV_SBS61_S\tBI_COMPOSITE_SNV_SBS62_S\tBI_COMPOSITE_SNV_SBS63_S\tBI_COMPOSITE_SNV_SBS64_P\tBI_COMPOSITE_SNV_SBS65_S\tBI_COMPOSITE_SNV_SBS66_S\tBI_COMPOSITE_SNV_SBS67_S\tBI_COMPOSITE_SNV_SBS68_P\tBI_COMPOSITE_SNV_SBS69_P\tBI_COMPOSITE_SNV_SBS70_P\tBI_COMPOSITE_SNV_SBS71_P\tBI_COMPOSITE_SNV_SBS72_P\tBI_COMPOSITE_SNV_SBS73_S\tBI_COMPOSITE_SNV_SBS74_S\tBI_COMPOSITE_SNV_SBS75_S\tBI_COMPOSITE_SNV_SBS76_S\tBI_COMPOSITE_SNV_SBS77_P\tBI_COMPOSITE_SNV_SBS78_S\tBI_COMPOSITE_SNV_SBS79_S\tBI_COMPOSITE_SNV_SBS80_P\tBI_COMPOSITE_SNV_SBS81_P\tBI_COMPOSITE_SNV_SBS82_P\tBI_COMPOSITE_SNV_SBS83_P\tn_pc"
# fourth, define eid, kmer_11, and center
#         if(center is C or T) keep kmer_11
#         else  take the reverse complement of kmer_11

awk '
function revcomp(arg) {
    o = ""
    for(i = length(arg); i > 0; i--)
        o = o c[substr(arg, i, 1)]
    return(o)
}

BEGIN{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; FS=OFS="\t"
{print "Chromosome","Start_position","End_position","Strand","Variant_Classification",
  "Variant_Type","Reference_Allele","maja","Tumor_Seq_Allele2","snpOverlap",
  "byFrequency","sample","normid","posAndType","ref_context","VAF",
  "BI_COMPOSITE_SNV_SBS1_P","BI_COMPOSITE_SNV_SBS2_P","BI_COMPOSITE_SNV_SBS3_P",
  "BI_COMPOSITE_SNV_SBS4_P","BI_COMPOSITE_SNV_SBS5_P","BI_COMPOSITE_SNV_SBS6_S",
  "BI_COMPOSITE_SNV_SBS7a_S","BI_COMPOSITE_SNV_SBS7b_S","BI_COMPOSITE_SNV_SBS7c_S",
  "BI_COMPOSITE_SNV_SBS8_P","BI_COMPOSITE_SNV_SBS9_P","BI_COMPOSITE_SNV_SBS10a_S",
  "BI_COMPOSITE_SNV_SBS11_S","BI_COMPOSITE_SNV_SBS12_P","BI_COMPOSITE_SNV_SBS13_P",
  "BI_COMPOSITE_SNV_SBS14_S","BI_COMPOSITE_SNV_SBS15_S","BI_COMPOSITE_SNV_SBS16_P",
  "BI_COMPOSITE_SNV_SBS17a_P","BI_COMPOSITE_SNV_SBS17b_P","BI_COMPOSITE_SNV_SBS18_P",
  "BI_COMPOSITE_SNV_SBS19_P","BI_COMPOSITE_SNV_SBS21_S","BI_COMPOSITE_SNV_SBS22_P",
  "BI_COMPOSITE_SNV_SBS26_S","BI_COMPOSITE_SNV_SBS28_P","BI_COMPOSITE_SNV_SBS30_P",
  "BI_COMPOSITE_SNV_SBS33_P","BI_COMPOSITE_SNV_SBS35_P","BI_COMPOSITE_SNV_SBS36_P",
  "BI_COMPOSITE_SNV_SBS37_P","BI_COMPOSITE_SNV_SBS38_S","BI_COMPOSITE_SNV_SBS39_P",
  "BI_COMPOSITE_SNV_SBS40_P","BI_COMPOSITE_SNV_SBS60_P","BI_COMPOSITE_SNV_SBS55_S",
  "BI_COMPOSITE_SNV_SBS44_S","BI_COMPOSITE_SNV_SBS61_S","BI_COMPOSITE_SNV_SBS62_S",
  "BI_COMPOSITE_SNV_SBS63_S","BI_COMPOSITE_SNV_SBS64_P","BI_COMPOSITE_SNV_SBS65_S",
  "BI_COMPOSITE_SNV_SBS66_S","BI_COMPOSITE_SNV_SBS67_S","BI_COMPOSITE_SNV_SBS68_P",
  "BI_COMPOSITE_SNV_SBS69_P","BI_COMPOSITE_SNV_SBS70_P","BI_COMPOSITE_SNV_SBS71_P",
  "BI_COMPOSITE_SNV_SBS72_P","BI_COMPOSITE_SNV_SBS73_S","BI_COMPOSITE_SNV_SBS74_S",
  "BI_COMPOSITE_SNV_SBS75_S","BI_COMPOSITE_SNV_SBS76_S","BI_COMPOSITE_SNV_SBS77_P",
  "BI_COMPOSITE_SNV_SBS78_S","BI_COMPOSITE_SNV_SBS79_S","BI_COMPOSITE_SNV_SBS80_P",
  "BI_COMPOSITE_SNV_SBS81_P","BI_COMPOSITE_SNV_SBS82_P","BI_COMPOSITE_SNV_SBS83_P",
  "n_pc","eid", "kmer_11"
}
}

{
eid=substr($1,4,length($1))"_"$3
kmer_11=toupper(substr($15,6,11))
center=substr(kmer_11,6,1)

if(center ~ /T|C/) print $0,eid,kmer_11
else print $0,eid,revcomp(kmer_11)
}' signature_snv_filtered_anno1.maf > signature_snv_filtered_annotated.maf

# add header is +1 line from last file
# wc -l signature_snv_filtered_annotated.maf
# 41,318,717

# replace empty entries in fields $10 and $11 with NAs
awk 'BEGIN {FS=OFS="\t"} {if($10 ~ /^ *$/) $10 = "NA"; if($11 ~ /^ *$/) $11 = "NA"}; 1' signature_snv_filtered_annotated.maf > signature_snv_filtered_annotated2.maf
# sort maf file by sample id (for joining purposes)
awk '{print $12,$0}' signature_snv_filtered_annotated2.maf | sort > signature_snv_filtered_annotated2sort.maf
# sort donor_labels.tsv by sample id (for joining purposes)
awk 'BEGIN{OFS="\t"} {print $4,$1,$2,$3,$5}' PCAWG_donorinfo/donor_labels.tsv | sort > temp_donor_labels.sort.tsv
# join maf file with donor_labels
join -1 1 -2 1 signature_snv_filtered_annotated2sort.maf temp_donor_labels.sort.tsv > signature_snv_filtered_annotated3.maf
# tab separate file 
awk 'BEGIN{FS=" ";OFS="\t"} {$1=$1;print $0}' signature_snv_filtered_annotated3.maf > signature_snv_filtered_annotated3.tsv
# exclude the duplicate sample id column ($1; keep the one in $12) and create a header
cut -f2-84 signature_snv_filtered_annotated3.tsv | \
awk 'BEGIN{FS=OFS="\t"; print "Chromosome","Start_position","End_position","Strand","Variant_Classification",
  "Variant_Type","Reference_Allele","maja","Tumor_Seq_Allele2","snpOverlap",
  "byFrequency","sample","normid","posAndType","ref_context","VAF",
  "BI_COMPOSITE_SNV_SBS1_P","BI_COMPOSITE_SNV_SBS2_P","BI_COMPOSITE_SNV_SBS3_P",
  "BI_COMPOSITE_SNV_SBS4_P","BI_COMPOSITE_SNV_SBS5_P","BI_COMPOSITE_SNV_SBS6_S",
  "BI_COMPOSITE_SNV_SBS7a_S","BI_COMPOSITE_SNV_SBS7b_S","BI_COMPOSITE_SNV_SBS7c_S",
  "BI_COMPOSITE_SNV_SBS8_P","BI_COMPOSITE_SNV_SBS9_P","BI_COMPOSITE_SNV_SBS10a_S",
  "BI_COMPOSITE_SNV_SBS11_S","BI_COMPOSITE_SNV_SBS12_P","BI_COMPOSITE_SNV_SBS13_P",
  "BI_COMPOSITE_SNV_SBS14_S","BI_COMPOSITE_SNV_SBS15_S","BI_COMPOSITE_SNV_SBS16_P",
  "BI_COMPOSITE_SNV_SBS17a_P","BI_COMPOSITE_SNV_SBS17b_P","BI_COMPOSITE_SNV_SBS18_P",
  "BI_COMPOSITE_SNV_SBS19_P","BI_COMPOSITE_SNV_SBS21_S","BI_COMPOSITE_SNV_SBS22_P",
  "BI_COMPOSITE_SNV_SBS26_S","BI_COMPOSITE_SNV_SBS28_P","BI_COMPOSITE_SNV_SBS30_P",
  "BI_COMPOSITE_SNV_SBS33_P","BI_COMPOSITE_SNV_SBS35_P","BI_COMPOSITE_SNV_SBS36_P",
  "BI_COMPOSITE_SNV_SBS37_P","BI_COMPOSITE_SNV_SBS38_S","BI_COMPOSITE_SNV_SBS39_P",
  "BI_COMPOSITE_SNV_SBS40_P","BI_COMPOSITE_SNV_SBS60_P","BI_COMPOSITE_SNV_SBS55_S",
  "BI_COMPOSITE_SNV_SBS44_S","BI_COMPOSITE_SNV_SBS61_S","BI_COMPOSITE_SNV_SBS62_S",
  "BI_COMPOSITE_SNV_SBS63_S","BI_COMPOSITE_SNV_SBS64_P","BI_COMPOSITE_SNV_SBS65_S",
  "BI_COMPOSITE_SNV_SBS66_S","BI_COMPOSITE_SNV_SBS67_S","BI_COMPOSITE_SNV_SBS68_P",
  "BI_COMPOSITE_SNV_SBS69_P","BI_COMPOSITE_SNV_SBS70_P","BI_COMPOSITE_SNV_SBS71_P",
  "BI_COMPOSITE_SNV_SBS72_P","BI_COMPOSITE_SNV_SBS73_S","BI_COMPOSITE_SNV_SBS74_S",
  "BI_COMPOSITE_SNV_SBS75_S","BI_COMPOSITE_SNV_SBS76_S","BI_COMPOSITE_SNV_SBS77_P",
  "BI_COMPOSITE_SNV_SBS78_S","BI_COMPOSITE_SNV_SBS79_S","BI_COMPOSITE_SNV_SBS80_P",
  "BI_COMPOSITE_SNV_SBS81_P","BI_COMPOSITE_SNV_SBS82_P","BI_COMPOSITE_SNV_SBS83_P",
  "n_pc","eid", "kmer_11","Donor_ID","Tumor_Type","colcode","project_code"
} {print $0}' > signature_snv_filtered_annotated4.tsv


# wc -l signature_snv_filtered_annotated4.tsv
# 41,318,717
head -n 1 signature_snv_filtered_annotated4.tsv | awk '{print NF}'
head -n 2 signature_snv_filtered_annotated4.tsv

# clean up (delete intermediate files)
rm signature_snv_autosome_noncds.maf signature_snv_filtered_anno1.maf
rm signature_snv_filtered_annotated.maf signature_snv_filtered_annotated2.maf
rm signature_snv_filtered_annotated2sort.maf signature_snv_filtered_annotated3.maf
rm signature_snv_filtered_annotated3.tsv temp_donor_labels.sort.tsv
