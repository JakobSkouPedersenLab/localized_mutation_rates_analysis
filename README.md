# localized_mutation_rates_analysis

Welcome to the GitHub repository for the analysis performed in “**Poulsgaard GA, Sørensen SG, Juul RI, Nielsen MM, Pedersen JS. Sequence dependencies and mutation rates of localized mutational processes in cancer. Genome Med. (2023)**”.

This repository contains the key scripts used for the analysis and to produce the results and figures of the paper. As our paper is not a dedicated Methods paper, we recognize that some adaptations are necessary for successful execution.

To navigate the analysis, please refer to the directory overview below.

## Directory overview
data_cleaning/	 	      # Filter genomic regions in PCAWG dataset and data cleaning
kmer_counting/          # Count 11-mers occurrences in the reference genome and mutations in the PCAWG dataset
kmer_logos/	          	# Run k-mer logo software with custom background 11-mer set
overlap_features/	      # Annotate SNVs with genomic features
simulate_rates/	        # Simulate mutation rates in 11-mer sets and evaluate significance
signature_analysis/	    # Analysis of mutational signatures
- hotspot_enrichment/		# Signature load enrichment in hotspots
- cosine_similarity/		# Compare similarities of mutational signature profiles
genome_stratification/	# Steps to stratify the genome into increasingly smaller parts
overall_rates_seqlogos/	# Plot mutation rate and sequence motif changes across genome stratifications
miscellaneous/	        # Additional analyses (APOBEC, repeat-elements, bimodal mutation rates)
