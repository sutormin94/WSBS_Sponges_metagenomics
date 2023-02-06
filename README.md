# WSBS_Sponges_metagenomics

An interplay between viruses and bacteria associated with White Sea sponges revealed by metagenomics

This repository contains a set of python scripts which have been used for metagenomic data analysis and visualization.
Original sequencing data were deposited in SRA (SRR17015069-SRR17015085). Metagenome assemblies were deposited in GeneBank under the BioProject PRJNA781598.

If you find this code useful and would like to use it in your own research, please, cite:
Rusanova A, Fedorchuk V, Toshchakov S, Dubiley S, Sutormin D. An Interplay between Viruses and Bacteria Associated with the White Sea Sponges Revealed by Metagenomics. Life (Basel). 2021 Dec 24;12(1):25. doi: 10.3390/life12010025. PMID: 35054418; PMCID: PMC8777954.

## Assembly_average_cov_depth.py

Calculates average coverage depth for de novo metagenome assemblies.

**Requirements:** Python 3

**Input:** Samtools depth output

**Output:** Average coverage depth of assemblies


## AntiPhage_systems_analysis.py

Quantification of antiphage defense systems in metagenomic data.

1) Script takes protein sequences from PADS Arsenal database and groups them by homology to prepare for alignment and preparation of HMM profiles. 
2) Analyses the output of HMM search performed in metagenomic ORFs with HMM profiles prepared for antiphage genes from PADS Arsenal.
3) Analyses the output of PADLOC and DefenseFinder.

**Supplementary data:** The HMM profiles constructed form PADS Arsenal database can be found in PADS_hmm_profiles directory.

**Requirements:** Python 3

**Input:** PADLOC output, DefenseFinder output, contigs coverage depth, hmmsearch output

**Output:** Antiphage defense genes and systems quantification heatmaps


## AntiPhage_systems_detection_in_viral_contigs.py

Quantification of antiphage defense systems in metagenomic data for viral assemblies.

1) Script takes protein sequences from PADS Arsenal database and groups them by homology.
2) Analyses the output of HMM search performed in metagenomic ORFs with HMM profiles prepared for antiphage genes from PADS Arsenal.
3) Analyses the output of PADLOC and DefenseFinder.

**Requirements:** Python 3

**Input:** PADLOC output, DefenseFinder output, contigs coverage depth, hmmsearch output

**Output:** Antiphage defense genes and systems quantification heatmaps


## CRISPR_spacers_and_protospacers_from_viral_contigs.py

Script reads CRISPR spacers blast output and groups it by spacers source and target.

**Requirements:** Python 3

**Input:** Blastn output for CRISPR protospacers search; CRISPR spacers clustering results prepared by CDHit

**Output:** Clustering of spacres by source and target


## CRISPR_spacers_Read_CDhit_MMseqs2_clustering_results.py

Analysis and comparison of CRISPR spacers repertoires from different samples.

1) Script parses MMseqs2 or CDHit clustering results, identifies clusters of CRISPR-Cas spacers.
2) Outputs spacer clusters and intersection of spacer sets derived from different samples.

**Requirements:** Python 3

**Input:** CRISPR spacers clustering results prepared by CDHit or MMseqs2

**Output:** Clustering of spacres, intersections of spacer sets


## Predicted_viral_contigs_data_comparison.py
Polishing of viral contigs prediction in metagenomes assembly.

1) Script takes as input data from CheckV, ViralComplete, and VirSorter2.
2) Filters data from CheckV.
3) Overlaps filtered contigs from CheckV, ViralComplete, and VirSorter2.
4) Returns fasta file with viral contigs of a high confidence (detected by at least 2 pipelines out of 3 used).

**Requirements:** Python 3

**Input:** Outputs of CheckV, ViralComplete, and VirSorter2; viral contig sequences

**Output:** Viral contigs of a high confidence, Venn-diagrams


## Predicted_viral_contigs_data_comparison.py
Comparison of viral contig sets.

1) Script compares the numbers of filtered predicted viral contigs across datasets.
2) Reads results of contigs clustering and draws heatmap representing the datasets overlappings.

**Requirements:** Python 3

**Input:** Viral contig sequences clustering results prepared by CDHit

**Output:** Heatmaps


## Spacers_blast_analysis.py
CRISPR spacers blastn analysis.

Script reads blastn output, filters it, and groups it by spacers clustering.

**Requirements:** Python 3

**Input:** CRISPR spacers sequences clustering results prepared by CDHit; Blastn output for CRISPR protospacers search

**Output:** Filtered and grouped blastn hits