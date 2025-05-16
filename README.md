# poaceae_tfbs
Evolution of cis-regulatory motifs across 589 grass species. Preprint here: https://www.biorxiv.org/content/10.1101/2025.04.23.650228v2

This repo contains code to reproduce all analyses in the manuscript. Here's what each subdirectory contains:
00_prepExternalData: Downloads publicly available data from JASPAR (motif PWMs) and NCBI SRA (raw reads used for Megahit genome assemblies).
01_shortReadAssembly: Assembles genomes from short reads with Megahit. Pipeline is adapted from Schulz et al 2023: https://doi.org/10.1101/2023.09.19.558246
02_orthgroup: Reconstructs ancestral protein sequences using a representative set of high-quality assemblies spanning diverse grasses, then queries against all 589 assemblies to identify orthologous regions.
03_phyloTreeConstruction: Constructs species tree and phylogenetic kinship matrix.
04_motifScanning: Scans orthologous upstream sequences for JASPAR PWMs, the collapses similar motifs into single intervals and counts occurrences per 500bp upstream region. Note: the current iteration uses only the UMR-enriched motifs from 05, so 05 should actually be run first.
05_motifEnrichment: Calculates motif enrichment within unmethylated regions and accessible chromatin regions relative to shuffled background sequence.
06_featureOverlap: Calculates proportion of maize motifs overlapping known regulatory features (CNS, MOA-seq, scATAC-seq, ChIP-seq).
07_motifTurnover: Quantifies conservation of motif instances across species and runs GO enrichments.
08_envirotyping: Pulls occurrence data and associated environmental data and constructs environmental features.
09_associationModeling: Runs motif-environment association models across species, calculate GO enrichments for top orthogroups, then plots top candidates.


