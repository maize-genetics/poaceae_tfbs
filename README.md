# poaceae_tfbs
Evolution of cis-regulatory motifs across 589 grass species. MBE paper here: https://doi.org/10.1093/molbev/msaf324

Questions? Contact Charlie Hale (chale295 AT gmail DOT com)

This repo contains code to reproduce all analyses in the manuscript. 

For running code, we recommend referring to the yaml files in the envs/ directory for dependency info. You can also load the environments directly with conda. Use mainEnv.yaml for command line operations unless otherwise noted (exceptions are for running Anchorwave, CrossMap, and asreml-R).

Here's what each subdirectory contains:

00_prepExternalData: Downloads publicly available data from JASPAR (motif PWMs) and NCBI SRA (raw reads used for the 57 contig-level genome assemblies generated for this study). Also pulls data for motif enrichment analyses.

01_shortReadAssembly: Assembles genomes from short reads with Megahit. Pipeline is adapted from Schulz et al 2023: https://doi.org/10.1101/2023.09.19.558246

02_orthogroup: Reconstructs ancestral protein sequences using a representative set of high-quality assemblies spanning diverse grasses, then queries against all 589 assemblies to identify orthologous regions.

03_phyloTreeConstruction: Constructs species tree and phylogenetic kinship matrix.

04_motifEnrichment: Calculates motif enrichment within unmethylated regions and accessible chromatin regions relative to shuffled background sequence.

05_motifScanning: Scans orthologous upstream sequences for JASPAR PWMs, the collapses similar motifs into single intervals and counts occurrences per 500bp upstream region. 

06_motifTurnover: Quantifies shared occupancy of motif instances across species and runs GO enrichments.

07_envirotyping: Pulls occurrence data and associated environmental data and constructs environmental features.

08_associationModeling: Runs motif-environment association models across species, calculate GO enrichments for top orthogroups, then plots top candidates.

