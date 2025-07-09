#!/bin/bash
# Scans UMR regions for motifs using FIMO and compares to dinucleotide shuffled background sequences across five species.
# Charlie Hale, 2025.07.09

COPIES=100 # number of copies for dinucleotide shuffle
THREADS=100 # number of threads for parallel

# Create key for UMRs
umr_key=(
  "Zea_mays Zea_mays.fa"
  "Brachypodium_distachyon Bdistachyon_314_v3.0.fa"
  "Hordeum_vulgare Hordeum_vulgare.fa"
  "Oryza_sativa Osativa_323_v7.0.fa"
  "Sorghum_bicolor Sbicolor_454_v3.0.1.fa"
)

mkdir -p output/motif_enrichment/umrs_crisp2020/fimo/Zea_mays 
bash src/05_motifEnrichment/ScanFeatureIntervals.sh \
    output/motif_enrichment/umrs_crisp2020/beds/Zea_mays.bed \
    data/genomes/motif_enrichment/Zea_mays.fa \
    data/homogenousNucFreqs.txt \
    output/JASPAR2024_CORE_plants_nr_PFM/pwms \
    output/motif_enrichment/umrs_crisp2020/fimo/Zea_mays \
    100 39


