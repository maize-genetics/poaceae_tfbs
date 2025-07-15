#!/bin/bash
# Scans UMR regions for motifs using FIMO and compares to dinucleotide shuffled background sequences across five species.
# Charlie Hale, 2025.07.09

COPIES=100 # number of copies for dinucleotide shuffle
THREADS=110 # number of threads for parallel

# Create key for UMRs
umr_key=(
  "Zea_mays Zea_mays.fa"
  "Brachypodium_distachyon Bdistachyon_314_v3.0.fa"
  "Hordeum_vulgare Hordeum_vulgare.fa"
  "Oryza_sativa Osativa_323_v7.0.fa"
  "Sorghum_bicolor Sbicolor_454_v3.0.1.fa"
)

for species in "${umr_key[@]}"; do
  read -r prefix fasta <<< "$species"
  echo "STARTING ${prefix}" 
  bash src/05_motifEnrichment/ScanFeatureIntervals.sh \
      output/motif_enrichment/umrs_crisp2020/beds/${prefix}.bed \
      data/genomes/motif_enrichment/${fasta} \
      dinuc \
      data/homogenousNucFreqs.txt \
      output/JASPAR2024_CORE_plants_nr_PFM/pwms \
      output/motif_enrichment/umrs_crisp2020/fimo \
      $COPIES $THREADS 100
  echo "FINISHED ${prefix}"
done

# Create key for UMRs
acr_key=(
  "Zea_mays Maize_7days_leaf_ACRs.bed Zea_mays.AGPv4.dna.toplevel.fa"
  "Brachypodium_distachyon Brachypodium_7days_leaf_ACRs.bed Bdistachyon_314_v3.0.fa"
  "Hordeum_vulgare Barley_7days_leaf_ACRs.bed Hordeum_vulgare.fa"
  "Oryza_sativa Rice_7days_leaf_ACRs.bed Osativa_323_v7.0.fa"
  "Sorghum_bicolor Sorghum_7days_leaf_ACRs.bed Sbicolor_454_v3.0.1.fa"
  "Setaria_viridis Setaria_7days_leaf_ACRs.bed Sviridis_311_v1.0.fa"
  "Arabidopsis_thaliana Arabidopsis_7days_leaf_ACRs.bed Athaliana_447_TAIR10.fa"
)

# for species in "${acr_key[@]}"; do
#   read -r prefix bed fasta <<< "$species"
#   echo "STARTING ${prefix}" 
#   bash src/05_motifEnrichment/ScanFeatureIntervals.sh \
#       output/motif_enrichment/umrs_crisp2020/beds/${bed} \
#       data/genomes/motif_enrichment/${fasta} \
#       data/homogenousNucFreqs.txt \
#       output/JASPAR2024_CORE_plants_nr_PFM/pwms \
#       output/motif_enrichment/umrs_crisp2020/fimo \
#       $COPIES $THREADS
#   echo "FINISHED ${prefix}"
# done



