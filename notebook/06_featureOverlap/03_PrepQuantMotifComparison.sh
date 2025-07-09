# Calculates motif enrichment in ChIP-seq regions vs di-nucleotide shuffled background sequences.
# Charlie Hale, 2025.07.07

#!/bin/bash
# Looking at TFs found both in ChIP-seq (Tu et al 2020) and in JASPAR motif dataset
# Define a list of TF names, their corresponding JASPAR motifs, and ChIP peak files
# Format: "TF_folder motif_id chip_rep1 chip_rep2"
# tfs=(
#   "ZmTF117_ereb17 MA1818.1 SRR8525082_region.bed SRR8525083_region.bed"
#   "ZmTFLG2_lg2   MA1835.2 SRR8525052_region.bed SRR8525051_region.bed"
#   "ZmTF58_glk53  MA1830.2 SRR8524997_region.bed SRR8524994_region.bed"
#   "ZmTF87_bhlh47 MA1834.2 SRR8525165_region.bed SRR8525162_region.bed"
#   "ZmTF174_hb34  MA1824.2 SRR8525010_region.bed SRR8525007_region.bed"
# )

# # Get consensus of ChIP_seq peaks for glk53
# mkdir -p output/quant_motif_comparison/chip/ZmTF58_glk53/motif
# bedtools intersect -a output/feature_overlap/chip/ZmTF58_glk53/ZmTF58_glk53_rep1_chip.bed \
#     -b output/feature_overlap/chip/ZmTF58_glk53/ZmTF58_glk53_rep2_chip.bed \
#     > output/quant_motif_comparison/chip/ZmTF58_glk53/ZmTF58_glk53_chip_consensus.bed

# #cat output/feature_overlap/chip/ZmTF58_glk53/ZmTF58_glk53_rep1_chip.bed \
#     # output/feature_overlap/chip/ZmTF58_glk53/ZmTF58_glk53_rep2_chip.bed | \
#     # sort -k1,1 -k2,2n | \
#     # bedtools merge -i - > output/quant_motif_comparison/chip/ZmTF58_glk53/ZmTF58_glk53_chip_union.bed

# # Copy the meme motif file for glk53
# cp output/JASPAR2024_CORE_plants_nr_PFM/pwms/MA1830.2.meme output/quant_motif_comparison/chip/ZmTF58_glk53/motif/

# # Scan for motifs with FIMO
# mkdir -p output/quant_motif_comparison/chip/ZmTF58_glk53/fimo

# bash src/06_featureOverlap/ScanFeatureIntervals.sh \
#   output/quant_motif_comparison/chip/ZmTF58_glk53/ZmTF58_glk53_chip_consensus.bed \
#   output/feature_overlap/maize_ref/Zm-B73-REFERENCE-NAM-5.0.fa \
#   data/homogenousNucFreqs.txt \
#   output/quant_motif_comparison/chip/ZmTF58_glk53/motif \
#   output/quant_motif_comparison/chip/ZmTF58_glk53/fimo \
#   100 39

# # Set useful paths
# chip_peak_dir="output/feature_overlap/chip/tf_peaks"



# # Convert upstream region bed file to standard coords 
# bash src/06_featureOverlap/Convert_ncbi_format.sh \
#   "$maize_upstream_regions" \
#   output/feature_overlap/maize_ref/B73v5_ncbi_key.tsv \
#   > output/feature_overlap/B73v5_500up_standardCoords.bed

# # Rename ChIP peak files for simplified parsing
# for tf in "${tfs[@]}"; do
#   read tf_name motif chip1 chip2 <<< "$tf"
#   cp "$chip_peak_dir/$chip1" "output/feature_overlap/chip/${tf_name}/${tf_name}_rep1_chip.bed"
#   cp "$chip_peak_dir/$chip2" "output/feature_overlap/chip/${tf_name}/${tf_name}_rep2_chip.bed"
# done

#!/usr/bin/env bash
set -euo pipefail

tfs=(
  "ZmTF117_ereb17 MA1818.1 SRR8525082_region.bed SRR8525083_region.bed"
  "ZmTFLG2_lg2   MA1835.2 SRR8525052_region.bed SRR8525051_region.bed"
  "ZmTF58_glk53  MA1830.2 SRR8524997_region.bed SRR8524994_region.bed"
  "ZmTF87_bhlh47 MA1834.2 SRR8525165_region.bed SRR8525162_region.bed"
  "ZmTF174_hb34  MA1824.2 SRR8525010_region.bed SRR8525007_region.bed"
)

for entry in "${tfs[@]}"; do
  # split fields
  read -r tf_name motif_id rep1_region rep2_region <<< "$entry"

  echo "Processing $tf_name (motif $motif_id)..."

  # Create outptut directories
  base=output/quant_motif_comparison/chip/$tf_name
  mkdir -p "$base"/{motif,fimo}

  # Get consensus peaks (intersect rep1 vs rep2)
  bedtools intersect \
    -a output/feature_overlap/chip/$tf_name/${tf_name}_rep1_chip.bed \
    -b output/feature_overlap/chip/$tf_name/${tf_name}_rep2_chip.bed \
    > "$base"/${tf_name}_chip_consensus.bed

  # Copy the corresponding PWM
  cp \
    output/JASPAR2024_CORE_plants_nr_PFM/pwms/${motif_id}.meme \
    "$base"/motif/

  # Run FIMO scan on consensus peaks
  bash src/06_featureOverlap/ScanFeatureIntervals.sh \
    "$base"/${tf_name}_chip_consensus.bed \
    output/feature_overlap/maize_ref/Zm-B73-REFERENCE-NAM-5.0.fa \
    data/homogenousNucFreqs.txt \
    "$base"/motif \
    "$base"/fimo \
    100 39

  echo
done
