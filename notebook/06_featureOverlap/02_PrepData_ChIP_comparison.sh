# Preps motif and ChIP-seq data for analysis of how motifs predict ChIP peaks in upstream regions
# Charlie Hale, 2025.03.31

#!/bin/bash

# Set useful paths
maize_upstream_regions="output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_500upstream_primaryAlignment/Zm-B73-REFERENCE-NAM-5.0.bed"
chip_peak_dir="output/feature_overlap/chip/tf_peaks"

# Looking at TFs found both in ChIP-seq (Tu et al 2020) and in JASPAR motif dataset
# Define a list of TF names, their corresponding JASPAR motifs, and ChIP peak files
# Format: "TF_folder motif_id chip_rep1 chip_rep2"
tfs=(
  "ZmTF117_ereb17 MA1818.1 SRR8525082_region.bed SRR8525083_region.bed"
  "ZmTFLG2_lg2   MA1835.2 SRR8525052_region.bed SRR8525051_region.bed"
  "ZmTF58_glk53  MA1830.2 SRR8524997_region.bed SRR8524994_region.bed"
  "ZmTF87_bhlh47 MA1834.2 SRR8525165_region.bed SRR8525162_region.bed"
  "ZmTF174_hb34  MA1824.2 SRR8525010_region.bed SRR8525007_region.bed"
)

# Convert upstream region bed file to standard coords 
bash src/06_featureOverlap/Convert_ncbi_format.sh \
  "$maize_upstream_regions" \
  output/feature_overlap/maize_ref/B73v5_ncbi_key.tsv \
  > output/feature_overlap/B73v5_500up_standardCoords.bed

# Rename ChIP peak files for simplified parsing
for tf in "${tfs[@]}"; do
  read tf_name motif chip1 chip2 <<< "$tf"
  cp "$chip_peak_dir/$chip1" "output/feature_overlap/chip/${tf_name}/${tf_name}_rep1_chip.bed"
  cp "$chip_peak_dir/$chip2" "output/feature_overlap/chip/${tf_name}/${tf_name}_rep2_chip.bed"
done

# Identify OGs with ChIP overlap (intersecting ChIP peak union)
for tf in "${tfs[@]}"; do
  read tf_name motif chip1 chip2 <<< "$tf"
  echo "Processing ChIP peaks for $tf_name"
  
  cat "output/feature_overlap/chip/${tf_name}/${tf_name}_rep1_chip.bed" \
      "output/feature_overlap/chip/${tf_name}/${tf_name}_rep2_chip.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    | bedtools intersect -a output/feature_overlap/B73v5_500up_standardCoords.bed -b stdin -wa \
    | cut -f 4 > "output/feature_overlap/chip/${tf_name}/og_overlap_chip_union.txt"
done