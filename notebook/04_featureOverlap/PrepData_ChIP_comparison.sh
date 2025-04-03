# Preps motif and ChIP-seq data for analysis of how motifs predict ChIP peaks in upstream regions
# Charlie Hale, 2025.03.31

# Looking at TFs found both in ChIP-seq (Tu et al 2020) and in JASPAR motif dataset
# 1. EREB_ZmTF117 (MA1818.1, ereb17)
# 2. bZIP_ZmTFLG2 (MA1835.2, lg2)
# 3. GLK_ZmTF58 (MA1830.2, glk53)
# 4. bHLH_ZmTF87 (MA1834.2, bhlh47)
# 5. HB_ZmTF174 (MA1824.2, hb34)

# Set some useful paths to reduce redundancy
full_fimo_file="output/motifOutput/fimo/uncollapsed_5kbUpstream/Zm-B73-REFERENCE-NAM-5.0.tsv"
intermed_dir="output/feature_overlap/chip/intermed_files"
maize_upstream_regions="output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_500upstream_primaryAlignment/Zm-B73-REFERENCE-NAM-5.0.bed"
chip_peak_dir="output/feature_overlap/chip/tf_peaks"

## 1. PREP MOTIFS
# Extract individual maize motifs from maize B73 FIMO output
grep MA1818.1 $full_fimo_file > $intermed_dir/ZmTF117_ereb17_motifs.tsv
grep MA1835.2 $full_fimo_file > $intermed_dir/ZmTFLG2_lg2_motifs.tsv
grep MA1830.2 $full_fimo_file > $intermed_dir/ZmTF58_glk53_motifs.tsv
grep MA1834.2 $full_fimo_file > $intermed_dir/ZmTF87_bhlh47_motifs.tsv
grep MA1824.2 $full_fimo_file > $intermed_dir/ZmTF174_hb34_motifs.tsv

# Convert files to the right format
for file in $intermed_dir/*.tsv; do
  tf_name=$(basename "$file" _motifs.tsv)
  mkdir -p output/feature_overlap/chip/$tf_name
  # convert to bed
  bash src/04_featureOverlap/FIMOtoBED.sh "$file" > "$intermed_dir/${tf_name}_motifs_tmp1.bed"
  # Filter to motifs 500bp upstream of primary alignments
  bedtools intersect \
    -a $intermed_dir/${tf_name}_motifs_tmp1.bed \
    -b $maize_upstream_regions \
    -wa \
    -wb \
    -f 1 \
    > "$intermed_dir/${tf_name}_motifs_tmp2.bed"
    # convert to normal chr numbers
    bash src/06_associationModeling/Convert_ncbi_format.sh \
        $intermed_dir/${tf_name}_motifs_tmp2.bed \
        output/feature_overlap/maize_ref/B73v5_ncbi_key.tsv \
        > output/feature_overlap/chip/$tf_name/${tf_name}_motifs.bed
    # remove tmp files
    rm $intermed_dir/${tf_name}_motifs_tmp*.bed 
done

## 2. Prep upstream regions
bash src/06_associationModeling/Convert_ncbi_format.sh \
  $maize_upstream_regions \
  output/feature_overlap/maize_ref/B73v5_ncbi_key.tsv \
  > output/feature_overlap/B73v5_500up_standardCoords.bed

## 3. Prep ChIP data (downloaded from MaizeGDB, Tu et al 2020 reads remapped to B73 v5)
# ereb17
cp $chip_peak_dir/SRR8525082_region.bed output/feature_overlap/chip/ZmTF117_ereb17/ZmTF117_ereb17_rep1_chip.bed
cp $chip_peak_dir/SRR8525083_region.bed output/feature_overlap/chip/ZmTF117_ereb17/ZmTF117_ereb17_rep2_chip.bed

# lg2
cp $chip_peak_dir/SRR8525052_region.bed output/feature_overlap/chip/ZmTFLG2_lg2/ZmTFLG2_lg2_rep1_chip.bed
cp $chip_peak_dir/SRR8525051_region.bed output/feature_overlap/chip/ZmTFLG2_lg2/ZmTFLG2_lg2_rep2_chip.bed

# glk53
cp $chip_peak_dir/SRR8524997_region.bed output/feature_overlap/chip/ZmTF58_glk53/ZmTF58_glk53_rep1_chip.bed
cp $chip_peak_dir/SRR8524994_region.bed output/feature_overlap/chip/ZmTF58_glk53/ZmTF58_glk53_rep2_chip.bed

# bhlh47
cp $chip_peak_dir/SRR8525165_region.bed output/feature_overlap/chip/ZmTF87_bhlh47/ZmTF87_bhlh47_rep1_chip.bed
cp $chip_peak_dir/SRR8525162_region.bed output/feature_overlap/chip/ZmTF87_bhlh47/ZmTF87_bhlh47_rep2_chip.bed

# hb34
cp $chip_peak_dir/SRR8525010_region.bed output/feature_overlap/chip/ZmTF174_hb34/ZmTF174_hb34_rep1_chip.bed
cp $chip_peak_dir/SRR8525007_region.bed output/feature_overlap/chip/ZmTF174_hb34/ZmTF174_hb34_rep2_chip.bed


# Get list of OGs intersecting ChIP peak union for glk53
cat output/feature_overlap/chip/ZmTF58_glk53/ZmTF58_glk53_rep1_chip.bed \
  output/feature_overlap/chip/ZmTF58_glk53/ZmTF58_glk53_rep2_chip.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  | bedtools intersect \
    -a output/feature_overlap/B73v5_500up_standardCoords.bed \
    -b stdin \
    -wa \
  | cut -f 4 > output/feature_overlap/chip/ZmTF58_glk53/og_overlap_chip_union.txt