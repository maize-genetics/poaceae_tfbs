# Preps motif and DAP-seq data for analysis of how motifs predict DAP peaks in upstream regions
# Charlie Hale, 2025.03.31

# Looking at TFs found both in DAP-seq (Ricci et al 2020) and in JASPAR motif dataset
# 1. ARF10 (MA1685.2, SRR7889790)
# 2. BZIP72 (MA2412.1, SRR7889811)
# 3. EREB127 (MA2414.1, SRR7889813)
# 4. IG1 (MA2418.1, SRR7889821)
# 5. SBP6 (MA2420.1, SRR7889825)

# Set some useful paths to reduce redundancy
full_fimo_file="output/motifOutput/fimo/uncollapsed_5kbUpstream/Zm-B73-REFERENCE-NAM-5.0.tsv"
intermed_dir="output/feature_overlap/dap/intermed_files"
maize_upstream_regions="output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_500upstream_primaryAlignment/Zm-B73-REFERENCE-NAM-5.0.bed"
dap_peak_dir="output/feature_overlap/dap/Ricci_2019_dapseq_peaks"
mkdir -p $intermed_dir
mkdir -p $dap_peak_dir
## 1. PREP MOTIFS
# Extract individual maize motifs from maize B73 FIMO output
grep MA1685.2 $full_fimo_file > $intermed_dir/ARF10_motifs.tsv
grep MA2412.1 $full_fimo_file > $intermed_dir/BZIP72_motifs.tsv
grep MA2414.1 $full_fimo_file > $intermed_dir/EREB127_motifs.tsv
grep MA2418.1 $full_fimo_file > $intermed_dir/IG1_motifs.tsv
grep MA2420.1 $full_fimo_file > $intermed_dir/SBP6_motifs.tsv

# Convert files to the right format
for file in $intermed_dir/*.tsv; do
  tf_name=$(basename "$file" _motifs.tsv)
  mkdir -p output/feature_overlap/dap/$tf_name
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
        > output/feature_overlap/dap/$tf_name/${tf_name}_motifs.bed
    # remove tmp files
    rm $intermed_dir/${tf_name}_motifs_tmp*.bed 
done

## 2. Prep upstream regions
bash src/06_associationModeling/Convert_ncbi_format.sh \
  $maize_upstream_regions \
  output/feature_overlap/maize_ref/B73v5_ncbi_key.tsv \
  > output/feature_overlap/B73v5_500up_standardCoords.bed

## 3. Prep ChIP data (downloaded from MaizeGDB, Tu et al 2020 reads remapped to B73 v5)
# ARF10
cp $dap_peak_dir/SRR7889790_GEM_peaks.bed output/feature_overlap/dap/ARF10/ARF10_dap.bed

# BZIP72
cp $dap_peak_dir/SRR7889811_GEM_peaks.bed output/feature_overlap/dap/BZIP72/BZIP72_dap.bed

# EREB127
cp $dap_peak_dir/SRR7889813_GEM_peaks.bed output/feature_overlap/dap/EREB127/EREB127_dap.bed

# IG1
cp $dap_peak_dir/SRR7889821_GEM_peaks.bed output/feature_overlap/dap/IG1/IG1_dap.bed

# SBP6
cp $dap_peak_dir/SRR7889825_GEM_peaks.bed output/feature_overlap/dap/SBP6/SBP6_dap.bed
