# Compute overlap of motifs with cis-regulatory features (ATAC-seq, MOA-seq, CNS)
# Charlie Hale
# 2024.12.23

mkdir -p output/feature_overlap/maize_ref
mkdir -p output/feature_overlap/maize_ncbi_coords
mkdir -p output/feature_overlap/atac # scATAC from Marand et al (2021)
mkdir -p output/feature_overlap/moa # MOA from Englehortn et al (2024)
mkdir -p output/feature_overlap/cns # CNS from Stitzer et al (in prep)


# First get 500bp upstream regions in maize
bedtools intersect \
-a output/motifOutput/fimo/collapsed_5kbUpstream/Zm-B73-REFERENCE-NAM-5.0.bed \
-b output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_500upstream_primaryAlignment/Zm-B73-REFERENCE-NAM-5.0.bed \
-f 1 \
> output/feature_overlap/maize_ncbi_coords/B73_500up_motifs_ncbiCoords.bed 
 
# Convert to standard coords
bash src/06_featureOverlap/Convert_ncbi_format.sh \
output/feature_overlap/maize_ncbi_coords/B73_500up_motifs_ncbiCoords.bed \
output/feature_overlap/maize_ref/B73v5_ncbi_key.tsv \
> output/feature_overlap/maize_ref/B73_500up_motifs.bed

# What % of motifs overlap MOA peaks?
bedtools intersect \
-a output/feature_overlap/maize_ref/B73_500up_motifs.bed \
-b output/feature_overlap/moa/MOA_all_peaks.merged.21_NAMs_new0624.bed \
-u \
> output/feature_overlap/moa/motifs_moa_overlap.bed # 184078 / 264,024 (69.7%) of motifs overlap with a MOA-seq peak

# WHat % of MOA peaks overlap motifs?
bedtools intersect \
-a output/feature_overlap/moa/MOA_all_peaks.merged.21_NAMs_new0624.bed \
-b output/feature_overlap/maize_ref/B73_500up_motifs.bed -u \
> output/feature_overlap/moa/moa_vs_motif_overlap.bed #14,487 / 701,701 (2.1%)
# scATAC overlap
# Download scATAC from Marand et al 2021
wget https://github.com/plantformatics/maize_single_cell_cis_regulatory_atlas/raw/refs/heads/master/all_ACRs.celltype_calls.overlap.bed.gz \
-O output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v4.bed.gz
gunzip output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v4.bed.gz
awk '{sub(/^chr/, "", $1); print $1, $2, $3}' \
	output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v4.bed \
> output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v4_tmp.bed # Change chr prefix to match chain
# Download chain file
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
 -O output/feature_overlap/ref/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain
 # Liftover to B73 v5 using CrossMap
export PYTHONPATH=/programs/CrossMap-0.7.3/lib64/python3.9/site-packages:/programs/CrossMap-0.7.3/lib/python3.9/site-packages
export PATH=/programs/CrossMap-0.7.3/bin:$PATH
python /programs/CrossMap-0.7.3/bin/CrossMap bed \
 output/feature_overlap/ref/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
 output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v4_tmp.bed \
    output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v5_tmp.bed
awk -v OFS='\t' '{print "chr"$1, $2, $3}' \
output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v5_tmp.bed \
> output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v5.bed 
# Intersect with motifs
bedtools intersect \
-a output/feature_overlap/maize_ref/B73_500up_motifs.bed \
-b output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v5.bed \
-u \
> output/feature_overlap/atac/motifs_atac_overlap.bed # 164,544 / 264,024 (62.3%) of maize motifs overlap with a scATAC peak

# WHat % of ATAC peaks overlap motifs?
bedtools intersect \
	-a output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v5.bed \
	-b output/feature_overlap/maize_ref/B73_500up_motifs.bed -u \
	> output/feature_overlap/atac/atac_vs_motif_overlap.bed #11,365 / 168009 (6.8%)
# PanAnd CNS overlap
bedtools intersect \
-a output/feature_overlap/maize_ref/B73_500up_motifs.bed \
-b output/feature_overlap/cns/panand_cns.bed \
-u \
> output/feature_overlap/cns/motifs_cns_overlap.bed # 71303 / 264,024  (27%) of motifs overlap with a PanAnd CNS

#What % of CNS overlap motifs?
bedtools intersect \
-a output/feature_overlap/cns/panand_cns.bed \
-b output/feature_overlap/maize_ref/B73_500up_motifs.bed -u \
> output/feature_overlap/cns/cns_vs_motif_overlap.bed #36,740 / 1,664,343 (2.2% of CNS overlap motifs)
# CNS vs ATAC overlap
bedtools intersect \
-a output/feature_overlap/cns/panand_cns.bed \
-b output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v5.bed -u \
> output/feature_overlap/cns/cns_atac_overlap.bed # 381,964 / 1,664,343 (22.9%) of CNS overlap with a scATAC peak

#ATAC vs CNS overlap
bedtools intersect \
-a output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v5.bed \
-b output/feature_overlap/cns/panand_cns.bed -u \
> output/feature_overlap/cns/cns_atac_overlap.bed # 84,706 / 168009 (50.4%) of ATAC peaks overlap with CNS
# CNS vs MOA overlap
bedtools intersect \
-a output/feature_overlap/cns/panand_cns.bed \
-b  output/feature_overlap/moa/MOA_all_peaks.merged.21_NAMs_new0624.bed \
-u \
> output/feature_overlap/cns/cns_moa_overlap.bed # 755,788 / 1,664,343 (45.4%) of CNS overlap with a MOA-seq peak

# MOA vs CNS overlap
bedtools intersect \
-a output/feature_overlap/moa/MOA_all_peaks.merged.21_NAMs_new0624.bed \
-b output/feature_overlap/cns/panand_cns.bed -u \
> output/feature_overlap/cns/cns_moa_overlap.bed # 196,897 / 701,701 (28.1%) of MOA peaks overlap with CNS

# Calculate 100 null overlaps for each feature type using shuffled feature intervals
# Generate maize genome file
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -O output/feature_overlap/maize_ref/Zm-B73-REFERENCE-NAM-5.0.fa.gz
gunzip output/feature_overlap/maize_ref/Zm-B73-REFERENCE-NAM-5.0.fa.gz
samtools faidx output/feature_overlap/maize_ref/Zm-B73-REFERENCE-NAM-5.0.fa
genomeFile="output/feature_overlap/maize_ref/Zm-B73-REFERENCE-NAM-5.0.fa.fai"

# MOA-seq
bash src/06_featureOverlap/CalculateNullFeatureOverlaps.sh \
--null_dir output/feature_overlap/moa/null_overlap/ \
--genome_file $genomeFile \
--motif_file output/feature_overlap/maize_ref/B73_500up_motifs.bed \
--feature_file output/feature_overlap/moa/MOA_all_peaks.merged.21_NAMs_new0624.bed

# scATAC-seq
bash src/06_featureOverlap/CalculateNullFeatureOverlaps.sh \
--null_dir output/feature_overlap/atac/null_overlap/ \
--genome_file $genomeFile \
--motif_file output/feature_overlap/maize_ref/B73_500up_motifs.bed \
--feature_file output/feature_overlap/atac/Marand_scATAC_allCellTypes_B73v5.bed 

# CNS
bash src/06_featureOverlap/CalculateNullFeatureOverlaps.sh \
--null_dir output/feature_overlap/cns/null_overlap/ \
--genome_file $genomeFile \
--motif_file output/feature_overlap/maize_ref/B73_500up_motifs.bed \
--feature_file output/feature_overlap/cns/panand_cns.bed