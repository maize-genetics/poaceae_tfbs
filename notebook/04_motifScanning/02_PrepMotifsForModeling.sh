# # Reads in collapsed motif intervals and outputs counts for each species per OG
threads="80"
OG_dir="output/motifOutput/motifs_by_orthogroup"
mkdir -p output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_500upstream_primaryAlignment 
mkdir -p $OG_dir/summarized_by_assembly
mkdir -p $OG_dir/summarized_by_OG
mkdir -p $OG_dir/filtered_OGs_200assemblies

echo "Step 1: Summarizing motifs across assemblies..."
# Summarize motifs by orthogroup, parallelizing over assemblies
parallel -j $threads "bash src/04_motifScanning/SummarizeMotifsByOrthogroup.sh {} \
   output/motifOutput/fimo/collapsed_5kbUpstream \
   output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG \
   output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_500upstream_primaryAlignment \
   500 \
   {//} \
   > $OG_dir/summarized_by_assembly/{/.}.txt" :::: lists/assembly_list.txt 

echo "Step 2: Merging all motif counts together..."
# # Join motif counts from all assemblies together
cat $OG_dir/summarized_by_assembly/*.txt > $OG_dir/mergedMotifCounts_Poaceae800_500upstream_primaryAlignment.txt

echo "Step 3: Summarizing counts by orthogroup..."
# Get counts for each assembly for each orthogroup
grep '^>' output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa | sed 's/>//' > lists/Helixer_OG_IDs_filtered.txt # Get list of OGs
parallel -j "$threads" --no-notice "grep -E '{}($| )' $OG_dir/mergedMotifCounts_Poaceae800_500upstream_primaryAlignment.txt > $OG_dir/summarized_by_OG/{}.txt" :::: lists/Helixer_OG_IDs_filtered.txt

# Remove empty files
find $OG_dir/summarized_by_OG/ -type f -size 0 -delete

echo "Step 4: Filtering out poorly-represented OGs..."

## Retain OGs with at least 200 taxa represented
# Define a function to filter out OGs with fewer than 200 taxa represented
filter_by_assembly_representation() {
    file=$1
    OG_dir=$2
    distinct_count=$(awk '{print $5}' "$file" | sort -u | wc -l)
    if [[ $distinct_count -ge 200 ]]; then
        cp "$file" "$OG_dir/filtered_OGs_200assemblies/"
    fi
}

# Export function and variable for use by parallel
export -f filter_by_assembly_representation
export OG_dir

# Process all files in the summarized_by_OG directory
find "$OG_dir/summarized_by_OG/" -name 'OG*.txt' | parallel -j "$threads" filter_by_assembly_representation {} $OG_dir