# Extracts motif hits from secondary alignments to allow comparision across paralogs
# Charlie Hale, 2025.08.11

#!/bin/bash
mkdir -p output/turnover_analysis/paralog_correlation

assemblyPath="data/genomes/longread/IWGSC_CS_RefSeq_v2.1.fa"
distance_upstream=500
motifDir="output/motifOutput/fimo/collapsed_5kbUpstream"
unfilt_interval_dir="output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG"
filt_interval_dir="output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_${distance_upstream}_allAlignments"
genomeFileDir="data/genomes/longread/"
assemblyID=$(basename "$assemblyPath" .fa)

mkdir -p "$filt_interval_dir"
# Extract motif hits from all alignments
temp_dir=$(mktemp -d)
temp_dir2=$(mktemp -d)
bash src/04_motifScanning/ExtractUpstreamCoords.sh $assemblyID $distance_upstream $genomeFileDir $unfilt_interval_dir $temp_dir T F
# Filter out intervals containing more than 5% N
bash src/04_motifScanning/FilterRegionsContainingNs.sh $temp_dir/$assemblyID.bed $assemblyPath $filt_interval_dir/$assemblyID.bed 5 $assemblyID

# Filter to only keep OGs that mapped to at least three locations
extract_og_ids() {
  local infile="$1"
  local query="$2"

  grep "$query" "$infile" | \
    awk -F'\t' '{
      match($4, /Target=([^;]+)/, m);
      if (m[1] != "") {
        split(m[1], parts, "_");
        print parts[1];
      }
    }' | \
    sort -u
}
rank1_ogs=$(extract_og_ids "$filt_interval_dir/$assemblyID.bed" "Rank=1")
rank2_ogs=$(extract_og_ids "$filt_interval_dir/$assemblyID.bed" "Rank=2")
rank3_ogs=$(extract_og_ids "$filt_interval_dir/$assemblyID.bed" "Rank=3")
rank4_ogs=$(extract_og_ids "$filt_interval_dir/$assemblyID.bed" "Rank=4")
# Get intersections of OGs across all ranks to recover triads

# Find OGs present in all three ranks
comm -12 <(comm -12 <(echo "$rank1_ogs" | sort -u) <(echo "$rank2_ogs" | sort -u)) <(echo "$rank3_ogs" | sort -u) > "$temp_dir2/triads.txt"

# Filter the intersection file to keep only triad entries
grep -F -f "$temp_dir2/triads.txt" "$filt_interval_dir/$assemblyID.bed" > "$temp_dir2/${assemblyID}_triads.bed" 

# Intersect motifs with upstream intervals, retaining only motifs completely contained within intervals
bedtools intersect -a "$motifDir/$assemblyID.bed" -b "$temp_dir2/${assemblyID}_triads.bed" -wa -wb -f 1 | 

# Clean up file
awk -v assemblyID="$assemblyID" -F'\t' '
BEGIN {
  OFS = "\t";
}

# Function to extract elements from the 9th field (rank, target)
function attr(field9, key,   i,n,parts,pair) {
  n = split(field9, parts, ";");
  for (i = 1; i <= n; i++) {
    split(parts[i], pair, "=");
    gsub(/^ +| +$/, "", pair[1]);
    gsub(/^ +| +$/, "", pair[2]);
    if (pair[1] == key) return pair[2];
  }
  return "";
}



{
  
  target = attr($9, "Target");
  split(target, t, "_"); og = t[1];

  rank = attr($9, "Rank"); if (rank == "") rank = "NA";
  og_rank = og "_Rank" rank;
  print $1, $2, $3, $5, og_rank, assemblyID;
}' | \
sort | \
uniq | \
# Get motif counts for each orthogroup 
awk '{count[$4" "$5]++} END {for (i in count) print i, count[i]}' \
 > output/turnover_analysis/paralog_correlation/${assemblyID}_homologCounts.txt

# Clean up temp file
rm -r $temp_dir
rm -r $temp_dir2