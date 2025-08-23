# Removes redundant motif calls by retaining only a single motif instance 
# for overlapping instances of the same JASPAR 2024 cluster
#
# Inputs: 1) Raw motif scanning output (.tsv from FIMO)
#         2) key file connecting motif Matrix IDs to Cluster IDs
#         3) Output directory


# Output: BED file containing non-redundant motif instances

# Charlie Hale, 2024.05.20 chale295@gmail.com

#!/bin/bash

# Accept command line arguments
raw_motifs=$1
motif_key=$2
workDir=$3

# Temporary file for sorted motifs
sorted_motifs=$(mktemp)

# Parse motifs and sort by motif ID
awk '
BEGIN { OFS = "\t" }
NR > 1 && $0 !~ /^#/ && $0 != "" {
  seqName = $3;
  start = $4; # FIMO default is to parse genomic coords and return 0-based coordinates
  stop = $5 + 1; # add 1 to end coordinate to match bed format
  motifID = $1;
  print seqName, start, stop, motifID;
}' "$raw_motifs" | sort -k4,4 > "$sorted_motifs"


# Temporary file for motifs with cluster info
motifs_with_cluster_info=$(mktemp)

key_to_join=$(mktemp)
awk -F '\t' -v  OFS='\t' 'NR > 1 {print $1, $10}' "$motif_key" | sort -k1,1 > "$key_to_join"
join -1 4 -2 1 "$sorted_motifs" "$key_to_join" | \
awk -v OFS='\t' '{print $2, $3, $4,  $1, $5, $6}' > "$motifs_with_cluster_info"
rm "$key_to_join"

# Merge overlapping motif intervals for the same cluster (if needed)
collapsed_motifs=$(mktemp)
fileID=$(basename "$raw_motifs" .tsv)
mkdir -p $workDir/$fileID
if [ -f "$workDir/$fileID/*.bed" ]; then
    rm "$workDir/$fileID/*.bed"
fi

bedtools sort -i $motifs_with_cluster_info |
# Split motif output by motif cluster, generating a separate bed file for each cluster
awk -v workDir="$workDir" -v fileID="$fileID" '{ print > workDir "/" fileID "/" $5 ".bed" }'
for file in "$workDir/$fileID"/*.bed; do
  bedtools cluster -i "$file" -d -1 |
  bedtools groupby -g 1,6 -c 2,3,4,5 -o min,max,collapse,first > "${file%.bed}_clustered_grouped.bed"
done
cat $workDir/$fileID/*_clustered_grouped.bed |
awk -F '\t' -v OFS='\t' '{print $1, $3, $4, $5, $6}' > "$collapsed_motifs"
cat $collapsed_motifs

# Clean up temporary files
rm "$sorted_motifs" "$motifs_with_cluster_info" "$collapsed_motifs"
rm $workDir/$fileID/*.bed
rmdir $workDir/$fileID 