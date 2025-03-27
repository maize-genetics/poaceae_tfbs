# Converts NCBI GenBank contig IDs from bed file to standard format (chr1,etc.)
# Then join with motif metatadata for IGV visualization
# Charlie Hale. 2024.05.23

# Parse command line arguments
bed_file=$1
genbank_key=$2
#motif_key=$3

# Prepare key file for joining
sorted_genbank_key=$(mktemp)
tail -n +2 "$genbank_key" | awk -F'\t' '{print $6, $14}' | sort -k1,1 > "$sorted_genbank_key"

# Join with bed file
motifs_fixed_contig=$(mktemp)
join -1 1 -2 1 <(sort -k1,1 $bed_file) $sorted_genbank_key |\
awk 'BEGIN { OFS = "\t" }
{
  print $(NF), $2, $3, $5
}' > "$motifs_fixed_contig"

cat $motifs_fixed_contig
# Join with motif metadata
#sorted_motif_key=$(mktemp)
#awk -F '\t' -v OFS='\t' '
#NR > 1 && $0 !~ /^#/ && $0 != "" {
#  print $10, $11, $12
#}' $motif_key | sort -k 1,1 | uniq > "$sorted_motif_key"

#join -t $'\t' -1 4 -2 1 <(sort -k4,4 $motifs_fixed_contig) $sorted_motif_key | \
#awk -F '\t' -v OFS='\t' '
#{
#  print $2, $3, $4, $6, "0", "+", $3, $4, $5
#}' | sort -k1,1 -k2,2n | cat # sort bed file and output to console

# Clean up temp files
rm $sorted_motif_key $motifs_fixed_contig $sorted_genbank_key
