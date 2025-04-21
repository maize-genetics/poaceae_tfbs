# Converts NCBI GenBank contig IDs from bed file to standard format (chr1,etc.)
# Then join with motif metatadata for IGV visualization
# Charlie Hale. 2024.05.23

# Parse command line arguments
bed_file=$1
genbank_key=$2
motif_key=$3

# Prepare key file for joining
sorted_genbank_key=$(mktemp)
tail -n +2 "$genbank_key" | awk -F'\t' '{print $6, $14}' | sort -k1,1 > "$sorted_genbank_key"

# Join with bed file and extract the Target ID from the attribute field
motifs_fixed_contig=$(mktemp)
join -1 1 -2 1 <(sort -k1,1 "$bed_file") "$sorted_genbank_key" |\
awk 'BEGIN { OFS = "\t" }
{
  # Extract target ID from column 4 (e.g., from "Target=OG0006243_1_247_Zm-B73-REFERENCE-NAM-5.0")
  match($4, /Target=([^_]+)/, arr);
  target = arr[1];
  if (target == "") target = ".";
  # Replace the first column with the extracted Target ID (or append it as needed)
  print $(NF), $2, $3, target, $5, $6 
}' > "$motifs_fixed_contig"

#cat "$motifs_fixed_contig"
# Join with motif metadata
join -t $'\t' -1 5 -2 2 <(sort -k5,5 $motifs_fixed_contig) $motif_key | \
awk -F '\t' -v OFS='\t' '
{
 print $2, $3, $4, $8, "0", "+", $3, $4, $7
}' | sort -k1,1 -k2,2n | cat # sort bed file and output to console

# Clean up temp files
rm $sorted_motif_key $motifs_fixed_contig $sorted_genbank_key
