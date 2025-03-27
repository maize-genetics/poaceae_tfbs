# Filters out intervals containing a user-specified percentage of Ns
input_gff=$1
input_fasta=$2
output_bed=$3
pct_n=$4
assemblyID=$5

# Extract sequence from intervals to enable filtering
tmpDir=$(mktemp -d)
bedtools getfasta -fi $input_fasta -bed $input_gff -bedOut -s > $tmpDir/$assemblyID.bed

# Filter on N content
awk -v threshold="$pct_n" '{
  seq = $7;
  total_length = length(seq);
  n_count = gsub(/[nN]/, "", seq);
  n_percentage = (n_count / total_length) * 100;
  if (n_percentage < threshold) {
    print $0;
  }
}' $tmpDir/$assemblyID.bed > $output_bed

rm -r $tmpDir



