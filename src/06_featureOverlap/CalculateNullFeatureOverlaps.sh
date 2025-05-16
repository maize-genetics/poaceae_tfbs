#!/bin/bash
# Shuffles bed intervals 100 times to generate a null set of intervals to compute motif overlap with

# Usage info
usage() {
  echo "Usage: $0 --null_dir <null_dir> --genome_file <genome_file> --motif_file <motif_file> --feature_file <feature_file>"
  exit 1
}

# Check argument count (4 named arguments = 8 total parameters)
if [ "$#" -ne 8 ]; then
  usage
fi

# Parse named arguments
while [ "$1" != "" ]; do
  case $1 in
    --null_dir )    shift; null_dir=$1 ;;
    --genome_file ) shift; genome_file=$1 ;;
    --motif_file )  shift; motif_file=$1 ;;
    --feature_file )shift; feature_file=$1 ;;
    * )             usage ;;
  esac
  shift
done

shuffled_dir="${null_dir}/shuffled"
overlap_dir="${null_dir}/overlap"

mkdir -p "$null_dir" "$shuffled_dir" "$overlap_dir"

# Export required variables for parallel function
export genome_file motif_file feature_file shuffled_dir overlap_dir

# Define function to perform one shuffle and intersection
process_iteration() {
  i=$1
  # Create filenames using the basename (removing .bed extension)
  shuffled_bed="${shuffled_dir}/$(basename ${feature_file} .bed)_shuffled_${i}.bed"
  overlap_bed="${overlap_dir}/$(basename ${motif_file} .bed)_overlap_null_${i}.bed"

  # Shuffle bed intervals and perform intersection
  bedtools shuffle -i "$feature_file" -g "$genome_file" -noOverlapping > "$shuffled_bed"
  bedtools intersect -a "$motif_file" -b "$shuffled_bed" -u > "$overlap_bed"
  count=$(wc -l < "$overlap_bed") # Count how many overlaps
  echo -e "${i}\t${count}"
}
export -f process_iteration

# Create a summary file
summary_file="${null_dir}/null_overlap_counts.txt"
echo -e "Iteration\tOverlapCount" > "$summary_file"

# Run 100 iterations in parallel
seq 1 100 | parallel -j 40 -k process_iteration {} >> "$summary_file"