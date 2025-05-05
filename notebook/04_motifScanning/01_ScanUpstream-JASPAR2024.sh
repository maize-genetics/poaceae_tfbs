#!/bin/bash

# Scans upstream regions for JASPAR 2024 plant TF motifs, then bins into consensus motifs
# Charlie Hale, 2024.07.03

# Set number of threads
nthreads=100
input_fasta_dir="output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_5kbUpstream"

# Download PWMs from JASPAR
wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt -O data/JASPAR2024_meme.txt

### FILTER MOTIF LIST to retain motifs with conserved UMR enrichment
# Define input files
meme_file="data/JASPAR2024_meme.txt"
motif_list="output/umr_enrichment/umr_enriched_plant_motifs_JASPAR2024.txt"
filtered_meme_file="data/umr_enriched_plant_motifs_JASPAR2024.meme"

# Start writing to the filtered meme file
echo "MEME version 4" > $filtered_meme_file
echo "" >> $filtered_meme_file
echo "ALPHABET= ACGT" >> $filtered_meme_file
echo "" >> $filtered_meme_file
echo "strands: + -" >> $filtered_meme_file
echo "" >> $filtered_meme_file
echo "Background letter frequencies" >> $filtered_meme_file
echo "A 0.25 C 0.25 G 0.25 T 0.25" >> $filtered_meme_file
echo "" >> $filtered_meme_file

# Read the motif list into an array
mapfile -t motifs < $motif_list

# Check if a given motif ID is in the list
is_motif_in_list() {
  local motif_id="$1"
  for id in "${motifs[@]}"; do
    if [[ "$id" == "$motif_id" ]]; then
      return 0
    fi
  done
  return 1
}

# Process the meme file
in_motif_section=0
current_motif_id=""

while read -r line; do
  if [[ $line == MOTIF* ]]; then
    current_motif_id=$(echo $line | awk '{print $2}')
    in_motif_section=0
    if is_motif_in_list "$current_motif_id"; then
      in_motif_section=1
      echo "$line" >> $filtered_meme_file
    fi
  elif [[ $line == URL* ]]; then
    if ((in_motif_section)); then
      echo "$line" >> $filtered_meme_file
    fi
  elif [[ $line == letter-probability* ]]; then
    if ((in_motif_section)); then
      echo "$line" >> $filtered_meme_file
    fi
  elif [[ $line =~ ^[0-9] ]]; then
    if ((in_motif_section)); then
      echo "$line" >> $filtered_meme_file
    fi
  else
    if ((in_motif_section)); then
      echo "$line" >> $filtered_meme_file
    fi
  fi
done < "$meme_file"

### SCAN FOR MOTIFS
# Prepare list of files to scan
realpath $input_fasta_dir/*.fa > lists/motif_scan_input.txt
# Create output directories
mkdir -p output/motifOutput/fimo/uncollapsed_5kbUpstream
mkdir -p output/motifOutput/fimo/collapsed_5kbUpstream
mkdir -p output/motifOutput/fimo/tmp

# Get start time
start=`date +%s`

# Create log file
log_file="output/motif_scan_log.txt"
# Remove previous file if it exists
if [ -f $log_file ]; then
    rm $log_file
fi

echo "Starting motif scanning process" > $log_file

# Run motif scanning in parallel and log output, including errors
parallel -j $nthreads "
  {
    process_start=\$(date +%s);
    echo Scanning motifs in {/.}... >> $log_file 2>&1;
    fimo -bfile data/homogenousNucFreqs.txt \
      --max-strand \
      --no-qvalue \
      --skip-matched-sequence \
      --text \
      --o output/motifOutput/fimo/tmp/{/.} \
      --max-stored-scores 10000000 \
      $filtered_meme_file \
      {} > output/motifOutput/fimo/uncollapsed_5kbUpstream/{/.}.tsv;
    if [ \$? -ne 0 ]; then
      echo \"Error during FIMO scanning for {/.}\" >> $log_file;
      continue;
    fi
    
    echo Collapsing motifs in {/.}... >> $log_file 2>&1;
    bash src/04_motifScanning/CollapseMotifsByCluster.sh \
    output/motifOutput/fimo/uncollapsed_5kbUpstream/{/.}.tsv \
    data/JASPAR2024_CORE_plants_nr_PFM_key.txt \
    output/motifOutput/fimo/collapsed_5kbUpstream > output/motifOutput/fimo/collapsed_5kbUpstream/{/.}.bed
    if [ \$? -ne 0 ]; then
      echo \"Error during motif collapsing for {/.}\" >> $log_file;
      continue;
    fi
    process_end=\$(date +%s);
    echo Execution time for {/.} was \$(expr \$process_end - \$process_start) seconds >> $log_file
  } >> $log_file 2>&1
" :::: lists/motif_scan_input.txt

echo "Motif scanning process completed" >> $log_file
