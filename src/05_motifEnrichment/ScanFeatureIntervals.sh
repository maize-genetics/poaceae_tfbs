#!/bin/bash

# Scan genomic intervals of interest for motifs using FIMO. Also performs scanning of dinucleotide-shuffled intervals for background motif comparison.
# Charlie Hale, 2025.07.07

#!/usr/bin/env bash
set -eo pipefail
# Usage: ScanFeatureIntervals.sh <bed_file> <genome_fa> <bfile> <meme_dir> <outdir> <copies> <threads> [num_batches]

BED_FILE=$1        # e.g. regions.bed
GENOME_FA=$2       # e.g. genome.fa
BG_TYPE=${3:-"dinuc"}  # background type, "dinuc" for dinucleotide shuffle or "random" for random genomic intervals
BFILE=${4:-"data/homogenousNucFreqs.txt"} # background nucleotide freq
MOTIF_DIR=$5       # dir of .meme files
OUTDIR=$6          # base output directory
COPIES=$7          # total dinuc-shuffle copies
THREADS=$8         # total FIMO threads
NUM_BATCHES=${9:-1}  # number of batches to split shuffled sequences into (default: 1)
PREFIX=$(basename "${BED_FILE}" .bed)
BASE_OUT="${OUTDIR}/${PREFIX}"

# 0) make top‐level and subdirs
mkdir -p "${BASE_OUT}/empirical_fimo/logs" \
         "${BASE_OUT}/${BG_TYPE}_bg_fimo/logs"

# 1) extract sequences for bed intervals
echo "1) Extracting → ${BASE_OUT}/${PREFIX}.fa"
bedtools getfasta -fi "${GENOME_FA}" -bed "${BED_FILE}" \
    -fo "${BASE_OUT}/${PREFIX}.fa"

# 2) scan empirical sequences
echo "2) Scanning empirical with ${THREADS} threads"
find "${MOTIF_DIR}" -name '*.meme' | parallel -j "${THREADS}" \
  "fimo -bfile ${BFILE} --max-strand --no-qvalue --skip-matched-sequence \
        --max-stored-scores 100000000000 {} ${BASE_OUT}/${PREFIX}.fa \
    > ${BASE_OUT}/empirical_fimo/{/.}.fimo.txt \
    2> ${BASE_OUT}/empirical_fimo/logs/{/.}.log"

# 3) Generate background intervals and scan them for motifs
if [[ "${BG_TYPE}" == "random" ]]; then
  echo "3) Generating random genomic intervals"
  mkdir -p "${BASE_OUT}/bg_random"
  # Create genome file index if not already present
  if [[ ! -f "${GENOME_FA}.fai" ]]; then
    echo "Indexing genome file..."
    samtools faidx "${GENOME_FA}"
  fi
  # Generate random genomic intervals using bedtools shuffle
  echo "Generating ${COPIES} random genomic intervals from ${BED_FILE}"
  seq 1 "${NUM_BATCHES}" | parallel -j "${THREADS}" \
    "bedtools shuffle -i ${BED_FILE} -g ${GENOME_FA}.fai  \
      > ${BASE_OUT}/bg_random/${PREFIX}_batch{}.bed; \
      bedtools getfasta -fi ${GENOME_FA} -bed ${BASE_OUT}/bg_random/${PREFIX}_batch{}.bed \
      -fo ${BASE_OUT}/bg_random/${PREFIX}_batch{}.fa"

  #  sanitize headers to enable parsing by FIMO
  echo "Sanitizing names"
  sed -i 's/:/_/g' "${BASE_OUT}/${PREFIX}.fa" \
              "${BASE_OUT}/bg_random/${PREFIX}_batch"*.fa

  # scan random sequences
  echo "Scanning random background sequences with ${THREADS} threads"
  parallel -j "${THREADS}" \
    "fimo -bfile ${BFILE} --max-strand --no-qvalue --skip-matched-sequence \
          --max-stored-scores 100000000000 {1} {2} \
      > ${BASE_OUT}/${BG_TYPE}_bg_fimo/{1/.}_{2/.}.fimo.txt \
      2> ${BASE_OUT}/${BG_TYPE}_bg_fimo/logs/{1/.}_{2/.}.log" \
    ::: "${MOTIF_DIR}"/*.meme \
    ::: "${BASE_OUT}/bg_random/${PREFIX}_batch"*.fa

# elif [[ "${BG_TYPE}" == "dinuc" ]]; then
elif [[ "${BG_TYPE}" == "dinuc" ]]; then
  echo " Using dinucleotide shuffle"
  mkdir -p "${BASE_OUT}/shuffled_batches"
  #shuffle intervals, preserving dinucleotide frequencies
  echo "Shuffling into ${NUM_BATCHES} batches of $((COPIES/NUM_BATCHES)) copies"
  seq 1 "${NUM_BATCHES}" | parallel -j "${THREADS}" \
    "/programs/meme-5.5.2/bin2/fasta-shuffle-letters \
      -kmer 2 -copies $((COPIES/NUM_BATCHES)) -dna ${BASE_OUT}/${PREFIX}.fa \
    > ${BASE_OUT}/shuffled_batches/${PREFIX}_batch{}.fa"

  # 4) sanitize headers to enable parsing by FIMO
  echo "Sanitizing names"
  sed -i 's/:/_/g' "${BASE_OUT}/${PREFIX}.fa" \
              "${BASE_OUT}/shuffled_batches/${PREFIX}_batch"*.fa
  # scan shuffled sequences
  echo "Scanning shuffled batches with ${THREADS} threads"
  parallel -j "${THREADS}" \
    "fimo -bfile ${BFILE} --max-strand --no-qvalue --skip-matched-sequence \
          --max-stored-scores 100000000000 {1} {2} \
      > ${BASE_OUT}/bg_fimo/{1/.}_{2/.}.fimo.txt \
      2> ${BASE_OUT}/bg_fimo/logs/{1/.}_{2/.}.log" \
    ::: "${MOTIF_DIR}"/*.meme \
    ::: "${BASE_OUT}/shuffled_batches/${PREFIX}_batch"*.fa

else
  echo "Error: Invalid background type '${BG_TYPE}'. Use 'dinuc' or 'random'."
  exit 1
fi


# # 6) Filter out negative scores
# echo "6) Filtering out negative scores"
# find "${BASE_OUT}/empirical_fimo" -name '*.fimo.txt' | parallel -j "${THREADS}" \
#   "awk 'NR==1||\$7+0>=1' {} > ${BASE_OUT}/empirical_fimo/filt_neg/{/}"

# find "${BASE_OUT}/bg_fimo" -name '*.fimo.txt' | parallel -j "${THREADS}" \
#   "awk 'NR==1||\$7+0>=1' {} > ${BASE_OUT}/bg_fimo/filt_neg/{/}"
# echo "Done. Outputs in ${BASE_OUT}/"