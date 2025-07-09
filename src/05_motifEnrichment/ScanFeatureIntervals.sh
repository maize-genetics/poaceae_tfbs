#!/bin/bash

# Scan genomic intervals of interest for motifs using FIMO. Also performs scanning of dinucleotide-shuffled intervals for background motif comparison.
# Charlie Hale, 2025.07.07

BED_FILE=$1 # Path to the BED file containing regions of interest
GENOME_FA=$2 # Path to the genome FASTA file
BFILE=$3 # Background nucleotide frequency file
MOTIF_DIR=$4 # Directory containing JASPAR motif files in .meme format
OUTDIR=$5 # Output directory
COPIES=$6 # number of copies for dinucleotide shuffle
THREADS=$7 # number of threads for parallel
PREFIX=$(basename ${BED_FILE} .bed)

mkdir -p ${OUTDIR}

echo "Extracting sequences from ${BED_FILE} using genome ${GENOME_FA} and saving to ${OUTDIR}/${PREFIX}.fa"
# get fasta sequences with bedtools
bedtools getfasta -fi ${GENOME_FA} -bed ${BED_FILE} -fo ${OUTDIR}/${PREFIX}.fa

echo "Shuffling sequences in ${OUTDIR}/${PREFIX}.fa to create dinucleotide shuffled sequences with ${COPIES} copies"
# get di-nucleotide shuffled sequences
/programs/meme-5.5.2/bin2/fasta-shuffle-letters -kmer 2 -copies ${COPIES} -dna ${OUTDIR}/${PREFIX}.fa > ${OUTDIR}/${PREFIX}_dinucleotide_shuffled_${COPIES}_copies.fasta

# modify sequence names to enable parsing by motif scanner
sed -i 's/:/_/g' ${OUTDIR}/${PREFIX}.fa
sed -i 's/:/_/g' ${OUTDIR}/${PREFIX}_dinucleotide_shuffled_${COPIES}_copies.fasta

# motif scanning of empirical regions
echo "Scanning motifs in empirical regions using FIMO with ${THREADS} threads"
mkdir -p ${OUTDIR}/empirical_fimo/
find ${MOTIF_DIR} -name '*.meme' | parallel -j ${THREADS} fimo -bfile ${BFILE} --max-strand --no-qvalue --skip-matched-sequence --max-stored-scores 100000000000 {} ${OUTDIR}/${PREFIX}.fa '>' ${OUTDIR}/empirical_fimo/{/.}.fimo.txt

# motif scanning of background sequences
echo "Scanning motifs in dinucleotide shuffled background sequences using FIMO with ${THREADS} threads"
mkdir -p ${OUTDIR}/bg_fimo/
find ${MOTIF_DIR} -name '*.meme' | parallel -j ${THREADS} fimo -bfile ${BFILE}  --max-strand --no-qvalue --skip-matched-sequence --max-stored-scores 100000000000 {} ${OUTDIR}/${PREFIX}_dinucleotide_shuffled_${COPIES}_copies.fasta '>' ${OUTDIR}/bg_fimo/{/.}.fimo.txt

