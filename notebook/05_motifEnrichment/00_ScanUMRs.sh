#!/bin/bash

BED_DIR="Zea_mays_UMR.bed"
GENOME_FA_DIR="Zm-B73-REFERENCE-GRAMENE-4.0.fa"
COPIES=100 # number of copies for dinucleotide shuffle
THREADS=96 # number of threads for parallel
PREFIX=$(basename ${BED_DIR} .bed)

# get fasta sequences with bedtools
bedtools getfasta -fi ${GENOME_FA_DIR} -bed ${BED_DIR} -fo ${PREFIX}.fa

# get di-nucleotide shuffled sequences
/programs/meme-5.5.2/bin2/fasta-shuffle-letters -kmer 2 -copies ${COPIES} -dna ${PREFIX}.fa > ${PREFIX}_dinucleotide_shuffled_${COPIES}_copies.fasta

# motif scanning of UMR regions
mkdir -p UMR_fimo/${PREFIX}
find ../../../JASPAR2024 -name '*.meme' | parallel -j ${THREADS} fimo -bfile ../../../homogenousNucFreqs.txt --max-strand --no-qvalue --skip-matched-sequence --max-stored-scores 100000000 {} ${PREFIX}.fa '>' UMR_fimo/${PREFIX}/{/.}.fimo.txt

# motif scanning of background sequences
mkdir -p bg_fimo/${PREFIX}
find ../../../JASPAR2024 -name '*.meme' | parallel -j ${THREADS} fimo -bfile ../../../homogenousNucFreqs.txt --max-strand --no-qvalue --skip-matched-sequence --max-stored-scores 100000000 {} ${PREFIX}_dinucleotide_shuffled_${COPIES}_copies.fasta '>' bg_fimo/${PREFIX}/{/.}.fimo.txt

parallel -j ${THREADS} "prefix=\$(basename {} .fimo.txt); bgFile=bg_fimo/${PREFIX}/\${prefix}.fimo.txt; Rscript UMR_enrichment.R {} \${bgFile}" ::: UMR_fimo/${PREFIX}/*.fimo.txt
