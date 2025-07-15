# Scan other maize regulatory features to calculate motif enrichment, using random genomic intervals as background.
# Charlie Hale, 2025.07.15
# scATAC from Marand et al (2021), MOA from Englehorn et al (2024), CNS from Stitzer et al (2025)

#!/bin/bash
mkdir -p output/motif_enrichment/{scATAC,moa,cns} 

COPIES=100 # number of copies for dinucleotide shuffle
THREADS=100 # number of threads for parallel
maize_fasta="data/genomes/motif_enrichment/Zm-B73-REFERENCE-NAM-5.0.fa"

# Scan scATAC from Marand et al 2021
bash src/05_motifEnrichment/ScanFeatureIntervals.sh \
    output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v5.bed \
    ${maize_fasta} \
    "random" \
    data/homogenousNucFreqs.txt \
    output/JASPAR2024_CORE_plants_nr_PFM/pwms \
    output/motif_enrichment/scATAC/fimo \
    $COPIES $THREADS 100

# Scan MOA from Englehorn et al 2024
bash src/05_motifEnrichment/ScanFeatureIntervals.sh \
    output/motif_enrichment/moa/MOA_all_peaks.merged.21_NAMs_new0624.bed \
    ${maize_fasta} \
    "random" \
    data/homogenousNucFreqs.txt \
    output/JASPAR2024_CORE_plants_nr_PFM/pwms \
    output/motif_enrichment/moa/fimo \
    $COPIES $THREADS 100

# Scan CNS from Stitzer et al 2025
bash src/05_motifEnrichment/ScanFeatureIntervals.sh \
    output/motif_enrichment/cns/panand_cns.bed \
    ${maize_fasta} \
    "random" \
    data/homogenousNucFreqs.txt \
    output/JASPAR2024_CORE_plants_nr_PFM/pwms \
    output/motif_enrichment/cns/fimo \
    $COPIES $THREADS 100

# Scan UMRs from Crisp et al 2020 (using genomic background instead of dinucleotide shuffle)
bash src/05_motifEnrichment/ScanFeatureIntervals.sh \
    output/motif_enrichment/umrs_crisp2020/beds/Zea_mays_B73v5.bed \
    ${maize_fasta} \
    "random" \
    data/homogenousNucFreqs.txt \
    output/JASPAR2024_CORE_plants_nr_PFM/pwms \
    output/motif_enrichment/umrs_crisp2020/fimo \
    $COPIES $THREADS 100