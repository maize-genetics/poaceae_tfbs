#!/bin/bash
# Scans UMR regions for motifs using FIMO and compares to dinucleotide shuffled background sequences across five species.
# Charlie Hale, 2025.07.09

COPIES=100 # number of copies for dinucleotide shuffle
THREADS=39 # number of threads for parallel

mkdir -p output/motif_enrichment/{fastas,umrs_crisp2020,acrs_lu2019}/{fimo,beds}
mkdir -p data/genomes/motif_enrichment
# Download fasta files
if [ ! -f data/genomes/motif_enrichment/Zea_mays.fa ]; then
  curl "https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz" | gzip -cd | sed 's/^>Chr/>/' > data/genomes/motif_enrichment/Zea_mays.fa
fi

# S. bicolor v3.1.1
# O. sativa v7.0
# B. distachyon v3.1
# H vulgare L. (ensembl v42)
# Arabidopsis thaliana (TAIR10)
# Setaria viridis (v1.0)

# Download UMR intervals from supplementary data of Crisp et al 2020 (https://doi.org/10.1073/pnas.2010250117) 
if [ ! -f output/motif_enrichment/umrs_crisp2020/beds/Zea_mays.bed ]; then
  wget https://www.pnas.org/doi/suppl/10.1073/pnas.2010250117/suppl_file/pnas.2010250117.sd03.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > output/motif_enrichment/umrs_crisp2020/beds/Zea_mays.bed
fi
if [ ! -f output/motif_enrichment/umrs_crisp2020/beds/Brachypodium_distachyon.bed ]; then
  wget https://www.pnas.org/doi/suppl/10.1073/pnas.2010250117/suppl_file/pnas.2010250117.sd10.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > output/motif_enrichment/umrs_crisp2020/beds/Brachypodium_distachyon.bed
fi
if [ ! -f output/motif_enrichment/umrs_crisp2020/beds/Hordeum_vulgare.bed ]; then
  wget https://www.pnas.org/doi/suppl/10.1073/pnas.2010250117/suppl_file/pnas.2010250117.sd07.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > output/motif_enrichment/umrs_crisp2020/beds/Hordeum_vulgare.bed
fi
if [ ! -f output/motif_enrichment/umrs_crisp2020/beds/Oryza_sativa.bed ]; then
  wget https://www.pnas.org/doi/suppl/10.1073/pnas.2010250117/suppl_file/pnas.2010250117.sd09.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > output/motif_enrichment/umrs_crisp2020/beds/Oryza_sativa.bed
fi
if [ ! -f output/motif_enrichment/umrs_crisp2020/beds/Sorghum_bicolor.bed ]; then
  wget https://www.pnas.org/doi/suppl/10.1073/pnas.2010250117/suppl_file/pnas.2010250117.sd08.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > output/motif_enrichment/umrs_crisp2020/beds/Sorghum_bicolor.bed
fi

# Create key for UMRs
umr_key=(
  "ZmTF117_ereb17 MA1818.1 SRR8525082_region.bed SRR8525083_region.bed"
  "ZmTFLG2_lg2   MA1835.2 SRR8525052_region.bed SRR8525051_region.bed"
  "ZmTF58_glk53  MA1830.2 SRR8524997_region.bed SRR8524994_region.bed"
  "ZmTF87_bhlh47 MA1834.2 SRR8525165_region.bed SRR8525162_region.bed"
  "ZmTF174_hb34  MA1824.2 SRR8525010_region.bed SRR8525007_region.bed"
)
mkdir -p output/motif_enrichment/umrs_crisp2020/fimo/Zea_mays 
bash src/05_motifEnrichment/ScanFeatureIntervals.sh \
    output/motif_enrichment/umrs_crisp2020/beds/Zea_mays.bed \
    data/genomes/motif_enrichment/Zea_mays.fa \
    data/homogenousNucFreqs.txt \
    output/JASPAR2024_CORE_plants_nr_PFM/pwms \
    output/motif_enrichment/umrs_crisp2020/fimo/Zea_mays \
    100 39

# Download ACR bed files from Lu et al 2019
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE128434&format=file&file=GSE128434%5FBED%5Ffiles%2Etar%2Egz" \
  -O output/motif_enrichment/acrs_lu2019/beds/acrs_lu2019_beds.tar.gz
tar -xzvf output/motif_enrichment/acrs_lu2019/beds/acrs_lu2019_beds.tar.gz -C output/motif_enrichment/acrs_lu2019/beds/
gunzip output/motif_enrichment/acrs_lu2019/beds/*.gz
