# Fetches genome assembly and bed files necessary to run motif enrichment analysis.
# Charlie Hale, 2025.07.09

#!/bin/bash
mkdir -p output/motif_enrichment/ref 
mkdir -p output/motif_enrichment/{fastas,umrs_crisp2020,acrs_lu2019}/{fimo,beds}
fasta_dir="data/genomes/motif_enrichment"
mkdir -p $fasta_dir

# Download genome FASTA files for each species
function download_decmp() {
  local url="$1"
  local out_path="$2"
  downloaded=0

  if [ ! -f "$out_path" ]; then
    echo "Downloading $url"
    curl "$url" | gzip --stdout --decompress > "$out_path"
    downloaded=1
  else
    echo "$out_path already exists, skipping"
  fi
}

# Maize v4
if [ ! -f $fasta_dir/Zea_mays.fa ]; then
  curl "https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz" | gzip -cd | sed 's/^>Chr/>/' > $fasta_dir/Zea_mays.fa
fi

# Maize v5
if [ ! -f $fasta_dir/Zm-B73-REFERENCE-NAM-5.0.fa ]; then
  curl "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz" | gzip -cd > $fasta_dir/Zm-B73-REFERENCE-NAM-5.0.fa
fi

download_decmp "ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/hordeum_vulgare/dna/Hordeum_vulgare.IBSC_v2.dna.toplevel.fa.gz" "$fasta_dir/Hordeum_vulgare.fa"
download_decmp "ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz" "$fasta_dir/Zea_mays.AGPv4.dna.toplevel.fa"

#TODO: code to download other genomes from JGI
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


# Download ACR bed files from Lu et al 2019
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE128434&format=file&file=GSE128434%5FBED%5Ffiles%2Etar%2Egz" \
  -O output/motif_enrichment/acrs_lu2019/beds/acrs_lu2019_beds.tar.gz
tar -xzvf output/motif_enrichment/acrs_lu2019/beds/acrs_lu2019_beds.tar.gz -C output/motif_enrichment/acrs_lu2019/beds/
gunzip output/motif_enrichment/acrs_lu2019/beds/*.gz


# Download scATAC from Marand et al 2021
wget https://github.com/plantformatics/maize_single_cell_cis_regulatory_atlas/raw/refs/heads/master/all_ACRs.celltype_calls.overlap.bed.gz \
-O output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v4.bed.gz
gunzip output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v4.bed.gz
awk '{sub(/^chr/, "", $1); print $1, $2, $3}' \
	output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v4.bed \
> output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v4_tmp.bed # Change chr prefix to match chain
# Download chain file
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
 -O output/motif_enrichment/ref/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain
 # Liftover to B73 v5 using CrossMap
export PYTHONPATH=/programs/CrossMap-0.7.3/lib64/python3.9/site-packages:/programs/CrossMap-0.7.3/lib/python3.9/site-packages
export PATH=/programs/CrossMap-0.7.3/bin:$PATH
python /programs/CrossMap-0.7.3/bin/CrossMap bed \
 output/motif_enrichment/ref/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
 output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v4_tmp.bed \
    output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v5_tmp.bed
awk -v OFS='\t' '{print "chr"$1, $2, $3}' \
output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v5_tmp.bed \
> output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v5.bed
sed -i s/chrscaf/scaf/g \
output/motif_enrichment/scATAC/Marand_scATAC_allCellTypes_B73v5.bed

# Lift over UMRs to B73 v5
python /programs/CrossMap-0.7.3/bin/CrossMap bed \
 output/motif_enrichment/ref/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
 output/motif_enrichment/umrs_crisp2020/beds/Zea_mays.bed \
    output/motif_enrichment/umrs_crisp2020/beds/Zea_mays_B73v5_tmp.bed
awk -v OFS='\t' '{print "chr"$1, $2, $3}' \
output/motif_enrichment/umrs_crisp2020/beds/Zea_mays_B73v5_tmp.bed \
> output/motif_enrichment/umrs_crisp2020/beds/Zea_mays_B73v5.bed
sed -i s/chrscaf/scaf/g \
output/motif_enrichment/umrs_crisp2020/beds/Zea_mays_B73v5.bed