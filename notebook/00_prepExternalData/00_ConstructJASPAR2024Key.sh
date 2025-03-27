# Pulls and aggregates data from JASPAR 2024 webpage to create a key file
# matching JASPAR motif IDs to TF metadata, cluster info, etc
# Charlie Hale, 2024.05.14

# Scrape motif profile metadatafrom JASPAR (Travis's code)
mkdir -p output/JASPAR2024_CORE_plants_nr_PFM/meta
mkdir -p output/JASPAR2024_CORE_plants_nr_PFM/pwms
wget 'https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_meme.zip' -P output/JASPAR2024_CORE_plants_nr_PFM/
unzip output/JASPAR2024_CORE_plants_nr_PFM/JASPAR2024_CORE_plants_non-redundant_pfms_meme.zip -d output/JASPAR2024_CORE_plants_nr_PFM/pwms/

echo "Scraping motif metadata from JASPAR..."
parallel --jobs 6 --joblog output/JASPAR2024_CORE_plants_nr_PFM/JASPAR.meta.joblog --resume-failed --bar --halt 'soon,fail=1' 'sleep 1 && curl "https://jaspar.elixir.no/api/v1/matrix/{/.}/?format=json" > output/JASPAR2024_CORE_plants_nr_PFM/meta/{/.}.json' ::: output/JASPAR2024_CORE_plants_nr_PFM/pwms/*.meme

# Convert jsons to a single merged tab-separated file
#tar -xzvf output/jaspar2024/JASPAR2024_CORE_plants_nr_PFM_meta.tgz -C output/jaspar2024
echo "Aggregating motif metadata..."
output_file="output/JASPAR2024_CORE_plants_nr_PFM/JASPAR2024_CORE_plants_nr_PFM_aggregated.tsv"
rm $output_file # Remove file if it already exists
# Add header row
#echo -e "MatrixID\tName\tClass\tFamily\tType\tUniprotID\tPazarTFID\tSpeciesName\tBaseID\tTaxGroup" > "$output_file"
for json in output/JASPAR2024_CORE_plants_nr_PFM/meta/*.json; do
    jq -r '. | [
    .matrix_id,
    .name,
    (.class | join(", ")),
    (.family | join(", ")),
    .type,
    (.uniprot_ids | join(", ")),
    ([.species] | flatten | map(.name) | join(", ")),
    .base_id,
    .tax_group
] | @tsv' "$json" >> $output_file
done

# Download clustering info from JASPAR
wget 'https://jaspar.elixir.no/static/clustering/2024/plants/CORE/interactive_trees/clusters.tab' -P output/JASPAR2024_CORE_plants_nr_PFM/ 
# Create cluster mapping file
awk -F '\t' 'NR > 1 { # Skip header
    split($2, ids, ","); # Split IDs by comma
    for (i in ids) {
        split(ids[i], parts, "_"); # Split each ID to extract MatrixID
        print parts[4] "\t" $1; # Print MatrixID and cluster
    }
}' output/JASPAR2024_CORE_plants_nr_PFM/clusters.tab | sort -k 2,2 > output/JASPAR2024_CORE_plants_nr_PFM/motifID_to_cluster.txt

# Annotate clusters with unique color mappings and descriptive names
sort -k2,2 data/JASPAR2024_SummarizedClusters.txt -o data/JASPAR2024_SummarizedClusters.txt # Assigned colors/descriptors by hand
sort -k2,2 output/JASPAR2024_CORE_plants_nr_PFM/motifID_to_cluster.txt -o output/JASPAR2024_CORE_plants_nr_PFM/motifID_to_cluster.txt
join -1 2 -2 2 -t $'\t' -a 1 output/JASPAR2024_CORE_plants_nr_PFM/motifID_to_cluster.txt data/JASPAR2024_SummarizedClusters.txt | sort -k 2,2 > output/JASPAR2024_CORE_plants_nr_PFM/annotated_clusters.txt
# Join with motif metadata
join -1 1 -2 2 -t $'\t' -a 1 <(sort -k1,1 $output_file) output/JASPAR2024_CORE_plants_nr_PFM/annotated_clusters.txt > output/JASPAR2024_CORE_plants_nr_PFM/tmp_motif_join.txt
# Add headers
echo -e "MatrixID\tName\tClass\tFamily\tType\tUniprotID\tSpeciesName\tBaseID\tTaxGroup\tCluster\tClusterColor\tClusterDesc" | cat - output/JASPAR2024_CORE_plants_nr_PFM/tmp_motif_join.txt > data/JASPAR2024_CORE_plants_nr_PFM_key.txt