# Downloaded bulk TF list for maize V5 from Grassius on 2025.02.12
# Link: https://grassius.org/download/Maize_TFome_Bulk_data.xls
wget https://grassius.org/download_species_gene_list/Maize/v5 -O output/turnover_analysis/tf_turnover/grassius_maize_tfs_2025.02.12.csv
# grep TF OGs using Grassius list
awk -F, 'NR > 1 {print $4}' output/turnover_analysis/tf_turnover/grassius_maize_tfs_2025.02.12.csv | grep -v '^$' | grep -Ff - output/OGToZm_mapping_v2.txt > output/turnover_analysis/tf_turnover/TF_OGs.txt