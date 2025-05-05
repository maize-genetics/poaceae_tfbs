# Run Kraken on raw read subsamples to check for contamination
# Charlie Hale, 2024.04.29

mkdir -p output/kraken/concat_reads
mkdir -p output/kraken/subsampled_reads
mkdir -p output/kraken/kraken_results
mkdir -p output/kraken/kraken_PlusPFP
mkdir -p output/kraken/krona 

# Prepare DB for Kraken using plsPFP database
#cd output/kraken
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20240112.tar.gz -O output/kraken/k2protocol_db/k2_pluspfp_20240112.tar.gz
tar -xzvf output/kraken/k2protocol_db/k2_pluspfp_20240112.tar.gz -C output/kraken/kraken_PlusPFP/
# Update taxonomy for Krona
cd ~/miniconda3/envs/androMotifs/opt/krona/
./updateTaxonomy.sh
cd /workdir/coh22/andro_tfbs
# Run Kraken on each accession and generate krona visualization
while IFS= read -r accession; do
    Concatenate forward and reverse reads into a single file
   find data/sra/ -maxdepth 2 -type f -name "${accession}*" -exec cat {} + > output/${accession}.fastq.gz
   echo "Concatenated $accession."
   gunzip output/${accession}.fastq.gz 
    subsample reads to a depth of 10M
   reformat.sh in=output/${accession}.fastq out=output/subsampled_reads/${accession}_subsampled.fastq reads=10000000
    Run Kraken
    kraken2 --db output/kraken/kraken_PlusPFP output/kraken/subsampled_reads/${accession}_subsampled.fastq --output output/kraken/kraken_results/${accession}_kraken_output.txt --report output/kraken/kraken_results/${accession}_kraken_report.txt --thread 110
    Generate krona plots
    ktImportTaxonomy -t 5 output/kraken/kraken_results/${accession}_kraken_report.txt -o output/kraken/kraken_results/${accession}_krona_output.html
    /programs/KrakenTools/kreport2krona.py -r output/kraken/kraken_results/${accession}_kraken_report.txt -o output/kraken/krona/${accession}.report.krona
done < lists/kraken_list.txt 

# Generate HTML page with krona results from multiple samples
cd output/kraken/krona/
python3 ../../../src/01_shortReadAssembly/BuildKronaHTMLpage.py SRA
cd ../../../..