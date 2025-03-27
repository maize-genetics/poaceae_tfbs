# Query Poaceae short read WGS accessions from NCBI database
# CHarlie Hale, 2024.01.04

# Accept arguments for memory and thread usage
# Default values
default_memory=400 # Default memory in GB
default_threads=110 # Default number of threads

# Use provided arguments or default values
memory=${1:-$default_memory}
threads=${2:-$default_threads}

# Download files to user repository
vdb-config --prefetch-to-user-repo

# Load accession list to iterate over
list="lists/accessions_with_run_accessions.txt"

mkdir -p data/sra
# Iteratively download accessions
while IFS=, read -r exp_accession run_accession
do
	echo "Downloading .sra file for ${run_accession}..."
	prefetch $run_accession --max-size u -O data/sra/${run_accession}/
	# Ensure that the file is downloaded correctly
	if [ ! -f "data/sra/${run_accession}/${run_accession}.sra" ] && [ ! -f "data/sra/${run_accession}/${run_accession}/${run_accession}.sra" ]; then
		echo "File not found: data/sra/${run_accession}/${run_accession}.sra"
		continue
	fi
	echo "Extracting fastq files for ${run_accession}..."
	fasterq-dump $run_accession --outdir data/sra/ --mem "${memory}G" --threads $threads --bufsize 500M --curcache 1G

	echo "Compressing fastqs for ${run_accession}..."
	find data/sra -mindepth 1 -maxdepth 2 -type f -name "${run_accession}"*.fastq | xargs -n 1 -P $threads pigz
	rm -r data/sra/${run_accession}
done < "$list"

