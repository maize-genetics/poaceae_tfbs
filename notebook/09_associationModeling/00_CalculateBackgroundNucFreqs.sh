# Calculates background nucleotide and dinucleotide frequencies in assemblies, based on introns < 150bp
# Charlie Hale, 2024.07.12

threads="39"
gff_dir="output/miniProt_alignments/OGconsensus"

mkdir -p output/background_nuc_freqs/introns_filtered_150bp
mkdir -p output/background_nuc_freqs/nucFreqs


    # Pull intron intervals from miniprot alignments
parallel -j $threads "bash src/09_associationModeling/pullIntrons.sh $gff_dir/{}.gff | \
awk -v OFS='\t' '{if((\$3 - \$2) < 150) print \$0}' > output/background_nuc_freqs/introns_filtered_150bp/{}.bed" :::: lists/assembly_list.txt

parallel -j 2 "bedtools getfasta -fi {} -bed output/background_nuc_freqs/introns_filtered_150bp/{/.}.bed -s > output/background_nuc_freqs/introns_filtered_150bp/{/.}.fa" :::: lists/all_assembly_paths.txt
    # Calculate nucleotide frequecies for each assembly
parallel -j $threads "fasta-get-markov -m 1 output/background_nuc_freqs/introns_filtered_150bp/{}.fa > output/background_nuc_freqs/nucFreqs/{}.txt" :::: lists/assembly_list.txt

# Clean up and aggregate into a single CSV
Rscript src/09_associationModeling/NucFreqsToCSV.R output/background_nuc_freqs/nucFreqs/ output/modelData/nucFreqs_shortIntrons.csv
