# Align reconstructed ancestral protein sequences to assemblies to pull out sequence 5kb upstream of orthologous genes
# Charlie Hale, 2024.10.09
# chale295@gmail.com

threads="100"

mkdir -p output/miniProt_alignments/unfiltered
mkdir -p output/miniProt_alignments/filtered_mRNA_stop_frameshift
mkdir -p output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG
mkdir -p output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_5kbUpstream/

# Rename shortreads for cleaner processing
 for file in data/genomes/shortread/*.final.contigs.fa; do
 	    mv "$file" "${file/.final.contigs.fa/.fa}"
 done

# Generate path list
find data/genomes/longread -type f -name "*.fa"  > lists/longreads.txt
find data/genomes/shortread -type f -name "*.fa" > lists/shortreads.txt
cat lists/longreads.txt lists/shortreads.txt > lists/assembly_list.txt

nJobs=$((threads / 40))
echo "Running miniprot alignments..."
parallel -j $nJobs "miniprot --gff-only -t 40 {} output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa > output/miniProt_alignments/unfiltered/{/.}.gff" :::: lists/assembly_list.txt

#Filter to retain  alignments without stop codons or frameshifts
echo "Filtering miniprot alignments..."
parallel -j "$threads" "grep mRNA output/miniProt_alignments/unfiltered/{/.}.gff | grep -v -e 'StopCodon=1;' -e 'Frameshift=1;' | sort | uniq > output/miniProt_alignments/filtered_mRNA_stop_frameshift/{/.}.gff" :::: lists/assembly_list.txt

# Retain alignments starting with ATG
# First pull out sequence for aligned regions
echo "Retaining alignments starting with ATG..."
parallel -j 2 "bedtools getfasta -fi {} -bed output/miniProt_alignments/filtered_mRNA_stop_frameshift/{/.}.gff -s > output/miniProt_alignments/filtered_mRNA_stop_frameshift/{/.}.fa" :::: lists/assembly_list.txt

# Then filter by ATG
parallel -j $threads "bash src/FilterGFF_byATG.sh {/.} {} output/miniProt_alignments/filtered_mRNA_stop_frameshift output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG 4" :::: lists/assembly_list.txt

echo "Extracting upstream sequence..."
# 5kb upstream (no length filter)
echo "Extracting upstream sequence..."
# 5kb upstream (no length filter)
parallel -j $threads "bash src/ExtractUpstreamCoords.sh {/.} 5000 {//} output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_5kbUpstream F F" :::: lists/assembly_list.txt

# Extract sequence
parallel -j 2 "bedtools getfasta -fi {} -bed output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_5kbUpstream/{/.}.bed > output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_5kbUpstream/{/.}.fa" :::: lists/assembly_list.txt

sed -i 's/::/_/g' output/miniProt_alignments/filtered_mRNA_stop_frameshift_ATG_5kbUpstream/*.fa # Modify seq names to enable parsing by motif scanner
