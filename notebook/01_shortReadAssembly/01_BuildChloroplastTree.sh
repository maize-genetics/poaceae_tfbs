# Builds a phylogeny across Poaceae using matK chloroplast gene
# Charlie Hale, 2024.04.22

# Index assemblies for minimap2
mkdir -p output/chloroplastTree/intermed_alignments
outDir="output/chloroplastTree"

# Use Poaceae matK and rbcL sequences from BOLD database to compare against (after downloading)
# First query for Poaceae and download fasta file
grep -A 1 "matK" output/chloroplastTree/poaceae_bold.fas > output/chloroplastTree/bold_matk.fas
grep -A 1 "rbcL" output/chloroplastTree/poaceae_bold.fas > output/chloroplastTree/bold_rbcL.fas
# Pull out rbcL from representative taxa to query against short reads
grep -A 1 -Ff lists/rep_species_rbcL.txt output/chloroplastTree/bold_rbcL.fas > output/chloroplastTree/rbcL_queries.fa 
# Align Streptochaeta matK vs short reads
 parallel -j 8 "minimap2 -ax asm20 --eqx -I 100g -t 5 --secondary=no data/genomes/shortread/{} ${outDir}/streptochaeta_angustifolia_matK.fasta > ${outDir}/intermed_alignments/{}_matK.sam 2>> ${outDir}/error_log_miniProtMatK.txt && \
                 samtools view -b ${outDir}/intermed_alignments/{}_matK.sam > ${outDir}/intermed_alignments/{}_matK.bam 2>> ${outDir}/error_log_samtoolsView.txt && \
                 bedtools bamtobed -i ${outDir}/intermed_alignments/{}_matK.bam > ${outDir}/intermed_alignments/{}_matK.bed 2>> ${outDir}/error_log_bamtobed.txt" :::: lists/shortreads.txt
# Align rbcLs vs short reads and keep top-scoring alignment
parallel -j 8 "minimap2 -ax asm20 --eqx -I 100g -t 5 --secondary=no data/genomes/shortread/{}.fa ${outDir}/rbcL_queries.fa > ${outDir}/intermed_alignments/{}_rbcL.sam 2>> ${outDir}/error_log_miniProtRbcL.txt && \
samtools view -b ${outDir}/intermed_alignments/{}_rbcL.sam | samtools sort -o ${outDir}/intermed_alignments/{}_rbcL.bam 
" :::: lists/shortreads.txt

# Extract matK sequence from assemblies
parallel -j 4 "bedtools getfasta -fi data/genomes/shortread/{} -bed ${outDir}/intermed_alignments/{}_matK.bed -fo ${outDir}/intermed_alignments/{}_matK.fa -s && \
awk '/^>/{print \">\" \"{}_\" substr(\$0,2); next} 1' ${outDir}/intermed_alignments/{}_matK.fa > ${outDir}/intermed_alignments/{}_matK_modified.fa && \
mv ${outDir}/intermed_alignments/{}_matK_modified.fa ${outDir}/intermed_alignments/{}_matK.fa" :::: lists/shortreads.txt 

# Merge matK sequences into a single file for mafft
cat ${outDir}/intermed_alignments/*_matK.fa output/chloroplastTree/bold_matk.fas > ${outDir}/matK_combined_bold.fa

# Run multiple alignment of matK
mafft --auto --thread 39 ${outDir}/matK_combined_bold.fa > ${outDir}/matK_MSA_bold.fa
sed -i 's/[:()]/_/g' ${outDir}/matK_MSA_bold.fa
# Build tree with RaxML
cd ${outDir}
raxmlHPC-PTHREADS -s matK_MSA_bold.fa -n matk_bold_tree -m GTRGAMMA -p 12345 -T 38