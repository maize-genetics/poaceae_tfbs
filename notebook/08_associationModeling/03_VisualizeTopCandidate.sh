# Perform whole genome alignments and lift over motifs to visualize gain/loss at key genes
# Charlie Hale, 2025.03.03
# NOTE: Use the promoterAnnotation.yaml environment to run anchorwave, and the crossMap.yaml environment to run CrossMap.
mkdir -p output/candidate_visualization/

threads=39

# Download rice gff
wget https://ftp.gramene.org/oryza/release-8/gff3/oryza_sativazs97/OsZS97.gff \
   -O output/candidate_visualization/OsZS97.gff
wget https://ftp.gramene.org/oryza/release-8/fasta/oryza_sativazs97/dna/OsZS97.dna.toplevel.fa \
    -O output/candidate_visualization/OsZS97.fa 
# Compare rice and Zizania
refName="Osativa_ZS97"
altName="Zpalustris"
refFasta="output/candidate_visualization/OsZS97.fa"
refGFF="output/candidate_visualization/OsZS97.gff"
altFasta="data/genomes/longread/ASM1927943v1.fa"
altMotifs="output/motifOutput/fimo/collapsed_5kbUpstream/ASM1927943v1.bed"
refMotifs="output/motifOutput/fimo/collapsed_5kbUpstream/ZS97RS3.bed"

# Whole genome alignment to maize
bash src/08_associationModeling/run_anchorwave.sh \
    --reffa $refFasta \
    --refgff $refGFF \
    --altfa $altFasta \
    --R 1 \
    --Q 1 \
    --threads $threads \
    --ref $refName \
    --alt $altName

# Generate bam
python3 lib/maf-convert sam output/candidate_visualization/$refName-$altName.maf \
    | sed 's/[0-9]+H//g' > output/candidate_visualization/$refName-$altName.sam
cat output/candidate_visualization/$refName-$altName.sam | \
samtools view \
    -O BAM \
    --threads $threads \
    --reference $refFasta - | \
samtools sort - > output/candidate_visualization/$refName-$altName.bam
samtools index output/candidate_visualization/$refName-$altName.bam

## Lift over motifs
# Generate chain file for liftover
# NOTE: MAFtoChain is designed to get around an Anchorwave bug with negative stranded alignments that has been fixed since v 1.2.3
g++ -o MAFInvert -fopenmp src/08_associationModeling/MAFInvert.cpp
./MAFInvert -t $threads output/candidate_visualization/$refName-$altName.maf > output/candidate_visualization/$altName-$refName.maf
rm MAFInvert # Remove compiled script

g++ -o MAFtoChain -fopenmp src/08_associationModeling/MAFtoChain.cpp
./MAFtoChain -t $threads output/candidate_visualization/$altName-$refName.maf  > output/candidate_visualization/$altName-$refName.chain
rm MAFtoChain # Remove compiled script

# Lift over motifs
export PYTHONPATH=/programs/CrossMap-0.7.3/lib64/python3.9/site-packages:/programs/CrossMap-0.7.3/lib/python3.9/site-packages
export PATH=/programs/CrossMap-0.7.3/bin:$PATH
python /programs/CrossMap-0.7.3/bin/CrossMap bed \
 output/candidate_visualization/$altName-$refName.chain \
 $altMotifs \
 output/candidate_visualization/$altName-motifs-$refName-Liftover_tmp.bed

# Annotate with colors and motif IDs for visualization
awk 'NR==FNR { colColor[$2]=$1; colCode[$2]=$3; next } \
    { print $1 "\t" $2 "\t" $3 "\t" colCode[$5] "\t0\t+\t" $2 "\t" $3 "\t" colColor[$5] }' \
    data/JASPAR2024_SummarizedClusters.txt \
    output/candidate_visualization/$altName-motifs-$refName-Liftover_tmp.bed > \
    output/candidate_visualization/$altName-motifs-$refName-Liftover.bed

rm  output/candidate_visualization/$altName-motifs-$refName-Liftover_tmp.bed

# Convert ref motifs to proper format for visualization
bash src/06_motifTurnover/Convert_ncbi_format.sh \
    $refMotifs \
    output/candidate_visualization/Os_Zs97_key.txt.tsv \
    > output/candidate_visualization/${refName}_motifs5kb_tmp.bed

awk 'NR==FNR { colColor[$2]=$1; colCode[$2]=$3; next } \
    { print $1 "\t" $2 "\t" $3 "\t" colCode[$4] "\t0\t+\t" $2 "\t" $3 "\t" colColor[$4] }' \
    data/JASPAR2024_SummarizedClusters.txt \
    output/candidate_visualization/${refName}_motifs5kb_tmp.bed > \
    output/candidate_visualization/${refName}_motifs5kb.bed

rm  output/candidate_visualization/${refName}_motifs5kb_tmp.bed



