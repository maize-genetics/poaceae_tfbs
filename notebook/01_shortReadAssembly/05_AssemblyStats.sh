# Calculates assembly statistics
# Charlie Hale, 2025.08.14

threads=110
assemblyDir="data/genomes/"
# Unzip files if necessary
pigz -d -p $threads ${assemblyDir}/shortread*.gz

# Calculate and extract assembly stats using assembly-stats 
mkdir -p results
outFile="results/assembly_stats.csv"

echo "assemblyID,assemblySize,nContig,avgContigSize,largestContigSize,N50,N90" > "$outFile"

find $assemblyDir -maxdepth 2 -type f \( -name "*.fa" \) -print0 |
parallel -0 -k -j "2" '
  asm_id="{/.}"
  assembly-stats {} |
  awk -v id="$asm_id" "BEGIN{FS=\"[ =,]+\"}
       /^sum =/ {sum=\$2; n=\$4; ave=\$6; largest=\$8}
       /^N50 =/ {N50=\$2}
       /^N90 =/ {N90=\$2}
       END      {printf \"%s,%s,%s,%s,%s,%s,%s\n\", id, sum, n, ave, largest, N50, N90}"
' >> "$outFile"

