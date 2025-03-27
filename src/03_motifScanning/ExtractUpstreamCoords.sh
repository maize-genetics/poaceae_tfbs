# Given a genome file and input gff, extracts coordinates upstream of gff intervals, a given number of bp upstream
# Returns a bed file containing  upstream intervals and seq metadata
# Charlie Hale, 2024.03.13

speciesID=$1
distUpstream=$2
genomeFileDir=$3
gffDir=$4
outDir=$5
filter_shorter=$6 # Indicate whether to filter out intervals shorter than designated distUpstream ("T" or "F")
primary_alignment=$7  # Indicate whether to filter out seconday alignments

# Make temp files
tmp_out_1=$(mktemp)
tmp_out_2=$(mktemp)
tmp_out_3=$(mktemp)

# Extract upstream intervals
bedtools flank -i "${gffDir}/${speciesID}.gff"  -g "${genomeFileDir}/${speciesID}.fa.fai" -s -l $distUpstream -r 0 > "$tmp_out_1"

# Re-arrange file to return in bed format, filtering out intervals shorter than upstream distance
if [ "$filter_shorter" = "T" ]; then
  awk -v speciesID="$speciesID" -v distUpstream="$distUpstream" 'BEGIN{ OFS=FS="\t" }; { $9 = $9 "_" speciesID; gsub(" ", "_", $9); t = $9; $9 = $3; $3 = t}; {if (($5 - $4 + 1) >= distUpstream) print $1, $4 - 1, $5, $3, ".", $7}' "$tmp_out_1" > "$tmp_out_2"
else
  awk -v speciesID="$speciesID" -v distUpstream="$distUpstream" 'BEGIN{ OFS=FS="\t" }; { $9 = $9 "_" speciesID; gsub(" ", "_", $9); t = $9; $9 = $3; $3 = t}; {print $1, $4 - 1, $5, $3, ".", $7}' "$tmp_out_1" > "$tmp_out_2"
fi

if [ "$primary_alignment" = "T" ]; then
  grep "Rank=1;" "$tmp_out_2" > "$tmp_out_3"
else 
  cp "$tmp_out_2" "$tmp_out_3"
fi

cat "$tmp_out_3" > "${outDir}/${speciesID}.bed"

rm $tmp_out_1
rm $tmp_out_2
rm $tmp_out_3
