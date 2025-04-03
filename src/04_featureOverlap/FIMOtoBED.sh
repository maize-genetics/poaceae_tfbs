# Converts FIMO .tsv format to BED
# Charlie Hale, 2024.06.07
# Usage: bash FIMOtoBED.sh <TSVpath>

raw_motifs=$1
awk '
BEGIN { OFS = "\t" }
NR > 1 && $0 !~ /^#/ && $0 != "" {
  seqName = $3;
    start = $4;
      stop = $5 + 1;
        motifID = $1;
	  print seqName, start, stop, motifID;
  }' "$raw_motifs" | sort -k4,4 | cat
