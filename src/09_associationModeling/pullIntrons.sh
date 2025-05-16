#!/bin/bash
# Script to extract intron sequences from miniProt GFFs. Output is in bed format.
# Original: Charlie Hale, 5/9/22. Modified 2025-02-19 to fix off-by-one error.

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.gff"
    exit 1
fi

input_gff="$1"

# Convert GFF to GTF and filter for CDS entries, extracting transcript_id.
gffread -T -o- "$input_gff" | \
awk '$3 == "CDS" {
    if (match($9, /transcript_id "([^"]+)"/, a)) {
        print $1, $4, $5, $7, a[1];
    }
}' | \
sort -k5,5 -k1,1 -k2,2n | \
awk 'BEGIN {OFS="\t"}
{
    if (prev_transcript == $5 && prev_chrom == $1 && prev_strand == $4) {
        if (prev_end + 1 < $2) {
            # Intron in GTF: [prev_end+1, $2-1] (1-based, inclusive)
            # Convert to BED: subtract 1 from the start to get 0-based coordinates.
            bed_start = prev_end;
            bed_end = $2 - 1;
            print $1, bed_start, bed_end, prev_transcript, ".", $4;
        }
    }
    prev_transcript = $5;
    prev_chrom = $1;
    prev_strand = $4;
    prev_end = $3;
}'