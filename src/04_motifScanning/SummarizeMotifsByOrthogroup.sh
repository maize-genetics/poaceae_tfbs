#!/bin/bash
# Processes motif bed files, subsetting to include only motifs within a desired upstream distance
assemblyPath=$1
motifDir=$2
unfilt_interval_dir=$3
filt_interval_dir=$4
distance_upstream=$5
genomeFileDir=$6

assemblyID=$(basename "$assemblyPath" .fa)
echo $assemblyID

temp_dir=$(mktemp -d)
# Extract upstream coordinates
bash src/04_motifScanning/ExtractUpstreamCoords.sh $assemblyID $distance_upstream $genomeFileDir $unfilt_interval_dir $temp_dir T T

# Filter out intervals containing more than 5% N
bash src/04_motifScanning/FilterRegionsContainingNs.sh $temp_dir/$assemblyID.bed $assemblyPath $filt_interval_dir/$assemblyID.bed 5 $assemblyID

# Intersect motifs with upstream intervals, retaining only motifs completely contained within intervals
bedtools intersect -a "$motifDir/$assemblyID.bed" -b "$filt_interval_dir/$assemblyID.bed" -wa -wb -f 1 | \
# Clean up file
awk -v assemblyID="$assemblyID" '
{
    # Split the 9th field by semicolons
    split($9, fields, ";");

    # Extract the OG identifier
    target = "";
    for (i in fields) {
        if (index(fields[i], "Target=") == 1) {
            # Remove "Target=" prefix
            target = substr(fields[i], 8);
            break;
        }
    }

    # Split the target field by underscores
    split(target, target_fields, "_");

    # Output the desired fields
    print $1, $2, $3, $5, target_fields[1], assemblyID
}' | \
sort | \
uniq | \
# Get motif counts for each orthogroup 
awk '{count[$1" "$4" "$5]++} END {for (i in count) print i, count[i], $6}'

# Clean up temp file
rm -r $temp_dir
