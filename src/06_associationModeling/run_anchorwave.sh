#!/bin/bash
# Written by Allen Gelfond and Charlie Hale (chale295@gmail.com)
# Runs pairwise anchorwave alignments and outputs maf and gvcf files

# Usage: run_anchorwave.sh --reffa <Reference fasta> --refgff <ref gff> --altfa <Alt fasta> --R <R_value> --Q <Q_value> --threads <threads> [--ref <ref short name>] [--alt <alt short name>]

# Function to extract the basename without extension
get_basename() {
    local full_name="$1"
    local name_without_extension
    name_without_extension="$(basename "$full_name" .fa)"
    echo "$name_without_extension"
}

# Parse named arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --reffa)
            reffa="$2"
            shift 2;;
        --refgff)
            refgff="$2"
            shift 2;;
        --altfa)
            altfa="$2"
            shift 2;;
        --R)
            R="$2"
            shift 2;;
        --Q)
            Q="$2"
            shift 2;;
        --threads)
            threads="$2"
            shift 2;;
        --ref)
            ref="$2"
            shift 2;;
        --alt)
            alt="$2"
            shift 2;;
        *)
            echo "Unknown argument: $1"
            exit 1;;
    esac
done

# Check for mandatory arguments
if [ -z "$reffa" ] || [ -z "$refgff" ] || [ -z "$altfa" ] || [ -z "$R" ] || [ -z "$Q" ] || [ -z "$threads" ]; then
    echo "Usage: $0 --reffa <Reference fasta> --refgff <ref gff> --altfa <Alt fasta> --R <R_value> --Q <Q_value> --threads <threads> [--ref <ref short name>] [--alt <alt short name>]"
    exit 1
fi

# Set default values for optional arguments if not provided
if [ -z "$ref" ]; then
    ref=$(get_basename "$reffa")
fi

if [ -z "$alt" ]; then
    alt=$(get_basename "$altfa")
fi

# Create output directories
mkdir -p output/candidate_visualization/tmp

## Set up gene anchors
echo "Setting up gene anchors for $ref using $reffa and $refgff"
anchorwave gff2seq -r "$reffa" -i "$refgff" -o output/candidate_visualization/tmp/"$ref".CDS.fa
minimap2 \
    -x splice \
    -t "$threads" \
    -k 12 \
    -a \
    -p 0.4 \
    -N 21 \
    "$reffa" output/candidate_visualization/tmp/"$ref".CDS.fa > output/candidate_visualization/tmp/"$ref".cds.sam

echo "Mapping anchors in $alt using $altfa..." 
if [ ! -f output/candidate_visualization/tmp/"${alt}".cds.sam ]; then
    minimap2 \
        -x splice \
        -t "$threads" \
        -k 12 \
        -a \
        -p 0.4 \
        -N 20 \
        "$altfa" output/candidate_visualization/tmp/"$ref".CDS.fa > output/candidate_visualization/tmp/"${alt}".cds.sam
fi

# Run anchorwave alignment
echo "Running whole-genome alignment for $ref vs $alt..."
if [ ! -f output/candidate_visualization/"$ref"-"$alt".maf ]; then
    /usr/bin/time -v anchorwave proali \
        -t 8 \
        -i "$refgff" \
        -as output/candidate_visualization/tmp/"$ref".CDS.fa \
        -r "$reffa" \
        -a output/candidate_visualization/tmp/"${alt}".cds.sam \
        -ar output/candidate_visualization/tmp/"$ref".cds.sam \
        -s "$altfa" \
        -n "$ref"-"$alt" \
        -R "$R" \
        -Q "$Q" \
        -o output/candidate_visualization/"$ref"-"$alt".maf \
        -f output/candidate_visualization/"$ref"-"$alt".f.maf
fi