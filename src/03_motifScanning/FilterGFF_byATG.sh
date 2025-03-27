# Given a GFF file, filter to retain intervals starting with ATG
# Charlie Hale, 2024.03.12
speciesID=$1
fasta=$2
unfilteredDir=$3
filteredDir=$4
atg_window=$5

echo $speciesID

# First extract sequence for unfiltered intervals
bedtools getfasta -fi $fasta -bed "${unfilteredDir}/${speciesID}.gff" -s > "${unfilteredDir}/${speciesID}.fa"


# Filter out sequences that don't start with ATG
awk -v atg_window="$atg_window" 'BEGIN {RS=">"; ORS=""} $0!="" {sequence = $0; sub("\n", "@", sequence); sub("\n", "", sequence); sub("@", "\n", sequence); if (substr(sequence, index(sequence, "\n"), atg_window) ~ /ATG/) {print ">" sequence "\n"}}' "${unfilteredDir}/${speciesID}.fa" > "${filteredDir}/${speciesID}.fa"

awk '/^>/{split(substr($0,2), a, "[:\\-()]"); print a[1]"\t"a[2]"\t"a[3]}' "${filteredDir}/${speciesID}.fa" > "${filteredDir}/${speciesID}.gff.tmp"

# Filter out the corresponding intervals in the GFF file
awk -v OFS="\t" 'NR==FNR {key=$1"\t"$2+1"\t"$3; arr[key]; next}
	{
		key=$1"\t"$4"\t"$5;
		if (key in arr) print $0
	}' "${filteredDir}/${speciesID}.gff.tmp" "${unfilteredDir}/${speciesID}.gff" > "${filteredDir}/${speciesID}.gff"

# Remove temp file
rm "${filteredDir}/${speciesID}.gff.tmp"
