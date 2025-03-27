# Usage: bash script.sh <rawReadDir> <outputFile>
# Define the parent directory and output file
parent_dir=$1
output_file=$2

# Initialize the output file with headers
echo -e "Sample\tPath1\tPath2" > "$output_file"

# Loop through all direct subdirectories of the parent directory
find "$parent_dir" -mindepth 1 -maxdepth 1 -type d | while read dir; do
  # Extract the final directory name
  final_dir_name=$(basename "$dir")
  
# Generate the list of files, comma-separated
  path1=$(cd "$dir" && realpath *_1.fastq.gz | paste -sd ",")
  path2=$(cd "$dir" && realpath *_2.fastq.gz | paste -sd ",")
  # Check if file_list is not empty
  if [ -n "$path1" ] && [ -n "$path2" ]; then
	echo -e "${final_dir_name}\t${path1}\t${path2}" >> "$output_file"
  fi
done

