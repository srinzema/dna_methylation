#!/bin/bash

# Check if enough arguments were provided
if [ $# -lt 3 ]; then
  echo "Usage: $0 <summary_output_file> <input_file1> <input_file2> ... <input_fileN>"
  exit 1
fi

# Assign the first argument to the output file
output_file=$1
shift

for file in "$@"; do
    sample_name="$(basename $file .coverage.summary.txt)"
    echo $sample_name
    awk -v sample="$sample_name" 'NR==1 {next} {print $1, $4}' "$file" | 
    awk -v fname="$file" '{print $1, $2 "\t" fname}' fname="$sample_name" > "${output_file}"
done