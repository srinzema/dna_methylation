#!/bin/bash

# Check if enough arguments were provided
if [ $# -lt 3 ]; then
  echo "Usage: $0 <summary_output_file> <input_file1> <input_file2> ... <input_fileN>"
  exit 1
fi

# Assign the first argument to the summary output file
summary_output_file=$1
shift

# Write header
echo "# id: 'Output from coverage_multiqc.sh'" > "$summary_output_file"
echo "# section_name: 'Coverage summary'" >> "$summary_output_file"
echo "# description: 'Methylated / unmethylated per sample'" >> "$summary_output_file"
echo "# format: 'tsv'" >> "$summary_output_file"
echo "# plot_type: 'bargraph'" >> "$summary_output_file"
echo -e "sample\tmethylated\tunmethylated" >> "$summary_output_file"


for summary in "$@"; do
    sample_name="$(basename $summary .coverage.summary.txt)"
    awk -v sample="$sample_name" '
    FNR == 2 { 
        printf "%s\t%d\t%d\n", sample, $2, $3
    }' "$summary" >> $summary_output_file
done
