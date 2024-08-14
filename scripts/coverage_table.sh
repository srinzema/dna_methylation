#!/bin/bash

# Check if enough arguments were provided
if [ $# -lt 2 ]; then
  echo "Usage: $0 <input_file> <output_file>"
  exit 1
fi

# Assign the first and second arguments to variables
input_file=$1
output_file=$2

# Process the input file with awk and write to the output file
awk -v outfile="$output_file" '
BEGIN { 
  # Open the output file for writing
  print "chromosome\tmethylated\tunmethylated\tpercentage" > outfile;
}

{ 
  # Accumulate totals for all chromosomes
  total_methylated += $4;
  total_unmethylated += $5;

  # Collect sums for each unique chromosome
  methylated_by_chr[$1] += $4;
  unmethylated_by_chr[$1] += $5;
}

END {
  total = total_methylated + total_unmethylated;
  
  # Print totals for all chromosomes
  printf "Total\t%d\t%d\t%.2f%%\n", total_methylated, total_unmethylated, (total_methylated / total) * 100 >> outfile;

  # Print results per chromosome in alphabetical order
  PROCINFO["sorted_in"] = "@ind_str_asc";
  for (chr in methylated_by_chr) {
    total_by_chr = methylated_by_chr[chr] + unmethylated_by_chr[chr];
    if (total_by_chr > 0) {
      printf "%s\t%d\t%d\t%.2f%%\n", chr, methylated_by_chr[chr], unmethylated_by_chr[chr], (methylated_by_chr[chr] / total_by_chr) * 100 >> outfile;
    } else {
      printf "%s\t%d\t%d\t0.0\n", chr, methylated_by_chr[chr], unmethylated_by_chr[chr] >> outfile;
    }
  }
}' "$input_file"