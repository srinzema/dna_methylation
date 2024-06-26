#!/bin/bash

# Check if an argument was provided
if [ -z "$1" ]; then
  echo "No argument provided. Please provide an input file."
  exit 1
fi

# Assign the first argument to a variable
input_file=$1

# Process the input file with awk
awk '
BEGIN { 
  # Initialize totals and arrays
  total_methylated = 0;
  total_unmethylated = 0;
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
  # Print TSV header
  print "chromosome\tmethylated\tunmethylated\tpercentage";

  # Print total
  total = total_methylated + total_unmethylated;
  printf "Total\t%d\t%d\t%.2f\n", total_methylated, total_unmethylated, (total_methylated / total) * 100;

  PROCINFO["sorted_in"] = "@ind_str_asc";
  for (chr in methylated_by_chr) {
    total_by_chr = methylated_by_chr[chr] + unmethylated_by_chr[chr];
    if (total_by_chr > 0) {
      printf "%s\t%d\t%d\t%.2f\n", chr, methylated_by_chr[chr], unmethylated_by_chr[chr], (methylated_by_chr[chr] / total_by_chr) * 100;
    } else {
      printf "%s\t%d\t%d\t0.0\n", chr, methylated_by_chr[chr], unmethylated_by_chr[chr];
    }
  }
}' "$input_file"
