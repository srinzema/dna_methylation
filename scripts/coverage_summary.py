#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path


def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Process some files.")

    # Add arguments
    parser.add_argument('output_file', type=Path, help="Output file path")
    parser.add_argument('input_files', type=Path, nargs='+', help="List of input files")

    # Parse the arguments
    args = parser.parse_args()

    # Access the output file and input files
    output_file = args.output_file
    input_files = args.input_files

    # Print out the results (or replace with actual processing)
    print(f"Output file: {output_file}")
    data_frames = [get_percentage(file) for file in input_files]

    merged_df = pd.concat(data_frames, axis=1)
    merged_df.to_csv(output_file, sep="\t")
        
    

def get_percentage(tsv_file: Path):
    sample_name = tsv_file.stem.replace(".coverage.summary", "")
    print(f"Parsing: {sample_name} \t {tsv_file}")

    # Read the TSV file into a pandas DataFrame with the first column as the index
    df = pd.read_csv(tsv_file, sep='\t', index_col=0)
    
    # Keep only the column named "percentage"
    df = df[['percentage']]

    # Remove any non-numeric characters (e.g., %) and convert the column to float
    df['percentage'] = df['percentage'].replace(r'[^\d.]+', '', regex=True).astype(float)

    df.columns = [sample_name]
    return df


if __name__ == "__main__":
    main()
