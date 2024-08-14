#!/usr/bin/env python

import argparse
import shutil
from pathlib import Path

parser = argparse.ArgumentParser("Trackhub generator")

parser.add_argument(
    "--output-dir",
    type=Path,
    required=True,
    help="Directory where the trackhub is made. Default is the current working directory."
)

parser.add_argument(
    "--assembly",
    type=str,
    required=True,
    help="Assembly name used in the trackhub. This argument is required."
)

parser.add_argument(
    "--email",
    type=str,
    required=True,
    help="Your email. This argument is required."
)

parser.add_argument(
    "--trackfiles",
    type=Path,
    nargs="+",
    required=True,
    help="List of files to host. This argument is required."
)

# Parse arguments
args = parser.parse_args()

# Make output dir
print(args.output_dir)
output_dir = Path(args.output_dir)
output_dir.mkdir(exist_ok=True)

# Make hub.txt
email = args.email
with open(output_dir / "hub.txt", "w") as hub_file:
    print(hub_file.name)
    hub_file.write("hub custom_trackhub\n")
    hub_file.write("shortLabel Custom Trackhub\n")
    hub_file.write("longLabel Custom Trackhub\n")
    hub_file.write("genomesFile genomes.txt\n")
    hub_file.write(f"email {email}\n")

# Make tracks
tracks = args.trackfiles
track_directory = output_dir / "tracks"
track_directory.mkdir(exist_ok=True)
with open(track_directory / "trackDb.txt", "w") as track_file:
    print(track_file.name)
    for input_track in tracks:
        track = track_directory / input_track.name
        shutil.copy(input_track, track)

        track_name = track.name.replace('.deduplicated.bedGraph.gz', '')
        track_file.write(f"track {track_name}\n")
        track_file.write("type bedGraph\n")
        track_file.write(f"shortLabel {track_name.replace('_', '')}\n")
        track_file.write(f"longLabel {track_name}\n")
        track_file.write("visibility full\n")
        track_file.write("color 0,0,255")


# Make genomes.txt
assembly = args.assembly
with open(output_dir / "genomes.txt", "w") as genome_file:
    print(genome_file.name)
    genome_file.write(f"genome {assembly}\n")
    genome_file.write(f"trackDb tracks/trackDb.txt\n")
  