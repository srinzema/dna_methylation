#!/usr/bin/env python

import argparse
import os
import gzip
import shutil

from pathlib import Path

parser = argparse.ArgumentParser("Trackhub generator")

parser.add_argument(
    "--output-dir",
    type=Path,
    required=True,
    help="Directory where the trackhub is made. Default is the current working directory.",
)

parser.add_argument(
    "--assembly",
    type=str,
    required=True,
    help="Assembly name used in the trackhub. This argument is required.",
)

parser.add_argument(
    "--email", type=str, required=True, help="Your email. This argument is required."
)

parser.add_argument(
    "--trackfiles",
    type=Path,
    nargs="+",
    required=True,
    help="List of files to host. This argument is required.",
)

# Parse arguments
args = parser.parse_args()

trackhub_root: Path = args.output_dir
trackhub_root.mkdir(parents=True, exist_ok=True)

assembly_dir: Path = trackhub_root / args.assembly
assembly_dir.mkdir(parents=True, exist_ok=True)


hub: Path = trackhub_root / "hub.txt"
trackdb: Path = assembly_dir / "trackDb.txt"
genomes: Path = trackhub_root / "genomes.txt"

# Make hub.txt
print(f"Making hub.txt at: {hub}")
with open(hub, "w") as file:
    file.write("hub MethylationHub\n")
    file.write("shortLabel Methylation hub\n")
    file.write("longLabel Methylation hub\n")
    file.write("genomesFile genomes.txt\n")
    file.write(f"email {args.email}\n")


# make genomes.txt
print(f"Making genomes.txt at: {genomes}")
with open(genomes, "w") as file:
    file.write(f"genome {args.assembly}\n")
    file.write(f"trackDb {args.assembly}/trackDb.txt\n")

# make trackDb.txt
print(f"Making trackDb.txt at: {trackdb}")
with open(trackdb, "w") as file:
    for trackfile in args.trackfiles:
        track_name = trackfile.name.replace("_", "").split(".")[0]
        print(f" - Adding {trackfile.name} as {track_name}")

        file.write(f"track {track_name}\n")
        file.write("type bedGraph\n")
        file.write(f"bigDataUrl {trackfile.name.replace('.gz', '')}\n")
        file.write(f"shortLabel {track_name}\n")
        file.write(f"longLabel {track_name}\n")
        file.write("visibility full\n")
        file.write("\n")

        print(f" - Unpacking {trackfile.name} to {assembly_dir}")
        new_location = assembly_dir / trackfile.name.replace(".gz", "")
        with gzip.open(trackfile, "rb") as f_in:
            with open(new_location, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
