#!/usr/bin/env python3
import sys
import math
import argparse
import os
from collections import namedtuple
from typing import List, Iterator
import pybigtools
# import numpy as np

# Define the named tuple.
# The fourth field is annotated as Union[int, str] to allow either type.
BedRecord = namedtuple("BedRecord", ["chrom", "start", "end", "value"])


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def iter_bed_file(filepath: str) -> Iterator[BedRecord]:
    """Load a BED-like file into memory as a list of BedRecord objects."""

    with open(filepath, "r") as f:
        try_convert = 10
        for line in f:
            # Skip empty lines or lines that are comments.
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"Line has insufficient columns: {line}")

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])

            # Handle the fourth field if it exists; try converting it to int.
            if len(parts) >= 4:
                fourth = parts[3:6]
                if try_convert > 0:
                    try:
                        fourth = float(parts[3])
                    except ValueError:
                        # If conversion fails, leave it as a string.
                        try_convert -= 1
            else:
                fourth = ""

            yield BedRecord(chrom, start, end, fourth)


def data_in_range(bw: pybigtools.BBIRead, chrom, start, end, missing=math.nan):
    if end < start:
        return bw.values(chrom, end, start, missing=missing)[::-1]
    else:
        return bw.values(chrom, start, end, missing=missing)


def compute_matrix(
    bw: pybigtools.BBIRead,
    coords: List[BedRecord],
    before=500,
    after=1000,
    handle=sys.stdout,
):
    names = set()
    for rec in coords:
        name = rec.value[0]
        # see: https://github.com/deeptools/deeptools_intervals/blob/master/deeptoolsintervals/parse.py
        # for ensuring non-duplicate names
        i = 0
        while name in names:
            i += 1
            name = f"{rec.value[0]}_r{i}"
        names.add(name)

        strand = rec.value[2]
        # TODO: might be some fencepost errors here
        if strand == "+":
            start = rec.start - before
            end = rec.start + after
        else:
            start = rec.end + before
            end = rec.end - after

        row_info = [rec.chrom, rec.start, rec.end, name, rec.value[1], strand]
        handle.write("\t".join(map(str, row_info)))
        values = data_in_range(bw, rec.chrom, start, end)
        handle.write("\t")
        # Replacing this with np.char.mod("%.6f", values) causes a huge performance regression
        # np.array2string also sucks
        handle.write("\t".join(f"{x:.6f}" for x in values))
        handle.write("\n")


parser = argparse.ArgumentParser(
    description="Process a BEDGRAPH and REGION BEDFILE with optional parameters."
)
# Two required positional arguments
parser.add_argument("bigwig", help="Input bigwig file")
parser.add_argument("region_bedfile", help="Input REGION BEDFILE")

# Optional integer parameters -a and -b
parser.add_argument(
    "-a",
    "--after",
    type=int,
    default=1500,
    help="Optional integer parameter a (default: 500)",
)
parser.add_argument(
    "-b",
    "--before",
    type=int,
    default=500,
    help="Optional integer parameter b (default: 1000)",
)

# Optional output file; if not provided, defaults to stdout.
parser.add_argument(
    "-o",
    "--output",
    nargs="?",
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="Output filename (defaults to stdout if not provided)",
)


def main():
    args = parser.parse_args()

    os.makedirs("output", exist_ok=True)

    bw = pybigtools.open(args.bigwig)
    bed_data = iter_bed_file(args.region_bedfile)

    eprint("computing matrix")
    compute_matrix(bw, bed_data, args.before, args.after, args.output)


if __name__ == "__main__":
    main()
