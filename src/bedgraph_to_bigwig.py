#!/usr/bin/env python3
"""Convert a sorted MACS3 bedGraph to BigWig without loading it into memory."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyBigWig


def read_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected chromosome and size")
            sizes[fields[0]] = int(fields[1])
    if not sizes:
        raise ValueError(f"Chromosome sizes file is empty: {path}")
    return sizes


def convert(source: Path, destination: Path, sizes_path: Path) -> None:
    sizes = read_sizes(sizes_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    temporary.unlink(missing_ok=True)
    writer = pyBigWig.open(str(temporary), "w")
    try:
        writer.addHeader(sorted(sizes.items()))
        chroms: list[str] = []
        starts: list[int] = []
        ends: list[int] = []
        values: list[float] = []
        previous_chrom = ""
        previous_start = -1

        def flush() -> None:
            if chroms:
                writer.addEntries(chroms, starts, ends=ends, values=values)
                chroms.clear()
                starts.clear()
                ends.clear()
                values.clear()

        with source.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith(("track", "browser", "#")):
                    continue
                chrom, start_text, end_text, value_text = line.rstrip("\n").split("\t")[:4]
                start, end = int(start_text), int(end_text)
                if chrom not in sizes:
                    raise ValueError(f"{source}:{line_number}: unknown chromosome {chrom}")
                if end <= start or end > sizes[chrom]:
                    raise ValueError(f"{source}:{line_number}: invalid interval")
                if chrom < previous_chrom or (
                    chrom == previous_chrom and start < previous_start
                ):
                    raise ValueError(f"{source}:{line_number}: input is not sorted")
                previous_chrom, previous_start = chrom, start
                chroms.append(chrom)
                starts.append(start)
                ends.append(end)
                values.append(float(value_text))
                if len(chroms) >= 100_000:
                    flush()
        flush()
    except Exception:
        writer.close()
        temporary.unlink(missing_ok=True)
        raise
    else:
        writer.close()
        temporary.replace(destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bedgraph", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    convert(args.bedgraph, args.output, args.chrom_sizes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
