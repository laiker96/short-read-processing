#!/usr/bin/env python3
"""Regenerate the deterministic paired-end workflow runtime fixture."""

import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent / "fixtures" / "workflow"
SEED = 1729
READ_LENGTH = 50


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def wrapped_fasta(name: str, sequence: str) -> str:
    lines = [f">{name}"]
    lines.extend(sequence[index : index + 80] for index in range(0, len(sequence), 80))
    return "\n".join(lines) + "\n"


def main() -> None:
    rng = random.Random(SEED)
    chromosome = "".join(rng.choices("ACGT", k=20_000))
    mitochondrion = "".join(rng.choices("ACGT", k=1_000))
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT / "tiny.fa").write_text(
        wrapped_fasta("chr1", chromosome) + wrapped_fasta("chrM", mitochondrion),
        encoding="utf-8",
    )
    (ROOT / "tiny.blacklist.bed").write_text("chr1\t19500\t19600\n", encoding="utf-8")
    (ROOT / "tiny.tss.bed").write_text("chr1\t10000\t10001\tgene1\t0\t+\n", encoding="utf-8")
    (ROOT / "tiny.autosomes.txt").write_text("chr1\n", encoding="utf-8")

    fragments = [(9_500 + index, 150 + index % 41) for index in range(400)]
    fragments.extend((200 + index * 75, 140 + index % 31) for index in range(100))
    r1_records = []
    r2_records = []
    for index, (start, length) in enumerate(fragments, start=1):
        end = start + length
        r1 = chromosome[start : start + READ_LENGTH]
        r2 = reverse_complement(chromosome[end - READ_LENGTH : end])
        quality = "I" * READ_LENGTH
        r1_records.append(f"@fragment{index}/1\n{r1}\n+\n{quality}\n")
        r2_records.append(f"@fragment{index}/2\n{r2}\n+\n{quality}\n")
    (ROOT / "sample_R1.fastq").write_text("".join(r1_records), encoding="utf-8")
    (ROOT / "sample_R2.fastq").write_text("".join(r2_records), encoding="utf-8")


if __name__ == "__main__":
    main()
