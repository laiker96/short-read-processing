import gzip
import runpy
from pathlib import Path
from types import SimpleNamespace


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_snakemake_scripts_do_not_use_future_imports():
    for script in (REPO_ROOT / "workflow/scripts").glob("*.py"):
        assert "from __future__ import" not in script.read_text(), script


def test_extract_tss_handles_both_strands_and_filters_unknown_contigs(tmp_path):
    annotation = tmp_path / "genes.gtf.gz"
    with gzip.open(annotation, "wt") as handle:
        handle.write(
            'chr2L\tRefSeq\ttranscript\t10\t20\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr2L\tRefSeq\ttranscript\t30\t50\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n'
            'chrUnknown\tRefSeq\ttranscript\t1\t5\t.\t+\t.\tgene_id "g3";\n'
        )
    fai = tmp_path / "genome.fa.fai"
    fai.write_text("chr2L\t100\t0\t0\t0\n")
    output = tmp_path / "tss.bed"
    log = tmp_path / "tss.log"
    snakemake = SimpleNamespace(
        input=SimpleNamespace(annotation=str(annotation), fai=str(fai)),
        output=SimpleNamespace(bed=str(output)),
        log=[str(log)],
    )

    runpy.run_path(
        REPO_ROOT / "workflow/scripts/extract_tss.py",
        init_globals={"snakemake": snakemake},
    )

    assert output.read_text().splitlines() == [
        "chr2L\t9\t10\tt1\t0\t+",
        "chr2L\t49\t50\tt2\t0\t-",
    ]


def test_write_autosomes_checks_fasta_contigs(tmp_path):
    fai = tmp_path / "genome.fa.fai"
    fai.write_text("chr2L\t100\t0\t0\t0\nchr2R\t100\t0\t0\t0\n")
    output = tmp_path / "autosomes.txt"
    log = tmp_path / "autosomes.log"
    snakemake = SimpleNamespace(
        input=SimpleNamespace(fai=str(fai)),
        output=SimpleNamespace(contigs=str(output)),
        params=SimpleNamespace(autosomes=["chr2L", "chr2R"]),
        log=[str(log)],
    )

    runpy.run_path(
        REPO_ROOT / "workflow/scripts/write_autosomes.py",
        init_globals={"snakemake": snakemake},
    )

    assert output.read_text() == "chr2L\nchr2R\n"
