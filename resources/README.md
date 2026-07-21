# Pipeline resources

This directory contains only versioned inputs needed to reproduce the dm6 atlas
runs and the adapter sequences used by the workflow.

## Atlas inputs

- `atlas_atac_seq_metadata.tsv`: curated ATAC source metadata and selection
  evidence.
- `atlas_atac_selected.sample_sheet.tsv`: canonical pipeline input containing
  the 23 selected runs representing 22 biological ATAC libraries.
- `atlas_atac_conditions.tsv`: explicit assignment of those 22 biological
  libraries to nine tissue/stage conditions for the optional atlas stage.
- `atlas_h3k27ac_metadata_ip_only.tsv`: curated H3K27ac source metadata and
  selection evidence.
- `atlas_h3k27ac_ip_only.sample_sheet.tsv`: canonical pipeline input containing
  the 15 selected H3K27ac IP-only runs. No external input-DNA libraries are
  included.

The two sample sheets are deterministic derivatives of the metadata tables:

```bash
python src/write_atlas_atac_sample_sheet.py \
  resources/atlas_atac_seq_metadata.tsv \
  --output resources/atlas_atac_selected.sample_sheet.tsv

python src/write_atlas_h3k27ac_sample_sheet.py \
  resources/atlas_h3k27ac_metadata_ip_only.tsv \
  --output resources/atlas_h3k27ac_ip_only.sample_sheet.tsv
```

## Adapter sequences

`adapters.fa` is copied unchanged from
`BBMap_39.96.tar.gz:bbmap/resources/adapters.fa` in the parent workspace.

- BBMap archive SHA-256: `e173bdd0d3ca047f378c71dad568a148596c1690bf36abca93e918569c9fb382`
- Extracted adapter FASTA SHA-256: `85abe9d3e40dc37c968f7e4c1227e05976a4ed0583d1dd442d375aa7516f13a9`
- Records: 158

The complete file contains TruSeq, Nextera, and other sequences. It must not be
passed wholesale to Cutadapt as though all 158 entries were equivalent 3-prime
adapters.

The reviewed `adapters/nextera.fa` and `adapters/truseq.fa` subsets contain the
read-through sequences associated with the named workflow presets. The full
BBMap collection is retained for provenance and custom adapter review.
