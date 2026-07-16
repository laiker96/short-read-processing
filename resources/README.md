# Adapter resource

`adapters.fa` is copied unchanged from
`BBMap_39.96.tar.gz:bbmap/resources/adapters.fa` in the parent workspace.

- BBMap archive SHA-256: `e173bdd0d3ca047f378c71dad568a148596c1690bf36abca93e918569c9fb382`
- Extracted adapter FASTA SHA-256: `85abe9d3e40dc37c968f7e4c1227e05976a4ed0583d1dd442d375aa7516f13a9`
- Records: 158

The complete file contains TruSeq, Nextera, and other sequences. The workflow
should select a small assay/library-specific subset before passing adapters to
Cutadapt; it should not blindly treat all 158 records as equivalent 3-prime
adapters.

The reviewed `adapters/nextera.fa` and `adapters/truseq.fa` subsets contain the
specific read-through sequences used by the workflow presets; the full BBMap
collection remains available for provenance and custom review.

## Sample-sheet examples

`sample_sheet.example.tsv` demonstrates the canonical accession-first input,
including default paired-end ATAC/HMMRATAC, an ATAC MACS3 shift/extsize
override, and matched TF ChIP treatment/control rows. The machine-readable
contract is `../schemas/sample-sheet.schema.yaml`.

`random_two.sample_sheet.tsv` is the validated pipeline input for the two
downloaded smoke-test accessions. `random_download_selection.tsv` preserves the
original seeded row selection and source metadata.

`abc_map_dm6_stage5_atac.sample_sheet.tsv` is a coherent bulk ATAC real-data
set derived from `../ABC-map-dm6/metadata/atac_seq_metadata.tsv`. It keeps the
two technical runs for stage-5 replicate 1 under one `sample_id` and the
recommended biological replicate as replicate 2.

`abc_map_dm6_stage5_h3k27ac.sample_sheet.tsv` uses the curated
`alternative_GSE140539_stage5_with_inputs` option from
`../ABC-map-dm6/metadata/h3k27ac_replicates_inputs_metadata.tsv`: two H3K27ac
IP replicates with their replicate-specific input controls.
