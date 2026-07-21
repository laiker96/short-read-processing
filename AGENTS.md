# AGENTS.md

Repository-specific guidance for coding agents and maintainers.

## Priorities

1. Scientific and computational correctness
2. Reproducibility and idempotency
3. Minimal focused changes
4. Verification
5. Clear user-facing behavior
6. Security and data hygiene

Do not silently change scientific defaults, thresholds, output formats, file
names, or sample relationships.

## Repository contract

- The public data input is an accession CSV/TSV conforming to
  `schemas/sample-sheet.schema.yaml`.
- Users do not provide raw FASTQ paths to the primary entry point.
- `src/run_pipeline.py` is the canonical accession-to-results command.
- `workflow/Snakefile` is the canonical processing workflow. Put shared rules
  under `workflow/rules/`; avoid creating parallel experimental Snakefiles for
  production behavior.
- `dm6` and `hg38` reference preparation belongs in the workflow DAG.
- Raw FASTQs are immutable inputs. Never overwrite or modify them in place.
- Technical runs share a `sample_id`, run lane-level QC/alignment separately,
  and merge before biological-library duplicate marking.
- Biological replicates have distinct `sample_id` values and remain separate
  through this processing pipeline.

## Scientific invariants

- Paired-end ATAC defaults to HMMRATAC for primary peaks.
- HMMRATAC is invalid for single-end data.
- The canonical paired-end ATAC refinement branch retains proper pairs with
  `0 < abs(TLEN) < 150`, creates a Tn5-shifted CPM BigWig, calls lenient MACS3
  candidates with `-f BAM --nomodel --shift -75 --extsize 150 -q 0.10`, and
  refines 50-400 bp intervals at mean CPM >= 2. CPM thresholds are evaluated
  from high to low. Modes separated by at least 25% relative saddle prominence
  retain their last separate qualifying boundaries; shallower subdivisions
  merge and continue expanding as one component.
- CPM-refined peaks are signal-derived and must not be described as q-values or
  independently FDR-controlled calls.
- The optional ATAC atlas uses an explicit condition map. Within conditions,
  pooled refined peaks require configured biological-replicate coverage.
  Across conditions, fixed-width summit windows undergo iterative overlap
  removal; never replace this with raw interval union/`bedtools merge` chaining.
- Cross-condition atlas presence, peak coverage, and CPM are separate matrices;
  a tissue-specific peak does not require support from another tissue.
- Contributor-aware signal shaping is downstream of DHS-seed grouping and must
  not change membership. Use at most one pooled profile per contributing
  condition, normalize profiles locally, and weight contributing conditions
  equally; absent conditions do not contribute zeros. Constrain each shaped
  summit to a contributing source DHS so neighboring elements cannot capture it.
- The fixed-atlas FWHM model uses all condition consensus DHSs but unions them
  within each condition before counting support. FWHM is the connected
  half-maximum condition-support component associated with the fixed anchor;
  do not substitute raw DHS-record counts or pooled CPM signal.
- Keep the narrow-source-first atlas as a comparison branch. It prioritizes
  original DHS width before signal when selecting non-overlapping 250-bp
  anchors and may annotate a broad bridging source to multiple retained anchors;
  do not silently replace the canonical signal-prioritized fixed atlas.
- Keep center-mode half-prominence boundaries separate from ordinary FWHM:
  select the local support mode nearest the anchor center and use its valley
  prominence rather than silently replacing the highest-support result.
- TF ChIP defaults to narrow MACS3 peaks; histone ChIP defaults to broad peaks.
- ChIP controls are explicit sample-sheet relationships. IP-only ChIP is valid
  but must be labeled as such and runs MACS3 without `-c`.
- MACS3 `callpeak` outputs include `-B --SPMR` treatment-pileup and
  control-lambda bedGraphs.
- FRiP denominators must be documented and fragment-aware for paired-end data.

If a task changes any invariant, update the schema, validators, tests, README,
and provenance together and call out the behavior change explicitly.

## Environment and dependencies

- Keep the orchestration environment at repository-root `.venv`.
- Set `MAMBA_ROOT_PREFIX="$PWD/.micromamba"` when practical.
- Keep Snakemake rule environments below `.snakemake/conda` using the local
  profile.
- Add orchestration dependencies to `environment.yml` only when required by
  repository commands or tests.
- Add bioinformatics dependencies to the smallest applicable file under
  `workflow/envs/`.
- R/Bioconductor dependencies belong only in `workflow/envs/chip_qc.yaml` unless
  a new R-dependent rule is explicitly introduced.

Do not create named environments outside the repository.

## Idempotency and file safety

- Every command must be safe to rerun with identical inputs and parameters.
- Prefer temporary/staging outputs followed by atomic replacement.
- Do not replace unchanged manifests or resolved configs merely to update a
  timestamp.
- Preserve aria2 resume state and validate reported checksums.
- Promote SRA-derived FASTQs only after conversion and compression complete.
- Let Snakemake own output completeness and keep `rerun-incomplete: true` in
  execution profiles.
- Never write a new parameter attempt into an existing scientific run
  namespace; use a new `run_id`.

## Parallelism and resources

- Preserve independent lane and sample jobs so Snakemake can schedule them
  concurrently.
- Declare realistic `threads` and `mem_mb` on every compute-heavy rule.
- Avoid nested oversubscription: one rule's subprocesses share its declared
  thread allocation.
- Make aggregate concurrency a CLI/profile concern (`--cores`, `--jobs`, and
  resource limits), not a hard-coded workflow limit.

## SLURM

- All site-specific SLURM scripts and profiles belong under ignored `slurm/`.
- Do not commit cluster hostnames, accounts, partitions, credentials, or home
  directory paths.
- Never run download, alignment, peak calling, or environment installation on
  the `cranex` login node. Submit compute work through SLURM.
- Long cluster runs must write logs under `logs/` and be safe to resubmit.

## Data and Git hygiene

Do not commit:

- `.venv/`, `.micromamba/`, `.snakemake/`, caches, or locks;
- raw FASTQs, references, BAM/CRAM files, BigWigs, bedGraphs, or generated
  results;
- SLURM scripts or site profiles;
- secrets, credentials, private metadata, or large temporary fixtures.

The versioned `resources/` directory is intentionally limited to adapter FASTAs
and the curated ATAC/H3K27ac atlas metadata and sample sheets. Do not add random
download selections or generated smoke-test data there.

Preserve unrelated user changes in a dirty worktree. Use `rg`/`rg --files` for
search and `apply_patch` for focused file edits.

## Verification

Run the smallest relevant checks during development, then the full practical
suite before handoff:

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
export XDG_CACHE_HOME="$PWD/.cache"

mamba run --prefix "$PWD/.venv" pytest -q

mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --lint

mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --cores 8 --dry-run

git diff --check
```

The loopback aria2 integration test may require network/socket permission in a
sandbox. Report skipped tests and environmental failures honestly.

For schema/configuration changes, include positive and negative validation
tests. For workflow changes, include a Snakemake dry-run assertion for every
new rule branch. For refinement changes, test interval boundaries, threshold
semantics, empty inputs, and deterministic output ordering.

## Documentation

Update `README.md` whenever setup, inputs, commands, defaults, dependencies,
outputs, or restart behavior changes. Update `resources/README.md` when atlas
or adapter provenance changes. Regenerate `docs/workflow-dag.svg` after
changing rule dependencies.

Final handoffs must state what changed, what was tested, what passed or failed,
what was not tested, and any remaining scientific uncertainty.
