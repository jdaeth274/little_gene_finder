# little_gene_finder

A lightweight pipeline for finding pneumococcal capsule (`cps`) genes in many
isolate assemblies, extracting the best nucleotide hit(s), and creating one
alignment per gene.

The intended input layout is:

```text
inputs/
  assemblies/
    ERR001.fasta
    ERR002.fasta
  refs/
    cps_genes.fna
    cps_genes.faa
```

## Why BLAST/tBLASTn rather than Prokka-only?

Prokka annotations are useful for intact genes, but disrupted capsule genes can
be missed, split into multiple CDS features, or assigned a different name. This
pipeline searches the raw assemblies directly with `blastn` and/or `tblastn`,
then merges non-overlapping HSPs per gene/isolate. That means a frameshifted,
truncated, or contig-broken gene can still be classified as `partial` or
`fragmented` and extracted for inspection/alignment.

A sensible default strategy is:

1. Use `blastn` against `cps_genes.fna` for close, intact nucleotide matches.
2. Use `tblastn` against `cps_genes.faa` to rescue frameshifted or divergent
   coding sequences from the assembly DNA.
3. Keep the best coverage/source per gene and isolate.
4. Classify hits as:
   - `complete`: high query coverage in one non-fragmented hit set.
   - `fragmented`: high query coverage split across HSPs and/or contigs.
   - `partial`: lower-but-useful query coverage.
   - `absent`: no passing BLAST HSPs.

## Requirements

Install BLAST+ and MAFFT on your `PATH`:

- `makeblastdb`
- `blastn`
- `tblastn`
- `mafft` (unless using `--skip-align`)

The Python script uses only the Python standard library. The supplied
`environment.yml` installs Python, BLAST+, and MAFFT from conda-forge/bioconda.

With mamba:

```bash
mamba env create -f environment.yml
mamba activate little_gene_finder
```

With pixi, import the same conda environment file into a new workspace:

```bash
pixi init --import environment.yml
pixi shell
```

If you already have a pixi workspace, run `pixi import environment.yml` instead.

## Usage

```bash
scripts/cps_gene_pipeline.py \
  --assemblies inputs/assemblies \
  --nucl-refs inputs/refs/cps_genes.fna \
  --prot-refs inputs/refs/cps_genes.faa \
  --outdir results/cps \
  --mode both \
  --threads 8
```

For a quick run without alignments:

```bash
scripts/cps_gene_pipeline.py --skip-align --threads 8
```

Important tunable thresholds:

- `--complete-qcov 0.90`: query coverage needed for `complete` or
  `fragmented` calls.
- `--partial-qcov 0.20`: minimum query coverage for `partial` calls.
- `--min-identity 70`: minimum HSP identity.
- `--min-hsp-len 30`: minimum HSP length (nucleotides for `blastn`, amino acids
  for `tblastn`).
- `--fragmented-gap 30`: query gap that marks a multi-HSP hit as fragmented.

## Outputs

The pipeline writes these files below `results/cps/` by default:

```text
blastdb/                 per-isolate nucleotide BLAST databases
blast/                   raw blastn/tblastn tabular output
summary.tsv              one row per isolate/gene with status and metrics
cps_hits.bed             BED6 coordinates for selected HSPs
sequences/               extracted hit FASTA per isolate
per_gene/                unaligned extracted FASTA split by gene
alignments/              MAFFT alignment per gene
```

`summary.tsv` includes status, selected search source, query coverage, mean
identity, HSP count, contig count, and the largest gap between selected HSPs on
the query.

For fragmented calls, selected HSP nucleotide sequences are oriented to the gene,
ordered by query coordinate, and concatenated with a configurable run of `N`s
(default: `--gap-n 20`). The original per-HSP assembly coordinates remain in
`cps_hits.bed`.

## Notes and caveats

- FASTA record IDs in `cps_genes.fna` and `cps_genes.faa` should use matching
  gene names if you run `--mode both`; the first whitespace-delimited token is
  used as the gene ID.
- `tblastn` coordinates are nucleotide coordinates on the assembly, but query
  coverage is measured over the protein query.
- The script chooses a non-overlapping set of HSPs greedily by bit score, then
  chooses the better of `blastn` and `tblastn` by status, coverage, score, and
  identity.
- Review `fragmented` and `partial` calls manually if they are biologically
  important, especially in repetitive capsule loci.
