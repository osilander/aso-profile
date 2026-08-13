# Common SNP ASO targets for dominant disease genes

## Overview

Antisense oligonucleotide (ASO) therapies for dominant genetic diseases are typically designed to target the causal pathogenic variant directly, limiting each therapy to patients who carry that specific allele. An alternative is to exploit the fact that a pathogenic allele is always inherited on a haplotype carrying multiple linked heterozygous variants. Any transcribed heterozygous variant in cis—common intronic SNPs, exonic SNPs, UTR variants—becomes a candidate ASO target for allele-selective knockdown. High-MAF SNPs (≥ 0.1) are particularly attractive because a single pre-validated ASO design would apply to every patient heterozygous at that site.

This repository builds gene-level catalogues of candidate ASO targets for 20 dominant disease genes: pathogenic/likely pathogenic missense SNVs from ClinVar and common SNPs (MAF ≥ 0.1) from 1000 Genomes, scored for ASO suitability using 17-mer sequence windows, RNAfold secondary-structure predictions, allele-selectivity ΔΔG (RNAduplex), and pre-mRNA off-target search (Bowtie2). Coverage of the 1000 Genomes cohort by greedy sets of common-SNP ASO targets is computed per gene, before and after specificity filtering.

## Focal genes

```
SCN2A  SCN8A  KCNQ2  KIF1A  SPTLC1  SPTAN1
COL1A1 COL1A2 COL2A1 KCNJ11 ABCC8   PIK3CA
AKT1   MTOR   STAT1  STAT3  NLRP3   CARD11
RELA   KIF5B
```

## Dependencies

| Tool | Purpose |
|------|---------|
| bcftools ≥ 1.22 | VCF extraction and querying |
| bedtools ≥ 2.31 | BED interval operations |
| ViennaRNA ≥ 2.4 (RNAfold, RNAduplex) | RNA structure and allele selectivity |
| Bowtie2 ≥ 2.5 | Off-target alignment |
| Python ≥ 3.8 | Off-target query generation and parsing |
| R ≥ 4.2 | Coverage analysis and figure generation |

R packages: `ggplot2`, `ggrepel`, `cowplot`

## Data sources

**GENCODE v48** (basic annotation) — place at `data/gencode.v48.basic.annotation.gtf.gz`:
```bash
wget -P data/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gtf.gz
```

**ClinVar** (GRCh38):
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz{,.tbi}
wget -P data/ https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
```

**1000 Genomes high-coverage GRCh38 VCFs** — place per-chromosome VCFs under `data/1000G_highcov_GRCh38/`.

**GRCh38 reference FASTA** — place at `refs/GRCh38.primary_assembly.genome.fa` (with `.fai` index).

## Running the pipeline

### 1. Build gene BED file

```bash
./make_gene_bed.sh
# outputs: gene_bed/genes.standard.bed
sed 's/^chr//' gene_bed/genes.standard.bed > gene_bed/genes.standard.nochr.bed
```

### 2. Build variant tables and sequence windows

```bash
./build_aso_variant_tables.sh
# key output: results/common_gene_body_snps.af0.1.with_gene.rates.tsv
```

### 3. Add transcript-context annotation

```bash
./make_common_snp_context.sh
```

### 4. Score 17-mer ASO windows

```bash
./score_17mer_aso_windows.sh
```

### 5. Run RNAfold

Local:
```bash
./run_rnafold_and_join_aso_scores.sh
```

On SLURM:
```bash
sbatch run_rnafold_201bp.sl   # edit --account before submitting
```

### 6. Add allele-selectivity ΔΔG and sequence motif scores

```bash
./add_extended_aso_scores.sh
# output: results/aso_17mer_scores_extended.tsv
```

On HPC, generate duplex inputs locally, run `run_rduplex_aso.sl` on cluster, then join:
```bash
SKIP_RDUPLEX=1 ./add_extended_aso_scores.sh
```

### 7. Build gene plot inputs

```bash
./make_all_gene_variant_plot_inputs.sh genes.txt
# outputs: gene_plot_inputs/
```

### 8. Off-target search

```bash
python scripts/01_make_aso_queries.py   # build allele-specific 17-mer query FASTAs
bash scripts/02_build_indices.sh        # build Bowtie2 pre-mRNA index (needs ref FASTA)
sbatch scripts/03_search_offtargets.sl  # align queries (SLURM; edit --account)
python scripts/04_parse_offtargets.py   # parse SAM → offtarget_summary.tsv
bash scripts/05_annotate_repeats.sh     # annotate target windows with RepeatMasker
# output: results/offtarget/offtarget_summary.tsv
```

### 9. Coverage curves

```bash
Rscript coverage_curves.R
# outputs: figures/coverage_specificity.csv
#          figures/pop_greedy_specificity.csv  (requires vcf_cache/ RDS files)
```

### 10. Generate figures and tables

```bash
Rscript plot_fig1.R       # Fig 1: example locus (STAT1/PIK3CA)
Rscript plot_fig2.R       # Fig 2: RNA accessibility and 4-mer context
Rscript plot_fig3.R       # Fig 3: allele selectivity ECDFs
Rscript plot_fig4.R       # Fig 4: coverage curves and specificity filtering
Rscript plot_figS2_figS4.R  # Fig S2: ancestry-stratified filtered; Fig S4: random vs greedy
Rscript plot_figS3.R      # Fig S3: unfiltered vs filtered coverage per gene
Rscript plot_figS5.R      # Fig S5: pre-mRNA off-target summary
Rscript make_tables.R     # Tables S1–S2
```

Final PDFs land in `submission/`.

## Repository structure

### Pipeline scripts

```
make_gene_bed.sh                    gene body BED from GENCODE
build_aso_variant_tables.sh         ClinVar P/LP SNVs + common 1000G SNPs (master pipeline)
make_common_snp_context.sh          add CDS/exon/intron context to common SNP table
score_17mer_aso_windows.sh          17-mer sequence scoring
run_rnafold_and_join_aso_scores.sh  RNAfold 201 bp + join to score table
run_rnafold_201bp.sl                SLURM job for RNAfold (edit --account)
add_extended_aso_scores.sh          allele-selectivity ΔΔG + motif scores
run_rduplex_aso.sl                  SLURM job for RNAduplex (edit --account)
make_all_gene_variant_plot_inputs.sh  transcript models and variant tables
run_pipeline_af0.1.sh               HPC SLURM pipeline driver (AF ≥ 0.1 analysis)
run_coverage_curves.sl              SLURM wrapper for coverage_curves.R
coverage_curves.R                   greedy coverage curves from 1000G het matrices
```

### Off-target search scripts (`scripts/`)

```
01_make_aso_queries.py    build allele-specific 17-mer query FASTAs
02_build_indices.sh       build Bowtie2 transcriptome and pre-mRNA indices
03_search_offtargets.sl   Bowtie2 alignment (SLURM)
04_parse_offtargets.py    parse SAM → off-target summary TSV
05_annotate_repeats.sh    RepeatMasker overlap annotation
```

### Figure scripts

```
plot_fig1.R          Fig 1: example locus variant facets
plot_fig2.R          Fig 2: RNA accessibility and 4-mer context
plot_fig3.R          Fig 3: allele-selectivity ECDFs and ΔΔG asymmetry
plot_fig4.R          Fig 4: coverage curves, filtering effect, off-target breakdown
plot_figS2_figS4.R   Fig S2 (ancestry-stratified filtered) + Fig S4 (random vs greedy)
plot_figS3.R         Fig S3: unfiltered vs filtered per-gene coverage
plot_figS5.R         Fig S5: pre-mRNA off-target ECDFs and bar summaries
make_tables.R        Tables S1–S2: coverage thresholds and candidate counts
```

### Committed data snapshots

These are reproducible pipeline outputs committed for convenience; regenerating them requires HPC cluster jobs or large reference files.

```
results/aso_17mer_scores_extended.tsv     scored 17-mer candidates (all genes)
results/offtarget/offtarget_summary.tsv   per-query off-target counts and filter flags
results/offtarget/query_meta.tsv          query-to-variant metadata mapping
results/fourmer_scores.tsv                4-mer motif scores
gene_plot_inputs/                         transcript models and variant tables (all genes)
figures/coverage_specificity.csv          greedy coverage curves before/after filtering
figures/pop_greedy_specificity.csv        per-superpopulation coverage curves
vcf_cache/                                per-gene heterozygosity matrices (1000G)
integrated_call_samples_v3.20130502.ALL.panel  1000G sample–population assignments
```

### Curated inputs

```
selected_STAT1_PIK3CA_contexts.tsv   hand-selected ASO candidates for STAT1/PIK3CA
genes.txt                            list of 20 focal genes
gene_bed/genes.standard.bed          gene body BED (chr-prefixed)
gene_bed/genes.standard.nochr.bed    no-chr version for ClinVar intersection
```

### Submission

```
submission/    final figures (fig1–4, figS1–5) and supplementary tables (S1–S2)
```

## License

MIT
