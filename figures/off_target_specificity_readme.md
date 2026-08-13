# Off-target and repeat-specificity triage for reusable SNP ASOs

This is a proposed analysis layer to address the reviewer-facing concern that ASO
efficacy features are not sufficient: reusable ASOs also need sequence-specificity
triage. The OligoAI/ASO Atlas-derived 4-mer score is useful as an efficacy/context
proxy, but it is not an off-target model. Specificity should therefore be handled
as a separate filter.

## Core question

For candidate 17-mer ASO target sites in the same disease genes, do common SNP
targets have a lower, similar, or higher predicted off-target burden than
pathogenic missense SNV targets?

The useful comparison is not whether a gene has paralogs. It is whether the local
17-mer target sequence has near-complementary hits outside the intended gene.
Common SNP candidates are mostly intronic, so they may have fewer paralog/coding
sequence hits than pathogenic missense sites, but they may also have more repeat
or low-complexity liabilities.

## Inputs

- `aso_17mer_scores_extended.tsv`
  - candidate common SNP and pathogenic missense 17-mers
  - includes gene, variant class, variant ID, ref/alt, strand, context, basic
    sequence filters, and allele-specific scoring fields
- Reference genome FASTA, e.g. GRCh38
- Transcriptome FASTA
  - mature transcriptome for mRNA off-targets
  - optional pre-mRNA/gene-body FASTA for nuclear RNase H gapmer off-targets
- Repeat annotation BED, e.g. RepeatMasker
- Optional paralog annotation table, e.g. Ensembl paralogs or HGNC gene families
- Gene coordinates BED for mapping intended gene versus off-gene hits

## Recommended outputs

- `results/offtarget/aso_candidates.allele_specific.tsv`
- `results/offtarget/aso_candidates.repeat_annotated.tsv`
- `results/offtarget/offtarget_hits.unmasked.tsv`
- `results/offtarget/offtarget_counts.unmasked.tsv`
- `results/offtarget/offtarget_counts.repeat_aware.tsv`
- `results/offtarget/common_snp_specificity_filtered.tsv`
- `figures/offtarget_specificity.pdf`
- `figures/coverage_specificity_filtered.pdf`

## Analysis outline

### 1. Build allele-specific ASO candidate table

For each candidate variant, emit one row per therapeutically relevant allele.

```text
for each row in aso_17mer_scores_extended.tsv:
  if variant_class == common_snp:
    emit candidate targeting REF allele
    emit candidate targeting ALT allele

  if variant_class == pathogenic_missense:
    emit candidate targeting pathogenic ALT allele
    optionally emit reciprocal WT/REF candidate for comparison only

for each emitted candidate:
  target_17mer = allele-specific target RNA/DNA sequence
  aso_seq = reverse_complement(target_17mer)
  candidate_id = gene | variant_id | targeted_allele
```

Keep the reciprocal WT candidate for pathogenic SNVs separate from the
therapeutically relevant pathogenic-allele candidate. It is useful for selectivity
asymmetry, but it is not usually the candidate one would treat with.

### 2. Annotate repeat overlap at the intended target site

Do not only hard-mask the candidate sequences. First annotate whether the intended
17-mer overlaps a repeat.

```text
intersect candidate 17-mer genomic intervals with RepeatMasker BED

add columns:
  repeat_overlap_any
  repeat_class
  repeat_family
  repeat_fraction_of_17mer
  repeat_overlap_variant_base
```

This allows three useful summaries:

- all candidates
- candidates excluding repeat-overlapping intended sites
- off-target hits separated into repeat and non-repeat hits

### 3. Search for near matches without masking

Align candidate target sequences against the unmasked reference. The exact query
orientation depends on the aligner setup; document whether the query is the
target 17-mer or the reverse-complement ASO sequence.

Suggested searches:

```text
exact full-length matches
full-length matches with <=1 mismatch
full-length matches with <=2 mismatches
hits with long perfect contiguous complement stretches, e.g. >=12, >=13, >=14 nt
```

Possible tools:

```text
bowtie -v 0 / -v 1 / -v 2
bowtie2 with short-read sensitive settings
blastn-short
custom exact/near-exact scanner for transcript FASTA
```

Search targets:

```text
genome.fa
mature_transcriptome.fa
premrna_or_gene_body.fa
```

### 4. Classify each hit

For every hit, annotate whether it is the intended site, another site in the
intended gene, or an off-gene hit.

```text
for each alignment hit:
  classify as:
    intended_site
    same_gene_other_site
    outside_intended_gene

  annotate:
    mismatch_count
    longest_perfect_match
    hit_gene
    hit_transcript
    hit_region_context: CDS / UTR / intron / intergenic / noncoding RNA
    repeat_overlap
    paralog_of_intended_gene
```

For the main specificity analysis, the most conservative counts are usually
outside-intended-gene hits. Same-gene secondary hits may still matter, but they
are conceptually different from paralog or genome-wide off-targets.

### 5. Summarize off-target burden per candidate

Create one candidate-level table with counts such as:

```text
exact_hits_outside_gene
one_mismatch_hits_outside_gene
two_mismatch_hits_outside_gene
max_contiguous_match_outside_gene
hits_in_paralog_genes_0mm
hits_in_paralog_genes_1mm
hits_in_paralog_genes_2mm
hits_in_repeats
hits_nonrepeat_outside_gene
hits_coding_outside_gene
hits_premrna_outside_gene
```

Then compare groups:

```text
common_snp_intronic
common_snp_exonic_nonCDS
common_snp_CDS
pathogenic_missense
```

Stratify or block by gene where possible so genes with many candidates do not
dominate the result.

### 6. Define a conservative specificity-pass flag

A first-pass filter could be:

```text
pass_specificity =
  repeat_overlap_any == FALSE
  and exact_hits_outside_gene == 0
  and one_mismatch_hits_outside_gene == 0
  and hits_in_paralog_genes_2mm == 0
```

This is intentionally conservative. It should be presented as a triage filter,
not as a complete toxicity prediction.

Also retain softer thresholds for sensitivity analysis:

```text
strict:
  no repeat overlap
  no off-gene exact or <=1 mismatch hits
  no paralog hits <=2 mismatches

moderate:
  no repeat overlap
  no off-gene exact hits
  no paralog exact or 1 mismatch hits

permissive:
  allow repeat-overlap flag but count/report it
  exclude only off-gene exact hits
```

### 7. Recalculate coverage after specificity filtering

Use the existing greedy heterozygosity coverage machinery, but restrict the SNP
pool to common SNPs with at least one allele-specific candidate passing the chosen
specificity filter.

```text
for each gene:
  raw_candidates = all common SNP candidates
  basic_candidates = pass_basic == 1
  specificity_candidates = basic_candidates and pass_specificity == TRUE

  rerun greedy coverage using:
    raw SNP set
    basic-filtered SNP set
    specificity-filtered SNP set

  report coverage at k = 1, 5, 10
  report k required for 80%, 90%, 95% where reached
```

This gives the most useful reviewer-facing result: how much of the apparent
population coverage survives when specificity is added as a real library-design
filter.

## Suggested figures

### Main or supplement figure: off-target specificity

Panels:

```text
A. Fraction of candidate 17-mers overlapping RepeatMasker regions
B. Distribution of off-gene <=1 mismatch hits
C. Distribution of off-gene <=2 mismatch hits
D. Paralog hit burden, grouped by candidate class
```

Compare:

```text
common SNP intronic
common SNP exonic
common SNP CDS
pathogenic missense
```

### Main or supplement figure: coverage shrinkage

Panels:

```text
A. Raw greedy heterozygosity coverage
B. Coverage after basic ASO sequence/context filters
C. Coverage after repeat/off-target specificity filters
```

Or as a compact triage plot:

```text
all common SNP targets
  -> pass basic sequence filters
  -> pass accessibility / 4-mer filters
  -> pass repeat / off-target filters
  -> retained k=5 or k=10 population coverage
```

## Manuscript framing

Useful wording:

```text
Because common-SNP libraries provide many possible allele-discriminating sites
per gene, sequence specificity can be treated as a design constraint rather than
an afterthought. We therefore incorporated repeat overlap, near-perfect off-gene
matches, and paralog cross-reactivity as downstream triage filters before
estimating retained library coverage.
```

Avoid claiming that in silico specificity screens prove safety. They only reduce
obvious hybridization-dependent liabilities and prioritize candidates for
experimental testing.

