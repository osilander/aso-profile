#!/usr/bin/env bash
set -euo pipefail

# Build and score 17-mers centred on:
#   1. ClinVar pathogenic missense SNVs
#   2. common gene-body SNPs with MAF >= 0.2
#
# The 17-mer is:
#   8 bp upstream + variant/SNP base + 8 bp downstream
#
# Sequence is extracted in gene/transcript orientation using bedtools getfasta -s.
#
# Outputs:
#   results/aso_17mer.pathogenic_missense.bed
#   results/aso_17mer.common_snps.bed
#   results/aso_17mer.pathogenic_missense.fa
#   results/aso_17mer.common_snps.fa
#   results/aso_17mer_scores.combined.tsv
#   results/aso_17mer_scores.top_common_per_gene.tsv

REF_FASTA="${1:-refs/GRCh38.primary_assembly.genome.fa}"
OUTDIR="results"

PATHOGENIC_TSV="${OUTDIR}/clinvar_pathogenic_snvs.missense.with_gene.tsv"
COMMON_TSV="${OUTDIR}/common_gene_body_snps.af0.2.with_gene.rates.context.tsv"

mkdir -p "$OUTDIR"

echo "== Inputs =="
echo "Reference FASTA: $REF_FASTA"
echo "Pathogenic TSV:  $PATHOGENIC_TSV"
echo "Common SNP TSV:  $COMMON_TSV"
echo

if [[ ! -s "$REF_FASTA" ]]; then
  echo "ERROR: missing reference FASTA: $REF_FASTA" >&2
  exit 1
fi

if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo "Indexing reference FASTA..."
  samtools faidx "$REF_FASTA"
fi

if [[ ! -s "$PATHOGENIC_TSV" ]]; then
  echo "ERROR: missing pathogenic TSV: $PATHOGENIC_TSV" >&2
  exit 1
fi

if [[ ! -s "$COMMON_TSV" ]]; then
  echo "ERROR: missing common SNP TSV: $COMMON_TSV" >&2
  exit 1
fi

###############################################################################
# Step 1: make 17-mer BED for pathogenic missense variants
###############################################################################

echo "== Step 1: make pathogenic missense 17-mer BED =="

# Pathogenic with-gene columns:
#  1 VAR_CHROM
#  2 VAR_START0
#  3 VAR_END
#  4 CLINVAR_ID
#  5 REF
#  6 ALT
#  7 GENEINFO
#  8 CLNSIG
#  9 CLNREVSTAT
# 10 MC
# 11 CLNHGVS
# 12 CLNDN
# 13 GENE_CHROM
# 14 GENE_START
# 15 GENE_END
# 16 GENE
# 17 SCORE
# 18 STRAND

awk -F'\t' 'BEGIN{OFS="\t"}
  {
    chrom=$1
    var_start=$2
    var_end=$3
    clinvar_id=$4
    ref=$5
    alt=$6
    clnsig=$8
    gene=$16
    strand=$18

    win_start=var_start-8
    win_end=var_end+8
    if (win_start < 0) win_start=0

    # Metadata fields:
    # variant_class|gene|id|ref|alt|clnsig|maf|het_rate|strand|context
    name="pathogenic_missense|"gene"|"clinvar_id"|"ref"|"alt"|"clnsig"|NA|NA|"strand"|pathogenic_missense"

    print chrom,win_start,win_end,name,".",strand
  }
' "$PATHOGENIC_TSV" \
  > "$OUTDIR/aso_17mer.pathogenic_missense.bed"

wc -l "$OUTDIR/aso_17mer.pathogenic_missense.bed"
head "$OUTDIR/aso_17mer.pathogenic_missense.bed" | column -t
echo

###############################################################################
# Step 2: make 17-mer BED for common SNPs
###############################################################################

echo "== Step 2: make common SNP 17-mer BED =="

# Common context columns:
#  1 SNP_CHROM
#  2 SNP_START0
#  3 SNP_END
#  4 SNP_ID
#  5 REF
#  6 ALT
#  7 AF
#  8 AC
#  9 AN
# 10 AC_HET
# 11 AC_HOM
# 12 MAF
# 13 HET_RATE
# 14 HOM_ALT_RATE
# 15 GENE
# 16 STRAND
# 17 TRANSCRIPT_CONTEXT

awk -F'\t' 'BEGIN{OFS="\t"}
  {
    chrom=$1
    snp_start=$2
    snp_end=$3
    snp_id=$4
    ref=$5
    alt=$6
    maf=$12
    het=$13
    gene=$15
    strand=$16
    context=$17

    win_start=snp_start-8
    win_end=snp_end+8
    if (win_start < 0) win_start=0

    # Metadata fields:
    # variant_class|gene|id|ref|alt|clnsig|maf|het_rate|strand|context
    name="common_snp|"gene"|"snp_id"|"ref"|"alt"|NA|"maf"|"het"|"strand"|"context

    print chrom,win_start,win_end,name,".",strand
  }
' "$COMMON_TSV" \
  > "$OUTDIR/aso_17mer.common_snps.bed"

wc -l "$OUTDIR/aso_17mer.common_snps.bed"
head "$OUTDIR/aso_17mer.common_snps.bed" | column -t
echo

###############################################################################
# Step 3: extract stranded FASTA
###############################################################################

echo "== Step 3: extract stranded FASTA =="

bedtools getfasta \
  -fi "$REF_FASTA" \
  -bed "$OUTDIR/aso_17mer.pathogenic_missense.bed" \
  -name \
  -s \
  -fo "$OUTDIR/aso_17mer.pathogenic_missense.fa"

bedtools getfasta \
  -fi "$REF_FASTA" \
  -bed "$OUTDIR/aso_17mer.common_snps.bed" \
  -name \
  -s \
  -fo "$OUTDIR/aso_17mer.common_snps.fa"

echo "Pathogenic FASTA records:"
grep -c '^>' "$OUTDIR/aso_17mer.pathogenic_missense.fa"

echo "Common SNP FASTA records:"
grep -c '^>' "$OUTDIR/aso_17mer.common_snps.fa"

echo "First pathogenic FASTA:"
head "$OUTDIR/aso_17mer.pathogenic_missense.fa"
echo

echo "First common FASTA:"
head "$OUTDIR/aso_17mer.common_snps.fa"
echo

###############################################################################
# Step 4: score 17-mers
###############################################################################

echo "== Step 4: score 17-mers =="

python3 - <<'PY'
from pathlib import Path
import re
import csv
import math

outdir = Path("results")

fastas = [
    outdir / "aso_17mer.pathogenic_missense.fa",
    outdir / "aso_17mer.common_snps.fa",
]

out_combined = outdir / "aso_17mer_scores.combined.tsv"
out_top_common = outdir / "aso_17mer_scores.top_common_per_gene.tsv"

def read_fasta(path):
    name = None
    seqs = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seqs).upper()
                name = line[1:]
                seqs = []
            else:
                seqs.append(line.strip())
        if name is not None:
            yield name, "".join(seqs).upper()

def longest_homopolymer(seq):
    best = 0
    current = 0
    last = None
    for b in seq:
        if b == last:
            current += 1
        else:
            current = 1
            last = b
        best = max(best, current)
    return best

def gc_fraction(seq):
    if len(seq) == 0:
        return float("nan")
    return (seq.count("G") + seq.count("C")) / len(seq)

def score_17mer(seq):
    seq = seq.upper()
    length = len(seq)

    has_n = any(b not in "ACGT" for b in seq)

    if length != 17 or has_n:
        return {
            "pass_basic": 0,
            "score_total": 0,
            "score_gc": 0,
            "score_homopolymer": 0,
            "score_dinuc_diversity": 0,
            "score_cpg": 0,
            "score_central9_gc": 0,
            "gc_frac": gc_fraction(seq),
            "central9_gc_frac": float("nan"),
            "max_homopolymer": longest_homopolymer(seq),
            "distinct_dinucs": 0,
            "cpg_count": seq.count("CG"),
            "fail_reason": "bad_length_or_N"
        }

    gc = gc_fraction(seq)
    central9 = seq[4:13]  # 0-based positions 5-13; centred on base 9
    central9_gc = gc_fraction(central9)
    max_run = longest_homopolymer(seq)
    dinucs = {seq[i:i+2] for i in range(len(seq)-1)}
    distinct_dinucs = len(dinucs)
    cpg_count = seq.count("CG")

    # GC score
    if 0.40 <= gc <= 0.60:
        score_gc = 3
    elif 0.30 <= gc <= 0.70:
        score_gc = 1
    else:
        score_gc = 0

    # Homopolymer score
    if max_run <= 3:
        score_homopolymer = 2
    elif max_run == 4:
        score_homopolymer = 1
    else:
        score_homopolymer = 0

    # Dinucleotide diversity score
    if distinct_dinucs >= 6:
        score_dinuc = 2
    elif distinct_dinucs >= 4:
        score_dinuc = 1
    else:
        score_dinuc = 0

    # CpG score
    score_cpg = 1 if cpg_count <= 2 else 0

    # Central 9-mer GC score
    if 0.33 <= central9_gc <= 0.67:
        score_central9 = 2
    elif 0.22 <= central9_gc <= 0.78:
        score_central9 = 1
    else:
        score_central9 = 0

    total = score_gc + score_homopolymer + score_dinuc + score_cpg + score_central9

    return {
        "pass_basic": 1,
        "score_total": total,
        "score_gc": score_gc,
        "score_homopolymer": score_homopolymer,
        "score_dinuc_diversity": score_dinuc,
        "score_cpg": score_cpg,
        "score_central9_gc": score_central9,
        "gc_frac": gc,
        "central9_gc_frac": central9_gc,
        "max_homopolymer": max_run,
        "distinct_dinucs": distinct_dinucs,
        "cpg_count": cpg_count,
        "fail_reason": "pass"
    }

def parse_header(header):
    # bedtools getfasta -name -s turns:
    # name::chr:start-end(strand)
    # into the FASTA header.
    # We only want the metadata before "::".
    meta = header.split("::")[0]
    parts = meta.split("|")

    # Expected:
    # variant_class|gene|id|ref|alt|clnsig|maf|het_rate|strand|context
    if len(parts) < 10:
        parts = parts + ["NA"] * (10 - len(parts))

    return {
        "variant_class": parts[0],
        "gene": parts[1],
        "variant_id": parts[2],
        "ref": parts[3],
        "alt": parts[4],
        "clnsig": parts[5],
        "maf": parts[6],
        "het_rate": parts[7],
        "strand": parts[8],
        "context": parts[9],
        "full_header": header,
    }

rows = []

for fasta in fastas:
    for header, seq in read_fasta(fasta):
        meta = parse_header(header)
        scores = score_17mer(seq)

        row = {}
        row.update(meta)
        row["seq_17mer"] = seq
        row.update(scores)
        rows.append(row)

fieldnames = [
    "variant_class",
    "gene",
    "variant_id",
    "ref",
    "alt",
    "clnsig",
    "maf",
    "het_rate",
    "strand",
    "context",
    "seq_17mer",
    "pass_basic",
    "score_total",
    "score_gc",
    "score_homopolymer",
    "score_dinuc_diversity",
    "score_cpg",
    "score_central9_gc",
    "gc_frac",
    "central9_gc_frac",
    "max_homopolymer",
    "distinct_dinucs",
    "cpg_count",
    "fail_reason",
    "full_header",
]

with open(out_combined, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

# Top common SNP per gene:
# Sort by score_total, then het_rate, then maf.
common = [r for r in rows if r["variant_class"] == "common_snp" and r["pass_basic"] == 1]

def to_float(x):
    try:
        return float(x)
    except Exception:
        return float("-inf")

best = {}
for r in common:
    g = r["gene"]
    key = (
        int(r["score_total"]),
        to_float(r["het_rate"]),
        to_float(r["maf"]),
    )
    if g not in best:
        best[g] = (key, r)
    elif key > best[g][0]:
        best[g] = (key, r)

with open(out_top_common, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for g in sorted(best):
        writer.writerow(best[g][1])

print(f"Wrote {out_combined}")
print(f"Wrote {out_top_common}")
print(f"Total scored rows: {len(rows)}")
print(f"Pathogenic rows: {sum(r['variant_class']=='pathogenic_missense' for r in rows)}")
print(f"Common SNP rows: {sum(r['variant_class']=='common_snp' for r in rows)}")
PY

echo

###############################################################################
# Step 5: sanity checks
###############################################################################

echo "== Step 5: sanity checks =="

echo "Combined score table:"
wc -l "$OUTDIR/aso_17mer_scores.combined.tsv"
head "$OUTDIR/aso_17mer_scores.combined.tsv" | column -t
echo

echo "Column counts:"
awk -F'\t' '{print NF}' "$OUTDIR/aso_17mer_scores.combined.tsv" | sort | uniq -c
echo

echo "Counts by variant class:"
awk -F'\t' 'NR>1 {count[$1]++} END{for (k in count) print k, count[k]}' \
  "$OUTDIR/aso_17mer_scores.combined.tsv" \
  | sort
echo

echo "Score distribution:"
awk -F'\t' 'NR>1 {count[$13]++} END{for (s in count) print s, count[s]}' \
  "$OUTDIR/aso_17mer_scores.combined.tsv" \
  | sort -n
echo

echo "Top common SNP per gene:"
column -t "$OUTDIR/aso_17mer_scores.top_common_per_gene.tsv"
echo

echo "== Done =="
echo "Main output for R:"
echo "  $OUTDIR/aso_17mer_scores.combined.tsv"
echo
echo "Top common SNPs for larger points:"
echo "  $OUTDIR/aso_17mer_scores.top_common_per_gene.tsv"
