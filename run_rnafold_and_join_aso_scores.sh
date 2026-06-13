#!/usr/bin/env bash
set -euo pipefail

# Build 201 bp RNAfold windows around pathogenic missense variants
# and common SNPs, run RNAfold, parse MFE, and join MFE values
# to the existing 17-mer ASO score table.
#
# Required existing inputs:
#   results/clinvar_pathogenic_snvs.missense.with_gene.tsv
#   results/common_gene_body_snps.af0.2.with_gene.rates.context.tsv
#   results/aso_17mer_scores.combined.tsv
#   refs/GRCh38.primary_assembly.genome.fa
#
# Outputs:
#   results/rnafold/pathogenic_missense.plusminus100.bed
#   results/rnafold/common_snps.plusminus100.bed
#   results/rnafold/pathogenic_missense.plusminus100.rna.fa
#   results/rnafold/common_snps.plusminus100.rna.fa
#   results/rnafold/pathogenic_missense.plusminus100.rnafold.txt
#   results/rnafold/common_snps.plusminus100.rnafold.txt
#   results/rnafold/rnafold_plusminus100_mfe.tsv
#   results/aso_17mer_scores_with_rnafold.tsv
#
# Usage:
#   bash run_rnafold_and_join_aso_scores.sh
#
# Optional:
#   bash run_rnafold_and_join_aso_scores.sh refs/GRCh38.primary_assembly.genome.fa

REF_FASTA="${1:-refs/GRCh38.primary_assembly.genome.fa}"

OUTDIR="results"
RNAFOLD_DIR="${OUTDIR}/rnafold"
LOGDIR="logs"

PATHOGENIC_TSV="${OUTDIR}/clinvar_pathogenic_snvs.missense.with_gene.tsv"
COMMON_TSV="${OUTDIR}/common_gene_body_snps.af0.2.with_gene.rates.context.tsv"
SCORE_TSV="${OUTDIR}/aso_17mer_scores.combined.tsv"

mkdir -p "$RNAFOLD_DIR" "$LOGDIR"

echo "###############################################################################"
echo "RNAfold local folding pipeline"
echo "###############################################################################"
echo
echo "Reference FASTA: $REF_FASTA"
echo "Pathogenic TSV:  $PATHOGENIC_TSV"
echo "Common SNP TSV:  $COMMON_TSV"
echo "17-mer scores:   $SCORE_TSV"
echo

for f in "$REF_FASTA" "$PATHOGENIC_TSV" "$COMMON_TSV" "$SCORE_TSV"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing required file: $f" >&2
    exit 1
  fi
done

echo "== Check FASTA index =="
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo "Indexing reference FASTA..."
  samtools faidx "$REF_FASTA"
fi

echo "Reference contigs:"
cut -f1 "${REF_FASTA}.fai" | head
echo

###############################################################################
# Step 1: build ±100 bp BED windows
###############################################################################

echo "###############################################################################"
echo "Step 1: Build 201 bp BED windows"
echo "###############################################################################"

echo "== Pathogenic missense windows =="

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

    win_start=var_start-100
    win_end=var_end+100
    if (win_start < 0) win_start=0

    name="pathogenic_missense|"gene"|"clinvar_id"|"ref"|"alt"|"clnsig"|NA|NA|"strand"|pathogenic_missense"

    print chrom,win_start,win_end,name,".",strand
  }
' "$PATHOGENIC_TSV" \
  > "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.bed"

echo "Pathogenic BED rows:"
wc -l "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.bed"
head "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.bed" | column -t
echo

echo "== Common SNP windows =="

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

    win_start=snp_start-100
    win_end=snp_end+100
    if (win_start < 0) win_start=0

    name="common_snp|"gene"|"snp_id"|"ref"|"alt"|NA|"maf"|"het"|"strand"|"context

    print chrom,win_start,win_end,name,".",strand
  }
' "$COMMON_TSV" \
  > "${RNAFOLD_DIR}/common_snps.plusminus100.bed"

echo "Common SNP BED rows:"
wc -l "${RNAFOLD_DIR}/common_snps.plusminus100.bed"
head "${RNAFOLD_DIR}/common_snps.plusminus100.bed" | column -t
echo

echo "Combined BED row counts:"
wc -l "${RNAFOLD_DIR}"/*.plusminus100.bed
echo

###############################################################################
# Step 2: extract stranded DNA FASTA
###############################################################################

echo "###############################################################################"
echo "Step 2: Extract stranded FASTA"
echo "###############################################################################"

bedtools getfasta \
  -fi "$REF_FASTA" \
  -bed "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.bed" \
  -name \
  -s \
  -fo "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.fa"

bedtools getfasta \
  -fi "$REF_FASTA" \
  -bed "${RNAFOLD_DIR}/common_snps.plusminus100.bed" \
  -name \
  -s \
  -fo "${RNAFOLD_DIR}/common_snps.plusminus100.fa"

echo "Pathogenic FASTA records:"
grep -c '^>' "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.fa"

echo "Common SNP FASTA records:"
grep -c '^>' "${RNAFOLD_DIR}/common_snps.plusminus100.fa"
echo

echo "First common FASTA records:"
head "${RNAFOLD_DIR}/common_snps.plusminus100.fa"
echo

###############################################################################
# Step 3: convert DNA FASTA to RNA FASTA
###############################################################################

echo "###############################################################################"
echo "Step 3: Convert DNA FASTA to RNA FASTA"
echo "###############################################################################"

awk '
  /^>/ {print; next}
  {gsub(/T/,"U"); print}
' "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.fa" \
  > "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.rna.fa"

awk '
  /^>/ {print; next}
  {gsub(/T/,"U"); print}
' "${RNAFOLD_DIR}/common_snps.plusminus100.fa" \
  > "${RNAFOLD_DIR}/common_snps.plusminus100.rna.fa"

echo "Check that common RNA FASTA has no T in first sequence:"
grep -v '^>' "${RNAFOLD_DIR}/common_snps.plusminus100.rna.fa" \
  | head -1 \
  | grep T || true
echo

###############################################################################
# Step 4: run RNAfold
###############################################################################

echo "###############################################################################"
echo "Step 4: Run RNAfold"
echo "###############################################################################"

# Load ViennaRNA if available as an environment module.
# If already loaded, this is harmless on most systems.
if command -v module >/dev/null 2>&1; then
  module load ViennaRNA/2.4.17-gimkl-2020a || true
fi

if ! command -v RNAfold >/dev/null 2>&1; then
  echo "ERROR: RNAfold not found. Try: module load ViennaRNA/2.4.17-gimkl-2020a" >&2
  exit 1
fi

echo "RNAfold path:"
which RNAfold
echo "RNAfold version:"
RNAfold --version || true
echo

echo "Running RNAfold: pathogenic missense..."
RNAfold --noPS < "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.rna.fa" \
  > "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.rnafold.txt"

echo "Running RNAfold: common SNPs..."
RNAfold --noPS < "${RNAFOLD_DIR}/common_snps.plusminus100.rna.fa" \
  > "${RNAFOLD_DIR}/common_snps.plusminus100.rnafold.txt"

echo "RNAfold record counts:"
grep -c '^>' "${RNAFOLD_DIR}/pathogenic_missense.plusminus100.rnafold.txt"
grep -c '^>' "${RNAFOLD_DIR}/common_snps.plusminus100.rnafold.txt"
echo

echo "First common RNAfold records:"
head -6 "${RNAFOLD_DIR}/common_snps.plusminus100.rnafold.txt"
echo

###############################################################################
# Step 5: parse RNAfold MFE
###############################################################################

echo "###############################################################################"
echo "Step 5: Parse RNAfold MFE"
echo "###############################################################################"

python3 - <<'PY'
from pathlib import Path
import re
import csv

files = [
    ("pathogenic_missense", Path("results/rnafold/pathogenic_missense.plusminus100.rnafold.txt")),
    ("common_snp", Path("results/rnafold/common_snps.plusminus100.rnafold.txt")),
]

out = Path("results/rnafold/rnafold_plusminus100_mfe.tsv")
rows = []

for expected_class, path in files:
    print(f"Parsing {path}")

    with open(path) as f:
        lines = [x.rstrip("\n") for x in f if x.strip()]

    if len(lines) % 3 != 0:
        print(f"WARNING: line count not divisible by 3 for {path}: {len(lines)}")

    for i in range(0, len(lines), 3):
        header = lines[i].lstrip(">")
        seq = lines[i + 1]
        struct_line = lines[i + 2]

        # Handles:
        #   (-54.60)
        #   ( -2.60)
        #   (  0.00)
        m = re.search(r"\(\s*([-+]?\d+(?:\.\d+)?)\s*\)\s*$", struct_line)
        mfe = float(m.group(1)) if m else None

        structure = re.sub(
            r"\s+\(\s*[-+]?\d+(?:\.\d+)?\s*\)\s*$",
            "",
            struct_line
        )

        meta = header.split("::")[0]
        parts = meta.split("|")
        while len(parts) < 10:
            parts.append("NA")

        rows.append({
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
            "window_length": len(seq),
            "mfe": mfe if mfe is not None else "NA",
            "mfe_per_nt": (mfe / len(seq)) if mfe is not None else "NA",
            "sequence": seq,
            "structure": structure,
            "full_header": header,
        })

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
    "window_length",
    "mfe",
    "mfe_per_nt",
    "sequence",
    "structure",
    "full_header",
]

with open(out, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(f"Wrote {out}")
print(f"Rows: {len(rows)}")
print(f"Rows with missing MFE: {sum(r['mfe'] == 'NA' for r in rows)}")
PY

echo

echo "Parsed RNAfold table check:"
head "${RNAFOLD_DIR}/rnafold_plusminus100_mfe.tsv" | column -t
wc -l "${RNAFOLD_DIR}/rnafold_plusminus100_mfe.tsv"
awk -F'\t' '{print NF}' "${RNAFOLD_DIR}/rnafold_plusminus100_mfe.tsv" | sort | uniq -c
echo

###############################################################################
# Step 6: join RNAfold MFE to 17-mer score table
###############################################################################

echo "###############################################################################"
echo "Step 6: Join RNAfold MFE to 17-mer score table"
echo "###############################################################################"

python3 - <<'PY'
from pathlib import Path
import csv

score_file = Path("results/aso_17mer_scores.combined.tsv")
fold_file = Path("results/rnafold/rnafold_plusminus100_mfe.tsv")
out_file = Path("results/aso_17mer_scores_with_rnafold.tsv")

def score_key(row):
    return row["full_header"].split("::")[0]

def fold_key(row):
    return row["full_header"].split("::")[0]

fold = {}

with open(fold_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    duplicate = 0
    for row in reader:
        k = fold_key(row)
        if k in fold:
            duplicate += 1
        fold[k] = row

with open(score_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    score_rows = list(reader)
    score_fields = reader.fieldnames

extra_fields = ["rnafold_window_length", "mfe", "mfe_per_nt"]
missing = 0

with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=score_fields + extra_fields,
        delimiter="\t"
    )
    writer.writeheader()

    for row in score_rows:
        k = score_key(row)
        frow = fold.get(k)

        if frow is None:
            missing += 1
            row["rnafold_window_length"] = "NA"
            row["mfe"] = "NA"
            row["mfe_per_nt"] = "NA"
        else:
            row["rnafold_window_length"] = frow["window_length"]
            row["mfe"] = frow["mfe"]
            row["mfe_per_nt"] = frow["mfe_per_nt"]

        writer.writerow(row)

print(f"Wrote {out_file}")
print(f"Score rows: {len(score_rows)}")
print(f"RNAfold unique keys: {len(fold)}")
print(f"Duplicate RNAfold keys: {duplicate}")
print(f"Missing RNAfold matches: {missing}")
PY

echo

###############################################################################
# Step 7: final checks
###############################################################################

echo "###############################################################################"
echo "Step 7: Final checks"
echo "###############################################################################"

FINAL="${OUTDIR}/aso_17mer_scores_with_rnafold.tsv"

echo "Final joined table:"
wc -l "$FINAL"
awk -F'\t' '{print NF}' "$FINAL" | sort | uniq -c
echo

echo "Counts by variant class:"
awk -F'\t' 'NR>1 {count[$1]++} END{for (k in count) print k, count[k]}' "$FINAL" | sort
echo

echo "Missing RNAfold values:"
awk -F'\t' 'NR>1 && ($27=="NA" || $28=="NA") {n++} END{print n+0}' "$FINAL"
echo

echo "First rows:"
head "$FINAL" | column -t
echo

echo "== Done =="
echo "Main output:"
echo "  $FINAL"
