#!/usr/bin/env bash
set -euo pipefail

# Add transcript-context annotation (CDS / exon_nonCDS / intron) to the
# common-SNP rates table produced by build_aso_variant_tables.sh.
#
# Input:  results/common_gene_body_snps.af0.2.with_gene.rates.tsv  (16 cols, no header)
# Output: results/common_gene_body_snps.af0.2.with_gene.rates.context.tsv  (17 cols, no header)
#
# Usage:
#   ./make_common_snp_context.sh
#   MIN_MAF=0.2 GTF=data/gencode.v48.basic.annotation.gtf.gz ./make_common_snp_context.sh

MIN_MAF="${MIN_MAF:-0.2}"
RATES="${1:-results/common_gene_body_snps.af${MIN_MAF}.with_gene.rates.tsv}"
GTF="${GTF:-data/gencode.v48.basic.annotation.gtf.gz}"
OUT="results/common_gene_body_snps.af${MIN_MAF}.with_gene.rates.context.tsv"

mkdir -p results

for f in "$RATES" "$GTF"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing required file: $f" >&2
    exit 1
  fi
done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "== Extracting exon and CDS union BEDs from GTF =="
zcat "$GTF" \
  | awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "$TMP/exons.bed"

zcat "$GTF" \
  | awk 'BEGIN{OFS="\t"} $3=="CDS"{print $1,$4-1,$5}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "$TMP/cds.bed"

echo "Exon union intervals: $(wc -l < "$TMP/exons.bed")"
echo "CDS union intervals:  $(wc -l < "$TMP/cds.bed")"
echo

echo "== Finding CDS and exon hits =="
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR}' "$RATES" > "$TMP/snps.bed"

bedtools intersect -a "$TMP/snps.bed" -b "$TMP/cds.bed"  -wa | cut -f4 | sort -u > "$TMP/cds_rows.txt"
bedtools intersect -a "$TMP/snps.bed" -b "$TMP/exons.bed" -wa | cut -f4 | sort -u > "$TMP/exon_rows.txt"

echo "CDS hits:  $(wc -l < "$TMP/cds_rows.txt")"
echo "Exon hits: $(wc -l < "$TMP/exon_rows.txt")"
echo

echo "== Writing context column =="
python3 - "$RATES" "$TMP/cds_rows.txt" "$TMP/exon_rows.txt" "$OUT" <<'PY'
import sys

rates_path = sys.argv[1]
cds_path   = sys.argv[2]
exon_path  = sys.argv[3]
out_path   = sys.argv[4]

cds_rows  = {int(x.strip()) for x in open(cds_path)  if x.strip()}
exon_rows = {int(x.strip()) for x in open(exon_path) if x.strip()}

written = 0
with open(rates_path) as fin, open(out_path, "w") as fout:
    for i, line in enumerate(fin, start=1):
        ctx = "CDS" if i in cds_rows else ("exon_nonCDS" if i in exon_rows else "intron")
        fout.write(line.rstrip("\n") + "\t" + ctx + "\n")
        written += 1

print(f"Wrote {written} rows to {out_path}")
PY

echo
echo "Context distribution:"
cut -f17 "$OUT" | sort | uniq -c | sort -nr
echo
echo "Output: $OUT"
