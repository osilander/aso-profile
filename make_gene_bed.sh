#!/usr/bin/env bash
set -euo pipefail

GENES_TXT="${1:-genes.txt}"
GTF_GZ="${2:-gencode.v48.basic.annotation.gtf.gz}"
OUTDIR="debug_gene_bed_v3"

mkdir -p "$OUTDIR"

echo "== Inputs =="
echo "Genes: $GENES_TXT"
echo "GTF:   $GTF_GZ"
echo

if [[ ! -s "$GENES_TXT" ]]; then
    echo "ERROR: missing or empty gene list: $GENES_TXT" >&2
    exit 1
fi

if [[ ! -s "$GTF_GZ" ]]; then
    echo "ERROR: missing or empty GTF: $GTF_GZ" >&2
    exit 1
fi

echo "== Normalising requested genes =="
sed 's/\r$//' "$GENES_TXT" \
  | awk '{$1=$1; print}' \
  | grep -v '^[[:space:]]*$' \
  | sort -u \
  > "$OUTDIR/requested_genes.txt"

wc -l "$OUTDIR/requested_genes.txt"
cat "$OUTDIR/requested_genes.txt"
echo

echo "== Checking gzip =="
gzip -t "$GTF_GZ"
echo "GTF gzip OK"
echo

echo "== Extracting all GTF gene rows first =="
zcat "$GTF_GZ" \
  | awk 'BEGIN{FS=OFS="\t"} $3=="gene"{print}' \
  > "$OUTDIR/gtf_gene_rows.tsv"

echo "Gene rows:"
wc -l "$OUTDIR/gtf_gene_rows.tsv"
echo "First 3 gene rows:"
head -3 "$OUTDIR/gtf_gene_rows.tsv"
echo

echo "== Testing gene_name extraction =="
awk 'BEGIN{FS=OFS="\t"}
    {
      attr=$9
      gene=attr
      sub(/^.*gene_name "/, "", gene)
      sub(/".*$/, "", gene)
      print $1, $4, $5, gene, attr
    }' "$OUTDIR/gtf_gene_rows.tsv" \
  > "$OUTDIR/gene_name_extraction_test.tsv"

echo "First 10 extracted names:"
head "$OUTDIR/gene_name_extraction_test.tsv"
echo

echo "== Extracting all gene coordinates =="
awk 'BEGIN{FS=OFS="\t"}
    {
      attr=$9
      gene=attr
      gene_id=attr
      gene_type=attr

      sub(/^.*gene_name "/, "", gene)
      sub(/".*$/, "", gene)

      sub(/^.*gene_id "/, "", gene_id)
      sub(/".*$/, "", gene_id)

      sub(/^.*gene_type "/, "", gene_type)
      sub(/".*$/, "", gene_type)

      if (gene != attr && gene != "") {
        print $1, $4-1, $5, gene, gene_id, gene_type, ".", $7
      }
    }' "$OUTDIR/gtf_gene_rows.tsv" \
  > "$OUTDIR/all_gene_coords.tsv"

echo "All extracted gene coordinate rows:"
wc -l "$OUTDIR/all_gene_coords.tsv"
echo "First 10:"
head "$OUTDIR/all_gene_coords.tsv"
echo

echo "== Checking requested genes are present =="
cut -f4 "$OUTDIR/all_gene_coords.tsv" | sort -u > "$OUTDIR/all_gene_names.txt"

grep -wFf "$OUTDIR/requested_genes.txt" "$OUTDIR/all_gene_names.txt" \
  > "$OUTDIR/found_genes.txt" || true

echo "Found:"
wc -l "$OUTDIR/found_genes.txt"
cat "$OUTDIR/found_genes.txt"
echo

echo "Missing:"
comm -23 "$OUTDIR/requested_genes.txt" "$OUTDIR/found_genes.txt" \
  > "$OUTDIR/missing_genes.txt" || true
cat "$OUTDIR/missing_genes.txt"
echo

echo "== Pulling coordinates for requested genes =="
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {wanted[$1]=1; next}
     ($4 in wanted) {print}
    ' "$OUTDIR/requested_genes.txt" \
      "$OUTDIR/all_gene_coords.tsv" \
  > "$OUTDIR/matched_gene_coords.tsv"

echo "Matched rows:"
wc -l "$OUTDIR/matched_gene_coords.tsv"
cat "$OUTDIR/matched_gene_coords.tsv"
echo

echo "== Writing BED6 =="
# all_gene_coords columns:
# 1 chrom, 2 start0, 3 end, 4 gene, 5 gene_id, 6 gene_type, 7 score, 8 strand
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$7,$8}' \
  "$OUTDIR/matched_gene_coords.tsv" \
  | sort -k1,1 -k2,2n \
  > "$OUTDIR/genes.bed"

echo "genes.bed:"
wc -l "$OUTDIR/genes.bed"
cat "$OUTDIR/genes.bed"
echo

echo "== Standard chromosomes only =="
grep -E '^chr([0-9]+|X|Y|M)[[:space:]]' "$OUTDIR/genes.bed" \
  > "$OUTDIR/genes.standard.bed" || true

echo "genes.standard.bed:"
wc -l "$OUTDIR/genes.standard.bed"
cat "$OUTDIR/genes.standard.bed"
echo

echo "== Done =="
echo "Use:"
echo "  $OUTDIR/genes.standard.bed"
