#!/usr/bin/env bash
set -euo pipefail

# Extract common 1000 Genomes SNPs from selected gene bodies.
#
# Assumes:
#   - genes.standard.bed uses chr-prefixed chromosome names, e.g. chr1
#   - 1000G VCFs also use chr-prefixed names, e.g. chr1
#   - 1000G VCFs contain INFO/AF
#
# Usage:
#   bash make_common_1000g_gene_snps.sh
#
# Optional usage:
#   bash make_common_1000g_gene_snps.sh genes.standard.bed data/1000G_highcov_GRCh38 results/common_snps 0.2

GENES_BED="${1:-genes.standard.bed}"
VCF_DIR="${2:-data/1000G_highcov_GRCh38}"
OUTDIR="${3:-results/common_snps}"
MIN_AF="${4:-0.2}"

# For MAF >= 0.2 using alternate allele frequency, keep AF between 0.2 and 0.8.
MAX_AF=$(awk -v min_af="$MIN_AF" 'BEGIN{print 1-min_af}')

mkdir -p "$OUTDIR"
mkdir -p results
mkdir -p logs

echo "== Inputs =="
echo "Genes BED: $GENES_BED"
echo "VCF dir:   $VCF_DIR"
echo "Outdir:    $OUTDIR"
echo "MIN_AF:    $MIN_AF"
echo "MAX_AF:    $MAX_AF"
echo

if [[ ! -s "$GENES_BED" ]]; then
  echo "ERROR: missing genes BED: $GENES_BED" >&2
  exit 1
fi

echo "== BED chromosomes =="
cut -f1 "$GENES_BED" | sort -u
echo

echo "== Check input 1000G VCFs and indexes =="
for c in {1..22}; do
  VCF="${VCF_DIR}/1kGP_high_coverage_Illumina.chr${c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

  if [[ ! -s "$VCF" ]]; then
    echo "WARNING: missing VCF for chr${c}: $VCF" >&2
    continue
  fi

  if [[ ! -s "${VCF}.tbi" ]]; then
    echo "ERROR: missing index for chr${c}: ${VCF}.tbi" >&2
    exit 1
  fi
done

echo "== Extract common biallelic SNPs in gene bodies =="
rm -f "${OUTDIR}"/chr*.genes.af"${MIN_AF}".snps.vcf.gz
rm -f "${OUTDIR}"/chr*.genes.af"${MIN_AF}".snps.vcf.gz.tbi
rm -f "${OUTDIR}/per_chrom_counts.tsv"

printf "CHROM\tN_SNPS\n" > "${OUTDIR}/per_chrom_counts.tsv"

for c in {1..22}; do
  VCF="${VCF_DIR}/1kGP_high_coverage_Illumina.chr${c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  OUT="${OUTDIR}/chr${c}.genes.af${MIN_AF}.snps.vcf.gz"
  LOG="logs/chr${c}.common_snps.af${MIN_AF}.log"

  echo "== chr${c} =="

  if [[ ! -s "$VCF" ]]; then
    echo "Missing VCF; writing count 0"
    printf "chr%s\t0\n" "$c" >> "${OUTDIR}/per_chrom_counts.tsv"
    continue
  fi

  bcftools view -R "$GENES_BED" "$VCF" \
    -m2 -M2 -v snps \
    -i "INFO/AF>=${MIN_AF} && INFO/AF<=${MAX_AF}" \
    -Oz -o "$OUT" 2> "$LOG"

  bcftools index -f "$OUT"

  N=$(bcftools view -H "$OUT" | wc -l)
  printf "chr%s\t%s\n" "$c" "$N" >> "${OUTDIR}/per_chrom_counts.tsv"

  echo "Variants kept: $N"

  if [[ -s "$LOG" ]]; then
    echo "Log:"
    cat "$LOG"
  fi

  echo
done

echo "== Per-chromosome counts =="
cat "${OUTDIR}/per_chrom_counts.tsv"
echo

echo "== Make list of non-empty VCFs to concatenate =="
rm -f "${OUTDIR}/files_to_concat.txt"

for c in {1..22}; do
  f="${OUTDIR}/chr${c}.genes.af${MIN_AF}.snps.vcf.gz"

  if [[ ! -s "$f" ]]; then
    continue
  fi

  n=$(bcftools view -H "$f" | wc -l)

  if [[ "$n" -gt 0 ]]; then
    echo "$f" >> "${OUTDIR}/files_to_concat.txt"
  fi
done

echo "Files to concatenate:"
cat "${OUTDIR}/files_to_concat.txt"
echo

if [[ ! -s "${OUTDIR}/files_to_concat.txt" ]]; then
  echo "ERROR: no non-empty VCFs found. Check chromosome naming and AF field." >&2
  exit 1
fi

echo "== Concatenate =="
FINAL_VCF="results/common_gene_body_snps.af${MIN_AF}.vcf.gz"

bcftools concat \
  -f "${OUTDIR}/files_to_concat.txt" \
  -Oz -o "$FINAL_VCF"

bcftools index -f "$FINAL_VCF"

FINAL_N=$(bcftools view -H "$FINAL_VCF" | wc -l)

echo "Final VCF: $FINAL_VCF"
echo "Final SNP count: $FINAL_N"
echo

echo "== Done =="
echo "Main outputs:"
echo "  ${OUTDIR}/per_chrom_counts.tsv"
echo "  ${OUTDIR}/files_to_concat.txt"
echo "  ${FINAL_VCF}"
echo "  ${FINAL_VCF}.tbi"
