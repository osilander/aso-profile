#!/usr/bin/env bash
set -euo pipefail

# Make ±20 bp stranded FASTA windows around ClinVar pathogenic SNVs.
#
# Inputs:
#   1. ClinVar pathogenic TSV with header
#   2. chr-prefixed gene BED
#   3. chr-prefixed reference FASTA
#
# Expected ClinVar TSV columns:
#   CHROM POS ID REF ALT GENEINFO CLNSIG CLNREVSTAT MC CLNHGVS CLNDN
#
# Example usage:
#   bash make_clinvar_pathogenic_snv_windows.sh \
#     results/clinvar_20genes_pathogenic.header.tsv \
#     genes.standard.bed \
#     refs/GRCh38.primary_assembly.genome.fa
#
# Outputs:
#   results/clinvar_pathogenic_snvs.tsv
#   results/clinvar_pathogenic_snvs.bedlike.tsv
#   results/clinvar_pathogenic_snvs.with_gene.tsv
#   results/clinvar_pathogenic_snvs.plusminus20.bed
#   results/clinvar_pathogenic_snvs.plusminus20.stranded.fa

CLINVAR_TSV="${1:-results/clinvar_20genes_pathogenic.header.tsv}"
GENES_BED="${2:-genes.standard.bed}"
REF_FASTA="${3:-refs/GRCh38.primary_assembly.genome.fa}"

OUTDIR="results"
mkdir -p "$OUTDIR"
mkdir -p logs

echo "== Inputs =="
echo "ClinVar TSV: $CLINVAR_TSV"
echo "Genes BED:   $GENES_BED"
echo "Reference:   $REF_FASTA"
echo

if [[ ! -s "$CLINVAR_TSV" ]]; then
  echo "ERROR: ClinVar TSV missing or empty: $CLINVAR_TSV" >&2
  exit 1
fi

if [[ ! -s "$GENES_BED" ]]; then
  echo "ERROR: gene BED missing or empty: $GENES_BED" >&2
  exit 1
fi

if [[ ! -s "$REF_FASTA" ]]; then
  echo "ERROR: reference FASTA missing or empty: $REF_FASTA" >&2
  exit 1
fi

echo "== Check FASTA index =="
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo "Indexing FASTA..."
  samtools faidx "$REF_FASTA"
fi

echo "Reference contigs:"
cut -f1 "${REF_FASTA}.fai" | head
echo

echo "Gene BED contigs:"
cut -f1 "$GENES_BED" | sort -u | head
echo

echo "ClinVar chromosomes from TSV:"
awk -F'\t' 'NR>1 {print $1}' "$CLINVAR_TSV" | sort -u | head
echo

echo "== Step 1: filter ClinVar TSV to SNVs only =="
# Keep header + rows where REF and ALT are both single bases.
# This keeps all pathogenic/likely pathogenic SNVs from the TSV.
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1 {
    print
    next
  }
  length($4)==1 && length($5)==1 {
    print
  }
' "$CLINVAR_TSV" \
  > "$OUTDIR/clinvar_pathogenic_snvs.tsv"

echo "ClinVar pathogenic SNVs:"
awk 'NR>1' "$OUTDIR/clinvar_pathogenic_snvs.tsv" | wc -l
echo

echo "Preview:"
head "$OUTDIR/clinvar_pathogenic_snvs.tsv" | column -t
echo

echo "== Step 2: convert ClinVar SNVs to BED-like coordinates =="
# ClinVar TSV chromosomes are non-chr, e.g. 1, 2, 11.
# GENCODE BED/FASTA are chr-prefixed, e.g. chr1, chr2, chr11.
#
# ClinVar POS is 1-based.
# BED start = POS - 1
# BED end = POS
#
# BED-like output columns:
#   CHROM START0 END ID REF ALT GENEINFO CLNSIG CLNREVSTAT MC CLNHGVS CLNDN
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1 {next}
  {
    chrom=$1
    if (chrom !~ /^chr/) chrom="chr" chrom

    pos=$2
    start=pos-1
    end=pos

    print chrom,start,end,$3,$4,$5,$6,$7,$8,$9,$10,$11
  }
' "$OUTDIR/clinvar_pathogenic_snvs.tsv" \
  > "$OUTDIR/clinvar_pathogenic_snvs.bedlike.tsv"

echo "BED-like ClinVar SNVs:"
wc -l "$OUTDIR/clinvar_pathogenic_snvs.bedlike.tsv"
head "$OUTDIR/clinvar_pathogenic_snvs.bedlike.tsv" | column -t
echo

echo "== Step 3: intersect ClinVar SNVs with gene BED to recover gene and strand =="
# -wa -wb keeps both variant columns and gene BED columns.
#
# Variant columns:
#   1  VAR_CHROM
#   2  VAR_START0
#   3  VAR_END
#   4  CLINVAR_ID
#   5  REF
#   6  ALT
#   7  GENEINFO
#   8  CLNSIG
#   9  CLNREVSTAT
#   10 MC
#   11 CLNHGVS
#   12 CLNDN
#
# Gene BED columns:
#   13 GENE_CHROM
#   14 GENE_START
#   15 GENE_END
#   16 GENE
#   17 SCORE
#   18 STRAND
bedtools intersect \
  -a "$OUTDIR/clinvar_pathogenic_snvs.bedlike.tsv" \
  -b "$GENES_BED" \
  -wa -wb \
  > "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv"

echo "ClinVar SNVs intersected to genes:"
wc -l "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv"
head "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv" | column -t
echo

echo "Counts by intersected gene:"
cut -f16 "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv" \
  | sort | uniq -c | sort -nr
echo

echo "== Step 4: make ±20 bp BED windows =="
# For a variant at BED start0/end:
#   variant base is start0..end
#   ±20 window is start0-20 .. end+20
# Total length = 41 bp.
#
# BED name contains:
#   gene|ClinVarID|chrom:pos:REF>ALT|CLNSIG=...|strand=...
#
# Use gene strand in column 6 so bedtools getfasta -s returns transcript-oriented sequence.
awk -F'\t' 'BEGIN{OFS="\t"}
  {
    var_chrom=$1
    var_start=$2
    var_end=$3
    clinvar_id=$4
    ref=$5
    alt=$6
    geneinfo=$7
    clnsig=$8
    revstat=$9
    mc=$10
    hgvs=$11
    disease=$12

    gene=$16
    strand=$18

    # Convert BED start back to 1-based position for readable name
    pos=var_start+1

    win_start=var_start-20
    win_end=var_end+20
    if (win_start < 0) win_start=0

    name=gene"|ClinVarID="clinvar_id"|"var_chrom":"pos":"ref">"alt"|CLNSIG="clnsig"|strand="strand

    print var_chrom,win_start,win_end,name,".",strand
  }
' "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv" \
  | sort -k1,1 -k2,2n \
  > "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.bed"

echo "ClinVar ±20 BED windows:"
wc -l "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.bed"
head "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.bed" | column -t
echo

echo "== Step 5: extract stranded FASTA =="
bedtools getfasta \
  -fi "$REF_FASTA" \
  -bed "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.bed" \
  -name \
  -s \
  -fo "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.stranded.fa"

echo "FASTA records:"
grep -c '^>' "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.stranded.fa"
echo

echo "First FASTA records:"
head "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.stranded.fa"
echo

echo "== Step 6: sanity checks =="
N_SNV=$(awk 'NR>1' "$OUTDIR/clinvar_pathogenic_snvs.tsv" | wc -l)
N_INTERSECT=$(wc -l < "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv")
N_BED=$(wc -l < "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.bed")
N_FASTA=$(grep -c '^>' "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.stranded.fa")

echo "SNVs in TSV:          $N_SNV"
echo "SNVs with gene hit:   $N_INTERSECT"
echo "BED windows:          $N_BED"
echo "FASTA records:        $N_FASTA"

if [[ "$N_INTERSECT" -ne "$N_BED" || "$N_BED" -ne "$N_FASTA" ]]; then
  echo "WARNING: counts differ. This may indicate variants overlapping multiple genes or missing FASTA contigs." >&2
fi

echo
echo "== Done =="
echo "Outputs:"
echo "  $OUTDIR/clinvar_pathogenic_snvs.tsv"
echo "  $OUTDIR/clinvar_pathogenic_snvs.bedlike.tsv"
echo "  $OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv"
echo "  $OUTDIR/clinvar_pathogenic_snvs.plusminus20.bed"
echo "  $OUTDIR/clinvar_pathogenic_snvs.plusminus20.stranded.fa"
