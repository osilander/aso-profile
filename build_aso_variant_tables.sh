#!/usr/bin/env bash
set -euo pipefail

# Build ASO variant target-site tables for:
#   1. ClinVar pathogenic / likely pathogenic missense SNVs
#   2. Common 1000G gene-body SNPs with MAF >= threshold
#   3. ±20 bp stranded sequence windows for both sets
#   4. Unified comparison table
#
# Assumptions:
#   - genes.standard.bed uses chr-prefixed contigs, e.g. chr1
#   - 1000G high-coverage VCFs use chr-prefixed contigs, e.g. chr1
#   - GENCODE/GRCh38 FASTA uses chr-prefixed contigs
#   - ClinVar TSV uses non-chr chromosomes in CHROM column; script adds chr
#
# Required inputs:
#   results/clinvar_20genes_pathogenic.header.tsv
#   genes.standard.bed
#   refs/GRCh38.primary_assembly.genome.fa
#   data/1000G_highcov_GRCh38/
#
# Usage:
#   bash build_aso_variant_tables.sh
#
# Optional:
#   bash build_aso_variant_tables.sh \
#     results/clinvar_20genes_pathogenic.header.tsv \
#     genes.standard.bed \
#     refs/GRCh38.primary_assembly.genome.fa \
#     data/1000G_highcov_GRCh38 \
#     0.2

CLINVAR_TSV="${1:-results/clinvar_20genes_pathogenic.header.tsv}"
GENES_BED="${2:-genes.standard.bed}"
REF_FASTA="${3:-refs/GRCh38.primary_assembly.genome.fa}"
VCF_DIR="${4:-data/1000G_highcov_GRCh38}"
MIN_MAF="${5:-0.2}"

MAX_AF=$(awk -v m="$MIN_MAF" 'BEGIN{print 1-m}')

OUTDIR="results"
COMMON_DIR="${OUTDIR}/common_snps"
LOGDIR="logs"

mkdir -p "$OUTDIR" "$COMMON_DIR" "$LOGDIR"

echo "== Inputs =="
echo "ClinVar TSV: $CLINVAR_TSV"
echo "Genes BED:   $GENES_BED"
echo "Reference:   $REF_FASTA"
echo "1000G dir:   $VCF_DIR"
echo "MIN_MAF:     $MIN_MAF"
echo "AF range:    ${MIN_MAF} to ${MAX_AF}"
echo

for f in "$CLINVAR_TSV" "$GENES_BED" "$REF_FASTA"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing required file: $f" >&2
    exit 1
  fi
done

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

###############################################################################
# PART A: ClinVar pathogenic SNVs and missense subset
###############################################################################

echo "###############################################################################"
echo "PART A: ClinVar pathogenic SNVs and missense subset"
echo "###############################################################################"
echo

echo "== A1: Filter ClinVar pathogenic TSV to SNVs only =="
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1 {print; next}
  length($4)==1 && length($5)==1 {print}
' "$CLINVAR_TSV" \
  > "$OUTDIR/clinvar_pathogenic_snvs.tsv"

N_CLINVAR_SNVS=$(awk 'NR>1' "$OUTDIR/clinvar_pathogenic_snvs.tsv" | wc -l)
echo "ClinVar pathogenic/likely pathogenic SNVs: $N_CLINVAR_SNVS"
echo

echo "== A2: Convert ClinVar SNVs to BED-like coordinates =="
# Input ClinVar TSV columns:
# CHROM POS ID REF ALT GENEINFO CLNSIG CLNREVSTAT MC CLNHGVS CLNDN
#
# Output:
# VAR_CHROM VAR_START0 VAR_END CLINVAR_ID REF ALT GENEINFO CLNSIG CLNREVSTAT MC CLNHGVS CLNDN
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

echo "ClinVar SNV BED-like rows:"
wc -l "$OUTDIR/clinvar_pathogenic_snvs.bedlike.tsv"
echo

echo "== A3: Intersect ClinVar SNVs with genes to recover gene and strand =="
# Variant columns:
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
#
# Gene BED columns:
# 13 GENE_CHROM
# 14 GENE_START
# 15 GENE_END
# 16 GENE
# 17 SCORE
# 18 STRAND
bedtools intersect \
  -a "$OUTDIR/clinvar_pathogenic_snvs.bedlike.tsv" \
  -b "$GENES_BED" \
  -wa -wb \
  > "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv"

N_CLINVAR_WITH_GENE=$(wc -l < "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv")
echo "ClinVar SNVs with gene hit: $N_CLINVAR_WITH_GENE"
echo

echo "Counts by gene, all pathogenic SNVs:"
cut -f16 "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv" \
  | sort | uniq -c | sort -nr
echo

echo "== A4: Make ClinVar missense-only subset =="
awk -F'\t' 'BEGIN{OFS="\t"}
  $10 ~ /missense_variant/ {print}
' "$OUTDIR/clinvar_pathogenic_snvs.with_gene.tsv" \
  > "$OUTDIR/clinvar_pathogenic_snvs.missense.with_gene.tsv"

N_MISSENSE=$(wc -l < "$OUTDIR/clinvar_pathogenic_snvs.missense.with_gene.tsv")
echo "ClinVar pathogenic/likely pathogenic missense SNVs: $N_MISSENSE"
echo

echo "Counts by gene, pathogenic missense:"
cut -f16 "$OUTDIR/clinvar_pathogenic_snvs.missense.with_gene.tsv" \
  | sort | uniq -c | sort -nr
echo

echo "== A5: Make ±20 bp BED windows for all ClinVar pathogenic SNVs =="
awk -F'\t' 'BEGIN{OFS="\t"}
  {
    var_chrom=$1
    var_start=$2
    var_end=$3
    clinvar_id=$4
    ref=$5
    alt=$6
    clnsig=$8
    gene=$16
    strand=$18

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

bedtools getfasta \
  -fi "$REF_FASTA" \
  -bed "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.bed" \
  -name \
  -s \
  -fo "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.stranded.fa"

echo "ClinVar pathogenic SNV FASTA records:"
grep -c '^>' "$OUTDIR/clinvar_pathogenic_snvs.plusminus20.stranded.fa"
echo

echo "== A6: Make ±20 bp BED/FASTA windows for ClinVar missense subset =="
awk -F'\t' 'BEGIN{OFS="\t"}
  {
    var_chrom=$1
    var_start=$2
    var_end=$3
    clinvar_id=$4
    ref=$5
    alt=$6
    clnsig=$8
    gene=$16
    strand=$18

    pos=var_start+1
    win_start=var_start-20
    win_end=var_end+20
    if (win_start < 0) win_start=0

    name=gene"|ClinVarID="clinvar_id"|"var_chrom":"pos":"ref">"alt"|CLNSIG="clnsig"|strand="strand
    print var_chrom,win_start,win_end,name,".",strand
  }
' "$OUTDIR/clinvar_pathogenic_snvs.missense.with_gene.tsv" \
  | sort -k1,1 -k2,2n \
  > "$OUTDIR/clinvar_pathogenic_snvs.missense.plusminus20.bed"

bedtools getfasta \
  -fi "$REF_FASTA" \
  -bed "$OUTDIR/clinvar_pathogenic_snvs.missense.plusminus20.bed" \
  -name \
  -s \
  -fo "$OUTDIR/clinvar_pathogenic_snvs.missense.plusminus20.stranded.fa"

echo "ClinVar missense FASTA records:"
grep -c '^>' "$OUTDIR/clinvar_pathogenic_snvs.missense.plusminus20.stranded.fa"
echo

###############################################################################
# PART B: Common 1000G gene-body SNPs
###############################################################################

echo "###############################################################################"
echo "PART B: Common 1000G gene-body SNPs"
echo "###############################################################################"
echo

echo "== B1: Extract common biallelic SNPs from gene bodies, per chromosome =="
rm -f "${COMMON_DIR}"/chr*.genes.af"${MIN_MAF}".snps.vcf.gz
rm -f "${COMMON_DIR}"/chr*.genes.af"${MIN_MAF}".snps.vcf.gz.tbi
rm -f "${COMMON_DIR}/per_chrom_counts.tsv"

printf "CHROM\tN_SNPS\n" > "${COMMON_DIR}/per_chrom_counts.tsv"

for c in {1..22}; do
  VCF="${VCF_DIR}/1kGP_high_coverage_Illumina.chr${c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  OUT="${COMMON_DIR}/chr${c}.genes.af${MIN_MAF}.snps.vcf.gz"
  LOG="${LOGDIR}/chr${c}.common_snps.af${MIN_MAF}.log"

  echo "== chr${c} =="

  if [[ ! -s "$VCF" ]]; then
    echo "WARNING: missing VCF: $VCF" >&2
    printf "chr%s\t0\n" "$c" >> "${COMMON_DIR}/per_chrom_counts.tsv"
    continue
  fi

  if [[ ! -s "${VCF}.tbi" ]]; then
    echo "ERROR: missing VCF index: ${VCF}.tbi" >&2
    exit 1
  fi

  bcftools view -R "$GENES_BED" "$VCF" \
    -m2 -M2 -v snps \
    -i "INFO/AF>=${MIN_MAF} && INFO/AF<=${MAX_AF}" \
    -Oz -o "$OUT" 2> "$LOG"

  bcftools index -f "$OUT"

  N=$(bcftools view -H "$OUT" | wc -l)
  printf "chr%s\t%s\n" "$c" "$N" >> "${COMMON_DIR}/per_chrom_counts.tsv"
  echo "Variants kept: $N"

  if [[ -s "$LOG" ]]; then
    cat "$LOG"
  fi
  echo
done

echo "Per-chromosome common SNP counts:"
cat "${COMMON_DIR}/per_chrom_counts.tsv"
echo

echo "== B2: Make list of non-empty common-SNP VCFs =="
rm -f "${COMMON_DIR}/files_to_concat.txt"

for c in {1..22}; do
  f="${COMMON_DIR}/chr${c}.genes.af${MIN_MAF}.snps.vcf.gz"
  if [[ ! -s "$f" ]]; then
    continue
  fi

  n=$(bcftools view -H "$f" | wc -l)
  if [[ "$n" -gt 0 ]]; then
    echo "$f" >> "${COMMON_DIR}/files_to_concat.txt"
  fi
done

echo "Files to concatenate:"
cat "${COMMON_DIR}/files_to_concat.txt"
echo

if [[ ! -s "${COMMON_DIR}/files_to_concat.txt" ]]; then
  echo "ERROR: no non-empty common SNP VCFs found." >&2
  exit 1
fi

echo "== B3: Concatenate common SNP VCFs =="
COMMON_VCF="${OUTDIR}/common_gene_body_snps.af${MIN_MAF}.vcf.gz"

bcftools concat \
  -f "${COMMON_DIR}/files_to_concat.txt" \
  -Oz -o "$COMMON_VCF"

bcftools index -f "$COMMON_VCF"

N_COMMON=$(bcftools view -H "$COMMON_VCF" | wc -l)
echo "Common SNP total: $N_COMMON"
echo

echo "== B4: Write common SNP TSV with AF and genotype counts =="
bcftools query \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%INFO/AC_Het\t%INFO/AC_Hom\n' \
  "$COMMON_VCF" \
  > "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.tsv"

printf "CHROM\tPOS\tID\tREF\tALT\tAF\tAC\tAN\tAC_Het\tAC_Hom\n" \
  | cat - "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.tsv" \
  > "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.header.tsv"

echo "Common SNP TSV rows:"
wc -l "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.tsv"
echo

echo "== B5: Intersect common SNPs with genes to recover gene and strand =="
bcftools query \
  -f '%CHROM\t%POS0\t%END\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%INFO/AC_Het\t%INFO/AC_Hom\n' \
  "$COMMON_VCF" \
  > "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.bedlike.tsv"

bedtools intersect \
  -a "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.bedlike.tsv" \
  -b "$GENES_BED" \
  -wa -wb \
  > "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.tsv"

printf "SNP_CHROM\tSNP_START0\tSNP_END\tSNP_ID\tREF\tALT\tAF\tAC\tAN\tAC_Het\tAC_Hom\tGENE_CHROM\tGENE_START\tGENE_END\tGENE\tSCORE\tSTRAND\n" \
  | cat - "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.tsv" \
  > "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.header.tsv"

N_COMMON_WITH_GENE=$(wc -l < "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.tsv")
echo "Common SNPs with gene hit: $N_COMMON_WITH_GENE"
echo

echo "Counts by gene, common SNPs:"
cut -f15 "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.tsv" \
  | sort | uniq -c | sort -nr
echo

echo "== B6: Add MAF, heterozygosity rate, homozygous-alt rate, gene, strand =="
awk -F'\t' 'BEGIN{OFS="\t"}
  {
    af=$7
    ac=$8
    an=$9
    ac_het=$10
    ac_hom=$11

    split(af, a_af, ","); af=a_af[1]
    split(an, a_an, ","); an=a_an[1]
    split(ac_het, a_het, ","); ac_het=a_het[1]
    split(ac_hom, a_hom, ","); ac_hom=a_hom[1]

    n_ind = an / 2
    maf = (af <= 0.5 ? af : 1 - af)
    het_rate = ac_het / n_ind

    # AC_Hom is allele count in homozygous-alt individuals.
    hom_alt_rate = (ac_hom / 2) / n_ind

    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,maf,het_rate,hom_alt_rate,$15,$17
  }
' "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.tsv" \
  > "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.rates.tsv"

printf "SNP_CHROM\tSNP_START0\tSNP_END\tSNP_ID\tREF\tALT\tAF\tAC\tAN\tAC_HET\tAC_HOM\tMAF\tHET_RATE\tHOM_ALT_RATE\tGENE\tSTRAND\n" \
  | cat - "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.rates.tsv" \
  > "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.rates.header.tsv"

echo "Checking common SNP rate-table column counts:"
awk -F'\t' '{print NF}' "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.rates.header.tsv" | sort | uniq -c
echo

echo "== B7: Make ±20 bp BED/FASTA windows for common SNPs =="
awk -F'\t' 'BEGIN{OFS="\t"}
  {
    chrom=$1
    start=$2 - 20
    end=$3 + 20
    snp_id=$4
    ref=$5
    alt=$6
    maf=$12
    het=$13
    gene=$15
    strand=$16

    if (start < 0) start=0

    name=gene"|"snp_id"|"ref">"alt"|MAF="maf"|HET="het"|strand="strand
    print chrom,start,end,name,".",strand
  }
' "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.rates.tsv" \
  > "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.plusminus20.bed"

bedtools getfasta \
  -fi "$REF_FASTA" \
  -bed "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.plusminus20.bed" \
  -name \
  -s \
  -fo "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.plusminus20.stranded.fa"

echo "Common SNP FASTA records:"
grep -c '^>' "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.plusminus20.stranded.fa"
echo

###############################################################################
# PART C: Unified comparison tables
###############################################################################

echo "###############################################################################"
echo "PART C: Unified comparison tables"
echo "###############################################################################"
echo

echo "== C1: Unified pathogenic missense table =="
awk -F'\t' 'BEGIN{
  OFS="\t"
  print "VARIANT_CLASS","GENE","CHROM","POS","START0","END","ID","REF","ALT","CLNSIG","MAF","HET_RATE","STRAND","MC","CLNHGVS","CLNDN"
}
{
  chrom=$1
  start0=$2
  end=$3
  id=$4
  ref=$5
  alt=$6
  clnsig=$8
  mc=$10
  hgvs=$11
  clndn=$12
  gene=$16
  strand=$18
  pos=start0+1

  print "pathogenic_missense",gene,chrom,pos,start0,end,id,ref,alt,clnsig,"NA","NA",strand,mc,hgvs,clndn
}' "$OUTDIR/clinvar_pathogenic_snvs.missense.with_gene.tsv" \
  > "$OUTDIR/unified_pathogenic_missense.tsv"

echo "Column counts:"
awk -F'\t' '{print NF}' "$OUTDIR/unified_pathogenic_missense.tsv" | sort | uniq -c
echo

echo "== C2: Unified common SNP table =="
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1 {
    print "VARIANT_CLASS","GENE","CHROM","POS","START0","END","ID","REF","ALT","CLNSIG","MAF","HET_RATE","STRAND","MC","CLNHGVS","CLNDN"
    next
  }
  {
    chrom=$1
    start0=$2
    end=$3
    id=$4
    ref=$5
    alt=$6
    maf=$12
    het=$13
    gene=$15
    strand=$16
    pos=start0+1

    print "common_snp",gene,chrom,pos,start0,end,id,ref,alt,"NA",maf,het,strand,"NA","NA","NA"
  }
' "$OUTDIR/common_gene_body_snps.af${MIN_MAF}.with_gene.rates.header.tsv" \
  > "$OUTDIR/unified_common_snps.af${MIN_MAF}.tsv"

echo "Column counts:"
awk -F'\t' '{print NF}' "$OUTDIR/unified_common_snps.af${MIN_MAF}.tsv" | sort | uniq -c
echo

echo "== C3: Combine pathogenic missense and common SNPs =="
COMBINED="$OUTDIR/unified_pathogenic_missense_and_common_snps.af${MIN_MAF}.tsv"

head -1 "$OUTDIR/unified_pathogenic_missense.tsv" > "$COMBINED"
tail -n +2 "$OUTDIR/unified_pathogenic_missense.tsv" >> "$COMBINED"
tail -n +2 "$OUTDIR/unified_common_snps.af${MIN_MAF}.tsv" >> "$COMBINED"

echo "Combined table:"
wc -l "$COMBINED"
echo "Column counts:"
awk -F'\t' '{print NF}' "$COMBINED" | sort | uniq -c
echo

echo "== C4: Summary by gene =="
awk -F'\t' '
  NR>1 {
    count[$2 FS $1]++
    genes[$2]=1
  }
  END {
    print "GENE\tpathogenic_missense\tcommon_snp"
    for (g in genes) {
      print g "\t" count[g FS "pathogenic_missense"]+0 "\t" count[g FS "common_snp"]+0
    }
  }
' "$COMBINED" \
  | sort \
  > "$OUTDIR/summary_by_gene.pathogenic_missense_vs_common_snp.af${MIN_MAF}.tsv"

column -t "$OUTDIR/summary_by_gene.pathogenic_missense_vs_common_snp.af${MIN_MAF}.tsv"
echo

###############################################################################
# PART D: optional example-candidate tables
###############################################################################

echo "###############################################################################"
echo "PART D: Example candidate tables"
echo "###############################################################################"
echo

echo "== D1: Top common SNPs per gene by heterozygosity rate =="
# This keeps top 5 common SNPs per gene by HET_RATE.
# It is useful for illustrative examples, not final ASO selection.
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1 {next}
  {
    print $0
  }
' "$OUTDIR/unified_common_snps.af${MIN_MAF}.tsv" \
  | sort -t $'\t' -k2,2 -k12,12nr \
  | awk -F'\t' 'BEGIN{OFS="\t"}
      {
        gene=$2
        n[gene]++
        if (n[gene] <= 5) print
      }
    ' \
  > "$OUTDIR/example_top5_common_snps_per_gene.tsv"

head -1 "$OUTDIR/unified_common_snps.af${MIN_MAF}.tsv" \
  | cat - "$OUTDIR/example_top5_common_snps_per_gene.tsv" \
  > "$OUTDIR/example_top5_common_snps_per_gene.header.tsv"

echo "Example top-5 common SNPs per gene:"
head "$OUTDIR/example_top5_common_snps_per_gene.header.tsv" | column -t
echo

echo "== D2: One pathogenic missense example per gene =="
# Simple example selector:
#   Pathogenic over Likely_pathogenic if sorted text happens to group;
#   this is not a clinical ranking, just a compact example table.
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1 {next}
  {
    print
  }
' "$OUTDIR/unified_pathogenic_missense.tsv" \
  | sort -t $'\t' -k2,2 -k10,10 \
  | awk -F'\t' 'BEGIN{OFS="\t"}
      {
        gene=$2
        if (!(gene in seen)) {
          print
          seen[gene]=1
        }
      }
    ' \
  > "$OUTDIR/example_one_pathogenic_missense_per_gene.tsv"

head -1 "$OUTDIR/unified_pathogenic_missense.tsv" \
  | cat - "$OUTDIR/example_one_pathogenic_missense_per_gene.tsv" \
  > "$OUTDIR/example_one_pathogenic_missense_per_gene.header.tsv"

echo "Example one pathogenic missense per gene:"
head "$OUTDIR/example_one_pathogenic_missense_per_gene.header.tsv" | column -t
echo

###############################################################################
# Final checks
###############################################################################

echo "###############################################################################"
echo "FINAL SUMMARY"
echo "###############################################################################"

echo "ClinVar pathogenic/likely pathogenic SNVs:       $N_CLINVAR_SNVS"
echo "ClinVar pathogenic/likely pathogenic missense:  $N_MISSENSE"
echo "Common gene-body SNPs with MAF >= ${MIN_MAF}:        $N_COMMON"
echo
echo "Main outputs:"
echo "  $OUTDIR/clinvar_pathogenic_snvs.plusminus20.stranded.fa"
echo "  $OUTDIR/clinvar_pathogenic_snvs.missense.plusminus20.stranded.fa"
echo "  $OUTDIR/common_gene_body_snps.af${MIN_MAF}.plusminus20.stranded.fa"
echo "  $OUTDIR/unified_pathogenic_missense.tsv"
echo "  $OUTDIR/unified_common_snps.af${MIN_MAF}.tsv"
echo "  $COMBINED"
echo "  $OUTDIR/summary_by_gene.pathogenic_missense_vs_common_snp.af${MIN_MAF}.tsv"
echo "  $OUTDIR/example_top5_common_snps_per_gene.header.tsv"
echo "  $OUTDIR/example_one_pathogenic_missense_per_gene.header.tsv"
echo
echo "Done."
