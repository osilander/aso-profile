# ASO Handle Profiling for Dominant Disease Genes

This repository contains a first-pass computational workflow for identifying candidate antisense oligonucleotide (ASO) target sites in dominant genetic disease genes.

The central idea is that ASO design should not be limited to the causal pathogenic variant. Once a pathogenic allele is phased, any transcribed heterozygous variant in cis with that allele may serve as an ASO handle. These handles may include the causal SNV itself, linked exonic SNPs, intronic SNPs in pre-mRNA, UTR SNPs, retained intron variants, pseudoexon sequences, or other transcript-specific features.

The immediate workflow:

1. Define a list of dominant disease genes.
2. Extract gene coordinates from GENCODE.
3. Pull ClinVar pathogenic/likely pathogenic variants in those genes.
4. Filter for ASO-relevant SNVs and splice-region variants.
5. Later: compare causal-variant ASO suitability with linked high-MAF SNP ASO suitability.

---

## Conceptual rationale

ASOs are increasingly being used as therapies for dominant gain-of-function and dominant-negative genetic diseases. Most current designs focus on the causal mutation itself, producing N-of-1 or N-of-few therapies.

This is unnecessarily restrictive. ASO efficacy depends on local sequence, accessibility, mismatch position, chemistry, off-target potential, transcript context, and toxicity. A causal SNV may be a poor ASO target, even if it is the disease-causing mutation. In contrast, a pathogenic allele often carries many linked heterozygous variants. If phase is known, these linked variants provide additional candidate handles for allele-selective or haplotype-selective ASO design.

This motivates gene-level ASO-handle libraries: catalogues of targetable high-MAF transcribed variants, ideally with reciprocal reference- and alternate-allele ASOs, that can be selected after phased diagnosis.

---

## Current status

The current workflow has successfully:

- generated a BED file for 20 dominant disease genes from GENCODE v48;
- extracted ClinVar variants overlapping those gene bodies;
- filtered ClinVar variants for pathogenic/likely pathogenic calls;
- produced an ASO-relevant SNV/splice-region candidate table.

Current ASO-relevant ClinVar variant count:

```text
4045 clinvar_20genes_pathogenic_snv_aso_relevant.tsv

```bash
bcftools
htslib/bgzip/tabix
awk
grep
sort
zcat
wget
```

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gtf.gz
```

```bash
gzip -t gencode.v48.basic.annotation.gtf.gz
zcat gencode.v48.basic.annotation.gtf.gz | grep -v '^#' | head
```

```bash
cat > genes.txt <<'EOF'
SCN2A
SCN8A
KCNQ2
KIF1A
SPTLC1
SPTAN1
COL1A1
COL1A2
COL2A1
KCNJ11
ABCC8
PIK3CA
AKT1
MTOR
STAT1
STAT3
NLRP3
CARD11
RELA
KIF5B
EOF
```

```bash
chmod +x make_gene_bed.sh
./make_gene_bed.sh genes.txt gencode.v48.basic.annotation.gtf.gz
```

```bash
sed 's/^chr//' genes.standard.bed > genes.standard.nochr.bed
```

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

bcftools view -R genes.standard.nochr.bed clinvar.vcf.gz \
  -Oz -o clinvar_20genes.vcf.gz

bcftools index -f clinvar_20genes.vcf.gz
```

```bash
(
  bcftools view -h clinvar_20genes.vcf.gz

  bcftools view -H clinvar_20genes.vcf.gz \
    | awk -F'\t' '
      $8 ~ /(^|;)CLNSIG=([^;]*[|,])?(Pathogenic|Likely_pathogenic|Pathogenic\/Likely_pathogenic|Likely_pathogenic\/Pathogenic)([|,;]|$)/
    '
) | bgzip -c > clinvar_20genes_pathogenic.vcf.gz

bcftools index -f clinvar_20genes_pathogenic.vcf.gz
```

```bash
bcftools query \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/GENEINFO\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/MC\t%INFO/CLNHGVS\t%INFO/CLNDN\n' \
  clinvar_20genes_pathogenic.vcf.gz \
  > clinvar_20genes_pathogenic.tsv
printf "CHROM\tPOS\tID\tREF\tALT\tGENEINFO\tCLNSIG\tCLNREVSTAT\tMC\tCLNHGVS\tCLNDN\n" \
  | cat - clinvar_20genes_pathogenic.tsv \
  > clinvar_20genes_pathogenic.header.tsv
```

```bash
awk -F'\t' 'BEGIN{OFS="\t"}
  length($4)==1 && length($5)==1 &&
  $9 ~ /missense_variant|splice_region_variant|splice_donor_variant|splice_acceptor_variant|inframe/
' clinvar_20genes_pathogenic.tsv \
  > clinvar_20genes_pathogenic_snv_aso_relevant.tsv

cut -f6 clinvar_20genes_pathogenic_snv_aso_relevant.tsv \
  | cut -d: -f1 \
  | sort | uniq -c | sort -nr

## result 4045 clinvar_20genes_pathogenic_snv_aso_relevant.tsv
```

### Counts by gene
```bash
715 COL1A1
650 COL2A1
626 COL1A2
449 SCN2A
431 KCNQ2
289 SCN8A
216 ABCC8
155 KIF1A
 96 STAT3
 87 PIK3CA
 75 STAT1
 62 NLRP3
 51 MTOR
 44 KCNJ11
 43 SPTAN1
 26 CARD11
 12 SPTLC1
  6 GLS
  5 RELA
  4 KIF5B
  3 AKT1
```

## Generate sequence windows around pathogenic variants

For each pathogenic SNV, extract local genomic sequence around the variant from GRCh38.

Score causal-variant ASO suitability

For each pathogenic SNV, generate all 16–22 nt windows containing the variant and score:

```bash
ASO length;
GC content;
SNP position within ASO;
homopolymer/repeat content;
local sequence complexity;
predicted Tm;
uniqueness/off-target risk;
predicted RNA accessibility;
exonic/intronic/UTR/splice context.
```

## Pull common SNPs in the same genes

Next, use gnomAD or 1000 Genomes to pull common SNPs across the same gene bodies.

## Score high-MAF SNP ASO suitability

For each high-MAF SNP, design reciprocal reference- and alternate-targeting ASOs where possible.

Score the same features as causal variants.

## Compare causal-variant versus linked-SNP targetability

For each pathogenic variant or gene, ask:

```bash
Is the causal variant itself a good ASO target?
Are there high-MAF linked SNPs with better ASO target scores?
Are the best candidate linked SNPs intronic, exonic, UTR, or splice-associated?
How many candidate ASO handles exist per gene?
How often would a phased patient have at least one usable linked ASO handle?
```

## Evaluate combination therapy candidates

After scoring individual handles, identify cases where two or more linked handles could be used together.

```bash
Potential combinations:

causal variant ASO + linked SNP ASO;
two linked SNP ASOs;
intronic pre-mRNA ASO + mature mRNA ASO;
RNase H degradation ASO + splice/isoform-targeting ASO.
