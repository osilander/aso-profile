#!/usr/bin/env bash
set -euo pipefail

GENES="${1:-genes.txt}"

GTF="${GTF:-data/gencode.v48.basic.annotation.gtf.gz}"
PATHOGENIC="${PATHOGENIC:-results/clinvar_pathogenic_snvs.missense.with_gene.tsv}"
COMMON="${COMMON:-results/common_gene_body_snps.af0.2.with_gene.rates.context.tsv}"
VAR_SUMMARY="${VAR_SUMMARY:-data/variant_summary.txt.gz}"

# Optional dbSNP/1000G VCF with rs IDs in the ID column.
# Example:
#   DBSNP_VCF=data/dbsnp/00-All.vcf.gz ./make_all_gene_variant_plot_inputs.sh
DBSNP_VCF="${DBSNP_VCF:-}"

OUTDIR="gene_plot_inputs"
mkdir -p "$OUTDIR"

for f in "$GENES" "$GTF" "$PATHOGENIC" "$COMMON"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing required file: $f" >&2
    exit 1
  fi
done

echo "Inputs:"
echo "  Genes:        $GENES"
echo "  GTF:          $GTF"
echo "  Pathogenic:   $PATHOGENIC"
echo "  Common SNPs:  $COMMON"
echo "  ClinVar summary: ${VAR_SUMMARY:-none}"
echo "  dbSNP VCF:    ${DBSNP_VCF:-none}"
echo "  Output dir:   $OUTDIR"
echo

###############################################################################
# Step 1: choose MANE Select / Ensembl canonical transcript per gene
###############################################################################

echo "###############################################################################"
echo "Step 1: Select MANE/canonical transcript per gene"
echo "###############################################################################"

python3 - "$GENES" "$GTF" "$OUTDIR/all_genes.transcripts.tsv" <<'PY'
import sys, gzip, re, csv
from pathlib import Path

genes_file = Path(sys.argv[1])
gtf_file = Path(sys.argv[2])
out_file = Path(sys.argv[3])

genes = [x.strip() for x in genes_file.read_text().splitlines() if x.strip()]
gene_set = set(genes)

def attr_value(attr, key):
    m = re.search(rf'{key} "([^"]+)"', attr)
    return m.group(1) if m else ""

def has_tag(attr, tag):
    return f'tag "{tag}"' in attr

records = []

with gzip.open(gtf_file, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue

        chrom, source, feature, start, end, score, strand, frame, attr = fields

        if feature != "transcript":
            continue

        gene_name = attr_value(attr, "gene_name")
        if gene_name not in gene_set:
            continue

        transcript_id = attr_value(attr, "transcript_id")
        transcript_name = attr_value(attr, "transcript_name")
        transcript_type = attr_value(attr, "transcript_type")
        gene_id = attr_value(attr, "gene_id")
        gene_type = attr_value(attr, "gene_type")

        start = int(start)
        end = int(end)
        span = end - start + 1

        mane = has_tag(attr, "MANE_Select")
        canonical = has_tag(attr, "Ensembl_canonical")
        appris = "appris_principal" in attr
        basic = has_tag(attr, "basic")

        # Priority:
        # 1. MANE_Select
        # 2. Ensembl_canonical
        # 3. APPRIS principal
        # 4. protein_coding
        # 5. basic
        # 6. longest transcript span
        priority = (
            1 if mane else 0,
            1 if canonical else 0,
            1 if appris else 0,
            1 if transcript_type == "protein_coding" else 0,
            1 if basic else 0,
            span
        )

        records.append({
            "GENE": gene_name,
            "CHROM": chrom,
            "TRANSCRIPT_START0": start - 1,
            "TRANSCRIPT_END": end,
            "STRAND": strand,
            "GENE_ID": gene_id,
            "GENE_TYPE": gene_type,
            "TRANSCRIPT_ID": transcript_id,
            "TRANSCRIPT_NAME": transcript_name,
            "TRANSCRIPT_TYPE": transcript_type,
            "TRANSCRIPT_SPAN_BP": span,
            "MANE_SELECT": "yes" if mane else "no",
            "ENSEMBL_CANONICAL": "yes" if canonical else "no",
            "APPRIS": "yes" if appris else "no",
            "BASIC": "yes" if basic else "no",
            "priority": priority
        })

best = {}
for r in records:
    g = r["GENE"]
    if g not in best or r["priority"] > best[g]["priority"]:
        best[g] = r

fields = [
    "GENE",
    "CHROM",
    "TRANSCRIPT_START0",
    "TRANSCRIPT_END",
    "STRAND",
    "GENE_ID",
    "GENE_TYPE",
    "TRANSCRIPT_ID",
    "TRANSCRIPT_NAME",
    "TRANSCRIPT_TYPE",
    "TRANSCRIPT_SPAN_BP",
    "MANE_SELECT",
    "ENSEMBL_CANONICAL",
    "APPRIS",
    "BASIC"
]

with open(out_file, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    w.writeheader()

    for g in genes:
        if g not in best:
            w.writerow({
                "GENE": g,
                "CHROM": "NA",
                "TRANSCRIPT_START0": "NA",
                "TRANSCRIPT_END": "NA",
                "STRAND": "NA",
                "GENE_ID": "NA",
                "GENE_TYPE": "NA",
                "TRANSCRIPT_ID": "NA",
                "TRANSCRIPT_NAME": "NA",
                "TRANSCRIPT_TYPE": "NA",
                "TRANSCRIPT_SPAN_BP": "NA",
                "MANE_SELECT": "no",
                "ENSEMBL_CANONICAL": "no",
                "APPRIS": "no",
                "BASIC": "no"
            })
        else:
            w.writerow({k: best[g][k] for k in fields})

print(f"Wrote {out_file}")
PY

column -t "$OUTDIR/all_genes.transcripts.tsv"
echo

###############################################################################
# Step 2: extract exon/CDS intervals from the selected transcripts
###############################################################################

echo "###############################################################################"
echo "Step 2: Extract selected-transcript exon/CDS intervals"
echo "###############################################################################"

python3 - "$OUTDIR/all_genes.transcripts.tsv" "$GTF" "$OUTDIR/all_genes.exons_for_plot.tsv" "$OUTDIR/all_genes.cds_for_plot.tsv" <<'PY'
import sys, gzip, re, csv
from pathlib import Path

tx_file = Path(sys.argv[1])
gtf_file = Path(sys.argv[2])
exon_out = Path(sys.argv[3])
cds_out = Path(sys.argv[4])

def attr_value(attr, key):
    m = re.search(rf'{key} "([^"]+)"', attr)
    return m.group(1) if m else None

tx = {}
tx_ids = set()

with open(tx_file) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row["TRANSCRIPT_ID"] == "NA":
            continue
        tx[row["TRANSCRIPT_ID"]] = row
        tx_ids.add(row["TRANSCRIPT_ID"])

def tx_interval_coords(s0, e0, tx_start0, tx_end, strand):
    """
    Convert genomic half-open interval [s0, e0) to transcript-oriented coordinates
    from 5' end of selected transcript.
    """
    if strand == "+":
        a = s0 - tx_start0
        b = e0 - tx_start0
    else:
        a = tx_end - e0
        b = tx_end - s0
    return min(a, b), max(a, b)

exon_rows = []
cds_rows = []

with gzip.open(gtf_file, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue

        chrom, source, feature, start, end, score, strand, frame, attr = fields

        if feature not in {"exon", "CDS"}:
            continue

        transcript_id = attr_value(attr, "transcript_id")
        if transcript_id not in tx_ids:
            continue

        row = tx[transcript_id]
        gene = row["GENE"]
        tx_start0 = int(row["TRANSCRIPT_START0"])
        tx_end = int(row["TRANSCRIPT_END"])
        tx_strand = row["STRAND"]

        s0 = int(start) - 1
        e0 = int(end)

        txs, txe = tx_interval_coords(s0, e0, tx_start0, tx_end, tx_strand)

        out_row = {
            "GENE": gene,
            "TRANSCRIPT_ID": transcript_id,
            "FEATURE": feature,
            "TX_START": txs,
            "TX_END": txe,
            "GENOMIC_CHROM": chrom,
            "GENOMIC_START0": s0,
            "GENOMIC_END": e0,
            "STRAND": tx_strand,
            "LENGTH": txe - txs
        }

        if feature == "exon":
            exon_rows.append(out_row)
        else:
            cds_rows.append(out_row)

fields = [
    "GENE",
    "TRANSCRIPT_ID",
    "FEATURE",
    "TX_START",
    "TX_END",
    "GENOMIC_CHROM",
    "GENOMIC_START0",
    "GENOMIC_END",
    "STRAND",
    "LENGTH"
]

for rows, outfile in [(exon_rows, exon_out), (cds_rows, cds_out)]:
    rows = sorted(rows, key=lambda x: (x["GENE"], int(x["TX_START"]), int(x["TX_END"])))

    with open(outfile, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote {outfile}: {len(rows)} rows")
PY

head "$OUTDIR/all_genes.exons_for_plot.tsv" | column -t
echo

###############################################################################
# Step 3: optional rsID lookup from dbSNP/1000G VCF
###############################################################################

echo "###############################################################################"
echo "Step 3: Optional rsID lookup"
echo "###############################################################################"

RS_LOOKUP="$OUTDIR/rs_lookup.tsv"

if [[ -n "$DBSNP_VCF" && -s "$DBSNP_VCF" ]]; then
  echo "Using VCF for rsID lookup: $DBSNP_VCF"

  awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $2!="NA" {
    print $2,$3,$4
  }' "$OUTDIR/all_genes.transcripts.tsv" \
    > "$OUTDIR/transcript_regions.bed"

  bcftools query \
    -R "$OUTDIR/transcript_regions.bed" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
    "$DBSNP_VCF" \
    | awk -F'\t' 'BEGIN{OFS="\t"}
        $5!="." && $5!="NA" {
          split($4, alts, ",")
          for (i in alts) {
            print $1":"$2":"$3":"alts[i], $5
          }
        }
      ' \
    > "$RS_LOOKUP"

  echo "rs lookup rows:"
  wc -l "$RS_LOOKUP"
else
  echo "No dbSNP/1000G VCF supplied. RS_ID will be NA unless present upstream."
  : > "$RS_LOOKUP"
fi

echo

###############################################################################
# Step 4: build variant tables
###############################################################################

echo "###############################################################################"
echo "Step 4: Build variant tables"
echo "###############################################################################"

python3 - "$OUTDIR/all_genes.transcripts.tsv" \
           "$OUTDIR/all_genes.exons_for_plot.tsv" \
           "$OUTDIR/all_genes.cds_for_plot.tsv" \
           "$PATHOGENIC" \
           "$COMMON" \
           "$VAR_SUMMARY" \
           "$RS_LOOKUP" \
           "$OUTDIR" <<'PY'
import sys, csv, gzip, re, math
from pathlib import Path

tx_file = Path(sys.argv[1])
exon_file = Path(sys.argv[2])
cds_file = Path(sys.argv[3])
pathogenic_file = Path(sys.argv[4])
common_file = Path(sys.argv[5])
variant_summary_file = Path(sys.argv[6])
rs_lookup_file = Path(sys.argv[7])
outdir = Path(sys.argv[8])

###############################################################################
# Helpers
###############################################################################

def to_int(x):
    try:
        return int(x)
    except Exception:
        return None

def to_float(x):
    try:
        if x in {"", "NA", None}:
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]

def pos_tx_coord(pos0, tx_start0, tx_end, strand):
    """
    Single-base coordinate from 5' end of selected transcript.
    For plus strand, first base of transcript is 0.
    For minus strand, last genomic base of transcript is 0.
    """
    if strand == "+":
        return pos0 - tx_start0
    else:
        return tx_end - (pos0 + 1)

def in_any_interval(pos0, intervals):
    return any(s <= pos0 < e for s, e in intervals)

def extract_hgvs_bits(name):
    cdna = "NA"
    protein = "NA"
    nuc = "NA"

    m = re.search(r":(c\.[A-Za-z0-9_\+\-\*?>=<\.]+(?:delins|del|dup|ins)?[A-Za-z0-9_\+\-\*?>=<\.]*)", name)
    if m:
        cdna = m.group(1)

    m = re.search(r":(g\.[A-Za-z0-9_\+\-\*?>=<\.]+(?:delins|del|dup|ins)?[A-Za-z0-9_\+\-\*?>=<\.]*)", name)
    if m:
        nuc = m.group(1)

    m = re.search(r"\((p\.[^)]+)\)", name)
    if m:
        protein = m.group(1)

    return nuc, cdna, protein

###############################################################################
# Load transcript models
###############################################################################

tx_by_gene = {}
genes = []

with open(tx_file) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        gene = row["GENE"]
        genes.append(gene)
        if row["TRANSCRIPT_ID"] == "NA":
            continue
        row["TRANSCRIPT_START0"] = int(row["TRANSCRIPT_START0"])
        row["TRANSCRIPT_END"] = int(row["TRANSCRIPT_END"])
        row["TRANSCRIPT_SPAN_BP"] = int(row["TRANSCRIPT_SPAN_BP"])
        tx_by_gene[gene] = row

exons_by_gene = {g: [] for g in genes}
cds_by_gene = {g: [] for g in genes}

with open(exon_file) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        g = row["GENE"]
        exons_by_gene.setdefault(g, []).append((int(row["GENOMIC_START0"]), int(row["GENOMIC_END"])))

with open(cds_file) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        g = row["GENE"]
        cds_by_gene.setdefault(g, []).append((int(row["GENOMIC_START0"]), int(row["GENOMIC_END"])))

for g in genes:
    exons_by_gene[g] = merge_intervals(exons_by_gene.get(g, []))
    cds_by_gene[g] = merge_intervals(cds_by_gene.get(g, []))

###############################################################################
# ClinVar summary lookup
###############################################################################

clinvar_name = {}
clinvar_review = {}
clinvar_type = {}

if variant_summary_file.exists() and variant_summary_file.stat().st_size > 0:
    with gzip.open(variant_summary_file, "rt") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            vid = row.get("VariationID", "")
            if not vid:
                continue
            clinvar_name[vid] = row.get("Name", "")
            clinvar_review[vid] = row.get("ReviewStatus", "")
            clinvar_type[vid] = row.get("Type", "")

###############################################################################
# rs lookup
###############################################################################

rs = {}
if rs_lookup_file.exists() and rs_lookup_file.stat().st_size > 0:
    with open(rs_lookup_file) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                rs[fields[0]] = fields[1]

def rs_id(chrom, pos1, ref, alt):
    return rs.get(f"{chrom}:{pos1}:{ref}:{alt}", "NA")

###############################################################################
# Build pathogenic table
###############################################################################

path_rows_all = []

with open(pathogenic_file) as f:
    r = csv.reader(f, delimiter="\t")

    for row in r:
        if len(row) < 18:
            continue

        gene = row[15]
        if gene not in tx_by_gene:
            continue

        tx = tx_by_gene[gene]
        chrom = tx["CHROM"]
        tx_start0 = tx["TRANSCRIPT_START0"]
        tx_end = tx["TRANSCRIPT_END"]
        strand = tx["STRAND"]
        tx_len = tx["TRANSCRIPT_SPAN_BP"]

        var_chrom = row[0]
        start0 = int(row[1])
        end = int(row[2])
        clinvar_id = row[3]
        ref = row[4]
        alt = row[5]
        clnsig = row[7]
        mc = row[9]
        genomic_hgvs = row[10]
        phenotype = row[11]

        if var_chrom != chrom:
            continue

        inside_tx = tx_start0 <= start0 < tx_end

        if inside_tx:
            tx_pos = pos_tx_coord(start0, tx_start0, tx_end, strand)
            rel_pos = tx_pos / tx_len if tx_len else "NA"

            if in_any_interval(start0, cds_by_gene.get(gene, [])):
                context = "CDS"
            elif in_any_interval(start0, exons_by_gene.get(gene, [])):
                context = "exon_nonCDS"
            else:
                context = "intron"
        else:
            tx_pos = "NA"
            rel_pos = "NA"
            context = "outside_selected_transcript"

        name = clinvar_name.get(clinvar_id, "")
        nuc_change, cdna_change, aa_change = extract_hgvs_bits(name)

        label = aa_change
        if label == "NA":
            label = cdna_change
        if label == "NA":
            label = genomic_hgvs
        if label == "NA" or label == "":
            label = clinvar_id

        path_rows_all.append({
            "GENE": gene,
            "VARIANT_CLASS": "pathogenic_missense",
            "CHROM": var_chrom,
            "POS": start0 + 1,
            "START0": start0,
            "END": end,
            "TX_POS_FROM_5P": tx_pos,
            "REL_POS_FROM_5P": rel_pos,
            "REF": ref,
            "ALT": alt,
            "ID": clinvar_id,
            "RS_ID": "NA",
            "LABEL": label,
            "CLNSIG": clnsig,
            "CONTEXT": context,
            "MAF": "NA",
            "HET_RATE": "NA",
            "GENOMIC_HGVS": genomic_hgvs,
            "NUC_CHANGE_FROM_CLINVAR_NAME": nuc_change,
            "CDNA_CHANGE": cdna_change,
            "AA_CHANGE": aa_change,
            "MOLECULAR_CONSEQUENCE": mc,
            "PHENOTYPE": phenotype,
            "CLINVAR_NAME": name if name else "NA",
            "CLINVAR_REVIEW_STATUS": clinvar_review.get(clinvar_id, "NA"),
            "TRANSCRIPT_ID": tx["TRANSCRIPT_ID"],
            "TRANSCRIPT_NAME": tx["TRANSCRIPT_NAME"],
            "TRANSCRIPT_SPAN_BP": tx_len,
            "STRAND": strand,
            "PLOT_ELIGIBLE": 1 if inside_tx else 0
        })

###############################################################################
# Build common SNP table
###############################################################################

common_rows_all = []

with open(common_file) as f:
    r = csv.reader(f, delimiter="\t")

    for row in r:
        if len(row) < 17:
            continue

        gene = row[14]
        if gene not in tx_by_gene:
            continue

        tx = tx_by_gene[gene]
        chrom = tx["CHROM"]
        tx_start0 = tx["TRANSCRIPT_START0"]
        tx_end = tx["TRANSCRIPT_END"]
        strand = tx["STRAND"]
        tx_len = tx["TRANSCRIPT_SPAN_BP"]

        snp_chrom = row[0]
        start0 = int(row[1])
        end = int(row[2])
        snp_id = row[3]
        ref = row[4]
        alt = row[5]
        af = row[6]
        maf = row[11]
        het_rate = row[12]
        old_context = row[16]

        if snp_chrom != chrom:
            continue

        inside_tx = tx_start0 <= start0 < tx_end

        if inside_tx:
            tx_pos = pos_tx_coord(start0, tx_start0, tx_end, strand)
            rel_pos = tx_pos / tx_len if tx_len else "NA"

            if in_any_interval(start0, cds_by_gene.get(gene, [])):
                context = "CDS"
            elif in_any_interval(start0, exons_by_gene.get(gene, [])):
                context = "exon_nonCDS"
            else:
                context = "intron"
        else:
            tx_pos = "NA"
            rel_pos = "NA"
            context = "outside_selected_transcript"

        rid = rs_id(snp_chrom, start0 + 1, ref, alt)
        label = rid if rid != "NA" else snp_id

        common_rows_all.append({
            "GENE": gene,
            "VARIANT_CLASS": "common_snp",
            "CHROM": snp_chrom,
            "POS": start0 + 1,
            "START0": start0,
            "END": end,
            "TX_POS_FROM_5P": tx_pos,
            "REL_POS_FROM_5P": rel_pos,
            "REF": ref,
            "ALT": alt,
            "ID": snp_id,
            "RS_ID": rid,
            "LABEL": label,
            "CLNSIG": "NA",
            "CONTEXT": context,
            "CONTEXT_OLD_GENE_BODY_UNION": old_context,
            "AF": af,
            "MAF": maf,
            "HET_RATE": het_rate,
            "GENOMIC_HGVS": "NA",
            "NUC_CHANGE_FROM_CLINVAR_NAME": "NA",
            "CDNA_CHANGE": "NA",
            "AA_CHANGE": "NA",
            "MOLECULAR_CONSEQUENCE": "NA",
            "PHENOTYPE": "NA",
            "CLINVAR_NAME": "NA",
            "CLINVAR_REVIEW_STATUS": "NA",
            "TRANSCRIPT_ID": tx["TRANSCRIPT_ID"],
            "TRANSCRIPT_NAME": tx["TRANSCRIPT_NAME"],
            "TRANSCRIPT_SPAN_BP": tx_len,
            "STRAND": strand,
            "PLOT_ELIGIBLE": 1 if inside_tx else 0
        })

###############################################################################
# Write outputs
###############################################################################

path_fields = [
    "GENE","VARIANT_CLASS","CHROM","POS","START0","END",
    "TX_POS_FROM_5P","REL_POS_FROM_5P",
    "REF","ALT","ID","RS_ID","LABEL","CLNSIG","CONTEXT",
    "GENOMIC_HGVS","NUC_CHANGE_FROM_CLINVAR_NAME","CDNA_CHANGE","AA_CHANGE",
    "MOLECULAR_CONSEQUENCE","PHENOTYPE","CLINVAR_NAME","CLINVAR_REVIEW_STATUS",
    "TRANSCRIPT_ID","TRANSCRIPT_NAME","TRANSCRIPT_SPAN_BP","STRAND","PLOT_ELIGIBLE"
]

common_fields = [
    "GENE","VARIANT_CLASS","CHROM","POS","START0","END",
    "TX_POS_FROM_5P","REL_POS_FROM_5P",
    "REF","ALT","ID","RS_ID","LABEL","CONTEXT",
    "CONTEXT_OLD_GENE_BODY_UNION","AF","MAF","HET_RATE",
    "TRANSCRIPT_ID","TRANSCRIPT_NAME","TRANSCRIPT_SPAN_BP","STRAND","PLOT_ELIGIBLE"
]

plot_fields = [
    "GENE","VARIANT_CLASS","CHROM","POS","START0","END",
    "TX_POS_FROM_5P","REL_POS_FROM_5P",
    "REF","ALT","ID","RS_ID","LABEL",
    "CLNSIG","CONTEXT","MAF","HET_RATE",
    "GENOMIC_HGVS","CDNA_CHANGE","AA_CHANGE",
    "MOLECULAR_CONSEQUENCE","PHENOTYPE",
    "TRANSCRIPT_ID","TRANSCRIPT_NAME","TRANSCRIPT_SPAN_BP","STRAND"
]

def numeric_sort_key(r):
    txp = r.get("TX_POS_FROM_5P", "NA")
    try:
        txp = float(txp)
    except Exception:
        txp = 1e30
    return (r["GENE"], txp, r["VARIANT_CLASS"], str(r["ID"]))

path_rows_all = sorted(path_rows_all, key=numeric_sort_key)
common_rows_all = sorted(common_rows_all, key=numeric_sort_key)

with open(outdir / "all_genes.pathogenic_table.tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=path_fields, delimiter="\t")
    w.writeheader()
    for r in path_rows_all:
        w.writerow({k: r.get(k, "NA") for k in path_fields})

with open(outdir / "all_genes.common_snps_table.tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=common_fields, delimiter="\t")
    w.writeheader()
    for r in common_rows_all:
        w.writerow({k: r.get(k, "NA") for k in common_fields})

plot_rows = []

for r in path_rows_all:
    if r["PLOT_ELIGIBLE"] == 1:
        plot_rows.append({k: r.get(k, "NA") for k in plot_fields})

for r in common_rows_all:
    if r["PLOT_ELIGIBLE"] == 1:
        plot_rows.append({k: r.get(k, "NA") for k in plot_fields})

plot_rows = sorted(plot_rows, key=numeric_sort_key)

with open(outdir / "all_genes.variants_for_plot.tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=plot_fields, delimiter="\t")
    w.writeheader()
    for r in plot_rows:
        w.writerow(r)

###############################################################################
# Summary
###############################################################################

summary_fields = [
    "GENE",
    "TRANSCRIPT_ID",
    "TRANSCRIPT_NAME",
    "TRANSCRIPT_SPAN_BP",
    "STRAND",
    "PATHOGENIC_TOTAL",
    "PATHOGENIC_IN_SELECTED_TRANSCRIPT",
    "COMMON_TOTAL_FROM_GENE_BODY_TABLE",
    "COMMON_IN_SELECTED_TRANSCRIPT",
    "COMMON_INTRON",
    "COMMON_EXON_NONCDS",
    "COMMON_CDS",
    "COMMON_OUTSIDE_SELECTED_TRANSCRIPT",
    "COMMON_WITH_RSID"
]

with open(outdir / "all_genes.summary.tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=summary_fields, delimiter="\t")
    w.writeheader()

    for gene in genes:
        tx = tx_by_gene.get(gene, {})

        p_all = [r for r in path_rows_all if r["GENE"] == gene]
        c_all = [r for r in common_rows_all if r["GENE"] == gene]
        c_in = [r for r in c_all if r["PLOT_ELIGIBLE"] == 1]

        row = {
            "GENE": gene,
            "TRANSCRIPT_ID": tx.get("TRANSCRIPT_ID", "NA"),
            "TRANSCRIPT_NAME": tx.get("TRANSCRIPT_NAME", "NA"),
            "TRANSCRIPT_SPAN_BP": tx.get("TRANSCRIPT_SPAN_BP", "NA"),
            "STRAND": tx.get("STRAND", "NA"),
            "PATHOGENIC_TOTAL": len(p_all),
            "PATHOGENIC_IN_SELECTED_TRANSCRIPT": sum(1 for r in p_all if r["PLOT_ELIGIBLE"] == 1),
            "COMMON_TOTAL_FROM_GENE_BODY_TABLE": len(c_all),
            "COMMON_IN_SELECTED_TRANSCRIPT": len(c_in),
            "COMMON_INTRON": sum(1 for r in c_in if r["CONTEXT"] == "intron"),
            "COMMON_EXON_NONCDS": sum(1 for r in c_in if r["CONTEXT"] == "exon_nonCDS"),
            "COMMON_CDS": sum(1 for r in c_in if r["CONTEXT"] == "CDS"),
            "COMMON_OUTSIDE_SELECTED_TRANSCRIPT": sum(1 for r in c_all if r["PLOT_ELIGIBLE"] == 0),
            "COMMON_WITH_RSID": sum(1 for r in c_all if r.get("RS_ID", "NA") != "NA")
        }
        w.writerow(row)

print(f"Wrote outputs to {outdir}")
print(f"Pathogenic rows total: {len(path_rows_all)}")
print(f"Common SNP rows total: {len(common_rows_all)}")
print(f"Plot rows: {len(plot_rows)}")
PY

echo
echo "Summary:"
column -t "$OUTDIR/all_genes.summary.tsv"
echo

echo "Output files:"
ls -lh "$OUTDIR"
echo

echo "Plot table preview:"
head "$OUTDIR/all_genes.variants_for_plot.tsv" | column -t
