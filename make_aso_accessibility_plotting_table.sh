#!/usr/bin/env bash
set -euo pipefail

# Make local accessibility and final R plotting table.
#
# Starting inputs from previous scripts:
#   results/rnafold/rnafold_plusminus100_mfe.tsv
#   results/aso_17mer_scores_with_rnafold.tsv
#
# Outputs:
#   results/rnafold/rnafold_plusminus100_local_accessibility.tsv
#   results/aso_17mer_scores_with_rnafold_accessibility.tsv
#   results/aso_17mer_plotting_table.tsv
#   results/top_common_snp_per_gene_by_composite_score.tsv
#
# The final table used in R is:
#   results/aso_17mer_plotting_table.tsv

RNAFOLD_TSV="results/rnafold/rnafold_plusminus100_mfe.tsv"
SCORE_TSV="results/aso_17mer_scores_with_rnafold.tsv"

ACCESS_TSV="results/rnafold/rnafold_plusminus100_local_accessibility.tsv"
JOINED_TSV="results/aso_17mer_scores_with_rnafold_accessibility.tsv"
PLOT_TSV="results/aso_17mer_plotting_table.tsv"
TOP_COMMON_TSV="results/top_common_snp_per_gene_by_composite_score.tsv"

echo "== Inputs =="
echo "RNAfold TSV: $RNAFOLD_TSV"
echo "Score TSV:   $SCORE_TSV"
echo

for f in "$RNAFOLD_TSV" "$SCORE_TSV"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing input file: $f" >&2
    exit 1
  fi
done

###############################################################################
# Step 1: parse local accessibility from RNAfold dot-bracket structures
###############################################################################

echo "== Step 1: parse central 17 nt accessibility =="

python3 - <<'PY'
from pathlib import Path
import csv

infile = Path("results/rnafold/rnafold_plusminus100_mfe.tsv")
outfile = Path("results/rnafold/rnafold_plusminus100_local_accessibility.tsv")

# 201 nt window centred on variant:
# central base is index 100 in 0-based coordinates.
# 17-mer target is 8 bases either side: 92:109.
TARGET_START_0 = 92
TARGET_END_0 = 109
CENTER_0 = 100

def clean_structure(s):
    return s.strip().split()[0]

with open(infile) as f, open(outfile, "w", newline="") as out:
    reader = csv.DictReader(f, delimiter="\t")

    fieldnames = reader.fieldnames + [
        "target17_structure",
        "target17_unpaired_bases",
        "target17_paired_bases",
        "target17_unpaired_fraction",
        "central_base_structure",
        "central_base_unpaired",
        "target17_accessibility_score"
    ]

    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    n = 0
    bad = 0

    for row in reader:
        n += 1

        seq = row["sequence"]
        struct = clean_structure(row["structure"])

        if len(seq) < TARGET_END_0 or len(struct) < TARGET_END_0:
            bad += 1
            row["target17_structure"] = "NA"
            row["target17_unpaired_bases"] = "NA"
            row["target17_paired_bases"] = "NA"
            row["target17_unpaired_fraction"] = "NA"
            row["central_base_structure"] = "NA"
            row["central_base_unpaired"] = "NA"
            row["target17_accessibility_score"] = "NA"
            writer.writerow(row)
            continue

        target_struct = struct[TARGET_START_0:TARGET_END_0]
        unpaired = target_struct.count(".")
        paired = len(target_struct) - unpaired
        unpaired_fraction = unpaired / len(target_struct)

        central_base_structure = struct[CENTER_0]
        central_base_unpaired = 1 if central_base_structure == "." else 0
        accessibility_score = 100 * unpaired_fraction

        row["target17_structure"] = target_struct
        row["target17_unpaired_bases"] = unpaired
        row["target17_paired_bases"] = paired
        row["target17_unpaired_fraction"] = unpaired_fraction
        row["central_base_structure"] = central_base_structure
        row["central_base_unpaired"] = central_base_unpaired
        row["target17_accessibility_score"] = accessibility_score

        writer.writerow(row)

print(f"Wrote {outfile}")
print(f"Rows parsed: {n}")
print(f"Rows with short sequence/structure: {bad}")
PY

echo "Accessibility table check:"
wc -l "$ACCESS_TSV"
awk -F'\t' '{print NF}' "$ACCESS_TSV" | sort | uniq -c
echo

###############################################################################
# Step 2: join local accessibility to 17-mer score table
###############################################################################

echo "== Step 2: join local accessibility to 17-mer score table =="

python3 - <<'PY'
from pathlib import Path
import csv

score_file = Path("results/aso_17mer_scores_with_rnafold.tsv")
access_file = Path("results/rnafold/rnafold_plusminus100_local_accessibility.tsv")
out_file = Path("results/aso_17mer_scores_with_rnafold_accessibility.tsv")

def key(row):
    return row["full_header"].split("::")[0]

access = {}

with open(access_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    duplicate = 0

    for row in reader:
        k = key(row)
        if k in access:
            duplicate += 1
        access[k] = row

extra_fields = [
    "target17_structure",
    "target17_unpaired_bases",
    "target17_paired_bases",
    "target17_unpaired_fraction",
    "central_base_structure",
    "central_base_unpaired",
    "target17_accessibility_score"
]

with open(score_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    fields = reader.fieldnames

missing = 0

with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields + extra_fields, delimiter="\t")
    writer.writeheader()

    for row in rows:
        k = key(row)
        arow = access.get(k)

        if arow is None:
            missing += 1
            for x in extra_fields:
                row[x] = "NA"
        else:
            for x in extra_fields:
                row[x] = arow[x]

        writer.writerow(row)

print(f"Wrote {out_file}")
print(f"Score rows: {len(rows)}")
print(f"Accessibility rows: {len(access)}")
print(f"Duplicate accessibility keys: {duplicate}")
print(f"Missing accessibility matches: {missing}")
PY

echo "Joined table check:"
wc -l "$JOINED_TSV"
awk -F'\t' '{print NF}' "$JOINED_TSV" | sort | uniq -c
echo

###############################################################################
# Step 3: calculate continuous sequence score and composite score
###############################################################################

echo "== Step 3: calculate sequence and composite scores =="

python3 - <<'PY'
from pathlib import Path
import csv
import math

infile = Path("results/aso_17mer_scores_with_rnafold_accessibility.tsv")
outfile = Path("results/aso_17mer_plotting_table.tsv")
topfile = Path("results/top_common_snp_per_gene_by_composite_score.tsv")

def to_float(x):
    try:
        if x == "NA" or x == "":
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")

rows = []

with open(infile) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        gc = to_float(row["gc_frac"])
        central9_gc = to_float(row["central9_gc_frac"])
        max_homopolymer = to_float(row["max_homopolymer"])
        distinct_dinucs = to_float(row["distinct_dinucs"])
        cpg_count = to_float(row["cpg_count"])
        target_unpaired = to_float(row["target17_unpaired_fraction"])

        if any(math.isnan(x) for x in [
            gc, central9_gc, max_homopolymer,
            distinct_dinucs, cpg_count, target_unpaired
        ]):
            sequence_score = "NA"
            accessibility_score = "NA"
            composite_score = "NA"
        else:
            gc_penalty = abs(gc - 0.50) * 120
            central9_penalty = abs(central9_gc - 0.50) * 80
            homopolymer_penalty = max(0, max_homopolymer - 3) ** 2 * 8
            dinuc_penalty = max(0, 8 - distinct_dinucs) * 3
            cpg_penalty = max(0, cpg_count - 1) * 4
            extreme_gc_penalty = 15 if (gc < 0.25 or gc > 0.75) else 0

            total_penalty = (
                gc_penalty +
                central9_penalty +
                homopolymer_penalty +
                dinuc_penalty +
                cpg_penalty +
                extreme_gc_penalty
            )

            sequence_score = max(0, 100 - total_penalty)
            accessibility_score = 100 * target_unpaired
            composite_score = 0.70 * sequence_score + 0.30 * accessibility_score

        row["sequence_score"] = sequence_score
        row["accessibility_score"] = accessibility_score
        row["composite_score"] = composite_score

        rows.append(row)

keep_fields = [
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
    "score_total",
    "sequence_score",
    "target17_unpaired_fraction",
    "target17_accessibility_score",
    "accessibility_score",
    "composite_score",
    "central_base_unpaired",
    "gc_frac",
    "central9_gc_frac",
    "max_homopolymer",
    "distinct_dinucs",
    "cpg_count",
    "mfe",
    "mfe_per_nt",
    "target17_structure",
    "full_header"
]

with open(outfile, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=keep_fields, delimiter="\t")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: row.get(k, "NA") for k in keep_fields})

# Top common SNP per gene by composite score
common = [
    r for r in rows
    if r["variant_class"] == "common_snp" and r["composite_score"] != "NA"
]

best = {}
for r in common:
    g = r["gene"]
    score = float(r["composite_score"])
    het = to_float(r["het_rate"])
    maf = to_float(r["maf"])
    key = (score, het, maf)

    if g not in best or key > best[g][0]:
        best[g] = (key, r)

top_fields = [
    "gene",
    "variant_id",
    "context",
    "maf",
    "het_rate",
    "seq_17mer",
    "sequence_score",
    "accessibility_score",
    "composite_score",
    "target17_unpaired_fraction",
    "central_base_unpaired",
    "target17_structure",
    "full_header"
]

with open(topfile, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=top_fields, delimiter="\t")
    writer.writeheader()

    for g in sorted(best):
        r = best[g][1]
        writer.writerow({k: r.get(k, "NA") for k in top_fields})

print(f"Wrote {outfile}")
print(f"Wrote {topfile}")
print(f"Rows written: {len(rows)}")
print(f"Top common SNP genes: {len(best)}")
PY

echo "Plotting table check:"
wc -l "$PLOT_TSV"
awk -F'\t' '{print NF}' "$PLOT_TSV" | sort | uniq -c
head "$PLOT_TSV" | column -t
echo

echo "Top common SNP table check:"
wc -l "$TOP_COMMON_TSV"
head "$TOP_COMMON_TSV" | column -t
echo

###############################################################################
# Step 4: summary statistics
###############################################################################

echo "== Step 4: quick summaries =="

echo "Counts by variant class:"
awk -F'\t' 'NR>1 {count[$1]++} END{for (k in count) print k, count[k]}' "$PLOT_TSV" | sort
echo

echo "Mean target17 unpaired fraction by variant class:"
awk -F'\t' 'NR>1 {
  class=$1
  sum[class]+=$14
  n[class]++
}
END {
  for (c in n) print c, n[c], sum[c]/n[c]
}' "$PLOT_TSV" | sort
echo

echo "High-accessibility target sites by variant class:"
awk -F'\t' 'NR>1 {
  class=$1
  n[class]++
  if ($14 >= 0.75) high[class]++
  if ($18 == 1) central[class]++
}
END {
  for (c in n) {
    print c, "N="n[c], "target17_unpaired>=0.75="high[c], "central_base_unpaired="central[c]
  }
}' "$PLOT_TSV" | sort
echo

echo "== Done =="
echo "Download for R:"
echo "  $PLOT_TSV"
echo
echo "Optional top-site table:"
echo "  $TOP_COMMON_TSV"
