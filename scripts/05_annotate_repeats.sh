#!/usr/bin/env bash
# Annotate each 17-mer candidate window with RepeatMasker overlap.
#
# Input:  results/aso_17mer_scores_extended.tsv  (full_header col has chr:start-end)
# Output: results/offtarget/repeat_overlap.tsv
#
# Requires: bedtools, python3 (pandas), wget
# module load BEDTools Python  # adjust module names for your HPC

set -euo pipefail

module load BEDTools 2>/dev/null || true
module load Python 2>/dev/null || true

RMSK_BED="data/rmsk.hg38.bed.gz"
OUTD="results/offtarget"
TSV="results/aso_17mer_scores_extended.tsv"
mkdir -p "${OUTD}" data

# ------------------------------------------------------------------ #
# 1. Download RepeatMasker table from UCSC if not present
# ------------------------------------------------------------------ #
# rmsk.txt.gz columns (1-indexed):
#   6=genoName(chrom)  7=genoStart  8=genoEnd  11=repName
#   10=strand  12=repClass  13=repFamily
if [ ! -f "${RMSK_BED}" ]; then
    echo "Downloading UCSC hg38 RepeatMasker table (~200 MB) ..."
    RMSK_TMP="data/rmsk.txt.gz"
    wget -q --show-progress -O "${RMSK_TMP}" \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
    echo "Converting to BED ..."
    zcat "${RMSK_TMP}" \
        | awk 'BEGIN{OFS="\t"} {print $6, $7, $8, $11, int($2), $10, $12, $13}' \
        | gzip -c > "${RMSK_BED}"
    rm "${RMSK_TMP}"
    echo "Written: ${RMSK_BED}"
fi

# ------------------------------------------------------------------ #
# 2. Build BED file of 17-mer candidate windows from full_header col
# ------------------------------------------------------------------ #
echo "Extracting 17-mer intervals from ${TSV} ..."
python3 - <<'PYEOF'
import pandas as pd, re, sys

df = pd.read_csv("results/aso_17mer_scores_extended.tsv", sep="\t")

records = []
for _, row in df.iterrows():
    fh = str(row.get("full_header", ""))
    # e.g. "intron::chr1:11108463-11108480(-)"
    m = re.search(r'::(chr[^:]+):(\d+)-(\d+)\([+-]\)', fh)
    if not m:
        sys.stderr.write(f"WARNING: cannot parse coords from: {fh}\n")
        continue
    chrom = m.group(1)
    start = int(m.group(2))
    end   = int(m.group(3))
    records.append({
        "chrom":         chrom,
        "start":         start,
        "end":           end,
        "variant_id":    str(row["variant_id"]),
        "gene":          row["gene"],
        "variant_class": row["variant_class"],
    })

bed = pd.DataFrame(records)
bed[["chrom","start","end","variant_id","gene","variant_class"]].to_csv(
    "results/offtarget/candidates_17mer.bed",
    sep="\t", header=False, index=False)
print(f"Written {len(bed)} intervals → results/offtarget/candidates_17mer.bed",
      flush=True)
PYEOF

# ------------------------------------------------------------------ #
# 3. Bedtools intersect: each 17-mer window vs RepeatMasker
#    -wao: report each A row even with no B overlap (overlap=0 col added)
# ------------------------------------------------------------------ #
echo "Running bedtools intersect ..."

TMP_PY=$(mktemp /tmp/parse_rmsk_XXXXXX.py)
cat > "${TMP_PY}" << 'PYEOF'
import sys, pandas as pd

# Candidate cols:  0=chrom 1=start 2=end 3=variant_id 4=gene 5=variant_class
# RepeatMasker:    6=chrom 7=start 8=end 9=repName 10=score 11=strand 12=repClass 13=repFamily
# Overlap col:    14=overlap_bp

rows = []
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    win_len    = int(f[2]) - int(f[1])
    variant_id = f[3]
    gene       = f[4]
    var_class  = f[5]
    no_hit     = (f[6] == ".")

    if no_hit:
        rows.append({
            "variant_id":         variant_id,
            "gene":               gene,
            "variant_class":      var_class,
            "repeat_overlap_any": False,
            "repeat_class":       None,
            "repeat_family":      None,
            "repeat_fraction":    0.0,
        })
    else:
        rep_class  = f[12]
        rep_family = f[13]
        overlap_bp = int(f[14]) if f[14].lstrip("-").isdigit() else 0
        rows.append({
            "variant_id":         variant_id,
            "gene":               gene,
            "variant_class":      var_class,
            "repeat_overlap_any": True,
            "repeat_class":       rep_class,
            "repeat_family":      rep_family,
            "repeat_fraction":    overlap_bp / win_len if win_len > 0 else 0.0,
        })

df = pd.DataFrame(rows)

def agg_fn(grp):
    if grp["repeat_overlap_any"].any():
        hit = grp[grp["repeat_overlap_any"]].sort_values(
            "repeat_fraction", ascending=False).iloc[0]
        return pd.Series({
            "gene":               hit["gene"],
            "variant_class":      hit["variant_class"],
            "repeat_overlap_any": True,
            "repeat_class":       hit["repeat_class"],
            "repeat_family":      hit["repeat_family"],
            "repeat_fraction":    hit["repeat_fraction"],
        })
    else:
        row = grp.iloc[0]
        return pd.Series({
            "gene":               row["gene"],
            "variant_class":      row["variant_class"],
            "repeat_overlap_any": False,
            "repeat_class":       None,
            "repeat_family":      None,
            "repeat_fraction":    0.0,
        })

out = df.groupby("variant_id", sort=False).apply(agg_fn).reset_index()
out.to_csv("results/offtarget/repeat_overlap.tsv", sep="\t", index=False)
print(f"Written {len(out)} rows → results/offtarget/repeat_overlap.tsv", flush=True)
n_rep = out["repeat_overlap_any"].sum()
print(f"  {n_rep}/{len(out)} candidates overlap a repeat element "
      f"({100*n_rep/max(len(out),1):.1f}%)")
PYEOF

bedtools intersect \
    -a "${OUTD}/candidates_17mer.bed" \
    -b "${RMSK_BED}" \
    -wao \
    | python3 "${TMP_PY}"

rm -f "${TMP_PY}"

echo "Done."
