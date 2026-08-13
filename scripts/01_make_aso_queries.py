#!/usr/bin/env python3
"""
Generate query FASTA files for ASO off-target search.

seq_17mer = mRNA target sense sequence (17 nt, variant at centre index 8).
ASO sequence = revcomp(seq_17mer).

For common SNPs  → two queries: ASO-REF (targets REF mRNA) and ASO-ALT (targets ALT mRNA)
For pathogenic   → one query:  ASO-ALT only (must target the pathogenic allele)

Output: results/offtarget/queries_common_ref.fa
        results/offtarget/queries_common_alt.fa
        results/offtarget/queries_pathogenic_alt.fa
        results/offtarget/query_meta.tsv   (ID → gene, class, allele, context)
"""

import os, sys
import pandas as pd

TSV   = "aso_17mer_scores_extended.tsv"
OUTD  = "results/offtarget"
os.makedirs(OUTD, exist_ok=True)

COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq):
    return seq.translate(COMP)[::-1]

def mRNA_alt_base(ref_genomic, alt_genomic, strand):
    """Return the ALT mRNA base at the variant position."""
    if strand == "+":
        return alt_genomic.upper()
    else:
        return alt_genomic.upper().translate(COMP)

df = pd.read_csv(TSV, sep="\t")
df = df[df["pass_basic"] == 1].copy()

CENTRE = 8  # 0-indexed centre of 17-mer

# Sanity check: centre base of seq_17mer should match the mRNA REF allele
def expected_mrna_ref(row):
    if row["strand"] == "+":
        return row["ref"].upper()
    else:
        return row["ref"].upper().translate(COMP)

df["expected_centre"] = df.apply(expected_mrna_ref, axis=1)
df["actual_centre"]   = df["seq_17mer"].str[CENTRE].str.upper()
mismatch = df[df["expected_centre"] != df["actual_centre"]]
if len(mismatch):
    print(f"WARNING: {len(mismatch)} rows where seq_17mer[8] ≠ expected mRNA REF base")
    print(mismatch[["gene","variant_id","strand","ref","expected_centre","actual_centre"]].head(10).to_string())
else:
    print(f"Sanity check passed: seq_17mer[8] matches mRNA REF base in all {len(df)} rows")

records = {"common_ref": [], "common_alt": [],
           "pathogenic_alt": [], "pathogenic_ref": []}
meta    = []

for _, row in df.iterrows():
    vc      = row["variant_class"]
    gene    = row["gene"]
    vid     = row["variant_id"]
    strand  = row["strand"]
    ref_g   = row["ref"].upper()
    alt_g   = row["alt"].upper()
    seq     = row["seq_17mer"].upper()
    context = row.get("context", "")

    uid = f"{vc}|{gene}|{vid}"

    # REF mRNA target sequence = seq_17mer as-is
    # ALT mRNA target: swap centre base
    alt_centre = mRNA_alt_base(ref_g, alt_g, strand)
    seq_alt    = seq[:CENTRE] + alt_centre + seq[CENTRE+1:]

    if vc == "common_snp":
        records["common_ref"].append((f"{uid}|REF", seq))
        records["common_alt"].append((f"{uid}|ALT", seq_alt))
        meta.append((f"{uid}|REF", gene, "common_snp", "REF", context))
        meta.append((f"{uid}|ALT", gene, "common_snp", "ALT", context))
    else:  # pathogenic_missense
        records["pathogenic_alt"].append((f"{uid}|ALT", seq_alt))
        meta.append((f"{uid}|ALT", gene, "pathogenic_missense", "ALT", context))
        # REF/WT allele — comparison only, not a therapeutic candidate
        records["pathogenic_ref"].append((f"{uid}|REF", seq))
        meta.append((f"{uid}|REF", gene, "pathogenic_missense", "REF", context))

for key, recs in records.items():
    fa = os.path.join(OUTD, f"queries_{key}.fa")
    with open(fa, "w") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n{seq}\n")
    print(f"Written {len(recs):>5} sequences → {fa}")

meta_df = pd.DataFrame(meta, columns=["query_id","gene","variant_class","allele","context"])
meta_tsv = os.path.join(OUTD, "query_meta.tsv")
meta_df.to_csv(meta_tsv, sep="\t", index=False)
print(f"Written metadata → {meta_tsv}")
