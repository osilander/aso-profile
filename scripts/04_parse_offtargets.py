#!/usr/bin/env python3
"""
Parse bowtie SAM output → per-ASO off-target summary.

Metrics per query × index:
  n_hits_0mm  : exact off-target hits (outside intended gene)
  n_hits_1mm  : 1-mismatch off-target hits
  n_hits_2mm  : 2-mismatch off-target hits
  max_perfect_stretch : longest contiguous perfect-match run in any hit
  n_hits_stretch_12+  : hits with longest perfect stretch ≥ 12 nt
  n_hits_stretch_13+  : ≥ 13 nt
  n_hits_stretch_14+  : ≥ 14 nt
  paralog_hits        : hits to known paralogous genes (if paralog file provided)

  If --repeats is provided (from 05_annotate_repeats.sh), adds:
  repeat_overlap_any  : whether 17-mer overlaps a RepeatMasker element
  repeat_class        : repeat class (LINE, SINE, etc.)
  repeat_family       : repeat family
  repeat_fraction     : fraction of 17-mer overlapping a repeat

  And three pass_specificity tiers:
  pass_specificity_strict     : no repeat + no 0mm/1mm off-target + no paralog
  pass_specificity_moderate   : no repeat + no 0mm off-target + no paralog exact
  pass_specificity_permissive : no 0mm off-target (allows repeat and 1mm)

Usage:
  python3 scripts/04_parse_offtargets.py \
      --meta     results/offtarget/query_meta.tsv \
      --samdir   results/offtarget \
      --out      results/offtarget/offtarget_summary.tsv \
      [--paralogs data/paralog_pairs.tsv] \
      [--repeats  results/offtarget/repeat_overlap.tsv]

Paralog file format (tab-separated, no header):
  GENE_A  GENE_B
"""

import argparse, os, re, sys
from collections import defaultdict
import pandas as pd

def longest_perfect_stretch(cigar_str, n_mm, seq_len=17):
    """
    Given the MD tag (mismatch descriptor) from a SAM record and
    number of mismatches, return the longest run of consecutive
    matched bases.

    MD tag examples:
      "17"     → perfect match
      "8A8"    → 1mm at pos 8
      "5T2A8"  → 2mm at pos 5 and 8
    """
    if n_mm == 0:
        return seq_len

    # Parse MD: alternating numbers (matches) and letters (mismatches/deletions)
    tokens = re.findall(r'(\d+|[A-Z^]+)', cigar_str)
    runs = []
    pos  = 0
    for tok in tokens:
        if tok.isdigit():
            runs.append(int(tok))
        else:
            # mismatch base(s) — count as 0 matches
            runs.append(0)
    return max(runs) if runs else seq_len

def parse_sam_hits(sam_file, intended_gene, paralogs_of_gene):
    """
    Read a SAM file and return per-hit records excluding the intended gene.
    Returns list of dicts with: ref_name, n_mm, longest_stretch, is_paralog
    """
    hits = []
    with open(sam_file) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag    = int(fields[1])
            if flag & 4:          # unmapped
                continue
            ref_name = fields[2]

            # Extract gene name from reference name
            # GENCODE transcriptome: ENST... (transcript ID) — need GTF lookup
            # gene_body index: gene_name from BED name field
            # Simple heuristic: split on '|' or '::' and take gene field
            ref_gene = extract_gene(ref_name)

            if ref_gene == intended_gene:
                continue          # skip on-target hits

            # Count mismatches from NM tag
            nm = 0
            md = ""
            for opt in fields[11:]:
                if opt.startswith("NM:i:"):
                    nm = int(opt[5:])
                if opt.startswith("MD:Z:"):
                    md = opt[5:]

            stretch = longest_perfect_stretch(md, nm)
            is_para = ref_gene in paralogs_of_gene

            hits.append({
                "ref_name":       ref_name,
                "ref_gene":       ref_gene,
                "n_mm":           nm,
                "longest_stretch": stretch,
                "is_paralog":     is_para,
            })
    return hits

def extract_gene(ref_name):
    """
    Extract gene name from bowtie reference sequence name.
    Handles:
      GENCODE cDNA:   ENST00000... (need gene_name from attributes — use pipe-separated FASTA header)
                      e.g. ENST000|ENSG000|OTTHUMG|OTTHUMT|gene_name|transcript_name|length|...
      gene_body BED:  GENENAME::chr1:1000-2000(+)
    """
    # gene_body: name::coords
    if "::" in ref_name:
        return ref_name.split("::")[0]
    # GENCODE pipe-separated (pipe fields: transcript|gene|HAVANA_t|HAVANA_g|gene_name|...)
    parts = ref_name.split("|")
    if len(parts) >= 6:
        return parts[5]   # gene_name field
    return ref_name

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta",     required=True)
    ap.add_argument("--samdir",   required=True)
    ap.add_argument("--out",      required=True)
    ap.add_argument("--paralogs", default=None)
    ap.add_argument("--repeats",  default=None,
                    help="repeat_overlap.tsv from 05_annotate_repeats.sh")
    args = ap.parse_args()

    meta = pd.read_csv(args.meta, sep="\t")

    # Load paralogs if provided
    paralog_map = defaultdict(set)   # gene → set of paralog gene names
    if args.paralogs and os.path.exists(args.paralogs):
        for line in open(args.paralogs):
            a, b = line.strip().split("\t")[:2]
            paralog_map[a].add(b)
            paralog_map[b].add(a)
        print(f"Loaded paralogs for {len(paralog_map)} genes")

    # Find all SAM files
    sam_files = [f for f in os.listdir(args.samdir) if f.endswith(".bwt")]

    # Build meta lookup: query_id → (gene, vc, allele, ctx, variant_id, paralogs)
    meta_by_qid = {}
    for _, row in meta.iterrows():
        qid = row["query_id"]
        vid_parts  = qid.split("|")
        variant_id = vid_parts[2] if len(vid_parts) >= 3 else ""
        meta_by_qid[qid] = {
            "gene":       row["gene"],
            "vc":         row["variant_class"],
            "allele":     row["allele"],
            "ctx":        row["context"],
            "variant_id": variant_id,
            "paras":      paralog_map.get(row["gene"], set()),
        }

    records = []
    for sam_name in sam_files:
        index_type = sam_name.replace(".bwt", "").split("_vs_")[-1]
        sam_path   = os.path.join(args.samdir, sam_name)
        print(f"Reading {sam_name} ...", flush=True)

        # Read entire SAM once → dict: query_id → list of hit dicts
        hits_by_qid = defaultdict(list)
        seen_qids   = set()
        with open(sam_path) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue
                qid  = fields[0]
                seen_qids.add(qid)
                if int(fields[1]) & 4:   # unmapped
                    continue
                if qid not in meta_by_qid:
                    continue
                ref_gene = extract_gene(fields[2])
                intended = meta_by_qid[qid]["gene"]
                if ref_gene == intended:
                    continue
                nm, md = 0, ""
                for opt in fields[11:]:
                    if opt.startswith("NM:i:"): nm = int(opt[5:])
                    if opt.startswith("MD:Z:"): md = opt[5:]
                stretch  = longest_perfect_stretch(md, nm)
                is_para  = ref_gene in meta_by_qid[qid]["paras"]
                hits_by_qid[qid].append({
                    "n_mm":            nm,
                    "longest_stretch": stretch,
                    "is_paralog":      is_para,
                })

        # Summarise per query
        for qid, m in meta_by_qid.items():
            if qid not in seen_qids:
                continue
            hits = hits_by_qid.get(qid, [])
            records.append({
                "query_id":            qid,
                "variant_id":          m["variant_id"],
                "gene":                m["gene"],
                "variant_class":       m["vc"],
                "allele":              m["allele"],
                "context":             m["ctx"],
                "index":               index_type,
                "n_offtarget_0mm":     sum(1 for h in hits if h["n_mm"] == 0),
                "n_offtarget_1mm":     sum(1 for h in hits if h["n_mm"] == 1),
                "n_offtarget_2mm":     sum(1 for h in hits if h["n_mm"] == 2),
                "max_perfect_stretch": max((h["longest_stretch"] for h in hits), default=0),
                "n_stretch_ge12":      sum(1 for h in hits if h["longest_stretch"] >= 12),
                "n_stretch_ge13":      sum(1 for h in hits if h["longest_stretch"] >= 13),
                "n_stretch_ge14":      sum(1 for h in hits if h["longest_stretch"] >= 14),
                "n_paralog_hits":      sum(1 for h in hits if h["is_paralog"]),
            })
        print(f"  → {len(seen_qids)} queries seen, {len(hits_by_qid)} with off-target hits",
              flush=True)

    out_df = pd.DataFrame(records)

    # Optionally join repeat annotation and compute pass_specificity tiers
    if args.repeats and os.path.exists(args.repeats):
        rep = pd.read_csv(args.repeats, sep="\t")
        rep_cols = ["variant_id", "repeat_overlap_any",
                    "repeat_class", "repeat_family", "repeat_fraction"]
        out_df = out_df.merge(rep[rep_cols], on="variant_id", how="left")
        out_df["repeat_overlap_any"] = out_df["repeat_overlap_any"].fillna(False)

        # strict: no repeat + no 0/1mm off-target + no paralog 2mm
        out_df["pass_specificity_strict"] = (
            ~out_df["repeat_overlap_any"] &
            (out_df["n_offtarget_0mm"] == 0) &
            (out_df["n_offtarget_1mm"] == 0) &
            (out_df["n_paralog_hits"]  == 0)
        )
        # moderate: no repeat + no exact off-target + no exact paralog
        out_df["pass_specificity_moderate"] = (
            ~out_df["repeat_overlap_any"] &
            (out_df["n_offtarget_0mm"] == 0) &
            (out_df["n_paralog_hits"]  == 0)
        )
        # permissive: just no exact off-target (allow repeat + 1mm)
        out_df["pass_specificity_permissive"] = (
            out_df["n_offtarget_0mm"] == 0
        )
        n_strict = out_df.groupby("query_id")["pass_specificity_strict"].all().sum()
        print(f"Loaded repeat annotation; "
              f"{n_strict}/{out_df['query_id'].nunique()} candidates pass strict filter")
    else:
        print("No repeat annotation provided; skipping pass_specificity columns.")

    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"Written: {args.out}  ({len(out_df)} rows)")

def parse_sam_hits_for_query(sam_path, query_id, intended_gene, paralogs):
    """Return list of hit dicts for a specific query_id in a SAM file."""
    hits = []
    found = False
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] != query_id:
                continue
            found = True
            flag = int(fields[1])
            if flag & 4:
                continue
            ref_name = fields[2]
            ref_gene = extract_gene(ref_name)
            if ref_gene == intended_gene:
                continue
            nm, md = 0, ""
            for opt in fields[11:]:
                if opt.startswith("NM:i:"): nm = int(opt[5:])
                if opt.startswith("MD:Z:"): md = opt[5:]
            stretch = longest_perfect_stretch(md, nm)
            hits.append({
                "ref_gene":        ref_gene,
                "n_mm":            nm,
                "longest_stretch": stretch,
                "is_paralog":      ref_gene in paralogs,
            })
    return hits if found else None

if __name__ == "__main__":
    main()
