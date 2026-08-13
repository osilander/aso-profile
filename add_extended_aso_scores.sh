#!/usr/bin/env bash
set -euo pipefail

# Add extended ASO scores for both common SNP handles AND pathogenic SNVs.
#
# For every 17-mer, two ASOs are scored:
#   ASO-REF: designed on REF allele; match=REF, mismatch=ALT → duplex_ddg_aso_ref
#   ASO-ALT: designed on ALT allele; match=ALT, mismatch=REF → duplex_ddg_aso_alt
# ΔΔG = dg_mismatch - dg_match; larger positive = more allele-selective.
# Both alleles are scored for all variants (common SNPs and pathogenic SNVs)
# because phasing is not known at this stage.
#
# Also adds 9 binary 4-mer motif columns (Hill et al. 2026):
#   negative: GGGG AAAA TAAA CTAA CCTA
#   positive: GCGT TTGT GTAT CGTA GTCG
#
# Required inputs:
#   results/aso_17mer.common_snps.fa           from score_17mer_aso_windows.sh
#   results/aso_17mer.pathogenic_missense.fa   from score_17mer_aso_windows.sh
#   results/aso_17mer_scores_with_rnafold.tsv  from run_rnafold_and_join_aso_scores.sh
#
# RNAduplex inputs (4 per variant class, generated in step 1):
#   results/rduplex/<class>.aso_ref.match.in
#   results/rduplex/<class>.aso_ref.mismatch.in
#   results/rduplex/<class>.aso_alt.match.in
#   results/rduplex/<class>.aso_alt.mismatch.in
#
# Output:
#   results/aso_17mer_scores_extended.tsv
#   New columns: duplex_dg_aso_ref_match  duplex_dg_aso_ref_mismatch  duplex_ddg_aso_ref
#                duplex_dg_aso_alt_match  duplex_dg_aso_alt_mismatch  duplex_ddg_aso_alt
#                variant_pos_in_17mer
#                motif_GGGG  motif_AAAA  motif_TAAA  motif_CTAA  motif_CCTA
#                motif_GCGT  motif_TTGT  motif_GTAT  motif_CGTA  motif_GTCG
#
# Usage:
#   ./add_extended_aso_scores.sh                   full run (prep + RNAduplex + join)
#   SKIP_RDUPLEX=1 ./add_extended_aso_scores.sh    join only (RNAduplex already run on HPC)

COMMON_FA="results/aso_17mer.common_snps.fa"
PATHOGENIC_FA="results/aso_17mer.pathogenic_missense.fa"
SCORE_TSV="results/aso_17mer_scores_with_rnafold.tsv"
RDUPLEX_DIR="results/rduplex"
OUT_TSV="results/aso_17mer_scores_extended.tsv"

SKIP_RDUPLEX="${SKIP_RDUPLEX:-0}"

mkdir -p "$RDUPLEX_DIR" logs

for f in "$COMMON_FA" "$PATHOGENIC_FA" "$SCORE_TSV"; do
  [[ -s "$f" ]] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done

echo "###############################################################################"
echo "Extended ASO scoring pipeline (common SNPs + pathogenic SNVs, both alleles)"
echo "###############################################################################"
echo

###############################################################################
# Step 1: generate RNAduplex inputs and motif score tables
###############################################################################

echo "###############################################################################"
echo "Step 1: generate duplex inputs and motif scores"
echo "###############################################################################"

python3 - "$COMMON_FA" "$PATHOGENIC_FA" "$RDUPLEX_DIR" <<'PY'
import sys
import csv
from pathlib import Path

common_fa     = sys.argv[1]
pathogenic_fa = sys.argv[2]
outdir        = Path(sys.argv[3])

RNA_COMP = str.maketrans('ACGUacgu', 'UGCAugca')
DNA_COMP = str.maketrans('ACGTacgt', 'TGCAtgca')

def revcomp_rna(seq):
    return seq.translate(RNA_COMP)[::-1]

def comp_dna(base):
    return base.upper().translate(DNA_COMP)

NEGATIVE_MOTIFS = ['GGGG', 'AAAA', 'TAAA', 'CTAA', 'CCTA']
POSITIVE_MOTIFS = ['GCGT', 'TTGT', 'GTAT', 'CGTA', 'GTCG']
ALL_MOTIFS  = NEGATIVE_MOTIFS + POSITIVE_MOTIFS
MOTIF_COLS  = [f'motif_{m}' for m in ALL_MOTIFS]

def read_fasta(path):
    header = None
    parts  = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts).upper()
                header = line[1:]
                parts  = []
            else:
                parts.append(line)
    if header is not None:
        yield header, ''.join(parts).upper()

def allele_in_transcript(allele_genomic, strand):
    """Genomic allele → transcript-strand base."""
    a = allele_genomic.upper()
    return a if strand == '+' else comp_dna(a)

def process_fasta(fasta_path, label):
    """
    For every valid 17-mer return four RNAduplex pairs:
      ASO-REF vs REF  (perfect match)
      ASO-REF vs ALT  (one mismatch)
      ASO-ALT vs ALT  (perfect match)
      ASO-ALT vs REF  (one mismatch)
    The variant/SNP is always at position 8 (0-indexed) = position 9 (1-indexed).
    """
    keys              = []
    aso_ref_match     = []   # ASO-REF vs REF
    aso_ref_mismatch  = []   # ASO-REF vs ALT
    aso_alt_match     = []   # ASO-ALT vs ALT
    aso_alt_mismatch  = []   # ASO-ALT vs REF
    motif_rows        = []
    skipped           = 0

    for header, seq in read_fasta(fasta_path):
        if len(seq) != 17 or any(b not in 'ACGT' for b in seq):
            skipped += 1
            continue

        meta  = header.split('::')[0]
        parts = meta.split('|')
        while len(parts) < 10:
            parts.append('NA')

        strand  = parts[8]
        ref_gen = parts[3]
        alt_gen = parts[4]

        ref_trans = allele_in_transcript(ref_gen, strand)
        alt_trans = allele_in_transcript(alt_gen, strand)

        # REF and ALT 17-mer sequences in transcript orientation (RNA)
        ref_rna = seq.replace('T', 'U')
        alt_seq = seq[:8] + alt_trans.replace('U', 'T') + seq[9:]
        alt_rna = alt_seq.replace('T', 'U')

        aso_ref = revcomp_rna(ref_rna)   # ASO complementary to REF allele
        aso_alt = revcomp_rna(alt_rna)   # ASO complementary to ALT allele

        keys.append(meta)
        aso_ref_match.append((aso_ref, ref_rna))
        aso_ref_mismatch.append((aso_ref, alt_rna))
        aso_alt_match.append((aso_alt, alt_rna))
        aso_alt_mismatch.append((aso_alt, ref_rna))

        flags = {f'motif_{m}': int(m in seq) for m in ALL_MOTIFS}
        motif_rows.append({'key': meta, 'variant_pos_in_17mer': 9, **flags})

    print(f"{label}: {len(keys)} valid, {skipped} skipped")
    return keys, aso_ref_match, aso_ref_mismatch, aso_alt_match, aso_alt_mismatch, motif_rows

def write_pairs(path, pairs):
    with open(path, 'w') as f:
        for aso, target in pairs:
            f.write(aso + '\n' + target + '\n')

def write_outputs(label, keys, ref_match, ref_mis, alt_match, alt_mis, motif_rows, outdir):
    write_pairs(outdir / f'{label}.aso_ref.match.in',    ref_match)
    write_pairs(outdir / f'{label}.aso_ref.mismatch.in', ref_mis)
    write_pairs(outdir / f'{label}.aso_alt.match.in',    alt_match)
    write_pairs(outdir / f'{label}.aso_alt.mismatch.in', alt_mis)

    with open(outdir / f'{label}.duplex_keys.txt', 'w') as f:
        for k in keys:
            f.write(k + '\n')

    with open(outdir / f'{label}.motif_scores.tsv', 'w', newline='') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=['key', 'variant_pos_in_17mer'] + MOTIF_COLS,
            delimiter='\t'
        )
        writer.writeheader()
        writer.writerows(motif_rows)

    print(f"  Wrote 4 duplex input files + keys + motif scores for {label}")

for fasta, label in [(common_fa, 'common_snps'), (pathogenic_fa, 'pathogenic_missense')]:
    data = process_fasta(fasta, label)
    write_outputs(label, *data, outdir)
PY

echo

###############################################################################
# Step 2: run RNAduplex (4 runs per variant class = 8 total)
###############################################################################

echo "###############################################################################"
echo "Step 2: run RNAduplex"
echo "###############################################################################"

ARMS="aso_ref.match aso_ref.mismatch aso_alt.match aso_alt.mismatch"

if [[ "$SKIP_RDUPLEX" == "1" ]]; then
  echo "SKIP_RDUPLEX=1 — skipping (expecting outputs already in $RDUPLEX_DIR)"
  for CLASS in common_snps pathogenic_missense; do
    for ARM in $ARMS; do
      f="$RDUPLEX_DIR/${CLASS}.${ARM}.out"
      [[ -s "$f" ]] || { echo "ERROR: missing expected output: $f" >&2; exit 1; }
    done
  done
else
  if command -v module >/dev/null 2>&1; then
    module load ViennaRNA/2.4.17-gimkl-2020a || true
  fi
  if ! command -v RNAduplex >/dev/null 2>&1; then
    echo "ERROR: RNAduplex not found." >&2
    echo "       Load ViennaRNA module or set SKIP_RDUPLEX=1 and use run_rduplex_aso.sl" >&2
    exit 1
  fi
  echo "RNAduplex: $(which RNAduplex)"

  for CLASS in common_snps pathogenic_missense; do
    for ARM in $ARMS; do
      echo "  ${CLASS}.${ARM} ..."
      RNAduplex < "$RDUPLEX_DIR/${CLASS}.${ARM}.in" \
        > "$RDUPLEX_DIR/${CLASS}.${ARM}.out"
    done
    echo "  ${CLASS}: done"
  done
fi

echo

###############################################################################
# Step 3: parse RNAduplex output and join all scores
###############################################################################

echo "###############################################################################"
echo "Step 3: parse RNAduplex and join"
echo "###############################################################################"

python3 - "$RDUPLEX_DIR" "$SCORE_TSV" "$OUT_TSV" <<'PY'
import sys
import re
import csv
from pathlib import Path

rduplex_dir = Path(sys.argv[1])
score_tsv   = Path(sys.argv[2])
out_tsv     = Path(sys.argv[3])

MFE_RE = re.compile(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)')

def parse_rduplex(path):
    """ΔG per pair; RNAduplex outputs 1 line per pair: structure  (-xx.xx)"""
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    dg = []
    for line in lines:
        m = MFE_RE.search(line)
        dg.append(float(m.group(1)) if m else None)
    return dg

def fmt(v):
    return f'{v:.2f}' if v is not None else 'NA'

def ddg(mis, mat):
    return (mis - mat) if (mis is not None and mat is not None) else None

CLASSES = ['common_snps', 'pathogenic_missense']

DUPLEX_COLS = [
    'duplex_dg_aso_ref_match',   'duplex_dg_aso_ref_mismatch',  'duplex_ddg_aso_ref',
    'duplex_dg_aso_alt_match',   'duplex_dg_aso_alt_mismatch',  'duplex_ddg_aso_alt',
]
MOTIF_COLS = [
    'variant_pos_in_17mer',
    'motif_GGGG', 'motif_AAAA', 'motif_TAAA', 'motif_CTAA', 'motif_CCTA',
    'motif_GCGT', 'motif_TTGT', 'motif_GTAT', 'motif_CGTA', 'motif_GTCG',
]

duplex = {}
motif  = {}

for cls in CLASSES:
    keys    = [l.rstrip('\n') for l in open(rduplex_dir / f'{cls}.duplex_keys.txt') if l.strip()]
    dg_rm   = parse_rduplex(rduplex_dir / f'{cls}.aso_ref.match.out')
    dg_rmis = parse_rduplex(rduplex_dir / f'{cls}.aso_ref.mismatch.out')
    dg_am   = parse_rduplex(rduplex_dir / f'{cls}.aso_alt.match.out')
    dg_amis = parse_rduplex(rduplex_dir / f'{cls}.aso_alt.mismatch.out')

    for lst, name in [(dg_rm,'ref.match'),(dg_rmis,'ref.mismatch'),(dg_am,'alt.match'),(dg_amis,'alt.mismatch')]:
        if len(lst) != len(keys):
            raise ValueError(f"{cls}.{name}: {len(lst)} dg values but {len(keys)} keys")

    for k, rm, rmis, am, amis in zip(keys, dg_rm, dg_rmis, dg_am, dg_amis):
        duplex[k] = {
            'duplex_dg_aso_ref_match':    fmt(rm),
            'duplex_dg_aso_ref_mismatch': fmt(rmis),
            'duplex_ddg_aso_ref':         fmt(ddg(rmis, rm)),
            'duplex_dg_aso_alt_match':    fmt(am),
            'duplex_dg_aso_alt_mismatch': fmt(amis),
            'duplex_ddg_aso_alt':         fmt(ddg(amis, am)),
        }

    with open(rduplex_dir / f'{cls}.motif_scores.tsv') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            motif[row['key']] = {c: row[c] for c in MOTIF_COLS}

    print(f"{cls}: {len(keys)} entries loaded")

# Join to score table
with open(score_tsv) as f:
    reader      = csv.DictReader(f, delimiter='\t')
    score_rows  = list(reader)
    orig_fields = list(reader.fieldnames)

all_fields = orig_fields + DUPLEX_COLS + MOTIF_COLS
missing_d  = 0
missing_m  = 0

with open(out_tsv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=all_fields, delimiter='\t')
    writer.writeheader()
    for row in score_rows:
        k = row['full_header'].split('::')[0]

        d = duplex.get(k)
        if d:
            row.update(d)
        else:
            missing_d += 1
            for c in DUPLEX_COLS:
                row[c] = 'NA'

        m = motif.get(k)
        if m:
            row.update(m)
        else:
            missing_m += 1
            for c in MOTIF_COLS:
                row[c] = 'NA'

        writer.writerow(row)

print(f"\nWrote {out_tsv}")
print(f"Total rows:     {len(score_rows)}")
print(f"Missing duplex: {missing_d}")
print(f"Missing motif:  {missing_m}")
PY

echo

###############################################################################
# Step 4: summary checks
###############################################################################

echo "###############################################################################"
echo "Step 4: checks"
echo "###############################################################################"

echo "Output rows and column count:"
wc -l "$OUT_TSV"
awk -F'\t' '{print NF}' "$OUT_TSV" | sort | uniq -c
echo

echo "ΔΔG summary by variant class and ASO allele:"
python3 - "$OUT_TSV" <<'PY'
import csv, sys
from collections import defaultdict

path = sys.argv[1]
stats = defaultdict(lambda: defaultdict(list))

with open(path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cls = row['variant_class']
        for col in ('duplex_ddg_aso_ref', 'duplex_ddg_aso_alt'):
            v = row.get(col, 'NA')
            if v != 'NA':
                stats[cls][col].append(float(v))

for cls in sorted(stats):
    print(f"  {cls}:")
    for col, vals in sorted(stats[cls].items()):
        mean = sum(vals) / len(vals)
        print(f"    {col}: n={len(vals)}  mean={mean:.2f}  "
              f"range [{min(vals):.2f}, {max(vals):.2f}]")
PY

echo
echo "Motif prevalence by variant class:"
python3 - "$OUT_TSV" <<'PY'
import csv, sys
from collections import defaultdict

path = sys.argv[1]
MOTIFS = ['motif_GGGG','motif_AAAA','motif_TAAA','motif_CTAA','motif_CCTA',
          'motif_GCGT','motif_TTGT','motif_GTAT','motif_CGTA','motif_GTCG']
counts = defaultdict(lambda: defaultdict(int))
totals = defaultdict(int)

with open(path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cls = row['variant_class']
        totals[cls] += 1
        for m in MOTIFS:
            if row.get(m) == '1':
                counts[cls][m] += 1

for cls in sorted(totals):
    n = totals[cls]
    print(f"  {cls} (n={n}):")
    for m in MOTIFS:
        c = counts[cls][m]
        print(f"    {m}: {c} ({100*c/n:.1f}%)")
PY

echo
echo "== Done =="
echo "Extended score table: $OUT_TSV"
echo
echo "NOTE: miRNA binding site overlap not yet scored."
echo "      Run add_mirna_overlap.sh (TargetScan BED + bedtools intersect) to add mirna_overlap column."
