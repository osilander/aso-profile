#!/usr/bin/env bash
# Build bowtie indices for off-target search.
# Run once on an HPC cluster (or locally if you have bowtie + the reference files).
#
# module load Bowtie2/2.5.1-GCC-11.3.0 BEDTools/2.31.0-GCC-11.3.0   (or similar)
# Locally:  conda install -c bioconda bowtie2 bedtools
module load Bowtie2 2>/dev/null || true
module load BEDTools 2>/dev/null || true
#
# Reference files needed:
#   TRANSCRIPTOME : GENCODE v48 cDNA + lncRNA (spliced)
#   PREMRNA       : gene-body sequences (extracted from genome using GTF)
#   GENOME        : GRCh38 genome FASTA (optional; large)

set -euo pipefail

OUTD="data/bowtie_indices"
mkdir -p "${OUTD}"

# ------------------------------------------------------------------ #
# 1. Spliced transcriptome (GENCODE v48)
# ------------------------------------------------------------------ #
TRANSCRIPT_FA="${OUTD}/gencode.v48.transcripts.fa"
if [ ! -f "${TRANSCRIPT_FA}" ]; then
    echo "Downloading GENCODE v48 transcriptome..."
    wget -q -O "${TRANSCRIPT_FA}.gz" \
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz"
    gunzip "${TRANSCRIPT_FA}.gz"
fi
echo "Building transcriptome bowtie index..."
bowtie2-build --threads 4 "${TRANSCRIPT_FA}" "${OUTD}/transcriptome"
echo "Done: ${OUTD}/transcriptome.*"

# ------------------------------------------------------------------ #
# 2. Pre-mRNA gene bodies (genome FASTA + GTF → gene-body FASTA)
# ------------------------------------------------------------------ #
# Requires: bedtools, genome FASTA, GTF
# Adjust GENOME_FA to your local path or HPC filesystem.
# Typical HPC example path:
# /path/to/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
GENOME_FA="${GENOME_FA:-data/GRCh38.primary_assembly.genome.fa}"
GTF="data/gencode.v48.basic.annotation.gtf.gz"
GENE_BED="${OUTD}/gene_bodies.bed"
GENE_FA="${OUTD}/gene_bodies.fa"

if [ -f "${GENOME_FA}" ]; then
    echo "Extracting gene body sequences..."
    # Extract gene-level intervals (+ 500 bp padding each side for splice site context)
    zcat "${GTF}" \
      | awk '$3 == "gene"' \
      | awk '{
          split($0, a, "\t");
          chrom=a[1]; start=a[4]-1; end=a[5]; strand=a[7];
          start=(start-500<0)?0:start-500; end=end+500;
          # gene_name from attributes
          match($0, /gene_name "([^"]+)"/, m);
          name=m[1];
          print chrom "\t" start "\t" end "\t" name "\t0\t" strand
        }' > "${GENE_BED}"
    bedtools getfasta -fi "${GENOME_FA}" -bed "${GENE_BED}" \
        -s -name -fo "${GENE_FA}"
    bowtie2-build --threads 4 "${GENE_FA}" "${OUTD}/premrna"
    echo "Done: ${OUTD}/premrna.*"
else
    echo "GENOME_FA not found (${GENOME_FA}); skipping pre-mRNA index."
    echo "Set GENOME_FA=/path/to/genome.fa and re-run to build it."
fi

# ------------------------------------------------------------------ #
# 3. Whole genome (optional; large — ~3 GB index)
# ------------------------------------------------------------------ #
if [ "${BUILD_GENOME_INDEX:-0}" = "1" ] && [ -f "${GENOME_FA}" ]; then
    echo "Building whole-genome bowtie index (slow, ~30 min)..."
    bowtie2-build --threads 8 "${GENOME_FA}" "${OUTD}/genome"
    echo "Done: ${OUTD}/genome.*"
fi

echo "All indices built in ${OUTD}/"
