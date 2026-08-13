#!/usr/bin/env bash
# run_pipeline_af0.1.sh
#
# Re-runs the full ASO pipeline at MAF >= 10% on an HPC cluster (SLURM).
# All outputs land in run_af0.1/ so the 20% results are untouched.
#
# Usage (from the project root on the cluster login node):
#   bash run_pipeline_af0.1.sh

set -euo pipefail

# ── paths ─────────────────────────────────────────────────────────────────────

PROJ=${PROJ:-$(pwd)}           # set PROJ or defaults to current directory
SCRIPTS=${PROJ}
RUNDIR=${PROJ}/run_af0.1
SLDIR=${RUNDIR}/slurm_scripts

ACCOUNT=${ACCOUNT:?Set ACCOUNT to your HPC scheduler account (e.g. export ACCOUNT=myproject)}
MAF=0.1

REF_FASTA=${REF_FASTA:-${PROJ}/refs/GRCh38.primary_assembly.genome.fa}
DATA_DIR=${DATA_DIR:-${PROJ}/data}               # full data dir (VCFs + GTF)
VCF_DIR=${DATA_DIR}/1000G_highcov_GRCh38
GTF=${DATA_DIR}/gencode.v48.basic.annotation.gtf.gz
GENES_BED=${PROJ}/genes.standard.bed
CLINVAR_TSV=${PROJ}/results/clinvar_20genes_pathogenic.header.tsv
PANEL_FILE=${PROJ}/integrated_call_samples_v3.20130502.ALL.panel

# ── set up run directory ──────────────────────────────────────────────────────

echo "Setting up ${RUNDIR} ..."
mkdir -p "${RUNDIR}"/{results,vcf_cache,figures,gene_plot_inputs,logs}
mkdir -p "${SLDIR}"

# Symlink shared resources (data/ is the full directory, not just the VCF subdir)
ln -snf "${PROJ}/refs" "${RUNDIR}/refs"                                                           2>/dev/null || true
ln -snf "${DATA_DIR}"  "${RUNDIR}/data"                                                           2>/dev/null || true
ln -snf "${GENES_BED}" "${RUNDIR}/genes.standard.bed"                                             2>/dev/null || true
ln -snf "${PANEL_FILE}" "${RUNDIR}/integrated_call_samples_v3.20130502.ALL.panel"                 2>/dev/null || true
ln -snf "${CLINVAR_TSV}" "${RUNDIR}/results/clinvar_20genes_pathogenic.header.tsv"                2>/dev/null || true

# ── write SLURM scripts ───────────────────────────────────────────────────────

# Step 1: build variant tables (ClinVar + common SNPs at MAF 0.1)
cat > "${SLDIR}/step1_build.sl" <<SLURM
#!/bin/bash
#SBATCH --job-name=aso_build_af0.1
#SBATCH --account=${ACCOUNT}
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=${RUNDIR}/logs/step1_%j.out
#SBATCH --error=${RUNDIR}/logs/step1_%j.err

set -euo pipefail
cd "${RUNDIR}"
module purge
module load BCFtools BEDTools

bash "${SCRIPTS}/build_aso_variant_tables.sh" \
    "${CLINVAR_TSV}" "${GENES_BED}" "${REF_FASTA}" "${VCF_DIR}" "${MAF}"
SLURM

# Step 2: add CDS/exon/intron context column — produces the .context.tsv that
# score_17mer_aso_windows.sh requires. This step was missing from the original chain.
cat > "${SLDIR}/step2_context.sl" <<SLURM
#!/bin/bash
#SBATCH --job-name=aso_context_af0.1
#SBATCH --account=${ACCOUNT}
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=${RUNDIR}/logs/step2_%j.out
#SBATCH --error=${RUNDIR}/logs/step2_%j.err

set -euo pipefail
cd "${RUNDIR}"
module purge
module load BEDTools Python

MIN_MAF=${MAF} GTF="${GTF}" bash "${SCRIPTS}/make_common_snp_context.sh"
SLURM

# Step 3: score 17-mer ASO windows
# score_17mer_aso_windows.sh hard-codes 'af0.2'; patch with sed | bash
cat > "${SLDIR}/step3_score17mer.sl" <<SLURM
#!/bin/bash
#SBATCH --job-name=aso_score17_af0.1
#SBATCH --account=${ACCOUNT}
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=${RUNDIR}/logs/step3_%j.out
#SBATCH --error=${RUNDIR}/logs/step3_%j.err

set -euo pipefail
cd "${RUNDIR}"
module purge
module load BEDTools

sed 's/af0\.2/af0.1/g' "${SCRIPTS}/score_17mer_aso_windows.sh" \
    | bash -s -- "${REF_FASTA}"
SLURM

# Step 4: RNAfold + join scores
# run_rnafold_and_join_aso_scores.sh also hard-codes 'af0.2'; patch same way
cat > "${SLDIR}/step4_rnafold.sl" <<SLURM
#!/bin/bash
#SBATCH --job-name=aso_rnafold_af0.1
#SBATCH --account=${ACCOUNT}
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=${RUNDIR}/logs/step4_%j.out
#SBATCH --error=${RUNDIR}/logs/step4_%j.err

set -euo pipefail
cd "${RUNDIR}"
module purge
module load ViennaRNA BEDTools

sed 's/af0\.2/af0.1/g' "${SCRIPTS}/run_rnafold_and_join_aso_scores.sh" \
    | bash -s -- "${REF_FASTA}"
SLURM

# Step 5: add extended ASO scores (RNAduplex ΔΔG + 4-mer motifs)
cat > "${SLDIR}/step5_extended.sl" <<SLURM
#!/bin/bash
#SBATCH --job-name=aso_extended_af0.1
#SBATCH --account=${ACCOUNT}
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=${RUNDIR}/logs/step5_%j.out
#SBATCH --error=${RUNDIR}/logs/step5_%j.err

set -euo pipefail
cd "${RUNDIR}"
module purge
module load ViennaRNA Python

bash "${SCRIPTS}/add_extended_aso_scores.sh"
SLURM

# Step 6: coverage_curves.R — query VCFs, build het matrices, greedy curves
# Needs BCFtools because the R script calls bcftools via system()
cat > "${SLDIR}/step6_coverage.sl" <<SLURM
#!/bin/bash
#SBATCH --job-name=aso_coverage_af0.1
#SBATCH --account=${ACCOUNT}
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=${RUNDIR}/logs/step6_%j.out
#SBATCH --error=${RUNDIR}/logs/step6_%j.err

set -euo pipefail
cd "${RUNDIR}"
module purge
module load BCFtools R

Rscript "${SCRIPTS}/coverage_curves.R"
SLURM

# Step 7: plot_ancestry_coverage.R — ancestry-stratified figure
cat > "${SLDIR}/step7_ancestry.sl" <<SLURM
#!/bin/bash
#SBATCH --job-name=aso_ancestry_af0.1
#SBATCH --account=${ACCOUNT}
#SBATCH --time=06:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --output=${RUNDIR}/logs/step7_%j.out
#SBATCH --error=${RUNDIR}/logs/step7_%j.err

set -euo pipefail
cd "${RUNDIR}"
module purge
module load R

Rscript "${SCRIPTS}/plot_ancestry_coverage.R"
SLURM

# ── submit with afterok chain ─────────────────────────────────────────────────

echo "Submitting SLURM jobs ..."

JOB1=$(sbatch --parsable "${SLDIR}/step1_build.sl")
echo "  Job 1 (build):    ${JOB1}"

JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} "${SLDIR}/step2_context.sl")
echo "  Job 2 (context):  ${JOB2}"

JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} "${SLDIR}/step3_score17mer.sl")
echo "  Job 3 (score17):  ${JOB3}"

JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} "${SLDIR}/step4_rnafold.sl")
echo "  Job 4 (rnafold):  ${JOB4}"

JOB5=$(sbatch --parsable --dependency=afterok:${JOB4} "${SLDIR}/step5_extended.sl")
echo "  Job 5 (extended): ${JOB5}"

JOB6=$(sbatch --parsable --dependency=afterok:${JOB5} "${SLDIR}/step6_coverage.sl")
echo "  Job 6 (coverage): ${JOB6}"

JOB7=$(sbatch --parsable --dependency=afterok:${JOB6} "${SLDIR}/step7_ancestry.sl")
echo "  Job 7 (ancestry): ${JOB7}"

echo ""
echo "Dependency chain: ${JOB1} → ${JOB2} → ${JOB3} → ${JOB4} → ${JOB5} → ${JOB6} → ${JOB7}"
echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Logs in:       ${RUNDIR}/logs/"
echo "Final figure:  ${RUNDIR}/figures/ancestry_coverage.pdf"
