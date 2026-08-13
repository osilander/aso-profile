#!/bin/bash
#SBATCH --job-name=coverage_curves
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/coverage_curves_%j.out
#SBATCH --error=logs/coverage_curves_%j.err

mkdir -p logs

module load R/4.3.1-gimkl-2022a   # adjust to whatever R module is available
module load BCFtools               # adjust to available module name

cd "$(dirname "$0")"
Rscript coverage_curves.R
