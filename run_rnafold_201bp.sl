#!/bin/bash -e
#SBATCH --job-name=rnafold_aso
#SBATCH --account=uoa04097
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/rnafold_aso.%j.out
#SBATCH --error=logs/rnafold_aso.%j.err

module purge
module load ViennaRNA/2.4.17-gimkl-2020a

mkdir -p results/rnafold logs

echo "RNAfold:"
which RNAfold
RNAfold --version

echo "Running pathogenic missense RNAfold..."
RNAfold --noPS < results/rnafold/pathogenic_missense.plusminus100.rna.fa \
  > results/rnafold/pathogenic_missense.plusminus100.rnafold.txt

echo "Running common SNP RNAfold..."
RNAfold --noPS < results/rnafold/common_snps.plusminus100.rna.fa \
  > results/rnafold/common_snps.plusminus100.rnafold.txt

echo "Checking counts..."
grep -c '^>' results/rnafold/pathogenic_missense.plusminus100.rnafold.txt
grep -c '^>' results/rnafold/common_snps.plusminus100.rnafold.txt

echo "Done."
