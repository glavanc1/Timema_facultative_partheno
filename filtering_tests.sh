#!/bin/bash

#SBATCH --account=tschwand_default
#SBATCH -p cpu
#SBATCH --job-name=test_filtering
#SBATCH --mem=12G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load gcc/9.3.0
module load vcftools/0.1.14

for DP in {3..20}
do
for MISS in $(seq -f "%f" 0.05 0.05 1)
do
vcftools --vcf populations.snps.vcf --minDP $DP --max-missing $MISS --mac 3 --stdout &>> vcftools_all.out
done
done

grep 'minDP' vcftools_all.out | cut -f 2 | cut -f 2 -d ' ' > DP
grep 'max-missing' vcftools_all.out | cut -f 2 | cut -f 2 -d ' ' > miss
grep 'out of a possible' vcftools_all.out | cut -f 4 -d ' ' > kept

paste DP miss kept > summary_filtering.txt

rm DP
rm miss
rm kept
