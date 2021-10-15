#!/bin/bash
#SBATCH --account=tschwand_default
#SBATCH -p cpu
#SBATCH --job-name=gstacks
#SBATCH --mem=1G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load gcc/9.3.0
module load vcftools/0.1.14
module load stacks/2.53

## Convert the vcf file into structure format
populations -V DP11_miss75_mac3_1snp_adults.recode.vcf -O . --structure
tail -n +3 DP11_miss75_mac3_1snp_adults.recode.p.structure | sed 's/0/-9/g' > DP11_miss75_mac3_1snp_adults.cleaned.str

## run structure with K ranging from 1 to 10
for i in {1..10}
do
echo "K = ${i}"
/software/EcologyEvolution/fastStructure/1.0/bin/structure.py -K $i --input=.cleaned --output=structure_DP11_miss75_mac3_1snp_adults_K --seed=100 --prior=simple --format=str
done

## run chooseK.py to select the best value of K
chooseK.py --input=structure_DP11_miss75_mac3_1snp_adults_K* > distruct_results.txt
