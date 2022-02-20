# Timema_facultative_partheno

This is the repository for the scripts used in this study:
>Larose* C., Lavanchy* G., Freitas S., Parker D.J., Schwander T. TITLE. BioRxiv. DOI

## Demultiplexed read acquisition
The SRA accession numbers for the demultiplexed reads used in this study will be available in Table S1 in the supplementary materials of the article.
They were demultiplexed using the `process_radtags` command of [stacks](https://catchenlab.life.illinois.edu/stacks/) v2.3e with the `-c -q -r -t 92 --filter_illumina`. options.

## Mapping to the reference genome of *Timema douglasi*.
The reference genome used in this study is available under BioSample SAMEA5384775. We used the `mem` algorithm of `bwa` version 0.7.17 and converted to bam format using `samtools` version 1.4. The script `bwa_mem_SLURM_script_maker.sh` was used to produce one SLURM script for each individual.

## SNP calling
We called SNPs using `stacks` v2.3e and the popmap.txt provided on this repository. We ran `gstacks` and `populations` using the following commands:

`gstacks -t 16 --phasing-dont-prune-hets --ignore-pe-reads -M /work/FAC/FBM/DEE/tschwand/default/glavanc1/stacks/popmap.txt -I /work/FAC/FBM/DEE/tschwand/default/glavanc1/mapped -O /work/FAC/FBM/DEE/tschwand/default/glavanc1/stacks`

`populations -t 16 -M /work/FAC/FBM/DEE/tschwand/default/glavanc1/Timema_facultative_partheno/popmap.txt -P /work/FAC/FBM/DEE/tschwand/default/glavanc1/Timema_facultative_partheno/stacks -O /work/FAC/FBM/DEE/tschwand/default/glavanc1/Timema_facultative_partheno/out_1snp --write-single-snp --vcf
`
## Filtering
We tested the effect of different filtering criteria on the number of retained SNPs using the script `filtering_tests.sh` We then visualised the results in `R` using the script `Filtering.R` to produce this graph (where different colours represent different values of minimum depth (`--minDP` in `vcftools`):
### INSERT PLOT

We generated the final vcf with the following command in `vcftools` version 0.1.14:

`vcftools --vcf populations.snps.vcf --minDP 11 --mac 3 --max-missing 0.75 --recode --out DP11_miss75_mac3_1snp_all`

and the subset dataset consisting only of adults with the keep_adults.txt file available on this repository using:

`vcftools --vcf populations.snps.vcf --keep keep_adults.txt --minDP 11 --mac 3 --max-missing 0.75 --recode --out DP11_miss75_mac3_1snp_adults`

## Population structure
We ran `fastSTRUCTURE` version 1.0 using the script `runfaststructure.sh`.

## Identifying the X chromosome of *Timema douglasi*.

* **ID_Tdi_X.sh** | Script to map male (BioProject: PRJNA808673) and female (BioProject: PRJNA670663) whole-genome reads to reference genome.  
* **Tim_cov_analysis.R** | Coverage analyses to ID the X.

## Further analyses
All subsequent analyses were run and plots were produced in `R` version 4.0.3 using the script `Analyses.R`. Plots were edited in Illustrator to improve display.
