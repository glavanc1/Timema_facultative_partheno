# ID_Tdi_X.sh
# I assume all code is run from the github repo directory

DIR='Timema_facultative_partheno'

#########################################################################################################################################################
#### download reference Tdi genome (accession number: PRJEB31411)

mkdir -p $DIR/Genomes/REFS/
cd $DIR/Genomes/REFS/
wget https://zenodo.org/record/5636226/files/Tdi_b3v08.fasta


#########################################################################################################################################################
#### download Tdi reads for males and females from NCBI's SRA to into $DIR/READS/rawreads/ 
#### Female reads from bioproject PRJNA670663

## ReSeq_Di02 
1_Tdi
Tdi_01
SRS7637469
SRR12928239, SRR12928250, SRR12928261, SRR12928272, SRR12928283, SRR12928294, SRR12928305, SRR12928961, SRR12928972, SRR12928983, SRR12929022, SRR12929034, SRR12929045
ReSeq_Di04
1_Tdi
Tdi_02
SRS7637489
SRR12928865-SRR12928870, SRR12928872, SRR12928878, SRR12928889, SRR12928928, SRR12928939, SRR12928950
ReSeq_Di06
1_Tdi
Tdi_03
SRS7637497
SRR12928806-SRR12928814, SRR12928861-SRR12928864
ReSeq_Di08
1_Tdi
Tdi_04
SRS7638280
SRR12928792, SRR12928794-SRR12928803, SRR12928805
ReSeq_Di10
1_Tdi
Tdi_05
SRS7638278
SRR12928779-SRR12928781, SRR12928783-SRR12928791


## Male reads from bioproject PRJNA808673






#########################################################################################################################################################
####  trim reads with trimmomatic v. 0.36

mkdir -p $DIR/READS/trimmed_reads_by_RG
cd       $DIR/READS/trimmed_reads_by_RG

cp       $DIR/data/AllIllumina-PEadapters.fa .

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/trimmomatic/0.36

for i in $DIR/READS/rawreads/*_R1_CC.fq.gz ; do
        foo1=`echo $i`
		basename=`echo $foo1 | sed 's/_R1_CC.fq.gz.*//' | sed 's/.*\///'`
        infileR1=`echo $foo1`
        infileR2=`echo $foo1 | sed 's/_R1_CC.fq.gz/_R2_CC.fq.gz/'`
        outfileR1=`echo "./"$basename"_R1_qtrimmed.fq.gz"`
        outfileR2=`echo "./"$basename"_R2_qtrimmed.fq.gz"`
        outfileR1_UP=`echo "./"$basename"_R1_qtrimmed_UNPAIRED.fq.gz"`
        outfileR2_UP=`echo "./"$basename"_R2_qtrimmed_UNPAIRED.fq.gz"`
        
        echo $infileR1
        echo $infileR2
        echo $outfileR1
        echo $outfileR1_UP
        echo $outfileR2
        echo $outfileR2_UP
		
        trimmomatic PE -threads 30 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:90
done



#########################################################################################################################################################
####  mapping to v8 genome with BWA v. 0.7.15

### prep ref

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3

for i in $DIR/Genomes/REFS/*fasta; do
	bwa index $i
done



#########################################################################################################################################################
#### Mapping stratagy 

## I want bams for coverage analysis

### Mapping as paired-end
#### remove duplicate reads, filter by mapq, and make sorted bam
### BWA MEM does not give an XT:A:U tag in MEM, but XA tags are still there.
### HOWEVER XA tags are only there if there are 1-3 alternative alignments.
### However, this should not be a problem since multiple alignments would result in a low mapping_v8 quality (and can thus be removed on this basis).
# A simple grep would work, but that could leave me with broken flags (filter out one read in a pair and not the other without adjusting the flag).
# SO need to ID multi reads the remove them (I do not need to do this for single end)
# note ONLY WORKS with this samtools 1.3 or higher

# I also filter out Supplemental alignments. These are chimeric alignments
# A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
# So it is a bit of a multimapper. Some of these supplemental alignments have very high mapQ values (60)
# They have an SA:Z tag. If these are paired reads I also remove the other pair, as I do with multi-mappers


# I also add read groups
# to ADD read groups - easiest at BWA stage - but can add / change with picard if I need to
#ID = Read group identifier This tag identifies which read group each read belongs to, so each read group's ID must be unique.
#PU = Platform Unit The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}. 
#SM = Sample The name of the sample sequenced in this read group. 
#PL = Platform/technology used to produce the read. Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.
#LB = DNA preparation library identifier 

# for BWA
# -R "@RG\tID:S1L6\tSM:S1\tPL:ILLUMINA\tLB:FC-140-1086"
# in my case lib = sample

# I then merge bams together by sample and remove PCR duplicates
 # With Picard
 
### I then do coverage analysis on these BAMs


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### map reads as paired reads with BWA

mkdir -p $DIR/mapping_v8/BWA_out/
mkdir -p $DIR/mapping_v8/BWA_out/mapped_as_paired
mkdir -p $DIR/mapping_v8/BWA_out/flagstat_out_paired

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3

read_dir=$DIR"/READS/trimmed_reads_by_RG"
ref_dir=$DIR"/Genomes/REFS/"
map_out_dir=$DIR"mapping_v8/BWA_out/mapped_as_paired"
flag_out_dir=$DIR"mapping_v8/BWA_out/flagstat_out_paired"
mapqfilt="30"

for i in $read_dir/*_R1_qtrimmed.fq.gz; do
		read_file_name_R1=`echo $i`
		read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
        base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
        base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
		badnames=`echo $base_read_name"_badnames.txt"`
        infile=`echo $i`
		sp=`echo $base_read_name | sed 's/_.*//'`
        ref_fa=`echo $ref_dir"/"$sp"_b3v08.fasta"`
        outsam=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA.sam"`
		outbam=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA_mapqfilt_"$mapqfilt".bam"`
		outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA_mapqfilt_"$mapqfilt"_sorted.bam"`
		flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA_flagstat_out.txt"`
		flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA_mapqfilt_"$mapqfilt"_flagstat_out.txt"`
		IFS='_' read -r -a sp_want_list <<< "$base_read_name"
		readgroup=`echo ${sp_want_list[-1]}`
		readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
		echo $read_file_name_R1
		echo $read_file_name_R2
        echo $base_read_name
		echo $base_read_name3
        echo $ref_fa
        echo $outsam
		echo $outbam
		echo $outbam_sorted
		echo $flagstat_out_sam
		echo $flagstat_out_bam
		echo $readgroup
		echo $readgroup_txt
		echo $badnames
		echo ""

		## map
        bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
		
		#flagstat
		samtools flagstat $outsam > $flagstat_out_sam
		
		# filter ## filter both reads out to avoid broken flags
		samtools view -S  $outsam | fgrep XA | cut -f 1 > $badnames
		samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
		# sort bam
		samtools sort $outbam > $outbam_sorted
		
		#flagstat
		samtools flagstat $outbam_sorted > $flagstat_out_bam
		
		#tidyup
		rm $outsam
		rm $outbam
		
done

rm *_badnames.txt


#########################################################################################################################
##### remove supp reads

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3


for i in $DIR/mapping_v8/BWA_out/mapped_as_paired/*_mapqfilt_30_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30_sorted.bam/_mapqfilt_30a_sorted.bam/'`
    badnames=`echo $i"_badnames.txt"`
    flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`

	echo $i
	echo $outbam
	echo $flagstat_out_bam
	
    samtools view $i | fgrep SA:Z: | cut -f 1 > $badnames
    samtools view -h $i | fgrep -vf $badnames | samtools view -b > $outbam
    samtools flagstat $outbam > $flagstat_out_bam
done

## remove original bams and _badnames.txt

rm $DIR/mapping_v8/BWA_out/mapped_as_paired/*badnames.txt
rm $DIR/mapping_v8/BWA_out/mapped_as_paired/*_mapqfilt_30_sorted.bam


#########################################################################################################################
##### merge bams

mkdir -p $DIR/mapping_v8/BWA_out/mapped_as_paired_merged
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tdi_F_ReSeq_Di02_to_Tdi_v8_pe_BWA_mapqfilt_30_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tdi_F_ReSeq_Di02*_mapqfilt_30_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tdi_F_ReSeq_Di04_to_Tdi_v8_pe_BWA_mapqfilt_30_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tdi_F_ReSeq_Di04*_mapqfilt_30_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tdi_F_ReSeq_Di06_to_Tdi_v8_pe_BWA_mapqfilt_30_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tdi_F_ReSeq_Di06*_mapqfilt_30_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tdi_F_ReSeq_Di08_to_Tdi_v8_pe_BWA_mapqfilt_30_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tdi_F_ReSeq_Di08*_mapqfilt_30_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tdi_F_ReSeq_Di10_to_Tdi_v8_pe_BWA_mapqfilt_30_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tdi_F_ReSeq_Di10*_mapqfilt_30_sorted.bam

Tdi_M_18-3997_to_Tdi_v8_pe_BWA_mapqfilt_30aDR_sorted.bam
Tdi_M_18-3998_to_Tdi_v8_pe_BWA_mapqfilt_30aDR_sorted.bam



#########################################################################################################################
##### remove PCR duplicates with picard 2.9.0 

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/picard-tools/2.9.0 
module load UHTS/Analysis/samtools/1.3

for i in $DIR/mapping_v8/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30a_sorted.bam/_mapqfilt_30aDR_sorted.bam/'`
	flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
	metric_file=`echo $outbam | sed 's/.bam/_metric.txt/'`

	echo $i
	echo $outbam
	echo $metric_file
	echo $flagstat_out_bam

	picard-tools MarkDuplicates REMOVE_DUPLICATES=true \
	INPUT=$i \
    OUTPUT=$outbam \
    METRICS_FILE=$metric_file
	
	samtools flagstat $outbam > $flagstat_out_bam

done

## tidy - remove original bams

rm $DIR/mapping_v8/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam 
rm $DIR/mapping_v8/BWA_out/mapped_as_paired_merged/*_metric.txt



############################################################################################
###### calc coverage

module add  Bioinformatics/Software/vital-it
module load UHTS/Analysis/BEDTools/2.26.0
full_ref_dir=$DIR/"Genomes/REFS"

for i in $DIR/mapping_v8/BWA_out/mapped_as_paired_merged/*_sorted.bam; do
    sp=`echo $i | sed 's/.*\///' | sed 's/_.*//'`
	inname=`echo $i`
    basename=`echo $i | sed 's/_sorted.bam//'`
    out_file=`echo $basename"_coverage.out"`
    ref_fa=`echo $full_ref_dir"/"$sp"_b3v08.fasta"`
    echo $sp
	echo $inname
    echo $basename
    echo $out_file
    echo $ref_fa
	echo ""
    genomeCoverageBed -ibam $inname -g $ref_fa > $out_file
	
done

mkdir $DIR/mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR
mv    $DIR/mapping_v8/BWA_out/*_merged*/*_coverage.out $DIR/mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR



#######################################################################################################################
#######################################################################################################################
## ALSO Calc cov for coverage analyses for scafs

mkdir $DIR/mapping_v8/mappingcoverage_ests_BWA_out_aDR
cp    $DIRmapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/*_coverage.out $DIR/mapping_v8/mappingcoverage_ests_BWA_out_aDR ## moving to sep folder for ease

for i in ./mapping_v8/mappingcoverage_ests_BWA_out_aDR/*_coverage.out; do
	in_name=`echo $i`
    basename=`echo $i | sed 's/_coverage.out//'`	
	echo $in_name
    echo $basename
	python ~/Timema_sex_chr_evol_code/accessory_scripts/genomeCoverageBed_tidier.py $in_name max $basename
done


    
#######################################################################################################################
## join together by sp with contig lengths

cd $DIR/mapping_v8/mappingcoverage_ests_BWA_out_aDR

## Filter down to smallest contig being 1000 bp 
python3 $DIR/accessory_scripts/Timema_cov_tidier.py -i $DIR/mapping_v8/mappingcoverage_ests_BWA_out_aDR/ -f $DIR/Genomes/REFS/Tdi_b3v08.fasta -m 1000 -o Tdi -e 30aDR_contig_cov.txt -P

mkdir $DIR/mapping_v8/v8aDR_cov_contig
mv    *_contig_cov.txt $DIR/mapping_v8/v8aDR_cov_contig

## This coverage file (Tdi_pairedcov_minlen=1000_contig_cov.txt) is provided in data/coverage for convenience 

#######################################################################################################################
### Then analyse coverage with Tim_cov_analysis.R
## stored the output (Tdi_v8_MFcov_filt_1000_pairedcov_.csv, sex_chr_peaks_1000_pairedcov_.csv) in data/coverage for convenience 


#######################################################################################################################
### Then class scafs as X linked
python3 accessory_scripts/sex_chr_cov_class.py  -c data/coverage/Tdi_v8_MFcov_filt_1000_pairedcov_.csv -p data/coverage/sex_chr_peaks_1000_pairedcov_.csv -o data/coverage/Tdi_1000_peakadj_2

## stored the output (Tdi_1000_peakadj_2) in data/coverage for convenience 
