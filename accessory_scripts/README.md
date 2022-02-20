# Timema_sex_chr_evol | accessory_scripts

* **genomeCoverageBed_tidier_wholegenomecov.py** | Tidy output of genomeCoverageBed for whole genome - see *2_genome_read_mapping.sh*
* **genomeCoverageBed_tidier_select_scafs.py** | Tidy output of genomeCoverageBed for selected scaffolds - see *2_genome_read_mapping.sh*
* **genomeCoverageBed_tidier.py** | Tidy output of genomeCoverageBed per scaffolds - see *2_genome_read_mapping.sh*
* **plot_genome_cov.R** | plot genome coverage - see *2_genome_read_mapping.sh*
* **Timema_cov_tidier.py**  | join coverage estimates together for each sample - see *2_genome_read_mapping.sh*
* **class_scafs_to_lg.py**  | add linkage group infomation to coverage table - see *2_genome_read_mapping.sh*
* **fasta_select_by_len.py** | select fasta sequences by length - see *2_genome_read_mapping.sh*
* **add_hetero_info.py** | add angsD outputs to coverage info - see *2_genome_read_mapping.sh*
* **Maker_gff_to_HTseq_gff.py** | Convert downloaded gff from maker to a form HTseq-count can use - see *3_RNAseq_mapping.sh*
* **HTSeq_to_edgeR.py** | Join read counts from HTseq-count for each sample together - see *3_RNAseq_mapping.sh*
* **gff_feature_lengths.py**  | get feature lengths from a gff - see *3_RNAseq_mapping.sh*
* **class_genes_to_lg.py**  | add linkage group infomation to read count table - see *3_RNAseq_mapping.sh*
* **sex_chr_cov_readcounts_tidier.py** | add ortholog infomation to read count table - see *3_RNAseq_mapping.sh*
* **sex_chr_cov_readcounts_orth_tidier.py** | get read counts for orthologs in each species - see *3_RNAseq_mapping.sh*
* **add_gene_name_to_sel_data.py** | add gene name to selection data - see *3_RNAseq_mapping.sh*
* **alignments_to_seq.py** | turn alignments into sequences - see  *4_selection_analyses.sh*
* **add_GC_contigs.py** | add GC ests to coverage file - see  *4_selection_analyses.sh*
