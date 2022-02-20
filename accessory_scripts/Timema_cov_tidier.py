import sys
import os
import getopt
import re
import decimal

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:f:m:e:s:l:rPh')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_dir_name = None
seqF1 = None
min_seq_len = 0
out_file_name_base = "COV_OUT"
cov_ext = "30a_contig_cov.txt"
paired_id = "_pe_"
single_id = "_se_"
unpaired_id = "_UNPAIRED_"
paired_only   = False

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** Timema_cov_tidier.py | Written by DJP, 29/10/19 in Python 3.5 in Lausanne, Swiss ****\n")
		print("\n**** USAGE example ****\n")
		print("python3 Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out/Tbi_cont_cov -f 4_Tbi_b3v08.fasta -m 1000 -o Tbi\n\n\n")

		sys.exit(2)
	
	elif opt in ('-i'):
		in_dir_name = arg
	elif opt in ('-o'):
		out_file_name_base = arg
	elif opt in ('-f'):
		seqF1 = arg
	elif opt in ('-m'):
		min_seq_len = arg
	elif opt in ('-l'):
		min_len = int(arg)
	elif opt in ('-e'):
		cov_ext = arg
	elif opt in ('-P'):
		paired_only  = True
	else:
		print("i dont know")
		sys.exit(2)





#############################################################################################################################################
#### first get scaf sizes into a dict


##### FIRST unwrap fasta - precautionary will be necessary for some files 
### note making a temp unwrapped fasta file  - removed at end
output_fasta_name = seqF1 + ".TEMP_extract_fasta_file" 

output_file = open(output_fasta_name, "w")
count = 0
in_file = open(seqF1)
for line in in_file:
	count = count + 1
	line = line.rstrip("\n")
	if line.startswith(">") and count == 1:
		output_file.write(line + "\n")
	elif line.startswith(">") and count > 1:
		output_file.write("\n" + line + "\n")
	else: 
		output_file.write(line)	

output_file.close()


seq_len_dict = {}

### add seqs to dictionary
name_list = []
seq_list = []
seq_dict = {}

done = 0
seq_file_1 = open(output_fasta_name)
for line in seq_file_1:
	lineA = line.rstrip("\n")
	if lineA.startswith(">"):
		lineB = lineA.replace(">", "")
		name_list.append(lineB)
	else:
		seq_list.append(lineA)
		done = done + 1

for element in range(0,len(name_list)):
	name1 = name_list[element]
	seq1 = seq_list[element]
	seq_dict[name1] = seq1


for el in name_list:
	a1 = len(seq_dict.get(el))
	seq_len_dict[el] = a1 

## tidyup
seq_file_1.close()
os.remove(output_fasta_name)

########################################################################################################
### sum coverage

print("min seq length: " + str(min_seq_len))

single_cov_dict = {}
unpaired_cov_dict = {}
paired_cov_dict = {}

sample_set = set()

paired_sample_set = set()
single_sample_set = set()
unpaired_sample_set = set()

all_contigs_set = set()

N_files_paired = 0
N_files_single = 0
N_files_unpaired = 0

path = in_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith(cov_ext):
			#print (os.path.join(path, name))
			#print(name)
			sample = re.split(r'_G[a-zA-Z0-9]*L[0-9]_', name)[0]
			sample_set.add(sample)
			
			sample_type = None
			if paired_id in name:
				sample_type = "paired"
				N_files_paired = N_files_paired + 1
			elif single_id in name:
				if unpaired_id in name:
					sample_type = "unpaired"
					N_files_unpaired = N_files_unpaired + 1
				else:
					sample_type = "single"
					N_files_single = N_files_single + 1
			else:
				sample_type = "ERROR"
				print("ERROR")
				sys.exit(2)
			
			#print(sample_type)
			
			line_N  = 0
			curr_file = open(os.path.join(path, name))
			for line in curr_file:
				line_N  = line_N  + 1
				if line_N > 1:
					line = line.rstrip().split(",")
					scaf_ID = line[0]
					all_contigs_set.add(scaf_ID)
					cov = decimal.Decimal(line[1])
					
					scaf_sample = (scaf_ID, sample)
					
					if sample_type == "paired":
						
						if scaf_sample not in paired_sample_set:
							paired_sample_set.add(scaf_sample)
							paired_cov_dict[scaf_sample] = cov
						else:
							rec = paired_cov_dict.get(scaf_sample)
							new_rec = rec + cov
							paired_cov_dict[scaf_sample] = new_rec
							
					
					elif sample_type == "single":
						
						if scaf_sample not in single_sample_set:
							single_sample_set.add(scaf_sample)
							single_cov_dict[scaf_sample] = cov
						else:
							rec = single_cov_dict.get(scaf_sample)
							new_rec = rec + cov
							single_cov_dict[scaf_sample] = new_rec
							
					elif sample_type == "unpaired":
						if scaf_sample not in unpaired_sample_set:
							unpaired_sample_set.add(scaf_sample)
							unpaired_cov_dict[scaf_sample] = cov
						else:
							rec = unpaired_cov_dict.get(scaf_sample)
							new_rec = rec + cov
							unpaired_cov_dict[scaf_sample] = new_rec														
					else:
						print("ERROR")
						sys.exit(2)					
							
					

all_contigs_list = list(all_contigs_set)
all_contigs_list_sorted = sorted(all_contigs_list)

sample_list = list(sample_set)
sample_list_sorted = sorted(sample_list)

print("\nSample names found:")
print(sample_list_sorted)		

print("\nNumber of paired files used:" + str(N_files_paired))
print("Number of single files used:" + str(N_files_single))
print("Number of unpaired files used:" + str(N_files_unpaired))

print("**Note** 'single files' = files mapped single end, excluding unpaired reads from trimming")

### output

out_pe_filename    = out_file_name_base + "_pairedcov_minlen="            + str(min_seq_len) + "_contig_cov.txt"
out_peUP_filename  = out_file_name_base + "_pairedandunpairedcov_minlen=" + str(min_seq_len) + "_contig_cov.txt"
out_seUP_filename  = out_file_name_base + "_singleandunpairedcov_minlen=" + str(min_seq_len) + "_contig_cov.txt"

out_pe_file    = open(out_pe_filename,   "w")
if paired_only != True:	
	out_peUP_file  = open(out_peUP_filename, "w")
	out_seUP_file  = open(out_seUP_filename, "w")

header = "contig,length"

for el in sample_list_sorted:
	header = header + "," + el

out_pe_file.write(header + "\n")
if paired_only != True:	
	out_peUP_file.write(header + "\n")
	out_seUP_file.write(header + "\n")

for cont in all_contigs_list_sorted:
	
	cont_len = seq_len_dict.get(cont)
	if cont_len >= int(min_seq_len):
		out_line_paired              = ""
		out_line_paired_and_unpaired = ""
		out_line_single_and_unpaired = ""
		
		for s in sample_list_sorted:
			rec_paired   = paired_cov_dict.get((cont,s))
			if paired_only != True:		
				rec_single   = single_cov_dict.get((cont,s))
				rec_unpaired = unpaired_cov_dict.get((cont,s))
			
			paired_cov              = rec_paired
			
			if paired_only != True:					
				paired_cov_and_unpaired = rec_paired + rec_unpaired
				single_cov_and_unpaired = rec_single + rec_unpaired
			
			out_line_paired = out_line_paired + "," + str(rec_paired)
			if paired_only != True:		
				out_line_paired_and_unpaired = out_line_paired_and_unpaired + "," + str(paired_cov_and_unpaired)
				out_line_single_and_unpaired = out_line_single_and_unpaired + "," + str(single_cov_and_unpaired)	
			
		out_pe_file.write(cont + "," + str(cont_len) + out_line_paired + "\n")
		if paired_only != True:		
			out_peUP_file.write(cont + "," + str(cont_len) + out_line_paired_and_unpaired+ "\n")
			out_seUP_file.write(cont + "," + str(cont_len) + out_line_single_and_unpaired + "\n")
		
	

print("\n\nFinished, Victor Frankenstein\n\n\n")







