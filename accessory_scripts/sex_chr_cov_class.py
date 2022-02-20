### sex_chr_cov_class.py

import sys
import os
import getopt
import decimal
from decimal import *

try:
	opts, args = getopt.getopt(sys.argv[1:], 'o:r:l:c:p:d:g:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


out_prefix = "NOTHINGSET"
M_F_cov_filename = None
M_F_cov_peak_adjust_filename = None
dist_from_peak = decimal.Decimal(0.1)
getcontext().prec = 7

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** sex_chr_cov_class.py | Written by DJP, 14/01/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("classes scafs as X or A\nSee 4_Sex_chr_expression_analyses.sh for actual usage")
		
		print("\n**** USAGE **** \n")
		print("python3 sex_chr_cov_class.py -c [M/F coverage file] -o [out prefix] [options] \n")
		
		print("\n**** OPTIONS ****\n")
		print("-o\toutput prefix")
		print("-c\tMale/female coverage file - output from male_tim_cov_v7.R")
		print("-p\tpeak adjustment. If want to adjust the value for the peak of X-linked scaffolds, specify the peak value file (from output from male_tim_cov_v7.R) here. Default: OFF")
		print("-d\tdistance around X-chr cov peak to keep - Default: 0.1")

		sys.exit(2)

	elif opt in ('-o'):
		out_prefix = arg
	elif opt in ('-c'):
		M_F_cov_filename = arg
	elif opt in ('-p'):
		M_F_cov_peak_adjust_filename = arg
	elif opt in ('-d'):
		dist_from_peak = decimal.Decimal(arg)
	else:
		print("i dont know")
		sys.exit(2)


################################################################################################################################################
######## peak adjustment

Sex_chr_upper = 0
Sex_chr_lower = 0

if M_F_cov_peak_adjust_filename != None:
	print("\nDoing peak adjustment\n")
	#out_file_stats.write("\nDoing peak adjustment\n")
	species_using = M_F_cov_filename.rstrip("/").split("/")[-1].split("_")[0]
	
	print("Peak value using:")
	#out_file_stats.write("Peak value using: ")
	peak_val = 0
	
	M_F_cov_peak_adjust_file = open(M_F_cov_peak_adjust_filename)
	for line in M_F_cov_peak_adjust_file:
		line = line.replace('"', '').rstrip("\n")
		
		if line.startswith(species_using):
			print(line)
			#out_file_stats.write(line + "\n")
			peak_val=decimal.Decimal(line.split(",")[1])
		
	M_F_cov_peak_adjust_file.close()
	
	Sex_chr_upper = decimal.Decimal(peak_val) + decimal.Decimal(dist_from_peak)
	Sex_chr_lower = decimal.Decimal(peak_val) - decimal.Decimal(dist_from_peak)	
		
else:
	print("\nNo peak adjustment\n")
	#out_file_stats.write("\nNo peak adjustment\n")
	Sex_chr_upper = decimal.Decimal(-1) + decimal.Decimal(dist_from_peak)
	Sex_chr_lower = decimal.Decimal(-1) - decimal.Decimal(dist_from_peak)

	
print("X chromosome contigs classed as having log2 M/F ratios between: " + str(Sex_chr_lower) + " and " + str(Sex_chr_upper) + "\n")
# out_file_stats.write("X chromosome contigs classed as having log2 M/F ratios between: " + str(Sex_chr_lower) + " and " + str(Sex_chr_upper) + "\n")


################################################################################################################################################
######## read M/F coverage in and class scaffs as X or A 

print(peak_val)

M_F_cov_file = open(M_F_cov_filename)

M_F_dict = {}


A_scafs_N = 0
X_scafs_N = 0
A_len = 0
X_len = 0

line_N = 0
for line in M_F_cov_file:
	line_N = line_N + 1
	if line_N > 1:
		line = line.rstrip("\n").replace('"', '').split(",")
		contig_name = line[0]
		M_F = decimal.Decimal(line[5])
		M_F_peak_adj = M_F + (-1 - peak_val)
		chr_len = int(line[6])
		chr_type = ""
		
		if Sex_chr_lower <= M_F <= Sex_chr_upper:
			chr_type = "X"
			X_len = X_len + chr_len
			X_scafs_N = X_scafs_N + 1
			# print(line)
			# print(M_F)
		else:
			chr_type = "A"
			A_len = A_len + chr_len
			A_scafs_N = A_scafs_N + 1
		
		M_F_dict[contig_name] = [chr_type, chr_len, M_F, M_F_peak_adj]

M_F_cov_file.close()


print("N Autosomal scafs: " + str(A_scafs_N))
print("N X-linked scafs: " + str(X_scafs_N))
print("Number of bases Autosomal: " + str(A_len))
print("Number of bases X-linked:  " + str(X_len))
print("Percentage of genome X-linked: " + str((X_len / (X_len + A_len)) * 100))

### output

out_file = open(out_prefix + "_chr_class.csv", "w")

out_file.write("scaf_name" + "," + "chr_type,MF,MF_adj" + "\n")

for el in M_F_dict:
	chr_type = M_F_dict.get(el)[0]
	MF_adj = M_F_dict.get(el)[3]
	MF = M_F_dict.get(el)[2]
	out_file.write(el + "," + chr_type +  "," + str(MF) + "," +  str(MF_adj) +  "\n")
	




# 
print("\n\nFinished, Walker.\n\n")	
# 		