# genomeCoverageBed_tidier.py

#####
import sys
import os
import decimal



args = sys.argv
arg_len = len(args)


if arg_len <2:
	print("genomeCoverageBed_tidier.py")
	
	print("\n**** genomeCoverageBed_tidier.py | Written by DJP, 13/08/16 in Python 3.5 in Lausanne ****\n")
	print("This program takes the default histogram output from genomeCoverageBed (e.g. genomeCoverageBed -ibam mybam.bam -g myref.fa > mybam_coverage.out) ") 
	print("It outputs: mean coverage for each contig and, in a seperate file, the mean coverage for the whole genome")
	print("\n**** USAGE **** \n")
	print("genomeCoverageBed_tidier.py [name of coverage file (e.g. mybam_coverage.out)] [max coverage] [output_basename] [\n")
	print("max coverage\t bases with coverage higher than [max coverage] will be set to the specified value. If no max coverage cutoff is required specify 'max' \n")
	
else:
	IN_bed_file_name = args[1] ###
	max_coverage = args[2] 
	output_file_name = args[3] ### output
	
	print("\nMax coverage used: " + str(max_coverage) + "\n")
	
	by_cont_dict = {}
	done = set()
	bed_file = open(IN_bed_file_name)
	for line in bed_file:
		line = line.rstrip("\n")
		line = line.split("\t")
		cont = line[0]
		depth = int(line[1])
		freq = int(line[2])
		length = int(line[3])
		
		if max_coverage != "max":
			if cont not in done:
				done.add(cont)
				if depth <= int(max_coverage):
					psud_cov = depth * freq ### could divide by length here, but could create issues with rounding from floating points
					by_cont_dict[cont] = [psud_cov, length]
				else:
					psud_cov = int(max_coverage) * freq ### set to max coverage
					by_cont_dict[cont] = [psud_cov, length]
			else:
				if depth <= int(max_coverage):
					psud_cov = depth * freq
					a1a = by_cont_dict.get(cont)
					a1 = int(a1a[0])
					a2 = a1 + psud_cov
					by_cont_dict[cont] = [a2, length]
				else:
					psud_cov = int(max_coverage) * freq
					a1a = by_cont_dict.get(cont)
					a1 = int(a1a[0])
					a2 = a1 + psud_cov
					by_cont_dict[cont] = [a2, length]
					
			
		else: 
			if cont not in done:
				psud_cov = depth * freq
				done.add(cont)
				by_cont_dict[cont] = [psud_cov, length]
			else:
				psud_cov = depth * freq
				a1a = by_cont_dict.get(cont)
				a1 = int(a1a[0])
				a2 = a1 + psud_cov
				by_cont_dict[cont] = [a2, length]
				
				
	
	#print(by_cont_dict)
	
	
	output_file_name_cont = output_file_name + "_contig_cov.txt"
	output_file_name_genome = output_file_name + "_genome_cov.txt"
	
	out_file_C = open(output_file_name_cont, "w")
	out_file_G = open(output_file_name_genome, "w")
	
	out_file_C.write("contig_name" + "," + "mean_coverage" + "\n")

	for el in by_cont_dict:
		a1 = by_cont_dict.get(el)
		tot_C = decimal.Decimal(a1[0])
		lent = decimal.Decimal(a1[1])
		cov = decimal.Decimal(tot_C / lent) ## using decimal module
		if el == "genome":
			out_file_G.write("Genome_coverage = " + "\t" + str(cov) + "\n")
			print("\nMean coverage of all contigs (genome): " + str(cov) + "\n")
		else:
			out_file_C.write(el + "," + str(cov) + "\n")
			
	
	print("\nFinished, Frank Cauldhame\n")
			
			
			
			
			
		

	
