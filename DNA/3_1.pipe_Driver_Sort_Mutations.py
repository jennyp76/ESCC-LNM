#!/usr/bin/python3

""" Sort mutations according to mutation type if gene is listed in the driver list """

# Import Modules
import os
import sys
import glob

# Input Parameters
Input = sys.argv[1]
Force_calling_Path = sys.argv[2]

Name =  Input.rstrip('.vcf')
Driver = open( "TCGA_IntOGen_1.0_CGC_tier1.txt" ,"r")
Data = open("{}/{}".format(Force_calling_Path,Input), 'r')
Potent = open("{}/{}.potent.vcf".format(Force_calling_Path,Name),'w')         
Indel = open("{}/{}.potent.indel.vcf".format(Force_calling_Path,Name),'w')
Missense = open("{}/{}.potent.missense.vcf".format(Force_calling_Path,Name),'w')
Nonsense = open("{}/{}.potent.nonsense.vcf".format(Force_calling_Path,Name),'w')

Database = set()

for line in Driver:
	line = line.rstrip('\n')
	gene = line.split('\t')[0]
	Database.add(gene)

for line in Data:
	line = line.rstrip('\n')
	if line.startswith('##'):
		continue
	elif line.startswith('#C'):
		Potent.write( line + '\n')
		Indel.write( line + '\n')
		Missense.write( line + '\n')
		Nonsense.write( line + '\n')
	else:
		col = line.split('\t')
		ref= col[3]
		alt = col[4]
		info = col[7]
		gene = info.split('|')[3]
		type = info.split('|')[1]

		if gene in Database:
			Potent.write(line +'\n')
			if len(ref) == 1 and len(alt) == 1:
				if type == 'missense_variant':
					Missense.write( line + '\n')
				else:
					Nonsense.write( line + '\n')
			else:
				Indel.write( line + '\n')
		else:
			continue 
 
