#!/usr/bin/python3

""" Collect all indels detected by Strelka Calling within each patients """

# Import Modules
import os
import glob
import sys

# Input Parameters
PATIENT = sys.argv[1]
Strelka_Path=sys.argv[2]

files = glob.glob("{}/{}_*indels.PASS.vcf".format(Strelka_Path,PATIENT))

def sort_chr(chr) :
        if chr[3:] == 'X':
                return 23
        elif chr[3:] == 'Y':
                return 24
        elif chr[3:] == 'M':
                return 0
        else :
                return int(chr[3:])

def sort_Chr(dict) :
        sorted_key = sorted(dict.keys() , key=lambda x : int(x.split('\t')[1]))
        sorted_key = sorted(sorted_key, key=lambda x : sort_chr(x.split('\t')[0]))
        return sorted_key

dict = {}
name_list =[]

for file in files:
	name = os.path.basename(file).split('.')[0]
	name_list.append(name)

	data = open( file, 'r')
	for line in data:
		line = line.rstrip('\n')
		if line.startswith('#'):
			continue
		else:
			col = line.split('\t')
			chr = col[0]
			pos = col[1]
			ref = col[3]
			alt = col[4]
			key = chr + '\t'+ pos + '\t' + ref + '\t' + alt
	
			if key not in dict.keys():
				dict[key]=[]
			else:
				pass
	
			dict[key].append(name)

sorted_key  = sort_Chr(dict)
output = open("{}/{}.somatic.indels.PASS.total.txt".format(Strelka_Path,PATIENT),"w")
output.write('chr\tpos\tref\talt\t'+ '\t'.join(name_list)+'\n')

for i in sorted_key :
	output_list = [] 
	for n in name_list:
		if n in dict[i]:
			output_list.append('O')
		else:
			output_list.append('-')

	output.write( i + '\t'+ '\t'.join(output_list)+ '\n' )

