#!/usr/bin/python3

""" Additional prediction for missense mutations. 
Functionally damaging if 2 or more methods consider the mutation deleterious """


# Import Modules
import os
import sys
import glob

# Input Parameters
Driver_Path = sys.argv[1]
Sample = sys.argv[2]
Name = sys.argv[3]

MUT = set()
INFO = {}

Anno = open( "{}/{}.missense.vcf".format(Driver_Path, Name), "r" )
for line in Anno:
	line = line.rstrip('\n')
	if line.startswith('#'):
		continue
	else:
		col = line.split('\t')
		chr = col[0]
		pos = col[1]
		id = col[2]
		ref = col[3]
		alt = col[4]
		info = col[5:]

		key = chr +'\t'+ pos + '\t' + id + '\t' + ref + '\t' + alt 
		INFO[key] = info

Data = open( "{}/{}.missense.dbnsfp.txt".format(Driver_Path,Name) , 'r')
Output= open("{}/{}.missense.dbnsfp.damaging.vcf".format(Driver_Path,Name),'w')

for line in Data:
        line = line.rstrip('\n')
        col = line.split('\t')
        if line.startswith('Chr'):
                output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'+ Sample + '\n')
        else:
                chr = 'chr' + col[0]
                pos = col[1]
                ref = col[3]
                alt = col[4]
		id = '.'
		key = chr +'\t'+ pos + '\t' + id + '\t' + ref + '\t' + alt
		
                sift = col[5]
		provean = col[10]
                revel = col[13]
		score = [ sift, provean ]
		cnt = 0 

                for i in score:
                        if 'D' in i:
				cnt = cnt +1
			else:
				continue

		if float(revel) > 0.5:
			cnt = cnt + 1
		else:
			pass
	
		if cnt >= 2 :
			MUT.add(key)
		else:
			continue

def sort_chr(chr):
        if chr[3:] == 'X':
                return 23
        elif chr[3:] == 'Y':
                return 24
        elif chr[3:] == 'M':
                return 0
        else :
                return int(chr[3:])

MUT = sorted(MUT, key=lambda x : int(x.split('\t')[1]))
MUT = sorted(MUT, key=lambda x : sort_chr(x.split('\t')[0]))

for i in MUT:
	if i in INFO:
		Output.write( i + '\t' + '\t'.join(INFO[i]) + '\n')
	else: 
		print("Not in info")
			
