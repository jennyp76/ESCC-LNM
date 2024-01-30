#!/usr/bin/python3

""" After force calling (manually check through IGV), divide indels (including rescued) into individual sample file """

# Import Modules
import glob
import os
import sys

# Input Parameters
PATIENT = sys.argv[1]
Force_calling_path = sys.argv[2]

Result = {}
Name = {}

data =  open('{}/{}.somatic.indels.PASS.total.rescue.txt'.format(Force_calling_path,PATIENT),'r')
for line in data :
	line = line.rstrip('\r\n')
	col = line.split('\t')
	if line.startswith('#chr'):
		samples = col[4:-2]
		for n in range(len(samples)):
			Result[samples[n]]= []  
			Name[str(n)]=samples[n]
	else: 
		stat = col[-1]
		if stat == '-':
			chr = col[0]
			pos = col[1]
			ref = col[2]
			alt = col[3]
			key = chr + '\t' + pos + '\t' + ref + '\t' + alt
			result = col[4:-2]
			for i, value in enumerate(result) :
				if value  == '-' :
					continue
				else: 
					Result[Name[str(i)]].append(key)
		else:
			continue   #exclue status == "insertion/deletion, soft-masked, normal"
			
for sample in Result.keys():
	output = open('{}/{}.somatic.indels.PASS.rescue.txt'.format(Force_calling_path, sample),'w')
	output.write('#CHR\tPOS\tREF\tALT\n')
	for i in Result[sample]:
		output.write( i + '\n')

			
	
