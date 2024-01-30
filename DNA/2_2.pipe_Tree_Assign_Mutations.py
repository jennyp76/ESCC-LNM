#!/usr/bin/python3

""" Sort mutations into trunk, shared and private for each patient
Trunk: Detected in all regions
Shared: Detected in not all but some regions
Private: Detected in only one region
Pre-metastasis (Trunk), Post-metastasis (Shared, Private)"""


# Import Modules
import glob
import sys
import os

# Input Parameters
Force_calling_Path = sys.argv[1]
PATIENT = sys.argv[2]

dic = dict()
trunk = []
shared = []
private = {}

def collect_mutation(file):
	data = open( file,"r")
	name = os.path.basename(file).split('.')[0]
	for line in data:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		else:
			col = line.split('\t')
			chr = col[0]
			pos = col[1]
			ref = col[2]
			alt = col[3]

			key = chr+'\t'+pos+'\t'+ref+'\t'+alt 
	
			if key in dic.keys():
				pass
			else:
				dic[key]=[]
			dic[key].append(name)

files = glob.glob("{}/{}_*.rescue.txt".format(Force_calling_Path, PATIENT)) 
n = len(files)

for file in files:
	collect_mutation(file)	

for i in dic.keys():
	if len(dic[i]) == n :
		trunk.append(i)
	elif len(dic[i]) == 1:
		if dic[i][0] in private.keys():
			pass
		else:
			private[dic[i][0]]=[]
			
		private[dic[i][0]].append(i)
	else:
		shared.append(i)

for i in private.keys():
	Output = open("{}/{}.RGadded.marked.fixed.merge.rescue.private.txt".format(Force_calling_Path,PATIENT),'w')
	Output.write("#chr\tpos\tref\talt\n")
	for j in private[i]:
		Output.write( j + '\n' )

Output = open("{}/{}_trunk.RGadded.marked.fixed.merge.rescue.trunk.txt".format(Force_calling_Path,PATIENT),'w')
Output.write('#chr\tpos\tref\talt\n')
for i in trunk:
	Output.write( i + '\n')

if len(shared) == 0 : 
	print ('No shared variant')
else:
	Output = open("{}/{}_shared.RGadded.marked.fixed.merge.rescue.shared.txt".format(Force_calling_Path,PATIENT),'w')
	Output.write("#chr\tpos\tref\talt\n")
	for i in shared :
		Output.write( i + '\n')
