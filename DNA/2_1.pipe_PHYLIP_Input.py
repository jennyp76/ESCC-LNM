#!/usr/bin/python3

""" Binary input for PHYLIP that includes SNVs and indels  """

# Import Modules
import os
import sys
import glob
import natsort


# Input Parameters
PATIENT = sys.argv[1]
Force_calling_Path =  sys.argv[3]
Phylip_Path = sys.argv[4]

files= natsort.natsorted(glob.glob('{}/{}_*.rescue.txt'.format(Force_calling_Path,PATIENT)),reverse = True)

variant_dict = {}
variant_set = set()
NAME = []
LIST=[]  

NAME.append('Normal    ')

for file in files:
	name = os.path.basename(file) 
	name = name.split('.')[0]
	if 'Primary'  in name:
		name = name.split('_')[-1]
	else:
		name = name.split('_')[-1]
		if len(name) > 7 :
			name = name[0:5]
		else:
			pass
		name = 'N'+ name 
	if len(name) != 10 :
		name = name + ' '*(10-len(name))

	NAME.append(name)
	variant_dict[name]=[]
	
	Data = open (file,'r')
	for line in Data :
		line = line.rstrip()
		if line.startswith('#CHR') :
			continue 
		else:
			col = line.split('\t')
			chr = col[0]
			pos = col[1]
			ref = col[2]
			alt = col[3]
			key = chr+'.'+pos+ref+'>'+alt 

			variant_set.add(key)
			variant_dict[name].append(key)

for i in range(0,len(NAME)) :
       temp = []
       for j in range(0, len(variant_set)) :
	       temp.append(str(0))
       LIST.append(temp)

for j, var in enumerate(variant_set) :
	for i , n in enumerate(NAME):
		if n == 'Normal    ':
			continue
		else:
			if var in variant_dict[n]:
				LIST[i][j] = str(1)
			else:
				LIST[i][j] =str(0) 

Output = open('{}/{}.phylip.txt'.format(Phylip_Path,PATIENT),"w")
Output.write('     '+str(len(NAME))+'    '+str(len(variant_set))+'\n')

for idx, i in enumerate(NAME):
	Output.write(i+''.join(LIST[idx])+'\n')
