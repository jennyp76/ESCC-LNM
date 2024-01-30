#!/usr/bin/python3

""" Rescue SNVs using multi-regional samples within each patients"""

# Import Modules
import os
import sys
import glob

# Input Parameters
PATIENT = sys.argv[1]
VCF_PATH = sys.argv[2]
OUTPUT_PATH = sys.argv[3]

file = open('{}/{}.RGadded.marked.fixed.mutect2.FILTER.PASS.snv.pysam.txt'.format(OUTPUT_PATH,PATIENT))
vcf_files=glob.glob(("{}/{}_*.RGadded.marked.fixed.mutect2.FILTER.PASS.snv.vcf").format(VCF_PATH,PATIENT))

def read_VCF(file, vcf_dict):
        data = open(file, 'r')
        NAME = os.path.basename(file).split('.')[0]
        for line in data :
                if line.startswith('##'):
                        continue
                elif line.startswith('#C'):
                        col = line.split('\t')
                        if 'Normal' in col[9] :
                                NORMAL = 9
                                TUMOR = 10
                        else:
                                NORMAL = 10
                                TUMOR = 9
                else:
                        line = line.rstrip()
                        col = line.split('\t')
                        CHR = col[0]
                        POS = col[1]
                        REF = col[3]
                        ALT = col[4]
                        key = CHR+'\t'+POS+'\t'+REF+'\t'+ALT
                        DP = col[TUMOR].split(':')[3]
                        nALT = col[TUMOR].split(':')[1].split(',')[1]

                        if key not in vcf_dict.keys():
                                vcf_dict[key]={}
                        else:
                                pass

                        vcf_dict[key][NAME]= str(nALT)+'/'+str(DP)
        return vcf_dict

def iterate_VCF(FILES, vcf_dict):
        for file in FILES:
                Mutect2_dict = read_VCF(file, vcf_dict)
        return Mutect2_dict

output_dict = {}
name = {}
vcf_dict = {}

Mutect2_dict = iterate_VCF(vcf_files, vcf_dict)

for line in file :
	line = line.rstrip()
	col = line.split('\t')
	if line.startswith('#'):
		for i in range(5,len(col)):
			name[i]=col[i]  	
		for i in range(5,len(col)):
			output_dict[col[i]]=[]  	
	else:
		k_key = col[0]+'\t'+col[1]+'\t'+ col[2]+'\t'+col[3]
		
		if col[4] == '-' or int(col[4].split('/')[0]) < 3 : 
			if k_key in Mutect2_dict.keys():
				for i in range(5,len(col)):
					if col[i] == '-':
						continue
					else:
						if name[i] in Mutect2_dict[k_key].keys() :
							output_dict[name[i]].append( [k_key,Mutect2_dict[k_key][name[i]]])
						else:
							nalt = int(col[i].split('/')[0])
							DP = int(col[i].split('/')[1])
							if nalt >= 3: 
								output_dict[name[i]].append([k_key, col[i]])
							else:
								continue
			else:
				print( 'Error with key' , k_key )	
		else:
			if k_key in Mutect2_dict.keys():
				for k in Mutect2_dict[k_key].keys():
					output_dict[k].append( [k_key,Mutect2_dict[k_key][k]])
			else:
				print( 'Error with key' , k_key )

for i in range(0,len(name)):
	output_file = open(('{}/{}.RGadded.marked.fixed.snv.rescue.txt').format(OUTPUT_PATH,name[i+5]),'w')
	output_file.write('#CHR\tPOS\tREF\tALT\tCount\n')
	for k in output_dict[name[i+5]] :
		output_file.write('\t'.join(k)+'\n')
