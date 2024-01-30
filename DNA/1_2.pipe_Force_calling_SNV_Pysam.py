#!/usr/bin/python3

""" Pysam: Count read for all SNV detected by Mutect2 within each patient """

# Import Modules
import pysam 
import glob
import sys
import os

# Input Parameters
PATIENT = sys.argv[1]
BAM_PATH =sys.argv[2]
VCF_PATH = sys.argv[3]
OUTPUT_PATH = sys.argv[4]

vcf_files=glob.glob(("{}/{}_*.RGadded.marked.fixed.mutect2.FILTER.PASS.snv.vcf").format(VCF_PATH,PATIENT))
bam_files=glob.glob(("{}/{}_*.RGadded.marked.fixed.MAPQ_30.bam").format(BAM_PATH,PATIENT))

def read_VCF(file, vcf_set ):
	data = open(file, 'r')
	NAME = os.path.basename(file).split('.')[0]
	print(NAME)
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
			key = CHR + '\t' + POS + '\t' + REF + '\t' + ALT
			if key not in vcf_set:
				vcf_set.add(key)
			else:
				continue
	return vcf_set

def read_BAM(bam_file, final_vcf_set, dict, sample ):
	bam = pysam.AlignmentFile('%s'%bam_file,'rb')
	for key in final_vcf_set :
		CHR=key.split('\t')[0]
		POS=key.split('\t')[1]
		REF=key.split('\t')[2]
		ALT=key.split('\t')[3]
		counts = bam.count_coverage(CHR,int(POS)-1, int(POS), quality_threshold = 20)
		count = {}
		count['A']=int(counts[0][0])
		count['C']=int(counts[1][0])
		count['G']=int(counts[2][0])
		count['T']=int(counts[3][0])
		DP =  count['A'] + count['C'] + count['G'] + count['T']

		if int(DP) >= 10:
			for alt,cnt in count.items():
				if alt == ALT and cnt != 0:
					if key in dict.keys():
						pass
					else:
						dict[key] = {}
					dict[key][sample]=(cnt,DP)
				else:
					continue
		else:
			continue

	return dict

def iterate_VCF(FILES, vcf_set):
	for file in FILES:
		vcf_set = read_VCF(file, vcf_set)
	return vcf_set

def iterate_BAM(FILES, final_vcf_set, dict):
	name = []
	final_dict={}
	for file in FILES:
		sample = os.path.basename(file).split('.')[0]
		name.append(sample)
		print(('Running....{}').format(file))
		dict = read_BAM( file, final_vcf_set, dict, sample)
		print(('DONE.....{}').format(file))
	final_dict = dict
	return final_dict, name

def sort_chr(chr):
        if chr[3:] == 'X':
                return 23
        elif chr[3:] == 'Y':
                return 24
        elif chr[3:] == 'M':
                return 0
        else :
                return int(chr[3:])
def sort_Chr(dict) :
        sorted_key = sorted(dict.keys(), key=lambda x : int(x.split('\t')[1]))
        sorted_key = sorted(sorted_key, key = lambda x : sort_chr(x.split('\t')[0]))
        return sorted_key

def write_output(final_dict, name):
	output=open("{}/{}.RGadded.marked.fixed.mutect2.FILTER.PASS.snv.pysam.txt".format(OUTPUT_PATH,PATIENT),"w")
	output.write("#CHR\tPOS\tREF\tALT\t"+'\t'.join(name)+'\n')
	print('Writing Header.......')
	sorted_key = sort_Chr(final_dict)
	for k_key in sorted_key :
		temp =[]
		for n_name in name :
			if n_name in final_dict[k_key].keys():
				nALT = final_dict[k_key][n_name][0]
				DP = final_dict[k_key][n_name][1]
				temp.append(str(nALT)+'/'+str(DP))
			else:
				temp.append('-')

		if len(set(temp)) == 1 :   # exclude when all = '-'
			print ('None counted in all samples:',k_key)
		else:
			output.write(k_key+'\t'+'\t'.join(temp)+'\n')
	
BAM_FILES = []
dict={}
vcf_set = set()

for i in range(0,2):
	BAM_FILES.append('')
for file in bam_files:
	if 'Normal' in file:
		BAM_FILES[0] = file		
	elif 'Primary' in file:
		BAM_FILES[1] = file
	else:
		BAM_FILES.append(file)

print(BAM_FILES)
print('STARTING..')
final_vcf_set = iterate_VCF(vcf_files,vcf_set)
print('MUTATION COLLECTING..')
final_dict, name  = iterate_BAM(BAM_FILES, final_vcf_set, dict)
print('WRITING OUTPUT...')
write_output(final_dict, name)	

