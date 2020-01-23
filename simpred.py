# SNPEff parser

import sys
import gzip
import binascii
import re
import numpy as np
import pandas as pd

def readAnnotationMatrices():
	data = {}
	data['aai7'] = pd.read_csv('data/idx_7aaI.csv', index_col = ['Src'], skiprows = 4)
	data['exchgb'] = pd.read_csv('data/idx_exchangeability.csv', index_col = ['Src'], skiprows = 4)
	data['sneath'] = pd.read_csv('data/idx_sneath_dissimilarity.csv', index_col = ['name'], skiprows = 4)
	return(data)
	
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

def load_aa_data():
	keys = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
	values = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	aas = dict(zip(keys, values))
	return(aas)
	
def aa3_to_aa1(aaa):
	global aas
	if aaa != '':
		code = aas[aaa]
	else:
		code = ''
	return(code)
	
def parseSNPEffFile(f):
	for line in f:
		if not line.startswith('#'):						# handle comments
			if 'missense_variant' in line:					# fish out missense variants only
				tmp = line.strip().split('\t')
				scaffold = tmp[0]
				coord = tmp[1]
				ref = tmp[3]
				var = tmp[4]
				annotation = tmp[7]
				#print(annotation)
				pattern = r'([^\|]+_variant\|)'
				new_string = re.sub(pattern, '[cuthere]\\1', annotation)
				#print(new_string)
				ann = new_string.split('[cuthere]')
				ann = ann[1:]
				#print(ann)
				for item in ann:
					vals = parseSNPEffAnn(item).values()
					parsed = ''
					for i in vals:
						parsed = parsed + ' ' + i
					if "MODIFIER" not in parsed:
						print(scaffold, coord, ref, var, parsed.strip())
        
def parseSNPEffAnn(ann):
#['missense_variant|MODERATE|Sc9M7eS_1763_HRSCAF_2674_28679|gene07994|transcript|mRNA07994|protein_coding|4/14|c.3640G>A|p.Glu1214Lys|3640/5781|3640/5781|1214/1926||,T|', 'missense_variant|MODERATE|Sc9M7eS_1763_HRSCAF_2674_28679|gene07996|transcript|mRNA07996|protein_coding|4/13|c.3640G>A|p.Glu1214Lys|3640/5610|3640/5610|1214/1869||,T|', 'intron_variant|MODIFIER|Sc9M7eS_1763_HRSCAF_2674_28679|gene07995|transcript|mRNA07995|protein_coding|3/9|c.455-4448G>A||||||,T|', 'intron_variant|MODIFIER|Sc9M7eS_1763_HRSCAF_2674_28679|gene07997|transcript|mRNA07997|protein_coding|3/10|c.455-4448G>A||||||,T|', 'intron_variant|MODIFIER|Sc9M7eS_1763_HRSCAF_2674_28679|gene07998|transcript|mRNA07998|protein_coding|3/8|c.455-4448G>A||||||;AN=26;AC=7']
	global ann_exchgb
	d = {}
	#print('Parsing ' + ann)
	tmp = ann.strip().split('|')
	d['type'] = tmp[0]
	d['eff'] = tmp[1]
	#d['scaff'] = tmp[1]
	#d['gene'] = tmp[2]
	d['tr'] = tmp[5]
	#d['bases'] = tmp[7].replace('c.', '')
	aa = tmp[9].replace('p.', '')
	d['ref'] = aa[:3]
	d['aa'] = aa[3:-3]
	d['var'] = aa[-3:]
	d['aa1'] = aa3_to_aa1(d['ref'])
	d['aa2'] = aa3_to_aa1(d['var'])
	if (d['ref']) != '':
		d['aai7'] = str(round((ann_aai7[d['aa1']][d['aa2']]) / ann_aai7_max, 2))
		d['exchgb1'] = str(round((ann_exchgb_max - ann_exchgb[d['aa1']][d['aa2']]) / ann_exchgb_max, 2))
		d['exchgb2'] = str(round((ann_exchgb_max - ann_exchgb[d['aa2']][d['aa1']]) / ann_exchgb_max, 2))
		d['sneath'] = str(round((ann_sneath[d['aa2']][d['aa1']]) / ann_sneath_max, 2))
	#print(d)
	return(d)
		
##################### Main ##################
ann_data = readAnnotationMatrices()
ann_aai7 = ann_data['aai7']
ann_exchgb = ann_data['exchgb']
ann_sneath = ann_data['sneath']
aas = load_aa_data()
ann_aai7_max = np.nanmax(ann_data['aai7'].values)
ann_exchgb_max = np.nanmax(ann_data['exchgb'].values)
ann_sneath_max = np.nanmax(ann_data['sneath'].values)

print('Scaffold Coord Ref Var Type Effect Transcript Ref_aa Coord_aa Var_aa Ref_aa_abbrev Var_aa_abbrev aaI7 exchgb_ref_var exchgb_var_ref sneath_dissim')
if is_gz_file(sys.argv[1]):
    with gzip.open(sys.argv[1], 'rt') as f:
    	parseSNPEffFile(f)
else:
    with open(sys.argv[1], 'r') as f:
    	parseSNPEffFile(f)



