# SNPEff parser
import os
import sys
import gzip
import binascii
import re
import numpy as np
import pandas as pd

def get_csv_files(directory):
    """Return a list of all CSV files in the specified directory."""
    try:
    	csv_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv')]
    	return csv_files 
    except Exception as e:
    	print(f"Problem finding descriptor files in {directory}")

def extract_name_from_comment(filepath):
    """Extract the name from the first comment line in the CSV file."""
    with open(filepath, 'r') as file:
        line = file.readline().strip() 
        if line.startswith('#'):  # Look for the comment line
                return line[1:].split(' ')[0].strip()  # Remove the # and any surrounding whitespace from the first word after the comment
    raise ValueError(f"First line invalid: {filepath}")

def readAnnotationMatrices():
	data = {}
	data_files = get_csv_files('data')
		
	for filename in data_files:
			print(f"Reading data from {filename}...")
			try:
				descriptor = extract_name_from_comment(filename)
				data[descriptor] = pd.read_csv(filename, index_col = 0, comment = '#', skip_blank_lines = True)
			except Exception as e:
				print(e)
	return(data)
	
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

def load_aa_data():
	keys = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
	values = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	aas = dict(zip(keys, values))
	return(aas)
	
def aa3_to_aa1(aaa, aa_dict):
	if aaa != '':
		code = aa_dict[aaa]
	else:
		code = ''
	return(code)
	
def parseSNPEffFile(f, aa_dict):
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
				    # Handle only variants annotated as not 'MODIFIER'
				    if not 'MODIFIER' in item:
				        try:
				            parsed_vals = parseSNPEffAnn(item, aa_dict)
				            parsed_str = ' '.join(parsed_vals.values()).strip()
				            print(scaffold, coord, ref, var, parsed_str)
				        except Exception as e:
				        	print(f"Error parsing SNPEff annotation: {e}")


def parseSNPEffAnn(ann, aa_dict):
#['missense_variant|MODERATE|Sc9M7eS_1763_HRSCAF_2674_28679|gene07994|transcript|mRNA07994|protein_coding|4/14|c.3640G>A|p.Glu1214Lys|3640/5781|3640/5781|1214/1926||,T|', 'missense_variant|MODERATE|Sc9M7eS_1763_HRSCAF_2674_28679|gene07996|transcript|mRNA07996|protein_coding|4/13|c.3640G>A|p.Glu1214Lys|3640/5610|3640/5610|1214/1869||,T|', 'intron_variant|MODIFIER|Sc9M7eS_1763_HRSCAF_2674_28679|gene07995|transcript|mRNA07995|protein_coding|3/9|c.455-4448G>A||||||,T|', 'intron_variant|MODIFIER|Sc9M7eS_1763_HRSCAF_2674_28679|gene07997|transcript|mRNA07997|protein_coding|3/10|c.455-4448G>A||||||,T|', 'intron_variant|MODIFIER|Sc9M7eS_1763_HRSCAF_2674_28679|gene07998|transcript|mRNA07998|protein_coding|3/8|c.455-4448G>A||||||;AN=26;AC=7']
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
	d['aa1'] = aa3_to_aa1(d['ref'], aa_dict)
	d['aa2'] = aa3_to_aa1(d['var'], aa_dict)
	if (d['ref']) != '':
		for key in ann_data:
			d[key + 'RV'] = str(round((max_vals[key] - ann_data[key][d['aa1']][d['aa2']]) / max_vals[key], 2))
			d[key + 'VR'] = str(round((max_vals[key] - ann_data[key][d['aa2']][d['aa1']]) / max_vals[key], 2))
	return(d)

def create_key_suffix_string(my_dict):
    """Create a string with each key appearing twice, with suffixes 'RV' and 'VR'."""
    result_list = []
    
    for key in my_dict:
        result_list.append(f"{key}_RefVar")  # Append the key with 'RV' suffix
        result_list.append(f"{key}_VarRef")  # Append the key with 'VR' suffix
    
    # Join the list elements into a single string, separated by spaces (or any delimiter you prefer)
    result_string = ' '.join(result_list)
    
    return result_string
    
def check_symmetry(df):
    """Check if a DataFrame is symmetric, allowing for floating-point precision and handling NaN values."""
    # First, check if the DataFrame is square
    if df.shape[0] != df.shape[1]:
        return False
    # Replace NaN values with 0 for comparison
    df_filled = np.nan_to_num(df)  # Fill NaN values with 0 or any suitable value
    # Use numpy's allclose for comparison to handle floating-point precision
    return np.allclose(df_filled, df_filled.T)
    
##################### Main ##################
aa_dict = load_aa_data()
ann_data = readAnnotationMatrices()
#print(ann_data)
max_vals = {}
for key in ann_data:
	max_vals[key] = np.nanmax(ann_data[key].values)
	
is_symmetric = {}
for key in ann_data:
	is_symmetric[key] = check_symmetry(ann_data[key].values)
#print(max_vals)

header1 = 'Scaffold Coord Ref Var Type Effect Transcript Ref_aa Coord_aa Var_aa Ref_aa_abbrev Var_aa_abbrev'
header2 = create_key_suffix_string(ann_data)
print(header1 + ' ' + header2)
if is_gz_file(sys.argv[1]):
    with gzip.open(sys.argv[1], 'rt') as f:
    	parseSNPEffFile(f, aa_dict)
else:
    with open(sys.argv[1], 'r') as f:
    	parseSNPEffFile(f, aa_dict)

print(is_symmetric)

