# SNPEff Parser
import sys
import gzip
import binascii
import re
import numpy as np
import pandas as pd

def read_annotation_matrices():
    """Reads and loads annotation matrices from CSV files into a dictionary."""
    try:
        data = {
            'aai7': pd.read_csv('data/idx_7aaI.csv', index_col='Src', skiprows=4),
            'exchgb': pd.read_csv('data/idx_exchangeability.csv', index_col='Src', skiprows=4),
            'sneath': pd.read_csv('data/idx_sneath_dissimilarity.csv', index_col='name', skiprows=4)
        }
        return data
    except Exception as e:
        print(f"Error reading annotation matrices: {e}")
        sys.exit(1)

def is_gz_file(filepath):
    """Checks if the provided file is a gzipped file."""
    try:
        with open(filepath, 'rb') as f:
            return binascii.hexlify(f.read(2)) == b'1f8b'
    except Exception as e:
        print(f"Error checking if file is gzipped: {e}")
        sys.exit(1)

def load_amino_acid_data():
    """Loads amino acid mapping from 3-letter to 1-letter abbreviations."""
    keys = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
    values = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    return dict(zip(keys, values))

def aa3_to_aa1(aa3, aa_map):
    """Converts a 3-letter amino acid code to a 1-letter code using the provided map."""
    try:
        return aa_map[aa3]
    except KeyError:
        return ''  # Return empty string if 3-letter code is not found

def parse_snpeff_file(file, aa_map, ann_data, ann_aai7_max, ann_exchgb_max, ann_sneath_max):
    """Parses the SNPEff file and prints out relevant information for missense variants."""
    for line in file:
        if not line.startswith('#') and 'missense_variant' in line:  # Process only missense variants
            fields = line.strip().split('\t')
            scaffold, coord, ref, var, annotation = fields[0], fields[1], fields[3], fields[4], fields[7]
            
            annotations = re.split(r'[^\|]+_variant\|', annotation)[1:]  # Split annotations by variant type
            for item in annotations:
                try:
                    parsed_vals = parse_snpeff_annotation(item, aa_map, ann_data, ann_aai7_max, ann_exchgb_max, ann_sneath_max)
                    parsed_str = ' '.join(parsed_vals.values()).strip()
                    if "MODIFIER" not in parsed_str:
                        print(scaffold, coord, ref, var, parsed_str)
                except Exception as e:
                    print(f"Error parsing SNPEff annotation: {e}")

def parse_snpeff_annotation(ann, aa_map, ann_data, ann_aai7_max, ann_exchgb_max, ann_sneath_max):
    """Parses a single SNPEff annotation line into a structured format."""
    ann_aai7, ann_exchgb, ann_sneath = ann_data['aai7'], ann_data['exchgb'], ann_data['sneath']
    parsed = {}
    
    fields = ann.strip().split('|')
    parsed['type'] = fields[0]
    parsed['eff'] = fields[1]
    parsed['tr'] = fields[5]
    
    aa_change = fields[9].replace('p.', '')
    parsed['ref'] = aa_change[:3]
    parsed['aa'] = aa_change[3:-3]
    parsed['var'] = aa_change[-3:]
    
    parsed['aa1'] = aa3_to_aa1(parsed['ref'], aa_map)
    parsed['aa2'] = aa3_to_aa1(parsed['var'], aa_map)
    
    if parsed['ref']:
        try:
            parsed['aai7'] = str(round(ann_aai7[parsed['aa1']][parsed['aa2']] / ann_aai7_max, 2))
            parsed['exchgb1'] = str(round((ann_exchgb_max - ann_exchgb[parsed['aa1']][parsed['aa2']]) / ann_exchgb_max, 2))
            parsed['exchgb2'] = str(round((ann_exchgb_max - ann_exchgb[parsed['aa2']][parsed['aa1']]) / ann_exchgb_max, 2))
            parsed['sneath'] = str(round(ann_sneath[parsed['aa2']][parsed['aa1']] / ann_sneath_max, 2))
        except KeyError as e:
            print(f"Error in accessing matrix data: {e}")
        except ZeroDivisionError:
            print("Division by zero occurred in matrix calculations.")
            parsed['aai7'], parsed['exchgb1'], parsed['exchgb2'], parsed['sneath'] = 'NaN', 'NaN', 'NaN', 'NaN'
    
    return parsed

##################### Main ##################
if __name__ == '__main__':
    # Load annotation matrices and amino acid data
    try:
        ann_data = read_annotation_matrices()
        aa_map = load_amino_acid_data()
        
        ann_aai7_max = np.nanmax(ann_data['aai7'].values)
        ann_exchgb_max = np.nanmax(ann_data['exchgb'].values)
        ann_sneath_max = np.nanmax(ann_data['sneath'].values)
        
        print('Scaffold Coord Ref Var Type Effect Transcript Ref_aa Coord_aa Var_aa Ref_aa_abbrev Var_aa_abbrev aaI7 exchgb_ref_var exchgb_var_ref sneath_dissim')

        # Open the input file (gzipped or plain text) and parse the SNPEff data
        input_file = sys.argv[1]
        if is_gz_file(input_file):
            with gzip.open(input_file, 'rt') as f:
                parse_snpeff_file(f, aa_map, ann_data, ann_aai7_max, ann_exchgb_max, ann_sneath_max)
        else:
            with open(input_file, 'r') as f:
                parse_snpeff_file(f, aa_map, ann_data, ann_aai7_max, ann_exchgb_max, ann_sneath_max)
    except Exception as e:
        print(f"An error occurred during file processing: {e}")
        sys.exit(1)
