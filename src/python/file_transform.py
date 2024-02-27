import argparse
from functools import reduce

import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description='Process count files and sample names to merge them.')

# Adding arguments
parser.add_argument('-f', '--count-file', required=True, help='Count files separated by spaces')
parser.add_argument('-s', '--samples', nargs='+', required=True, help='Sample names separated by spaces')
parser.add_argument('-o', '--output', required=True, help='Output folder path')
parser.add_argument('-n', '--name', required=True, help='Name for the output file')

# Parse arguments
args = parser.parse_args()


dfs = pd.read_table(args.count_file, usecols=['sgRNA'] + args.samples, sep=',', low_memory=False)


def extract_gene(s):
    # Check if 'negative control' is in the string
    if 'negative_control' in s:
        return s  # Return the original string
    else:
        # Split the string by underscore and return the first part
        return s.split('_')[0]

df_merge = dfs
df_merge['Gene'] = df_merge['sgRNA'].apply(extract_gene)
cols = df_merge.columns.tolist()
cols = [cols[0]] + [cols[-1]] + cols[1:-1]
df_merge = df_merge[cols]
df_merge.to_csv(args.output + '/%s_merged_counts.txt' % args.name, sep='\t', index=False)