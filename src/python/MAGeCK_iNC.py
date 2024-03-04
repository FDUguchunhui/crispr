## Originally developed by Ruilin Tian, Kampmann lab, UCSF; 01/03/2019
# updated by Chunhui Gu, 02/20/2024
## Requirement:
#   MAGeCK--can be downloaded at https://sourceforge.net/p/mageck/wiki/Home/
#   Python packages: pandas, numpy, scipy
# Usage: python MAGeCK+u.py
# The script takes the counts file as input and generates following output files:
#   1. MAGeck outputs: a PDF summary, sgRNA_summary.txt, gene_summary.txt   
#   2. U test outputs: phenotypes-pvalues for all genes (_all_genes.csv), phenotypes-pvalues for FDR hits (_hits.csv), a volcano plot

#  The read count file should contain both groups for comparison and list the names of the sgRNA, the gene it is targeting, followed by the read counts in each sample. 
#  Each item should be separated by the tab ('\t'). A header line is optional. 


import argparse
import random
import shlex
import subprocess
from random import shuffle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def execute(command):
    print(command)
    subprocess.call(shlex.split(command))


def product_threshold_fdr(df, fdr=0.05):
    maxi = abs(df['product']).max()
    for pro in np.arange(0, maxi, 0.01):
        df_thres = df[abs(df['product']) > pro]
        if (1.0 * len(df_thres[df_thres['index'].str.contains('NTC')]) / len(df_thres)) < fdr:
            break
    return pro, df_thres


def rank_test(df):
    df = df[df['treat_mean'] > 20] # only keep sgRNAs with mean read counts > threshold
    df_ntc = df[df['Gene'].str.contains(args.negative_control_keyword)]
    df_targeting = df[~df['Gene'].str.contains(args.negative_control_keyword)]
    ntc_sgRNA_p = list(df_ntc['p.twosided'])
    ntc_sgRNA_p_lfc = zip(list(df_ntc['p.twosided']), list(df_ntc['LFC']))
    genes = df_targeting['Gene'].unique()
    num_of_genes = len(genes)
    gene_lfc_p = {}
    for gene in genes:
        df_gene = df_targeting[df_targeting['Gene'] == gene].sort_values('p.twosided')
        lfc = df_gene.iloc[:5]['LFC'].mean() # the original study has 3 sgRNAs per gene, so we use 5 here
        x, pvalue = mannwhitneyu(list(df_gene['p.twosided'])[:5], ntc_sgRNA_p, alternative='two-sided')
        gene_lfc_p[gene] = [lfc, pvalue]

    random.seed(10)
    for j in range(num_of_genes):
        ntc_sgRNA_p_lfc = list(ntc_sgRNA_p_lfc)
        shuffle(ntc_sgRNA_p_lfc)
        ntc_selected = ntc_sgRNA_p_lfc[:5]
        ntc_selected_p = [i[0] for i in ntc_selected]
        ntc_lfc = np.mean([i[1] for i in sorted(ntc_selected, key=lambda x: x[0])][:5])
        x, ntc_pvalue = mannwhitneyu(ntc_selected_p, ntc_sgRNA_p, alternative='two-sided')
        gene_lfc_p['NTC' + str(j)] = [ntc_lfc, ntc_pvalue]
    return gene_lfc_p



if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Process command line arguments for the script.")

    # Add arguments
    parser.add_argument("--fdr", type=float, required=True, help="False Discovery Rate threshold.")
    parser.add_argument("--counts_file", type=str, required=True, help="Path to the counts file.")
    parser.add_argument("--control_group", type=str, nargs='+', required=True,
                        help="Control group label(s), separated by commas if more than one.")
    parser.add_argument("--counts_thres_control", type=float, default=-1,
                        help="Counts threshold that needed to be greater than for a sgRNA to keep control groups.")
    parser.add_argument("--treatment_group", type=str, nargs='+', required=True,
                        help="Treatment group label(s), separated by commas if more than one.")
    parser.add_argument("--counts_thres_treatment", type=float, default=-1,
                        help="Counts threshold for treatment groups.")
    parser.add_argument("--negative_control_keyword", type=str, required=True, help="sgRNA(s) contain the keyword "
                                                                                    "will be considered as negative control")
    parser.add_argument("--control_sgrna", type=str, required=True, help="Path to the negative control file.")
    parser.add_argument("--output_folder", type=str, required=True, help="Output folder path.")
    parser.add_argument("--output_name", type=str, required=True, help="Comparison name for the output files.")
    parser.add_argument("--gene_list_path", type=str, default='0',
                        help="Path of gene list to label, enter 0 for not labeling genes.")
    parser.add_argument("--label_all_sig_genes", type=int, default=0, choices=[0, 1],
                        help="Label all hit genes? 0 for no, 1 for yes.")
    parser.add_argument("--make_plot", type=int, choices=[0, 1], required=True, help="Generate plot? 0 for no, 1 for yes.")

    args = parser.parse_args()

    # thresholding
    fdr = args.fdr
    counts_file = args.counts_file
    control_groups = args.control_group
    control_thres = args.counts_thres_control
    treatment_groups = args.treatment_group
    treatment_thres = args.counts_thres_treatment
    output_folder = args.output_folder
    output_name = args.output_name
    gene_list_path = args.gene_list_path
    label_all_sig_genes = args.label_all_sig_genes
    make_plot = args.make_plot

    if str(gene_list_path) == '0':
        genes_to_label = []
    else:
        with open(gene_list_path.strip(), 'r') as f:
            genes_to_label = [i.strip() for i in f.readlines()]

    df = pd.read_table(args.counts_file)
    df_thres = df[(df[control_groups] > control_thres).all(axis=1) & (df[treatment_groups] > treatment_thres).all(axis=1)]
    df_thres.to_csv(args.output_folder + '/%s_thresholded_counts.txt' % args.output_name, sep='\t', index=False)

    # MAGeCK
    print("running MAGeCK.....")
    # !mageck test -k data/mageck_input_merged_counts.txt -t sample_het_1 -c sample_wt_1 --control-sgrna data/negative_controls.txt --norm-method control -n report/sample_1/sample_1 --pdf-report --normcounts-to-file
    execute("mageck test -k " + output_folder + '/%s_thresholded_counts.txt' % output_name + " -t " +
            ','.join(treatment_groups) + " -c " + ','.join(control_groups) +
            " -n " + output_folder + "/" + output_name +
            # " --control-sgrna " + args.control_sgrna +
            " --pdf-report --norm-method none")
    # execute("mageck test -k " + output_folder + '/%s_thresholded_counts.txt' % output_name + " -t " +
    #         ','.join(treatment_groups) + " -c " + ','.join(control_groups) +
    #         " -n " + output_folder + "/" + output_name +
    #         " --pdf-report --norm-method none")
    print("running u test.....")
    # u test
    df_mageck = pd.read_table(output_folder + "/" + output_name + '.sgrna_summary.txt')
    df = pd.DataFrame(rank_test(df_mageck)).T
    df.columns = ['epsilon', 'pvalue']
    df.reset_index(inplace=True)
    df['gene'] = df['index'].apply(lambda x: x.split('_')[0])
    df['epsilon'] = -df['epsilon']
    df_ntc = df[df['gene'].str.contains('NTC')]
    df['epsilon'] = df['epsilon'] - df_ntc['epsilon'].median()
    df['product'] = df['epsilon'] * (-np.log10(df['pvalue']))
    df.sort_values('product', ascending=False)

    # FDR
    thres, df_hits = product_threshold_fdr(df, fdr)
    df.sort_values('product', ascending=False).to_csv(output_folder + '/' + output_name + '_all_genes.csv', index=False)
    df_hits.sort_values('product', ascending=False).to_csv(
        output_folder + '/' + output_name + '_fdr%s_product%s_hits.csv' % (fdr, thres), index=False)
    df_ntc = df[df['index'].str.contains('NTC')]

    ########### volcano plot
    if make_plot == 1:

        print("plotting....")
        npg = ["#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2",
               "#7E6148B2"]
        plt.figure(figsize=[10, 8])
        df_pos = df_hits[df_hits['epsilon'] > 0]
        df_neg = df_hits[df_hits['epsilon'] < 0]
        plt.scatter(df_pos['epsilon'], -np.log10(df_pos['pvalue']), c="#DC0000FF", s=5, label='Positive hits')
        plt.scatter(df_neg['epsilon'], -np.log10(df_neg['pvalue']), c="#3C5488FF", s=5, label='Negative hits')
        plt.scatter(df['epsilon'], -np.log10(df['pvalue']), c='#F39B7FFF', s=5, label='Other genes')
        plt.scatter(df_pos['epsilon'], -np.log10(df_pos['pvalue']), c="#DC0000FF", s=5, label=None)
        plt.scatter(df_neg['epsilon'], -np.log10(df_neg['pvalue']), c="#3C5488FF", s=5, label=None)
        plt.scatter(df_ntc['epsilon'], -np.log10(df_ntc['pvalue']), c='grey', s=5, label="Negative control")
        genes = list(df['gene'])
        phenotype = list(df['epsilon'])
        p_value = list(-np.log10(df['pvalue']))
        i = 0
        for x, y, s in zip(phenotype, p_value, genes):
            if s in genes_to_label:
                plt.annotate(s, (x, y), fontsize=16)
                if i == 0:
                    plt.scatter(x, y, c='darkgreen', s=20, label='Genes of interest')
                else:
                    plt.scatter(x, y, c='darkgreen', s=20)
                i = 1
        if str(label_all_sig_genes) == '1':
            genes = list(df_hits['gene'])
            phenotype = list(df_hits['epsilon'])
            p_value = list(-np.log10(df_hits['pvalue']))
            for x, y, s in zip(phenotype, p_value, genes):
                plt.annotate(s, (x, y), fontsize=16)

        x = np.arange(0.01, 10, 0.01)
        y = [thres / i for i in x]
        plt.plot(x, y, '--', c='k', )
        plt.plot(-x, y, '--', c='k', label='FDR = %s' % fdr)
        lim = max(abs(df['epsilon'].min() - 1), (df['epsilon'].max() + 2))
        plt.xlim(-lim, lim)
        plt.ylim(0, -np.log10(df['pvalue'].min()) + 0.5)
        plt.legend(loc=1, fontsize='large', fancybox=True)
        plt.xlabel('Phenotype', fontsize=14)
        plt.ylabel('-log10 P', fontsize=14)
        plt.title(output_name, fontsize=18)
        plt.savefig(output_folder + "/" + output_name + '_volcano_plot.pdf')
        plt.show()
