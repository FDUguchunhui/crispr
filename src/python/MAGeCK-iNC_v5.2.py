## Developed by Ruilin Tian, Kampmann lab, UCSF; 01/03/2019
## Requirement:
#   MAGeCK--can be downloaded at https://sourceforge.net/p/mageck/wiki/Home/
#   Python packages: pandas, numpy, scipy
# Usage: python MAGeCK+u.py
# The script takes the counts file as input and generates following output files:
#   1. MAGeck outputs: a PDF summary, sgRNA_summary.txt, gene_summary.txt   
#   2. U test outputs: phenotypes-pvalues for all genes (_all_genes.csv), phenotypes-pvalues for FDR hits (_hits.csv), a volcano plot

#  The read count file should contain both groups for comparison and list the names of the sgRNA, the gene it is targeting, followed by the read counts in each sample. 
#  Each item should be separated by the tab ('\t'). A header line is optional. 


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import random
from random import shuffle
import subprocess
import shlex
import sys

def execute(command):
    print command
    subprocess.call(shlex.split(command))

def product_threshold_fdr(df,fdr = 0.05):
    maxi = abs(df['product']).max()
    for pro in np.arange(0,maxi,0.1):
        df_thres = df[abs(df['product'])>pro]
        if (1.0 * len(df_thres[df_thres['index'].str.contains('NTC')]) / len(df_thres)) < fdr:
            break
    return pro, df_thres

def rank_test(df):
    df = df[df['treat_mean'] > 20]
    df_ntc = df[df['Gene'].str.contains('non-targeting')]
    df_targeting =  df[~df['Gene'].str.contains('non-targeting')]  
    ntc_sgRNA_p = list(df_ntc['p.twosided'])
    ntc_sgRNA_p_lfc = zip(list(df_ntc['p.twosided']),list(df_ntc['LFC']))
    genes = df_targeting['Gene'].unique()
    num_of_genes = len(genes)
    gene_lfc_p = {}
    for gene in genes:
        df_gene = df_targeting[df_targeting['Gene'] == gene].sort_values('p.twosided')
        lfc = df_gene.iloc[:3]['LFC'].mean()
        x, pvalue = mannwhitneyu(list(df_gene['p.twosided'])[:3],ntc_sgRNA_p,alternative='two-sided')
        gene_lfc_p[gene] = [lfc,pvalue]

    random.seed(10)
    for j in range(num_of_genes):
        shuffle(ntc_sgRNA_p_lfc)
        ntc_selected = ntc_sgRNA_p_lfc[:5]
        ntc_selected_p = [i[0] for i in ntc_selected]
        ntc_lfc = np.mean([i[1] for i in sorted(ntc_selected, key=lambda x: x[0])][:3])
        x, ntc_pvalue = mannwhitneyu(ntc_selected_p,ntc_sgRNA_p,alternative='two-sided')
        gene_lfc_p['NTC_' + str(j)] = [ntc_lfc, ntc_pvalue]
    return gene_lfc_p

print ('fdr: ')
fdr = float(raw_input('-->'))
print ('counts file: ')
counts_file = raw_input('-->')

print ('control_group label (if more than 1, seperate by ,): ')
control_group = raw_input('-->')

print ('counts threshold_control_groups: ')
counts_thres_control = float(raw_input('-->'))

print ('treatment_group label (if more than 1, seperate by ,): ')
treatment_group = raw_input('-->')

print ('counts threshold_treatment_groups: ')
counts_thres_treatment = float(raw_input('-->'))

print ('output folder: ')
output_folder = raw_input('-->')
print ('comparison name: ')
output_name = raw_input('-->')
print ('path of gene list to label: (enter 0 for not labeling genes) ')
gene_list_path = raw_input('-->')
print ('label all hit genes? (0 or 1) ')
label_all_sig_genes = raw_input('-->')
print ('plot( 0 or 1 )?: ')
make_plot = eval(raw_input('-->'))
if str(gene_list_path) == '0':
    genes_to_label = []
else:
    with open(gene_list_path.strip(),'r') as f:
        genes_to_label = [i.strip() for i in f.readlines()]

# thresholding
control_groups = str(control_group).split(',')
control_thres = counts_thres_control
treatment_groups = str(treatment_group).split(',')
treatment_thres = counts_thres_treatment
df = pd.read_table(counts_file)
df_thres = df[(df[control_groups] > control_thres).all(axis=1)&(df[treatment_groups] > treatment_thres).all(axis=1)]
df_thres.to_csv(output_folder+'/%s_thresholded_counts.txt'%output_name,sep='\t',index=False)

# MAGeCK
print "running MAGeCK....."
execute("mageck test -k " + output_folder+'/%s_thresholded_counts.txt'%output_name + " -t "+ treatment_group + " -c " + control_group + " -n " + output_folder + "/" + output_name + " --pdf-report")
print "running u test....."
# u test
df_mageck = pd.read_table(output_folder + "/" + output_name + '.sgrna_summary.txt')
df = pd.DataFrame(rank_test(df_mageck)).T
df.columns = ['epsilon','pvalue']
df.reset_index(inplace=True)
df['gene'] = df['index'].apply(lambda x:x.split('_')[0])
df['epsilon'] = -df['epsilon']
df_ntc = df[df['gene']=='NTC']
df['epsilon'] = df['epsilon'] - df_ntc['epsilon'].median()
df['product'] = df['epsilon'] * (-np.log10(df['pvalue']))
df.sort_values('product',ascending=False)


# FDR
thres, df_hits = product_threshold_fdr(df,fdr)
df.sort_values('product',ascending=False).to_csv(output_folder + '/' + output_name + '_all_genes.csv',index=False)
df_hits.sort_values('product',ascending=False).to_csv(output_folder + '/' + output_name + '_fdr%s_product%s_hits.csv'%(fdr,thres),index=False)
df_ntc = df[df['index'].str.contains('NTC')]

########### volcano plot
if make_plot == 1:

    print "plotting...."
    npg = ["#E64B35B2","#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2"]
    plt.figure(figsize=[10,8])
    df_pos = df_hits[df_hits['epsilon']>0]
    df_neg = df_hits[df_hits['epsilon']<0]
    plt.scatter(df_pos['epsilon'],-np.log10(df_pos['pvalue']),c="#DC0000FF",s=5,label='Positive hits')
    plt.scatter(df_neg['epsilon'],-np.log10(df_neg['pvalue']),c="#3C5488FF",s=5,label='Negative hits')
    plt.scatter(df['epsilon'],-np.log10(df['pvalue']),c='#F39B7FFF',s=5,label='Other genes')
    plt.scatter(df_pos['epsilon'],-np.log10(df_pos['pvalue']),c="#DC0000FF",s=5,label=None)
    plt.scatter(df_neg['epsilon'],-np.log10(df_neg['pvalue']),c="#3C5488FF",s=5,label=None)
    plt.scatter(df_ntc['epsilon'],-np.log10(df_ntc['pvalue']),c='grey',s=5,label="Negative control")
    genes=list(df['gene'])
    phenotype = list(df['epsilon'])
    p_value = list(-np.log10(df['pvalue']))
    i = 0
    for x, y, s in zip(phenotype, p_value,genes):
        if s in genes_to_label:
            plt.annotate(s,(x,y),fontsize=16)
            if i == 0:
                plt.scatter(x,y,c='darkgreen',s=20,label = 'Genes of interest')
            else:
                 plt.scatter(x,y,c='darkgreen',s=20)
            i = 1
    if str(label_all_sig_genes) == '1':
        genes=list(df_hits['gene'])
        phenotype = list(df_hits['epsilon'])
        p_value = list(-np.log10(df_hits['pvalue']))
        for x, y, s in zip(phenotype, p_value,genes):
            plt.annotate(s,(x,y),fontsize=16)
        
    x = np.arange(0.01,10,0.01)
    y = [thres / i for i in x]
    plt.plot(x,y,'--', c= 'k',)
    plt.plot(-x,y,'--',c = 'k',label='FDR = %s'%fdr)
    lim = max(abs(df['epsilon'].min() - 1),(df['epsilon'].max() + 2))
    plt.xlim(-lim,lim)
    plt.ylim(0,-np.log10(df['pvalue'].min()) + 0.5)
    plt.legend(loc = 1,fontsize='large',fancybox=True)
    plt.xlabel('Phenotype',fontsize = 14)
    plt.ylabel('-log10 P',fontsize = 14)
    plt.title(output_name,fontsize=18)
    plt.savefig(output_folder + "/" + output_name + '_volcano_plot.pdf')
    plt.show()

