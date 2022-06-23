############################################################################################
# Objective:
#   Perform non-parametric hypothesis testing between abundance values in multiple groups\
#       (Control vs GWI) and (Chow vs WD)

# Author: Vitalii Stebliankion (vsteb002@fiu.edu)
#   Bioinformatics Research Group (BioRG)
############################################################################################


import pandas as pd
from scipy.stats import mannwhitneyu
import PyPluMA

#ABR_df_abundance = pd.read_csv("./ABR_counts/ABR_profile_normalized.csv")

def get_significance_genes(metadata_df, df, group1, group2, group1a, group2a):
    abr_genes = list(df.columns)
    abr_genes.pop(0)
    df = metadata_df.merge(df, how='left', right_on='Samples', left_on='Samples#')

    df_CC = df[df['Group']==group1]
    df_GWI = df[df['Group']==group2]
    df_WD_CC=df[df['Group']==group1a]
    df_WD_GWI=df[df['Group']==group2a]

    significant_genes = set()
    for abr in abr_genes:
        p, p_WD = 1,1
        if df_CC[abr].mean()!=0 and df_GWI[abr].mean()!=0:
            stat, p = mannwhitneyu(df_CC[abr], df_GWI[abr])
        if df_WD_CC[abr].mean() != 0 and df_WD_GWI[abr].mean() != 0:
            stat_WD, p_WD = mannwhitneyu(df_WD_CC[abr], df_WD_GWI[abr])
        if p<0.05 or p_WD<0.05:
            print('{}; '+group1a+' p-val: {}; '+group2a+' p-val: {}'.format(abr, p, p_WD))
            significant_genes.add(abr)
    print(significant_genes)
    significant_genes2 = set()
    #print("Chow vs WD")
    for abr in abr_genes:
        p_CC, p_GWI = 1,1
        if df_CC[abr].mean()!=0 and df_WD_CC[abr].mean()!=0:
            stat, p_CC = mannwhitneyu(df_CC[abr], df_WD_CC[abr])
        if df_GWI[abr].mean() != 0 and df_WD_GWI[abr].mean() != 0:
            stat_GWI, p_GWI = mannwhitneyu(df_GWI[abr], df_WD_GWI[abr])
        if p_CC<0.05 or p_GWI<0.05:
            print('{}; '+group1+' p-val: {}; '+group2+' p-val: {}'.format(abr, p_CC, p_GWI))
            significant_genes2.add(abr)
    #print(significant_genes)

    return significant_genes, significant_genes2

class MannWhitneyPValuePlugin:
    def input(self, inputfile):
       self.params=dict()
       infile = open(inputfile, 'r')
       for line in infile:
           contents = line.strip().split('\t')
           self.params[contents[0]] = contents[1]
       self.metadata_df = pd.read_csv(PyPluMA.prefix()+"/"+self.params['metadata'])
       # Check which species are statistically significant between Control and GWI
       #metadata_df = pd.read_csv('../../data/metadata/samples_metagen_metadata.csv')
       self.abundance_species = pd.read_csv(PyPluMA.prefix()+"/"+self.params['abundance'])
       self.group1 = self.params['group1']
       self.group2 = self.params['group2']
       self.group1a = self.params['group1a']
       self.group2a = self.params['group2a']
       #abundaunce_species = pd.read_csv("./out_abundance/abundance_species_transposed.csv")

    def run(self):
       self.sig_genes, self.sig_genes2 = get_significance_genes(self.metadata_df, self.abundance_species, self.group1, self.group2, self.group1a, self.group2a)

    def output(self, outputfile):
        outfile1 = open(outputfile+".coarse.txt", 'w')
        for gene in self.sig_genes:
            outfile1.write(gene+"\n")
        outfile2 = open(outputfile+".fine.txt", 'w')
        for gene in self.sig_genes2:
            outfile2.write(gene+"\n")
#print("Species:")
#get_significance_genes(metadata_df, abundaunce_species)

# print()
# print("Abundance:")
# get_significance_genes(metadata_df, ABR_df_abundance)
