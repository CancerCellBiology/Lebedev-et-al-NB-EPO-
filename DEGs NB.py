import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os

def diff_genes(data, cl, status):
    if status== 'withMYCN':
        df= data 
    elif status== 'noMYCN':
        df= data[data['#mycn']!=2]
    df1= df[df['cluster']==cl].drop(['Dataset', 'Dataset_n', '#mycn', '#stage', 'cluster', 0, 1], axis=1)
    df2= df[df['cluster']!=cl].drop(['Dataset', 'Dataset_n','#mycn', '#stage', 'cluster', 0, 1], axis=1)
    genes= set(df1.drop('Patient', axis=1).columns)
    genes_set= set(dict.fromkeys(df_set['Gene'].dropna()))
    df1.set_index('Patient', inplace=True)
    df2.set_index('Patient', inplace=True)
    gene_list= list(genes & genes_set)
    df1_t = df1.transpose()
    df2_t = df2.transpose()
    df1_mean= df1_t.mean(axis=1)
    df2_mean= df2_t.mean(axis=1)
    p_list=list()
    dif_list=list()
    mean_1_list=list()
    mean_2_list=list()
    for gene in gene_list:
        stat, p= mannwhitneyu(df1_t.loc[gene], df2_t.loc[gene], alternative='two-sided')
        mean_1= df1_mean.loc[gene]
        mean_2= df2_mean.loc[gene]
        dif= mean_1 / mean_2
        p_list.append(p)
        dif_list.append(dif)
        mean_1_list.append(mean_1)
        mean_2_list.append(mean_2)
    p_adj= multipletests(pvals=p_list, alpha=0.01, method='fdr_tsbky')
    df_out= pd.DataFrame({'Gene': gene_list, 'cl_mean': mean_1_list, 'other_mean': mean_2_list, 'p-value': p_list, 'p-adj': p_adj[1],
                          'ratio': dif_list, 'significant': p_adj[0]})
    df_out.sort_values('p-adj', inplace=True)
    df_out.to_excel(input_file.split('.')[0]+ '_Diff_genes_cl'+str(cl)+'_'+status+'.xlsx')
    return



os.chdir('C:\\Work_dir\\')
input_file='UMAP_NB_v2_noSEQC_correlation_50_v2_no-outliners.xlsx'
df_umap= pd.read_excel(input_file)
file_name= input_file.split('.')[0].split('_')
df_set= pd.read_excel('Genes for UMAP Elastic GF only signaling Erbb.xlsx')
for cluster in list(dict.fromkeys(df_umap['cluster'])):
    diff_genes(df_umap, cluster, 'noMYCN')
    diff_genes(df_umap, cluster, 'withMYCN')




