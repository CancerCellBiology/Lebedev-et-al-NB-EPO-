import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os


'''set working folder name'''
os.chdir('work folder') 
'''select file with clustering results from working folder'''
input_file='NB_tumors_clusters.xlsx' 
'''select cluster'''
cluster=1 
df= pd.read_excel(input_file)
genes = list(df.drop(['Patient','Dataset', 'mycn', 'stage', 'cluster', 0, 1], axis=1).columns) #exclude annotations
df1= df[df['cluster']==cluster][genes]
df2= df[df['cluster']!=cluster][genes]
df1_t = df1.transpose()
df2_t = df2.transpose()
df1_mean= df1_t.mean(axis=1)
df2_mean= df2_t.mean(axis=1)
p_list=list()
dif_list=list()
mean_1_list=list()
mean_2_list=list()
'''perform t-tests and calculate difference (ratio)'''
for gene in genes:
    stat, p= mannwhitneyu(df1_t.loc[gene], df2_t.loc[gene], alternative='two-sided')
    mean_1= df1_mean.loc[gene]
    mean_2= df2_mean.loc[gene]
    dif= mean_1 / mean_2
    p_list.append(p)
    dif_list.append(dif)
    mean_1_list.append(mean_1)
    mean_2_list.append(mean_2)
'''perform FDR correction'''
p_adj= multipletests(pvals=p_list, alpha=0.01, method='fdr_tsbky') 
df_out= pd.DataFrame({'Gene': genes, 'cl_mean': mean_1_list, 'other_mean': mean_2_list, 'p-value': p_list, 'p-adj': p_adj[1],
                      'ratio': dif_list, 'significant': p_adj[0]})
output_file= input_file.split('.')[0]+ '_Diff_genes_cl'+str(cluster)+'vs_all.xlsx'
df_out.to_excel(output_file) #save results
