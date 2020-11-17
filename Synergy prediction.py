# -*- coding: utf-8 -*-
"""
Code developed in the laboratory of Cancer Cell Biology,
Engelhardt Institute of Molecular Biology, Moscow, Russia.
If you have any questions about the code or its use contact:
lebedevtd@gmail.com
"""

import pandas as pd
import os


os.chdir('C:\\work folder') #set work folder
input_file='dataset.xlsx'
df= pd.read_excel(input_file)

param= {'Das': {'file': 'Das weights.xlsx', 'cell1': [-1.735, 12.7], 'cell2': [-0.906, 28.1]},
        'ERK': {'file': 'ERK weights.xlsx', 'cell1': [1.614, 7.8], 'cell2': [1.879, 18.5]},
        'IM': {'file': 'IM weights.xlsx', 'cell1': [-7.693, -4.5], 'cell2': [-5.593, 24.4]},
        'Vin': {'file': 'Vin weights.xlsx', 'cell1': [-2.604, 0.1], 'cell2': [-1.846, 13.7]}}

df_out= df
for score in param.keys():
    df_weights= pd.read_excel(param[score]['file'])
    genes_elastic= set(dict.fromkeys(df_weights.Gene.dropna()))
    df_weights.set_index('Gene', inplace=True)
    genes_umap= set(df.columns)
    genes= list(genes_elastic & genes_umap)
    df1= df[genes]
    synergy_list=list()
    for i in df1.index:
        synergy=0
        df2=df1.loc[i]
        for gene in genes:
            synergy += df2[gene]*df_weights.at[gene, 'score']
        synergy_list.append(synergy)
    scale= param[score]
    k= (scale['cell2'][1]-scale['cell1'][1])/(scale['cell2'][0]-scale['cell1'][0])
    y0= scale['cell2'][1]-k*scale['cell2'][0]
    synergy_rescale= [k*x+y0 for x in synergy_list]
    synergy_rescale= [100 if x>100 else x for x in synergy_rescale]
    synergy_rescale= [-100 if x<-100 else x for x in synergy_rescale]
    col_name= 'Synergy_'+score
    df_out[col_name]= synergy_rescale

output_file= input_file.split('.')[0]+'_Synergies_rescale.xlsx'
df_out.to_excel(output_file)


