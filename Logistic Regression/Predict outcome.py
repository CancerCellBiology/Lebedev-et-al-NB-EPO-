import numpy as np
import pandas as pd
import os
import math

def predict_survival (df):
    """Calculates event probability (death) based on gene expression from input dataframe df
    and precalculated models (nomycn_model, and mycn_model)"""
    p_list= list()
    for ind1, row1 in df.iterrows():
        """probability for mycn amplified tumors"""
        if df.at[ind1, '#mycn_status']=='yes':
            score=0.075 #intercept parameter
            for ind2, row2 in mycn_model.iterrows():
                gene= mycn_model.at[ind2, 'Gene']
                if gene in df.columns:
                    score= score+ mycn_model.at[ind2,'coef']*df.at[ind1, gene]
            p= 1/(1+math.exp(-score))
        else: 
            """probability for mycn non-amplified tumors"""
            score=0.523 #intercept parameter
            for ind2, row2 in nomycn_model.iterrows():
                gene= nomycn_model.at[ind2, 'Gene']
                if gene in df.columns:
                    score= score+ nomycn_model.at[ind2,'coef']*df.at[ind1, gene]
            p= 1/(1+math.exp(-score))
        p_list.append(p)
    df['p']=p_list
    df['group']=''
    """assing risk groups based on predicted probabilities"""
    for ind, row in df.iterrows():
        if df.at[ind, 'p']>=0.6:
            df.at[ind, 'group']=5
        elif df.at[ind, 'p']>=0.4:
            df.at[ind, 'group']=4
        elif df.at[ind, 'p']>=0.15:
            df.at[ind, 'group']=3
        elif df.at[ind, 'p']>=0.05:
            df.at[ind, 'group']=2
        else:
            df.at[ind, 'group']=1
    """returns df with calculated event probabilities and risk groups"""
    return df 

"""Set work folder with input files"""

os.chdir('C:\\Work_dir\\')
"""list of datasets for survival prediction"""
datasets= {'Cangelosi': 'Tumor Neuroblastoma integrated platforms - Cangelosi - 786_GFgenes_norm.xlsx',
           'NRC': 'NRC_GFgenes_norm.xlsx',
           'Westermann': 'Tumor Neuroblastoma ALT - Westermann - 144 - tpm - gencode19_GFgenes_norm.xlsx',
           'Versteeg': 'Versteeg_GFgenes_norm.xlsx',
           'SEQC': 'SEQC RPM_GFgenes_norm.xlsx'}
nomycn_model= pd.read_excel('LogReg model MYCN_no.xlsx') #load model parameters for mycn non-amplified tumors
mycn_model= pd.read_excel('LogReg model MYCN_yes.xlsx') #load model parameters for mycn amplified tumors
"""calculates survival probability for each dataset"""
for dataset in list(datasets.keys()):
    data=pd.read_excel(datasets[dataset])
    result= predict_survival(data)
    result.to_excel(dataset+'_scores.xlsx') #saves results to excel file


        
