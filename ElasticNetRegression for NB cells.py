# -*- coding: utf-8 -*-
"""
Code developed in the laboratory of Cancer Cell Biology,
Engelhardt Institute of Molecular Biology, Moscow, Russia.
If you have any questions about the code or its use contact:
lebedevtd@gmail.com
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error
import os
from scipy.stats import zscore

os.chdir('C:\\work folder') #set work folder
files={'test': 'test cells exp with score.xlsx'} #dictionary of dataset files with synergy scores and gene expressions
score='IM_score' #select score for modeling
dataset='test' #select dataset from files dict
df= pd.read_excel(files[dataset])
df_genes= pd.read_excel('Receptors kinases ligands oncogenes unique.xlsx') #select geneset
data_genes= set(dict.fromkeys(df.Cell_line.dropna()))
input_genes= set(dict.fromkeys(df_genes.Gene.dropna()))
df.set_index('Cell_line', inplace=True)
genes=list(input_genes & data_genes)
X=df.loc[genes].transpose()
y=df.transpose()[score]

#used for cross-validation
"""
elastic=ElasticNet(normalize=True, tol=0.005)
search=GridSearchCV(estimator=elastic,param_grid={'alpha':[0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 1],'l1_ratio':[.1,.2,.4,.6,.8]},scoring='neg_mean_squared_error',n_jobs=1,refit=True,cv=5)
search.fit(X,y)
print(search.best_params_)
print(abs(search.best_score_))
"""

a=0.001
l1=0.6
elastic=ElasticNet(normalize=True,alpha=a,l1_ratio=l1, tol=0.005, selection='random')
df_elastic=pd.DataFrame()
for x in range(0, 1000):
    elastic.fit(X,y)
    coef_dict_baseline = {}
    for coef, feat in zip(elastic.coef_,X.columns):
        coef_dict_baseline[feat] = coef
    df_temp = pd.DataFrame.from_dict(coef_dict_baseline, orient='index')
    df_elastic = pd.concat([df_elastic, df_temp], axis=1)
df_out = pd.DataFrame(index=df_elastic.index)
df_out['Result'] = df_elastic.mean(axis=1)
values = np.array(df_out['Result'].dropna())
df_out['z_score']= np.abs(zscore(values))
df_out.sort_values(by='z_score', ascending=False, inplace=True)
if not os.path.exists('ElasticNet results'):
    os.mkdir('ElasticNet results')
df_out.to_excel('ElasticNet results\\Elastic_'+dataset+'_'+'_'+score+'_a'+str(a)+'_l1'+str(l1)+'.xlsx')


