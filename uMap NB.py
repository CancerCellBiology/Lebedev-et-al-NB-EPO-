import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import os
import hdbscan
import matplotlib



'''set working folder name'''
os.chdir('work folder path')
'''select file with clustering results from working folder'''
input_file='NB_tumors_normalized_expression.xlsx' 
df= pd.read_excel(input_file)
file_name= input_file.split('.')[0].split('_')
'''set UMAP parameters'''
param= {'n': 50, 'method': 'correlation'}
df = pd.read_excel(file_name+'.xlsx')
genes = set(df.drop(['Patient', 'Dataset', 'mycn', 'stage'], axis=1).columns) #exclude annotations
'''perform UMAP embedding and HDBSCAN'''
embedding = umap.UMAP(n_neighbors=param['n'], min_dist=0.0, metric=param['method']).fit_transform(df[genes]) 
hdbscan_labels = hdbscan.HDBSCAN(min_samples=50, min_cluster_size=100).fit_predict(embedding) 
hdbscan_labels = hdbscan_labels+1
df1 = pd.DataFrame(embedding)
df1['cluster'] = hdbscan_labels
df_out= pd.concat([df, df1], axis=1)
df_out.set_index('Patient', inplace=True)
df_out.to_excel(file_name+'_UMAP_'+'n='+str(param['n'])+'_'+param['method']+'.xlsx') #save UMAP results
'''plot results with cluster annotation'''
sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})
custom_color=['#323232', '#a0cf4c', '#c8def9','#da225f', '#782c8c', '#ff6000']
cmap= matplotlib.colors.ListedColormap(custom_color)
plt.scatter(df1[0], df1[1], c=df1.cluster, cmap=cmap, s=20)
n= len(set(dict.fromkeys(df1.cluster)))
plt.colorbar(boundaries=np.arange(n+1)-0.5, ticks=np.arange(n))

