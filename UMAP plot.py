# -*- coding: utf-8 -*-
"""
Code developed in the laboratory of Cancer Cell Biology,
Engelhardt Institute of Molecular Biology, Moscow, Russia.
If you have any questions about the code or its use contact:
lebedevtd@gmail.com
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import matplotlib



os.chdir('C:\\work folder') #set work folder
file_name='UMAP with synergies' #excel file for plotting
df = pd.read_excel(file_name+'.xlsx')
sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})
matplotlib.rcParams['axes.linewidth']= 5
matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)
custom_color=['#323232', '#ff6000', '#c8def9','#da225f']
cmap= matplotlib.colors.ListedColormap(custom_color)
labels={'cluster': ['na','cl1','cl2', 'cl3'], 'Dataset': ['Kocak','Versteeg', 'NRC']}

#for plotting synergy scores
plt.scatter(df[0], df[1], c=df.Synergy_Vin, cmap='RdYlBu_r', s=20, vmin=-100, vmax=100)
plt.colorbar()

"""
#for plotting clusters/dataset data
data= 'Dataset'
plt.scatter(df[0], df[1], c=df[data], cmap=cmap, s=20)
n= len(set(dict.fromkeys(df[data])))
plt.colorbar(boundaries=np.arange(n+1)-0.5, ticks=np.arange(n)).set_ticklabels(labels[data])
"""

name='Vincristine Synergy' #name the figure
plt.title(name, fontsize=24)
if not os.path.exists('UMAP figures'):
    os.mkdir('UMAP figures')
plt.savefig('UMAP figures//'+file_name+'_'+ name+'.png')