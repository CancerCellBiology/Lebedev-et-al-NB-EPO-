import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
import matplotlib.colors as colors



os.chdir('C:\\Lab\\Python\\NB\\UMAP')
file_name='UMAP_NB_v2_noSEQC_correlation_50_v2_no-outliners'
df = pd.read_excel(file_name+'.xlsx')
fig, ax = plt.subplots(figsize=(10, 10))
matplotlib.rcParams['axes.linewidth']= 5
custom_color=['#323232', '#9F6DAC', '#FBC570','#8AC7B9', '#ab2b6f']
custom_color1=['#504f4f', '#7d2a57', '#c8def9','#a0cf4c'] #clusters 
custom_color2=['#504f4f', '#7d2a57', '#c8def9','#a0cf4c', '#da225f']
custom_color3=['#323232', '#c8def9', '#da225f']
custom_color4=['#323232', '#a0cf4c', '#c8def9','#da225f', '#782c8c', '#ff6000']
cmap= matplotlib.colors.ListedColormap(custom_color3)

labels= {'cluster': ['na','cl1','cl2', 'cl3'], 'Dataset_n': ['Kocak', 'Maris', 'NRC', 'Versteeg', 'Westermann'],
         '#mycn':['na','non-amp','amp'], '#stage': ['na', 'st 1', 'st 2', 'st 3', 'st 4', 'st 4S']}
label= '#mycn'
"""Plot discrete values"""
sc= ax.scatter(df[0], df[1], c=df[label], cmap=cmap, s=20)
n= len(set(dict.fromkeys(df[label])))
#plt.colorbar(sc, boundaries=np.arange(n+1)-0.5, ticks=np.arange(n)).set_ticklabels(labels[label])
"""Plot expression"""
#df['log_label']= np.log2(df[label])
#sc= ax.scatter(df[0], df[1], c=df['log_label'], cmap='RdYlBu_r', s=20, vmin=-2, vmax=2)
plt.xticks([], [])
plt.yticks([], [])
ax.set_xlabel('UMAP1', fontsize=18, weight='bold')
ax.set_ylabel('UMAP2', fontsize=18, weight='bold')

# legend for expression plots
"""
cbar= fig.colorbar(sc, shrink=0.5, aspect=8)
cbar.ax.tick_params(labelsize=18)
cbar.set_label('log(FC)', fontsize=18, rotation=90)
"""

#plt.colorbar(boundaries=np.arange(7)-0.5, ticks=np.arange(6))
name='UMAP NB v2 n0 SEQC'+label+' all samples'
#plot_name= ''
#plt.title(plot_name, fontsize=24, weight='bold')
#save_name=input('Figure name:')
plt.savefig('new//'+name+'.png', bbox_inches = 'tight')