import pandas as pd
import os
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import numpy as np
from statistics import mean
from scipy.stats import zscore
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import matplotlib

os.chdir('C:\\Lab\\Python\\NB\\LogReg\\')

datasets= {'SEQC': ['SEQC_scores_v4.xlsx', '#323232'],
           'Versteeg': ['Versteeg_scores_v4.xlsx', '#5baf06'],
           'NRC': ['NRC_scores_v4.xlsx', '#dddd73'],
           'Westermann': ['Westermann_scores_v4.xlsx','#790079'],
           }
matplotlib.rcParams['axes.linewidth']= 3
fig, ax = plt.subplots(figsize=(10, 10))
for data in datasets.keys():
    df= pd.read_excel(datasets[data][0])
    df= df[df['#event_overall']!='na']
    df['event']= df['#event_overall'].apply(lambda x: 1 if x== 'yes' else 0)
    df1=df
    #df1= df[df['#mycn_status']=='yes']
    event= df1['event']
    log_p= df1['p']
    # calculate ROC
    ns_p = [0 for _ in range(len(event))]
    ns_auc = roc_auc_score(event, ns_p)
    log_auc = roc_auc_score(event, log_p)
    ns_fpr, ns_tpr, _ = roc_curve(event, ns_p)
    lr_fpr, lr_tpr, _ = roc_curve(event, log_p)
    # plot AUC
    ax.plot(lr_fpr, lr_tpr, marker='.', linewidth=6, color=datasets[data][1], label=data+' AUC= %.2f' % log_auc)
ax.plot(ns_fpr, ns_tpr, linestyle='--', linewidth=6, color='#4b4b4b')
ax.set_xlabel('False Positive Rate', fontsize=24, weight='bold')
ax.set_ylabel('True Positive Rate', rotation=90, fontsize=24, weight='bold')
ax.tick_params(labelsize=18)
# add parameters (F1 or AUC)
#ax.text(0, 0.95, 'AUC= %.2f' % log_auc, fontsize=24, weight='bold')
# show the legend
plt.legend(prop={'size': 24}, loc=(0.32, 0.0))
#plt.legend(prop={'size': 24}, loc=(0.36, -0.5))
plt.savefig('LogReg_datasets_ROC_AUC_v4.png', bbox_inches = 'tight')