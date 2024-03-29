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

def find_param(df1, genes):
    """Finds optimal parameters for logistic regression by going through l1 and c parameters
    and then compares them using mean for random 5 different datasets splits
    df1- dataframe with input parameters, genes- list of genes used for log regression"""
    X= df1[genes]
    y= df1['event']
    l1_list=list()
    c_list=list()
    score_list=list()
    for l1 in np.arange(0,1, 0.2):
        for c in [1,10,100,1000]:
            test_score_list=list()
            for state in np.arange(0,5):
                X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=split, random_state=state)
                log_model = LogisticRegression(C=c, solver='saga', random_state=0, penalty='elasticnet', l1_ratio=l1, max_iter= 10000).fit(X_train,y_train)
                test_score = log_model.score(X_test,y_test)
                test_score_list.append(test_score)
            l1_list.append(l1)
            c_list.append(c)
            score_list.append(mean(test_score_list))
    df_param= pd.DataFrame({'l1':l1_list, 'c':c_list, 'score':score_list})
    df_param.sort_values('score', inplace=True, ascending=False)
    df_param.reset_index(inplace=True)
    """returns l1 and c parameters for log regression"""
    return df_param.at[0,'l1'], df_param.at[0,'c']

def log_reg (df1, l1, c, genes):
    """performs logistic regression for data from input dataframe (df1) using precalculated l1 and c parameters 
    for a subset from imput genes list (genes)"""
    X= df1[genes]
    y= df1['event']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=split, random_state=0)
    log_model = LogisticRegression(C=c, solver='saga', random_state=0, penalty='elasticnet', l1_ratio=l1, max_iter= 10000).fit(X_train, y_train)
    df_coef= pd.DataFrame({'Gene': list(X.columns), 'coef':log_model.coef_[0], 'z-score': np.abs(zscore(log_model.coef_[0]))}) 
    df_coef.sort_values('z-score', inplace=True, ascending=False)
    """returns log model, intercept parameter and model parameters (gene coefficients) in df_coef"""
    return log_model, log_model.intercept_, df_coef
    
"""Set work folder and load dataset with survival data and genes used for analysis""" 

os.chdir('C:\\Work_dir\\')
data = pd.read_excel('Tumor Neuroblastoma integrated platforms - Cangelosi - 786_GFgenes_norm.xlsx')
df_set = pd.read_excel('Diff genes.xlsx')
genes_valid= list(set(data.columns) & set(df_set['Gene']))
data['event']= data['#event_overall'].apply(lambda x: 1 if x== 'yes' else 0) #transforms events into 1- dead, and 0- alive
"""Set parameters for dataset splitting and MYCN status"""
status='no' #yes for MYCN amp model
split=0.3 #select ratio for dataset split (30% Test dataset and 70% Training dataset)
df= data[data['#mycn_status']==status]
"""Calculates optimal parameters for initial model, containing all genes""" 
L10, C0= find_param(df, genes_valid)
print('Initial model parameters for mycn_status='+status+':')
print('L1_ratio:'+str(L10)+' C:'+str(C0))
"""Applies parameters to create the first model"""
model0, inter0, model0_table= log_reg(df, L10, C0, genes_valid)
"""Selects top 20% genes, based on their absolute score in the first model
   to change the percentage change 0.2 to any number between 0 and 1"""
genes_select= list(model0_table.head(round(len(genes_valid)*0.2))['Gene'])
print('Genes selected:')
print(genes_select)

L1, C= find_param(df, genes_valid)
model, inter, model_table= log_reg(df, L1, C, genes_select)
print('Final model parameters for mycn_status='+status+':')
print('L1_ratio='+str(L1)+' C='+str(C)+' intercept='+str(inter[0]))
model_table.to_excel('LogReg model MYCN_'+status+'.xlsx')
"""Calculates optimal parameters for second model, containing only top genes, and creates final model""" 
Exp= df[genes_select]
event= df['event']
Exp_train, Exp_test, event_train, event_test = train_test_split(Exp, event, test_size=split, random_state=0)
print('Train:'+classification_report(event_train, model.predict(Exp_train)))
print('Test:'+classification_report(event_test, model.predict(Exp_test)))
"""Tests created model Precision/Recall for Test dataset"""
log_p = model.predict_proba(Exp_test)
log_p = log_p[:, 1]
yhat = model.predict(Exp_test)
lr_precision, lr_recall, _ = precision_recall_curve(event_test, log_p)
lr_f1, lr_auc = f1_score(event_test, yhat), auc(lr_recall, lr_precision)
no_skill = len(event_test[event_test==1]) / len(event_test)
"""Plot results for Precision/Recall"""
matplotlib.rcParams['axes.linewidth']= 3
fig, ax = plt.subplots(figsize=(10, 10))
plt.plot([0, 1], [no_skill, no_skill], linestyle='--', linewidth=6, color='k', label='No Skill')
plt.plot(lr_recall, lr_precision, marker='.', linewidth=6, color='#7fc355', label='MYCN '+status+'-amp model')
ax.set_xlabel('Recall', fontsize=24, weight='bold')
ax.set_ylabel('Precision', rotation=90, fontsize=24, weight='bold')
ax.tick_params(labelsize=18)
plt.legend(prop={'size': 24}) # show the legend
plt.savefig('LogReg_MYCN_'+status+'_test_PreRec.png', bbox_inches = 'tight')




