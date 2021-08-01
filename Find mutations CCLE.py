import pandas as pd
import os
import glob

os.chdir('C:\\Lab\\Python\\NB\\Mutations')
path= 'C:\\Lab\\Python\\NB\\NB CCLE'
files= glob.glob(path+'//*.csv')
mut_list= list()
for file in files:
    df= pd.read_csv(file)
    for ind, row in df.iterrows():
        df.at[ind,'Mutation']= df.at[ind,'Gene']+'_'+df.at[ind, 'Variant Type']
    mut_list+= list(set(df[df['Variant Annotation']!= 'silent']['Gene']))
mutations= list(set(mut_list))
counts= list()
for mut in mutations:
    counts.append(mut_list.count(mut))
df_out= pd.DataFrame({'Mutation': mutations, 'Occurance': counts})
df_out.sort_values('Occurance', ascending=False, inplace=True)
df_out.to_excel('NB CCLE non silent mutations by gene.xlsx')



path= 'C:\\Lab\\Python\\NB\\Mutations\\NB CCLE\\selected'
files= glob.glob(path+'//*.csv')
df_rate= pd.read_excel('TARGET vs CCLE mutation rate.xlsx')
selected= list(df_rate[df_rate['Target_rate']>=0.02]['Gene'])
selected+= ['NRAS', 'NF1', 'KRAS', 'TP53', 'RET', 'ABCA13', 'ABCB1', 'MET', 'EGFR']
df_cell= pd.DataFrame()
mut_list= list()
for gene in selected:
    cell_list= list()
    for file in files:
        df= pd.read_csv(file)
        cell_line= file.split('\\')[-1].split(' ')[0]
        mut_list= list(set(df[df['Variant Annotation']!= 'silent']['Gene']))
        if gene in mut_list:
            cell_list.append(cell_line)
    df_temp= pd.DataFrame({gene: cell_list})
    df_cell= pd.concat([df_cell, df_temp], axis=1)
df_cell.to_excel('NB mutations selected cells expanded.xlsx')
