# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 10:21:31 2019

@author: bberg
"""
#%% import packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%% GOI choice
title = "ECM components"
#defnine genes of interest(GOIs) and their annotation.
gois={"comp25324_c0_seq1":"HCol-6|comp25324_c0_seq1",
           "comp27401_c0_seq1":"Laminin, alpha|comp27401_c0_seq1",
           "comp25312_c0_seq1":"Laminin, beta|comp25312_c0_seq1",
           "comp28519_c0_seq1":"HCol-1|comp28519_c0_seq1",
           "comp18600_c0_seq1":"HCol-5|comp18600_c0_seq1",
           "comp24511_c0_seq2":"HCol-3|comp24511_c0_seq2",
           "comp26030_c0_seq1":"Laminin, gamma|comp26030_c0_seq1",
           "comp20380_c0_seq1":"Col IV|comp20380_c0_seq1",
           "comp29438_c0_seq1":"Col VI-like|comp29438_c0_seq1",
           "comp24178_c0_seq1":"HmTSP|comp24178_c0_seq1",
           "comp17112_c0_seq1":"Hemicentin|comp17112_c0_seq1",
           "comp28013_c0_seq1":"Protease Inhibitor|comp28013_c0_seq1",
           "comp24822_c0_seq3":"Dystrophin|comp24822_c0_seq3",
           "comp25485_c0_seq15":"Fibrillin|comp25485_c0_seq15",
           }
#%%plot style
colors = "PRGn_r"

#%%
#read data from file
data =  pd.read_csv("transcriptomic_data.txt",  sep="\t", index_col=0)
#%%
#extract the gene of interet rows from the dataset
df = data.loc[gois.keys()]
#df= data.copy().dropna(axis="rows") #optional cleaning step

#%%
#selects only the samples which have a significant value in one regeneration steps
df_filtered = df.loc[(df["X3h.qvalue"] < 0.05)  | 
        (df["X0.5h.qvalue"] < 0.05)             |
        (df["X6h.qvalue"] < 0.05)               |
        (df["X12h.qvalue"] < 0.05)              |
        (df["X24h.qvalue"] < 0.05)              |
        (df["X48h.qvalue"] < 0.05)].copy()

#Select only the relevant columns part 1       
df_log2=df_filtered[['X0.5h.log2FC','X3h.log2FC','X6h.log2FC','X12h.log2FC',
                    'X24h.log2FC','X48h.log2FC']].copy()
df_log2.columns=["0.5h", "3h","6h","12h","24h","48h"] #Fixes the column names

#replace the IDs with the annotation
new_id=[gois.get(index,index) for index in df_log2.index ]
df_log2["Genes"]=new_id
df_log2=df_log2.set_index("Genes")

#Select only the relevant columns part 2
#df_qvalues=df_filtered[['X0.5h.qvalue','X3h.qvalue','X6h.qvalue','X12h.qvalue',
#                    'X24h.qvalue','X48h.qvalue']].copy()


#%% Build the Plot
g = sns.clustermap(df_log2,col_cluster=False,cmap=colors,center=False)

#styling of the plot
g.fig.suptitle(title)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.subplots_adjust(top=0.985,
            bottom=0.135,
            left=0.04,
            right=0.655,
            hspace=0.2,
            wspace=0.2)
