import os
import sys
import pandas as pd

import numpy as np
import matplotlib as mplt
import matplotlib.pyplot as plt
import seaborn as sns

config_txt = sys.argv[1]
out_png = sys.argv[2]

config_dict = dict()
lines = open(config_txt).readlines()
for l in lines:
    splits = l.strip().split()
    print(splits)
    values = splits[1].split(',')
    values.reverse()
    config_dict[splits[0]] = values

gene_counts = []
gene_df = {'SharedNumV' : [], 'SpecificNumV' : []}
for csv in config_dict['GeneFiles']:
    df = pd.read_csv(csv)
    df_v = df.loc[(df['GeneType'] == 'V') & (df['Productive'])]
    df_v_ape = df_v.loc[df_v['IsSpeciesSpecific']]
    num_vs = len(df_v)
    gene_counts.append(num_vs)
    gene_df['SharedNumV'].append(num_vs - len(df_v_ape))
    gene_df['SpecificNumV'].append(len(df_v_ape))

print(np.mean(gene_counts))
gene_df = pd.DataFrame(gene_df, index = config_dict['GenomeNames'])

x = np.array(range(len(gene_counts)))
plt.figure(figsize = (12, 4))

gene_df.plot(stacked=True, kind="barh", color = ['#377DB8', '#E31A1B'], figsize = (4, 12), legend=False)
plt.yticks([], [])
plt.savefig(out_png, dpi = 300)
plt.clf()
