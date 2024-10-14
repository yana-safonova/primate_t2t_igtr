import os
import sys
import pandas as pd
import random as rd
from pygenomeviz import GenomeViz

import numpy as np
import matplotlib as mplt
import matplotlib.pyplot as plt
import seaborn as sns

def GetColorByNormalizedValue(cmap_name, norm_value):
    if norm_value < 0 or norm_value > 1:
        print("ERROR: value " + str(norm_value) + ' does not belong to [0, 1]')
    cmap = plt.cm.get_cmap(cmap_name)
    color = cmap(norm_value)
    return mplt.colors.rgb2hex(color[:3])

def GetPercentIdentityColor(pi):
    pi = pi / 100
    if pi < 0.9:
        return '#96C6D3'
    if pi < 0.95:
        return '#A5CDD6'
    if pi < 0.99:
        return '#B5D4DB'
    if pi < 0.995:
        return '#D1D9B4'
    if pi < 0.996:
        return '#EEDF8E'
    if pi < 0.997:
        return '#EBD783'
    if pi < 0.998:
        return '#E9D079'
    if pi < 0.999:
        return '#EDAB79'
    return '#F28679'

def AdjustPos(pos, strand, l):
    if strand == 1:
        return pos
    else:
        return l - pos - 1

config_txt = sys.argv[1]
out_png = sys.argv[2]

config_dict = dict()
lines = open(config_txt).readlines()
for l in lines:
    splits = l.strip().split()
    print(splits)
    config_dict[splits[0]] = splits[1].split(',')

alignment_files = config_dict['RangeFiles']
dfs = [pd.read_csv(f, sep = '\t') for f in alignment_files]

genomes = config_dict['GenomeNames']
lens = [int(l) for l in config_dict['Lengths']]

gene_files = config_dict['GeneFiles']
gene_dfs  = [pd.read_csv(f) for f in gene_files]

#bp_files = config_dict['BPFiles']
#bp_dfs = [pd.read_csv(f, sep = '\t') for f in bp_files]

unit_files = []
unit_dfs = [pd.DataFrame({})] * len(genomes)
if 'Units' in config_dict:
    unit_files = config_dict['Units']
    unit_dfs = [pd.read_csv(f, sep = '\t') for f in unit_files]

strands = [int(s) for s in config_dict['Strands']]

gv = GenomeViz(tick_style="axis")
#gv.set_scale_xticks()

tracks = []
gene_colors = {'V' : '#1F77B4', 'D' : '#FF7F0F', 'J' : '#2AA02B'}
strand_dict = {'+' : 1, '-' : -1}
gene_width = max(lens) / 500
for g, l, gene_df, locus_strand in zip(genomes, lens, gene_dfs, strands):
    # adding genes
    track = gv.add_feature_track(g, l)
    tracks.append(track)

    for i in range(len(gene_df)):
        pos = AdjustPos(gene_df['LocusPos'][i], locus_strand, l)
        end_pos = AdjustPos(min(l, gene_df['LocusPos'][i] + gene_width), locus_strand, l)
        species_specific = False
        if 'IsSpeciesSpecific' in gene_df.columns:
            species_specific = gene_df['IsSpeciesSpecific'][i]
        color = gene_colors[gene_df['GeneType'][i]]
        if species_specific:
            color = 'red'
        track.add_feature(min(pos, end_pos), max(pos, end_pos), strand = 1, plotstyle = 'box', facecolor = color)

min_pi = 100
max_pi = 0
for df in dfs:
    min_pi = min(min_pi, min(df['PI']))
    max_pi = max(max_pi, max(df['PI']))

print(min_pi, max_pi)
min_pi = 90

for idx, df in enumerate(dfs):
    l1 = genomes[idx]
    l2 = genomes[idx + 1]
    s1 = strands[idx]
    s2 = strands[idx + 1]
    len1 = lens[idx]
    len2 = lens[idx + 1]
    print('==== ' + l1 + ' ' + l2)
    for i in range(len(df)):
        p11, p12 = AdjustPos(df['G1_start'][i] - 1, s1, len1), AdjustPos(df['G1_end'][i] - 1, s1, len1)
        p21, p22 = AdjustPos(df['G2_start'][i] - 1, s2, len2), AdjustPos(df['G2_end'][i] - 1, s2, len2)
        pi = max(min_pi, df['PI'][i])
        pi_color = GetPercentIdentityColor(pi) #GetColorByNormalizedValue('coolwarm', (pi - min_pi) / (max_pi - min_pi))
        print(p11, p12, p21, p22, pi)
        gv.add_link((l1, p11, p12), (l2, p21, p22), normal_color=pi_color, inverted_color=pi_color, curve=True)

max_cov = 8 #max([max(df['AvgCov']) for df in bp_dfs])
fig = gv.plotfig()

width_pos = 0.5
width_neg = 0.5
for g, l, t, unit_df, strand in zip(genomes, lens, tracks, unit_dfs, strands):
    t = gv.get_track(g)
    # genome box
    x, y = (0, 0, l, l), (width_pos * 2, 0, 0, width_pos * 2)
    t.ax.fill(x, y, fc="#C7C7C7", linewidth=1, ec = 'black', alpha=0.5, zorder=-10)
    # break point track
    #for i in range(len(bp_df)):
    #    continue
    #    start = AdjustPos(bp_df['Start'][i], strand, l)
    #    end = AdjustPos(max(bp_df['End'][i], start + 1000), strand, l)
    #    x, y = (start, start, end, end), (-width_neg, 0, 0, -width_neg)
    #    color = GetColorByNormalizedValue('Greens', min(0.75, bp_df['AvgCov'][i] / max_cov))
    #    t.ax.fill(x, y, fc = color, linewidth = 0, ec = 'black', alpha=1, zorder=-10)
    for i in range(len(unit_df)):
        start = AdjustPos(unit_df['Start'][i], strand, l)
        end = AdjustPos(unit_df['End'][i], strand, l)
        color = unit_df['Color'][i]
        x, y = (start, start, end, end), (-width_neg, 0, 0, -width_neg)
        t.ax.fill(x, y, fc = color, linewidth = 1, ec = color, alpha=1, zorder=-10)

fig.savefig(out_png)
