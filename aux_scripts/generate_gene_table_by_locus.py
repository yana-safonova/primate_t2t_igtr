import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def IsGeneSpeciesSpecific(gene_seq, species_specific_genes):
    for g in species_specific_genes:
        gene_dir = g
        gene_rev = str(Seq(g).reverse_complement())
        if gene_seq.find(gene_dir) != -1 or gene_seq.find(gene_rev) != -1:
            return True
    return False    

genome_fasta = sys.argv[1]
bed_fname = sys.argv[2]
locus = sys.argv[3]
gene_csv = sys.argv[4]

gene_types = ['V', 'D', 'J']
contig_name = ''
lines = open(bed_fname).readlines()
genes = []
for l in lines:
    splits = l.strip().split()
    contig_name = splits[0]
    start_pos = int(splits[1]) - 1
    end_pos = int(splits[2]) - 1
    gene_name = splits[3]
    type_of_gene = ''
    for gene_type in gene_types:
        gene_prefix = locus + gene_type
        if gene_name.find(gene_prefix) != -1 and len(gene_name) > 4:
            type_of_gene = gene_type
            break
    if type_of_gene == '':
        continue
    genes.append((gene_name, type_of_gene, start_pos, end_pos))
    print(genes[-1])

contig = ''
for r in SeqIO.parse(genome_fasta, 'fasta'):
    if r.id == contig_name:
        contig = str(r.seq).upper()
        break

min_pos = min([p[2] for p in genes])
max_pos = max([p[2] for p in genes])

offset = 20000
adj_min_pos = max(0, min_pos - offset)
adj_max_pos = min(len(contig), max_pos + offset)

locus_seq = contig[adj_min_pos : adj_max_pos]

gene_df = pd.read_csv(gene_csv)
species_specific_genes = list(gene_df.loc[gene_df['IsSpeciesSpecific']]['Sequence'])
species_specific_pos = list(gene_df.loc[gene_df['IsSpeciesSpecific']]['Pos'])
for gene in species_specific_genes:
    continue
    gene_rev = str(Seq(gene).reverse_complement())
    pos_dir = locus_seq.find(gene)
    pos_rev = locus_seq.find(gene_rev)
    print(gene, pos_dir, pos_rev)
    if max(pos_dir, pos_rev) == -1:
        continue
    species_specific_pos.append(max(pos_dir, pos_rev))
print(len(species_specific_genes), len(species_specific_pos))    

df = {'GeneType' : [], 'Contig' : [], 'Pos' : [], 'Strand' : [], 'Sequence' : [], 'Productive' : [], 'Locus' : [], 'LocusPos' : [], 'IsSpeciesSpecific' : []}
for gene_name, gene_type, start_pos, end_pos in genes:
    if gene_type == 'V':
        continue
    df['GeneType'].append(gene_type)
    df['Contig'].append(contig_name)
    df['Pos'].append(start_pos)
    df['Strand'].append('NA')
    sequence = contig[start_pos : end_pos]
    df['Sequence'].append(sequence)
    df['Productive'].append(True)
    df['Locus'].append(locus)
    df['LocusPos'].append(start_pos - adj_min_pos)
    is_species_specific = IsGeneSpeciesSpecific(sequence, species_specific_genes) #len([pos for pos in species_specific_pos if pos >= (start_pos - adj_min_pos) and pos <= (end_pos - adj_min_pos)]) != 0
    df['IsSpeciesSpecific'].append(is_species_specific)

df = pd.DataFrame(df)
v_df = gene_df.loc[gene_df['GeneType'] == 'V']
df = pd.concat([df, v_df])
df = df.sort_values(by = 'Pos').reset_index()

min_pos = min(df['Pos'])
max_pos = max(df['Pos'])

adj_min_pos = max(0, min_pos - offset)
adj_max_pos = min(len(contig), max_pos + offset)

df['LocusPos'] = [df['Pos'][i] - adj_min_pos for i in range(len(df))]
df.to_csv('human_' + locus + '.csv', index = False, columns = ['GeneType', 'Contig', 'Pos', 'Strand', 'Sequence', 'Productive', 'Locus', 'LocusPos', 'IsSpeciesSpecific'])

fh = open('human_' + locus + '.fasta', 'w')
fh.write('>' + contig_name + '|START_POS:' + str(adj_min_pos) + '|END_POS:' + str(adj_max_pos) + '\n' + contig[adj_min_pos : adj_max_pos] + '\n')
fh.close()
