import os
import sys
import pandas as pd
from Bio import SeqIO

locus_dir = sys.argv[1]
locus = sys.argv[2]
locus_csv = sys.argv[3]
output_dir = sys.argv[4]
config_txt = sys.argv[5]

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

df = pd.read_csv(locus_csv)

locus_df = df.loc[df['Locus'] == locus].reset_index()
locus_df['Label'] = [locus_df['Species'][i] + '_' + locus_df['Haplotype'][i] for i in range(len(locus_df))]

order = ['mPanPan1_mat', 'mPanPan1_pat', 'human_T2T', 'mGorGor1_mat', 'mGorGor1_pat', 'mPonAbe1_hap1', 'mPonAbe1_hap2', 'mPonPyg2_hap1', 'mPonPyg2_hap2']

human_dir = '../human_t2t_annotation'
human_fasta = os.path.join(human_dir, 'human_' + locus + '.fasta')

fasta_files = []
for i in range(2):
    fasta_files.append(os.path.join(locus_dir, locus_df.loc[locus_df['Label'] == order[i]].reset_index()['Fname'][0]))

fasta_files.append(human_fasta)

for i in range(3, len(order)):
    fasta_files.append(os.path.join(locus_dir, locus_df.loc[locus_df['Label'] == order[i]].reset_index()['Fname'][0]))

lengths = []
for fasta in fasta_files:
    lengths.append([len(r.seq) for r in SeqIO.parse(fasta, 'fasta')][0])

range_files = []
min_length = 30000

for i in range(len(order) - 1):
    l1 = order[i]
    l2 = order[i + 1]
    fasta1 = fasta_files[i]
    fasta2 = fasta_files[i + 1]
    yass_fname = os.path.join(output_dir, l1 + '_' + l2 + '_yass.txt')
    print(fasta1, fasta2)
    os.system('yass -d 2 -o ' + yass_fname + ' ' + fasta1 + ' ' + fasta2)
    range_txt = os.path.join(output_dir, l1 + '_' + l2 + '_range.txt')
    os.system('python process_alignment.py ' + yass_fname + ' ' + str(min_length) + ' ' + range_txt)
    range_files.append(range_txt)

fh = open(config_txt, 'w')
fh.write('RangeFiles\t' + ','.join(range_files) + '\n')
fh.write('GenomeNames\t' + ','.join(order) + '\n')
fh.write('Lengths\t' + ','.join([str(l) for l in lengths]))
fh.close()
