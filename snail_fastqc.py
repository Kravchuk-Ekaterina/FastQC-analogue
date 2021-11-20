# Importing modules

import argparse
import os
import csv
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Parser for command line options
from argparse import ArgumentParser

parser = ArgumentParser(description="FastQC analogue")
parser.add_argument("-i", "--input", dest="filename", required=True, help="input .fastqc file", metavar="FILE")
parser.add_argument("-o", "--outdir", dest="outdir", required=True, help="input the directory for the output file", metavar = 'OUTDIR')
args = parser.parse_args()

input_file = args.filename
outdir = args.outdir
if not os.path.isdir(outdir):
     os.mkdir(outdir)

# Reading the input file

lines =[]
with open(input_file) as inp:
    for line in inp:
        lines.append(line)
# Descarding metadata if there is metadata

n = 0
while True:
    seq = lines[n]
    if seq[0] == '@':
        break
    n += 1

file_len = len(lines)

seq_lines = lines[1::4]
qual_lines = lines[3::4]

# Setting the working directory

os.chdir(outdir)

# Basic statistics

f = input_file.split('/')
Filename = (f[-1])[0:len(f[-1])-6]

File_type = 'Conventional base calls'

"""
for Illumina 1.3-1.8 quality starts whith ASCII 64
for Solexa / Illumina 1.0 quality starts whith ASCII 59
for Sanger / Illumina 1.9 quality starts whith ASCII 33
"""

Encoding = 'Illumina 1.3-1.8'

for line in qual_lines:
    if '?' in line:
        Encoding = 'Solexa / Illumina 1.0'
        break

for line in qual_lines:
    if ':' in line:
        Encoding = 'Sanger / Illumina 1.9'
        break

Total_Sequences = round(file_len/4)

len_list = []
for line in seq_lines:
    len_list.append(len(line)-1)

max_len = max(len_list)
min_len = min(len_list)

Sequence_Length = 'From {min_len} to {max_len}'

if min_len == max_len:
    Sequence_Length = min_len

gc = 0
for line in seq_lines:
    gc += line.count('G')
    gc += line.count('C')
GC = round(100*(gc/sum(len_list)))

with open(Filename + '_Basic_Statistics.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Filename', Filename + '.fastq'])
    tsv_writer.writerow(['File type', File_type])
    tsv_writer.writerow(['Encoding', Encoding])
    tsv_writer.writerow(['Total Sequences', Total_Sequences])
    tsv_writer.writerow(['Sequence Length', Sequence_Length])
    tsv_writer.writerow(['%GC', GC])

print('Calculated Basic Statistics...')

# Per Sequence Quality Scores

def average_score(lst):
    plus = 0
    qual_list = []
    if Encoding == 'Solexa / Illumina 1.0':
        plus = 26
    if Encoding == 'Illumina 1.3-1.8':
        plus = 31
    delta = 33 + plus
    for line in lst:
        quality = line.strip()
        sum_q = 0
        for j in quality:
            sum_q += (ord(j) - delta)
        mean_quality = sum_q/len(quality)
        qual_list.append(mean_quality)
    return qual_list


qual_list = average_score(qual_lines)

med = np.median(qual_list)
status = 'Green (entirely normal)'
col ='green'
if med < 27:
    status = 'Orange (slightly abnormal)'
    col = 'orange'
if med < 20:
    status = 'Red (very unusual)'
    col = 'red'

x, y = np.unique(qual_list, return_counts = True)
plt.plot(x, y, color = col)
plt.title('Per Sequence Quality Scores Destribution')
plt.suptitle(status)
plt.savefig(Filename + '_per_seq_qscore.png')

print('Calculated Per Sequence Quality Scores... ')

# Overrepresented sequences

count_seqs = {}
unique_seqs = list(set(seq_lines))
for line in unique_seqs[0:200]:
    n = seq_lines.count(line)
    count_seqs[line] = n

status = 'Green (entirely normal)'
over_seqs = {}
len_seq_lines = len(seq_lines)
for key in count_seqs:
    if count_seqs[key]/len_seq_lines > 0.001:
        status = 'Orange (slightly abnormal)'
        over_seqs[key] = count_seqs[key]
    if n/len_seq_lines > 0.01:
        status = 'Red (very unusual)'

with open(Filename + '_Overrepresented_sequences.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Status', status])
    if status != 'Green (entirely normal)':
        tsv_writer.writerow(['Sequence', 'Count', 'Percentage'])
        for key in over_seqs:
            tsv_writer.writerow([key, over_seqs[key], 100*over_seqs[key]/len_seq_lines])

print('Finished searching for overrepresented sequences...')

