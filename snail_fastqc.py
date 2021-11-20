# Importing modules

import argparse
import os
import csv

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

print('Calculated basic Statistics...')

