# Importing modules

import argparse
import os

# Parser for command line options
from argparse import ArgumentParser

parser = ArgumentParser(description="FastQC analogue")
parser.add_argument("-i", "--input", dest="filename", required=True, help="input .fastqc file", metavar="FILE")
parser.add_argument("-o", "--outdir", dest="outdir", required=True, help="input the directory for the output file")
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

# Setting the working directory
os.chdir(outdir)

