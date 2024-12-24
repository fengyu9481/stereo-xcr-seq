import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--raw', type=str, required=True, help='Path to the cyclone raw FASTQ file')
    
parser.add_argument('--step1', type=str,required=True, help='Path to the step 1 R2 FASTQ file')
    
parser.add_argument('--out', type=str, required=True, help='Path to the output FASTA file')

args = parser.parse_args()

cyclone_raw_fastq = args.raw
step1_fastq = args.step1
filename = args.out

import pandas as pd
import gzip

def read_fastq(fastq):
    rnames = list()
    seqs = list()
    qs = list()
    if fastq.endswith('.gz'):
        f = gzip.open(fastq, 'rt')
    else:
        f = open(fastq, 'r')
    for idx, l in enumerate(f):
        if idx % 4 == 0:
            rnames.append(l.split()[0][1:])
        if idx % 4 == 1:
            seqs.append(l.rstrip())
        if idx % 4 == 3:
            qs.append(l.rstrip())
    return rnames, seqs,qs

def rev_compl(s):  
    rev_compl_l = [chr(i) for i in range(128)]
    rev_compl_l[ord('A')] = 'T'
    rev_compl_l[ord('C')] = 'G'
    rev_compl_l[ord('G')] = 'C'
    rev_compl_l[ord('T')] = 'A'
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

rawrn,rawsq,rawqs = read_fastq(cyclone_raw_fastq)
fixrn,fixsq,fixqs = read_fastq(step1_fastq)
fixrn = [ i.split('|||')[0] for i in fixrn]
rawdict = dict(zip(rawrn,rawsq))
fixdict = dict(zip(fixrn,fixsq))
with gzip.open(filename, 'wt') as file:
    for i in range(len(fixrn)):
        raw_fastq = rawdict[fixrn[i]]
        fix_fastq = fixdict[fixrn[i]]
        idx = raw_fastq.index(fix_fastq)
        if idx == 0:
            rest = raw_fastq[len(fix_fastq):]
            rest = rev_compl(rest)
        else:
            rest = raw_fastq[:idx]
        file.write(f">{fixrn[i]}\n")
        file.write(f"{rest}\n")


