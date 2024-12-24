import pandas as pd
import gzip
from glob import glob
import os
import sys 

def read_fastq(cln,fastq):
    rnames = list()
    if fastq.endswith('.gz'):
        f = gzip.open(fastq, 'rt')
    else:
        f = open(fastq, 'r')
    for idx, l in enumerate(f):
        if idx % 4 == 0:
            rnames.append([cln,l.split()[0][1:]])
    return rnames

sampleid = sys.argv[1]
mixcr_align = sys.argv[2]  # "./mc38_0304/mc38_0304.align.tsv"
fastq_need_subset = sys.argv[3]    #"/storage/liuyi/09.ma_tcr/p1/split_fix/MC38-TRBC_mixcr_CDR3.r2.fastq.gz"
clone_info = sys.argv[4] #"./mc38_0304/mc38_0304.contigs.tsv"
outdir = sys.argv[5]
barcode_file = sys.argv[6]

out_fastq = f"{outdir}/{sampleid}.align.r2.fastq.gz"
out_barcode = f"{outdir}/{sampleid}.align.barcode.csv.gz"
out_read_clone = f"{outdir}/{sampleid}.ReadsCDR3.csv"

raw_mixcr_align = pd.read_csv(mixcr_align,sep = '\t')
clone_info = pd.read_csv(clone_info,sep = '\t')
barcode = pd.read_csv(barcode_file,compression='gzip',header=None,names = ['rname','barcode'])
del raw_mixcr_align['readId']
mixcr_align = raw_mixcr_align[raw_mixcr_align['cloneId'] != -1].copy()
subset_barcode = barcode[barcode['rname'].isin(mixcr_align['descrsR1'])].copy()

# write
align_dict = dict(zip(subset_barcode.index,[1] * subset_barcode.shape[0]))
if fastq_need_subset.endswith('.gz'):
    f = gzip.open(fastq_need_subset, 'rt')
else:
    f = open(fastq_need_subset, 'r')

with gzip.open(out_fastq,'wb') as fq_f:
    for idx, l in enumerate(f):
        
        if align_dict.get(idx//4,False):
            if idx % 4 == 0:
                rnames = l.split()[0][1:]
                fq_f.write(f'@{rnames}\n'.encode('utf-8'))
            if idx % 4 == 2:
                fq_f.write('+\n'.encode('utf-8'))
            if idx % 4 == 1:
                seqs = l.rstrip()
                fq_f.write(f'{seqs}\n'.encode('utf-8'))
            if idx % 4 == 3:
                qs = l.rstrip()
                fq_f.write(f'{qs}\n'.encode('utf-8'))
                
subset_barcode.to_csv(out_barcode,compression = 'gzip',header=None,index=None)
cloneseq_dict = dict(zip(clone_info['cloneId'],clone_info['nSeqCDR3']))
cloneaa_dict = dict(zip(clone_info['cloneId'],clone_info['aaSeqCDR3']))
total_read = barcode.shape[0]
align_read = raw_mixcr_align.shape[0]
clone_read = subset_barcode.shape[0]
mixcr_align['CDR3_aa'] = mixcr_align['cloneId'].map(cloneaa_dict)
mixcr_align['CDR3_n'] = mixcr_align['cloneId'].map(cloneseq_dict)
with open(f'{outdir}/step2.Mixcr.report','w') as f:
    f.write(f'## MIXCR Reads ##\nTotalReads\t{total_read}\t100.00%\nAlignReads\t{align_read}\t{align_read/total_read*100:.2f}%\nCloneReads\t{clone_read}\t{clone_read/total_read*100:.2f}% <<-\n\n')
mixcr_align.to_csv(out_read_clone,index=None)
