######
# sample_proxy = "7-ZXH_LC"
# mixcr_align_file = '/storage/fengyu/liuzhong/Mixcr_align/PE150_data_analyze/7.01/7-mix-V350247322/mixcr/L01/7-mix-V350247322.align.tsv'  #MIXCR-NGS-align
# mixcr_clone_file = "/storage/liuyi/09.ma_tcr/tmp/7-mix-V350247322.contigs.tsv"  #MIXCR-NGS-contig
# barcodemap_fastq_file = "/storage/liuyi/09.ma_tcr/data_precoess/p7/erdai/V350247322_L01_read_2.barcode.fq.gz"   #NGS-barcode_map-fastq
# cyclone_result_dir = "./cyclone_result"
# cyclone_match_dir = "./cyclone_match"
# cyclone_isotype = ["IGHA",'IGHG','IGHM','TRAC','TRBC','IGK','IGL']
# outmeta_file = "out.meta.gz"
######
import argparse
import sys


parser = argparse.ArgumentParser(description="Build meta")

parser.add_argument('--sample_proxy', type=str, required=True, 
                    help="Sample proxy same in cyclone analysis")
parser.add_argument('--mixcr_align_file', type=str, required=True, 
                    help="MIXCR-NGS-align file path")
parser.add_argument('--mixcr_clone_file', type=str, required=True, 
                    help="MIXCR-NGS-contig file path")
parser.add_argument('--barcodemap_fastq_file', type=str, required=True, 
                    help="NGS-barcode_map-fastq file path")
parser.add_argument('--cyclone_result_dir', type=str, required=True, 
                    help="Directory to cyclone result")
parser.add_argument('--cyclone_match_dir', type=str, required=True, 
                    help="Directory to cyclone barcode match results")
parser.add_argument('--outmeta_file', type=str, required=True, 
                    help="Output meta file")
parser.add_argument('--cyclone_isotype', type=str, nargs='+', default=["IGHA", "IGHG", "IGHM", "TRAC", "TRBC", "IGK", "IGL"], 
                    help="List of cyclone isotypes (default: ['IGHA','IGHG','IGHM','TRAC','TRBC','IGK','IGL'])")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
sample_proxy = args.sample_proxy
mixcr_align_file = args.mixcr_align_file
mixcr_clone_file = args.mixcr_clone_file
barcodemap_fastq_file = args.barcodemap_fastq_file
cyclone_result_dir = args.cyclone_result_dir
cyclone_match_dir = args.cyclone_match_dir
outmeta_file = args.outmeta_file
cyclone_isotype = args.cyclone_isotype
print(cyclone_isotype)

tmp_meta_file = outmeta_file[:-3] + '.raw.gz'

import gzip
import pandas as pd
import numpy as np

def read_fastq(fastq):
    rnames = list()
    seqs = list()
    if fastq.endswith('.gz'):
        f = gzip.open(fastq, 'rt')
    else:
        f = open(fastq, 'r')
    for idx, l in enumerate(f):
        if idx % 4 == 0:
            rnames.append(l.split()[0][1:])
        if idx % 4 == 1:
            seqs.append(l.rstrip())
    return rnames, seqs

def align_iterator(file_path):

    with open(file_path, 'r') as file:
        headers = next(file).strip().split('\t')
        for line in file:
            values = line.strip().split('\t')
            yield dict(zip(headers, values))


r2_rnames,r2_seqs = read_fastq(barcodemap_fastq_file)
locdict = {}
# ciddict = {}
# umidict = {}
cid_umidict = {}

for idx,rname in enumerate(r2_rnames):
    readid = rname.split('|||')[0]
    pos = rname.split('|||')[1].split(':')[2]
    locdict[readid] = pos
    cidwithumi = r2_seqs[idx]
    cid_umidict[readid] = cidwithumi
    
    # ciddict[readid] = cidwithumi[:25]
    # umidict[readid] = cidwithumi[25:]

import gc
del r2_rnames
del r2_seqs
gc.collect()

clone = pd.read_csv(mixcr_clone_file,sep = '\t')
topchaindict = dict(zip(clone['cloneId'],clone['topChains']))
cloneaadict = dict(zip(clone['cloneId'],clone['aaSeqCDR3']))
clonentdict = dict(zip(clone['cloneId'],clone['nSeqCDR3']))
clonev = dict(zip(clone['cloneId'],clone['allVHitsWithScore']))
cloned = dict(zip(clone['cloneId'],clone['allDHitsWithScore']))
clonej = dict(zip(clone['cloneId'],clone['allJHitsWithScore']))
clonec = dict(zip(clone['cloneId'],clone['allCHitsWithScore']))
isotypePrimary = dict(zip(clone['cloneId'],clone['isotypePrimary'].map({'IgA':'IGHA','IgG':'IGHG','IgM':'IGHM','IgD':'IGHD','IgE':'IGHE'})))

            
with gzip.open(tmp_meta_file, 'wb') as file:
    for i in align_iterator(mixcr_align_file):
        cloneId = int(i['cloneId'])
        readname = i['descrsR1'].split('/')[0]
        loc = locdict.get(readname,False)
        
        if cloneId != -1 and loc:
            AlignVHits = i['allVHitsWithScore']
            AlignDHits = i['allDHitsWithScore']
            AlignJHits = i['allJHitsWithScore']
            AlignCHits = i['allCHitsWithScore']
            incloneID = 'SR_' + str(cloneId)
            cdr3aa = cloneaadict[cloneId]
            cdr3nt = clonentdict[cloneId]
            isotype = topchaindict[cloneId]
            cid_umi_seq = cid_umidict[readname]
            
            cid = cid_umi_seq[:25]
            umi = cid_umi_seq[25:]
            
            # cid = ciddict[readname]
            # umi = umidict[readname]
            x,y = loc.split('_')
            
            priisotype = isotypePrimary[cloneId]
            
            if isinstance(priisotype,float):
                priisotype = isotype
            
            line = f'{readname}\t{AlignVHits}\t{AlignDHits}\t{AlignJHits}\t{AlignCHits}\t{incloneID}\t{cdr3aa}\t{cdr3nt}\t{isotype}\t{priisotype}\t{cid}\t{umi}\t{x}\t{y}\tSR\n'
            file.write(line.encode('utf-8'))
    print('erdai OK')
    # cyclone
    for priisotype in cyclone_isotype:
        tmpdict = {"IGHA":'IGHA','IGHG':'IGHG','IGHM':'IGHM','TRAC':'TRA','TRBC':'TRB','IGK':'IGK','IGL':'IGL'}
        cyclone_result = pd.read_csv(f'{cyclone_result_dir}/{sample_proxy}_{priisotype}.tsv',sep = '\t')
        match = pd.read_csv(f'{cyclone_match_dir}/{sample_proxy}_{priisotype}.match.csv')
        
        isotype = priisotype[:3]
        
        priisotype = tmpdict[priisotype]
        
        
        matchdict = dict(zip(match['rname'],match['mask']))
        for idx,row in cyclone_result.iterrows():
            # try:
            read_name = row['descrsR1']
            AlignVHits = row['allVHitsWithScore']
            AlignDHits = row['allDHitsWithScore']
            AlignJHits = row['allJHitsWithScore']
            AlignCHits = row['allCHitsWithScore']
            cloneId = row['cloneId']
            incloneId = 'LR_'+isotype+"_"+ str(cloneId)
            cdr3aa = row['CDR3_aa']
            cdr3nt = row['CDR3_n']
            if matchdict.get(read_name,False):
                cid = matchdict[read_name]
                umi = read_name.split('|||UB:')[1]
                x = row['loc'].split('_')[0]
                y = row['loc'].split('_')[1]

                line = f'{read_name}\t{AlignVHits}\t{AlignDHits}\t{AlignJHits}\t{AlignCHits}\t{incloneId}\t{cdr3aa}\t{cdr3nt}\t{isotype}\t{priisotype}\t{cid}\t{umi}\t{x}\t{y}\tLR\n'
                file.write(line.encode('utf-8'))
            # except:
            #     pass
        print(f'>>> {priisotype} done')


# meta2adata
meta = pd.read_csv(tmp_meta_file,sep = '\t',compression='gzip',header=None,names=['readname','AlignVHits','AlignDHits','AlignJHits','AlignCHits','cloneID','cdr3aa','cdr3nt','isotype','priisotype','cid','umi','x','y','readtype'])
meta['loc'] = meta['x'].map(str)+'_'+meta['y'].map(str)
meta['func'] = meta['cdr3aa'].map(lambda x: 'None'  if '_' in x or '*' in x  else 'Func')
meta['cdr3name'] = meta['cdr3aa']+'@'+meta['isotype']
clone = pd.read_csv(mixcr_clone_file,sep = '\t')
clone['cdr3name'] = clone['aaSeqCDR3'] + '@' + clone['topChains']
erset = set(clone['cdr3name'])
sanset = set(meta[meta['readtype'] == 'LR']['cdr3name'])
overlapmeta = meta[meta['cdr3name'].isin(erset)].copy()



check_list = [ x for x in cyclone_isotype if x in ['IGHA', 'IGHG', 'IGHM']]
demeta = overlapmeta.drop_duplicates(['cdr3aa','priisotype','loc','readtype'])
ighdf = demeta[demeta['isotype'] == 'IGH']


if len(check_list) != 0 and ighdf.shape[0] != 0:
    addsr = []
    for loc,tmp in ighdf.groupby('loc'):
        readtypes = tmp['readtype'].unique().tolist()
        if len(readtypes) == 2:
            sr = tmp[(tmp['readtype'] == 'SR') & (tmp['priisotype'] == 'IGH')].copy()
            lr = tmp[tmp['readtype'] == 'LR']
            if len(set(sr['cdr3aa']).intersection(set(lr['cdr3aa']))) > 0:
                lrdict = dict(zip(lr['cdr3aa'],lr['priisotype']))
                sr['addisotype'] = sr['cdr3aa'].map(lrdict)
                sr = sr[~sr['addisotype'].isna()]
                if sr.shape[0] != 0:
                    addsr.append(sr)
    addsr = pd.concat(addsr)
    adddict = dict(zip(addsr.index,addsr['addisotype']))
    overlapmeta['addinfo'] = overlapmeta.index.map(adddict)
    overlapmeta['CombineIsotype'] = overlapmeta.apply(lambda row: row['addinfo'] if pd.notna(row['addinfo']) else row['priisotype'], axis=1)
    overlapmeta['cdr3name'] = overlapmeta['cdr3aa'] + '@' + overlapmeta['CombineIsotype']
else:
    overlapmeta['addinfo'] = np.nan
    overlapmeta['CombineIsotype'] = overlapmeta['isotype']
    overlapmeta['cdr3name'] = overlapmeta['cdr3aa'] + '@' + overlapmeta['CombineIsotype']
    
if outmeta_file.endswith('gz'):
    overlapmeta.to_csv(outmeta_file,index=None,compression='gzip')
else:
    overlapmeta.to_csv(outmeta_file,index=None)

