import pandas as pd
import numpy as np
from collections import Counter
from tqdm import tqdm
import matplotlib.pyplot as plt
import re
import sys
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
import Levenshtein

def cluster_strings(series, max_distance=2):
    n = len(series)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dist = Levenshtein.distance(series[i], series[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    Z = linkage(squareform(dist_matrix), 'complete')
    clusters = fcluster(Z, max_distance, criterion='distance')
    tmp = pd.DataFrame({'umi':series,'cluster':clusters}).sort_values('cluster')
    return dict(zip(tmp['umi'],tmp['cluster']))


def get_Hits(x):
    if isinstance(x,float):
        return 'NA'
    return re.sub(r'\(.*?\)', '', x)

def find_kmers(input_list, k):
    kmers = []
    for sequence in input_list:
        sequence = sequence.replace('_', '').replace('*', '')
        for i in range(0, len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmers.append(kmer)
    kmer_counts = Counter(kmers)
    return kmer_counts

def score_sequences(input_list, k):
    kmer_counts = find_kmers(input_list,k = k)
    scores = []
    for sequence in input_list:
        clean_sequence = sequence.replace('_', '').replace('*', '')
        sequence_score = 0
        for i in range(len(clean_sequence) - k + 1):
            kmer = clean_sequence[i:i+k]
            sequence_score += kmer_counts[kmer]
        scores.append([sequence,sequence_score])
    return pd.DataFrame(scores,columns=['name','score']).sort_values('score',ascending=False).set_index('name')

def get_valid_cdr3(subdf,target = 'CDR3_aa',k = 5):
    cdr3aa_df = score_sequences(subdf[target].unique(),k = k)
    cdr3aa_df['count'] = subdf[target].value_counts()
    cdr3aa_df['total_count'] = cdr3aa_df['score'] * cdr3aa_df['count']
    cdf3aa = cdr3aa_df['total_count'].sort_values(ascending=False).index[0]
    return cdf3aa


match_file = sys.argv[1]
align_file = sys.argv[2]
result_profix = sys.argv[3]

match = pd.read_csv(match_file)
align = pd.read_csv(align_file)

match['loc'] = match['bam_x'].map(str) +'_' + match['bam_y'].map(str)
locdict = dict(zip(match['rname'],match['loc']))
align['loc'] = align['descrsR1'].map(locdict)
cdr3 = align[~align['loc'].isna()].reset_index(drop=True)
cdr3 = cdr3[~cdr3['CDR3_n'].isna()].copy()
cdr3['umi'] = cdr3['descrsR1'].map(lambda x : x.split('|||')[1].split('UB:')[1])
cdr3 = cdr3[cdr3['CDR3_aa'].map(len).apply(lambda x  :  x > 5 and x < 30)].copy()
cdr3['func'] = cdr3['CDR3_aa'].apply(lambda x : 'nofunc' if '_' in x or '*' in x  else 'func')

vc = cdr3['loc'].value_counts()
weak_evidence = cdr3[cdr3['loc'].isin(vc[vc == 1].index)].copy()
pass_loc = vc[vc > 1].index
pass_filter_v1 = cdr3[cdr3['loc'].isin(vc[vc > 1].index)].copy()

m = 0
result_list = []
for idx,tmp in tqdm(pass_filter_v1.groupby('loc'),total = len(pass_loc) ):
    umi_list = tmp['umi'].drop_duplicates().tolist()
    if len(umi_list) != 1:
        umi_dict = cluster_strings(umi_list,max_distance = 2)
        tmp['cluster_umi'] = tmp['umi'].map(umi_dict)
    else:
        tmp['cluster_umi'] = 1
    result_list.append(tmp)
umi_cdr3 = pd.concat(result_list)
umi_cdr3['stat'] = umi_cdr3['loc'] + '@' + umi_cdr3['cluster_umi'].map(str)
for typ in ['V','D','J','C']:
    umi_cdr3[f'{typ}Hits'] = umi_cdr3[f'all{typ}HitsWithScore'].apply(get_Hits)

result_list = []
m = 0
for idx,tmp in tqdm(umi_cdr3.groupby('loc'),total = umi_cdr3['loc'].drop_duplicates().shape[0]):
    tmp_dict = {}
    subdf_list = []
    for umi,subdf in tmp.groupby('cluster_umi'):
        if subdf.shape[0] == 1:
            subdf['umiCDR3'] = subdf['CDR3_aa']
            subdf['umistate'] = 'single'
        else:
            subdf['vj'] = subdf['VHits']+'@'+subdf['JHits']
            vc = subdf['vj'].value_counts()
            if vc[vc>1].shape[0] == 0:
                subdf['umiCDR3']= subdf['CDR3_aa']
                subdf['umistate'] = 'weakevi'
            else:
                cdf3aa = get_valid_cdr3(subdf,target = 'CDR3_aa', k =5)
                subdf['umiCDR3'] = cdf3aa
                subdf['umistate'] = 'pass'
        subdf_list.append(subdf)
    loc_umi_modify = pd.concat(subdf_list)
    result_list.append(loc_umi_modify)

modify_umi_cdf3 = pd.concat(result_list)
modify_umi_cdf3.to_csv(f'{result_profix}.tsv',sep = '\t',index=None)
weak_evidence.to_csv(f'{result_profix}.weakevidence.tsv',sep = '\t',index=None)
