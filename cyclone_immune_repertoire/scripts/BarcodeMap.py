import sys

sample_name = sys.argv[1]
r1_bf = sys.argv[2]
bwl = sys.argv[3]#'./erdai_list.csv.gz'
pkl = sys.argv[4] #'./erdai.pkl'
outdir = sys.argv[5]
threads = int(sys.argv[6])

global ann,k
ann = int(sys.argv[7])
k = int(sys.argv[8])
chunk_size = int(sys.argv[9])


result_file = f'{outdir}/{sample_name}.match.csv'
report_file = f'{outdir}/step3.BarcodeMap.report'


import datasketch as ds
import edlib
import pandas as pd
import numpy as np  
from tqdm import tqdm
from multiprocess import Pool
import multiprocessing
from datetime import datetime
import pickle
import gc

def time():
    now = datetime.now()
    formatted_now = now.strftime("%Y-%m-%d %H:%M:%S")
    return formatted_now


def read_csv_or_gzip(file_path):
    if file_path.endswith('.gz'):
        df = pd.read_csv(file_path, header=None,compression='gzip')
    else:
        df = pd.read_csv(file_path, header=None)
    return df

def load_info(r1_bf,bwl):
    barcode_candidate = read_csv_or_gzip(r1_bf)
    candidate_list = barcode_candidate.values.tolist()
    barcodes_info = pd.read_csv(bwl,header=None,names=['bam_x','bam_y','barcode'])
    barcode_dict = {idx: value for idx, value in enumerate(barcodes_info['barcode'])}
    return candidate_list,barcode_dict

def _lsh_search(query_list):
    query_name = query_list[0]
    query_sequence = query_list[1]    
    query_minhash = ds.MinHash(num_perm = 128)
    for i in range(len(query_sequence) - k + 1):
        kmer = query_sequence[i:i + k]
        query_minhash.update(kmer.encode('utf8'))
    query_results = lsh_forest.query(query_minhash,ann)
    return [query_name,query_sequence,query_results]

def _lw_match(match_tmp):
    rname,seq,read_idx = match_tmp
    min_distance = 10
    final_match_seq = "NA"
    for idx in read_idx:
        result = barcode_dict.get(idx,None)
        e_r = edlib.align(result, seq, 'HW', 'locations',k = 5)
        if e_r['editDistance'] < 5 and e_r['editDistance'] != -1:
            if  e_r['editDistance'] < min_distance:
                min_distance = e_r['editDistance']
                final_match_seq = result
    if min_distance != 10:
        return [rname,seq,final_match_seq,min_distance]
    else:
        return 'NA'

from Levenshtein import distance
def lw_match_25bp(match_tmp):
    rname,seq,read_idx = match_tmp
    min_distance = 100
    mapping_list = []
    for idx in read_idx:
        mask_barcode = barcode_dict.get(idx,None)
        dis = distance(mask_barcode,seq)
        if dis < min_distance:
            mapping_list = []
            mapping_list.append(mask_barcode)
            min_distance = dis
        elif dis == min_distance:
            mapping_list.append(mask_barcode)
    return [rname,seq,mapping_list,min_distance]
    
with open(pkl, 'rb') as file:
    lsh_forest = pickle.load(file)


candidate_list,barcode_dict = load_info(r1_bf,bwl)
sub_lists = [candidate_list[i:i + chunk_size] for i in range(0, len(candidate_list), chunk_size)]
match_file_list = []
for idx,sub_list in enumerate(sub_lists):
    print(f'>> Round {idx}/{len(sub_lists)}')
    match_list = []
    with tqdm(total = len(sub_list)) as pbar:
        def update(v):
            match_list.append(v)
            pbar.update()
        pool = multiprocessing.Pool(threads)
        for i in sub_list:
            pool.apply_async(_lsh_search,(i,),callback = update)
        pool.close()
        pool.join()
    
    match_file_list.append(f'{outdir}/{sample_name}.{idx}.pkl')
    with open(f'{outdir}/{sample_name}.{idx}.pkl', 'wb') as file:
        pickle.dump(match_list, file)


lw_match_list = []
for match_file in match_file_list:
    with open(match_file, 'rb') as file:
        match_list = pickle.load(file)
    
    # with tqdm(total = len(match_list)) as pbar:
    #     def update(v):
    #         lw_match_list.append(v)
    #         pbar.update()
    #     pool = multiprocessing.Pool(threads)
    #     for i in match_list:
    #         pool.apply_async(_lw_match,(i,),callback = update)
    #     pool.close()
    #     pool.join()
    with tqdm(total = len(match_list)) as pbar:
        def update(v):
            rname,seq,mapping_list,min_distance = v
            if len(mapping_list) == 1 and min_distance <= 4:
                lw_match_list.append([rname,seq,mapping_list[0],min_distance])
            pbar.update()
        pool = multiprocessing.Pool(threads)
        for i in match_list:
            pool.apply_async(lw_match_25bp,(i,),callback = update)
        pool.close()
        pool.join()


# result_list = [i for i in lw_match_list if i != 'NA']
result_list = lw_match_list
result = pd.DataFrame(result_list,columns = ['rname','split','mask','dis'])

barcodes_info = pd.read_csv(bwl,header = None,names = ['bam_x','bam_y','mask'])
result = result.merge(barcodes_info,how = 'left', on ='mask')
result.to_csv(result_file,index = None)

with open(report_file,'w') as f:
    total_counts = len(candidate_list)
    match_reads = result.shape[0]
    f.write(f"Total Count\t{total_counts}\n## Barcode Map ##\nMapped Reads\t{match_reads}\t{match_reads/total_counts*100:.2f}%\n")
