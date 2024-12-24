import argparse
parser = argparse.ArgumentParser(description="Barcode Index")

parser.add_argument("--sample", type=str, required=True, help="Name of the sample")
parser.add_argument("--bwl", type=str, required=True, help="Path to the barcode white list file")
parser.add_argument("--kmer", type=int, default = 5, help="K-mer size")
parser.add_argument("--threads", type=int, default = 20, help="Number of threads to use")
parser.add_argument("--pkl", type=str, required=True, help="Path to output pickle file")
args = parser.parse_args()

sample_name = args.sample
barcode_white_list_file = args.bwl
pkl = args.pkl
kmer = args.kmer
threads = args.threads

import datasketch as ds
import pyfastx 
import edlib
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocess import Pool
from datetime import datetime
import argparse
import pickle
import gzip


def time():
    now = datetime.now()
    formatted_now = now.strftime("%Y-%m-%d %H:%M:%S")
    return formatted_now

class BarcodeMatcher:
    
    def __init__(self, barcode_white_list_file, kmer = 5,threads = 20):
        print(">>>>>>",time())
        self.WL = barcode_white_list_file
        self.LSH_kmer = kmer
        self.threads = threads
        self.barcode_white_list = self.load_info()
        self.lsh_forest = ds.MinHashLSHForest(num_perm = 128)
        self.bulit_LSH_index()
        print(">>>>>>",time())
        
    def time(self):
        now = datetime.now()
        formatted_now = now.strftime("%Y-%m-%d %H:%M:%S")
        return formatted_now
        
    def load_info(self):

        if self.WL.endswith('.gz'):
            compression_option = 'gzip'
        else:
            compression_option = None
        barcodes_info = pd.read_csv(self.WL,header=None,names=['bam_x','bam_y','barcode'],compression=compression_option)
        return barcodes_info['barcode']
    
    def bulit_LSH_index(self):
        barcode_white_list = self.barcode_white_list
        minhashes = []
        print('Buliding Minhash Vector...')
        with Pool(self.threads) as p:
            for m in tqdm(p.imap(_get_minhash,barcode_white_list,chunksize = 10000),total=len(barcode_white_list)):
                minhashes.append(m)
        print('Buliding LSH Forest...')
        for i, minhash in tqdm(enumerate(minhashes),total=len(minhashes)):
            self.lsh_forest.add(i, minhash)
        self.lsh_forest.index()

def _get_minhash(seq):
    m = ds.MinHash(num_perm = 128)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        m.update(kmer.encode('utf8'))
    return m


def main(barcode_white_list_file, sample_name, kmer, threads, pkl):
    global k,barcode_white_list,lsh_forest
    k = kmer
    barcodematcher = BarcodeMatcher(barcode_white_list_file,kmer = kmer,threads = threads)
    barcode_white_list = barcodematcher.barcode_white_list
    lsh_forest = barcodematcher.lsh_forest
    data = pickle.dumps(lsh_forest)
    if pkl == None:
        with open(f'{sample_name}.barcode.pkl', 'wb') as file:
            file.write(data)
    else:
        with open(pkl, 'wb') as file:
            file.write(data)

main(barcode_white_list_file,sample_name, kmer = kmer,threads = threads,pkl = pkl)

