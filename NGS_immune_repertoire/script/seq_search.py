import sys
from multiprocessing import Pool
import numpy as np
import pandas as pd
import gzip
from tqdm import tqdm
import gzip
import edlib
import os

fastq = sys.argv[1]
filename = sys.argv[2]
sequence = sys.argv[3]
logname = os.path.basename(filename)[:-3]+'.log'



def replace_ranges_with_X(input_str, ranges):
    result_str = list(input_str)
    for start, end in ranges:
        for i in range(start, end + 1):
            result_str[i] = 'X'
    return ''.join(result_str)

def rev_seq(s):  
        rev_compl_l = [chr(i) for i in range(128)]
        rev_compl_l[ord('A')] = 'T'
        rev_compl_l[ord('C')] = 'G'
        rev_compl_l[ord('G')] = 'C'
        rev_compl_l[ord('T')] = 'A'
        return ''.join(rev_compl_l[ord(c)] for c in reversed(s))
    
    
class new_finder:
    
    def __init__(self, fastq,threads=10):
        self.fastq = fastq
        self.threads = threads
        basename = os.path.basename(fastq)
    
        if 'fastq' in basename or 'fq' in basename:
            self._rnames,self._seqs,self._qs = self.read_fastq(fastq)
        else:
            self._rnames,self._seqs = self.read_fasta(fastq)
        self.aln = []

    def alns(self,df=True):
        if df:
            return pd.concat([pd.DataFrame(self._rnames,columns = ['readID']),pd.DataFrame(self.aln,columns = ['strand','distance','location','start'])],axis = 1)
        return self._rnames,self.aln
    
    def read_fasta(self, fasta):
        rnames = []
        seqs = []

        if fasta.endswith('.gz'):
            f = gzip.open(fasta, 'rt')
        else:
            f = open(fasta, 'r')

        for idx, line in enumerate(f):
            if idx % 2 == 0:  
                rnames.append(line.strip()[1:])
            else:           
                seqs.append(line.strip())

        f.close()  
        return rnames, seqs

    def read_fastq(self,fastq):
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
    
    def rev_compl(self,s):  
        rev_compl_l = [chr(i) for i in range(128)]
        rev_compl_l[ord('A')] = 'T'
        rev_compl_l[ord('C')] = 'G'
        rev_compl_l[ord('G')] = 'C'
        rev_compl_l[ord('T')] = 'A'
        return ''.join(rev_compl_l[ord(c)] for c in reversed(s))
    
    
    def get_query_alns(self,query,max_distance = 5):
        seqs = self._seqs
        rev_query = self.rev_compl(query)
        q_list = [query] * len(seqs)
        rq_list = [rev_query] * len(seqs)
        max_distance = [max_distance] *len(seqs)
        tmp_result = []
        with Pool(self.threads) as p:
            for single_pos_result,single_neg_result in tqdm(p.imap(self.search_all,zip(seqs,q_list,rq_list,max_distance),chunksize=10000),total = len(seqs)):
                tmp_result.append([single_pos_result,single_neg_result])
        self.alns_result = tmp_result
        return tmp_result
    
    def get_start(self,with_distance = False):
        num_stat = []
        result_pos = []
        
        for alns in self.alns_result:
            pos,neg = alns
            pos_dict,pos_num = self.dedup_start(pos)
            neg_dict,neg_num = self.dedup_start(neg)
            num_stat.append([pos_num,neg_num])
            
            
            if with_distance:
                result_pos.append([pos_dict,neg_dict])
            else:
                pos_list,neg_list = 'NA','NA'
                if len(pos_dict) != 0:
                    pos_list = []
                    for i in pos_dict.values():
                        pos_list.extend(i)
                if len(neg_dict) != 0:
                    neg_list = []
                    for i in neg_dict.values():
                        neg_list.extend(i)
                result_pos.append([pos_list,neg_list])
            
            
        num_stat_df = pd.DataFrame(num_stat,columns=['+','-'])
        num_stat_df['sum'] = num_stat_df['+'] + num_stat_df['-']
        return result_pos,num_stat_df
    
    @staticmethod
    def dedup_start(aln_dict):
        cadidate_dict = {}
        match_num = 0

        for distance,candidate_loci in aln_dict.items():
            start_list = []
            for start,end in candidate_loci:
                start_list.append(start)
            start_list = list(set(start_list))
            cadidate_dict[distance] = start_list
            match_num += len(start_list)
        return cadidate_dict,match_num
    
    @staticmethod
    def search_all(in_zip):
        raw_seq,query,r_query,max_distance = in_zip
        single_pos_result = {}
        search_type = True
        seq = raw_seq
        while search_type:
            search_result = edlib.align(query, seq, 'HW', 'locations', k = max_distance)
            if search_result['editDistance'] == -1:
                search_type = False
            else:
                seq = replace_ranges_with_X(seq,search_result['locations'])
                single_pos_result[search_result['editDistance']] = [(x[0], x[1] + 1) for x in search_result['locations']]

        single_neg_result = {}
        seq = raw_seq
        search_type = True
        while search_type:
            search_result = edlib.align(r_query, seq, 'HW', 'locations', k = max_distance)
            if search_result['editDistance'] == -1:
                search_type = False
            else:
                seq = replace_ranges_with_X(seq,search_result['locations'])
                single_neg_result[search_result['editDistance']] = [(x[0], x[1] + 1) for x in search_result['locations']]

        return single_pos_result,single_neg_result
    
    @staticmethod
    def add(x,num):
        return x+num

def getmindict(indict):
    keys = indict.keys()
    minkey = min(keys)
    return indict[minkey]


fastq_finder = new_finder(fastq,threads = 10)
r1 = sequence
r1_alns_result = fastq_finder.get_query_alns(r1,max_distance = 1)
r1_result_pos,r1_num_stat = fastq_finder.get_start(with_distance=True)

duetoambi = 0
fix_valid = 0
barcode_valid = 0
total = 0
with gzip.open(filename, 'wt') as file:
    for idx,(f,r) in tqdm(enumerate(r1_result_pos),total = len(r1_result_pos)):
        total += 1
        fpos = 'NA'
        rpos = 'NA'
        if len(f) !=0:        # forward    umi_fix_barcode
            if len(f) >1 :
                fpos = getmindict(f)[0]
            else:
                fpos = list(f.values())[0][0]

            fix_start = fpos
            fix_end = fpos + 15
            barcode_start = fix_end
            barcode_end = barcode_start + 25
            strand = '+'

        if len(r) !=0:
            if len(r) >1:
                rpos = getmindict(r)[0]
            else:
                rpos = list(r.values())[0][0]
            fix_start = rpos
            fix_end = rpos + 15
            
            barcode_end = fix_start
            barcode_start = barcode_end -25
            strand = '-'

        if rpos!='NA' and fpos !='NA':
            duetoambi += 1
            pass
        elif rpos!='NA' or fpos!='NA':
            seq = fastq_finder._seqs[idx]
            barcode_seq = seq[barcode_start:barcode_end]
            fix_seq = seq[fix_start:fix_end]
            if strand == '-':
                barcode_seq = rev_seq(barcode_seq)
                fix_seq = rev_seq(fix_seq)
            total_seq = barcode_seq 
            fix_valid += 1
            if len(total_seq) == 25:
                read_name = fastq_finder._rnames[idx]
                file.write(f"@{read_name}\n")
                file.write(rev_seq(barcode_seq) + "\n")  # reverse
                file.write("+\n")
                file.write("I" * len(barcode_seq) + "\n")
                barcode_valid +=1
        else:
            pass
        
with open(logname, 'w') as file:
    file.write(f'Total\t{total}\t\n')
    file.write(f'Duetoambi\t{duetoambi}\t{duetoambi/total*100:.2f}%\n')
    file.write(f'Fix_valid\t{fix_valid}\t{fix_valid/total*100:.2f}%\n')
    file.write(f'Barcode_valid\t{barcode_valid}\t{barcode_valid/total*100:.2f}%\n')
