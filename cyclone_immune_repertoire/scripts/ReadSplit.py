import sys

fastq = sys.argv[1]
outdir = sys.argv[2]
sampleId = sys.argv[3]
threads = int(sys.argv[4])

from multiprocessing import Pool
import numpy as np
import pandas as pd
import gzip
from tqdm import tqdm
import gzip
import edlib


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
        self._rnames,self._seqs,self._qs = self.read_fastq(fastq)
        self.aln = []

    def alns(self,df=True):
        if df:
            return pd.concat([pd.DataFrame(self._rnames,columns = ['readID']),pd.DataFrame(self.aln,columns = ['strand','distance','location','start'])],axis = 1)
        return self._rnames,self.aln

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
    
def get_state(row):
    state = 'NA'
    r1 = row['r1']
    r2 = row['r2']
    fix = row['fix']
    if r1+r2+fix == 0:
        state = 'MissALL'
    elif r1 == 0:
        if r2!=0 and fix !=0:
            state = 'MissR1'
        else:
            state = 'MissMut'
    elif r2 == 0 and r1!=0 and fix !=0:
        if r1!=0 and fix !=0:
            state = 'MissR2'
        else:
            state = 'MissMut'
    elif fix == 0 and r1!= 0 and r2 !=0:
        if r1!=0 and r2 !=0:
            state = 'MissFIX'
        else:
            state = 'MissMut'
    else:
        multiple =[]
        for idx,check in zip(['1','2','f'],[r1,r2,fix]):
            if check >= 3:
                multiple.append(idx)
        if len(multiple) == 3:
            state = 'MutALL'
        elif len(multiple) == 2:
            state = 'MutMut'
        elif len(multiple) == 1:
            match_dict = {'f':'MutFIX','2':'MutR2','1':'MutR1'}
            state = match_dict[multiple[0]]
        else:
            state = 'Simple'
    return state

def find_closest_triplet(a_list, b_list, c_list):
    for tmp in [a_list,b_list,c_list]:
        if tmp == 'NA' or len(tmp) == 0:
            return 'NA'
    a_list.sort(reverse=True)
    b_list.sort(reverse=True)
    c_list.sort(reverse=True)
    closest_a = a_list[0]
    closest_b = b_list[0]
    closest_c = c_list[0]
    closest_diff = float("inf")
    for a in a_list:
        for b in b_list:
            if a <= b:
                continue
            for c in c_list:
                if b <= c:
                    continue
                diff = a - c
                if diff < closest_diff:
                    closest_a, closest_b, closest_c = a, b, c
                    closest_diff = diff
    if closest_a > closest_b and closest_b > closest_c:
        return [closest_a, closest_b, closest_c]
    else:
        return 'NA'

def find_all_closest_triplet(a,b,c,reverse = False):
    result_list = []
    tmp_sign = True
    while tmp_sign:
        position = find_closest_triplet(a, b, c)
        if position != 'NA':
            a = [num for num in a if num != position[0]]
            b = [num for num in b if num != position[1]]
            c = [num for num in c if num != position[2]]
            result_list.append(position)
        else:
            tmp_sign = False
    if reverse:
        new_list = []
        for result in result_list :
            new_list.append(result[::-1])
        return new_list
    return result_list

def find_shortest_tuple(num, tmp):
    candidate_tuples = []
    for key, value in tmp.items():
        for tup in value:
            if tup[0] == num:
                candidate_tuples.append(tup)
    if candidate_tuples:
        return min(candidate_tuples, key=lambda x: x[1] - x[0])
    else:
        return None 
    
def get_split_tuple(row):
    idx = row['idx']
    fix = row['fix']
    r1 = row['r1']
    r2 = row['r2']
    if fix >= r1:
        fix_tuple = find_shortest_tuple(fix,fix_alns_result[idx][0])
        r1_tuple = find_shortest_tuple(r1,r1_alns_result[idx][0])
        r2_tuple = find_shortest_tuple(r2,r2_alns_result[idx][1])
        
        umi_split_tuple = (fix_tuple[1],fix_tuple[1] + 10)
        cid_split_tuple = (r1_tuple[1],fix_tuple[0])
        r2_split_tuple = (0,r2_tuple[0])
        strand = '+'
    else:
        fix_tuple = find_shortest_tuple(fix,fix_alns_result[idx][1])
        r1_tuple = find_shortest_tuple(r1,r1_alns_result[idx][1])
        r2_tuple = find_shortest_tuple(r2,r2_alns_result[idx][0])
        
        umi_split_tuple = (fix_tuple[0]-10,fix_tuple[0])
        cid_split_tuple = (fix_tuple[1],r1_tuple[0])
        r2_split_tuple = (r2_tuple[1],0)
        strand = '-'
    return fix_tuple[1],r1_tuple[1],r2_tuple[1],umi_split_tuple,cid_split_tuple,r2_split_tuple,strand

def r2_distance(row):
    r2_range = row['r2_range']
    idx = row['idx']
    if r2_range[0] == 0:
        distance = r2_range[1]
    else:
        start = r2_range[0]
        total = len(fastq_finder._seqs[idx])
        distance = total - start 
    return distance


fastq_finder = new_finder(fastq,threads = threads)

r1 = 'ATGGCGACCTTATCAG'
r1_alns_result = fastq_finder.get_query_alns(r1,max_distance = 3)
r1_result_pos,r1_num_stat = fastq_finder.get_start(with_distance=False)

r2 = 'GCCATGTCGTTCTGTGAGCCAAGGAGTT'
r2_alns_result = fastq_finder.get_query_alns(r2,max_distance = 5)
r2_result_pos,r2_num_stat = fastq_finder.get_start(with_distance=False)

fix = 'TTGTCTTCCTAAGAC'
fix_alns_result = fastq_finder.get_query_alns(fix,max_distance = 3)
fix_result_pos,fix_num_stat = fastq_finder.get_start(with_distance=False)

# pattern
stat = pd.concat([r1_num_stat,r2_num_stat,fix_num_stat],axis = 1)
stat.columns = ['r1+','r1-','r1','r2+','r2-','r2','fix+','fix-','fix']
stat = stat[['r1','r2','fix']].copy()
stat['state'] = stat.apply(lambda x :get_state(x),axis = 1)
position_list = []
for m in tqdm(stat.index,total = stat.shape[0]):
    a,b = fix_result_pos[m]
    c,d = r1_result_pos[m]
    e,f = r2_result_pos[m]
    position_1 = find_all_closest_triplet(a, c, f)
    position_2 = find_all_closest_triplet(e, d, b,reverse = True)        
    position_list.append([position_1,position_2])

stat_list = []
oneD1_list = []
oneD2_list = []
for com1,com2 in position_list:
    com1_sign,com2_sign = False,False
    if len(com1) == 0:
        com1_sign = True
    if len(com2) == 0:
        com2_sign = True
    
    if com1_sign and com2_sign:
        stat_list.append('No')
    elif not com1_sign and not com2_sign:
        stat_list.append('Ambi')            
    else:
        if not com1_sign:
            if len(com1) > 1:
                stat_list.append('1dxm')
            else:
                stat_list.append('1d1m')
        elif not com2_sign:
            if len(com2) > 1:
                stat_list.append('1dxm')
            else:
                stat_list.append('1d1m')

# get valid
t = pd.DataFrame(stat_list)
match_df = t[t[0] == '1d1m']
pos_stat_list = []
for i in match_df.index:
    position = position_list[i]
    pos_stat_list.append([i]+[pattern for pattern in position if len(pattern) != 0][0][0])
pos_stat = pd.DataFrame(pos_stat_list,columns=['idx','fix','r1','r2'])
pos_stat[['fix_end','r1_end','r2_end','umi_range','cid_range','r2_range','strand']] = pos_stat.apply(lambda row : get_split_tuple(row),axis = 1,result_type='expand')
pos_stat['cid_lenth'] = pos_stat['cid_range'].apply(lambda x : x[1] - x[0])
pos_stat['r2_lenth'] = pos_stat.apply(lambda x :r2_distance(x),axis = 1)
valid_reads_df = pos_stat[(pos_stat['cid_lenth'] >= 20) & (pos_stat['cid_lenth'] <= 30) & (pos_stat['r2_lenth'] <= 10000)].copy()

# write fastq
with gzip.open(f'{outdir}/{sampleId}.barcode.csv.gz','wb') as barcode_file ,gzip.open(f'{outdir}/{sampleId}.r2.fastq.gz','wb') as r2_fastq_file:
    for idx,row in tqdm(valid_reads_df.iterrows(),total = valid_reads_df.shape[0]):


        fastq_idx = row['idx']
        umi_start,umi_end = row['umi_range']
        cid_start,cid_end = row['cid_range']
        r2_start,r2_end = row['r2_range']
        strand = row['strand']
        
        
        seq = fastq_finder._seqs[fastq_idx]
        q = fastq_finder._qs[fastq_idx]
        rname = fastq_finder._rnames[fastq_idx]
        
        cid = seq[cid_start:cid_end]
        umi = seq[umi_start:umi_end]
        
        if r2_start ==0:
            r2_seq = seq[:r2_end]
            r2_seq_q = q[:r2_end]
        else:
            r2_seq = seq[r2_start:]
            r2_seq_q = q[r2_start:]

        if strand == '-':
                cid = rev_seq(cid)
                umi = rev_seq(umi)
        rname += f'|||UB:{umi}'
        barcode_file.write(f'{rname},{cid}\n'.encode('utf-8'))
        r2_fastq_file.write(f'@{rname}\n{r2_seq}\n+\n{r2_seq_q}\n'.encode('utf-8'))


# stat and plot
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = 'Arial'
sns.set_style('white', {'axes.grid' : False})

plot_state = stat['state'].value_counts().reset_index()
plt.figure(figsize=(5, 4))
ax = sns.barplot(data=plot_state, x='index', y='state')
plt.xticks(rotation=45)

total = float(len(plot_state)) 
for idx,p in enumerate(ax.patches):
    height = p.get_height()
    percentage = '{:.1f}%'.format(100 *  plot_state['state'].iloc[idx]/plot_state['state'].sum())
    ax.text(p.get_x() + p.get_width()/2., height, percentage, ha="center", va='bottom', rotation=0)
plt.savefig(f'{outdir}/fix_search.pdf', bbox_inches='tight')

pos_stat['r1_fix'] = np.abs(pos_stat['r1'] - pos_stat['fix'])
pos_stat['r1_r2'] = np.abs(pos_stat['r1'] - pos_stat['r2'])
r1_fix_df = pos_stat[pos_stat['r1_fix'] < 100]

plt.figure(figsize=(4, 3))
plt.hist(r1_fix_df['r1_fix'], bins=50, color='grey', edgecolor='black')
total = pos_stat.shape[0]
filter_r12fix = pos_stat[pos_stat['r1_fix'] < 100].shape[0]
plt.title(f'Fix Start-R1 Start Lenth\nTotal:{total} reads\n<100bp:{filter_r12fix} {filter_r12fix/total *100:.2f}%')
plt.savefig(f'{outdir}/r1_to_fix.lenth.pdf', bbox_inches='tight')

r1_r2_df = pos_stat[pos_stat['r1_r2'] < 100]
plt.figure(figsize=(4, 3))
plt.hist(r1_r2_df['r1_r2'], bins=50, color='grey', edgecolor='black')
filter_r12r2 = pos_stat[pos_stat['r1_r2'] < 100].shape[0]
plt.title(f'R2 Lenth\nTotal:{total} reads\n<100bp:{filter_r12r2} {filter_r12r2/total *100:.2f}%')
plt.savefig(f'{outdir}/r1_to_r2.lenth.pdf', bbox_inches='tight')


cid_len_df = pos_stat['cid_lenth']
plt.figure(figsize=(4, 3))
plt.hist(cid_len_df[cid_len_df<50],bins = 50,color='grey', edgecolor='black')
total_cid_len_reads = cid_len_df.shape[0]
filter_cid_len_reads = cid_len_df[(cid_len_df<=30) & (cid_len_df>=20)].shape[0]
plt.title(f'CID Lenth\nTotal:{total_cid_len_reads} reads\nBarcode20~30bp:{filter_cid_len_reads} {filter_cid_len_reads/total_cid_len_reads *100:.2f}%')
plt.savefig(f'{outdir}/barcode.lenth.pdf', bbox_inches='tight')


r2_len_df = pos_stat['r2_lenth']
plt.figure(figsize=(4, 3))
plt.hist(r2_len_df[r2_len_df<1000],bins = 50,color='grey', edgecolor='black')
total_r2_len_reads = r2_len_df.shape[0]
filter_r2_len_reads = r2_len_df[r2_len_df < 1000].shape[0]
plt.title(f'R2 Lenth\nTotal:{total_r2_len_reads} reads\n<1000bp:{filter_r2_len_reads} {filter_r2_len_reads/total_r2_len_reads *100:.2f}%')
plt.savefig(f'{outdir}/r2.lenth.pdf', bbox_inches='tight')



from matplotlib.colors import ListedColormap
t_counts = t[0].value_counts()
set3_palette = sns.color_palette("Set3", len(t_counts))
set3_cmap = ListedColormap(set3_palette.as_hex())
plt.figure(figsize=(5, 3.5))
plt.pie(
    t_counts,
    autopct='%1.1f%%',
    startangle=90,
    colors=set3_cmap.colors, 
    wedgeprops={'edgecolor': 'black'}, 
    pctdistance=0.5
)
plt.axis('equal')
plt.legend(t_counts.index, bbox_to_anchor=(0.8, 1), loc='upper left')
plt.savefig(f'{outdir}/patten.pie.pdf', bbox_inches='tight')




with open(f'{outdir}/step1.Split.report','w') as f:
    total_counts = len(r1_alns_result)
    match_r1 = r1_num_stat[r1_num_stat['sum']>0].shape[0]
    match_r2 = r2_num_stat[r2_num_stat['sum']>0].shape[0]
    match_fix = fix_num_stat[fix_num_stat['sum']>0].shape[0]
    f.write(f'Total Count\t{total_counts}\n## Fixed Sequence Search ##\nR1\t{match_r1}\t{match_r1/total_counts*100:.2f}%\nR2\t{match_r2}\t{match_r2/total_counts*100:.2f}%\nFix\t{match_fix}\t{match_fix/total_counts*100:.2f}%\n\n')
    state_dict = plot_state.set_index('index').to_dict()['state']
    simple = state_dict.get('Simple',0)
    miss = state_dict.get('MissALL',0) + state_dict.get('MissMut',0) + state_dict.get('MissR1',0) + state_dict.get('MissR2',0) + state_dict.get('MissFIX',0)
    mul = state_dict.get('MutALL',0) + state_dict.get('MutMut',0) + state_dict.get('MutR1',0) + state_dict.get('MutR2',0) + state_dict.get('MutFIX',0)
    f.write(f'## Fixed Sequcence Match ##\nSimple\t{simple}\t{simple/total_counts*100:.2f}% <<-\nMiss\t{miss}\t{miss/total_counts*100:.2f}%\nMul\t{mul}\t{mul/total_counts*100:.2f}%\n\n')
    no = t_counts.get('No',0)
    d1m1 = t_counts.get('1d1m',0)
    d1mx = t_counts.get('1dxm',0)
    ambi = t_counts.get('Ambi',0)
    f.write(f"## Pattern Search ##\nNo\t{no}\t{no/total_counts*100:.2f}%\n1d1m\t{d1m1}\t{d1m1/total_counts*100:.2f}% <<-\n1dxm\t{d1mx}\t{d1mx/total_counts*100:.2f}%\nambi\t{ambi}\t{ambi/total_counts*100:.2f}%\n\n")
    valid = valid_reads_df.shape[0]
    f.write(f"## Lenth Filter ##\nValid\t{valid}\t{valid/total_counts*100:.2f}% <<-\n")

valid_reads_df.to_csv(f'{outdir}/valid_read.info',index=None,sep = '\t')
