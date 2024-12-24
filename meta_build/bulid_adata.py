#######
# sample_proxy = "IBD"
# meta_file = "/storage/fengyu/liuzhong/Cyclone/IBD1-B03501C4/metadata/IBD1-B03501C4.meta.gz"
# cell_meta_file = "/home/zangyupeng/05.ExtendDisk/00.Data/10.ZH/pancancer/IBD_isotype.csv"
# rna_adata = '/home/zangyupeng/05.ExtendDisk/00.Data/07.NBT/ST/IBD/B03501C4_norm.h5ad'
# cellbin_file = "/storage/fengyu/liuzhong/Cellpose/20240725/IBD1-B03501C4/GUI/IBD1-B03501C4.registration.cellpose_fill_expand20.csv"
# outline = "/storage/fengyu/liuzhong/Cellpose/20240725/IBD1-B03501C4/GUI/IBD1-B03501C4.registration.cellpose_outlines_expand20.csv"
# outadatadir = "./ibd_umi"
# binsize = 50

#######
import argparse
import sys


parser = argparse.ArgumentParser(description="Build adata")

parser.add_argument('--sample_proxy', type=str, required=True, 
                    help="Sample proxy same in cyclone analysis")
parser.add_argument('--meta_file', type=str, required=True, 
                    help="Meta file path")
parser.add_argument('--cell_meta_file', type=str, required=True, 
                    help="Cell Meta file path")
parser.add_argument('--rna_adata', type=str, required=True, 
                    help="RNA_adata")
parser.add_argument('--cellbin_file', type=str, required=True, 
                    help="cellbin_file file path")
parser.add_argument('--outline', type=str, required=True, 
                    help="outline file path")
parser.add_argument('--outadatadir', type=str, required=True, 
                    help="Directory to out adata")
parser.add_argument('--bin', type=int, default=50, 
                    help="binsize")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
    
args = parser.parse_args()


sample_proxy = args.sample_proxy
meta_file = args.meta_file
cell_meta_file = args.cell_meta_file
rna_adata = args.rna_adata
cellbin_file = args.cellbin_file
outline = args.outline
outadatadir = args.outadatadir
binsize = args.binsize


import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad
import scanpy as sc
import os 
from anndata import AnnData
import concave_hull

def RestAdata(NeedAddAdata,TotalObs):

    NeedAddAdata = NeedAddAdata[NeedAddAdata.obs_names.intersection(TotalObs)]
    restcell_set = set(TotalObs) - set(NeedAddAdata.obs_names)
    restcell = AnnData(csr_matrix((len(restcell_set), NeedAddAdata.shape[1])))
    restcell.obs_names = restcell_set
    restcell.var_names = NeedAddAdata.var_names
    imudata = sc.concat([NeedAddAdata,restcell])
    return imudata

def ExtractAdata(obs,var,value,meta):
    id_codes, binid_uniques = pd.factorize(meta[obs])
    cdr3_codes, cdr3_uniques = pd.factorize(meta[var])
    row = id_codes
    col = cdr3_codes
    count_data = meta[value]
    sparse_matrix = csr_matrix((count_data, (row, col)), shape=(len(binid_uniques), len(cdr3_uniques)))
    adata = AnnData(X=sparse_matrix)
    adata.obs_names = binid_uniques
    adata.var_names = cdr3_uniques
    return adata

def CDR3info(adata):
    adata.var['CDR3'] = adata.var_names.map(lambda x : x.split('@')[0])
    adata.var['Chain'] = adata.var_names.map(lambda x : x.split('@')[1])
    adata.var['Func'] = adata.var['CDR3'].apply(lambda x : 'Non' if '*' in x or '_' in x else 'Func')
    return adata

def AddCDR3uns(adata,meta,obsid = 'binid'):
    Total_CDR3_to_ID = {}
    Total_ID_to_CDR3 = {}
    statdf = meta[['CDR3PriName',obsid,'c_priisotype',"c_isotype"]].drop_duplicates(['CDR3PriName',obsid,'c_priisotype',"c_isotype"])
    for tmpchain in statdf['c_isotype'].unique():
        tmp = statdf[statdf['c_isotype'] == tmpchain]
        cdr3_to_id = tmp.groupby('CDR3PriName')[obsid].apply(list).to_dict()
        id_to_cdr3 = tmp.groupby(obsid)['CDR3PriName'].apply(list).to_dict()
        Total_CDR3_to_ID[tmpchain] = cdr3_to_id
        Total_ID_to_CDR3[tmpchain] = id_to_cdr3

    HeavyChain_CDR3_to_ID = {}
    HeavyChain_ID_to_CDR3 = {}

    for tmpchain in statdf[statdf['c_priisotype'].map(len)==4]['c_priisotype'].unique():
        tmp = statdf[statdf['c_priisotype'] == tmpchain]
        cdr3_to_id = tmp.groupby('CDR3PriName')[obsid].apply(list).to_dict()
        id_to_cdr3 = tmp.groupby(obsid)['CDR3PriName'].apply(list).to_dict()
        HeavyChain_CDR3_to_ID[tmpchain] = cdr3_to_id
        HeavyChain_ID_to_CDR3[tmpchain] = id_to_cdr3
    adata.uns['XCR'] = {'c2i':Total_CDR3_to_ID,'i2c':Total_ID_to_CDR3,'Hc2i':HeavyChain_CDR3_to_ID,'Hi2c':HeavyChain_ID_to_CDR3}
    return adata

def BulidNGSadata(meta,rna,binsize = 50):
    overlapmeta = meta.copy()
    overlapmeta['binid'] = 'DNB_' + (overlapmeta['x']// binsize* binsize).map(str) + '_' + (overlapmeta['y']// binsize* binsize).map(str)
    overlapmeta = overlapmeta[overlapmeta['binid'].isin(rna.obs_names)].copy()
    overlapmeta['count'] = 1
    adata = ExtractAdata(obs = 'binid',var = 'CDR3PriName',value = 'count',meta = overlapmeta)
    adata = RestAdata(adata,rna.obs_names)
    adata = adata[rna.obs_names].copy()
    adata.obs['x'] = adata.obs_names.map(lambda x: int(x.split('_')[1]))
    adata.obs['y'] = adata.obs_names.map(lambda x: int(x.split('_')[2]))
    adata.obsm['spatial'] = adata.obs[['x','y']].values
    adata = CDR3info(adata)
    adata = AddCDR3uns(adata,meta = overlapmeta)
    return adata

def AddCellInfo(adata,rna,outline,tissuecellmask):
    outline = pd.read_csv(outline,index_col = 0)
    outline['cellID'] = 'Cell.'+outline['cellID'].map(str)
    tissueoutline = outline[outline['cellID'].isin(adata.obs_names)].copy()
    tissueoutline = tissueoutline.reset_index(drop=True)
    tissuecellmask = tissuecellmask.reset_index(drop=True)

    adata.uns['outline'] = tissueoutline
    adata.uns['cellpose'] = tissuecellmask[['cellID','x','y']]
    
    
    if 'x' in rna.obs.columns:
        points = rna.obs[['x','y']].values
    else:
        if 'X_spatial' in rna.obsm.keys():
            key = 'X_spatial'
        elif 'spatial' in rna.obsm.keys():
            key = 'spatial'

        points = rna.obsm[key]
        
    hull = concave_hull.concave_hull(points,concavity = 0.1)
    adata.uns['hull'] = hull
    return adata

def Bulidnpdata(meta,rna,cellbin_file,outline,binsize = 50):
    overlapmeta = meta.copy()
    cellmask = pd.read_csv(cellbin_file,index_col = 0)
    cellmask['cellID'] = 'Cell.' + cellmask['cellID'].map(str)
    cellmask['pos'] = cellmask['x'].map(str)+'_'+cellmask['y'].map(str)
    cellmask['binloc'] = 'DNB_' + (cellmask['x']//binsize*binsize).map(str) + '_' + (cellmask['y']//binsize*binsize).map(str)
    tissuecellmask = cellmask[cellmask['binloc'].isin(rna.obs_names)].copy()
    tissuecell = tissuecellmask['cellID'].drop_duplicates().tolist()
    celldict = dict(zip(cellmask['pos'],cellmask['cellID']))
    centroid = cellmask.groupby('cellID').agg({'x': 'mean', 'y': 'mean'}).reset_index()

    overlapmeta['id'] = overlapmeta['loc'].map(celldict)
    overlapmeta = overlapmeta[~overlapmeta['id'].isna()].copy()
    overlapmeta['count'] = 1
    
    adata = ExtractAdata(obs = 'id',var = 'CDR3PriName',value = 'count',meta = overlapmeta)
    adata = RestAdata(adata,tissuecell)
    adata.obs['x'] = adata.obs_names.map(dict(zip(centroid['cellID'],centroid['x'])))
    adata.obs['y'] = adata.obs_names.map(dict(zip(centroid['cellID'],centroid['y'])))
    adata.obsm['spatial'] = adata.obs[['x','y']].values
    adata = CDR3info(adata)
    adata = AddCellInfo(adata,rna,outline,tissuecellmask)
    adata = AddCDR3uns(adata,meta = overlapmeta,obsid = 'id')
    return adata

cell_meta = pd.read_csv(cell_meta_file)#compression='gzip')
cell_meta['CDR3PriName'] = cell_meta['cdr3@c_isotype']
cell_meta['final_func'] = cell_meta['CDR3PriName'].apply(lambda x : 'Non' if '*' in x or '_' in x else 'Func')
cell_meta = cell_meta[cell_meta['final_func'] == 'Func'].copy()


cellmask = pd.read_csv(cellbin_file,index_col = 0)
cellmask['cellID'] = 'Cell.' + cellmask['cellID'].map(str)
cellmask['pos'] = cellmask['x'].map(str)+'_'+cellmask['y'].map(str)
cellmask['binloc'] = 'DNB_' + (cellmask['x']//binsize*binsize).map(str) + '_' + (cellmask['y']//binsize*binsize).map(str)


tissuecellmask = cellmask[cellmask['binloc'].isin(rna.obs_names)].copy()
tissuecell = tissuecellmask['cellID'].drop_duplicates().tolist()
celldict = dict(zip(cellmask['pos'],cellmask['cellID']))


meta = pd.read_csv(meta_file,compression='gzip')
meta['c_isotype'] = meta['isotype']
meta['c_priisotype'] = meta['priisotype']
meta['cdr3@c_isotype'] = meta['cdr3name']
meta['CDR3PriName'] = meta['cdr3@c_isotype'].map(lambda x: x.split('@')[0])+ '@'+meta['c_isotype']
meta['final_func'] = meta['CDR3PriName'].apply(lambda x : 'Non' if '*' in x or '_' in x else 'Func')
meta = meta[meta['final_func'] == 'Func'].copy()
meta['id'] = meta['loc'].map(celldict)


# load RNA data for tissue cut

rna = sc.read_h5ad(rna_adata)
try:
    rna.obs_names = rna.obs_names.map(lambda x : 'DNB_' + x.split('@')[1])
except:
    pass

ngsadata = BulidNGSadata(meta,rna,binsize = binsize)
binadatafile = os.path.join(outadatadir,f'{sample_proxy}.bin{binsize}.h5ad')
ngsadata.write(binadatafile)


npdata = Bulidnpdata(cell_meta,rna,cellbin_file,outline,binsize = binsize)
npadatafile = os.path.join(outadatadir,f'{sample_proxy}.cellbin.h5ad')
npdata.write(npadatafile)