import pandas as pd
import numpy as np

import os


def create_count_matrix(dname):
    """creates count martix for alignments
    dname: path to alignmnets"""
    fnames = os.listdir(dname)
    # OS specific hidden file
    try:
        fnames.remove('.DS_Store')
    except ValueError:
        pass
    # sample names
    snames = pd.Series(fnames).apply(lambda x: x.split('.')[0]).values
    assert len(snames) == len(set(snames)), 'Not unique sample names!'
    # just for get row names
    data_temp = pd.read_csv(dname+'/'+fnames[0], sep='\t', header=None, skiprows=4,
                            index_col=0)
    count_mtx = pd.DataFrame(index=data_temp.index, columns=snames)
    for fname in fnames:
        data_temp = pd.read_csv(dname+'/'+fname, sep='\t', header=None,
                                skiprows=4, index_col=0)
        sname = fname.split('.')[0]
        count_mtx.loc[data_temp.index, sname] = data_temp[1]
    return count_mtx

def translate_gene_ids(data):
    """ replaces ensmbl gene ids with HGNC symbols
    removes non mappable genes"""
    gene_dict = pd.read_csv('../results/hgnc_symbol.csv', sep=',',
                           header=0, index_col=0)
    fil = ~(gene_dict['hgnc_symbol'].isna())
    gene_dict = gene_dict[fil]
    gene_dict = gene_dict.drop_duplicates('ensembl_gene_id')
    gene_dict = gene_dict.drop_duplicates('hgnc_symbol')
    gene_dict.index = gene_dict['ensembl_gene_id']
    gene_dict = gene_dict['hgnc_symbol']
    genes = list(set(data.index) & set(gene_dict.index))
    data = data.loc[genes]
    data.index = gene_dict[data.index]
    return data

def remove_0_genes(data):
    """ removes genes with 0 expression in all samples"""
    fil = data.sum(1) > 0
    return data[fil]
