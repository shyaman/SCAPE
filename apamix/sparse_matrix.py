from loguru import logger
import os
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse
import math


def convert_sparse_matrix(arg):
    df_pa, pseudocount = arg

    # TODO: optimize this
    df_pa_grouped = df_pa[['cb', 'umi', 'gene_id']].groupby(['gene_id','cb'])[['umi']]
    df_pa_grouped_cb_uniq_tot = df_pa_grouped.nunique()
    df_pa_grouped_onlypA = df_pa[df_pa['pa_len_arr']>0][['cb', 'umi', 'gene_id']].groupby(['gene_id','cb'])[['umi']]
    df_pa_grouped_cb_uniq = df_pa_grouped_onlypA.nunique()

    n_uniq_transcript_tot = df_pa_grouped_cb_uniq_tot['umi'].sum()
    n_uniq_transcript_pa = df_pa_grouped_cb_uniq['umi'].sum()
    logger.info(f'Sample pA pct: {n_uniq_transcript_pa*100/n_uniq_transcript_tot}')

    # if pseudocount == 'auto':
    #     pseudocount = math.floor(n_uniq_transcript_tot/n_uniq_transcript_pa)

    logger.info(f'Pseudocount: {pseudocount}')
    
    pct_df = df_pa_grouped_cb_uniq_tot.join(df_pa_grouped_cb_uniq, lsuffix='_totr', rsuffix='_npa').fillna(0)

    # for each cell, there should be atleast one read for a gene to be considered
    pct_df['pct_read'] = (pct_df['umi_npa'] + 1) / (pct_df['umi_totr'] + pseudocount)

    # convert to sparse matrix and then adata
    pa_df = pct_df[['pct_read']].reset_index()
    gene_id = pa_df['gene_id'].drop_duplicates().sort_values()
    cb = pa_df['cb'].drop_duplicates().sort_values()
    data = pa_df['pct_read'].tolist()
    row = pa_df['cb'].astype('category').cat.codes
    col = pa_df['gene_id'].astype('category').cat.codes
    adata_pa = ad.AnnData(sparse.csr_matrix((data, (row, col)), shape=(len(cb), len(gene_id))))
    adata_pa.obs.index = cb
    adata_pa.var.index = gene_id
    
    return adata_pa