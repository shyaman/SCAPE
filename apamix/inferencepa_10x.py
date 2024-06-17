import os
import time
import pysam
import pickle
import numpy as np
import pandas as pd
from loguru import logger
from collections import defaultdict
import anndata as ad
from scipy import sparse
import math

from apamix.apamix import APA
from utils.utils import dotdict, dict_list
from utils.bam import cigar_support_adj, collapse_intron, check_strand

def convert2sparse_10x(arg):
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

def run(arg):
    line, bamfile, cb_df, outdir, tag, verbose, n_max_apa, n_min_apa, LA_dis_arr, pmf_LA_dis_arr,mu_f,sigma_f = arg
    region_name = line
    chrom, left_site, right_site, strand, gene_id_fbed = line.split('\t')

    pickle_tmp = f'{outdir}/tmp/{chrom}_{left_site}_{right_site}_{strand}.pickle'

    cb_tag, umi_tag, gid_tag, ensid_tag= list(map(lambda x: x.strip(), tag.split(',')))

    t0 = time.time()
    if os.path.exists(pickle_tmp) and os.stat(pickle_tmp).st_size != 0:
        logger.info(f'cache found, return cache results. {pickle_tmp}')
        pickle_in = open(pickle_tmp, 'rb')
        try:
            res = pickle.load(pickle_in)
            res = res.drop('cb', axis=1).fillna(0).astype('int64')
            pickle_in.close()
            return res
        except pickle.UnpicklingError as e:
            pass

    bamfh = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)

    apa_reads = dotdict({
            # 'r2_len_arr': [],
            # 'r2_utr_st_arr': [],
            # 'r1_len_arr': [],
            # 'polya_len_arr': [],
            # 'pa_site_arr': [],
            'cb': [],
            'umi': [],
            'pa_len_arr': [],
            # 'pai_len_arr': [],
            'strand': [],
            'gene_id': [],
            'ensid': [],
            # 'sequence': [],
        })

    left_site = int(left_site)
    right_site = int(right_site)


    res = defaultdict(dict_list)

    for line in bamfh.fetch(chrom, int(left_site) - 1, int(right_site) + 1):

        # here to filter the score was not 255, just for STAR alignment methods.
        if line.mapping_quality != 255:
            continue

        if not check_strand(line, strand):
            continue

        start_site = line.reference_start + 1
        end_site = line.reference_end + 1
        query_name = line.query_name

        # ignore reads that are not in the region
        if start_site < left_site and end_site > right_site:
            continue

        if line.is_read1:
            res[query_name]['r1'].append(line)
        else:
            res[query_name]['r2'].append(line)

    for k, v in res.items():
        label = '_'.join(sorted(v.keys()))
        if label == 'r1':
            continue
        elif label == 'r2':
            # NOTE: 10x is only actually r2 only
            read2 = v[label]
            if len(read2) == 1:
                read2 = read2[0]
            else:
                # if there were multiple reads, continue
                continue

            r2_relative_start = collapse_intron(read2)
            if r2_relative_start > right_site and strand == '+':
                continue

            if r2_relative_start < left_site and strand == '-':
                continue

            pa_support, pa_len = cigar_support_adj(read2, strand)

            # try:
            #     dt = len(read2.get_tag('DT'))
            # except KeyError:
            #     dt = 0

            # try:
            #     pai = read2.get_tag(pa_tag)
            # except KeyError:
            #     pai = 0

            try:
                gid = read2.get_tag(gid_tag)
            except KeyError:
                gid = None
            
            try:
                ensid = read2.get_tag(ensid_tag)
            except KeyError:
                ensid = None

            # if strand == '+':
            #     r2_st = r2_relative_start - left_site +  1
            #     if pa_support == 'yes':
            #         pa_site = read2.reference_end + 1 - left_site
            #     else:
            #         pa_site = np.NaN
            # else:
            #     r2_st = right_site - r2_relative_start +  1
            #     if pa_support == 'yes':
            #         pa_site = right_site - read2.reference_start + 1
            #     else:
            #         pa_site = np.NaN
            # r1_dt = np.NaN
        else:
            r1 = v['r1']
            r2 = v['r2']
            if len(r1) != 1 or len(r2) != 1:
                continue
            read1 = r1[0]
            read2 = r2[0]

            # Actually this line was unused. because most of reads were single-end reads
            if max(read1.reference_end, read2.reference_end) - \
                min(read1.reference_start, read2.reference_start) > 1000:
                continue

            cigar_info = read1.cigar

            r2_relative_start = collapse_intron(read2)

            if r2_relative_start > right_site and strand == '+':
                continue

            if r2_relative_start < left_site and strand == '-':
                continue

            pa_support, pa_len = cigar_support_adj(read2, strand)

            # try:
            #     dt = len(read2.get_tag('DT'))
            #     if read2.is_reverse:
            #         if cigar_info[0][0] == 4:
            #             r1_dt = int(cigar_info[0][-1]) + dt
            #         else:
            #             r1_dt = dt
            #     else:
            #         if cigar_info[-1][0] == 4:
            #             r1_dt = int(cigar_info[-1][-1]) + dt
            #         else:
            #             r1_dt = dt

            # except KeyError:
            #     dt = 0
            
            # try:
            #     pai = read2.get_tag('pa')
            # except KeyError:
            #     pai = 0

            try:
                gid = read2.get_tag(gid_tag)
            except KeyError:
                gid = None
            
            try:
                ensid = read2.get_tag(ensid_tag)
            except KeyError:
                ensid = None

            # if strand == '+':
            #     pa_site = read1.reference_end + 1 - left_site
            #     r2_st = r2_relative_start - left_site +  1
            # else:
            #     pa_site = right_site - read1.reference_start + 1
            #     r2_st = right_site - r2_relative_start +  1
        try:
            # for 10X
            # append info if only read assigned to a gene
            if gid is None:
                continue
            cb = read2.get_tag(cb_tag)
            umi = read2.get_tag(umi_tag)
            # apa_reads.r2_utr_st_arr.append(r2_st)
            # apa_reads.r2_len_arr.append(read2.query_alignment_length)
            # apa_reads.r1_len_arr.append(dt)
            # apa_reads.polya_len_arr.append(r1_dt)
            # apa_reads.pa_site_arr.append(pa_site)
            apa_reads.cb.append(cb)
            apa_reads.umi.append(umi)
            apa_reads.pa_len_arr.append(pa_len)
            # apa_reads.pai_len_arr.append(pai)
            apa_reads.strand.append(strand)
            apa_reads.gene_id.append(gid)
            apa_reads.ensid.append(ensid)
        except KeyError:
            continue

    # if len(apa_reads.r2_utr_st_arr) < 50:
    # if len(apa_reads.pa_len_arr) < 50:
    #     return None
    # elif len(apa_reads.r2_utr_st_arr) > 1000000:
    # elif len(apa_reads.pa_len_arr) > 1000000:
    #     # skip the huge region
    #     # huge region was not the main reason for huge memory usage, so i commented this line
    #     print(f'Too many reads in {chrom}_{left_site}_{right_site}_{strand}_{gid}')
    #     huge_tmp = f'{outdir}/huge/{chrom}_{left_site}_{right_site}_{strand}'
    #     huge_in = open(huge_tmp, 'w')
    #     huge_in.write('\n')
    #     huge_in.close()
    #     return None

    # else:
    #     pass

    apa_reads = pd.DataFrame(apa_reads)

    # relative_offset = 150
    # apa_reads.r2_utr_st_arr = apa_reads.r2_utr_st_arr + relative_offset
    # apa_reads.pa_site_arr = apa_reads.pa_site_arr + relative_offset

    if cb_df is not None:
        cb_df = pd.read_csv(cb_df, names=['cb'])
        apa_reads = pd.merge(apa_reads, cb_df, on='cb', how='inner')

    return apa_reads

    # # here only add 100 bp for extending the utr length, just for decreasing time consumption.
    # utr_l = max([max(apa_reads.r2_utr_st_arr) + max(apa_reads.r2_len_arr) + 50, 1000])

    # test = APA(
    #     n_max_apa=int(n_max_apa),
    #     n_min_apa=int(n_min_apa),
    #     r1_utr_st_arr=apa_reads.r2_utr_st_arr,
    #     r1_len_arr=apa_reads.r2_len_arr,
    #     r2_len_arr=apa_reads.r1_len_arr,
    #     polya_len_arr=apa_reads.polya_len_arr,
    #     pa_site_arr=apa_reads.pa_site_arr,
    #     LA_dis_arr=LA_dis_arr,
    #     pmf_LA_dis_arr=pmf_LA_dis_arr,
    #     utr_len=utr_l,
    #     cb=apa_reads.cb,
    #     umi=apa_reads.umi,
    #     mu_f=mu_f,
    #     sigma_f=sigma_f,
    #     region_name=region_name,
    #     verbose=verbose
    # )

    # cb = apa_reads.cb
    # umi = apa_reads.umi
    # del apa_reads

    # logger.debug(f'Processing {chrom}, {left_site}, {right_site}, {strand}')

    # res = test.inference()
    # del test

    # df = pd.concat([pd.Series(res.label, name='label'), cb, umi], axis=1)
    # df.columns = ['label', 'cb', 'umi']

    # # remove last noise component
    # ignore_label = len(res.alpha_arr)

    # # Correct pA sites from relative position to genome position
    # if strand == '+':
    #     alpha_arr = left_site + res.alpha_arr - relative_offset - 1
    # else:
    #     alpha_arr = right_site - res.alpha_arr + relative_offset + 1
    # alpha_arr = [f'{chrom}:{int(alpha)}:{int(beta)}:{strand}' for alpha, beta in zip(alpha_arr, res.beta_arr.values.tolist())]

    # # final_res = defaultdict(int_dic)
    # df = df[df.label != ignore_label]
    # df = df.groupby(['cb', 'umi'])['label'].max().reset_index()

    # df = df.groupby(['cb', 'label'])['umi'].size().reset_index()# calculate umi counts
    # df = df.pivot_table(index=['cb'], columns='label', values='umi').reset_index()# make pivot table
    # rename_dict = dict(zip(range(len(alpha_arr)), alpha_arr))
    # df.rename(columns=rename_dict, inplace=True)
    # res = pd.merge(cb_df, df, on=['cb'], how='left')

    # pickle_out = open(pickle_tmp, 'wb')
    # pickle.dump(res, pickle_out)
    # pickle_out.close()

    # res = res.drop('cb', axis=1).fillna(0).astype('int64')

    # logger.debug(f'Done! {chrom}, {left_site}, {right_site}, {strand}')
    # return res



def wraper_process_10x(arg):
    return run(arg)
