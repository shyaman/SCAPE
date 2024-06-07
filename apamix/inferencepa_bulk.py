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
from utils.bam import cigar_support_bulk, collapse_intron, check_strand

def convert2table(arg):
    df_pa, pseudocount = arg

    # TODO: optimize this
    df_pa_grouped = df_pa.groupby(['gene_id'])
    df_pa_grouped_cb_uniq_tot = df_pa_grouped.size().to_frame('count')
    df_pa_grouped_onlypA = df_pa[df_pa['pa_len_arr']>0].groupby(['gene_id'])
    df_pa_grouped_cb_uniq = df_pa_grouped_onlypA.size().to_frame('count')

    n_uniq_transcript_tot = df_pa_grouped_cb_uniq_tot['count'].sum()
    n_uniq_transcript_pa = df_pa_grouped_cb_uniq['count'].sum()
    logger.info('Sample pA pct: {:.2f}%'.format(n_uniq_transcript_pa*100/n_uniq_transcript_tot))

    logger.info(f'Pseudocount: {pseudocount}')
    
    pct_df = df_pa_grouped_cb_uniq_tot.join(df_pa_grouped_cb_uniq, lsuffix='_totr', rsuffix='_npa').fillna(0)

    # for each cell, there should be atleast one read for a gene to be considered
    pct_df['pct_read'] = (pct_df['count_npa'] + 1) / (pct_df['count_totr'] + pseudocount)
    
    return pct_df

def run(arg):
    line, bamfile, outdir, verbose, n_max_apa, n_min_apa, LA_dis_arr, pmf_LA_dis_arr,mu_f,sigma_f = arg
    region_name = line
    chrom, left_site, right_site, strand, gene_id_fbed = line.split('\t')

    pickle_tmp = f'{outdir}/tmp/{chrom}_{left_site}_{right_site}_{strand}.pickle'

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
            'pa_len_arr': [],
            'strand': [],
            'gene_id': [],
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
            # read1 = v[label]
            # if len(read1) == 1:
            #     read1 = read1[0]
            # else:
            #     # if there were multiple reads, continue
            #     continue

            # r1_relative_start = collapse_intron(read1)
            # if r1_relative_start > right_site and strand == '+':
            #     continue

            # if r1_relative_start < left_site and strand == '-':
            #     continue

            # pa_support, pa_len = cigar_support_bulk(read1, strand)

            # if pa_support=='yes':
            #     print(pa_support, pa_len, k, gene_id_fbed)
        elif label == 'r2':
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

            pa_support, pa_len = cigar_support_bulk(read2, strand)

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

            pa_support_read2, pa_len_read2 = cigar_support_bulk(read2, strand)
            pa_support_read1, pa_len_read1 = cigar_support_bulk(read1, strand)

            if pa_support_read1 == 'yes' and pa_support_read2 == 'yes':
                pa_len = max(pa_len_read1, pa_len_read2)
            elif pa_support_read1 == 'yes':
                pa_len = pa_len_read1
            elif pa_support_read2 == 'yes':
                pa_len = pa_len_read2
            else:
                pa_len = 0

        try:
            apa_reads.pa_len_arr.append(pa_len)
            apa_reads.strand.append(strand)
            apa_reads.gene_id.append(gene_id_fbed)
        except KeyError:
            continue

    apa_reads = pd.DataFrame(apa_reads)

    return apa_reads

def wraper_process_bulk(arg):
    return run(arg)
