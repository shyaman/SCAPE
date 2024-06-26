import os
import sys
import click
import signal
import tqdm
import time
import pandas as pd
import scipy.io, scipy.sparse
from numpy.lib.arraysetops import isin
from pathlib import Path


from loguru import logger
from multiprocessing import Pool
from multiprocessing.pool import MaybeEncodingError
from multiprocessing import set_start_method

from utils.utils import dict_list
from apamix.inferencepa_10x import wraper_process_10x, convert2sparse_10x
from apamix.inferencepa_bulk import wraper_process_bulk, convert2table


logger.add('apamix.log',
            rotation='10 MB',
            colorize=True,
            level="DEBUG")


@click.command()
@click.option(
    '--bed',
    type=str,
    help='The target regions (bed format) used for mix inference',
    required=True
    )
@click.option(
    '--bam',
    type=str,
    help='The bam file (sorted and indexed)',
    required=True
    )
@click.option(
    '--out',
    '-o',
    type=str,
    help='The output path',
    required=True
    )
@click.option(
    '--cores',
    type=int,
    help='Num (cores) of region are infering at once',
    default=1
    )
@click.option(
    '--cb',
    type=str,
    help='The cell barcode file, one cb for one line.',
    default=None,
    )
@click.option(
    '--tag',
    type=str,
    default='CB,UB,GN,GX',
    help='The cell barcode and UMI tag, for 10X: CB,UB,GN,GX.',
    required=True
    )

@click.option(
    '--n_max_apa',
    type=int,
    default='5',
    help='The maximum number of pA sites. Default value is 5.',
    required=True
    )

@click.option(
    '--n_min_apa',
    type=int,
    default='1',
    help='The minimum number of pA sites. Default value is 1.',
    required=True
    )

@click.option(
    '--max_utr_len',
    type=int,
    default=6000,
    help='The maximum length of UTR. Default value is 6000.'
    )

@click.option(
    '--la_dis_arr',
    type=str,
    default=None,
    help='The distinct lengths of polyA lengths in the dataset. Default: np.arange(10, 150, 10). User counld pass \'[10, 30, 50, 70, 90, 110, 130]\''
    )

@click.option(
    '--pmf_la_dis_arr',
    type=str,
    default=None,
    help='The the number of reads for each distinct polyA length. Default: Unif(10, 150). \'[309912, 4107929, 802856, 518229, 188316, 263208, 101]\''
    )

@click.option(
    '--mu_f',
    type=int,
    default=300,
    help='The mean of insert size for illumina library. Default: 300'
    )
@click.option(
    '--sigma_f',
    type=int,
    default=50,
    help='The std of insert size for illumina library. Default: 50'
    )

@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose mode'
)

@click.option(
    '--pseudocount',
    type=int,
    default=10,
    help='The pseudocount for the sparse matrix. Default: 10'
    )

@click.option(
    '--mode',
    type=str,
    default='10x',
    help='The input data type. Default: 10x'
    )

@click.option(
    '--create-table',
    is_flag=True,
    default=False,
    help='Only create the table'
    )
    
@click.option(
    '--freq-table',
    is_flag=True,
    default=False,
    help='create the frequency table'
    )

def readpa(
    bed,
    bam,
    out,
    cb,
    tag,
    cores,
    n_max_apa,
    n_min_apa,
    max_utr_len,
    la_dis_arr,
    pmf_la_dis_arr,
    mu_f,
    sigma_f,
    pseudocount,
    verbose,
    mode,
    create_table,
    freq_table,
    ):
    if not all([bed, bam, out]):
        cli(['readpa', '--help'])
        sys.exit(1)

    if not create_table:
        if not os.path.exists(os.path.join(out, 'tmp')):
            os.makedirs(os.path.join(out, 'tmp'))
            # os.makedirs(os.path.join(out, 'huge'))
            # os.makedirs(os.path.join(out, 'TooLongRegion'))
            # os.makedirs(os.path.join(out, 'TimeConsulting'))

        target_region = open(bed, 'r')
        res_lst = []
        # cb_df.cb = list(map(lambda x : x.split('-')[0], cb_df.cb.values))

        if la_dis_arr and pmf_la_dis_arr:
            la_dis_arr = eval(la_dis_arr)
            pmf_la_dis_arr = eval(pmf_la_dis_arr)
            if not (isinstance(la_dis_arr, list) and isinstance(pmf_la_dis_arr, list)):
                sys.exit('la_dis_arr and pmf_la_dis_arr were not list format in python.')
        else:
            la_dis_arr = None
            pmf_la_dis_arr = None

        peak_lst = []

        for i in target_region:
            if not i:
                continue
            chrom, st, en, strand, gene_id_fbed = i.strip().split('\t')
            if int(en) - int(st) + 1 > max_utr_len:
                logger.info(f'Skip more than {max_utr_len} UTR, {chrom}_{st}_{en}_{strand}')
                continue

            peak_lst.append(i.strip())

        target_region.close()

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

        pool = Pool(processes=cores)
        signal.signal(signal.SIGINT, original_sigint_handler)

        res_lst = []
        try:
            for x in range(len(peak_lst)):
                if mode == '10x':
                    cb_df = pd.read_csv(cb, names=['cb']) if cb else None
                    arg = [peak_lst[x], bam, cb_df, out, tag, verbose, n_max_apa, n_min_apa, la_dis_arr, pmf_la_dis_arr, mu_f, sigma_f, freq_table]
                    res_lst.append(pool.apply_async(wraper_process_10x, (arg,)))
                elif mode == 'bulk':
                    arg = [peak_lst[x], bam, out, verbose, n_max_apa, n_min_apa, la_dis_arr, pmf_la_dis_arr, mu_f, sigma_f]
                    res_lst.append(pool.apply_async(wraper_process_bulk, (arg,))) 

        except KeyboardInterrupt:
            logger.info('Caught KeyboardInterrupt, terminating workers')
            pool.terminate()
            sys.exit('KeyboardInterrupt')

        pool.close()
        pool.join()

        logger.info('Concating your final sheet')

        md, hd='w', True
        for df in res_lst:
            try:
                df = df.get()
            except MaybeEncodingError:
                continue

            if not isinstance(df, pd.DataFrame):
                continue

            # df = df.transpose()
            # if hd:
            #     df.columns = cb_df.cb.values.tolist()
            # df.to_csv(f'{out}/pasite.csv.gz', mode=md, header=hd, compression='gzip')
            df.to_csv(f'{out}/pa_lengths.csv', mode=md, header=hd, index=False)
            md, hd='a', False

    logger.info('Converting to sparse matrix/table')
    df = pd.read_csv(os.path.join(out, 'pa_lengths.csv'))
    if mode == '10x':
        sparse_matrix_adata = convert2sparse_10x([df, pseudocount])
        sparse_matrix_adata.write(os.path.join(out, 'adata_pa_sparse.h5ad'))
    elif mode == 'bulk':
        table = convert2table([df, pseudocount])
        table.to_csv(os.path.join(out, 'pa_table.csv'))

    logger.info('All done')
