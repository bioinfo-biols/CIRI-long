#!/home/zhangjy/.virtualenvs/Benchmarking/bin/python
import os
import sys
import pickle
import argparse
import json
from collections import defaultdict


def main():
    from CIRI.version import __version__
    from CIRI.rofinder import find_ccs_reads, load_ccs_reads
    from CIRI.alignment import index_annotation, scan_ccs_reads, recover_ccs_reads
    from CIRI.alignment import scan_raw_reads
    from CIRI.utils import check_file, check_dir
    from CIRI.logger import get_logger

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in', dest='input', metavar='READS', default=None,
                        help='Input reads.fq.gz', )
    parser.add_argument('-o', '--out', dest='output', metavar='DIR', default=None,
                        help='Output directory, default: ./', )
    parser.add_argument('-r', '--ref', dest='reference', metavar='REF', default=None,
                        help='Reference genome FASTA file', )
    parser.add_argument('-p', '--prefix', dest='prefix', metavar='PREFIX', default="CIRI-long",
                        help='Output sample prefix, (default: %(default)s)', )
    parser.add_argument('-a', '--anno', dest='gtf', metavar='GTF', default=None,
                        help='Genome reference gtf', )
    parser.add_argument('--canonical', dest='canonical', default=True, action='store_true',
                        help='Use canonical splice signal (GT/AG) only, default: %(default)s)')
    parser.add_argument('-t', '--threads', dest='threads', metavar='INT', default=os.cpu_count(),
                        help='Number of threads', )
    parser.add_argument('--debug', dest='debug', default=False, action='store_true',
                        help='Run in debuggin mode, (default: %(default)s)', )
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()

    if args.input is None or args.output is None:
        sys.exit('Please provide input and output file, run CIRI-long using -h or --help for detailed information.')
    if args.reference is None:
        sys.exit('Please specific FASTA of reference genome')

    # Check parameters
    in_file = check_file(args.input)
    gtf_file = None if args.gtf is None else check_file(args.gtf)
    out_dir = check_dir(args.output)
    ref_fasta = check_file(args.reference)
    check_dir(out_dir + '/tmp')
    prefix = args.prefix
    threads = int(args.threads)
    debugging = args.debug
    is_canonical = args.canonical

    logger = get_logger('CIRI-long', fname='{}/{}.log'.format(out_dir, prefix), verbosity=debugging)
    logger.info('Input reads: ' + os.path.basename(in_file))
    logger.info('Output directory: ' + os.path.basename(out_dir))
    logger.info('Multi threads: {}'.format(args.threads))

    # Scan for repeats and CCS
    reads_count = defaultdict(int)
    if os.path.exists('{}/tmp/{}.ccs.fa'.format(out_dir, prefix)) and os.path.exists('{}/tmp/{}.raw.fa'.format(out_dir, prefix)):
        logger.info('Step 1 - Loading circRNA candidates in previous run')
        ccs_seq = load_ccs_reads(out_dir, prefix)
        reads_count['consensus'] = len(ccs_seq)
    else:
        logger.info('Step 1 - Scanning raw reads to find circRNA candidates')
        total_reads, ro_reads, ccs_seq = find_ccs_reads(in_file, out_dir, prefix, threads, debugging)
        reads_count['total'] = total_reads
        reads_count['consensus'] = ro_reads

    if 'total' in reads_count:
        logger.info('Total Reads: {}'.format(reads_count['total']))
    logger.info('Cyclic Consensus Reads: {}'.format(reads_count['consensus']))

    # generate index of splice site and annotation
    if gtf_file is None:
        logger.warn('No genome annotation provided, entering \'De novo\' mode')
        gtf_idx, ss_idx = None, None
    else:
        idx_file = out_dir + '/tmp/ss.idx'
        if os.path.exists(idx_file):
            logger.info('Loading pre-built splice site index from: {}'.format(idx_file))
            with open(idx_file, 'rb') as idx:
                gtf_idx, ss_idx = pickle.load(idx)
        else:
            gtf_idx, ss_idx = index_annotation(gtf_file)
            with open(idx_file, 'wb') as idx:
                pickle.dump([gtf_idx, ss_idx], idx)

    # Find circRNAs
    logger.info('Step 2.1 - Find circRNAs from CCS reads')
    tmp_cnt, short_seq = scan_ccs_reads(ccs_seq, ref_fasta, ss_idx, is_canonical, out_dir, prefix, threads)
    for key, value in tmp_cnt.items():
        reads_count[key] += value

    # Recover short reads
    logger.info('Step 2.2 - Recover short CCS reads')
    tmp_cnt = recover_ccs_reads(short_seq, ref_fasta, ss_idx, is_canonical, out_dir, prefix, threads)
    for key, value in tmp_cnt.items():
        reads_count[key] += value

    # Find BSJs
    logger.info('Step 3 - Find circRNAs with partial structure')
    tmp_cnt, short_seq = scan_raw_reads(in_file, ref_fasta, ss_idx, is_canonical, out_dir, prefix, threads)
    for key, value in tmp_cnt.items():
        reads_count[key] += value

    # logger.info('Step 3.2 - Second scanning for BSJs')
    # tmp_cnt = recover_raw_reads(short_seq, ref_fasta, ss_idx, is_canonical, out_dir, prefix, threads)

    logger.info('Raw unmapped: {}'.format(reads_count['raw_unmapped']))
    logger.info('CCS mapped: {}'.format(reads_count['ccs_mapped']))
    logger.info('BSJ: {}'.format(reads_count['bsj']))
    logger.info('Splice signal: {}'.format(reads_count['signal']))
    logger.info('Partial reads: {}'.format(reads_count['partial']))

    with open('{}/{}.json'.format(out_dir, prefix), 'w') as f:
        json.dump(reads_count, f)

    logger.info('All Finished!')


if __name__ == '__main__':
    main()
