#!/home/zhangjy/.virtualenvs/Benchmarking/bin/python
import os
import sys
import pickle
import json
from collections import defaultdict


def call(args):
    from CIRI_long.logger import get_logger
    from CIRI_long.utils import check_file, check_dir
    from CIRI_long.align import index_annotation, index_circ
    from CIRI_long.find_ccs import find_ccs_reads, load_ccs_reads
    from CIRI_long.find_bsj import scan_ccs_reads, recover_ccs_reads
    from CIRI_long.find_bsj import scan_raw_reads
    from CIRI_long.pipeline import run_ccs

    lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
    os.environ['PATH'] = lib_path + ':' + os.environ['PATH']
    os.chmod(lib_path + '/ccs', 0o755)

    if args.input is None or args.output is None:
        sys.exit('Please provide input and output file, run CIRI-long using -h or --help for detailed information.')
    if args.reference is None:
        sys.exit('Please specific FASTA of reference genome')

    # Check parameters
    in_file = check_file(args.input)
    gtf_file = None if args.gtf is None else check_file(args.gtf)
    circ_file = None if args.circ is None else check_file(args.circ)
    out_dir = check_dir(args.output)
    ref_fasta = check_file(args.reference)
    check_dir(out_dir + '/tmp')
    prefix = args.prefix
    threads = int(args.threads)
    debugging = args.debug
    is_canonical = args.canonical

    logger = get_logger('CIRI-long', fname='{}/{}.log'.format(out_dir, prefix), verbosity=debugging)
    logger.info('----------------- Input paramters ------------------')
    logger.info('Input reads: ' + os.path.basename(in_file))
    logger.info('Output directory: ' + os.path.basename(out_dir))
    logger.info('Multi threads: {}'.format(args.threads))
    logger.info('----------------- Calling circRNAs -----------------')

    # Scan for repeats and CCS
    reads_count = defaultdict(int)
    is_fast = 0

    if not debugging and os.path.exists('{}/tmp/{}.ccs.fa'.format(out_dir, prefix)) and os.path.exists('{}/tmp/{}.raw.fa'.format(out_dir, prefix)):
        logger.info('Step 1 - Loading circRNA candidates in previous run')
        ccs_seq = load_ccs_reads(out_dir, prefix)
        reads_count['consensus'] = len(ccs_seq)
    else:
        try:
            run_ccs(in_file, out_dir, prefix, threads, debugging)
            ccs_seq = load_ccs_reads(out_dir, prefix)
            reads_count['consensus'] = len(ccs_seq)
            is_fast = 1
        except Exception as e:
            logger.warn('Failed to run ccs command in fast mode, falling back to slower python version.')
            total_reads, ro_reads, ccs_seq = find_ccs_reads(in_file, out_dir, prefix, threads, debugging)
            reads_count['total'] = total_reads
            reads_count['consensus'] = ro_reads
            is_fast = 0

    if is_fast == 1 and reads_count['consensus'] == 0:
        logger.warn('Failed to run ccs command in fast mode, try again in slower python version.')
        total_reads, ro_reads, ccs_seq = find_ccs_reads(in_file, out_dir, prefix, threads, debugging)
        reads_count['total'] = total_reads
        reads_count['consensus'] = ro_reads

    if 'total' in reads_count:
        logger.info('Total Reads: {}'.format(reads_count['total']))
    logger.info('Cyclic Consensus Reads: {}'.format(reads_count['consensus']))

    # generate index of splice site and annotation
    if gtf_file is None and circ_file is None:
        logger.warn('No annotation provided, entering \'De novo\' mode')
        gtf_idx, ss_idx = None, None
    else:
        idx_file = out_dir + '/tmp/ss.idx'
        if os.path.exists(idx_file):
            logger.info('Loading pre-built splice site index from: {}'.format(idx_file))
            with open(idx_file, 'rb') as idx:
                gtf_idx, intron_idx, ss_idx = pickle.load(idx)
        else:
            if gtf_file is not None:
                gtf_idx, intron_idx, ss_idx = index_annotation(gtf_file)
            else:
                gtf_idx, intron_idx, ss_idx = None, None, None
            if circ_file is not None:
                ss_idx = index_circ(circ_file, ss_idx)

            with open(idx_file, 'wb') as idx:
                pickle.dump([gtf_idx, intron_idx, ss_idx], idx, -1)

    # Find circRNAs
    logger.info('Step 2.1 - Find circRNAs from CCS reads')
    tmp_cnt, short_seq = scan_ccs_reads(ccs_seq, ref_fasta, ss_idx, gtf_idx, intron_idx, is_canonical, out_dir, prefix, threads)
    for key, value in tmp_cnt.items():
        reads_count[key] += value

    # Recover short reads
    logger.info('Step 2.2 - Recover short CCS reads')
    tmp_cnt = recover_ccs_reads(short_seq, ref_fasta, ss_idx, gtf_idx, intron_idx, is_canonical, out_dir, prefix, threads)
    for key, value in tmp_cnt.items():
        reads_count[key] += value

    # Find BSJs
    logger.info('Step 3 - Find circRNAs with partial structure')
    tmp_cnt, short_seq = scan_raw_reads(in_file, ref_fasta, gtf_idx, intron_idx, ss_idx, is_canonical, out_dir, prefix, threads)
    for key, value in tmp_cnt.items():
        reads_count[key] += value

    logger.info('Raw unmapped: {}'.format(reads_count['raw_unmapped']))
    logger.info('CCS mapped: {}'.format(reads_count['ccs_mapped']))
    logger.info('BSJ: {}'.format(reads_count['bsj']))
    logger.info('Splice signal: {}'.format(reads_count['signal']))
    logger.info('Partial reads: {}'.format(reads_count['partial']))

    with open('{}/{}.json'.format(out_dir, prefix), 'w') as f:
        json.dump(reads_count, f)

    logger.info('Calling circRNAs finished!')


def collapse(args):
    from CIRI_long.logger import get_logger
    from CIRI_long.utils import check_file, check_dir
    from CIRI_long.align import index_annotation, index_circ
    from CIRI_long import collapse

    if args.input is None or args.output is None:
        sys.exit('Please provide input and output file, run CIRI-long using -h or --help for detailed information.')

    in_file = check_file(args.input)
    out_dir = check_dir(args.output)
    check_dir(out_dir + '/tmp')
    prefix = args.prefix

    gtf_file = None if args.gtf is None else check_file(args.gtf)
    circ_file = None if args.circ is None else check_file(args.circ)
    ref_fasta = check_file(args.reference)

    threads = int(args.threads)
    debugging = args.debug
    is_canonical = args.canonical

    logger = get_logger('CIRI-long', fname='{}/{}.log'.format(out_dir, prefix), verbosity=debugging)
    logger.info('----------------- Input paramters ------------------')
    logger.info('Input reads: ' + os.path.basename(in_file))
    logger.info('Output directory: ' + os.path.basename(out_dir))
    logger.info('Multi threads: {}'.format(args.threads))
    logger.info('-------------- Collapse circular reads -------------')

    # generate index of splice site and annotation
    if gtf_file is None and circ_file is None:
        logger.warn('No annotation provided, entering \'De novo\' mode')
        gtf_idx, ss_idx = None, None
    else:
        idx_file = out_dir + '/tmp/ss.idx'
        if os.path.exists(idx_file):
            logger.info('Loading pre-built splice site index from: {}'.format(idx_file))
            with open(idx_file, 'rb') as idx:
                gtf_idx, intron_idx, ss_idx = pickle.load(idx)
        else:
            if gtf_file is not None:
                gtf_idx, intron_idx, ss_idx = index_annotation(gtf_file)
            else:
                gtf_idx, intron_idx, ss_idx = None, None, None
            if circ_file is not None:
                ss_idx = index_circ(circ_file, ss_idx)

            with open(idx_file, 'wb') as idx:
                pickle.dump([gtf_idx, intron_idx, ss_idx], idx, -1)

    # Load reads
    cand_reads = collapse.load_cand_circ(in_file)

    # Consensus reads
    corrected_file = '{}/tmp/{}.corrected.pkl'.format(out_dir, prefix)
    if not debugging and os.path.exists(corrected_file):
        logger.info('Step 1 - Loading clustered circular reads in previous run')
        with open(corrected_file, 'rb') as pkl:
            circ_num, corrected_reads = pickle.load(pkl)
    else:
        logger.info('Step 1 - Clustering candidate circular reads')
        # Cluster reads
        reads_cluster = collapse.cluster_reads(cand_reads)
        logger.info('Circular reads clusters: {}'.format(len(reads_cluster)))

        # Generate consensus reads
        circ_num, corrected_reads = collapse.correct_reads(reads_cluster, ref_fasta, gtf_idx, intron_idx, ss_idx, threads)
        with open(corrected_file, 'wb') as pkl:
            pickle.dump([circ_num, corrected_reads], pkl, -1)
        logger.info('Corrected clusters: {}, {}/{}/{}/{} annotated/denovo/lariat/unknown'.format(
            len(corrected_reads), circ_num['Annotated'], circ_num['Denovo signal'],
            circ_num['High confidence lariat'], circ_num['Unknown signal']))

    logger.info('Step 2 - Calculating expression matrix')
    circ_cnt = collapse.cal_exp_mtx(cand_reads, corrected_reads, ref_fasta, gtf_idx, out_dir, prefix)
    logger.info('Final circRNAs: {}'.format(circ_cnt))

    # # Find circRNAs again
    # logger.info('Correct circRNAs from consensus reads!')
    # corrected_circ, short_reads = scan_corrected_reads(corrected_reads, ref_fasta, ss_idx, threads)
    #
    # # Recover short circRNAs
    # logger.info('Recover short circRNAs from consensus reads!')
    # corrected_circ += recover_corrected_reads(short_reads, ref_fasta, ss_idx)
    #
    # logger.info('Clustered circRNAs: {}'.format(len(corrected_reads)))
    # logger.info('Corrected circRNAs: {}'.format(len(corrected_circ)))
    #
    # with open('{}/{}_circ.pkl'.format(out_dir, prefix), 'wb') as idx:
    #     pickle.dump(corrected_circ, idx, -1)

    logger.info('Correction of Back-Spliced Junctions finished!')


def main():
    import argparse
    from CIRI_long.version import __version__

    # Init parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s v{version}'.format(version=__version__))

    # Init subparsers
    subparsers = parser.add_subparsers(help='commands')

    # Parsers for circRNA calling
    call_parser = subparsers.add_parser('call')
    call_parser.add_argument('-i', '--in', dest='input', metavar='READS', default=None,
                             help='Input reads.fq.gz', )
    call_parser.add_argument('-o', '--out', dest='output', metavar='DIR', default=None,
                             help='Output directory, default: ./', )
    call_parser.add_argument('-r', '--ref', dest='reference', metavar='REF', default=None,
                             help='Reference genome FASTA file', )
    call_parser.add_argument('-p', '--prefix', dest='prefix', metavar='PREFIX', default="CIRI-long",
                             help='Output sample prefix, (default: %(default)s)', )
    call_parser.add_argument('-a', '--anno', dest='gtf', metavar='GTF', default=None,
                             help='Genome reference gtf, (optional)', )
    call_parser.add_argument('-c', '--circ', dest='circ', metavar='CIRC', default=None,
                             help='Additional circRNA annotation in bed/gtf format, (optional)', )
    call_parser.add_argument('--canonical', dest='canonical', default=True, action='store_true',
                             help='Use canonical splice signal (GT/AG) only, default: %(default)s)')
    call_parser.add_argument('-t', '--threads', dest='threads', metavar='INT', default=os.cpu_count(),
                             help='Number of threads, (default: use all cores)', )
    call_parser.add_argument('--debug', dest='debug', default=False, action='store_true',
                             help='Run in debugging mode, (default: %(default)s)', )
    call_parser.set_defaults(func=call)

    # Parsers for concat circular reads
    collapse_parser = subparsers.add_parser('collapse')
    collapse_parser.add_argument('-i', '--in', dest='input', metavar='LIST', default=None,
                                 help='Input list of CIRI-long results', )
    collapse_parser.add_argument('-o', '--out', dest='output', metavar='DIR', default=None,
                                 help='Output directory, default: ./', )
    collapse_parser.add_argument('-p', '--prefix', dest='prefix', metavar='PREFIX', default="CIRI-long",
                                 help='Output sample prefix, (default: %(default)s)', )
    collapse_parser.add_argument('-r', '--ref', dest='reference', metavar='REF', default=None,
                                 help='Reference genome FASTA file', )
    collapse_parser.add_argument('-a', '--anno', dest='gtf', metavar='GTF', default=None,
                                 help='Genome reference gtf, (optional)', )
    collapse_parser.add_argument('-c', '--circ', dest='circ', metavar='CIRC', default=None,
                                 help='Additional circRNA annotation in bed/gtf format, (optional)', )
    collapse_parser.add_argument('--canonical', dest='canonical', default=True, action='store_true',
                                 help='Use canonical splice signal (GT/AG) only, default: %(default)s)')
    collapse_parser.add_argument('-t', '--threads', dest='threads', metavar='INT', default=os.cpu_count(),
                                 help='Number of threads, (default: use all cores)', )
    collapse_parser.add_argument('--debug', dest='debug', default=False, action='store_true',
                                 help='Run in debugging mode, (default: %(default)s)', )

    collapse_parser.set_defaults(func=collapse)

    # Parse parsers
    args = parser.parse_args()

    # Run function
    try:
        args.func(args)
    except AttributeError as e:
        parser.print_help()


if __name__ == '__main__':
    main()
