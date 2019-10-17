#!/home/zhangjy/.virtualenvs/Benchmarking/bin/python
import os
import sys
import argparse


def main():
    from CIRI.version import __version__
    from CIRI.rofinder import find_ccs_reads, load_ccs_reads
    from CIRI.alignment import filter_ccs_reads
    from CIRI.utils import check_file, check_dir
    from CIRI.logger import get_logger

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in', dest='input', metavar='READS', default=None,
                        help='Input reads.fq.gz', )
    parser.add_argument('-o', '--out', dest='output', metavar='DIR', default=None,
                        help='Output directory, default: ./', )
    parser.add_argument('-p', '--prefix', dest='prefix', metavar='PREFIX', default="CIRI-long",
                        help='Output sample prefix, default: CIRI-long', )
    parser.add_argument('-t', '--threads', dest='threads', metavar='INT', default=os.cpu_count(),
                        help='Number of threads', )
    parser.add_argument('-x', '--mmi', dest='mmi', metavar='MMI', default=None,
                        help='Minimap2 index of reference genome', )
    parser.add_argument('--debug', dest='debug', default=False, action='store_true',
                        help='Run in debuggin mode', )

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()

    if args.input is None or args.output is None:
        sys.exit('Please provide input and output file, run CIRI-long using -h or --help for detailed information.')
    if args.mmi is None:
        sys.exit('Please specific minimap2 index or reference fasta')

    in_file = check_file(args.input)
    out_dir = check_dir(args.output)
    minimap_index = check_file(args.mmi)
    check_dir(out_dir + '/tmp')
    prefix = args.prefix
    threads = int(args.threads)
    debugging = args.debug

    logger = get_logger('CIRI-long', fname='{}/{}.log'.format(out_dir, prefix), verbosity=debugging)
    logger.info('Input reads: ' + os.path.basename(in_file))
    logger.info('Output directory: ' + os.path.basename(out_dir))
    logger.info('Multi threads: {}'.format(args.threads))

    if os.path.exists('{}/{}.ccs.fa'.format(out_dir, prefix)) and os.path.exists('{}/{}.trimmed.seq'.format(out_dir, prefix)):
        logger.info('Step 1 - Loading circRNA candidates in previous run')
        ccs_seq = load_ccs_reads(out_dir, prefix)
        reads_count = {
            'Total': 'Unknown',
            'with_repeats': len(ccs_seq)
        }
    else:
        logger.info('Step 1 - Scanning raw reads to find circRNA candidates')
        total_reads, ro_reads, ccs_seq = find_ccs_reads(in_file, out_dir, prefix, threads, debugging)
        reads_count = {
            'Total': total_reads,
            'with_repeats': ro_reads,
        }

    logger.info('Step 2 - Candidate reads alignment & filter')
    tmp_cnt = filter_ccs_reads(ccs_seq, minimap_index, out_dir, prefix, threads, debugging)
    reads_count.update(tmp_cnt)

    logger.info('All Finished!')
    logger.debug('Summary:')
    logger.debug('Total Reads: {}'.format(reads_count['Total']))
    logger.debug('Cyclic Tandem Repeats: {}'.format(reads_count['with_repeats']))
    logger.debug('Consensus: {}'.format(reads_count['consensus']))
    logger.debug('Raw unmapped: {}'.format(reads_count['raw_unmapped']))
    logger.debug('CCS mapped: {}'.format(reads_count['ccs_mapped']))
    logger.debug('CCS & Raw aligned concordantly: {}'.format(reads_count['accordance']))
    logger.debug('BSJ: {}'.format(reads_count['bsj']))


if __name__ == '__main__':
    main()
