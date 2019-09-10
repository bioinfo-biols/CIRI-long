#!/home/zhangjy/.virtualenvs/Benchmarking/bin/python
import os
import sys
import argparse


def main():
    from CIRI.version import __version__
    from CIRI.rofinder import find_ccs_reads
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
    parser.add_argument('--debug', dest='debug', default=False, action='store_true',
                        help='Run in debuggin mode', )

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()

    if args.input is None or args.output is None:
        sys.exit('Please specific input and output files!')

    logger = get_logger('CIRI-long')
    logger.info('Input reads: ' + os.path.basename(args.input))
    logger.info('Output CCS reads: ' + os.path.basename(args.output))
    logger.info('Multi threads: {}'.format(args.threads))

    if args.input is None or args.output is None:
        sys.exit('Please provide input and output file, run CIRI-long using -h or --help for detailed information.')

    in_file = check_file(args.input)
    out_dir = check_dir(args.output)
    check_dir(out_dir + '/tmp')
    prefix = args.prefix
    threads = int(args.threads)
    debugging = args.debug

    total_reads, ro_reads = find_ccs_reads(in_file, out_dir, prefix, threads, debugging)

    logger.info('Total Reads: {}, RO reads: {}'.format(total_reads, ro_reads))


if __name__ == '__main__':
    main()
