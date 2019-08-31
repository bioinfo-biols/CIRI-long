#!/home/zhangjy/.virtualenvs/Benchmarking/bin/python
import os
import sys
import argparse
from logger import get_logger


def main():
    from rofinder import find_ccs_reads
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in', dest='input', metavar='READS', required=True,
                        help='Input reads.fq.gz', )
    parser.add_argument('-o', '--out', dest='output', metavar='CCS', required=True,
                        help='Output reads info', )
    parser.add_argument('-t', '--threads', dest='threads', metavar='INT', default=os.cpu_count(),
                        help='Number of threads', )
    args = parser.parse_args()

    logger = get_logger('CIRI-long')
    logger.info('Input reads: ' + os.path.basename(args.input))
    logger.info('Output CCS reads: ' + os.path.basename(args.output))
    logger.info('Multi threads: {}'.format(args.threads))

    in_file = args.input
    out_file = args.output
    threads = int(args.threads)

    total_reads, ro_reads = find_ccs_reads(in_file, out_file, threads)
    logger.info('Total Reads: {}, RO reads: {}'.format(total_reads, ro_reads))


if __name__ == '__main__':
    main()
