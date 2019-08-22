#!/home/zhangjy/.virtualenvs/Benchmarking/bin/python
import argparse
import sys

from CIRI import assemble


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    # parser.add_argument('-G', '--gap', type=int, default=-2, help='Gap penalty, default=-1')
    # parser.add_argument('-g', '--globalAlign', action='store_true', help='Global alignment (default: local)')
    args = parser.parse_args()

    fasta = args.infile

    fasta = '/home/zhangjy/Data/PROJECT/01.Benchmarking_Study/pacbio/m54148_190411_071934.ccs.fasta.gz'
    # fasta = '/home/zhangjy/git/CIRI-PacBio/test_data/m54148_190408_064610.ccs.fasta.gz'
    # fasta = '/home/zhangjy/git/CIRI-PacBio/test_data/BSJ.fa.gz'

    ccs_reads = assemble.circularize(fasta)
    print(len(ccs_reads))


if __name__ == '__main__':
    main()
