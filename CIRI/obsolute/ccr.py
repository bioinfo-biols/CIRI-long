#!/usr/bin/env python
import os
import sys
import pysam
import numpy as np
import mappy as mp

READ_ID = None
READ_SEQ = ''
CIRC_ID = None
CIRC_SEQUENCE = ''
SAMPLE = None


def load_sequence():
    global CIRC_SEQUENCE, READ_SEQ
    with open('/home/zhangjy/03.CIRIpacbio/1908_Nanopore/alignment/circRNA_pseudo.fa', 'r') as f:
        for line in f:
            header = line.rstrip().lstrip('>')
            seq = f.readline().rstrip()
            if header == CIRC_ID:
                CIRC_SEQUENCE = seq[:len(seq) // 2]
                break

    read = None
    sam = pysam.AlignmentFile('/home/zhangjy/03.CIRIpacbio/1908_Nanopore/alignment/{}_circRNA_pseudo.sorted.bam'.format(SAMPLE), 'rb')
    for x in sam.fetch(CIRC_ID, multiple_iterators=True):
        if x.query_name != READ_ID:
            continue
        if x.is_supplementary or x.is_secondary:
            continue
        read = x
        break

    assert read is not None
    READ_SEQ = read.query_sequence


def trim_sequence(seq, primer = 'AAGCAGTGGTATCAACGCAGAGTAC'):
    from skbio import DNA
    from skbio.alignment import local_pairwise_align_ssw
    l_alignments, l_score, l_positions = local_pairwise_align_ssw(DNA(seq[:100]), DNA(primer))
    l_trim = l_positions[0][1] + 1 if l_score >= 25 else 0

    r_alignments, r_score, r_positions = local_pairwise_align_ssw(DNA(seq[-100:]), DNA(mp.revcomp(primer)))
    r_trim = -100 + r_positions[0][0] if r_score >= 25 else -1
    return seq[l_trim:r_trim]


def find_consensus(seq):
    score_mtx = np.ndarray(shape=(len(seq), len(seq)))
    score_mtx[0] = 1


def main():
    global READ_ID
    global CIRC_ID
    global SAMPLE
    READ_ID = 'e16e46cb-7063-4353-83b4-10dab4a5625c'
    CIRC_ID = '11:94097218-94099152'
    SAMPLE = 'barcode20'

    load_sequence()
    print('Read ID: {}'.format(READ_ID))
    print('Read Length: {}'.format(len(READ_SEQ)))
    print('CIRC_ID: {}'.format(CIRC_ID))
    print('circRNA Length: {}'.format(len(CIRC_SEQUENCE)))

    seq = trim_sequence(READ_SEQ, 'AAGCAGTGGTATCAACGCAGAGTAC')
    print('Trimmed length: {}'.format(len(READ_SEQ) - len(seq)))

    ccs = find_consensus(seq[:20])


if __name__ == '__main__':
    main()
