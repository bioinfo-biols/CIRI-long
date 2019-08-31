#!/usr/bin/env python
import gzip

from bwapy import BwaAligner
# index = '/home/zhangjy/git/CIRI-PacBio/test_data/bwa/mm10_genome.fa'
# bwa_aligner = BwaAligner(index, options='-T 19')

import mappy as mp
# mmi = '/home/zhangjy/git/CIRI-PacBio/test_data/minimap/mm10_genome.mmi'
# mp_aligner = mp.Aligner(fn_idx_in=mmi, preset='splice')

from skbio import DNA
from skbio.alignment import StripedSmithWaterman, local_pairwise_align_ssw

from CIRI import poagraph
from CIRI import seqgraphalignment
from CIRI.utils import log
from CIRI.utils import to_str, to_bytes


class Query(object):
    def __init__(self, header='', sequence=''):
        self.header = header
        self.sequence = sequence


class Trimmed(object):
    def __init__(self, header='', sequence='', subseq=''):
        self.header = header
        self.sequence = sequence
        self.subseq = subseq


def read_fasta(infile):
    query = None
    with gzip.open(infile, 'rb') as f:
        for line in f:
            content = to_str(line.rstrip())
            if content.startswith('>'):
                if query:
                    yield query
                query = Query(header=content.lstrip('>'), sequence='')
            else:
                query.sequence += content
        if query:
            yield query


def circularize(infile):
    log('Input: {}'.format(infile))

    P5_primer = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC.T'
    P7_primer = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG'
    P5_query = StripedSmithWaterman(P5_primer)
    P7_query = StripedSmithWaterman(P7_primer)

    ccs_reads = []
    cnt = 0
    for query in read_fasta(infile):
        cnt += 1
        if cnt / 1000 == cnt // 1000:
            log('Loaded {} CCS reads'.format(cnt))

        if len(query.sequence) < 300:
            continue

        # Remove adapter sequence
        P5_alignment = P5_query(query.sequence)
        P7_alignment = P7_query(query.sequence)
        rev_P5_alignment = P5_query(mp.revcomp(query.sequence))
        rev_P7_alignment = P7_query(mp.revcomp(query.sequence))

        if P5_alignment.query_end - P5_alignment.query_begin >= len(P5_primer) - 10 \
           and P7_alignment.query_end - P7_alignment.query_begin >= len(P7_primer) - 10:
            # Forward strand
            trimmed_seq = query.sequence[len(P5_primer):-len(P7_primer)]
        elif rev_P5_alignment.query_end - rev_P5_alignment.query_begin >= len(P5_primer) - 10 \
           and rev_P7_alignment.query_end - rev_P7_alignment.query_begin >= len(P7_primer) - 10:
            # Reverse strand
            trimmed_seq = query.sequence[len(P7_primer):-len(P5_primer)]
        else:
            continue

        sub_seq = get_subsequence(trimmed_seq)
        if sub_seq:
            ccs_reads.append(Trimmed(query.header, trimmed_seq, sub_seq))
        else:
            ccs_reads.append(Trimmed(query.header, trimmed_seq, ''))

    return ccs_reads


def get_subsequence(sequence):
    # Seed matching
    junc = 100

    while True:
        alignments, score, positions = local_pairwise_align_ssw(DNA(sequence[:junc]), DNA(sequence[junc:]))
        (l_st, l_en), (r_st, r_en) = positions
        if r_st == 0:
            break
        else:
            junc += r_st

    # keep reads containing 5' overlap
    if l_st <= 10 and len(sequence) - 10 <= junc + r_en and l_en - l_st >= 19 and r_en - r_st >= 19:
        return sequence[l_st:junc]
    else:
        return None


def get_consensus(fasta):
    match = 1
    gap = -3
    mismatch = -1
    globalAlign = 0
    simple = 0

    graph = poagraph.POAGraph(fasta[0][1], fasta[0][0])
    for label, sequence in fasta[1:]:
        alignment = seqgraphalignment.SeqGraphAlignment(sequence, graph, fastMethod=not simple,
                                                        globalAlign=globalAlign,
                                                        matchscore=match, mismatchscore=mismatch,
                                                        gapscore=gap)
        graph.incorporateSeqAlignment(alignment, sequence, label)

    consensus = graph.allConsenses()
    return ''.join(consensus[1])


if __name__ == '__main__':
    print()
