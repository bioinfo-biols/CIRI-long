#!/usr/bin/env python
import os
import re
import sys
import gzip
import itertools
import operator

from datetime import datetime
from collections import defaultdict
from collections import namedtuple

import mappy as mp
from CIRI.utils import log

Hit = namedtuple('Hit', 'q_st q_en r_st r_en')


def alignment_to_hit(align):
    return Hit(align.q_st, align.q_en, align.r_st, align.r_en)


def hit_len(x):
    return x.q_en - x.q_st


def merge_alignment(a, b):
    q_st = min(a.q_st, b.q_st)
    q_en = max(a.q_en, b.q_en)
    r_st = min(a.r_st, b.r_st)
    r_en = max(a.r_en, b.r_en)
    return Hit(q_st=q_st, q_en=q_en, r_st=r_st, r_en=r_en)


def print_hit(query):
    for i in query:
        print(i)
    return 1


def main():
    # aligner = mp.Aligner(fn_idx_in='/home/zhangjy/git/CIRI/test_data/chr1.mmi', preset='splice', )
    aligner = mp.Aligner(fn_idx_in='/home/zhangjy/git/CIRI/test_data/28NT_pseudo.mmi', preset='splice')
    log('Loaded Minimap2 aligner')

    fasta = '/home/zhangjy/git/CIRI/test_data/m54148_190408_064610.ccs.fasta.gz'
    # fasta = '/home/zhangjy/git/CIRI/test_data/BSJ.fa.gz'

    ccs_reads = {}
    with gzip.open('/home/zhangjy/git/CIRI/test_data/circ_ccs.fasta.gz', 'wb') as out:
        for query in read_fasta(fasta):
            ccs_reads[query.header] = 0

            if len(ccs_reads) / 500 == len(ccs_reads) // 500:
                log('Loaded {} CCS reads'.format(len(ccs_reads)))

            match = empty_iter(aligner.map(query.sequence, cs=True))
            if match is None:
                continue

            # traverse hits
            tmp_hit = defaultdict(dict)
            for hit in match:
                tmp_hit[hit.ctg].setdefault(hit.strand, []).append(hit)

            for contig, contig_hit in tmp_hit.items():
                for strand, query_hit in contig_hit.items():
                    # skip uniquely aligned reads
                    # TODO: de novo align unaligned query sequence

                    if len(query_hit) == 1:
                        continue

                    query_hit.sort(key=lambda x: x.q_en * strand)
                    query_hit.sort(key=lambda x: x.q_st * strand)

                    sorted_hit = []
                    tmp_hit = alignment_to_hit(query_hit[0])
                    for _, hit in enumerate(query_hit, 2):
                        overlap = min(hit.q_en, tmp_hit.q_en) - max(hit.q_st, tmp_hit.q_st)
                        if overlap < max(hit_len(hit), hit_len(tmp_hit)) * 0.4:
                            sorted_hit.append(tmp_hit)
                            tmp_hit = alignment_to_hit(hit)
                        else:
                            tmp_hit = merge_alignment(tmp_hit, hit)

                    sorted_hit.append(tmp_hit)

                    r_hit = sorted_hit.copy()
                    r_hit.sort(key=lambda x: x.r_en)
                    r_hit.sort(key=lambda x: x.r_st)

                    # Drop reads mapped in correct order
                    if operator.eq(sorted_hit, r_hit):
                        continue

                    ccs_reads[query.header] = 1

                    r_sect = ['{}-{}'.format(i.r_st, i.r_en) for i in query_hit]
                    q_sect = ['{}-{}'.format(i.q_st, i.q_en) for i in query_hit]
                    tmp_header = to_bytes('>{}\t{}\t{}\t{}\t{}\n'.format(
                        query.header, strand, query_hit[0].ctg, ','.join(q_sect), ','.join(r_sect)))
                    out.write(tmp_header)

                    tmp_seq = to_bytes('{}\n'.format(query.sequence))
                    out.write(tmp_seq)

            # tmp_hit = defaultdict(dict)
            # for hit in match:
            #     tmp_hit[hit.ctg].setdefault(hit.strand, []).append(hit)
            #
            # for contig, contig_hit in tmp_hit.items():
            #     for strand, query_hit in contig_hit.items():
            #         if len(query_hit) == 1:
            #             continue
            #
            #         query_hit.sort(key=lambda x: x.q_en * strand)
            #         query_hit.sort(key=lambda x: x.q_st * strand)
            #
            #         if len(query_hit) == 2:
            #             # different query segment
            #             if strand > 0:
            #                 gap = query_hit[1].q_st - query_hit[0].q_en
            #             else:
            #                 gap = query_hit[1].q_en - query_hit[0].q_st
            #
            #             if gap < 0 and abs(gap) >= min(query_hit[0].blen, query_hit[1].blen) * 0.5:
            #                 continue
            #
            #             # on same strand
            #             if query_hit[0].strand != query_hit[1].strand:
            #                 continue
            #
            #             # Back spliced junction
            #             if query_hit[1].r_st - query_hit[0].r_en >= gap + 10:
            #                 continue
            #
            #             print_hit(query_hit)
            #             r_sect = ['{}-{}'.format(i.r_st, i.r_en) for i in query_hit]
            #             q_sect = ['{}-{}'.format(i.q_st, i.q_en) for i in query_hit]
            #             tmp_header = to_bytes('>{}\t{}\t{}\t{}\n'.format(
            #                 query_name, query_hit[0].ctg, ','.join(q_sect), ','.join(r_sect)))
            #             out.write(tmp_header)
            #
            #             tmp_seq = to_bytes('{}\n'.format(query_sequence))
            #             out.write(tmp_seq)
            #
            #             circ_reads += 1
            #             total_ccs[ccs_id] += 1

    log('Total CCS Reads: {}'.format(len(ccs_reads)))
    log('circRNA CCS Reads: {}'.format(len([i for i, j in ccs_reads.items() if j == 1])))


if __name__ == '__main__':
    main()
