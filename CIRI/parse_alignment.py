import os
import sys
import pysam
import mappy as mp
import time
from CIRI import poagraph
from CIRI import seqgraphalignment
from CIRI.logger import ProgressBar


def msa(fasta, verbose=0):
    match = 1
    gap = -2
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
    if verbose:
        alignments = graph.generateAlignmentStrings()
        for label, alignstring in alignments:
            print("{0:15s} {1:s}".format(label, alignstring))
        print()
    return [''.join(i[1]) for i in graph.allConsenses()]


def get_blocks(hit):
    r_start = hit.r_st
    r_end = hit.r_st
    r_block = []
    for length, operation in hit.cigar:
        if operation == 0:
            r_end += length
        elif operation == 1:
            pass
        elif operation == 2:
            r_end += length
        elif operation == 3:
            r_block.append((r_start, r_end, r_end - r_start + 1))
            r_start = r_end + length
            r_end = r_start
        elif operation == 4:
            pass
    if r_end > r_start:
        r_block.append((r_start, r_end, r_end - r_start + 1))
    return ','.join(['{}-{}|{}'.format(i[0], i[1], i[2]) for i in r_block])


def grouper(iterable, n, fillvalue=None):
    from itertools import zip_longest
    """
    Collect data info fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    """

    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=None)


def parse_chunk(chunk):
    total_cnt = 0
    unmapped_cnt = 0
    accordance_cnt = 0
    circ_cnt = 0
    ret = []
    for read_id, segments, ccs, raw in chunk:
        total_cnt += 1

        # Remove linear mapped reads
        is_raw_mapped = 0
        hit = None
        for hit in ALIGNER.map(raw):
            if not hit.is_primary:
                continue

            # Filter out linear reads
            if hit.q_en - hit.q_st > min(len(raw) * 0.8, len(raw) - 50):
                is_raw_mapped = 1
        if is_raw_mapped:
            continue
        unmapped_cnt += 1

        if hit is not None:
            raw_en, raw_st = hit.q_en, hit.q_st
        else:
            raw_en, raw_st = None, None

        # Remove reads that ccs doesn't mapped to genome
        ccs_hit = None
        for hit in ALIGNER.map(ccs):
            if hit.is_primary:
                ccs_hit = hit
        if ccs_hit is None or ccs_hit.ctg in ['MT']:
            continue

        # check the alignment results of non-repeat segments
        seg_st = int(segments.split(';')[0].split('-')[0])
        seg_en = int(segments.split(';')[-1].split('-')[1])

        seg_st_hit = None
        for hit in ALIGNER.map(raw[:seg_st]):
            if hit.is_primary:
                seg_st_hit = hit
        is_st_mapped = 1 if seg_st_hit is not None and seg_st_hit.mlen > seg_st * 0.6 else 0

        seg_en_hit = None
        for hit in ALIGNER.map(raw[seg_en:]):
            if hit.is_primary:
                seg_en_hit = hit
        is_en_mapped = 1 if seg_en_hit is not None and seg_en_hit.mlen > (len(raw) - seg_en) * 0.6 else 0

        is_accordance = 1
        if is_st_mapped and seg_st_hit.ctg != ccs_hit.ctg:
            is_accordance = 0
        if is_en_mapped and seg_en_hit.ctg != ccs_hit.ctg:
            is_accordance = 0
        if is_accordance == 0:
            continue

        # Remove other mapped region that intersect with ccs
        if raw_st is None or raw_en is None:
            pass
        elif raw_st + 50 < seg_st or seg_en + 50 < raw_en:
            continue
        elif raw_en - raw_st > 1.5 * len(ccs):
            continue
        else:
            pass

        # Remove CCS that not start from tailed end
        if seg_st > 25 or seg_en < len(raw) - 25:
            continue

        accordance_cnt += 1

        # Remove ccs which not that consensus
        fasta = []
        for pos in segments.split(';'):
            i_st, i_en = pos.split('-')
            fasta.append(('{}-{}'.format(i_st, i_en), raw[int(i_st):int(i_en)]))

        fasta_ccs = msa(fasta)
        if len(fasta_ccs) > 1:
            continue

        # Find junction site in CCS
        tmp = ccs
        tmp_junc = 0
        iter_num = 0
        while True:
            iter_num += 1
            tmp_hits = [hit for hit in ALIGNER.map(tmp) if hit.is_primary]
            if len(tmp_hits) == 0:
                break
            tmp_hit = sorted(tmp_hits, key=lambda x: x.mlen, reverse=True)[0]

            if tmp_hit.q_st < 10 < len(tmp) - tmp_hit.q_en:
                tmp_junc = (tmp_junc + tmp_hit.q_en) % len(tmp)
                tmp = ccs[tmp_junc:] + ccs[:tmp_junc]
            elif tmp_hit.q_st > 10 > len(tmp) - tmp_hit.q_en:
                tmp_junc = (tmp_junc + tmp_hit.q_st) % len(tmp)
                tmp = ccs[tmp_junc:] + ccs[:tmp_junc]
            else:
                break
            if iter_num >= 10:
                break

        if iter_num >= 10:
            continue

        tmp_hits = []
        for hit in ALIGNER.map(tmp):
            if hit.is_primary:
                tmp_hits.append(hit)
        if len(tmp_hits) == 0:
            continue
        tmp_hit = sorted(tmp_hits, key=lambda x: x.mlen, reverse=True)[0]
        if tmp_hit.r_en - tmp_hit.r_st < 0.9 * len(ccs):
            continue

        # Retrive circRNA positions
        ret.append((read_id, segments, tmp_junc, '{}:{}-{}'.format(tmp_hit.ctg, tmp_hit.r_st, tmp_hit.r_en - 1), tmp))
        circ_cnt += 1

    return total_cnt, unmapped_cnt, accordance_cnt, circ_cnt, ret


ALIGNER = None


def initializer(aligner):
    global ALIGNER
    ALIGNER = aligner


def parse_sample(sample, aligner):
    from collections import namedtuple
    from multiprocessing import Pool
    Read = namedtuple('Read', 'segments ccs raw')

    all_reads = {}
    with open('/home/zhangjy/03.CIRIpacbio/1908_Nanopore/CIRI-long/results/{}.ccs.fa'.format(sample), 'r') as f:
        for line in f:
            header = line.rstrip()
            content = header.split('\t')
            seq = f.readline().rstrip()
            all_reads[content[0].lstrip('>')] = [content[1], seq]

    with open('/home/zhangjy/03.CIRIpacbio/1908_Nanopore/CIRI-long/results/{}.trimmed.seq'.format(sample), 'r') as f:
        for line in f:
            header = line.rstrip().split('\t')[0].lstrip('>')
            seq = f.readline().rstrip()
            all_reads[header].append(seq)

    #     with open('../{}.cand_circ.fa'.format(sample), 'w') as out:
    total_cnt = 0
    unmapped_cnt = 0
    accordance_cnt = 0
    circ_cnt = 0
    chunk_size = 250
    chunk_cnt = 0
    jobs = []
    pool = Pool(12, initializer, (aligner,))
    for reads in grouper(list(all_reads), chunk_size):
        chunk = [[i, ] + all_reads[i] for i in reads if i is not None]
        chunk_cnt += 1
        jobs.append(pool.apply_async(parse_chunk, (chunk,)))
    pool.close()

    prog = ProgressBar()
    finished_cnt = 0
    with open('/home/zhangjy/git/CIRI-long/test_data/{}.cand_circ.fa'.format(sample), 'w') as out:
        for job in jobs:
            finished_cnt += 1
            tmp_total, tmp_unmapped, tmp_accordance, tmp_circ, ret = job.get()
            total_cnt += tmp_total
            unmapped_cnt += tmp_unmapped
            accordance_cnt += tmp_accordance
            circ_cnt += tmp_circ

            for read_id, segments, tmp_junc, tmp_circ, tmp in ret:
                out.write('>{}\t{}\t{}\t{}\n{}\n'.format(read_id, segments, tmp_junc, tmp_circ, tmp))
            prog.update(100 * finished_cnt / chunk_cnt)
    pool.join()
    prog.update(100)

    return total_cnt, unmapped_cnt, accordance_cnt, circ_cnt


def main():
    genome_aligner = mp.Aligner('/home/zhangjy/03.CIRIpacbio/1908_Nanopore/CIRI-long/results/genome_splice.mmi')

    print('Sample', 'Total', 'Unmapped', 'Accordance', 'Circ')
    for sample in 'barcode20', 'barcode21', 'barcode22':
        total_cnt, unmapped_cnt, accordance_cnt, circ_cnt = parse_sample(sample, genome_aligner)
        print(sample, total_cnt, unmapped_cnt, accordance_cnt, circ_cnt)


if __name__ == '__main__':
    main()
