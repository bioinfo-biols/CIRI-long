import logging
LOGGER = logging.getLogger('CIRI-long')

import numpy as np
import mappy as mp

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
    consensus_cnt = 0
    raw_unmapped_cnt = 0
    ccs_mapped_cnt = 0
    accordance_cnt = 0
    bsj_cnt = 0

    ret = []
    for read_id, segments, ccs, raw in chunk:
        # Filter 1 - Remove ccs which not that consensus
        fasta = []
        for pos in segments.split(';'):
            i_st, i_en = pos.split('-')
            fasta.append(('{}-{}'.format(i_st, i_en), raw[int(i_st):int(i_en)]))

        d_mean = np.mean([len(i) for _, i in fasta])
        d_delta = int(2.3 * np.sqrt(0.1 * d_mean))

        fasta_ccs = msa(fasta)
        if len(fasta_ccs) > 1:
            continue
        if 0.9 * (d_mean - d_delta) < len(fasta_ccs[0]) < 1.1 * (d_mean + d_delta):
            continue
        consensus_cnt += 1

        # Filter 2 - Remove linear mapped reads
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
        raw_unmapped_cnt += 1

        if hit is not None:
            raw_en, raw_st = hit.q_en, hit.q_st
        else:
            raw_en, raw_st = None, None

        # Filter 3 - Remove reads that ccs doesn't mapped to genome
        ccs_hit = None
        for hit in ALIGNER.map(ccs):
            if hit.is_primary:
                ccs_hit = hit
        if ccs_hit is None or ccs_hit.ctg in ['MT']:
            continue
        ccs_mapped_cnt += 1

        # Filter 4 - check the accordance alignment results of non-repeat segments and ccs
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

        # Filter 5 - Find BSJ site in CCS
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
        bsj_cnt += 1

    return consensus_cnt, raw_unmapped_cnt, ccs_mapped_cnt, accordance_cnt, ret


ALIGNER = None


def initializer(aligner):
    global ALIGNER
    ALIGNER = aligner


def filter_ccs_reads(ccs_seq, minimap_index, out_dir, prefix, threads, debugging):
    from collections import namedtuple
    from multiprocessing import Pool
    Read = namedtuple('Read', 'segments ccs raw')

    minimap_aligner = mp.Aligner(minimap_index)

    reads_count = {
        'consensus': 0,
        'raw_unmapped': 0,
        'ccs_mapped': 0,
        'accordance': 0,
        'bsj': 0,
    }

    chunk_size = 250
    chunk_cnt = 0
    jobs = []
    pool = Pool(threads, initializer, (minimap_aligner,))
    for reads in grouper(list(ccs_seq), chunk_size):
        chunk = [[i, ] + ccs_seq[i] for i in reads if i is not None]
        chunk_cnt += 1
        jobs.append(pool.apply_async(parse_chunk, (chunk,)))
    pool.close()

    prog = ProgressBar()
    finished_cnt = 0
    with open('{}/{}.cand_circ.fa'.format(out_dir, prefix), 'w') as out:
        for job in jobs:
            finished_cnt += 1

            tmp_consensus, tmp_raw_unmapped, tmp_ccs_mapped, tmp_accordance, ret = job.get()
            reads_count['consensus'] += tmp_consensus
            reads_count['raw_unmapped'] += tmp_raw_unmapped
            reads_count['ccs_mapped'] += tmp_ccs_mapped
            reads_count['accordance'] += tmp_accordance
            reads_count['bsj'] += len(ret)

            for read_id, segments, tmp_junc, tmp_circ, tmp in ret:
                out.write('>{}\t{}\t{}\t{}\n{}\n'.format(read_id, segments, tmp_junc, tmp_circ, tmp))
            prog.update(100 * finished_cnt / chunk_cnt)
    pool.join()
    prog.update(100)

    return reads_count

