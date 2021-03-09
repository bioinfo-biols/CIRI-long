import os
import sys
import re
import logging
from multiprocessing import Pool
from collections import defaultdict

import pysam
from CIRI_long import env
from CIRI_long.align import *
from CIRI_long.logger import ProgressBar
from CIRI_long.utils import grouper, revcomp

LOGGER = logging.getLogger('CIRI-long')


def search_splice_signal(contig, start, end, clip_base, search_length=10, shift_threshold=3):
    # Find free sliding region
    # start | real_start <-> end | real_end
    ds_free = 0
    for i in range(100):
        if end + i > env.CONTIG_LEN[contig]:
            break
        if env.GENOME.seq(contig, start, start + i) == env.GENOME.seq(contig, end, end + i):
            ds_free = i
        else:
            break

    us_free = 0
    for j in range(100):
        if start - j < 0:
            break
        if env.GENOME.seq(contig, start - j, start) == env.GENOME.seq(contig, end - j, end):
            us_free = j
        else:
            break

    if start - search_length - us_free - 2 < 0 or end + search_length + ds_free + 2 > env.CONTIG_LEN[contig]:
        return None, us_free, ds_free
    # Splice site: site_id, strand, us_shift, ds_shift, site_weight, altered_len, altered_total
    # First: Find flanking junction from annotation gtf
    if env.SS_INDEX is not None:
        anno_ss = []
        for strand in ['+', '-']:
            tmp_us_sites = []
            for us_shift in range(-search_length, search_length):
                us_pos = start + us_shift
                if contig in env.SS_INDEX and us_pos in env.SS_INDEX[contig] and strand in env.SS_INDEX[contig][us_pos]:
                    tmp_us_sites.append(us_shift - 1)

            tmp_ds_sites = []
            for ds_shift in range(-search_length, search_length):
                ds_pos = end + ds_shift
                if contig in env.SS_INDEX and ds_pos in env.SS_INDEX[contig] and strand in env.SS_INDEX[contig][ds_pos]:
                    tmp_ds_sites.append(ds_shift)

            if len(tmp_us_sites) == 0 or len(tmp_ds_sites) == 0:
                continue

            for i in tmp_us_sites:
                for j in tmp_ds_sites:
                    if abs(i - j) > shift_threshold + clip_base:
                        continue
                    us_ss = env.GENOME.seq(contig, start + i - 2, start + i)
                    ds_ss = env.GENOME.seq(contig, end + j, end + j + 2)
                    if strand == '-':
                        us_ss, ds_ss = revcomp(ds_ss), revcomp(us_ss)
                    ss_id = '{}-{}|{}-{}'.format(us_ss, ds_ss, i, j)
                    ss_weight = SPLICE_SIGNAL[(ds_ss, us_ss)] if (ds_ss, us_ss) in SPLICE_SIGNAL else 3

                    anno_ss.append((
                        ss_id, strand, i, j, ss_weight, *get_ss_altered_length(i, j, us_free, ds_free, clip_base)
                    ))

        if len(anno_ss) > 0:
            return sort_ss(anno_ss, us_free, ds_free, clip_base), us_free, ds_free

    # Second: Find Denovo BSJ using pre-defined splice signal
    us_search_length = search_length + us_free
    ds_search_length = search_length + ds_free
    us_seq = env.GENOME.seq(contig, start - us_search_length - 2, start + ds_search_length)
    ds_seq = env.GENOME.seq(contig, end - us_search_length, end + ds_search_length + 2)

    if us_seq is None or len(us_seq) < ds_search_length - us_search_length + 2:
        return None, us_free, ds_free
    if ds_seq is None or len(ds_seq) < ds_search_length - us_search_length + 2:
        return None, us_free, ds_free

    putative_ss = []
    for strand in ['+', '-']:
        for (tmp_ds_ss, tmp_us_ss), ss_weight in SPLICE_SIGNAL.items():
            if strand == '-':
                ds_ss, us_ss = revcomp(tmp_us_ss), revcomp(tmp_ds_ss)
            else:
                ds_ss, us_ss = tmp_ds_ss, tmp_us_ss

            # Find upstream signal
            tmp_us_start = 0
            tmp_us_sites = []
            while 1:
                tmp_us = us_seq.find(us_ss, tmp_us_start + 1)
                if tmp_us == -1:
                    break
                tmp_us_sites.append(tmp_us)
                tmp_us_start = tmp_us

            # Find downstream signal
            tmp_ds_start = 0
            tmp_ds_sites = []
            while 1:
                tmp_ds = ds_seq.find(ds_ss, tmp_ds_start + 1)
                if tmp_ds == -1:
                    break
                tmp_ds_sites.append(tmp_ds)
                tmp_ds_start = tmp_ds

            # Filter paired splice signal in concordance position
            if len(tmp_us_sites) == 0 or len(tmp_ds_sites) == 0:
                continue

            for i in tmp_us_sites:
                for j in tmp_ds_sites:
                    if abs(i - j) > clip_base + shift_threshold:
                        continue
                    us_shift = i - us_search_length
                    ds_shift = j - us_search_length
                    ss_id = '{}-{}*|{}-{}'.format(tmp_us_ss, tmp_ds_ss, us_shift, ds_shift)
                    putative_ss.append((
                        ss_id, strand, us_shift, ds_shift, ss_weight,
                        *get_ss_altered_length(us_shift, ds_shift, us_free, ds_free, clip_base)
                    ))

    if len(putative_ss) > 0:
        return sort_ss(putative_ss, us_free, ds_free, clip_base), us_free, ds_free

    return None, us_free, ds_free


def find_bsj(ccs):
    """
    Find junction using aligner
    """
    init_hit = get_primary_alignment(env.ALIGNER.map(ccs * 2))
    if init_hit is None:
        return None, None

    circ_junc = init_hit.q_st % len(ccs)
    circ = ccs[circ_junc:] + ccs[:circ_junc]

    last_junc = 0
    last_m = 0

    itered_junc = {}
    while True:
        circ_hit = get_primary_alignment(env.ALIGNER.map(circ))
        if circ_hit is None or circ_hit.mlen <= last_m:
            circ_junc = last_junc
            break
        last_m = circ_hit.mlen
        last_junc = circ_junc

        st_clip, en_clip = circ_hit.q_st, len(circ) - circ_hit.q_en
        if st_clip == 0 and en_clip == 0:
            break

        if st_clip >= en_clip:
            circ_junc = (circ_junc + st_clip) % len(circ)
        else:
            circ_junc = (circ_junc + circ_hit.q_en) % len(circ)

        if circ_junc in itered_junc:
            circ_junc = last_junc
            break

        circ = ccs[circ_junc:] + ccs[:circ_junc]
        itered_junc[circ_junc] = 1

    circ = ccs[circ_junc:] + ccs[:circ_junc]
    return circ, circ_junc


def align_clip_segments(circ, hit):
    """
    Align clip bases
    """
    from libs.striped_smith_waterman.ssw_wrap import Aligner
    from collections import Counter
    st_clip, en_clip = hit.q_st, len(circ) - hit.q_en
    clip_r_st, clip_r_en, clipped_circ = None, None, None

    if st_clip + en_clip >= 20:
        clip_seq = circ[hit.q_en:] + circ[:hit.q_st]
        if len(clip_seq) > 0.6 * len(circ):
            return None, None, None, None

        tmp_start = max(hit.r_st - 200000, 0)
        tmp_end = min(hit.r_en + 200000, env.CONTIG_LEN[hit.ctg])

        tmp_seq = env.GENOME.seq(hit.ctg, tmp_start, tmp_end)
        if Counter(tmp_seq)['N'] >= 0.3 * (tmp_end - tmp_start):
            return None, None, None, None

        if hit.strand > 0:
            ssw = Aligner(tmp_seq, match=1, mismatch=1, gap_open=1, gap_extend=1)
            align_res = ssw.align(clip_seq)
            clip_r_st, clip_r_en = tmp_start + align_res.ref_begin, tmp_start + align_res.ref_end
            if clip_r_st < hit.r_st:
                clipped_circ = clip_seq[align_res.query_begin:] + \
                               circ[hit.q_st:hit.q_en] + \
                               clip_seq[:align_res.query_begin]
            else:
                clipped_circ = circ[hit.q_st:] + circ[:hit.q_st]
        else:
            ssw = Aligner(revcomp(tmp_seq), match=1, mismatch=1, gap_open=1, gap_extend=1)
            align_res = ssw.align(clip_seq)
            clip_r_st, clip_r_en = tmp_end - align_res.ref_end, tmp_end - align_res.ref_begin
            if clip_r_en > hit.r_en:
                clipped_circ = clip_seq[align_res.query_begin:] + \
                               circ[hit.q_st:hit.q_en] + \
                               clip_seq[:align_res.query_begin]
            else:
                clipped_circ = circ[hit.q_st:] + circ[:hit.q_st]

        clip_base = hit.q_st + len(circ) - hit.q_en - (align_res.query_end - align_res.query_begin) + 1
        circ_start = min(hit.r_st, clip_r_st) - 1
        circ_end = max(hit.r_en, clip_r_en)
    else:
        clipped_circ = circ[hit.q_st:] + circ[:hit.q_st]
        clip_base = st_clip + en_clip
        circ_start = hit.r_st - 1
        circ_end = hit.r_en

    return clipped_circ, circ_start, circ_end, (clip_r_st, clip_r_en, clip_base)


def scan_ccs_chunk(chunk, is_canonical):
    reads_cnt = defaultdict(int)
    ret = []

    short_reads = []
    for read_id, segments, ccs, raw in chunk:
        # Filter 1 - Remove linear mapped reads
        raw_hit = get_primary_alignment(env.ALIGNER.map(raw))
        if raw_hit and raw_hit.mlen > max(len(raw) * 0.8, len(raw) - 200):
            continue
        if raw_hit and raw_hit.mlen > 1.5 * len(ccs):
            continue

        raw_st = raw_hit.q_st if raw_hit else None
        raw_en = raw_hit.q_en if raw_hit else None
        reads_cnt['raw_unmapped'] += 1

        # Filter 2 - Remove other mapped region that intersect with ccs
        seg_st = int(segments.split(';')[0].split('-')[0])
        seg_en = int(segments.split(';')[-1].split('-')[1])
        if raw_hit and (raw_en < seg_st or raw_st > seg_en):
            continue

        ccs_hit = get_primary_alignment(env.ALIGNER.map(ccs * 2))
        if ccs_hit is None and len(ccs) < 150:
            short_reads.append((read_id, segments, ccs, raw))
        if ccs_hit is None or seg_en - seg_st < ccs_hit.q_en - ccs_hit.q_st:
            continue

        reads_cnt['ccs_mapped'] += 1

        # Find back-spliced junction site
        circ, junc = find_bsj(ccs)

        # Candidate alignment situation, more than 85%
        circ_hit = get_primary_alignment(env.ALIGNER.map(circ))
        if circ_hit is None or circ_hit.mlen < 0.75 * len(circ):
            continue

        clipped_circ, circ_start, circ_end, clip_info = align_clip_segments(circ, circ_hit)
        if circ_start is None or circ_end is None:
            continue

        clip_base = clip_info[2]
        if clip_base > 0.15 * len(ccs) or clip_base > 20:
            continue

        reads_cnt['bsj'] += 1

        # Retrive circRNA positions, convert minimap2 position to real position
        host_strand = find_host_gene(circ_hit.ctg, circ_start, circ_end)
        ss_site, us_free, ds_free, tmp_signal = find_annotated_signal(circ_hit.ctg, circ_start, circ_end, clip_base, clip_base + 10)
        if ss_site is None:
            ss_site = find_denovo_signal(circ_hit.ctg, circ_start, circ_end, host_strand, tmp_signal,
                                         us_free, ds_free, clip_base, clip_base + 10, 3, True)

        if ss_site is None:
            ss_id = 'NA'
            strand = 'NA'
            correction_shift = 0
        else:
            reads_cnt['signal'] += 1
            ss_id, strand, us_shift, ds_shift = ss_site
            circ_start += us_shift
            circ_end += ds_shift
            correction_shift = min(max(us_shift, us_free), ds_free)

        circ_id = '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end)

        # Get Cirexons
        cir_exons = get_blocks(circ_hit)
        cir_exons = merge_clip_exon(cir_exons, clip_info)

        cir_exons[0][0] = circ_start
        cir_exons[-1][1] = circ_end

        cir_exon_tag = []
        for cir_exon_start, cir_exon_end, cir_exon_length in cir_exons:
            cir_exon_tag.append('{}-{}|{}'.format(cir_exon_start + 1, cir_exon_end, cir_exon_length))

        # BSJ correction for 5' prime region
        circ_seq = clipped_circ if circ_hit.strand > 0 else revcomp(clipped_circ)
        circ_seq = circ_seq[correction_shift:] + circ_seq[:correction_shift]

        ret.append((
            read_id, circ_id, strand, ','.join(cir_exon_tag), ss_id,
            '{}|{}-{}'.format(junc, clip_base, len(circ)), segments, circ_seq
        ))

    return reads_cnt, short_reads, ret


def scan_ccs_reads(ccs_seq, ref_fasta, ss_index, gtf_index, intron_index, is_canonical, out_dir, prefix, threads):
    import mappy as mp

    faidx = Faidx(ref_fasta)
    contig_len = faidx.contig_len
    faidx.close()

    # First scanning using minimap2
    minimap_aligner = mp.Aligner(ref_fasta, n_threads=threads, preset='splice')

    chunk_size = 250
    jobs = []
    pool = Pool(threads, env.initializer,
                (minimap_aligner, contig_len, minimap_aligner, gtf_index, intron_index, ss_index))

    for reads in grouper(list(ccs_seq), chunk_size):
        chunk = [[i, ] + ccs_seq[i] for i in reads if i is not None]
        jobs.append(pool.apply_async(scan_ccs_chunk, (chunk, is_canonical)))
    pool.close()

    prog = ProgressBar()
    finished_cnt = 0

    reads_count = defaultdict(int)
    short_reads = []
    with open('{}/{}.cand_circ.fa'.format(out_dir, prefix), 'w') as out:
        for job in jobs:
            finished_cnt += 1

            tmp_cnt, tmp_short, ret = job.get()
            for key, value in tmp_cnt.items():
                reads_count[key] += value

            short_reads += tmp_short

            for read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, segments, circ_seq in ret:
                out.write('>{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\n'.format(
                    read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, segments, circ_seq
                ))
            prog.update(100 * finished_cnt / len(jobs))

    pool.join()
    prog.update(100)

    return reads_count, short_reads


def recover_ccs_chunk(chunk, is_canonical):
    reads_cnt = defaultdict(int)
    ret = []

    for read_id, segments, ccs, raw in chunk:
        # Remove other mapped region that intersect with ccs
        seg_st = int(segments.split(';')[0].split('-')[0])
        seg_en = int(segments.split(';')[-1].split('-')[1])

        ccs_hit = get_primary_alignment(env.ALIGNER.map(ccs * 2))
        if ccs_hit is None or seg_en - seg_st < ccs_hit.q_en - ccs_hit.q_st:
            continue

        reads_cnt['ccs_mapped'] += 1

        # Find back-spliced junction site
        circ, junc = find_bsj(ccs)

        # Candidate alignment situation, more than 85%
        circ_hit = get_primary_alignment(env.ALIGNER.map(circ))
        if circ_hit is None:
            continue

        clipped_circ, circ_start, circ_end, clip_info = align_clip_segments(circ, circ_hit)
        if circ_start is None or circ_end is None:
            continue

        clip_base = clip_info[2]
        if clip_base > 0.15 * len(ccs) or clip_base > 20:
            continue

        reads_cnt['bsj'] += 1

        # Retrive circRNA positions, convert minimap2 position to real position
        host_strand = find_host_gene(circ_hit.ctg, circ_start, circ_end)
        ss_site, us_free, ds_free, tmp_signal = find_annotated_signal(circ_hit.ctg, circ_start, circ_end, clip_base, clip_base + 10)
        if ss_site is None:
            ss_site = find_denovo_signal(circ_hit.ctg, circ_start, circ_end, host_strand, tmp_signal,
                                         us_free, ds_free, clip_base, clip_base + 10, 3, True)

        if ss_site is None:
            ss_id = 'NA'
            strand = 'NA'
            correction_shift = 0
        else:
            reads_cnt['signal'] += 1
            ss_id, strand, us_shift, ds_shift = ss_site
            circ_start += us_shift
            circ_end += ds_shift
            correction_shift = min(max(us_shift, us_free), ds_free)

        circ_id = '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end)

        # Get Cirexons
        cir_exons = get_blocks(circ_hit)
        cir_exons = merge_clip_exon(cir_exons, clip_info)

        cir_exons[0][0] = circ_start
        cir_exons[-1][1] = circ_end

        cir_exon_tag = []
        for cir_exon_start, cir_exon_end, cir_exon_length in cir_exons:
            cir_exon_tag.append('{}-{}|{}'.format(cir_exon_start + 1, cir_exon_end, cir_exon_length))

        # BSJ correction for 5' prime region
        circ_seq = clipped_circ if circ_hit.strand > 0 else revcomp(clipped_circ)
        circ_seq = circ_seq[correction_shift:] + circ_seq[:correction_shift]

        ret.append((
            read_id, circ_id, strand, ','.join(cir_exon_tag), ss_id,
            '{}|{}-{}'.format(junc, clip_base, len(circ)), segments, circ_seq
        ))

    return reads_cnt, ret


def recover_ccs_reads(short_reads, ref_fasta, ss_index, gtf_index, intron_index, is_canonical, out_dir, prefix, threads):
    from bwapy import BwaAligner

    # Second scanning of short reads
    genome = Fasta(ref_fasta)

    options = '-x ont2d -T 19'
    bwa_aligner = Aligner(BwaAligner(ref_fasta, options=options))

    chunk_size = 250
    jobs = []
    pool = Pool(threads, env.initializer, (bwa_aligner, genome.contig_len, genome, gtf_index, intron_index, ss_index))
    for reads in grouper(short_reads, chunk_size):
        chunk = [i for i in reads if i is not None]
        jobs.append(pool.apply_async(recover_ccs_chunk, (chunk, is_canonical)))
    pool.close()

    prog = ProgressBar()
    prog.update(0)
    finished_cnt = 0

    reads_count = defaultdict(int)
    with open('{}/{}.cand_circ.fa'.format(out_dir, prefix), 'a') as out:
        for job in jobs:
            finished_cnt += 1

            tmp_cnt, ret = job.get()
            for key, value in tmp_cnt.items():
                reads_count[key] += value

            for read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, segments, circ_seq in ret:
                out.write('>{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\n'.format(
                    read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, segments, circ_seq
                ))
            prog.update(100 * finished_cnt / len(jobs))

    pool.join()
    prog.update(100)

    return reads_count


def check_read(segments, seq):
    from spoa import poa
    fasta = [seq[int(i.split('-')[0]):int(i.split('-')[1])] for i in segments.split(';')]
    poa(fasta, 1, True, -1, -1, -1, -1, -1)


def scan_raw_chunk(chunk, is_canonical, circ_reads):
    reads_cnt = defaultdict(int)

    ret = []
    short_reads = []

    for read_id, seq in chunk:
        if read_id in circ_reads:
            continue

        # API for short reads
        if len(seq) < 300:
            short_reads.append((read_id, seq))
            continue

        # Remove reads that have ambiguous mapping
        raw_hits = sorted([i for i in env.ALIGNER.map(seq) if i.is_primary], key=lambda x: [x.q_st, x.q_en])
        if len(raw_hits) == 0:
            continue
        elif len(raw_hits) == 1:
            raw_hit = remove_long_insert(raw_hits[0])
            if raw_hit.mlen < len(seq) * .45 or raw_hit.mlen > len(seq) - 50:
                continue
            if raw_hit.q_st < 50 and raw_hit.q_en > len(seq) - 50:
                continue
            circ, junc = find_bsj(seq)
            if junc is None:
                continue

        elif len(raw_hits) == 2:
            head, tail = remove_long_insert(raw_hits[0]), remove_long_insert(raw_hits[1])
            if head.ctg != tail.ctg:
                continue
            if not head.q_st + head.mlen * 0.45 < tail.q_st:
                continue
            if head.r_en - 20 < tail.r_st:
                continue
            if head.q_en < tail.q_st - 50:
                continue
            circ, junc = find_bsj(seq)
            if junc is None or junc < head.q_en - 10 or junc > tail.q_st + 10:
                continue
        else:
            continue

        circ_hits = sorted([remove_long_insert(i) for i in env.ALIGNER.map(circ) if i.is_primary], key=lambda x: [x.q_st, x.q_en])
        if len(circ_hits) == 0:
            continue
        elif len(circ_hits) == 1:
            circ_hit = circ_hits[0]
            if circ_hit.mlen <= max([i.mlen for i in raw_hits]):
                continue
            if min(junc, len(seq) - junc) < 30:
                continue
            if not junc + circ_hit.q_st < len(seq) < junc + circ_hit.q_en:
                continue
            circ_ctg, circ_start, circ_end, circ_strand = circ_hit.ctg, circ_hit.r_st, circ_hit.r_en, circ_hit.strand
            clip_base = circ_hit.q_st + len(seq) - circ_hit.q_en
            cir_exons = get_parital_blocks(circ_hit, len(seq) - junc)
        elif len(circ_hits) == 2:
            head, tail = circ_hits[0], circ_hits[1]
            if head.ctg != tail.ctg or head.strand != tail.strand:
                continue
            if not head.q_st + (head.q_en - head.q_st) * 0.5 < tail.q_st:
                continue
            if head.r_en - 20 < tail.r_st:
                continue
            if head.q_en < tail.q_st - 20:
                continue
            circ_ctg, circ_start, circ_end, circ_strand = head.ctg, tail.r_st, head.r_en, head.strand
            clip_base = abs(tail.q_st - head.q_en)

            head_exons = get_blocks(head)
            tail_exons = get_blocks(tail)

            cir_exons = merge_exons(tail_exons, head_exons)

            circ = circ[tail.q_st:] + circ[:tail.q_st]
        else:
            continue

        if clip_base > 20:
            continue

        # Retrive circRNA positions, convert minimap2 position to real position
        host_strand = find_host_gene(circ_ctg, circ_start, circ_end)
        try:
            ss_site, us_free, ds_free, tmp_signal = find_annotated_signal(circ_ctg, circ_start, circ_end, clip_base, clip_base + 10)
        except Exception as e:
            print(e)
        if ss_site is None:
            ss_site = find_denovo_signal(circ_ctg, circ_start, circ_end, host_strand, tmp_signal,
                                         us_free, ds_free, clip_base, clip_base + 10, 3, True)

        if ss_site is None:
            strand = 'NA'
            ss_id = 'NA'
            correction_shift = 0
        else:
            ss_id, strand, us_shift, ds_shift = ss_site
            circ_start += us_shift
            circ_end += ds_shift
            correction_shift = min(max(us_shift, -us_free), ds_free)

        circ_id = '{}:{}-{}'.format(circ_ctg, circ_start + 1, circ_end)
        cir_exons[0][0] = circ_start
        cir_exons[-1][1] = circ_end

        cir_exon_tag = []
        for cir_exon_start, cir_exon_end, cir_exon_len in cir_exons:
            cir_exon_tag.append('{}-{}|{}'.format(cir_exon_start, cir_exon_end, cir_exon_len))

        circ_seq = circ if circ_strand > 0 else revcomp(circ)
        circ_seq = circ_seq[correction_shift:] + circ_seq[:correction_shift]

        ret.append((
            read_id, circ_id, strand, ','.join(cir_exon_tag), ss_id, '{}|{}-NA'.format(junc, clip_base), 'partial', circ_seq
        ))

        reads_cnt['partial'] += 1

    return reads_cnt, ret, short_reads


def scan_raw_reads(in_file, ref_fasta, gtf_index, intron_index, ss_index, is_canonical, out_dir, prefix, threads):
    import gzip
    import mappy as mp
    from CIRI_long.utils import to_str

    circ_reads = {}
    with open('{}/{}.cand_circ.fa'.format(out_dir, prefix), 'r') as f:
        for line in f:
            read_id = line.rstrip().split()[0].lstrip('>')
            circ_reads[read_id] = 1
            f.readline()

    # Auto detect input format
    is_fastq = 1
    is_gz = 1
    if in_file.endswith('.fa') or in_file.endswith('.fasta'):
        is_fastq = 0
        is_gz = 0
        fq = open(in_file, 'r')
    elif in_file.endswith('.fa.gz') or in_file.endswith('.fasta.gz'):
        is_fastq = 0
        is_gz = 1
        fq = gzip.open(in_file, 'rb')
    elif in_file.endswith('.fq') or in_file.endswith('.fastq'):
        is_gz = 0
        fq = open(in_file, 'r')
    elif in_file.endswith('.fq.gz') or in_file.endswith('.fastq.gz'):
        fq = gzip.open(in_file, 'rb')
    else:
        sys.exit('Wrong format of input')

    # Prepare aligners
    faidx = Faidx(ref_fasta)
    contig_len = faidx.contig_len
    faidx.close()

    aligner = mp.Aligner(ref_fasta, n_threads=threads, preset='splice')

    jobs = []
    pool = Pool(threads, env.initializer, (aligner, contig_len, aligner, gtf_index, intron_index, ss_index))

    # Init jobs
    chunk = []
    chunk_size = 1000
    for line in fq:
        if is_gz:
            header = to_str(line).rstrip().split(' ')[0]
            seq = to_str(fq.readline()).rstrip()
        else:
            header = line.rstrip().split(' ')[0]
            seq = fq.readline().rstrip()

        if is_fastq:
            header = header.lstrip('@')
            fq.readline()
            fq.readline()
        else:
            header = header.lstrip('>')

        # Split into chunks
        chunk.append((header, seq))
        if len(chunk) == chunk_size:
            jobs.append(pool.apply_async(scan_raw_chunk, (chunk, is_canonical, circ_reads)))
            chunk = []

    if len(chunk) > 0:
        jobs.append(pool.apply_async(scan_raw_chunk, (chunk, is_canonical, circ_reads)))
    pool.close()

    # Receive results
    reads_cnt = defaultdict(int)

    prog = ProgressBar()
    prog.update(0)
    finished_job = 0

    short_reads = []
    with open('{}/{}.low_confidence.fa'.format(out_dir, prefix), 'w') as out:
        for job in jobs:
            tmp_cnt, tmp_ret, tmp_short = job.get()
            for key, value in tmp_cnt.items():
                reads_cnt[key] += value
            short_reads += tmp_short

            for read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, segments, circ_seq in tmp_ret:
                out.write('>{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\n'.format(
                    read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, segments, circ_seq
                ))
            finished_job += 1

            prog.update(100 * finished_job // len(jobs))

    prog.update(100)
    pool.join()

    return reads_cnt, short_reads


# def recover_raw_chunk(chunk, is_canonical):
#     reads_cnt = defaultdict(int)
#
#     ret = []
#     for read_id, seq in chunk:
#         raw_hits = ALIGNER.map(seq)
#         if raw_hits is None:
#             continue
#         raw_hits = sorted([i for i in raw_hits if i.is_primary], key=lambda x: [x.q_st, x.q_en])
#         if len(raw_hits) == 0:
#             continue
#
#         if len(raw_hits) == 1:
#             raw_hit = raw_hits[0]
#             if raw_hit.q_en - raw_hit.q_st < len(seq) * .45:
#                 continue
#             if raw_hit.q_en - raw_hit.q_st > len(seq) - 50:
#                 continue
#             if raw_hit.q_st < 50 and len(seq) - raw_hit.q_en < 50:
#                 continue
#
#             circ, junc = find_bsj(seq)
#             if junc is None:
#                 continue
#             circ_hits = ALIGNER.map(circ)
#             if circ_hits is None:
#                 continue
#             circ_hits = sorted([i for i in circ_hits if i.is_primary], key=lambda x: [20, x.q_en])
#             if len(circ_hits) == 0:
#                 continue
#
#             if len(circ_hits) == 1:
#                 circ_hit = circ_hits[0]
#                 if circ_hit.mlen <= raw_hit.mlen:
#                     continue
#                 if min(junc, len(seq) - junc) < 30:
#                     continue
#                 if not junc + circ_hit.q_st < len(seq) < junc + circ_hit.q_en:
#                     continue
#                 circ_ctg, circ_start, circ_end = circ_hit.ctg, circ_hit.r_st, circ_hit.r_en
#                 clip_base = circ_hit.q_st + len(seq) - circ_hit.q_en
#
#             elif len(circ_hits) == 2:
#                 head, tail = circ_hits[0], circ_hits[1]
#                 if head.ctg != tail.ctg:
#                     continue
#                 if not head.q_st + (head.q_en - head.q_st) * 0.5 < tail.q_st:
#                     continue
#                 if head.r_en - 20 < tail.r_st:
#                     continue
#                 if head.q_en < tail.q_st - 50:
#                     continue
#                 circ_ctg, circ_start, circ_end = head.ctg, tail.r_st, head.r_en
#                 clip_base = abs(tail.q_st - head.q_en)
#             else:
#                 continue
#
#         elif len(raw_hits) == 2:
#             head, tail = raw_hits[0], raw_hits[1]
#             if head.ctg != tail.ctg:
#                 continue
#             if not head.q_st + (head.q_en - head.q_st) * 0.5 < tail.q_st:
#                 continue
#             if head.r_en - 20 < tail.r_st:
#                 continue
#             if head.q_en < tail.q_st - 50:
#                 continue
#             circ_ctg, circ_start, circ_end = head.ctg, tail.r_st, head.r_en
#             clip_base = abs(tail.q_st - head.q_en)
#
#         else:
#             continue
#
#         # Retrive circRNA positions, convert minimap2 position to real position
#         ss_site, us_free, ds_free = search_splice_signal(circ_ctg, circ_start, circ_end, clip_base)
#         if ss_site is None:
#             continue
#
#         ss_id, strand, us_shift, ds_shift = ss_site
#         circ_start += us_shift
#         circ_end += ds_shift
#
#         # if is_canonical: keep canonical splice site only
#         ss = ss_id.split('|')[0]
#         if ss_site is None or (is_canonical and ss[-1] == '*'):
#             continue
#
#         circ_id = '{}:{}-{}'.format(circ_ctg, circ_start + 1, circ_end)
#
#     return reads_cnt, ret
#
#
# def recover_raw_reads(short_reads, ref_fasta, ss_index, is_canonical, out_dir, prefix, threads):
#     from bwapy import BwaAligner
#
#     faidx = Faidx(ref_fasta)
#     contig_len = faidx.contig_len
#
#     options = '-x ont2d -T 19'
#     bwa_aligner = Aligner(BwaAligner(ref_fasta, options=options))
#
#     chunk_size = 250
#     jobs = []
#     pool = Pool(threads, initializer, (bwa_aligner, faidx, ss_index, contig_len))
#     for reads in grouper(short_reads, chunk_size):
#         chunk = [i for i in reads if i is not None]
#         jobs.append(pool.apply_async(recover_raw_chunk, (chunk, is_canonical)))
#     pool.close()
#
#     prog = ProgressBar()
#     prog.update(0)
#     finished_job = 0
#
#     reads_cnt = defaultdict(int)
#     with open('{}/{}.low_confidence.fa'.format(out_dir, prefix), 'a') as out:
#         for job in jobs:
#             finished_job += 1
#
#             tmp_cnt, tmp_ret = job.get()
#             for key, value in tmp_cnt.items():
#                 reads_cnt[key] += value
#
#             prog.update(100 * finished_job // len(jobs))
#
#     pool.join()
#     prog.update(100)
#
#     faidx.close()
#
#     return reads_cnt
