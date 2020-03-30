import logging
import pandas as pd
import numpy as np
from multiprocessing import Pool
from collections import defaultdict, Counter, namedtuple

from CIRI.utils import tree, grouper
from CIRI.logger import ProgressBar

LOGGER = logging.getLogger('CIRI-long')
READ = namedtuple('Read', 'read_id circ_id sample seq type')
globals()['Read'] = READ


def load_cand_circ(in_file):
    sample_attr = {}
    with open(in_file, 'r') as f:
        for line in f:
            sample, fname = line.rstrip().split('\t')
            sample_attr[sample] = fname

    cand_reads = {}
    for sample, fname in sample_attr.items():
        with open(fname, 'r') as f:
            for line in f:
                content = line.rstrip().lstrip('>').split('\t')
                clip_base = int(content[5].split('|')[1].split('-')[0])
                seq = f.readline().rstrip()
                if clip_base > 20:
                    continue

                read_id, circ_id = content[0], content[1]
                read_type = 'partial' if content[6] == 'partial' else 'full'
                cand_reads[read_id] = READ(read_id, circ_id, sample, seq, read_type)

    return cand_reads


CIRC_READS = defaultdict(list)
CIRC_START = defaultdict(dict)
CIRC_END = defaultdict(dict)


def load_index(read):
    import re
    global CIRC_READS, CIRC_START, CIRC_END
    contig, start, end = re.split('[:-]', read.circ_id)
    start, end = int(start), int(end)

    # Store reads & BSJ sites
    CIRC_READS[contig].append((start, end, read.read_id))
    CIRC_START[contig].setdefault(start, []).append(read.read_id)
    CIRC_END[contig].setdefault(end, []).append(read.read_id)


def cluster_reads(cand_reads):
    from operator import itemgetter

    # Load circRNA junction sites
    for read_id, read in cand_reads.items():
        load_index(read)

    # Build index of circRNA junction sites
    reads_cluster = []
    for contig in CIRC_READS:
        circ_start_index = {}
        circ_end_index = {}

        # Cluster of start junction site
        tmp = [[], ]
        for x in sorted(CIRC_START[contig]):
            if not tmp[-1]:
                tmp[-1].append(x)
            elif x > tmp[-1][-1] + 20:
                tmp.append([x, ])
            else:
                tmp[-1].append(x)
        for x in tmp:
            for i in range(min(x) // 500, max(x) // 500 + 1):
                circ_start_index.setdefault(i, []).append(x)

        # Cluster of end junction site
        tmp = [[], ]
        for x in sorted(CIRC_END[contig]):
            if not tmp[-1]:
                tmp[-1].append(x)
            elif x > tmp[-1][-1] + 20:
                tmp.append([x, ])
            else:
                tmp[-1].append(x)
        for x in tmp:
            for i in range(min(x) // 500, max(x) // 500 + 1):
                circ_end_index.setdefault(i, []).append(x)

        # Cluster reads
        reads_itered = {}

        for (start, end, read_id) in sorted(CIRC_READS[contig], key=itemgetter(0, 1)):
            if read_id in reads_itered:
                continue

            tmp_reads = []
            p = [i for i in circ_start_index[start // 500] if start in i][0]
            q = [i for i in circ_end_index[end // 500] if end in i][0]

            for i in p:
                tmp_start = CIRC_START[contig][i]
                for j in q:
                    tmp_end = CIRC_END[contig][j]
                    tmp = set(tmp_start) & set(tmp_end)
                    if tmp:
                        tmp_reads += tmp

            for i in tmp_reads:
                reads_itered[i] = 1

            reads_cluster.append(sorted([cand_reads[i] for i in tmp_reads], key=lambda x: len(x.seq), reverse=True))
    return reads_cluster


def shift_base(seq):
    base = 0
    for i in seq:
        if i in 'atcgATGC':
            break
        else:
            base += 1
    return base


def transform_seq(seq, bsj):
    return seq[bsj:] + seq[:bsj]


def correct_chunk(chunk):
    from CIRI.poa import consensus
    from libs.striped_smith_waterman.ssw_wrap import Aligner

    cs_cluster = []
    for cluster in chunk:
        if cluster is None:
            continue
        if len(cluster) <= 1:
            continue

        ref = cluster[0]
        ssw = Aligner(ref.seq, match=10, mismatch=4, gap_open=8, gap_extend=2)

        head_pos = []
        for query in cluster[1:]:
            alignment = ssw.align(query.seq)

            head_pos.append(alignment.ref_begin)
            if alignment.ref_begin > 20:
                from spoa import poa
                _, msa = poa([ref.seq, query.seq], 0, True, 10, -4, -8, -2, -24, -4)
                print(msa)

        trans_reads = []
        template = transform_seq(ref.seq, max(head_pos))

        ssw = Aligner(template, match=10, mismatch=4, gap_open=8, gap_extend=2)
        trans_reads.append(READ(ref.read_id, ref.circ_id, ref.sample, template, ref.type))

        for query in cluster[1:]:
            alignment = ssw.align(query.seq)
            trans_reads.append(READ(query.read_id, query.circ_id, query.sample,
                                    transform_seq(query.seq, alignment.query_begin), query.type))

        cs_reads = consensus([(i.read_id, i.seq) for i in trans_reads], alignment_type=0,
                             match=10, mismatch=-4, gap=-8, extension=-2, gap_affine=-24, extension_affine=-4,
                             debug=0)

        # from spoa import poa
        # poa([i.seq for i in trans_reads], 0, True, 10, -4, -8, -2, -24, -4)

        cs_cluster.append(([i.read_id for i in trans_reads], cs_reads))

    return cs_cluster


def correct_reads(reads_cluster, threads):
    corrected_reads = []

    jobs = []
    pool = Pool(threads)
    for cluster in grouper(reads_cluster, 1000):
        jobs.append(pool.apply_async(correct_chunk, (cluster,)))
    pool.close()

    prog = ProgressBar()
    prog.update(0)
    cnt = 0
    for job in jobs:
        tmp_cluster = job.get()
        corrected_reads += tmp_cluster
        cnt += 1
        prog.update(100 * cnt // len(jobs))
    pool.join()
    prog.update(100)

    return corrected_reads


def scan_corrected_chunk(chunk):
    from CIRI import alignment
    from CIRI.preprocess import revcomp

    ret = []
    short_reads = []
    for reads, ccs in chunk:
        ccs_hit = alignment.get_primary_alignment(alignment.ALIGNER.map(ccs * 2))
        if ccs_hit is None and len(ccs) < 150:
            short_reads.append((reads, ccs))
            continue

        # Find back-spliced junction site
        circ, junc = alignment.find_bsj(ccs)

        # Candidate alignment situation, more than 85%
        circ_hit = alignment.get_primary_alignment(alignment.ALIGNER.map(circ))
        if circ_hit is None or circ_hit.mlen < 0.75 * len(circ):
            continue

        circ_start, circ_end, clip_info = alignment.align_clip_segments(circ, circ_hit)
        if circ_start is None or circ_end is None:
            continue

        clip_base = clip_info[2]
        if clip_base > 0.15 * len(ccs):
            continue

        # Retrive circRNA positions, convert minimap2 position to real position
        ss_site, us_free, ds_free = alignment.search_splice_signal(circ_hit.ctg, circ_start, circ_end, clip_base, clip_base + 10)
        if ss_site is None:
            circ_id = '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end)
            circ_seq = circ if circ_hit.strand > 0 else revcomp(circ)
            ret.append([reads, circ_id, circ_seq])
            continue

        ss_id, strand, us_shift, ds_shift = ss_site
        circ_start += us_shift
        circ_end += ds_shift

        # if is_canonical: keep canonical splice site only
        ss = ss_id.split('|')[0]
        # if is_canonical and ss[-1] == '*' and ss != 'AG-GT*':
        #     continue

        circ_id = '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end)

        # BSJ correction for 5' prime region
        correction_shift = min(max(us_shift, us_free), ds_free)
        circ_seq = circ if circ_hit.strand > 0 else revcomp(circ)
        circ_seq = circ_seq[correction_shift:] + circ_seq[:correction_shift]

        ret.append([reads, circ_id, circ_seq])

    return ret, short_reads


def scan_corrected_reads(corrected_reads, ref_fasta, ss_index, threads):
    import mappy as mp
    from CIRI import alignment

    faidx = alignment.Faidx(ref_fasta)
    contig_len = faidx.contig_len
    faidx.close()

    # First scanning using minimap2
    minimap_aligner = mp.Aligner(ref_fasta, n_threads=threads, preset='splice')

    chunk_size = 250
    jobs = []
    pool = Pool(threads, alignment.initializer, (minimap_aligner, minimap_aligner, ss_index, contig_len))
    for chunk in grouper(corrected_reads, chunk_size):
        jobs.append(pool.apply_async(scan_corrected_chunk, (chunk, )))
    pool.close()

    prog = ProgressBar()
    cnt = 0

    corrected_circ = []
    short_reads = []
    for job in jobs:
        tmp_res, tmp_short = job.get()
        corrected_circ += tmp_res
        short_reads += tmp_short
        cnt += 1
        prog.update(100 * cnt // len(jobs))
    pool.join()
    prog.update(100)

    return corrected_circ, short_reads


def recover_corrected_chunk(chunk, is_canonical):
    from CIRI import alignment
    from CIRI.preprocess import revcomp

    ret = []

    for reads, ccs in chunk:
        ccs_hit = alignment.get_primary_alignment(alignment.ALIGNER.map(ccs * 2))
        if ccs_hit is None:
            continue

        # Find back-spliced junction site
        circ, junc = alignment.find_bsj(ccs)

        # Candidate alignment situation, more than 85%
        circ_hit = alignment.get_primary_alignment(alignment.ALIGNER.map(circ))
        if circ_hit is None:
            continue

        circ_start, circ_end, clip_info = alignment.align_clip_segments(circ, circ_hit)
        if circ_start is None or circ_end is None:
            continue

        clip_base = clip_info[2]
        if clip_base > 0.15 * len(ccs):
            continue

        # Retrive circRNA positions, convert minimap2 position to real position
        ss_site, us_free, ds_free = alignment.search_splice_signal(circ_hit.ctg, circ_start, circ_end, clip_base, clip_base + 10)
        if ss_site is None:
            circ_id = '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end)
            circ_seq = circ if circ_hit.strand > 0 else revcomp(circ)
            ret.append([reads, circ_id, circ_seq])
            continue

        ss_id, strand, us_shift, ds_shift = ss_site
        circ_start += us_shift
        circ_end += ds_shift

        # if is_canonical: keep canonical splice site only
        ss = ss_id.split('|')[0]
        # if is_canonical and ss[-1] == '*' and ss != 'AG-GT*':
        #     continue

        circ_id = '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end)

        # BSJ correction for 5' prime region
        correction_shift = min(max(us_shift, us_free), ds_free)
        circ_seq = circ if circ_hit.strand > 0 else revcomp(circ)
        circ_seq = circ_seq[correction_shift:] + circ_seq[:correction_shift]

        ret.append([reads, circ_id, circ_seq])

    return ret


def recover_corrected_reads(short_reads, ref_fasta, ss_index):
    from bwapy import BwaAligner
    from CIRI import alignment

    # Second scanning of short reads
    faidx = alignment.Faidx(ref_fasta)
    contig_len = faidx.contig_len

    options = '-x ont2d -T 19'
    bwa_aligner = alignment.Aligner(BwaAligner(ref_fasta, options=options))

    chunk_size = 250
    jobs = []
    pool = Pool(1, alignment.initializer, (bwa_aligner, faidx, ss_index, contig_len))
    for reads in grouper(short_reads, chunk_size):
        chunk = [i for i in reads if i is not None]
        jobs.append(pool.apply_async(recover_corrected_chunk, (chunk, )))
    pool.close()

    prog = ProgressBar()
    prog.update(0)
    cnt = 0

    corrected_circ = []
    for job in jobs:
        tmp_ret = job.get()
        corrected_circ += tmp_ret
        cnt += 1
        prog.update(100 * cnt / len(jobs))

    pool.join()
    prog.update(100)

    faidx.close()

    return corrected_circ
