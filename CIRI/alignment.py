import os
import sys
import re
import logging
import subprocess
from collections import defaultdict

import pysam
import numpy as np

from CIRI.logger import ProgressBar
from CIRI.preprocess import revcomp

LOGGER = logging.getLogger('CIRI-long')


# Parsing GTF format
GENE_ID = re.compile(r'gene_id "(\S+)";')
GENE_BIOTYPE = re.compile(r'gene_biotype "(\S+)";')
GENE_TYPE = re.compile(r'gene_type "(\S+)";')
GENE_NAME = re.compile(r'gene_name "(\S+)";')
TSCP_ID = re.compile(r'transcript_id "(\S+)";')
TSCP_NAME = re.compile(r'transcript_name "(\S+)";')
TSCP_BIOTYPE = re.compile(r'transcript_biotype "(\S+)";')

ATTR = {
    # 'gene_id': GENE_ID,
    # 'gene_name': GENE_NAME,
    'gene_biotype': GENE_BIOTYPE,
    'gene_type': GENE_TYPE
    # 'transcript_id': TSCP_ID,
    # 'transcript_name': TSCP_NAME,
    # 'transcript_biotype': TSCP_BIOTYPE,
}


class GTFParser(object):
    """
    Class for parsing annotation gtf
    """

    def __init__(self, content):
        self.contig = content[0]
        self.source = content[1]
        self.type = content[2]
        self.start, self.end = int(content[3]), int(content[4])
        self.strand = content[6]
        self.attr_string = content[-1]

    def attr(self, key):
        """match attribute string for row"""
        try:
            return ATTR[key].search(self.attr_string).group(1)
        except AttributeError:
            return None


class Aligner(object):
    """
    API for unify bwapy to mappy::Aligner
    """
    def __init__(self, aligner):
        self.aligner = aligner

    def map(self, seq):
        alignments = self.aligner.align_seq(seq)
        if alignments:
            hits = [Hit(aln) for aln in alignments]
            hits[0].is_primary = 1
            return hits
        else:
            return None


OPERATION = {
    'M': 0,
    'I': 1,
    'D': 2,
    'N': 3,
    'S': 4,
    'H': 5,
    'P': 6,
    '=': 7,
    'X': 8,
}


class Hit(object):
    """
    API for bwapy alignment
    Alignment(rname='chr7', orient='-', pos=131329219, mapq=60, cigar='51S37M', NM=0)
    """

    def __init__(self, aln):
        self.ctg = aln.rname
        self.strand = 1 if aln.orient == '+' else -1
        self.cigar_string = aln.cigar
        self.r_st = aln.pos
        self.r_en, self.q_st, self.q_en = self.__parse_cigar()
        self.mlen = self.q_en - self.q_st
        self.is_primary = 0

    @property
    def cigar(self):
        return [(int(length), OPERATION[operation]) for length, operation in
                re.findall(r'(\d+)([MIDNSHP=X])', self.cigar_string)]

    def __parse_cigar(self):
        r_en = self.r_st
        q_st, q_en = 0, 0
        for length, operation in self.cigar:
            if operation == 0:
                q_en += length
                r_en += length
            elif operation == 1:
                q_en += length
            elif operation in [2, 3]:
                r_en += length
            elif operation in [4, 5]:
                if q_st == 0:
                    q_st = length
                    q_en = length
            else:
                pass
        return r_en, q_st, q_en

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}'.format(self.ctg, self.r_st, self.r_en, self.q_st, self.q_en, self.cigar_string)


class Faidx(object):
    """
    API for unify pysam::Fastafile and mappy2::Aligner
    """
    def __init__(self, infile):
        if not os.path.exists(infile + '.fai'):
            LOGGER.warn('Index of reference genome not found, generating ...')
        try:
            faidx = pysam.FastaFile(infile)
        except ValueError:
            LOGGER.error('Cannot generate index of reference genome, index file is missing')
            sys.exit()
        except IOError:
            LOGGER.error('Cannot generate index of reference genome, file could not be opened')
            sys.exit()

        self.faidx = faidx
        self.contig_len = {contig: faidx.get_reference_length(contig) for contig in faidx.references}

    def seq(self, contig, start, end):
        return self.faidx.fetch(contig, start, end)

    def close(self):
        self.faidx.close()


def index_annotation(gtf):
    """
    Generate binned index for element in gtf
    """
    from CIRI.utils import tree

    LOGGER.info('Loading annotation gtf ..')
    gtf_index = defaultdict(dict)
    splice_site_index = tree()
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            content = line.rstrip().split('\t')
            # only include gene and exon feature for now
            if content[2] not in ['gene', 'exon']:
                continue

            parser = GTFParser(content)
            # Extract splice site
            if content[2] == 'exon':
                splice_site_index[parser.contig][parser.start][parser.strand] = 1
                splice_site_index[parser.contig][parser.end][parser.strand] = 1

            # if parser.attr('gene_biotype') is None or parser.attr('gene_biotype') in ['lincRNA', 'pseudogene']:
            #     continue
            # if parser.attr('gene_type') is None or parser.attr('gene_type') in ['lincRNA', 'pseudogene']:
            #     continue

            # Binned index
            start_div, end_div = parser.start // 500, parser.end // 500
            for i in range(start_div, end_div + 1):
                gtf_index[parser.contig].setdefault(i, []).append(parser)

    return gtf_index, splice_site_index


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
            r_block.append([r_start, r_end, r_end - r_start + 1])
            r_start = r_end + length
            r_end = r_start
        elif operation == 4:
            pass
    if r_end > r_start:
        r_block.append([r_start, r_end, r_end - r_start + 1])
    return r_block


def merge_exons(exons, clip_info):
    clip_st, clip_en, clip_base = clip_info
    exon_st, exon_en = exons[0][0], exons[-1][1]

    if clip_st and clip_en:
        if clip_en < exon_st:
            exons = [[clip_st, clip_en, clip_base], ] + exons
        elif exon_en < clip_st:
            exons = exons + [[clip_st, clip_en, clip_base], ]
        elif clip_st < exon_st < clip_en:
            exons[0] = [clip_st, exons[0][1], exons[0][1] - clip_st + 1]
        elif clip_st < exon_en < clip_en:
            exons[-1] = [exons[-1][0], clip_en, clip_en - exons[-1][0] + 1]
        else:
            pass

    return exons

SPLICE_SIGNAL = {
    ('GT', 'AG'): 0, # U2-type
    ('GC', 'AG'): 1, # U2-type
    ('AT', 'AC'): 2, # U12-type
    # ('GT', 'TG'): 3, # non-canonical
    # ('AT', 'AG'): 3, # non-canonical
    # ('GA', 'AG'): 3, # non-canonical
    # ('GG', 'AG'): 3, # non-canonical
    # ('GT', 'GG'): 3, # non-canonical
    # ('GT', 'AT'): 4, # non-canonical
    # ('GT', 'AA'): 4, # non-canonical
}


def search_splice_signal(contig, start, end, clip_base, search_length=10, shift_threshold=3):
    # Find free sliding region
    # start | real_start <-> end | real_end
    ds_free = 0
    for i in range(100):
        if FAIDX.seq(contig, start, start + i) == FAIDX.seq(contig, end, end + i):
            ds_free = i
        else:
            break

    us_free = 0
    for j in range(100):
        if FAIDX.seq(contig, start - j, start) == FAIDX.seq(contig, end - j, end):
            us_free = j
        else:
            break

    # Splice site: site_id, strand, us_shift, ds_shift, site_weight, altered_len, altered_total
    # First: Find flanking junction from annotation gtf
    if SS_INDEX is not None:
        anno_ss = []
        for strand in ['+', '-']:
            tmp_us_sites = []
            for us_shift in range(-search_length, search_length):
                us_pos = start + us_shift
                if contig in SS_INDEX and us_pos in SS_INDEX[contig] and strand in SS_INDEX[contig][us_pos]:
                    tmp_us_sites.append(us_shift - 1)

            tmp_ds_sites = []
            for ds_shift in range(-search_length, search_length):
                ds_pos = end + ds_shift
                if contig in SS_INDEX and ds_pos in SS_INDEX[contig] and strand in SS_INDEX[contig][ds_pos]:
                    tmp_ds_sites.append(ds_shift)

            if len(tmp_us_sites) == 0 or len(tmp_ds_sites) == 0:
                continue

            for i in tmp_us_sites:
                for j in tmp_ds_sites:
                    if abs(i - j) > shift_threshold:
                        continue
                    us_ss = FAIDX.seq(contig, start + i - 2, start + i)
                    ds_ss = FAIDX.seq(contig, end + j, end + j + 2)
                    if strand == '-':
                        us_ss, ds_ss = revcomp(ds_ss), revcomp(us_ss)
                    ss_id = '{}-{}|{}-{}'.format(us_ss, ds_ss, i, j)
                    ss_weight = SPLICE_SIGNAL[(ds_ss, us_ss)] if (ds_ss, us_ss) in SPLICE_SIGNAL else 3
                    anno_ss.append((
                        ss_id, strand, i, j, ss_weight,
                        abs(i - j),
                        min(abs(i - ds_free), abs(i - us_free)) + min(abs(j - ds_free), abs(j - us_free)),
                    ))

        if len(anno_ss) > 0:
            return sort_splice_sites(anno_ss, us_free, ds_free, clip_base), us_free, ds_free

    # Second: Find Denovo BSJ using pre-defined splice signal
    us_search_length = search_length + us_free
    ds_search_length = search_length + ds_free
    us_seq = FAIDX.seq(contig, start - us_search_length - 2, start + ds_search_length)
    ds_seq = FAIDX.seq(contig, end - us_search_length, end + ds_search_length + 2)

    if us_seq is None or len(us_seq) < ds_search_length - us_search_length + 2:
        return None, us_free, ds_free
    if ds_seq is None or len(ds_seq) < ds_search_length - us_search_length + 2:
        return None, us_free, ds_free

    putative_ss = []
    for strand in ['+', '-']:
        for (ds_ss, us_ss), ss_weight in SPLICE_SIGNAL.items():
            ss_id = '{}-{}*'.format(us_ss, ds_ss)
            if strand == '-':
                ds_ss, us_ss = revcomp(us_ss), revcomp(ds_ss)

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
                    if abs(i - j) > shift_threshold:
                        continue
                    us_shift = i - us_search_length
                    ds_shift = j - us_search_length
                    putative_ss.append((
                        '{}|{}-{}'.format(ss_id, us_shift, ds_shift), strand, us_shift, ds_shift,
                        ss_weight,
                        abs(i - j),
                        min(abs(us_shift - ds_free), abs(us_shift - us_free)) + min(abs(ds_shift - ds_free), abs(ds_shift - us_free)),
                    ))

    if len(putative_ss) > 0:
        return sort_splice_sites(putative_ss, us_free, ds_free, clip_base), us_free, ds_free

    return None, us_free, ds_free


def sort_splice_sites(sites, us, ds, clip_base):
    # Splice site: site_id, strand, us_shift, ds_shift, site_weight, altered_len, altered_total
    from operator import itemgetter
    get_ss = itemgetter(0, 1, 2, 3)

    # Confidential splice sites
    confident_sites = [i for i in sites if -us <= i[2] <= ds and -us <= i[3] <= ds]
    if len(confident_sites) > 0:
        return get_ss(sorted(confident_sites, key=itemgetter(5, 4, 6))[0])

    # Ambiguous splice sites
    ambiguous_sites = [i for i in set(sites) - set(confident_sites) if -clip_base <= i[2] <= 0 <= i[3] <= clip_base]
    if len(ambiguous_sites) > 0:
        return get_ss(sorted(ambiguous_sites, key=itemgetter(4, 5, 6))[0])

    # Other sites
    normal_sites = set(sites) - set(confident_sites) - set(ambiguous_sites)
    if len(normal_sites) > 0:
        return get_ss(sorted(normal_sites, key=itemgetter(4, 5, 6))[0])
    return None


def get_primary_alignment(hits):
    if hits is None:
        return None

    for hit in hits:
        if hit.is_primary:
            return hit

    return None


def find_bsj(ccs):
    """
    Find junction using aligner
    """
    init_hit = get_primary_alignment(ALIGNER.map(ccs * 2))
    if init_hit is None:
        return None

    circ_junc = init_hit.q_st
    circ = ccs[circ_junc:] + ccs[:circ_junc]

    last_junc = 0
    last_m = 0

    last_junc = 0
    last_m = 0

    itered_junc = {}
    while True:
        circ_hit = get_primary_alignment(ALIGNER.map(circ))
        if circ_hit is None or circ_hit.mlen < last_m:
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
            circ_junc = (circ_junc - en_clip) % len(circ)

        if circ_junc in itered_junc:
            break

        circ = ccs[circ_junc:] + ccs[:circ_junc]
        itered_junc[circ_junc] = 1
    return circ


def align_clip_segments(circ, hit):
    """
    Align clip bases
    """
    from skbio import DNA, local_pairwise_align_ssw
    from collections import Counter
    st_clip, en_clip = hit.q_st, len(circ) - hit.q_en
    clip_r_st, clip_r_en = None, None

    if st_clip + en_clip >= 20:
        clip_seq = circ[hit.q_en:] + circ[:hit.q_st]
        if len(clip_seq) > 0.6 * len(circ):
            return None, None, None

        tmp_start = max(hit.r_st - 200000, 0)
        tmp_end = min(hit.r_en + 200000, CONTIG_LEN[hit.ctg])

        tmp_seq = FAIDX.seq(hit.ctg, tmp_start, tmp_end)
        if Counter(tmp_seq)['N'] >= 0.3 * (tmp_end - tmp_start):
            return None, None, None

        if hit.strand > 0:
            tmp_alignment, tmp_score, tmp_interval = local_pairwise_align_ssw(
                DNA(clip_seq), DNA(tmp_seq),
                gap_open_penalty=1, gap_extend_penalty=1, match_score=1, mismatch_score=-1
            )
            clip_r_st, clip_r_en = tmp_start + tmp_interval[1][0], tmp_start + tmp_interval[1][1]
        else:
            tmp_alignment, tmp_score, tmp_interval = local_pairwise_align_ssw(
                DNA(clip_seq), DNA(revcomp(tmp_seq)),
                gap_open_penalty=1, gap_extend_penalty=1, match_score=1, mismatch_score=-1
            )
            clip_r_st, clip_r_en = tmp_end - tmp_interval[1][1], tmp_end - tmp_interval[1][0]

        clip_base = hit.q_st + len(circ) - hit.q_en - (tmp_interval[1][1] - tmp_interval[1][0]) + 1
        circ_start = min(hit.r_st, clip_r_st) - 1
        circ_end = max(hit.r_en, clip_r_en)
    else:
        circ_start = hit.r_st - 1
        circ_end = hit.r_en
        clip_base = st_clip + en_clip

    return circ_start, circ_end, (clip_r_st, clip_r_en, clip_base)


FAIDX = None
ALIGNER = None
SS_INDEX = None
CONTIG_LEN = None
def initializer(aligner, faidx, splice_site_index, contig_len):
    global ALIGNER, FAIDX, SS_INDEX, CONTIG_LEN
    ALIGNER = aligner
    FAIDX = faidx
    SS_INDEX = splice_site_index
    CONTIG_LEN = contig_len


def scan_chunk(chunk, is_canonical):
    reads_cnt = defaultdict(int)
    ret = []

    short_reads = []
    for read_id, segments, ccs, raw in chunk:
        # Filter 1 - Remove linear mapped reads
        raw_hit = get_primary_alignment(ALIGNER.map(raw))
        if raw_hit and raw_hit.q_en - raw_hit.q_st > max(len(raw) * 0.8, len(raw) - 200):
            continue
        if raw_hit and raw_hit.q_en - raw_hit.q_st > 1.5 * len(ccs):
            continue

        raw_st = raw_hit.q_st if raw_hit else None
        raw_en = raw_hit.q_en if raw_hit else None
        reads_cnt['raw_unmapped'] += 1

        # Filter 2 - Remove other mapped region that intersect with ccs
        seg_st = int(segments.split(';')[0].split('-')[0])
        seg_en = int(segments.split(';')[-1].split('-')[1])
        if raw_hit and (raw_en < seg_st or raw_st > seg_en):
            continue

        ccs_hit = get_primary_alignment(ALIGNER.map(ccs * 2))
        if ccs_hit is None and len(ccs) < 150:
            short_reads.append((read_id, segments, ccs, raw))
        if ccs_hit is None or seg_en - seg_st < ccs_hit.q_en - ccs_hit.q_st:
            continue

        reads_cnt['ccs_mapped'] += 1

        # Find back-spliced junction site
        circ = find_bsj(ccs)

        # Candidate alignment situation, more than 85%
        circ_hit = get_primary_alignment(ALIGNER.map(circ))
        if circ_hit is None or circ_hit.mlen < 0.75 * len(circ):
            continue

        circ_start, circ_end, clip_info = align_clip_segments(circ, circ_hit)
        if circ_start is None or circ_end is None:
            continue

        clip_base = clip_info[2]
        if clip_base > 0.15 * len(ccs):
            continue

        reads_cnt['bsj'] += 1

        # Retrive circRNA positions, convert minimap2 position to real position
        ss_site, us_free, ds_free = search_splice_signal(circ_hit.ctg, circ_start, circ_end, clip_base)
        if ss_site is None:
            # Keep reads that seemed to generate from circRNAs, but cannot get splice signal
            if not is_canonical:
                ret.append((
                    read_id, '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end),
                    'NA', 'NA', 'NA', '{}-{}'.format(clip_base, len(ccs)), circ
                ))
            continue

        ss_id, strand, us_shift, ds_shift = ss_site
        circ_start += us_shift
        circ_end += ds_shift

        # if is_canonical: keep canonical splice site only
        ss = ss_id.split('|')[0]
        if is_canonical and ss[-1] == '*' and ss != 'AG-GT*':
            continue

        circ_id = '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end)

        # Get Cirexons
        cir_exons = get_blocks(circ_hit)
        cir_exons = merge_exons(cir_exons, clip_info)

        cir_exons[0][0] = circ_start
        cir_exons[-1][1] = circ_end

        cir_exon_tag = []
        for cir_exon_start, cir_exon_end, _ in cir_exons:
            cir_exon_tag.append('{}-{}'.format(cir_exon_start + 1, cir_exon_end))

        # BSJ correction for 5' prime region
        if strand == '+':
            correction_shift = min(max(us_shift, us_free), ds_free)
            circ_seq = circ[correction_shift:] + circ[:correction_shift]
        else:
            correction_shift = min(max(ds_shift, us_free), ds_free)
            circ_seq = revcomp(circ[correction_shift:] + circ[:correction_shift])

        ret.append((
            read_id, circ_id, strand, ','.join(cir_exon_tag), ss_id, '{}-{}'.format(clip_base, len(circ)), circ_seq
        ))
        reads_cnt['signal'] += 1

    return reads_cnt, short_reads, ret


def scan_ccs_reads(ccs_seq, ref_fasta, ss_index, is_canonical, out_dir, prefix, threads):
    from multiprocessing import Pool
    from CIRI.utils import grouper

    import mappy as mp

    faidx = Faidx(ref_fasta)
    contig_len = faidx.contig_len
    faidx.close()

    # First scanning using minimap2
    minimap_aligner = mp.Aligner(ref_fasta, n_threads=threads, preset='splice')

    chunk_size = 250
    chunk_cnt = 0
    jobs = []
    pool = Pool(threads, initializer, (minimap_aligner, minimap_aligner, ss_index, contig_len))
    for reads in grouper(list(ccs_seq), chunk_size):
        chunk = [[i, ] + ccs_seq[i] for i in reads if i is not None]
        chunk_cnt += 1
        jobs.append(pool.apply_async(scan_chunk, (chunk, is_canonical)))
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

            for read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, circ_seq in ret:
                out.write('>{}\t{}\t{}\t{}\t{}\t{}\n{}\n'.format(
                    read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, circ_seq
                ))
            prog.update(100 * finished_cnt / chunk_cnt)

    pool.join()
    prog.update(100)

    return reads_count, short_reads


def recover_chunk(chunk, is_canonical):
    reads_cnt = defaultdict(int)
    ret = []

    for read_id, segments, ccs, raw in chunk:
        # Remove other mapped region that intersect with ccs
        seg_st = int(segments.split(';')[0].split('-')[0])
        seg_en = int(segments.split(';')[-1].split('-')[1])

        ccs_hit = get_primary_alignment(ALIGNER.map(ccs * 2))
        if ccs_hit is None or seg_en - seg_st < ccs_hit.q_en - ccs_hit.q_st:
            continue

        reads_cnt['ccs_mapped'] += 1

        # Find back-spliced junction site
        circ = find_bsj(ccs)

        # Candidate alignment situation, more than 85%
        circ_hit = get_primary_alignment(ALIGNER.map(circ))
        if circ_hit is None:
            continue

        circ_start, circ_end, clip_info = align_clip_segments(circ, circ_hit)
        if circ_start is None or circ_end is None:
            continue

        clip_base = clip_info[2]
        if clip_base > 0.15 * len(ccs):
            continue

        reads_cnt['bsj'] += 1

        # Retrive circRNA positions, convert minimap2 position to real position
        ss_site, us_free, ds_free = search_splice_signal(circ_hit.ctg, circ_start, circ_end, clip_base)
        if ss_site is None:
            # Keep reads that seemed to generate from circRNAs, but cannot get splice signal
            if not is_canonical:
                ret.append((
                    read_id, '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end),
                    'NA', 'NA', 'NA', '{}-{}'.format(clip_base, len(ccs)), circ
                ))
            continue

        ss_id, strand, us_shift, ds_shift = ss_site
        circ_start += us_shift
        circ_end += ds_shift

        # if is_canonical: keep canonical splice site only
        ss = ss_id.split('|')[0]
        if is_canonical and ss[-1] == '*' and ss != 'AG-GT*':
            continue

        circ_id = '{}:{}-{}'.format(circ_hit.ctg, circ_start + 1, circ_end)

        # Get Cirexons
        cir_exons = get_blocks(circ_hit)
        cir_exons = merge_exons(cir_exons, clip_info)

        cir_exons[0][0] = circ_start
        cir_exons[-1][1] = circ_end

        cir_exon_tag = []
        for cir_exon_start, cir_exon_end, _ in cir_exons:
            cir_exon_tag.append('{}-{}'.format(cir_exon_start + 1, cir_exon_end))

        # BSJ correction for 5' prime region
        if strand == '+':
            correction_shift = min(max(us_shift, us_free), ds_free)
            circ_seq = circ[correction_shift:] + circ[:correction_shift]
        else:
            correction_shift = min(max(ds_shift, us_free), ds_free)
            circ_seq = revcomp(circ[correction_shift:] + circ[:correction_shift])

        ret.append((
            read_id, circ_id, strand, ','.join(cir_exon_tag), ss_id, '{}-{}'.format(clip_base, len(circ)), circ_seq
        ))
        reads_cnt['signal'] += 1

    return reads_cnt, ret


def recover_ccs_reads(short_reads, ref_fasta, ss_index, is_canonical, out_dir, prefix, threads):
    from multiprocessing import Pool
    from CIRI.utils import grouper

    from bwapy import BwaAligner

    # Second scanning of short reads
    faidx = Faidx(ref_fasta)
    contig_len = faidx.contig_len

    options = '-x ont2d -T 19'
    bwa_aligner = Aligner(BwaAligner(ref_fasta, options=options))

    chunk_size = 250
    chunk_cnt = 0
    jobs = []
    pool = Pool(1, initializer, (bwa_aligner, faidx, ss_index, contig_len))
    for reads in grouper(short_reads, chunk_size):
        chunk = [i for i in reads if i is not None]
        jobs.append(pool.apply_async(recover_chunk, (chunk, is_canonical)))
        chunk_cnt += 1
    pool.close()

    prog = ProgressBar()
    finished_cnt = 0

    reads_count = defaultdict(int)
    with open('{}/{}.cand_circ.fa'.format(out_dir, prefix), 'a') as out:
        for job in jobs:
            finished_cnt += 1

            tmp_cnt, ret = job.get()
            for key, value in tmp_cnt.items():
                reads_count[key] += value

            for read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, circ_seq in ret:
                out.write('>{}\t{}\t{}\t{}\t{}\t{}\n{}\n'.format(
                    read_id, circ_id, strand, cir_exon_tag, ss_id, clip_info, circ_seq
                ))
            prog.update(100 * finished_cnt / chunk_cnt)

    pool.join()
    prog.update(100)

    faidx.close()

    return reads_count
