import os
import sys
import re
import logging
from collections import defaultdict

import pysam
from CIRI_long import env
from CIRI_long.utils import *

LOGGER = logging.getLogger('CIRI-long')

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
    0: 'M',
    1: 'I',
    2: 'D',
    3: 'N',
    4: 'S',
    5: 'H',
    6: 'P',
    7: '=',
    9: 'X',
}

SPLICE_SIGNAL = {
    ('GT', 'AG'): 0,  # U2-type
    ('GC', 'AG'): 1,  # U2-type
    ('AT', 'AC'): 2,  # U12-type
    ('GT', 'AC'): 2,  # U12-type
    ('AT', 'AG'): 2,  # U12-type
    # ('GT', 'TG'): 3,  # non-canonical
    # ('AT', 'AG'): 3,  # non-canonical
    # ('GA', 'AG'): 3,  # non-canonical
    # ('GG', 'AG'): 3,  # non-canonical
    # ('GT', 'GG'): 3,  # non-canonical
    # ('GT', 'AT'): 4,  # non-canonical
    # ('GT', 'AA'): 4,  # non-canonical
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
        self.attr_string = content[8]

    @property
    def attr(self):
        """
        Parsing attribute column in gtf file
        """
        import re
        field = {}
        for attr_values in [re.split(r'\s+', i.strip()) for i in self.attr_string.split(';')[:-1]]:
            key, value = attr_values[0], attr_values[1:]
            field[key] = ' '.join(value).strip('"')
        return field


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
        return '\t'.join([str(x) for x in [self.q_st, self.q_en, self.ctg, self.r_st, self.r_en, self.mlen, self.blen,
                                           self.cigar_string]])


class SubHit(object):
    def __init__(self, hit, r_st, q_st, cigar):
        self.ctg = hit.ctg
        self.strand = hit.strand
        self.cigar = cigar
        self.r_st = r_st
        self.r_en, self.q_st, self.q_en = self.__parse_cigar(q_st)
        self.mlen, self.blen = self.__match()
        self.is_primary = 0

    def __parse_cigar(self, q_st):
        r_en = self.r_st
        q_en = q_st
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
                    q_st += length
                    q_en += length
            else:
                pass
        return r_en, q_st, q_en

    def __match(self):
        mlen, blen = 0, 0
        for l, o in self.cigar:
            if o in [0, 1]:
                mlen += l
            if o in [0, 1, 2]:
                blen += l
        return mlen, blen

    @property
    def cigar_string(self):
        return ''.join(['{}{}'.format(length, OPERATION[operation]) for length, operation in self.cigar])

    def __str__(self):
        return '\t'.join([str(x) for x in [self.q_st, self.q_en, self.ctg, self.r_st, self.r_en, self.mlen, self.blen,
                                           self.cigar_string]])


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


class Fasta(object):
    """
    load fasta file into memory
    """
    def __init__(self, infile):
        faidx = Faidx(infile)
        self.contig_len = faidx.contig_len
        self.genome = {ctg: faidx.seq(ctg, 0, size) for ctg, size in self.contig_len.items()}
        faidx.close()

    def seq(self, contig, start, end):
        if contig not in self.genome:
            return None
        return self.genome[contig][start:end]


def index_annotation(gtf):
    """
    Generate binned index for element in gtf
    """
    from CIRI_long.utils import tree

    LOGGER.info('Loading annotation gtf ..')
    gtf_index = defaultdict(dict)
    intron_index = defaultdict(dict)
    splice_site_index = tree()

    last_exon = None
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
                splice_site_index[parser.contig][parser.start][parser.strand]['start'] = 1
                splice_site_index[parser.contig][parser.end][parser.strand]['end'] = 1

                # Load intron
                if last_exon is not None and last_exon.attr['transcript_id'] == parser.attr['transcript_id']:
                    intron_start = last_exon.end if last_exon.strand == '+' else last_exon.start
                    intron_end = parser.start if parser.strand == '+' else parser.end
                    intron_strand = parser.strand

                    intron_start, intron_end = min(intron_start, intron_end), max(intron_start, intron_end)
                    start_div, end_div = intron_start // 500, intron_end // 500
                    for i in range(start_div, end_div + 1):
                        intron_index[parser.contig].setdefault(i, []).append((intron_start, intron_end, intron_strand))

                last_exon = parser

            # Binned index
            start_div, end_div = parser.start // 500, parser.end // 500
            for i in range(start_div, end_div + 1):
                gtf_index[parser.contig].setdefault(i, []).append(parser)

    return gtf_index, intron_index, splice_site_index


def index_circ(circ_file, circ_ss_idx):
    """
    Generate binned index for element in gtf
    """
    from CIRI_long.utils import tree
    from pathlib import Path

    circ_path = Path(circ_file)
    if circ_ss_idx is None:
        circ_ss_idx = tree()

    if circ_path.suffix == ".gtf":
        LOGGER.info('Loading additional circRNA gtf ..')
        with open(circ_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                content = line.rstrip().split('\t')
                parser = GTFParser(content)
                circ_ss_idx[parser.contig][parser.start][parser.strand]['start'] = 1
                circ_ss_idx[parser.contig][parser.end][parser.strand]['end'] = 1
    elif circ_path.suffix == ".bed":
        LOGGER.info('Loading additional circRNA bed ..')
        n_skip = 0
        with open(circ_path, 'r') as f:
            for line in f:
                content = line.rstrip().split('\t')
                contig = content[0]
                try:
                    start = int(content[1])
                    end = int(content[2])
                except ValueError:
                    n_skip += 1
                    continue
                strand = content[3]
                circ_ss_idx[contig][start][strand]['start'] = 1
                circ_ss_idx[contig][end][strand]['end'] = 1
        LOGGER.warn('Skipping {} lines in bed file'.format(n_skip))
    else:
        sys.exit('{} is not a valid bed/gtf file'.format(str(circ_path)))

    return circ_ss_idx


def get_blocks(hit):
    """
    Get blocks of aligned segments
    :param hit:
    :return:
    """
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


def get_exons(hit):
    r_start, r_end = hit.r_st, hit.r_st
    q_start, q_end = hit.q_st, hit.q_st
    r_block = []
    for length, operation in hit.cigar:
        if operation == 0:
            r_end += length
            q_end += length
        elif operation == 1:
            q_end += length
        elif operation == 2:
            r_end += length
        elif operation == 3:
            r_block.append([r_start, r_end, q_start, q_end])
            r_start = r_end + length
            r_end = r_start
            q_start = q_end
        elif operation == 4:
            pass

    if r_end > r_start:
        r_block.append([r_start, r_end, q_start, q_end])
    else:
        pass

    return r_block


def get_parital_blocks(hit, junc):
    exons = get_exons(hit)
    blocks = []
    for r_st, r_en, q_st, q_en in exons:
        if abs(q_st - junc) <= 10:
            blocks.append([r_st, r_en, '*-'])
        elif abs(q_en - junc) <= 10:
            blocks.append([r_st, r_en, '-*'])
        else:
            blocks.append([r_st, r_en, r_en - r_st + 1])
    return blocks


def merge_blocks(blocks):
    from operator import itemgetter
    tmp = sorted(blocks, key=itemgetter(0, 1))
    merged = []
    last_st, last_en = tmp[0][0], tmp[0][1]
    for st, en, length in tmp[1:]:
        if st <= last_en:
            last_en = max(en, last_en)
            last_st = min(st, last_st)
        else:
            merged.append([last_st, last_en, last_en - last_st + 1])
            last_st, last_en = st, en
    merged.append([last_st, last_en, last_en - last_st + 1])
    return merged


def merge_exons(tail_exons, head_exons):
    if head_exons[0][0] < tail_exons[-1][1]:
        return merge_blocks(tail_exons + head_exons)
    else:
        head_exons[0] = [head_exons[0][0], head_exons[0][1], '*-']
        tail_exons[-1] = [tail_exons[-1][0], tail_exons[-1][1], '-*']
        return tail_exons + head_exons


def merge_clip_exon(exons, clip_info):
    clip_st, clip_en = clip_info[0], clip_info[1]
    exon_st, exon_en = exons[0][0], exons[-1][1]

    if clip_st and clip_en:
        if clip_en < exon_st:
            exons = [[clip_st, clip_en, clip_en - clip_st + 1], ] + exons
        elif exon_en < clip_st:
            exons = exons + [[clip_st, clip_en, clip_en - clip_st + 1], ]
        elif clip_st < exon_st < clip_en:
            exons[0] = [clip_st, exons[0][1], exons[0][1] - clip_st + 1]
        elif clip_st < exon_en < clip_en:
            exons[-1] = [exons[-1][0], clip_en, clip_en - exons[-1][0] + 1]
        else:
            pass

    return exons


def remove_long_insert(hit):
    r_st, q_st = hit.r_st, hit.q_st
    last_r_st, last_q_st = r_st, q_st
    last_cigar = []
    sub_hits = []
    for length, operation in hit.cigar:
        if operation == 0:
            r_st += length
            q_st += length
        elif operation == 1:
            q_st += length
            if length > 20:
                sub_hits.append(SubHit(hit, last_r_st, last_q_st, last_cigar))
                last_cigar = []
                last_r_st, last_q_st = r_st, q_st
                continue
        elif operation in [2, 3]:
            r_st += length
        elif operation in [4, 5]:
            if q_st == hit.q_st:
                q_st += length
        else:
            pass
        last_cigar.append((length, operation))
    if last_cigar:
        sub_hits.append(SubHit(hit, last_r_st, last_q_st, last_cigar))
    primary_hit = sorted(sub_hits, key=lambda x: x.mlen, reverse=True)[0]
    primary_hit.is_primary = 1

    return primary_hit


def get_primary_alignment(hits):
    if not hits:
        return None

    for hit in hits:
        if hit.is_primary:
            return remove_long_insert(hit)

    return None


def find_annotated_signal(contig, start, end, clip_base, search_length=10, shift_threshold=3):
    tmp_annotated_signal = {}

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
        return None, us_free, ds_free, tmp_annotated_signal

    # First: find annotated splice signal
    if env.SS_INDEX is not None and contig in env.SS_INDEX:
        anno_ss = []
        for strand in ['+', '-']:
            # Upstream
            tmp_us_sites = []
            for us_shift in range(-search_length, search_length):
                us_pos = start + us_shift + 1
                if us_pos not in env.SS_INDEX[contig]:
                    continue
                if strand not in env.SS_INDEX[contig][us_pos]:
                    continue
                if 'start' not in env.SS_INDEX[contig][us_pos][strand]:
                    continue
                tmp_us_sites.append(us_shift)

            for us_shift in range(-search_length, search_length):
                us_pos = start + us_shift
                if us_pos not in env.SS_INDEX[contig]:
                    continue
                if strand not in env.SS_INDEX[contig][us_pos]:
                    continue
                if 'end' in env.SS_INDEX[contig][us_pos][strand]:
                    tmp_us_sites.append(us_shift)

            # Downstream
            tmp_ds_sites = []
            for ds_shift in range(-search_length, search_length):
                ds_pos = end + ds_shift + 1
                if ds_pos not in env.SS_INDEX[contig]:
                    continue
                if strand not in env.SS_INDEX[contig][ds_pos]:
                    continue
                if 'start' not in env.SS_INDEX[contig][ds_pos][strand]:
                    continue
                tmp_ds_sites.append(ds_shift)

            for ds_shift in range(-search_length, search_length):
                ds_pos = end + ds_shift
                if ds_pos not in env.SS_INDEX[contig]:
                    continue
                if strand not in env.SS_INDEX[contig][ds_pos]:
                    continue
                if 'end' not in env.SS_INDEX[contig][ds_pos][strand]:
                    continue
                tmp_ds_sites.append(ds_shift)

            tmp_annotated_signal[strand] = (tmp_us_sites, tmp_ds_sites)

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
            return sort_ss(anno_ss, us_free, ds_free, clip_base), us_free, ds_free, tmp_annotated_signal

    return None, us_free, ds_free, tmp_annotated_signal


def find_denovo_signal(contig, start, end, host_strand, tmp_signal, us_free, ds_free, clip_base,
                       search_length=10, shift_threshold=3, is_canonical=False):
    # Second: find GT-AG splice signal
    us_search_length = search_length + us_free
    ds_search_length = search_length + ds_free
    us_seq = env.GENOME.seq(contig, start - us_search_length - 2, start + ds_search_length)
    ds_seq = env.GENOME.seq(contig, end - us_search_length, end + ds_search_length + 2)

    if us_seq is None or len(us_seq) < ds_search_length - us_search_length + 2:
        return None
    if ds_seq is None or len(ds_seq) < ds_search_length - us_search_length + 2:
        return None

    # Same strand
    prior_ss = []
    if host_strand:
        prior_strand = set(list(host_strand))
        for strand in prior_strand:
            for (tmp_ds_ss, tmp_us_ss), ss_weight in SPLICE_SIGNAL.items():
                # Only search canonical GT-AG
                if is_canonical and ss_weight != 0:
                    continue
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
                    tmp_us_sites.append(tmp_us - us_search_length)
                    tmp_us_start = tmp_us

                # Find downstream signal
                tmp_ds_start = 0
                tmp_ds_sites = []
                while 1:
                    tmp_ds = ds_seq.find(ds_ss, tmp_ds_start + 1)
                    if tmp_ds == -1:
                        break
                    tmp_ds_sites.append(tmp_ds - us_search_length)
                    tmp_ds_start = tmp_ds

                if strand in tmp_signal:
                    tmp_us_signal, tmp_ds_signal = tmp_signal[strand]
                    tmp_us_sites = sorted(list(set(tmp_us_sites + tmp_us_signal)))
                    tmp_ds_sites = sorted(list(set(tmp_ds_sites + tmp_ds_signal)))

                # Filter paired splice signal in concordance position
                if len(tmp_us_sites) == 0 or len(tmp_ds_sites) == 0:
                    continue

                for i in tmp_us_sites:
                    for j in tmp_ds_sites:
                        if abs(i - j) > clip_base + shift_threshold:
                            continue
                        ss_id = '{}-{}*|{}-{}'.format(tmp_us_ss, tmp_ds_ss, i, j)
                        prior_ss.append((
                            ss_id, strand, i, j, ss_weight,
                            *get_ss_altered_length(i, j, us_free, ds_free, clip_base)
                        ))

    if len(prior_ss) > 0:
        return sort_ss(prior_ss, us_free, ds_free, clip_base)

    # Anti-sense
    other_ss = []
    other_strand = {'+', '-'} - set(list(host_strand)) if host_strand else {'+', '-'}
    if other_strand:
        for strand in other_strand:
            for (tmp_ds_ss, tmp_us_ss), ss_weight in SPLICE_SIGNAL.items():
                if is_canonical and ss_weight != 0:
                    continue
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
                    tmp_us_sites.append(tmp_us - us_search_length)
                    tmp_us_start = tmp_us

                # Find downstream signal
                tmp_ds_start = 0
                tmp_ds_sites = []
                while 1:
                    tmp_ds = ds_seq.find(ds_ss, tmp_ds_start + 1)
                    if tmp_ds == -1:
                        break
                    tmp_ds_sites.append(tmp_ds - us_search_length)
                    tmp_ds_start = tmp_ds

                if strand in tmp_signal:
                    tmp_us_signal, tmp_ds_signal = tmp_signal[strand]
                    tmp_us_sites = sorted(list(set(tmp_us_sites + tmp_us_signal)))
                    tmp_ds_sites = sorted(list(set(tmp_ds_sites + tmp_ds_signal)))

                # Filter paired splice signal in concordance position
                if len(tmp_us_sites) == 0 or len(tmp_ds_sites) == 0:
                    continue

                for i in tmp_us_sites:
                    for j in tmp_ds_sites:
                        if abs(i - j) > clip_base + shift_threshold:
                            continue
                        ss_id = '{}-{}*|{}-{}'.format(tmp_us_ss, tmp_ds_ss, i, j)
                        other_ss.append((
                            ss_id, strand, i, j, ss_weight,
                            *get_ss_altered_length(i, j, us_free, ds_free, clip_base)
                        ))

    if len(other_ss) > 0:
        return sort_ss(other_ss, us_free, ds_free, clip_base)

    return None


def get_ss_altered_length(i, j, us_free, ds_free, clip_base):
    clip_altered = min(abs(j - i - clip_base), abs(j - i + clip_base))
    us_altered = min(abs(i + us_free), abs(i - ds_free))
    ds_altered = min(abs(j + us_free), abs(j - ds_free))
    return abs(i-j), clip_altered, us_altered + ds_altered


def sort_ss(sites, us, ds, clip_base):
    # Splice site: site_id, strand, us_shift, ds_shift, site_weight, altered_len, clip_altered, altered_total
    from operator import itemgetter
    get_ss = itemgetter(0, 1, 2, 3)

    tmp_sites = set(sites)

    # Clipped sites
    clipped_sites = [i for i in tmp_sites if (-clip_base <= i[2] - i[3] <= clip_base)]
    if len(clipped_sites) > 0:
        return get_ss(sorted(clipped_sites, key=itemgetter(6, 5, 4, 7))[0])
    tmp_sites = set(sites) - set(clipped_sites)

    # Confidential splice sites
    confident_sites = [i for i in tmp_sites if -us <= i[2] <= ds and -us <= i[3] <= ds]
    if len(confident_sites) > 0:
        return get_ss(sorted(confident_sites, key=itemgetter(5, 4, 6, 7))[0])
    tmp_sites = tmp_sites - set(confident_sites)

    # Ambiguous splice sites
    ambiguous_sites = [i for i in tmp_sites if -clip_base <= i[2] <= 0 <= i[3] <= clip_base]
    if len(ambiguous_sites) > 0:
        return get_ss(sorted(ambiguous_sites, key=itemgetter(4, 5, 6, 7))[0])
    tmp_sites = tmp_sites - set(ambiguous_sites)

    # Other sites
    if len(tmp_sites) > 0:
        return get_ss(sorted(tmp_sites, key=itemgetter(4, 5, 6, 7))[0])
    return None


def find_host_gene(ctg, start, end):
    if ctg not in env.GTF_INDEX:
        return None
    start_div, end_div = start // 500, end // 500

    host_gene = {}
    for x in range(start_div, end_div + 1):
        if x not in env.GTF_INDEX[ctg]:
            continue
        for element in env.GTF_INDEX[ctg][x]:
            # start site
            if element.end < start or element.start > end:
                continue
            if element.start - 500 <= start <= element.end + 500 or element.start - 500 <= end <= element.end + 500:
                host_gene.setdefault(element.strand, []).append(element)

    if host_gene:
        return host_gene
    else:
        return None


def find_retained_introns(ctg, start, end):
    if ctg not in env.INTRON_INDEX:
        return None
    start_div, end_div = start // 500, end // 500

    host_gene = {}
    for x in range(start_div, end_div + 1):
        if x not in env.INTRON_INDEX[ctg]:
            continue
        for st, en, strand in env.INTRON_INDEX[ctg][x]:
            if st - 25 <= start and end <= en + 25:
                host_gene.setdefault(strand, []).append((st, en, strand))

    if host_gene:
        return host_gene
    else:
        return None


def find_overlap_exons(ctg, start, end):
    if ctg not in env.GTF_INDEX:
        return None
    start_div, end_div = start // 500, end // 500

    host_gene = {}
    for x in range(start_div, end_div + 1):
        if x not in env.GTF_INDEX[ctg]:
            continue
        for element in env.GTF_INDEX[ctg][x]:
            if element.type != 'exon':
                continue
            if element.end - 25 < start or end < element.start + 25:
                continue
            host_gene.setdefault(element.strand, []).append((element.start, element.end, element.strand))

    if host_gene:
        return host_gene
    else:
        return None


def convert_cigar_string(x):
    return [(int(l), OPERATION[op]) for l, op in re.findall(r'(\d+)([MIDNSHP=X])', x)]


def find_alignment_pos(alignment, pos):
    r_st, r_en = alignment.ref_begin, alignment.ref_begin
    q_st, q_en = alignment.query_begin, alignment.query_begin
    for l, op in convert_cigar_string(alignment.cigar_string):
        if op == 0:
            r_en += l
            q_en += l
        elif op == 1:
            q_en += l
        elif op in [2, ]:
            r_en += l
        elif op in [4, 5]:
            pass
        if r_st <= pos <= r_en:
            return q_st + pos - r_st
        r_st = r_en
        q_st = q_en
    return None
