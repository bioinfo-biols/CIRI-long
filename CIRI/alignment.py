import re
import logging
from collections import defaultdict

import numpy as np
import mappy as mp

from CIRI import poagraph
from CIRI import seqgraphalignment
from CIRI.logger import ProgressBar

LOGGER = logging.getLogger('CIRI-long')
ALIGNER = None
SS_INDEX = None

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

SPLICE_SIGNAL = {
    ('GT', 'AG'): 0,  # U2-type
    ('GC', 'AG'): 1,  # U2-type
    ('AT', 'AC'): 1,  # U12-type
    ('GA', 'AG'): 2,  # non-canonical
    ('GT', 'TG'): 2,  # non-canonical
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

            if parser.attr('gene_biotype') is None or parser.attr('gene_biotype') in ['lincRNA', 'pseudogene']:
                continue
            if parser.attr('gene_type') is None or parser.attr('gene_type') in ['lincRNA', 'pseudogene']:
                continue

            # Binned index
            start_div, end_div = parser.start // 500, parser.end // 500
            for i in range(start_div, end_div + 1):
                gtf_index[parser.contig].setdefault(i, []).append(parser)

    return gtf_index, splice_site_index


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
            r_block.append([r_start, r_end, r_end - r_start + 1])
            r_start = r_end + length
            r_end = r_start
        elif operation == 4:
            pass
    if r_end > r_start:
        r_block.append([r_start, r_end, r_end - r_start + 1])
    # return ','.join(['{}-{}|{}'.format(i[0], i[1], i[2]) for i in r_block])
    return r_block


def search_splice_signal(contig, start, end, search_length=10, shift_threshold=3):
    from operator import itemgetter
    us_seq = ALIGNER.seq(contig, start - search_length - 2, start + search_length)
    ds_seq = ALIGNER.seq(contig, end - search_length, end + search_length + 2)

    get_ss = itemgetter(0, 1, 2, 3)
    # First: Find flanking junction from annotation gtf
    anno_ss = []
    for strand in ['+', '-']:
        tmp_us_sites = []
        for us_shift in range(-search_length, search_length):
            us_pos = start + us_shift
            if contig in SS_INDEX and us_pos in SS_INDEX[contig] and strand in SS_INDEX[contig][us_pos]:
                tmp_us_sites.append(us_shift)

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
                us_ss = ALIGNER.seq(contig, start + i - 2 - 1, start + i - 1)
                ds_ss = ALIGNER.seq(contig, end + j, end + j + 2)
                if strand == '-':
                    us_ss, ds_ss = mp.revcomp(ds_ss), mp.revcomp(us_ss)
                ss_id = '{}-{}|{}-{}'.format(us_ss, ds_ss, i, j)
                ss_weight = SPLICE_SIGNAL[(ds_ss, us_ss)] if (ds_ss, us_ss) in SPLICE_SIGNAL else 3
                anno_ss.append((ss_id, strand, i, j, ss_weight, abs(i - j), abs(i) + abs(j)))

    if len(anno_ss) > 0:
        return get_ss(sorted(anno_ss, key=itemgetter(4, 5, 6))[0])

    # Second: Find Denovo BSJ using pre-defined splice signal
    putative_ss = []
    for strand in ['+', '-']:
        for (ds_ss, us_ss), ss_weight in SPLICE_SIGNAL.items():
            ss_id = '{}-{}*'.format(us_ss, ds_ss)
            if strand == '-':
                ds_ss, us_ss = mp.revcomp(us_ss), mp.revcomp(ds_ss)

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
                    us_shift = i - search_length
                    ds_shift = j - search_length
                    putative_ss.append(('{}|{}-{}'.format(ss_id, us_shift, ds_shift), strand, us_shift,
                                        ds_shift, ss_weight, abs(i - j), abs(us_shift) + abs(ds_shift)))

    if len(putative_ss) == 0:
        return None

    return get_ss(sorted(putative_ss, key=itemgetter(4, 5, 6))[0])


def initializer(aligner, splice_site_index):
    global ALIGNER, SPLICE_SIGNAL, SS_INDEX
    ALIGNER = aligner
    SS_INDEX = splice_site_index


def parse_chunk(chunk, is_canonical):
    from CIRI.preprocess import revcomp
    from Levenshtein import distance

    consensus_cnt = 0
    raw_unmapped_cnt = 0
    ccs_mapped_cnt = 0
    accordance_cnt = 0
    bsj_cnt = 0

    ret = []
    ccs_dis = []
    for read_id, segments, ccs, raw in chunk:
        # Filter 1 - Remove ccs with strange consensus length
        fasta = []
        for pos in segments.split(';'):
            i_st, i_en = pos.split('-')
            fasta.append(('{}-{}'.format(i_st, i_en), raw[int(i_st):int(i_en)]))

        if len(fasta) < 2:
            continue

        d_mean = np.mean([len(i) for _, i in fasta[:-1]])
        d_delta = int(2.3 * np.sqrt(0.1 * d_mean))
        if 0.9 * (d_mean - d_delta) > len(ccs) or 1.1 * (d_mean + d_delta) < len(ccs):
            continue

        tmp_dis = []
        for i in fasta[:-1]:
            tmp_dis.append(distance(i[1], ccs) / len(i[1]))
        tmp_dis.append(distance(i[1], ccs[:len(i[1])]) / len(i[1]))
        if np.max(tmp_dis) > 0.1:
            continue

        fasta_msa = msa(fasta)
        if len(fasta_msa) > 1:
            continue

        consensus_cnt += 1

        # Filter 2 - Remove linear mapped reads
        is_raw_mapped = 0
        hit = None
        for hit in ALIGNER.map(raw):
            if not hit.is_primary:
                continue

            # Filter out linear reads
            if hit.q_en - hit.q_st > max(len(raw) * 0.8, len(raw) - 50):
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
        if seg_st > 15 or seg_en < len(raw) - 15:
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

        # CCS aligned more than 90%
        tmp_hits = []
        for hit in ALIGNER.map(tmp):
            if hit.is_primary:
                tmp_hits.append(hit)
        if len(tmp_hits) == 0:
            continue
        tmp_hit = sorted(tmp_hits, key=lambda x: x.mlen, reverse=True)[0]
        if tmp_hit.r_en - tmp_hit.r_st < 0.9 * len(ccs):
            continue
        bsj_cnt += 1

        # Retrive circRNA positions, convert minimap2 position to real position
        ss_site = search_splice_signal(tmp_hit.ctg, tmp_hit.r_st, tmp_hit.r_en)
        if ss_site is None:
            continue
        ss_id, strand, us_shift, ds_shift = ss_site

        # if is_canonical: keep canonical splice site only
        ss = ss_id.split('|')[0]
        if is_canonical and ss[-1] == '*' and ss != 'AG-GT*':
            continue

        circ_ctg = tmp_hit.ctg
        circ_start = tmp_hit.r_st + us_shift
        circ_end = tmp_hit.r_en + ds_shift
        circ_id = '{}:{}-{}'.format(circ_ctg, circ_start, circ_end)

        # Get Cirexons
        cir_exons = get_blocks(tmp_hit)
        cir_exons[0][0] = circ_start - 1 # Trim 1bp for correct boundary of minimap2 alignment
        cir_exons[-1][-1] = circ_end

        cir_exon_tag = []
        for cir_exon_start, cir_exon_end, _ in cir_exons:
            cir_exon_tag.append('{}-{}'.format(cir_exon_start + 1, cir_exon_end))
            # cir_exon_ss = search_splice_signal(tmp_hit.ctg, cir_exon_start, cir_exon_end)
            # if cir_exon_ss is None:
            #     cir_exon_tag.append('{}-{}|?'.format(cir_exon_start, cir_exon_end))
            # else:
            #     ss_id, ss_strand, us_shift, ds_shift = cir_exon_ss
            #     cir_exon_tag.append('{}-{}|{}'.format(cir_exon_start))
            # cir_exon_tag.append('{}-{}|{}|{}'.format(cir_exon_start, cir_exon_end, strand, ss_id))

        circ_seq = tmp if strand == '+' else revcomp(tmp)

        ret.append((read_id, circ_id, strand, ','.join(cir_exon_tag), ss_id, circ_seq))

    from scipy.stats import describe
    print(describe(ccs_dis))

    return consensus_cnt, raw_unmapped_cnt, ccs_mapped_cnt, accordance_cnt, bsj_cnt, ret


def filter_ccs_reads(ccs_seq, minimap_index, ss_indx, is_canonical, out_dir, prefix, threads, debugging):
    from collections import namedtuple
    from multiprocessing import Pool
    from CIRI.utils import grouper
    # Read = namedtuple('Read', 'segments ccs raw')
    reads_count = {
        'consensus': 0,
        'raw_unmapped': 0,
        'ccs_mapped': 0,
        'accordance': 0,
        'bsj': 0,
        'splice_signal': 0,
    }

    minimap_aligner = mp.Aligner(minimap_index, preset='splice')

    chunk_size = 250
    chunk_cnt = 0
    jobs = []
    pool = Pool(threads, initializer, (minimap_aligner, ss_indx))
    for reads in grouper(list(ccs_seq), chunk_size):
        chunk = [[i, ] + ccs_seq[i] for i in reads if i is not None]
        chunk_cnt += 1
        jobs.append(pool.apply_async(parse_chunk, (chunk, is_canonical)))
    pool.close()

    prog = ProgressBar()
    finished_cnt = 0

    with open('{}/{}.cand_circ.fa'.format(out_dir, prefix), 'w') as out:
        for job in jobs:
            finished_cnt += 1

            tmp_consensus, tmp_raw_unmapped, tmp_ccs_mapped, tmp_accordance, tmp_bsj, ret = job.get()
            reads_count['consensus'] += tmp_consensus
            reads_count['raw_unmapped'] += tmp_raw_unmapped
            reads_count['ccs_mapped'] += tmp_ccs_mapped
            reads_count['accordance'] += tmp_accordance
            reads_count['bsj'] += tmp_bsj
            reads_count['splice_signal'] += len(ret)

            for read_id, circ_id, strand, cir_exon_tag, ss_id, circ_seq in ret:
                out.write('>{}\t{}\t{}\t{}\t{}\n{}\n'.format(read_id, circ_id, strand, cir_exon_tag, ss_id, circ_seq))
            prog.update(100 * finished_cnt / chunk_cnt)
    pool.join()
    prog.update(100)

    return reads_count
