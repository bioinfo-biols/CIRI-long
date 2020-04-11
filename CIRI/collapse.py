import logging
import pandas as pd
import numpy as np
from multiprocessing import Pool
from collections import defaultdict, Counter, namedtuple

from CIRI.utils import tree, grouper
from CIRI.logger import ProgressBar
from CIRI.alignment import SPLICE_SIGNAL

LOGGER = logging.getLogger('CIRI-long')
READ = namedtuple('Read', 'read_id circ_id strand cirexon ss clip segments seq sample type')
globals()['Read'] = READ


def load_cand_circ(in_file):
    from pathlib import Path
    sample_attr = {}
    with open(in_file, 'r') as f:
        for line in f:
            sample, fname = line.rstrip().split('\t')
            sample_attr[sample] = fname

    cand_reads = {}
    for sample, fname in sample_attr.items():
        # Full length reads
        cand_circ = Path(fname)
        with open(cand_circ, 'r') as f:
            for line in f:
                content = line.rstrip().lstrip('>').split('\t')
                clip_base = int(content[5].split('|')[1].split('-')[0])
                seq = f.readline().rstrip()
                if clip_base > 20:
                    continue
                cand_reads[content[0]] = READ(*content, seq, sample, 'full')

        # Partial reads
        prefix = cand_circ.name.split('.')[0]
        with open(cand_circ.parent / (prefix + '.low_confidence.fa')) as f:
            for line in f:
                content = line.rstrip().lstrip('>').split('\t')
                clip_base = int(content[5].split('|')[1].split('-')[0])
                seq = f.readline().rstrip()
                if clip_base > 20:
                    continue
                cand_reads[content[0]] = READ(*content, seq, sample, 'partial')

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


def junction_seq(seq, bsj, width=25):
    st, en = bsj - width, bsj + width
    if len(seq) <= 2 * width:
        return seq[bsj-len(seq)//2:] + seq[:bsj-len(seq)//2]

    if st < 0:
        if en < 0:
            return seq[st:en]
        else:
            return seq[st:] + seq[:en]
    elif en > len(seq):
        return seq[st:] + seq[:en-len(seq)]
    else:
        return seq[st:en]


CONTIG_LEN = None
GENOME = None
SS_INDEX = None


def initializer(genome, ss_index, contig_len):
    global CONTIG_LEN, GENOME, SS_INDEX
    CONTIG_LEN = contig_len
    GENOME = genome
    SS_INDEX = ss_index


def genome_junction_seq(contig, start, end, width=25):
    return GENOME[contig][end - width:end] + GENOME[contig][start:start + width]


def avg_score(alignment, ref, query):
    from Levenshtein import distance
    x = query[alignment.query_begin:alignment.query_end]
    return distance(ref, x) / len(ref)


def curate_junction(ctg, st, en, junc):
    from libs.striped_smith_waterman.ssw_wrap import Aligner
    from operator import itemgetter
    junc_scores = []
    for i in range(min(st) - 10, max(st) + 10):
        for j in range(min(en) - 10, max(en) + 10):
            tmp = genome_junction_seq(ctg, i, j, width=10)
            tmp_aligner = Aligner(tmp, match=10, mismatch=4, gap_open=8, gap_extend=2)
            tmp_scores = avg_score(tmp_aligner.align(junc), tmp, junc)
            junc_scores.append((i, j, np.mean(tmp_scores)))
    return sorted(junc_scores, key=itemgetter(2))


def equivalent_seq(genome, contig_len, contig, start, end, strand):
    from CIRI.preprocess import revcomp

    if strand is None:
        return 'Unknown'

    ds_seq = ''
    for i in range(100):
        if end + i > contig_len[contig]:
            break
        if genome[contig][start:start + i] == genome[contig][end:end + i]:
            ds_seq += genome[contig][start:start + i]
        else:
            break

    us_seq = ''
    for j in range(100):
        if start - j < 0:
            break
        if genome[contig][start - j:start] == genome[contig][end - j:end]:
            us_seq += genome[contig][start - j:start]
        else:
            break

    tmp = us_seq[::-1] + ds_seq
    if strand == '+':
        return tmp
    else:
        return revcomp(tmp)


def strict_splice_signal(contig, start, end, clip_base, search_length=10, shift_threshold=3):
    from CIRI.preprocess import revcomp
    from CIRI.alignment import ss_altered_length, sort_splice_sites

    ds_free = 0
    for i in range(100):
        if end + i > CONTIG_LEN[contig]:
            break
        if GENOME[contig][start:start + i] == GENOME[contig][end:end + i]:
            ds_free = i
        else:
            break

    us_free = 0
    for j in range(100):
        if start - j < 0:
            break
        if GENOME[contig][start - j:start] == GENOME[contig][end - j:end]:
            us_free = j
        else:
            break

    if start - search_length - us_free - 2 < 0 or end + search_length + ds_free + 2 > CONTIG_LEN[contig]:
        return None, us_free, ds_free
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
                    if abs(i - j) > shift_threshold + clip_base:
                        continue
                    us_ss = GENOME[contig][start + i - 2:start + i]
                    ds_ss = GENOME[contig][end + j:end + j + 2]
                    if strand == '-':
                        us_ss, ds_ss = revcomp(ds_ss), revcomp(us_ss)
                    ss_id = '{}-{}|{}-{}'.format(us_ss, ds_ss, i, j)
                    ss_weight = SPLICE_SIGNAL[(ds_ss, us_ss)] if (ds_ss, us_ss) in SPLICE_SIGNAL else 3

                    anno_ss.append((
                        ss_id, strand, i, j, ss_weight, *ss_altered_length(i, j, us_free, ds_free, clip_base)
                    ))

        if len(anno_ss) > 0:
            return sort_splice_sites(anno_ss, us_free, ds_free, clip_base), us_free, ds_free

    # Second: Find Denovo BSJ using pre-defined splice signal
    us_search_length = search_length + us_free
    ds_search_length = search_length + ds_free
    us_seq = GENOME[contig][start - us_search_length - 2:start + ds_search_length]
    ds_seq = GENOME[contig][end - us_search_length:end + ds_search_length + 2]

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
                        *ss_altered_length(us_shift, ds_shift, us_free, ds_free, clip_base)
                    ))

    if len(putative_ss) > 0:
        return sort_splice_sites(putative_ss, us_free, ds_free, clip_base), us_free, ds_free

    return None, us_free, ds_free


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
        ssw = Aligner(ref.seq[:50], match=10, mismatch=4, gap_open=8, gap_extend=2)

        head_pos = []
        for query in cluster[1:]:
            alignment = ssw.align(query.seq[:50])
            head_pos.append(alignment.ref_begin)

        template = transform_seq(ref.seq, max(head_pos))
        ssw = Aligner(template, match=10, mismatch=4, gap_open=8, gap_extend=2)
        junc_seqs = [junction_seq(template, -max(head_pos)//2, 25), ]

        for query in cluster[1:]:
            alignment = ssw.align(query.seq)
            tmp = transform_seq(query.seq, alignment.query_begin)
            junc_seqs.append(junction_seq(tmp, -max(head_pos)//2, 25))

        cs_junc = consensus(junc_seqs, alignment_type=2,
                            match=10, mismatch=-4, gap=-8, extension=-2, gap_affine=-24, extension_affine=-1,
                            debug=0)

        ctg = Counter([i.circ_id.split(':')[0] for i in cluster]).most_common()[0][0]
        tmp_st = [int(i.circ_id.split(':')[1].split('-')[0]) for i in cluster]
        tmp_en = [int(i.circ_id.split(':')[1].split('-')[1]) for i in cluster]

        scores = curate_junction(ctg, tmp_st, tmp_en, cs_junc)
        circ_start, circ_end, score = scores[0]

        for shift_threshold in [0, 3]:
            ss_site, us_free, ds_free = strict_splice_signal(ctg, circ_start, circ_end, 0, 5, shift_threshold)
            if ss_site is not None:
                ss_id, strand, us_shift, ds_shift = ss_site
                circ_start += us_shift
                circ_end += ds_shift
                break

        if ss_site is None:
            ss_id = 'None'
            strand = 'None'

        circ_id = '{}:{}-{}'.format(ctg, circ_start + 1, circ_end)

        # from spoa import poa
        # _, msa = poa([genome_junction_seq(ctg, circ_start, circ_end, width=10), cs_junc] + junc_seqs, 2, True,
        #              10, -4, -8, -2, -24, -1)
        # print(circ_id, ss_id)
        cs_cluster.append(([i.read_id for i in cluster], circ_id, strand, ss_id, us_free, ds_free))

    return cs_cluster


def correct_reads(reads_cluster, ref_fasta, ss_index, threads):
    from CIRI.alignment import Faidx
    faidx = Faidx(ref_fasta)
    contig_len = faidx.contig_len
    genome = {ctg: faidx.seq(ctg, 0, size) for ctg, size in contig_len.items()}
    faidx.close()

    corrected_reads = []
    jobs = []
    pool = Pool(threads, initializer, (genome, ss_index, contig_len, ))
    for cluster in grouper(reads_cluster, 250):
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


def circ_pos(x):
    ctg, pos = x.split(':')
    st, en = pos.split('-')
    return ctg, int(st), int(en)


def by_circ(x):
    ctg, pos = x.split(':')
    if ctg.startswith('chr'):
        ctg = ctg.lstrip('chr')

    try:
        idx = '{:02d}'.format(int(ctg))
    except ValueError as e:
        if ctg in ['X', 'x', 'Y', 'y']:
            idx = 'a'
        elif ctg in ['M', 'm']:
            idx = 'b'
        else:
            idx = 'c'

    st, en = pos.split('-')

    return idx, ctg, int(st), int(en)


def cal_exp_mtx(cand_reads, corrected_reads, ref_fasta, gtf_idx, out_dir, prefix):
    from collections import Counter
    from CIRI.alignment import Faidx
    faidx = Faidx(ref_fasta)
    contig_len = faidx.contig_len
    genome = {ctg: faidx.seq(ctg, 0, size) for ctg, size in contig_len.items()}
    faidx.close()

    exp_df = {}
    circ_info = {}
    reads_df = []

    for reads, circ_id, strand, ss_id, us_free, ds_free in corrected_reads:
        # circRNA information
        ctg, st, en = circ_pos(circ_id)
        field = circ_attr(gtf_idx, ctg, st, en, strand)

        tmp_attr = 'circ_id "{}"; splice_site "{}"; equivalent_seq "{}"; circ_type "{}";'.format(
            circ_id,
            ss_id,
            equivalent_seq(genome, contig_len, ctg, st, en, strand),
            field['circ_type'] if field else 'Unknown',
        )
        for key in 'gene_id', 'gene_name', 'gene_type':
            if key in field:
                tmp_attr += ' {} "{}";'.format(key, field[key])
        circ_info[circ_id] = [ctg, 'CIRI-long', 'circRNA', st, en, len(reads), strand, '.', tmp_attr, ]

        # Expression levels
        exp_df[circ_id] = Counter([cand_reads[i].sample for i in reads])

        # Corrected reads
        for read_id in reads:
            read = cand_reads[read_id]
            tmp = [read_id, read.circ_id, read.strand, read.ss, read.clip, read.segments, read.sample, circ_id]
            reads_df.append(tmp)

    reads_df = pd.DataFrame(
        reads_df, columns=['read_id', 'circ_id', 'strand', 'signal', 'alignment', 'segments', 'sample', 'corrected'])
    reads_df.to_csv('{}/{}.reads'.format(out_dir, prefix), sep="\t", index=False)

    sorted_circ = sorted(list(circ_info), key=by_circ)
    with open('{}/{}.info'.format(out_dir, prefix), 'w') as out:
        for circ_id in sorted_circ:
            tmp_line = circ_info[circ_id]
            out.write('\t'.join([str(x) for x in tmp_line]) + '\n')

    exp_df = pd.DataFrame.from_dict(exp_df).transpose().fillna(0).reindex(sorted_circ)
    exp_df.to_csv('{}/{}.expression'.format(out_dir, prefix), sep="\t", index_label='circ_ID')

    return 0


def circ_attr(gtf_index, ctg, start, end, strand):
    """
    annotate circRNA information
    """
    if ctg not in gtf_index:
        # LOGGER.warn('chrom of contig "{}" not in annotation gtf, please check'.format(ctg))
        return {}
    start_div, end_div = start // 500, end // 500

    host_gene = {}
    start_element = defaultdict(list)
    end_element = defaultdict(list)

    for x in range(start_div, end_div + 1):
        if x not in gtf_index[ctg]:
            continue
        for element in gtf_index[ctg][x]:
            # start site
            if element.start <= start <= element.end and element.strand == strand:
                start_element[element.type].append(element)
            # end site
            if element.start <= end <= element.end and element.strand == strand:
                end_element[element.type].append(element)
            # if element.type != 'gene':
            #     continue
            if element.end < start or end < element.start:
                continue
            if element.attr['gene_id'] not in host_gene:
                host_gene[element.attr['gene_id']] = element

    circ_type = {}
    forward_host_gene = []
    antisense_host_gene = []

    if len(host_gene) > 0:
        for gene_id in host_gene:
            if host_gene[gene_id].strand == strand:
                forward_host_gene.append(host_gene[gene_id])
                if 'exon' in start_element and 'exon' in end_element:
                    circ_type['exon'] = 1
                else:
                    circ_type['intron'] = 1
            else:
                antisense_host_gene.append(host_gene[gene_id])
                circ_type['antisense'] = 1
    else:
        circ_type['intergenic'] = 1

    if len(forward_host_gene) > 1:
        circ_type['gene_intergenic'] = 1

    field = {}
    if 'exon' in circ_type:
        field['circ_type'] = 'exon'
    elif 'intron' in circ_type:
        field['circ_type'] = 'intron'
    # elif 'gene_intergenic' in circ_type:
    #     field['circ_type'] = 'gene_intergenic'
    elif 'antisense' in circ_type:
        field['circ_type'] = 'antisense'
    else:
        field['circ_type'] = 'intergenic'

    if len(forward_host_gene) == 1:
        if 'gene_id' in forward_host_gene[0].attr:
            field['gene_id'] = forward_host_gene[0].attr['gene_id']
        if 'gene_name' in forward_host_gene[0].attr:
            field['gene_name'] = forward_host_gene[0].attr['gene_name']
        if 'gene_type' in forward_host_gene[0].attr:
            field['gene_type'] = forward_host_gene[0].attr['gene_type']
        elif 'gene_biotype' in forward_host_gene[0].attr:
            field['gene_type'] = forward_host_gene[0].attr['gene_biotype']
        else:
            pass
    elif len(forward_host_gene) > 1:
        tmp_gene_id = []
        tmp_gene_name = []
        tmp_gene_type = []
        for x in forward_host_gene:
            if 'gene_id' in x.attr:
                tmp_gene_id.append(x.attr['gene_id'])
            if 'gene_name' in x.attr:
                tmp_gene_name.append(x.attr['gene_name'])
            if 'gene_type' in x.attr:
                tmp_gene_type.append(x.attr['gene_type'])
            elif 'gene_biotype' in x.attr:
                tmp_gene_type.append(x.attr['gene_biotype'])
            else:
                pass
        if len(tmp_gene_id) > 0:
            field['gene_id'] = ','.join(tmp_gene_id)
        if len(tmp_gene_name) > 0:
            field['gene_name'] = ','.join(tmp_gene_name)
        if len(tmp_gene_type) > 0:
            field['gene_type'] = ','.join(tmp_gene_type)
    elif field['circ_type'] == 'antisense' and len(antisense_host_gene) > 0:
        tmp_gene_id = []
        tmp_gene_name = []
        tmp_gene_type = []
        for x in antisense_host_gene:
            if 'gene_id' in x.attr:
                tmp_gene_id.append(x.attr['gene_id'])
            if 'gene_name' in x.attr:
                tmp_gene_name.append(x.attr['gene_name'])
            if 'gene_type' in x.attr:
                tmp_gene_type.append(x.attr['gene_type'])
            elif 'gene_biotype' in x.attr:
                tmp_gene_type.append(x.attr['gene_biotype'])
            else:
                pass
        if len(tmp_gene_id) > 0:
            field['gene_id'] = ','.join(tmp_gene_id)
        if len(tmp_gene_name) > 0:
            field['gene_name'] = ','.join(tmp_gene_name)
        if len(tmp_gene_type) > 0:
            field['gene_type'] = ','.join(tmp_gene_type)
    else:
        pass

    return field


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
        if circ is None:
            continue

        # Candidate alignment situation, more than 85%
        circ_hit = alignment.get_primary_alignment(alignment.ALIGNER.map(circ))
        if circ_hit is None or circ_hit.mlen < 0.75 * len(circ):
            continue

        clipped_circ, circ_start, circ_end, clip_info = alignment.align_clip_segments(circ, circ_hit)
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
            ret.append([reads, circ_id, circ_seq, 'Unknown'])
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
        circ_seq = clipped_circ if circ_hit.strand > 0 else revcomp(clipped_circ)
        circ_seq = circ_seq[correction_shift:] + circ_seq[:correction_shift]

        ret.append([reads, circ_id, circ_seq, ss_id])

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
    for reads in grouper(corrected_reads, chunk_size):
        chunk = [i for i in reads if i is not None]
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


def recover_corrected_chunk(chunk):
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

        clipped_circ, circ_start, circ_end, clip_info = alignment.align_clip_segments(circ, circ_hit)
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
            ret.append([reads, circ_id, circ_seq, 'Unknown'])
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
        circ_seq = clipped_circ if circ_hit.strand > 0 else revcomp(clipped_circ)
        circ_seq = circ_seq[correction_shift:] + circ_seq[:correction_shift]

        ret.append([reads, circ_id, circ_seq, ss_id])

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

