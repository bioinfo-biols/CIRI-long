import logging
import pandas as pd
import numpy as np
from multiprocessing import Pool
from collections import defaultdict, Counter, namedtuple

from CIRI import env
from CIRI.align import *
from CIRI.utils import *
from CIRI.logger import ProgressBar


LOGGER = logging.getLogger('CIRI-long')
READ = namedtuple('Read', 'read_id circ_id strand cirexon ss clip segments seq sample type')
CIRC = namedtuple('Circ', 'contig start end strand')
globals()['Read'] = READ
globals()['Circ'] = CIRC


class Segment(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.last = None
        self.next = None

    def __str__(self):
        return '{}-{}'.format(self.start, self.end)


class Exon(Segment):
    def __init__(self, start, end):
        self.start = int(start)
        self.end = int(end)


def load_cand_circ(in_file):
    from pathlib import Path
    sample_attr = {}
    with open(in_file, 'r') as f:
        for line in f:
            content = line.rstrip().split()
            if content:
                sample, fname = content
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


def cluster_reads(cand_reads):
    import re
    from operator import itemgetter

    circ_reads = defaultdict(list)
    circ_start = defaultdict(dict)
    circ_end = defaultdict(dict)

    # Load circRNA junction sites
    for read_id, read in cand_reads.items():
        contig, start, end = re.split('[:-]', read.circ_id)
        start, end = int(start), int(end)
        if end - start > 200 * 1000:
            continue

        # Store reads & BSJ sites
        circ_reads[contig].append((start, end, read.read_id))
        circ_start[contig].setdefault(start, []).append(read.read_id)
        circ_end[contig].setdefault(end, []).append(read.read_id)

    # Build index of circRNA junction sites
    reads_cluster = []
    for contig in circ_reads:
        circ_start_index = {}
        circ_end_index = {}

        # Cluster of start junction site
        tmp = [[], ]
        for x in sorted(circ_start[contig]):
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
        for x in sorted(circ_end[contig]):
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
        for (start, end, read_id) in sorted(circ_reads[contig], key=itemgetter(0, 1)):
            if read_id in reads_itered:
                continue

            tmp_reads = []
            p = [i for i in circ_start_index[start // 500] if start in i][0]
            q = [i for i in circ_end_index[end // 500] if end in i][0]

            for i in p:
                tmp_start = circ_start[contig][i]
                for j in q:
                    tmp_end = circ_end[contig][j]
                    tmp = set(tmp_start) & set(tmp_end)
                    if tmp:
                        tmp_reads += tmp

            for i in tmp_reads:
                reads_itered[i] = 1

            reads_cluster.append(sorted([cand_reads[i] for i in tmp_reads], key=lambda x: len(x.seq), reverse=True))

    return reads_cluster


def genome_junction_seq(contig, start, end, width=25):
    return env.GENOME.seq(contig, end-width, end) + env.GENOME.seq(contig, start, start + width)


def avg_score(alignment, ref, query):
    from Levenshtein import distance
    x = query[alignment.query_begin:alignment.query_end]
    return distance(ref, x) / len(ref)


def curate_junction(ctg, st, en, junc):
    from libs.striped_smith_waterman.ssw_wrap import Aligner
    from operator import itemgetter
    junc_scores = []
    for i in range(max(0, min(st) - 25), max(st) + 25):
        for j in range(min(en) - 25, min(max(en) + 25, env.CONTIG_LEN[ctg])):
            tmp = genome_junction_seq(ctg, i, j, width=10)
            tmp_aligner = Aligner(tmp, match=10, mismatch=4, gap_open=8, gap_extend=2)
            tmp_score = avg_score(tmp_aligner.align(junc), tmp, junc)
            junc_scores.append((i, j, tmp_score))
    return sorted(junc_scores, key=itemgetter(2))


def annotated_hit(contig, scores):
    if contig not in env.SS_INDEX:
        return None

    weighted = []
    for st, en, score in scores:
        w = 0
        if st + 1 in env.SS_INDEX[contig]:
            tmp = set(flatten([p for _, p in env.SS_INDEX[contig][st + 1].items()]))
            if 'start' in tmp:
                w += 1
        elif st in env.SS_INDEX[contig]:
            tmp = set(flatten([p for _, p in env.SS_INDEX[contig][st].items()]))
            if 'end' in tmp:
                w += 1
        else:
            pass

        if en in env.SS_INDEX[contig]:
            tmp = set(flatten([p for _, p in env.SS_INDEX[contig][en].items()]))
            if 'end' in tmp:
                w += 1
        elif en + 1 in env.SS_INDEX[contig]:
            tmp = set(flatten([p for _, p in env.SS_INDEX[contig][en + 1].items()]))
            if 'start' in tmp:
                w += 1
        else:
            pass

        weighted.append([st, en, w])

    return min_sorted_items(weighted, 2, True)


def correct_chunk(chunk):
    from collections import Counter
    from CIRI.poa import consensus
    from libs.striped_smith_waterman.ssw_wrap import Aligner

    cs_cluster = []
    cnt = defaultdict(int)
    for cluster in chunk:
        if cluster is None:
            continue
        if len(cluster) <= 1:
            continue

        if 'full' not in set([i.type for i in cluster]):
            continue

        # if '13a24928-f025-4a52-ae90-b7fa98c3b3ee' not in [i.read_id for i in cluster]:
        #     continue

        counter = Counter([i.circ_id for i in cluster if i.type == 'full']).most_common(n=1)
        ref = sorted([i for i in cluster if i.circ_id == counter[0][0] and i.type == 'full'],
                     key=lambda x: len(x.seq), reverse=True)[0]
        ssw = Aligner(ref.seq[:50], match=10, mismatch=4, gap_open=8, gap_extend=2)

        head_pos = []
        for query in cluster[1:]:
            alignment = ssw.align(query.seq[:50])
            head_pos.append(alignment.ref_begin)

        template = transform_seq(ref.seq, max(head_pos))
        ssw = Aligner(template, match=10, mismatch=4, gap_open=8, gap_extend=2)
        junc_seqs = [get_junc_seq(template, -max(head_pos)//2, 25), ]

        for query in cluster[1:]:
            alignment = ssw.align(query.seq)
            tmp = transform_seq(query.seq, alignment.query_begin)
            junc_seqs.append(get_junc_seq(tmp, -max(head_pos)//2, 25))

        cs_junc = consensus(junc_seqs, alignment_type=2,
                            match=10, mismatch=-4, gap=-8, extension=-2, gap_affine=-24, extension_affine=-1,
                            debug=0)

        ctg = Counter([i.circ_id.split(':')[0] for i in cluster]).most_common()[0][0]
        tmp_st = [int(i.circ_id.split(':')[1].split('-')[0]) for i in cluster]
        tmp_en = [int(i.circ_id.split(':')[1].split('-')[1]) for i in cluster]

        # Curate junction sequence
        scores = curate_junction(ctg, tmp_st, tmp_en, cs_junc)
        aval_junc = min_sorted_items(scores, 2)

        # Annotated sites
        anno_junc = annotated_hit(ctg, aval_junc)
        if anno_junc:
            circ_start, circ_end, circ_score = anno_junc[0]
        else:
            circ_start, circ_end, circ_score = aval_junc[0]

        # Annotated sites
        for shift_threshold in [5, 10]:
            ss_site, us_free, ds_free = find_annotated_signal(ctg, circ_start, circ_end, 0, 10, shift_threshold)
            if ss_site is not None:
                ss_id, strand, us_shift, ds_shift = ss_site
                circ_start += us_shift
                circ_end += ds_shift
                break

        host_strand = find_host_gene(ctg, circ_start, circ_end)

        # Canonical sites
        if ss_site is None:
            for shift_threshold in [5, 10]:
                ss_site = find_denovo_signal(ctg, circ_start, circ_end, host_strand, us_free, ds_free,
                                             0, 10, shift_threshold, True)
                if ss_site is not None:
                    ss_id, strand, us_shift, ds_shift = ss_site
                    circ_start += us_shift
                    circ_end += ds_shift
                    cnt['Annotated'] += 1
                    break

        # Intronic circRNAs
        if ss_site is None:
            retained_introns = find_retained_introns(ctg, circ_start + 1, circ_end)
            overlap_exons = find_overlap_exons(ctg, circ_start + 1, circ_end)

            is_lariat = 0
            if retained_introns is not None and overlap_exons is None:
                is_lariat = 1
                # Find high-confidence ciRNAs
                retained_introns = set(sum([i for _, i in retained_introns.items()], []))
                retained_strand = set([i[2] for i in retained_introns])
                tmp_circ = []
                for intron_start, intron_end, intron_strand in retained_introns:
                    if abs(intron_start - circ_start) > 50 or abs(intron_end - circ_end) > 50:
                        continue
                    if intron_strand == '+':
                        tmp_site = [i for i in scores if i[0] == intron_start]
                    else:
                        tmp_site = [i for i in scores if i[1] == intron_end]
                    if tmp_site:
                        tmp_circ.append([*tmp_site[0], intron_strand])

                ss_id = 'lariat'
                if tmp_circ:
                    circ_start, circ_end, circ_score, strand = sorted(tmp_circ, key=lambda x: x[2])[0]
                    cnt['High confidence lariat'] += 1
                else:
                    # Lariat with recursive splicing branchpoint
                    is_lariat = 0
                    tmp_circ = []
                    for tmp_strand in retained_strand:
                        tmp_start, tmp_end, tmp_score = recursive_splice_site(scores, ctg, tmp_strand)
                        if tmp_score is not None:
                            tmp_circ.append([tmp_start, tmp_end, tmp_score, tmp_strand])
                    if tmp_circ:
                        circ_start, circ_end, circ_score, strand = sorted(tmp_circ, key=lambda x: x[2])[0]
                        # cnt['Recursive splicing lariat'] += 1
                    else:
                        # cnt['Unknown lariat'] += 1
                        strand = 'None'

            # Find denovo splice signal
            if is_lariat == 0:
                ss_site = find_denovo_signal(ctg, circ_start, circ_end, host_strand, us_free, ds_free,
                                             5, 10, 3, False)
                if ss_site is not None:
                    ss_id, strand, us_shift, ds_shift = ss_site
                    circ_start += us_shift
                    circ_end += ds_shift
                    cnt['Denovo signal'] += 1
                else:
                    ss_id = 'None'
                    strand = 'None'
                    cnt['Unknown signal'] += 1

        circ_id = '{}:{}-{}'.format(ctg, circ_start + 1, circ_end)

        cluster_seq = []
        circ_junc_seq = genome_junction_seq(ctg, circ_start, circ_end)
        ssw = Aligner(circ_junc_seq, match=10, mismatch=4, gap_open=8, gap_extend=2, report_cigar=True)
        for query in cluster:
            # if query.type != 'full':
            #     continue
            alignment = ssw.align(query.seq * 2)
            tmp_pos = find_alignment_pos(alignment, len(circ_junc_seq)//2) % len(query.seq)
            if tmp_pos is None:
                cluster_seq.append(query.seq)
            else:
                tmp_seq = transform_seq(query.seq, tmp_pos)
                cluster_seq.append(tmp_seq)

        circ = CIRC(ctg, circ_start + 1, circ_end, strand)

        # Cluster sequence
        # seq_freq = [(i.seq, 1) for i in cluster if i.type == 'full']
        # seq_cluster = cluster_sequence(seq_freq)
        #
        # n_iter = 0
        # while 1:
        #     n_cluster = cluster_sequence(seq_cluster)
        #     if len(n_cluster) == len(seq_cluster):
        #         break
        #     seq_cluster = n_cluster
        #     n_iter += 1
        #     if n_iter >= 10:
        #         LOGGER.warn('Sequence not consensus for circRNA: {}'.format(circ_id))
        #         break

        isoforms, circ_len = curate_cirexons(circ, cluster)

        # tmp_fasta = sorted([i.seq for i in cluster if i.segments != 'partial'], key=len)
        # from CIRI.poa import consensus
        # tmp_consensus = consensus(tmp_fasta, 2, 10, -4, -8, -2, -24, -1, 0)

        # field = circ_attr(GTF_INDEX, ctg, circ_start + 1, circ_end, strand)
        # if field and field['circ_type'] == 'exon' and ss_id == 'lariat':
        #     print(circ_id)
        # host_gene = transcript_strand(GTF_INDEX, ctg, circ_start + 1, circ_end)
        # field = circ_attr(GTF_INDEX, ctg, circ_start + 1, circ_end, strand)
        # nearby_st, nearby_en = nearby_ss(ctg, circ_start, circ_end)
        # if nearby_st or nearby_en and field['circ_type'] != 'exon':
        #     print(circ_id)

        cs_cluster.append(([i.read_id for i in cluster], cluster_seq, circ_id, strand,
                           ss_id, us_free, ds_free, circ_len, isoforms))

    return cs_cluster, cnt


def batch_cluster_sequence(cluster_seq):
    reads = [i[0] for i in x]
    freq = [i[1] for i in x]

    res = []
    cnt = 0
    for chunk in grouper(x, 10):
        tmp = [i for i in chunk if i is not None]
        res = cluster_sequence(res + tmp)
        cnt += len(tmp)
        if sum([i[1] for i in res]) > cnt:
            print(res)
    return res


def cluster_sequence(segment_freq):
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
    from CIRI.poa import consensus
    from CIRI.find_ccs import compress_seq
    from Levenshtein import distance

    if len(segment_freq) == 1:
        return [(segment_freq[0][0], segment_freq[0][1]), ]

    hpc_seq = [compress_seq(i[0]) for i in segment_freq]

    # Calculate distance between segments
    dist = np.zeros((len(segment_freq), len(segment_freq)))
    for i in range(len(segment_freq)):
        for j in range(len(segment_freq)):
            if i > j:
                continue
            dist[i][j] = distance(hpc_seq[i], hpc_seq[j]) / max(len(hpc_seq[i]), len(hpc_seq[j]))

    dist = dist + dist.T
    # Hierarchical Cluster
    if dist.sum() != 0:
        z = leaves_list(linkage(squareform(dist), "ward", optimal_ordering=True))
    else:
        z = list(range(len(segment_freq)))

    clusters = [[z[0], ]]
    for i, j in pairwise(z):
        if i > j:
            if dist[j][i] < 0.3:
                clusters[-1].append(j)
            else:
                clusters.append([j, ])
        else:
            if dist[i][j] < 0.3:
                clusters[-1].append(j)
            else:
                clusters.append([j, ])

    ccs_seq = []
    for cluster in clusters:
        if len(cluster) == 1:
            ccs_seq.append((segment_freq[cluster[0]]))
            continue

        cluster_seq = [segment_freq[i][0] for i in cluster]
        cluster_freq = [segment_freq[i][1] for i in cluster]

        # Generate consensus sequence
        ccs = consensus(cluster_seq, 2, 10, -4, -8, -2, -24, -1, 0)
        ccs_seq.append((ccs, sum(cluster_freq)))
    return ccs_seq


def find_intron_signal(contig, starts, ends):
    if env.SS_INDEX is not None and contig in env.SS_INDEX:
        return None

    for strand in ['+', '-']:
        tmp_us_sites = []
        for st in range(min(starts)-10, max(starts)+10):
            us_pos = st
            if us_pos not in env.SS_INDEX[contig]:
                continue
            if strand not in env.SS_INDEX[contig][us_pos]:
                continue
            if 'end' not in env.SS_INDEX[contig][us_pos][strand]:
                continue
            tmp_us_sites.append(st)

        tmp_ds_sites = []
        for en in range(min(ends)-10, max(ends)+10):
            ds_pos = en
            if ds_pos not in env.SS_INDEX[contig]:
                continue
            if strand not in env.SS_INDEX[contig][ds_pos]:
                continue
            if 'start' not in env.SS_INDEX[contig][ds_pos][strand]:
                continue
            tmp_ds_sites.append(en)

        if len(tmp_us_sites) == 0 and len(tmp_ds_sites) == 0:
            continue
        elif len(tmp_us_sites) > 0 and len(tmp_ds_sites) > 0:
            pass
        elif len(tmp_us_sites) > 0:
            pass
        elif len(tmp_ds_sites) > 0:
            pass
        else:
            continue


def recursive_splice_site(scores, ctg, strand):
    for st, en, scr in scores:
        if strand == '+' and (env.GENOME.seq(ctg, st-2, st) == 'AG' and env.GENOME.seq(ctg, st, st+2) == 'GT'):
            return st, en, scr
        if strand == '-' and (env.GENOME.seq(ctg, en, en+2) == 'CT' and env.GENOME.seq(ctg, en-2, en) == 'CA'):
            return st, en, scr
    return None, None, None


def curate_cirexons(circ, cluster):
    isoforms = {}
    starts = []
    ends = []
    for read in cluster:
        if read.cirexon == 'NA':
            continue
        if read.type == 'partial':
            continue
        exons = parse_cirexons(circ, read)
        if len(exons) == 0:
            continue
        isoforms[read.read_id] = exons
        for exon in exons:
            starts.append(exon.start)
            ends.append(exon.end)
    if len(isoforms) == 0:
        return [], None

    # Cluster junction sites
    tmp_starts = cluster_bins(starts, dis=10)
    tmp_ends = cluster_bins(ends, dis=10)

    convert_st = {}
    for tmp_st in tmp_starts:
        # circRNA start
        if circ.start in tmp_st:
            for i in tmp_st:
                convert_st[i] = circ.start

        # Canonical sites
        aval_st = []
        for i in set(tmp_st):
            i_ss = env.GENOME.seq(circ.contig, i-3, i-1)
            if circ.strand == '+' and i_ss == 'AG':
                aval_st.append(i)
            elif circ.strand == '-' and revcomp(i_ss) == 'GT':
                aval_st.append(i)
            else:
                pass

        tmp_st_counter = Counter(tmp_st)
        if aval_st:
            final_st = sorted(aval_st, key=lambda x: tmp_st_counter[x], reverse=True)[0]
        else:
            final_st = tmp_st_counter.most_common(n=1)[0][0]
        for i in tmp_st:
            convert_st[i] = final_st

    convert_en = {}
    for tmp_en in tmp_ends:
        # circRNA start
        if circ.end in tmp_en:
            for i in tmp_en:
                convert_en[i] = circ.end

        # Canonical sites
        aval_en = []
        for i in set(tmp_en):
            i_ss = env.GENOME.seq(circ.contig, i, i + 1)
            if circ.strand == '-' and i_ss == 'AG':
                aval_en.append(i)
            elif circ.strand == '+' and revcomp(i_ss) == 'GT':
                aval_en.append(i)
            else:
                pass

        tmp_en_counter = Counter(tmp_en)
        if aval_en:
            final_en = sorted(aval_en, key=lambda x: tmp_en_counter[x], reverse=True)[0]
        else:
            final_en = tmp_en_counter.most_common(n=1)[0][0]
        for i in tmp_en:
            convert_en[i] = final_en

    curated_isoform = []
    iso_len = {}
    for read_id, exons in isoforms.items():
        tmp_iso = []
        tmp_len = 0
        for exon in exons:
            tmp_st = convert_st[exon.start]
            tmp_en = convert_en[exon.end]
            tmp_iso.append('{}-{}'.format(tmp_st, tmp_en))
            tmp_len += tmp_en - tmp_st + 1

        iso_name = ','.join(tmp_iso)
        curated_isoform.append(iso_name)
        iso_len[iso_name] = tmp_len

    final_len = iso_len[Counter(curated_isoform).most_common(n=1)[0][0]]

    return set(curated_isoform), final_len


def parse_cirexons(circ, read):
    exon_str = read.cirexon.split(',')
    exons = []
    for x in exon_str:
        st, en = x.split('|')[0].split('-')
        exons.append(Exon(st, en))
    if len(exons) == 0:
        return []

    if exons[0].end <= circ.start:
        exons = exons[1:]
    else:
        exons[0].start = circ.start
    if len(exons) == 0:
        return []

    if exons[-1].start >= circ.end:
        exons = exons[:-1]
    else:
        exons[-1].end = circ.end

    return exons


def cluster_bins(pos, dis=10):
    clustered = []
    last_i = None
    for i in sorted(pos):
        if last_i is None:
            last_i = [i, ]
            continue
        if i > last_i[-1] + dis:
            clustered.append(last_i)
            last_i = [i, ]
        else:
            last_i.append(i)
    clustered.append(last_i)
    return clustered


def intron_splice_sites(contig, start, end, strand):
    us = env.GENOME.seq(contig, start, start + 2)
    ds = env.GENOME.seq(contig, end - 3, end - 1)
    if strand == '-':
        us, ds = revcomp(ds), revcomp(us)
    return signal_weight(us, ds)


def signal_weight(us, ds):
    if us == 'GT' and ds == 'AG':
        return 0
    else:
        return 1


def correct_reads(reads_cluster, ref_fasta, gtf_index, intron_index, ss_index, threads):
    # Load reference genome
    genome = Fasta(ref_fasta)

    corrected_reads = []
    jobs = []
    pool = Pool(threads, env.initializer, (None, genome.contig_len, genome, gtf_index, intron_index, ss_index, ))
    for cluster in grouper(reads_cluster, 250):
        jobs.append(pool.apply_async(correct_chunk, (cluster,)))
    pool.close()

    prog = ProgressBar()
    prog.update(0)
    cnt = 0
    circ_num = defaultdict(int)
    for job in jobs:
        tmp_cluster, tmp_num = job.get()
        corrected_reads += tmp_cluster
        for i in tmp_num:
            circ_num[i] += tmp_num[i]
        cnt += 1
        prog.update(100 * cnt // len(jobs))
    pool.join()
    prog.update(100)

    return circ_num, corrected_reads


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

    genome = Fasta(ref_fasta)

    circ_reads = defaultdict(list)
    circ_info = {}
    reads_df = []

    with open('{}/{}.fa'.format(out_dir, prefix), 'w') as fa:
        for reads, seqs, circ_id, strand, ss_id, us_free, ds_free, circ_len, isoforms in corrected_reads:
            for tmp_id, tmp_seq in zip(reads, seqs):
                fa.write('>{}\t{}\t{}\n{}\n'.format(tmp_id, circ_id, strand, tmp_seq))

            # circRNA information
            ctg, st, en = circ_pos(circ_id)
            field = circ_attr(gtf_idx, ctg, st, en, strand)

            tmp_attr = 'circ_id "{}"; splice_site "{}"; equivalent_seq "{}"; circ_type "{}"; circ_len "{}";'.format(
                circ_id,
                ss_id,
                equivalent_seq(genome, ctg, st, en, strand),
                field['circ_type'] if field else 'Unknown',
                circ_len,
            )
            if isoforms:
                tmp_attr += ' isoform "{}";'.format('|'.join(isoforms))

            for key in 'gene_id', 'gene_name', 'gene_type':
                if key in field:
                    tmp_attr += ' {} "{}";'.format(key, field[key])
            circ_info[circ_id] = [ctg, 'CIRI-long', 'circRNA', st, en, len(reads), strand, '.', tmp_attr, ]

            # Expression levels
            circ_reads[circ_id] += reads

            # Corrected reads
            for read_id in reads:
                read = cand_reads[read_id]
                tmp = [read_id, circ_id, read.circ_id, read.strand, read.cirexon,
                       read.ss, read.clip, read.segments, read.sample, read.type]
                reads_df.append(tmp)

    # circular reads
    reads_df = pd.DataFrame(
        reads_df, columns=['read_id', 'circ_id', 'tmp_id', 'strand', 'cirexons',
                           'signal', 'alignment', 'segments', 'sample', 'type'])
    reads_df.to_csv('{}/{}.reads'.format(out_dir, prefix), sep="\t", index=False)

    # circRNA information
    sorted_circ = sorted(list(circ_info), key=by_circ)
    with open('{}/{}.info'.format(out_dir, prefix), 'w') as out:
        for circ_id in sorted_circ:
            tmp_line = circ_info[circ_id]
            out.write('\t'.join([str(x) for x in tmp_line]) + '\n')

    # Expression level
    exp_df = {}
    for circ_id, reads in circ_reads.items():
        exp_df[circ_id] = Counter([cand_reads[i].sample for i in reads])
    exp_df = pd.DataFrame.from_dict(exp_df).transpose().fillna(0).reindex(sorted_circ)
    exp_df.to_csv('{}/{}.expression'.format(out_dir, prefix), sep="\t", index_label='circ_ID')

    return len(sorted_circ)


def equivalent_seq(genome, contig, start, end, strand):
    if strand is None:
        return 'Unknown'

    ds_seq = ''
    for i in range(100):
        if end + i > genome.contig_len[contig]:
            break
        if genome.seq(contig, start, start + i) == genome.seq(contig, end, end + i):
            ds_seq += genome.seq(contig, start, start + i)
        else:
            break

    us_seq = ''
    for j in range(100):
        if start - j < 0:
            break
        if genome.seq(contig, start - j, start) == genome.seq(contig, end - j, end):
            us_seq += genome.seq(contig, start - j, start)
        else:
            break

    tmp = us_seq[::-1] + ds_seq
    if strand == '+':
        return tmp
    else:
        return revcomp(tmp)


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
            if element.start <= start <= element.end and (element.strand == strand or strand is None):
                start_element[element.type].append(element)
            # end site
            if element.start <= end <= element.end and (element.strand == strand or strand is None):
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
            if strand is None or host_gene[gene_id].strand == strand:
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
