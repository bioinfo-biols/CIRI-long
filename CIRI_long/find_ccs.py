import os
import sys
import logging
from collections import namedtuple, defaultdict
LOGGER = logging.getLogger('CIRI-long')

import numpy as np


def collect_kmers(seq, k=8, use_hpc=False, is_circular=False):
    """
    split sequence into kmers
    :param seq: sequence
    :param k: k-mer size
    :param use_hpc: True if use homopolymer-compressed kmers
    :return: dict of kmers and occurence position
    """
    kmers = {}
    kmer_occ = defaultdict(list)

    split_func = split_hpc_kmers if use_hpc else split_kmers
    if is_circular:
        tmp_kmers, tmp_kmer_occ = collect_kmers(seq * 2, k, use_hpc, False)
        for i in tmp_kmers:
            if i < len(seq):
                continue
            kmers[i - len(seq)] = tmp_kmers[i]
            kmer_occ[tmp_kmers[i]].append(i - len(seq))
    else:
        for i, kmer in split_func(seq, k):
            kmers[i] = kmer
            kmer_occ[kmer].append(i)

    return kmers, kmer_occ


def all_kmers(k):
    import itertools
    all_bases = ['A', 'C', 'G', 'T']
    all_list = [all_bases, ] * k
    all_kmers = {''.join(kmer): idx for idx, kmer in enumerate(itertools.product(*all_list))}
    return all_kmers


def split_kmers(seq, k=8):
    for x in range(len(seq)):
        if x >= k - 1:
            yield x, seq[x-k+1:x+1]


def split_hpc_kmers(seq, k=8):
    hpc = [seq[0], ]
    for x, (i, j) in enumerate(zip(seq[:-1], seq[1:])):
        if i != j:
            hpc.append(j)
            if len(hpc) >= k:
                yield x + 1, ''.join(hpc[-k:])


def compress_seq(seq):
    hpc = [seq[0], ]
    for x, (i, j) in enumerate(zip(seq[:-1], seq[1:])):
        if i != j:
            hpc.append(j)
    return ''.join(hpc)


def estimate_distance(kmer_occ, p_indel, d_min):
    from scipy.stats.kde import gaussian_kde
    tuple_dis = []
    for k_occ in kmer_occ.values():
        for x1, x2 in zip(k_occ[:-1], k_occ[1:]):
            if x2 - x1 < d_min:
                continue
            tuple_dis.append(x2 - x1)
    if len(tuple_dis) <= 2:
        return None, None, 0

    dis_uniq = np.unique(tuple_dis)
    if len(dis_uniq) == 1:
        tuple_dis_mean = dis_uniq[0]
    else:
        # Kernel density estimation
        kernel = gaussian_kde(tuple_dis)
        tuple_dis_mean = dis_uniq[np.argmax(kernel(dis_uniq))]

    # Random walk model
    tuple_dis_delta = int(2.3 * np.sqrt(p_indel * tuple_dis_mean))
    tuple_dis_support = len([i for i in tuple_dis if abs(i - tuple_dis_mean) <= tuple_dis_delta])

    return tuple_dis_mean, tuple_dis_delta, tuple_dis_support


def optimal_hit(kmers, seed, s, e, use_hpc):
    """
    Return hit with smallest distance and nearest to optimal distance
    """
    from Levenshtein import distance
    from operator import itemgetter

    hits = []
    for i in range(s, e + 1):
        if i not in kmers:
            continue
        hits.append([i, kmers[i], distance(seed, kmers[i]), max(i - s, i - e)])

    # No hit
    if len(hits) == 0:
        return None, None, len(seed), None

    # For HPC hits
    if use_hpc and len(hits) > 1:
        hits = collapse_hits(hits)

    return sorted(hits, key=itemgetter(2, 3))[0]


def collapse_hits(hits):
    """
    For hits with same HPC k-mer content and distance, collapse into the last base
    """
    collapsed = []
    for i, j in zip(hits[:-1], hits[1:]):
        if j[0] == i[0] + 1 and j[1] == i[1]:
            continue
        else:
            collapsed.append(j)
    return collapsed


def circular_hits(kmers, k, tuple_dis_mean, tuple_dis_delta, p_match, use_hpc):
    pos = sorted(list(kmers))
    hits = {}
    for i in pos:
        j, _, score, _ = optimal_hit(kmers, kmers[i],
                                     i + tuple_dis_mean - tuple_dis_delta,
                                     i + tuple_dis_mean + tuple_dis_delta,
                                     use_hpc)
        if score > int(k * (1 - p_match)):
            continue
        hits[i] = (i, j, j - i, score)
    return hits


def optimal_chains(hits):
    """
    Intersect hits into chains
    """
    pos = sorted(list(hits))
    s0, e0, _, score0 = hits[pos[0]]

    chains = [[(s0, e0, score0), ], ]
    for i in pos[1:]:
        s, e, _, score = hits[i]
        if chains[-1][-1][0] < s < chains[-1][-1][1] <= e:
            chains[-1].append((s, e, score))
        elif chains[-1][-1][1] < s:
            chains.append([(s, e, score), ])
        else:
            pass

    # Sort chains by length
    chains = sorted(chains, key=lambda x: x[-1][1] - x[0][0], reverse=True)

    return chains[0]


def split_sequence(primary_chain, seq, p_match, p_indel):
    from Bio import pairwise2

    match = 10
    mismatch = -4
    gap = -8
    extension = -2

    final = primary_chain[-1][1]

    for s, e, score in primary_chain:
        if score == 0:
            break
    boundaries = [s, e, ]

    idx = 0

    hits = {hit[0]: hit for hit in primary_chain}
    hit_pos = [hit[0] for hit in primary_chain]

    last_x = e
    while last_x < final:
        if last_x in hits:
            idx = hit_pos.index(last_x)
            last_x = hits[last_x][1]
            boundaries.append(last_x)
        else:
            # Find nearest hits
            while hit_pos[idx] < last_x:
                if idx >= len(hit_pos) - 1 or hit_pos[idx + 1] >= last_x:
                    break
                idx += 1
            if idx >= len(hit_pos) - 1:
                break

            # local alignment
            q_s, r_s, _ = hits[hit_pos[idx]]
            q_e, r_e, _ = hits[hit_pos[idx + 1]]

            alignments = pairwise2.align.globalms(seq[q_s:q_e + 1], seq[r_s:r_e + 1],
                                                  match, mismatch, gap, extension)

            q_shift, r_shift = 0, 0
            q_seq, r_seq, score, _, _ = alignments[0]
            for i, j in zip(q_seq, r_seq):
                q_shift += int(i in 'atcgATCG')
                r_shift += int(j in 'atcgATCG')
                if q_shift == last_x - q_s:
                    break
            # TODO: Stop searching when the score is too slow
            last_x = r_s + r_shift
            boundaries.append(last_x)

    segments = [(i, j) for i, j in zip(boundaries[:-1], boundaries[1:])]
    if segments[-1][-1] < final:
        segments.append((segments[-1][-1], final))

    return segments


def circular_finder(read_id, seq, k=8, use_hpc=True, p_match=.85, p_indel=.1, d_min=40, support_min=2):
    assert 0 <= p_indel <= 1
    assert 0 <= p_match <= 1

    # Split k-mers
    kmers, kmer_occ = collect_kmers(seq, k, use_hpc)

    # Approximately distance of circular segments
    tuple_dis_mean, tuple_dis_delta, tuple_dis_support = estimate_distance(kmer_occ, p_indel, d_min)
    if tuple_dis_support < support_min:
        return None, None

    # Hits of kmers
    tuple_hits = circular_hits(kmers, k, tuple_dis_mean, tuple_dis_delta, p_match, use_hpc)
    if not tuple_hits:
        return None, None

    # Primary chain
    primary_chain = optimal_chains(tuple_hits)

    # Split into segments
    segments = split_sequence(primary_chain, seq, p_match, p_indel)

    is_circular = 1
    # Head and tail position
    if segments[0][0] > 100 or segments[-1][-1] < len(seq) - 100:
        is_circular = 0

    if len(segments) < 2:
        is_circular = 0

    if len(segments) == 2 and segments[-1][1] - segments[-1][0] < 30:
        is_circular = 0

    return segments, is_circular


def find_consensus(header, seq):
    from spoa import poa
    from Levenshtein import distance

    # if header not in 'ENSMUST00000021181|ENSMUSG00000020831|11:70236877-70237625|-|197_1126_aligned_43558_F_52_821_60':
    #     return None, None, None

    # Trim sequence
    if len(seq) <= 50:
        return None, None, None

    # Repeat segments
    chains = None
    is_circular = 0
    for k, is_hpc in [(11, False), (8, False), (11, True), (8, True)]:
        tmp_chain, tmp_circular = circular_finder(header, seq, k=k, use_hpc=is_hpc)
        if tmp_chain is None:
            continue
        if tmp_circular == 1:
            chains = tmp_chain
            is_circular = 1
            break
        elif chains is None:
            chains = tmp_chain
        else:
            pass
    if chains is None:
        return None, None, None

    # Chains
    if len(chains) < 2:
        return None, None, None

    fasta = [seq[s:e] for s, e in chains]
    ccs, _ = poa(fasta, 2, False, 10, -4, -8, -2, -24, -1)

    # Check segment similarity
    tail = fasta[-1]
    if len(fasta) == 2:
        dis_body = distance(fasta[0][:len(tail)], ccs[:len(tail)]) / len(tail)
    else:
        dis_body = max([distance(i, ccs) / len(ccs) for i in fasta[:-1]])
    dis_tail = distance(tail, ccs[:len(tail)]) / len(tail)

    if dis_body > 0.2 or dis_tail > 0.35:
        return None, None, None

    segments = ';'.join(['{}-{}'.format(s, e) for s, e in chains])

    return segments, ccs, is_circular


def worker(chunk):
    ret = []
    for header, seq in chunk:
        segments, ccs, is_circular = find_consensus(header, seq)
        if segments is None or ccs is None or is_circular == 0:
            continue
        ret.append((header, seq, segments, ccs))
    return len(chunk), ret


def find_ccs_reads(in_file, out_dir, prefix, threads, debugging):
    import gzip
    from multiprocessing import Pool
    from CIRI_long.utils import to_str
    from CIRI_long.logger import ProgressBar
    pool = Pool(threads)
    jobs = []

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

    chunk = []
    chunk_size = 250
    chunk_cnt = 0
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
            jobs.append(pool.apply_async(worker, (chunk, )))
            chunk = []
            chunk_cnt += 1

    if len(chunk) > 0:
        jobs.append(pool.apply_async(worker, (chunk, )))
        chunk_cnt += 1
    pool.close()
    fq.close()

    prog = ProgressBar()
    prog.update(0)
    finished_chunk = 0

    total_reads = 0
    ro_reads = 0

    ccs_seq = {}
    with open('{}/tmp/{}.ccs.fa'.format(out_dir, prefix), 'w') as out, \
            open('{}/tmp/{}.raw.fa'.format(out_dir, prefix), 'w') as trimmed:
        for job in jobs:
            tmp_cnt, ret = job.get()
            total_reads += tmp_cnt
            for header, seq, segments, ccs in ret:
                ro_reads += 1
                out.write('>{}\t{}\t{}\n{}\n'.format(header, segments, len(ccs), to_str(ccs)))
                trimmed.write('>{}\n{}\n'.format(header, seq))
                ccs_seq[header] = [segments, to_str(ccs), seq]

            finished_chunk += 1
            prog.update(100 * finished_chunk // chunk_cnt)
    prog.update(100)
    pool.join()

    return total_reads, ro_reads, ccs_seq


def load_ccs_reads(out_dir, prefix):
    ccs_seq = {}
    with open('{}/tmp/{}.ccs.fa'.format(out_dir, prefix), 'r') as f:
        for line in f:
            header = line.rstrip()
            content = header.split('\t')
            seq = f.readline().rstrip()
            ccs_seq[content[0].lstrip('>')] = [content[1], seq]

    with open('{}/tmp/{}.raw.fa'.format(out_dir, prefix), 'r') as f:
        for line in f:
            header = line.rstrip().split('\t')[0].lstrip('>')
            seq = f.readline().rstrip()
            ccs_seq[header].append(seq)
    return ccs_seq
