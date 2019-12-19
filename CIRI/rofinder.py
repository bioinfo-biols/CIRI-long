import os
import sys
import logging
from collections import namedtuple, defaultdict
LOGGER = logging.getLogger('CIRI-long')

import numpy as np
from skbio import DNA, local_pairwise_align_ssw


def find_consensus(header, seq, out_dir, debugging):
    from .preprocess import trim_primer
    from poa import consensus

    # Trim sequence
    if len(seq) <= 50:
        return None, None, None

    # Repeat segments
    chains = None
    is_circular = 0
    for k in [11, 9, 7]:
        tmp_chain, tmp_circular = ROF(header, seq, k=k)
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

    fasta = [('{}-{}'.format(s, e), seq[s:e]) for s, e in chains]
    ccs = consensus(fasta, alignment_type=1,
                    match=10, mismatch=-4, gap=-8, extension=-2, gap_affine=-24, extension_affine=-4,
                    debug=0)
    # if header == 'ENSMUST00000177117|ENSMUSG00000021546|13:58395270-58396899|-|211_544_aligned_43648_R_28_1692_15':
    #     print(header)

    segments = ';'.join(['{}-{}'.format(s, e) for s, e in chains])

    return segments, ccs, is_circular


def ROF(header, seq, k=11, p_match=0.8, p_indel=0.1, d_min=40, support_min=2):
    """
    Circular Finder

    Paramters
    ---------
    seq : str
        input sequence
    k : int, default 11
        k-mer length
    p_match : float, default 0.8
        probability of match in kmer
    p_indel : float, default 0.1
        probability of indel in kmer
    d_min : int, default 40
        minimum distance of repeat kmer
    support_min : int, default 19
        minimum support kmer of distance

    Returns
    -------
    junc : tuple
        position of junction site

    """

    assert 0 <= p_indel <= 1
    assert 0 <= p_match <= 1

    # Occurence of kmer
    kmer_occ = defaultdict(list)
    for i in range(len(seq) - k):
        kmer = seq[i:i + k]
        kmer_occ[kmer].append(i)

    # Distance of repeat kmers
    tuple_dis = defaultdict(int)
    for k_occ in kmer_occ.values():
        for x1, x2 in zip(k_occ[:-1], k_occ[1:]):
            # filter out repeat near than d_min
            if x2 - x1 >= d_min:
                tuple_dis[x2 - x1] += 1

    if len(tuple_dis) == 0:
        return None, None

    tuple_dis_mean = list(tuple_dis)[np.argmax(list(tuple_dis.values()))]

    # Minimum suppport kmer for circular segments
    if tuple_dis[tuple_dis_mean] <= support_min:
        return None, None

    # Random walk distribution for alignment length
    tuple_dis_delta = int(2.3 * np.sqrt(p_indel * tuple_dis_mean))

    # Chain of match kmers
    intervals = [(0, len(seq)), ]
    chains = [[], ]

    for i in range(len(seq) - k):
        j, score = best_hit(seq, seq[i:i+k], i + tuple_dis_mean - tuple_dis_delta, i + tuple_dis_mean + tuple_dis_delta)
        if score > int(0.25 * k):
            continue

        for x, (s, e) in enumerate(intervals):
            if s <= i < e and s <= j < e:
                intervals = intervals[:x] + [(s, j), (j, e)] + intervals[x+1:]
                chains[x].append(i)
                chains = chains[:x+1] + [[j],] + chains[x+1:]
                break
            elif s <= i <= e:
                chains[x].append(i)
            elif s <= j <= e:
                chains[x].append(j)
            else:
                pass

    if min([len(c) for c in chains]) == 0 or len(chains) <= 1:
        LOGGER.warn('Dropped strange repeat patterns for {}'.format(header))
        return None, None

    is_circular = 1

    # Over 75% similar for full length sequence
    if len(chains) > 2 and min([max(c) - min(c) for c in chains[:-1]]) < 0.75 * tuple_dis_mean:
        is_circular = 0

    # Drop repeat patterns that have too long head / tail region
    if chains[0][0] > 100 or chains[-1][-1] < len(seq) - 100:
        is_circular = 0

    # Extend chains to adjacent
    chains = [[min(c), max(c)] for c in chains]
    for i, (s, e) in enumerate(chains[:-1]):
        if e < chains[i+1][0]:
            chains[i] = (s, chains[i+1][0])
    chains[-1] = (chains[-1][0], chains[-1][1] + k)

    return chains, is_circular


def valid_junc(juncs, d_mean, d_delta):
    valid_s = []
    last_s = None
    for s in juncs:
        if last_s is not None and s != last_s:
            continue

        for x in range(s + d_mean - d_delta, s + d_mean + d_delta):
            if x in juncs:
                valid_s.append(s)
                last_s = x
                break

    if last_s is not None:
        valid_s.append(last_s)

    return valid_s


def best_hit(seq, seed, start, end):
    from Levenshtein import distance
    k = len(seed)
    pos = []
    dis = []
    for i in range(start, end):
        pos.append(i)
        dis.append(distance(seq[i:i + k], seed))

    min_x = np.argmin(dis)
    return pos[min_x], dis[min_x]


def worker(chunk, out_dir, debugging):
    ret = []
    for header, seq in chunk:
        segments, ccs, is_circular = find_consensus(header, seq, out_dir, debugging)
        if segments is None or ccs is None or is_circular == 0:
            continue
        ret.append((header, seq, segments, ccs))
    return len(chunk), ret


def find_ccs_reads(in_file, out_dir, prefix, threads, debugging):
    import gzip
    from multiprocessing import Pool
    from CIRI.utils import to_str
    from CIRI.logger import ProgressBar
    pool = Pool(threads)
    jobs = []

    # Auto detect input format
    is_fastq = 1
    is_gz = 1
    if in_file.endswith('.fa') or in_file.endswith('.fasta'):
        is_fastq = 0
        is_gz = 0
        fq = open(in_file, 'r')
    elif in_file.endswith('.fq') or in_file.endswith('.fastq'):
        is_gz = 0
        fq = open(in_file, 'r')
    elif in_file.endswith('.fq.gz') or in_file.endswith('.fastq.gz'):
        fq = gzip.open(in_file, 'rb')
    else:
        sys.exit('Wrong format of input')

    total_cnt = 0
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

        total_cnt += 1
        # Split into chunks
        chunk.append((header, seq))
        if len(chunk) == chunk_size:
            jobs.append(pool.apply_async(worker, (chunk, out_dir, debugging, )))
            # worker(chunk, out_dir, debugging)
            chunk = []
            chunk_cnt += 1

    if len(chunk) > 0:
        jobs.append(pool.apply_async(worker, (chunk, out_dir, debugging, )))
        chunk_cnt += 1
        # worker(chunk, out_dir, debugging)
    pool.close()

    prog = ProgressBar()
    prog.update(0)
    finished_chunk = 0

    total_reads = 0
    ro_reads = 0

    ccs_seq = {}
    with open('{}/{}.ccs.fa'.format(out_dir, prefix), 'w') as out, \
            open('{}/{}.raw.fa'.format(out_dir, prefix), 'w') as trimmed:
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
    with open('{}/{}.ccs.fa'.format(out_dir, prefix), 'r') as f:
        for line in f:
            header = line.rstrip()
            content = header.split('\t')
            seq = f.readline().rstrip()
            ccs_seq[content[0].lstrip('>')] = [content[1], seq]

    with open('{}/{}.raw.fa'.format(out_dir, prefix), 'r') as f:
        for line in f:
            header = line.rstrip().split('\t')[0].lstrip('>')
            seq = f.readline().rstrip()
            ccs_seq[header].append(seq)
    return ccs_seq
