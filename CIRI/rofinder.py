import numpy as np
from skbio import DNA, local_pairwise_align_ssw
from collections import namedtuple, defaultdict

import poa


def find_consensus(header, seq, out_dir, debugging):
    from .preprocess import trim_primer
    from .poa import consensus

    # Trim sequence
    trimmed_seq = trim_primer(seq)
    if len(trimmed_seq) <= 50:
        return None, None

    junc_sites = ROF(trimmed_seq)
    if junc_sites is None:
        return None, None

    fasta = [('{}-{}'.format(i, j), trimmed_seq[i:j]) for i, j in zip(junc_sites[:-1], junc_sites[1:])]

    print(header)
    ccs = consensus(fasta)

    # poa_graph = partial_order_alignment(fasta)
    # ccs_reads = [''.join(i[1]) for i in poa_graph.allConsenses()]
    # ccs = sorted(ccs_reads, key=lambda x: len(x), reverse=True)[0]

    # if debugging is True:
    #     alignments = poa_graph.generateAlignmentStrings()
    #     with open('{}/tmp/{}.msa'.format(out_dir, header), 'w') as out:
    #         for label, alignstring in alignments:
    #             out.write("{0:15s} {1:s}\n".format(label, alignstring))
    #     with open('{}/tmp/{}.fa'.format(out_dir, header), 'w') as out:
    #         for label, sequence in fasta:
    #             out.write('>{}\n{}\n'.format(label, sequence))

    segments = ';'.join(['{}-{}'.format(i, j) for i, j in zip(junc_sites[:-1], junc_sites[1:])])
    return segments, ccs


def ROF(seq, k=11, p_match=0.8, p_indel=0.1, d_min=40, support_min=10):
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
        return None

    tuple_dis_mean = list(tuple_dis)[np.argmax(list(tuple_dis.values()))]

    # Minimum suppport kmer for circular segments
    if tuple_dis[tuple_dis_mean] <= support_min:
        return None

    # Random walk distribution for alignment length
    tuple_dis_delta = int(2.3 * np.sqrt(p_indel * tuple_dis_mean))

    # Kmer with maximum occurence
    sorted_kmers = sorted(list(kmer_occ), key=lambda x: (len(kmer_occ[x]), -kmer_occ[x][0]), reverse=True)
    occ_max = len(kmer_occ[sorted_kmers[0]])

    cand_junc = []
    for kmer in sorted_kmers:
        if len(kmer_occ[kmer]) != occ_max:
            break
        valid_s = valid_junc(kmer_occ[kmer], tuple_dis_mean, tuple_dis_delta)
        if len(valid_s) == 0:
            continue
        cand_junc.append(valid_s)


    if len(cand_junc) == 0:
        return None

    final_junc = sorted(cand_junc, key=lambda x: (len(x), -x[0]), reverse=True)[0]

    while final_junc[0] >= 15:
        if tuple_dis_mean < final_junc[0]:
            x, score = best_hit(seq, seq[final_junc[0]:final_junc[0] + k],
                                max(0, final_junc[0] - tuple_dis_mean - tuple_dis_delta),
                                final_junc[0] - tuple_dis_mean + tuple_dis_delta)
            if score > 3:
                return None
            final_junc = [x, ] + final_junc
        else:
            # final_junc = [0, ] + final_junc
            break

    while final_junc[-1] <= len(seq) - 15:
        if final_junc[-1] + tuple_dis_mean < len(seq):
            x, score = best_hit(seq, seq[final_junc[-1]:final_junc[-1] + k],
                                final_junc[-1] + tuple_dis_mean - tuple_dis_delta,
                                min(final_junc[-1] + tuple_dis_mean + tuple_dis_delta, len(seq)))
            if score > 3:
                return None
            final_junc.append(x)
        else:
            final_junc.append(len(seq))

    return final_junc


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
        segments, ccs = find_consensus(header, seq, out_dir, debugging)
        ret.append((header, segments, ccs))
    return ret


def find_ccs_reads(in_file, out_dir, prefix, threads, debugging):
    import gzip
    from multiprocessing import Pool
    from .utils import to_str
    from .logger import ProgressBar
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

    with open('{}/{}.ccs.fa'.format(out_dir, prefix), 'w') as out:
        for job in jobs:
            ret = job.get()
            for header, segments, ccs in ret:
                total_reads += 1
                if segments is None and ccs is None:
                    continue
                ro_reads += 1
                out.write('>{}\t{}\t{}\n{}\n'.format(header, segments, len(ccs), ccs))

            finished_chunk += 1
            prog.update(100 * finished_chunk // chunk_cnt)
    prog.update(100)
    pool.join()

    return total_reads, ro_reads
