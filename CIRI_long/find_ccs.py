import os
import sys
import logging
from collections import namedtuple, defaultdict
LOGGER = logging.getLogger('CIRI-long')

import numpy as np
from pyccs import find_consensus


def worker(chunk):
    ret = []
    for header, seq in chunk:
        segments, ccs = find_consensus(seq)
        if segments is None or ccs is None:
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
            content = header.split()
            seq = f.readline().rstrip()
            ccs_seq[content[0].lstrip('>')] = [content[1], seq]

    with open('{}/tmp/{}.raw.fa'.format(out_dir, prefix), 'r') as f:
        for line in f:
            header = line.rstrip().split()[0].lstrip('>')
            seq = f.readline().rstrip()
            ccs_seq[header].append(seq)
    return ccs_seq
