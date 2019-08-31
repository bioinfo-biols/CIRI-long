import os
import sys
import time
import re
from collections import defaultdict
from multiprocessing import Pool

import pysam
from skbio import DNA
from skbio.alignment import local_pairwise_align_ssw

CIRC = {}


class ProgressBar(object):

    def __init__(self, width=50):
        self.last_x = -1
        self.width = width

    def update(self, x):
        assert 0 <= x <= 100
        if self.last_x == int(x):
            return
        self.last_x = int(x)
        p = int(self.width * (x / 100.0))
        time_stamp = time.strftime("[%a %Y-%m-%d %H:%M:%S]", time.localtime())
        sys.stderr.write('\r%s [%-5s] [%s]' % (time_stamp, str(int(x)) + '%', '#' * p + '.' * (self.width - p)))
        sys.stderr.flush()
        if x == 100:
            sys.stderr.write('\n')


def load_circ(circ_file):
    circ_info = defaultdict(list)
    with open(circ_file, 'r') as f:
        for line in f:
            contig, start, end = re.split(r'[:|]', line.rstrip())
            circ_info[contig].append((int(start), int(end)))
    return circ_info


def circular_reads(bam_file, contig):
    if contig not in CIRC:
        return None

    sam = pysam.AlignmentFile(bam_file, 'rb')
    circ_reads = 0
    ro_reads = 0
    content = []
    for (start, end) in CIRC[contig]:
        circ_id = '{}:{}-{}'.format(contig, start, end)
        circ_length = end - start + 1
        for read in sam.fetch(circ_id, multiple_iterators=True):
            # Remove supplementary and secondary alignment
            if read.is_supplementary or read.is_secondary:
                continue

            # Spanning Back-Spliced Junction
            is_bsj = 0
            for (block_start, block_end) in read.get_blocks():
                if block_start + 50 <= circ_length <= block_end - 50:
                    is_bsj = 1

            if is_bsj == 0:
                continue

            # Aligned more than 80%
            block_start = min(i[0] for i in read.get_blocks())
            block_end = max(i[1] for i in read.get_blocks())
            block_length = block_end - block_start
            aligned_length = sum([block_end - block_start for (block_start, block_end) in read.get_blocks()])
            read_length = read.infer_read_length()

            is_circ = 0
            if circ_length <= block_length or aligned_length >= 0.6 * read_length:
                circ_reads += 1
                is_circ = 1

            # Seed matching
            seed, is_ro = get_subsequence(read.query_sequence[15:])
            if seed is not None and is_ro is True:
                ro_reads += 1

            if is_ro or is_circ:
                tmp = (read.query_name, read_length, circ_id, circ_length,
                       is_circ, aligned_length, is_ro, len(seed) if seed else "NA", seed)
                content.append(tmp)
    return circ_reads, ro_reads, content


def get_subsequence(sequence, is_circ=False):
    if len(sequence) <= 200:
        return sequence, is_circ

    # Seed matching
    junc = 100
    while True:
        alignments, score, positions = local_pairwise_align_ssw(DNA(sequence[:junc]), DNA(sequence[junc:]))
        (l_st, l_en), (r_st, r_en) = positions
        if r_st == 0:
            break
        else:
            junc += r_st

    seed = sequence[:junc]
    # Filter unmapped reads
    if l_st > 10:
        return None, False

    # keep reads containing 5' RO
    if r_en >= len(sequence) - junc - 10:
        return seed, True

    # Full length RO
    if l_en >= junc - 10:
        if len(sequence) <= 2 * junc + 20:
            return seed, True

        is_circ = True
        sub_seed, sub_circ = get_subsequence(sequence[junc:], is_circ)
        if sub_circ is True and abs(len(sub_seed) - len(seed)) <= 20:
            return seed, True
    return None, False


def initializer(circ_info):
    global CIRC
    CIRC = circ_info


def count_bam(bam_file, circ_info, out_dir, sample):
    thread = 6

    initializer(circ_info)

    pool = Pool(thread, initializer, (circ_info, ), maxtasksperchild=4)

    jobs = []
    for contig in circ_info:
        jobs.append(pool.apply_async(circular_reads, (bam_file, contig, )))
    pool.close()

    prog = ProgressBar()
    prog.update(0)
    cnt = 0
    circ_reads, ro_reads = 0, 0

    with open('{}/{}_aligned_reads.list'.format(out_dir, sample), 'w') as out:
        for job in jobs:
            tmp_circ, tmp_ro, tmp_content = job.get()
            circ_reads += tmp_circ
            ro_reads += tmp_ro

            for line in tmp_content:
                out.write('\t'.join([str(x) for x in line]) + '\n')

            cnt += 1
            prog.update(100 * cnt // len(circ_info))
    pool.join()
    prog.update(100)

    return circ_reads, ro_reads


def main():
    os.chdir('/home/zhangjy/03.CIRIpacbio/1908_Nanopore/alignment')

    out_dir = '/home/zhangjy/03.CIRIpacbio/1908_Nanopore/alignment/RO_reads'
    circ_file = './all_mouse_liver.uniq.list'
    circ_info = load_circ(circ_file)

    for sample in ['barcode20', 'barcode21', 'barcode22']:
        bam_file = './{}_circRNA_pseudo.sorted.bam'.format(sample)
        print('Loading Minimap2 results for: ' + sample)
        circ_reads, ro_reads = count_bam(bam_file, circ_info, out_dir, sample)
        print('Circular reads for {}: {}'.format(sample, circ_reads))
        print('Circular reads for {}: {}'.format(sample, ro_reads))


if __name__ == '__main__':
    main()
