import os
import sys
import time
import numpy as np
from functools import cmp_to_key
from collections import namedtuple

import pysam
import poagraph
import seqgraphalignment
import mappy as mp
from Levenshtein import distance
from skbio import DNA
from skbio.alignment import local_pairwise_align_ssw

Node = namedtuple('path', 'row col score')
SCORE = {
    'A': {'A': 10, 'G': -5, 'C': -7, 'T': -7},
    'G': {'A': -5, 'G': 10, 'C': -7, 'T': -7},
    'C': {'A': -7, 'G': -7, 'C': 10, 'T': -5},
    'T': {'A': -7, 'G': -7, 'C': -5, 'T': 10},
}


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


class ROFinder(object):
    def __init__(self, sequence, match=10, mismatch=-3, gap=-5, local_alignment=True):
        self.sequence = sequence
        self._match = 10
        self._mismatch = -3
        self._gap = -5
        self._len = len(sequence)
        self.localAlignment = local_alignment
        # Init matrix
        self.init_mtx()

    def init_mtx(self):
        """Init scoring matrix"""
        self.mtx = np.zeros((self._len + 1, self._len + 1))
        for i in range(self._len + 1):
            self.mtx[0][i] = 0 if self.localAlignment is True else self._gap * i

        # For a symmetric matrix, only the upper half is computed
        for i in range(1, self._len + 1):
            for j in range(i, self._len + 1):
                self.mtx[i][j] = max(
                    self.mtx[i - 1][j - 1] + SCORE[self.sequence[i - 1]][self.sequence[j - 1]],
                    self.mtx[i - 1][j] + self._gap,
                    self.mtx[i][j - 1] + self._gap,
                )
                if self.localAlignment is True:
                    self.mtx[i][j] = max(self.mtx[i][j], 0)

    def symmetrize(self):
        """Symmetric the full matrix

        Return: symmetrized scoring matrix
        """
        return self.mtx + self.mtx.T - np.diag(self.mtx.diagonal())

    @staticmethod
    def cmp_parent(a, b):
        if a[0] == b[0]:
            # Proprity From small to large
            return b[3] - a[3]
        else:
            # Score from large to small
            return a[0] - b[0]

    def recursive_path(self, i, j):
        """Recursive find path
        Input:
            i: index of row, j: index of col
        Return:
            path of best alignment
        """
        if i is None and j is None:
            return []
        elif i == 0:
            return [Node(i, j, self.mtx[i][j]), ]
        else:
            pass

        parents = [
            (self.mtx[i - 1][j - 1] + SCORE[self.sequence[i - 1]][self.sequence[j - 1]], i - 1, j - 1, 1),
            (self.mtx[i - 1][j] + self._gap, i - 1, j, 2),
            (self.mtx[i][j - 1] + self._gap, i, j - 1, 3),
        ]
        if self.localAlignment:
            parents.append((0, None, None, 4))

        l_score, l_row, l_col, l_priority = sorted(parents, key=cmp_to_key(self.cmp_parent), reverse=True)[0]
        return self.recursive_path(l_row, l_col) + [Node(i, j, self.mtx[i][j]), ]

    def find_ro_paths(self):
        """Find Reverse Overlap Regions"""
        ro_paths = {}
        idx = self._len
        while idx > 0:
            tmp_path = self.recursive_path(idx, self._len)
            idx -= 1
            if len(tmp_path) < 15:
                continue

            s_node, e_node = tmp_path[0], tmp_path[-1]
            # Path start at head line
            if s_node.row != 0:
                continue
            if s_node.col not in ro_paths or ro_paths[s_node.col].score < e_node.score:
                ro_paths[s_node.col] = e_node

        return ro_paths

    @staticmethod
    def merge_ro_paths(ro_paths):
        """Merge RO paths

        Input:
            ro_paths: list

        Return:
            merged_paths: list
        """
        tmp_site = 0
        ro_sects = [[], ]
        for i in sorted(list(ro_paths)):
            if i > tmp_site + 20:
                ro_sects.append([i, ])
            else:
                ro_sects[-1].append(i)
            tmp_site = i

        merged_paths = {}
        for cols in ro_sects:
            col = sorted(cols, key=lambda x: ro_paths[x].score)[-1]
            merged_paths[col] = ro_paths[col]

        return merged_paths

    def filter_ro_paths(self, ro_paths):
        if len(ro_paths) <= 2:
            return ro_paths

        cols = sorted(list(ro_paths))
        kmers = [self.sequence[i:i + 11] for i in cols]

        valid_paths = []
        for col, kmer in zip(cols, kmers):
            if col == 0 or distance(kmer, kmers[0]) <= 3:
                valid_paths.append((col, ro_paths[col]))
        return valid_paths

    def circular_segments(self, ro_paths):
        segments = []
        divs = sorted(list(ro_paths))
        for i, x in enumerate(divs):
            if i < len(divs) - 1:
                segments.append((x[0], divs[i + 1][0]))
            else:
                segments.append((x[0], self._len - 1))
        return segments

    def get_seq(self, segment):
        return self.sequence[segment[0]:segment[1]]

    def consensus(self):
        """Generate consensus sequence using RO features"""
        ro_paths = self.find_ro_paths()
        merged_paths = self.merge_ro_paths(ro_paths)
        filtered_paths = self.filter_ro_paths(merged_paths)
        ro_segments = self.circular_segments(filtered_paths)
        ro_seqs = [('{}-{}'.format(*segment), self.get_seq(segment)) for segment in ro_segments]
        if len(ro_segments) < 2:
            return None
        else:
            # all_consenses = MSA(ro_seqs)
            # ccs = sorted([''.join(x[1]) for x in all_consenses], key=lambda x: len(x), reverse=True)[0]
            return ro_seqs


def trim_sequence(seq, primer = 'AAGCAGTGGTATCAACGCAGAGTAC'):
    l_alignments, l_score, l_positions = local_pairwise_align_ssw(DNA(seq[:100]), DNA(primer))
    l_trim = l_positions[0][1] + 1 if l_score >= 25 else 0

    r_alignments, r_score, r_positions = local_pairwise_align_ssw(DNA(seq[-100:]), DNA(mp.revcomp(primer)))
    r_trim = -100 + r_positions[0][0] if r_score >= 25 else -1
    return seq[l_trim:r_trim]


def MSA(fasta):
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

    return graph.allConsenses()


def generate_consensus(data_dir, sample, read_id, circ_id):
    sam = pysam.AlignmentFile('{}/{}_circRNA_pseudo.sorted.bam'.format(data_dir, sample))
    for x in sam.fetch(circ_id, multiple_iterators=True):
        if x.query_name != read_id:
            continue
        if x.is_supplementary or x.is_secondary:
            continue
        read = x
        break

    read_seq = read.query_sequence
    trimmed_seq = trim_sequence(read_seq)

    # ro_finder = ROFinder(trimmed_seq, match=10, mismatch=-3, gap=-5)
    # ccs = ro_finder.consensus()

    ccs = blast_search(trimmed_seq)
    return ccs


def blast_search(seq):
    def merge_hsp(x):
        tmp = x[0]
        merged = [tmp, ]
        for i in x[1:]:
            if i - tmp <= 35:
                tmp = i
            else:
                merged.append(i)
                tmp = i
        return merged

    def is_seed(pos):
        if len(pos) <= 2:
            return 1
        else:
            segments = ([val - pos[i - 1] for i, val in enumerate(pos) if i > 0])
            return 1 if np.max(segments) - np.min(segments) <= 20 else 0

    def extend_hsp(x, seq_len):
        extended = []
        for i, (start, end) in enumerate(x):
            start = min(0, start) if i == 0 else x[i - 1][1]
            extended.append([start, end])

        extended[-1][1] = seq_len
        return extended

    kmers = [seq[i:i + 11] for i in range(len(seq)) if i + 11 < len(seq)]
    hsp = {}
    for i, j in enumerate(kmers):
        hsp.setdefault(j, []).append(i)

    for i in hsp:
        hsp[i] = merge_hsp(hsp[i])

    cmp_hsp = lambda a, b: min(hsp[b]) - min(hsp[a]) if len(hsp[a]) == len(hsp[b]) else len(hsp[a]) - len(hsp[b])

    seed_num = max([len(j) for _, j in hsp.items()])
    seed_list = sorted([i for i, j in hsp.items() if len(j) == seed_num and is_seed(j)], key=cmp_to_key(cmp_hsp),
                       reverse=True)

    try:
        init_seed = [(start, end + 11) for start, end in zip(hsp[seed_list[0]], hsp[seed_list[-1]])]
        extended_seed = extend_hsp(init_seed, len(seq))

        if len(extended_seed) >= 2:
            seqs = [('{}-{}'.format(start, end), seq[start:end]) for (start, end) in extended_seed]
            # ccs = MSA(seqs)
            ccs = seqs
        else:
            ccs = None
    except IndexError as e:
        ccs = None

    return ccs


def worker1(seq):
    trimmed_seq = trim_sequence(seq)
    ro_finer = ROFinder(trimmed_seq)
    return ro_finer.consensus()


def worker2(seq):
    trimmed_seq = trim_sequence(seq)
    if len(trimmed_seq) == 0:
        return None
    return blast_search(trimmed_seq)


def main():
    import subprocess
    from multiprocessing import Pool
    data_dir = '/home/zhangjy/03.CIRIpacbio/1908_Nanopore/alignment/'
    fq_dir = '/home/zhangjy/03.CIRIpacbio/1908_Nanopore/fastq'
    for sample in 'barcode20', 'barcode21', 'barcode22':
        print('Loading {}'.format(sample))
        # wc_results = subprocess.getoutput('wc -l {}/RO_reads/{}_aligned_reads.list'.format(data_dir, sample))
        # total_cnt = int(wc_results.split(' ')[0])
        #
        # # Load reads id
        # cand_reads = {}
        # with open('{}/RO_reads/{}_aligned_reads.list'.format(data_dir, sample), 'r') as f:
        #     for line in f:
        #         content = line.rstrip().split('\t')
        #         q_name, q_len, circ_id, circ_length, is_circ, aligned_length, is_ro, seed_len, seed = content
        #         cand_reads.setdefault(circ_id, []).append(q_name)
        #
        # # Load read sequence
        # cand_seqs = {}
        # sam = pysam.AlignmentFile('{}/{}_circRNA_pseudo.sorted.bam'.format(data_dir, sample))
        # for circ_id in cand_reads:
        #     for x in sam.fetch(circ_id, multiple_iterators=True):
        #         if x.query_name not in cand_reads[circ_id]:
        #             continue
        #         if x.is_supplementary or x.is_secondary:
        #             continue
        #         cand_seqs[x.query_name] = x.query_sequence
        import gzip
        pool = Pool(4)
        jobs = []

        total_cnt = 0
        fastq = gzip.open('{}/{}.fq.gz'.format(fq_dir, sample))
        for line in fastq:
            total_cnt += 1
            header = line.rstrip()
            seq = fastq.readline().rstrip()
            separator = fastq.readline()
            qual = fastq.readline()

            jobs.append(pool.apply_async(worker1, (seq, )))
        pool.close()

        read_cnt = 0
        ro_cnt = 0
        prog = ProgressBar()
        prog.update(0)
        for job in jobs:
            ccs = job.get()
            if ccs is not None:
                ro_cnt += 1
            read_cnt += 1
            prog.update(100 * read_cnt // total_cnt)
        pool.join()
        prog.update(100)
        print('Sample: {}, Total reads: {}, RO reads: {}\n'.format(sample, total_cnt, ro_cnt))


if __name__ == '__main__':
    main()
