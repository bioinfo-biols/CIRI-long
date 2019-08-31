import numpy as np
from skbio import DNA
from skbio import local_pairwise_align_ssw
from collections import namedtuple
Node = namedtuple('path', 'row col score')


class TracerError(Exception):
    pass


class ROFinder(object):
    def __init__(self, sequence, match=1, mismatch=-1, gap=-2):
        self.sequence = sequence
        self._match = match
        self._mismatch = mismatch
        self._gap = gap
        self._len = len(sequence)
        self.mtx = None
        self.tracer = None
        # Init matrix
        self.init_mtx()

    def init_mtx(self):
        # Init scoring matrix
        self.mtx = np.zeros((self._len + 1, self._len + 1)).astype(int)
        for i in range(self._len + 1):
            self.mtx[0][i] = 0

        # Init tracer matrix
        self.tracer = np.zeros((self._len + 1, self._len + 1)).astype(int)

        # For a symmetric matrix, only the upper half is computed
        for i in range(1, self._len + 1):
            for j in range(i, self._len + 1):
                p = self.mtx[i - 1][j - 1] + self.kmer_score(self.sequence[i - 1], self.sequence[j - 1])
                q = self.mtx[i - 1][j] + self._gap
                r = self.mtx[i][j - 1] + self._gap
                s = 0

                # back trace
                self.mtx[i][j] = max(p, q, r, s)
                if self.mtx[i][j] == p:
                    self.tracer[i][j] = 1
                elif self.mtx[i][j] == q:
                    self.tracer[i][j] = 2
                elif self.mtx[i][j] == r:
                    self.tracer[i][j] = 3
                else:
                    self.tracer[i][j] = 4

    def kmer_score(self, a, b):
        from Levenshtein import distance
        if distance(a, b) < len(a) // 3:
            return self._match
        else:
            return self._mismatch

    def symmetrize(self):
        """Symmetric the full matrix

        Return: symmetrized scoring matrix
        """
        return self.mtx + self.mtx.T - np.diag(self.mtx.diagonal())

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

        if self.tracer[i][j] == 1:
            l_row, l_col = i - 1, j - 1
        elif self.tracer[i][j] == 2:
            l_row, l_col = i - 1, j
        elif self.tracer[i][j] == 3:
            l_row, l_col = i, j - 1
        elif self.tracer[i][j] == 4:
            l_row, l_col = None, None
        else:
            raise TracerError('Wrong tracer value for {} {}'.format(i, j))

        return self.recursive_path(l_row, l_col) + [Node(i, j, self.mtx[i][j]), ]

    def find_ro_paths(self):
        """Find Reverse Overlap Regions"""
        ro_paths = {}
        idx = self._len
        tmp_score = self.mtx[idx][self._len]

        while idx > 1:
            if self.mtx[idx][self._len] < tmp_score or self.mtx[idx - 1][self._len] > self.mtx[idx][self._len]:
                tmp_score = self.mtx[idx, self._len]
                idx -= 1
                continue

            tmp_path = self.recursive_path(idx, self._len)
            idx -= 1
            if self.mtx[idx][self._len] < 3:
                continue

            s_node, e_node = tmp_path[0], tmp_path[-1]
            # Path start at head line
            if s_node.row != 0:
                continue
            if s_node.col not in ro_paths or ro_paths[s_node.col].score < e_node.score:
                ro_paths[s_node.col] = e_node

        return ro_paths


def get_kmers(in_seq, k=31, step=7):
    kmer_seq = in_seq
    kmers = []
    for i in range(0, len(kmer_seq) - k, step):
        kmers.append(kmer_seq[i:min(i + k, len(kmer_seq))])
    return kmers


def multiple_sequence_alignment(fasta, match=1, mismatch=-1, gap=-2, globalAlign=0, simple=0):
    import poagraph
    import seqgraphalignment
    graph = poagraph.POAGraph(fasta[0][1], fasta[0][0])
    for label, sequence in fasta[1:]:
        alignment = seqgraphalignment.SeqGraphAlignment(sequence, graph, fastMethod=not simple,
                                                        globalAlign=globalAlign,
                                                        matchscore=match, mismatchscore=mismatch,
                                                        gapscore=gap)
        graph.incorporateSeqAlignment(alignment, sequence, label)

    return graph.allConsenses()


def find_ro(raw_seq, k=31):
    from preprocess import trim_primer
    # Trim sequence
    trimmed_seq = trim_primer(raw_seq)
    if len(trimmed_seq) <= 50:
        return None, None

    # Kmer graph
    step = k // 4
    kmers = get_kmers(trimmed_seq, k=k, step=step)
    ro_finder = ROFinder(kmers, match=1, mismatch=-1, gap=-2)
    ro_paths = ro_finder.find_ro_paths()

    # Adjust accurate junction site
    seed = ro_finder.sequence[0]
    junc_sites = [0]
    for start, end_node in ro_paths.items():
        if start == 0:
            continue
        align, score, pos = local_pairwise_align_ssw(DNA(seed),
                                                     DNA(trimmed_seq[start * step - 10: start * step + k - 10]))
        if score > 30:
            junc_sites.append(start * step - 10 + pos[1][0] - pos[0][0])

    # Split last tail if possible
    if len(junc_sites) == 1:
        return None, None

    mean_length = int(np.mean([j - i for i, j in zip(junc_sites[:-1], junc_sites[1:])]))
    tail_length = len(trimmed_seq) - junc_sites[-1]
    if tail_length > mean_length + 15:
        last_seq = trimmed_seq[junc_sites[-1] + mean_length - 15:]
        last_align, last_score, last_pos = local_pairwise_align_ssw(DNA(seed), DNA(last_seq))
        if last_score > 30:
            junc_sites.append(junc_sites[-1] + mean_length - 15 + last_pos[1][0] - last_pos[0][0])

    fasta = []
    for i, j in zip(junc_sites, junc_sites[1:] + [len(trimmed_seq), ]):
        fasta.append(('{}-{}'.format(i, j), trimmed_seq[i:j]))

    segments = ','.join([segment_id for segment_id, _ in fasta])
    ccs = ''.join(multiple_sequence_alignment(fasta)[0][1])
    return segments, ccs


def worker(chunk):
    ret = []
    for header, seq in chunk:
        segments, ccs = find_ro(seq)
        ret.append((header, segments, ccs))
    return ret


def find_ccs_reads(in_file, out_file, threads):
    import gzip
    from multiprocessing import Pool
    from utils import to_str
    from logger import ProgressBar
    pool = Pool(threads)
    jobs = []

    total_cnt = 0
    with gzip.open(in_file, 'rb') as fastq:
        chunk = []
        chunk_size = 250
        chunk_cnt = 0
        for line in fastq:
            total_cnt += 1
            header = to_str(line.rstrip())
            seq = to_str(fastq.readline().rstrip())
            separator = to_str(fastq.readline())
            qual = to_str(fastq.readline())

            # Split into chunks
            chunk.append((header, seq))
            chunk_cnt += 1
            if len(chunk) == chunk_size:
                jobs.append(pool.apply_async(worker, (chunk,)))
                chunk = []
    if len(chunk) > 0:
        jobs.append(pool.apply_async(worker, (chunk,)))
    pool.close()

    prog = ProgressBar()
    prog.update(0)
    finished_chunk = 0
    total_reads = 0
    ro_reads = 0

    with open(out_file, 'w') as out:
        for job in jobs:
            ret = job.get()
            for header, segments, ccs in ret:
                total_reads += 1
                if segments is None and ccs is None:
                    continue
                ro_reads += 1
                out.write('{}\t{}\t{}\n'.format(header, segments, ccs))
            finished_chunk += 1
            prog.update(100 * finished_chunk // chunk_cnt)
    prog.update(100)
    pool.join()

    return total_reads, ro_reads
