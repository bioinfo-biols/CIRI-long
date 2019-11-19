# cython: language_level=3, boundscheck=False
cimport cpoa


def consensus(fasta, alignment_type=1, match=10, mismatch=-4,
              gap=-8, extension=-2, gap_affine=-24, extension_affine=-1, debug=0):
    """
    Generate consensus sequence using libspoa (https://github.com/rvaser/spoa)
    :param fasta: tuples of input sequences
    :param alignment_type: alignment mode, 0:local, 1:global, 2:semi-global [1]
    :param match: score for matching bases [5]
    :param mismatch: score for mismatching bases [-4]
    :param gap: gap opening penalty (non-positive) [-8]
    :param extension: gap extension penalty (non-positive) [-6]
    :param gap_affline: gap opening penalty of the second affine function (non-positive) [-10]
    :param extension_affline: gap extension penalty of the second affine function (non-positive) [-4]
    :param debug: whether output multiple sequence alignment [0]
    :return: string of consensus sequence
    """

    sequences = [seq.encode('utf-8') for name, seq in sorted(fasta, key=lambda x: len(x[1]), reverse=True)]
    cdef int l = alignment_type
    cdef int m = match
    cdef int n = mismatch
    cdef int g = gap
    cdef int e = extension
    cdef int q = gap_affine
    cdef int c = extension_affine
    cdef int v = debug

    cdef char *s = cpoa.fasta_to_ccs(sequences, l, m, n, g, e, q, c, v)
    py_string = s.decode('utf-8')
    return py_string

