#!/usr/bin/env python
import sys

def test_ssw():
    with open('./tests/test.fa', 'r') as f:
        f.readline()
        seq1 = f.readline().rstrip()
        f.readline()
        seq2 = f.readline().rstrip()

    # from skbio import local_pairwise_align_ssw, DNA
    from lib.striped_smith_waterman.ssw_wrap import Aligner
    ssw = Aligner(str(seq1), match=1, mismatch=1, gap_open=1, gap_extend=1)
    res = ssw.align(str(seq2), 0, 19)
    sys.stderr.write(str(res))


if __name__ == '__main__':
    test_ssw()



# print('finished')