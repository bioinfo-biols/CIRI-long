def revcomp(seq):
    trantab = str.maketrans("ATCG", "TAGC")
    return seq.translate(trantab)[::-1]


def trim_primer(seq, primer='AAGCAGTGGTATCAACGCAGAGTAC'):
    from skbio import DNA
    from skbio.alignment import local_pairwise_align_ssw
    l_alignments, l_score, l_positions = local_pairwise_align_ssw(DNA(seq[:100]), DNA(primer))
    l_trim = l_positions[0][1] + 1 if l_score >= 25 else 0

    r_alignments, r_score, r_positions = local_pairwise_align_ssw(DNA(seq[-100:]), DNA(revcomp(primer)))
    r_trim = -100 + r_positions[0][0] if r_score >= 25 else -1
    return seq[l_trim:r_trim]
