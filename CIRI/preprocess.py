def revcomp(seq):
    trantab = str.maketrans("ATCG", "TAGC")
    return seq.translate(trantab)[::-1]

