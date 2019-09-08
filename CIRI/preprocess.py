def revcomp(seq):
    trantab = str.maketrans("ATCG", "TAGC")
    return seq.translate(trantab)[::-1]


def trim_primer(seq, primer='AAGCAGTGGTATCAACGCAGAGTAC', last_bases=150):
    from skbio import DNA
    from skbio.alignment import local_pairwise_align_ssw
    l_alignments, l_score, l_positions = local_pairwise_align_ssw(DNA(seq[:last_bases]), DNA(primer))
    l_trim = l_positions[0][1] + 1 if l_score >= 25 else 0

    r_alignments, r_score, r_positions = local_pairwise_align_ssw(DNA(seq[-last_bases:]), DNA(revcomp(primer)))
    r_trim = -last_bases + r_positions[0][0] if r_score >= 25 else len(seq)
    return seq[l_trim:r_trim]


def partial_order_alignment(fasta, match=2, mismatch=-3, gap=-2, globalAlign=0, simple=0):
    import poagraph
    import seqgraphalignment
    graph = poagraph.POAGraph(fasta[0][1], fasta[0][0])

    if len(fasta) > 2:
        for label, sequence in fasta[1:-1]:
            alignment = seqgraphalignment.SeqGraphAlignment(sequence, graph, fastMethod=not simple,
                                                            globalAlign=1,
                                                            matchscore=match, mismatchscore=mismatch,
                                                            gapscore=gap)
            graph.incorporateSeqAlignment(alignment, sequence, label)

    alignment = seqgraphalignment.SeqGraphAlignment(fasta[-1][1], graph, fastMethod=not simple,
                                                    globalAlign=0,
                                                    matchscore=match, mismatchscore=mismatch,
                                                    gapscore=gap)
    graph.incorporateSeqAlignment(alignment, fasta[-1][1], fasta[-1][0])
    return graph
