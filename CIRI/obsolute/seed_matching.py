#!/home/zhangjy/.virtualenvs/Benchmarking/bin/python
import poagraph
import seqgraphalignment
#
# from skbio.alignment import local_pairwise_align_ssw
# from skbio import TabularMSA
# from skbio import DNA
#
# query_sequence = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGTCCGTCCGTCCGTCCTGTGGGGATGGCAGAGACATCTTGACCACGGCGTCCGCTTCTACGTGTGTGCTTTCGTGGGTGACATCACTGTGGGCGGCGGCAGTGACCGGTGGACGCTGACGGCGTGCGCGCGTCCGTCCGTCCGTCCAGCGGGGCGAGGGTCGGTGGCACATACATCTTGACCGCGGCGGTGAACATCTTGACCGCGGCGGCTCGTCCGTGTCTGCGTGGCGGGTGCGGGTCCTGCGGATGCGTGTGCTGGCCTGCCTGCCGTGTGTGCCGTGTGTGCCGCGTGTGCCTGGGCGGGCGCGGGGCCCGCGTGTGTCCTGTGTGCTGTGGTGGGTGCGGCGGTGCGCTCGGCGCAGGCTTCCTACCTTGTGGAGTCCGTCCGTGCGTCCGTGCGTCCGTGCGTCCCGTGGGCATGGCAGAGACATCTTGACCGCGGCGTCCGCTTCTGCGCGTGTGGGTGACGTCGCTGGGGGCGGCGGCCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTTTGCCGTCTTCTGCTTG'
#
# print('Sequence Length: {}'.format(len(query_sequence)))
# print()
#
# # print('-- Adapter Alignment --')
# # adapter = 'ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT'
# # alignment, score, positions = local_pairwise_align_ssw(DNA(adapter), DNA(query_sequence))
# # print(alignment)
# # print(score)
# # print(positions)
# # print()
#
# print('-- Seed-and-Extend --')
# junc = 200
# while True:
#     alignment, score, positions = local_pairwise_align_ssw(DNA(query_sequence[:junc]), DNA(query_sequence[junc:]))
#     print(junc, positions)
#     (l_st, l_en), (r_st, r_en) = positions
#     if r_st != 0:
#         junc = junc + r_st
#     else:
#         fasta = list()
#         fasta.append(['seq1', query_sequence[l_st:junc]])
#         fasta.append(['seq2', query_sequence[junc:junc + r_en]])
#         break

fasta = [('My', 'GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'),
         ('7', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG'),
         ('8', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG'),
         ('9', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG'),
         ('10', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG'),
         ]
print('-- Generate Consensus --')
match = 1
gap = -1
mismatch = -1
globalAlign = 1
simple = 0

graph = poagraph.POAGraph(fasta[0][1], fasta[0][0])
for label, sequence in fasta[1:]:
    alignment = seqgraphalignment.SeqGraphAlignment(sequence, graph, fastMethod=not simple,
                                                    globalAlign=globalAlign,
                                                    matchscore=match, mismatchscore=mismatch,
                                                    gapscore=gap)
    graph.incorporateSeqAlignment(alignment, sequence, label)

alignments = graph.generateAlignmentStrings()
for label, alignstring in alignments:
    print("{0:15s} {1:s}".format(label, alignstring))
print()

# print(query_sequence[:l_st])
# print(query_sequence[junc + r_en:])
# print('-- Output Consensus Sequence --')
# consenses = graph.allConsenses()
# for i, consensus in enumerate(consenses):
#     print('Consensus'+str(i))
#     print(''.join(consensus[1]))
#
# print(query_sequence[55:junc])