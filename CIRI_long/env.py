ALIGNER = None
CONTIG_LEN = None
GENOME = None
GTF_INDEX = None
INTRON_INDEX = None
SS_INDEX = None


def initializer(aligner, contig_len, genome, gtf_index, intron_index, ss_index):
    global ALIGNER
    global CONTIG_LEN
    global GENOME
    global GTF_INDEX
    global INTRON_INDEX
    global SS_INDEX
    ALIGNER = aligner
    CONTIG_LEN = contig_len
    GENOME = genome
    GTF_INDEX = gtf_index
    INTRON_INDEX = intron_index
    SS_INDEX = ss_index


