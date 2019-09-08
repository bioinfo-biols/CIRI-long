cimport cpoa


def consensus(fasta):
    sequecnes = [seq.encode('utf-8') for name, seq in fasta]
    cdef char *s = cpoa.fasta_to_ccs(sequecnes)
    py_string = s
    return py_string

