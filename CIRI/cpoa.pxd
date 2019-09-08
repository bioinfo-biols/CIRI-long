from libcpp.vector cimport vector
from libcpp.string cimport string


cdef extern from "cpoa.h":
    char *fasta_to_ccs(vector[string] sequences);

