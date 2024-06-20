# cython: language_level=3
from tqdm import tqdm
cimport cython

@cython.boundscheck(False)  
@cython.wraparound(False)   
cpdef inline extract_kmer(str seq, int mer):
    return [seq[i:i+mer] for i in range(len(seq) - mer + 1)]

@cython.boundscheck(False)  
@cython.wraparound(False)   
cpdef dict fast_compute_edge_weight(list seqs, int mer):
    cdef dict weight_dict = dict()
    cdef list k_mers
    cdef int i
    for seq in seqs:
        k_mers = extract_kmer(seq, mer)
        for i in range(len(k_mers) - 1):
            if (k_mers[i], k_mers[i+1]) in weight_dict:
                weight_dict[(k_mers[i], k_mers[i+1])] += 1
            else:
                weight_dict[(k_mers[i], k_mers[i+1])] = 1
    return weight_dict
