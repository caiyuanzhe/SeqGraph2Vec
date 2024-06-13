# -*- coding: utf-8 -*-
import sys
sys.path.extend(['.', '..'])
import os
import re
from Bio import SeqIO

from util.faiss_getprecision import create_index
from util.faiss_getprecision import accuracy
from util.perf_tools import Tee


def parse_seq(path_to_input: str):
    """ Return a list containing DNA seqment(s) captured in fna file(s)."""
    seq_files = list()
    for input_file_dir in path_to_input:
        print(input_file_dir)
        for root, dirs, files in os.walk(input_file_dir):
            for file in files:
                if file.endswith('.fna'):
                    seq_files.append(os.path.join(root, file))
    seqs = list()
    for seq_file in seq_files:
        for seq_record in SeqIO.parse(seq_file, 'fasta'):
            seq = re.sub('[^ACGTacgt]+', '', str(seq_record.seq))
            seqs.append(seq.upper())

    print('There are ' + str(len(seqs)) + ' seqs')

    return seqs


class SequenceRetrieval:

    def __init__(
        self,
        ref_segment_file: str,
        ref_segment_vec_file: str,
        query_ref_segment_file: str,
        query_segment_vec_file: str,
        faiss_index_file: str,
        faiss_log: str,
        top_kn: int,
    ):
        """ Note: we use segments from which ramdomly selected subsegments 
            originally come to take the place of subsegments, in order to 
            make comparisons with segments in corpus . """
        self.ref_segment_file = ref_segment_file
        self.ref_segment_vec_file = ref_segment_vec_file
        self.query_ref_segment_file = query_ref_segment_file
        self.query_segment_vec_file = query_segment_vec_file

        self.faiss_index_file = faiss_index_file
        self.faiss_log = faiss_log
        self.top_kn = top_kn

    def train(
        self,
        dimension: int,
        index_method: str,
        vertex_connection: int,
        ef_search: int,
        ef_construction: int,
    ):
        """ Print the Top-K result for sequence retrieval task """
        logger = Tee(self.faiss_log)
        sys.stdout = logger

        if create_index(
            self.ref_segment_vec_file,
            self.faiss_index_file,
            dimension,
            index_method,
            vertex_connection,
            ef_search,
            ef_construction,
        ):
            accuracy(
                self.query_segment_vec_file,
                self.query_ref_segment_file,
                self.ref_segment_file,
                self.faiss_index_file,
                self.top_kn
            )


if __name__ == '__main__':
    work_dir = '../data_dir/input/2_data_for_seq_search/'

    clf = SequenceRetrieval(
        ref_segment_file=work_dir+'ref_segment.txt',  # default name
        ref_segment_vec_file=work_dir+'ref_segment_vecs.txt',  # default name
        query_ref_segment_file=work_dir+'query_ref_segment.txt',  # default name
        query_segment_vec_file=work_dir+'query_segment_vecs.txt',  # default name
        faiss_index_file=work_dir+'faiss-index-file',  # default name
        faiss_log=work_dir+'faiss-result.log',  # default name
        top_kn=2,
    )
    clf.train(
        dimension=128,  # set this, consistent with dimensionality of sequence embedding
        index_method='HNSW',  # default 
        vertex_connection=100,  # default 
        ef_search=2000,  # default 
        ef_construction=128,  # default 
    )
