# -*- coding: utf-8 -*-
import sys
sys.path.extend(['.', '..'])
import argparse
from util.faiss_getprecision import create_index
from util.faiss_getprecision import accuracy
from util.perf_tools import Tee


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
    parser = argparse.ArgumentParser()
    parser.add_argument('--work_dir', type=str, default='../data_dir/input/2_data_for_seq_search/', 
                        help='the output dir where exists reference/query segments and their embedding vectors. Search result is also produced in this dir.')  
    parser.add_argument('--vertex_connection', type=int, default=100, help='number of connections each vertex will have') 
    parser.add_argument('--ef_search', type=int, default=2000, help='depth of layers explored during search') 
    parser.add_argument('--ef_construction', type=int, default=128, help='depth of layers explored during index construction') 
    args = parser.parse_args()
    print(args)
    work_dir = args.work_dir

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
        index_method='HNSW', 
        vertex_connection=args.vertex_connection,  
        ef_search=args.ef_search, 
        ef_construction=args.ef_construction, 
    )
