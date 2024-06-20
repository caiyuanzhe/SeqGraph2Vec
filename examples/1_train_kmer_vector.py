# -*- coding: utf-8 -*-
import sys
sys.path.extend(['.', '..'])
import os
import re
import gc
import multiprocessing
from collections import Counter
import networkx as nx
from Bio import SeqIO

from src.seqgraph2vec import SeqGraph2Vec
from src.cython_function import fast_compute_edge_weight



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


class KMerEmbeddings:

    def __init__(
        self,
        p: float,
        q: float,
        dimensions: int,  
        workers: int,
        path_to_edg_list_file: str,
        kmer_vec_output_dir: str,
        pgr: dict,
    ):
        self.p = p
        self.q = q
        self.workers = workers
        self.dimensions = dimensions
        self.path_to_edg_list_file = path_to_edg_list_file
        self.kmer_vec_output_dir = kmer_vec_output_dir
        self.pgr = pgr

    def train(self):
        """ Obtain the k-mer embedding. """
        print(self.path_to_edg_list_file)
        clf = SeqGraph2Vec(p=self.p, q=self.q, workers=self.workers, dimensions=self.dimensions, pgr=self.pgr)
        clf.fit(
            path_to_edg_list_file=self.path_to_edg_list_file,
            path_to_embeddings_file=self.kmer_vec_output_dir + f"kmer-node2vec-embedding.txt",
        )


def compute_edge_weight(seqs, global_weight_dict, mer, process_id):
    weight_dict = fast_compute_edge_weight(seqs, mer)
    global_weight_dict[process_id] = weight_dict


def split_list(lst, n):
    k, m = divmod(len(lst), n)
    return list(lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


if __name__ == '__main__':

    ###### Step1: train the kmer embedding, you can skip this step if you have pre-trained kmer embedding. ######    
    work_dir = '../data_dir/input/1_data_for_kmer_vector_training/small_data/'  # set this, all input data sets should be placed in one directory. An example of input data is given in "data_dir/input/small_data/toy_seqs.fna".
    edge_list_path = work_dir + 'networkfile.edg'
    mer = 8
    dataprocess_workers = 8
    seq_file_num_to_load = 10  # could set a large value for parallel loading, if possible 
    alpha_in_pagerank = 0.85

    seq_files = list()
    for input_file_dir in [work_dir]:
        print('input_file_dir', input_file_dir)
        for root, dirs, files in os.walk(input_file_dir):
            for file in files:
                if file.endswith('.fna'):
                    seq_files.append(os.path.join(root, file))
    

    all_edge_weight_dict = dict()
    for i in range(0, len(seq_files), seq_file_num_to_load):

        # Load 'seq_file_num_to_load' fna files in memory each time
        files = seq_files[i:i+seq_file_num_to_load]
        seqs = list()
        for f in files:
            for seq_record in SeqIO.parse(f, 'fasta'):
                seq = re.sub('[^ACGTacgt]+', '', str(seq_record.seq))
                seqs.append(seq.upper())

        manager = multiprocessing.Manager()  # multiprocess for computing edge weights
        edge_weight_dict = manager.dict()
        processes = []
        for i,partition_seqs in enumerate(split_list(seqs, n=dataprocess_workers)):
            p = multiprocessing.Process(
                target=compute_edge_weight,
                args=(partition_seqs, edge_weight_dict, mer, i)
            )
            processes.append(p)
            p.start()
        for p in processes:
            p.join()

        tmp = Counter({})  # merge multiple edge weight dictionaries into single one
        for i in range(0, dataprocess_workers):
            tmp.update(Counter(edge_weight_dict[i]))
        all_edge_weight_dict.update(dict(tmp))

        # memory recycling
        del tmp
        del edge_weight_dict
        del seqs
        del p
        del partition_seqs
        del processes
        del manager
        leak = gc.collect()
        print("gc.collect() returned", leak)

    edge_list = [(nodes[0], nodes[1], weight) for nodes, weight in all_edge_weight_dict.items()]
    with open(edge_list_path, 'w', encoding='utf-8') as edge_list_file:
        for edge_pair in edge_list:
            write_content = str(edge_pair[0]) + '\t' + str(edge_pair[1]) + '\t' + str(edge_pair[2]) + '\n'
            edge_list_file.write(write_content)
    
    graph = nx.DiGraph()
    for nodes, weight in all_edge_weight_dict.items():
        graph.add_edge(
            nodes[0], nodes[1],
            weight=weight
        )
    # print(graph.edges.data("weight"))
    pgr = nx.pagerank(graph, alpha=alpha_in_pagerank)

    """ Use random walk sampling and Skip-Gram to learn kmer embedding """
    clf = KMerEmbeddings(
        kmer_vec_output_dir=work_dir,
        path_to_edg_list_file=edge_list_path,  # default setting
        p=1.0,
        q=0.001,
        dimensions=128,
        workers=8,
        pgr=pgr
    )
    clf.train()
 
