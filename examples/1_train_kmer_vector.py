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
import argparse

from src.seqgraph2vec import SeqGraph2Vec
from src.cython_function import fast_compute_edge_weight


def compute_edge_weight(seqs, global_weight_dict, mer, process_id):
    weight_dict = fast_compute_edge_weight(seqs, mer)
    global_weight_dict[process_id] = weight_dict


def split_list(lst, n):
    k, m = divmod(len(lst), n)
    return list(lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


if __name__ == '__main__':

    ###### Step1: train the kmer embedding, you can skip this step if you have pre-trained kmer embedding. ######    
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_data_dir', type=str, default='../data_dir/input/1_data_for_kmer_vector_training/small_data/',
                        help='all input data sets should be placed in this directory. An example of input data file is "data_dir/input/small_data/toy_seqs.fna')
    parser.add_argument('--path_to_kmer_embedding_file', type=str, default='../data_dir/input/1_data_for_kmer_vector_training/small_data/'+"kmer-embedding.txt",
                        help='the k-mer embedding vector file trained by SeqGraph2Vec')
    parser.add_argument('--kmer_size', type=int, default=8, help='the length of k-mers split from DNA sequences') 
    parser.add_argument('--dataprocess_workers', type=int, default=8, help='cpu-cores for multiprocessing in reading input data sets') 
    parser.add_argument('--seq_file_num_to_load', type=int, default=8, help='could set a large value for parallel loading, if possible ') 
    parser.add_argument('--pagerank_damping_factor', type=float, default=0.85, help='the damping factor in PageRank for computing node equilibrium distribution probability') 
    parser.add_argument('--p', type=float, default=1.0, help='2nd random walks') 
    parser.add_argument('--q', type=float, default=0.001, help='2nd random walks') 
    parser.add_argument('--damping_factor_for_teleportation', type=float, default=0.99,
                        help='the damping factor that controls random walks to teleport to a random node with a probability of (1-damping_factor_for_teleportation)')  
    parser.add_argument('--num_walks', type=int, default=40, 
                        help='The number of random walks each node will perform is: max((4^k * num_walks * walks_length * node_equilibrium_distribution_prob), 1)')
    parser.add_argument('--walks_length', type=int, default=150,
                        help='The number of random walks each node will perform is: max((4^k * num_walks * walks_length * node_equilibrium_distribution_prob), 1)')
    parser.add_argument('--kmer_vec_dimension', type=int, default=128, help='k-mer embedding vector dimension') 
    parser.add_argument('--skip_gram_workers', type=int, default=8, help='accelerating Skip-Gram model training') 
    args = parser.parse_args()
    print(args)


    ########### start: de Bruijn sum graph construction ###########
    edge_list_path = args.input_data_dir + 'networkfile.edg' # TMP file

    seq_files = list()
    for input_file_dir in [args.input_data_dir]:
        print('input_file_dir', input_file_dir)
        for root, dirs, files in os.walk(input_file_dir):
            for file in files:
                if file.endswith('.fna'):
                    seq_files.append(os.path.join(root, file))
    

    all_edge_weight_dict = dict()
    for i in range(0, len(seq_files), args.seq_file_num_to_load):
        # Load 'seq_file_num_to_load' fna files in memory each time
        files = seq_files[i:i+args.seq_file_num_to_load]
        seqs = list()
        for f in files:
            for seq_record in SeqIO.parse(f, 'fasta'):
                seq = re.sub('[^ACGTacgt]+', '', str(seq_record.seq))
                seqs.append(seq.upper())

        manager = multiprocessing.Manager()  # multiprocess for computing edge weights
        edge_weight_dict = manager.dict()
        processes = []
        for i,partition_seqs in enumerate(split_list(seqs, n=args.dataprocess_workers)):
            p_ = multiprocessing.Process(
                target=compute_edge_weight,
                args=(partition_seqs, edge_weight_dict, args.kmer_size, i)
            )
            processes.append(p_)
            p_.start()
        for p_ in processes:
            p_.join()

        tmp = Counter({})  # merge multiple edge weight dictionaries into single one
        for i in range(0, args.dataprocess_workers):
            tmp.update(Counter(edge_weight_dict[i]))
        all_edge_weight_dict.update(dict(tmp))

        # memory recycling
        del tmp
        del edge_weight_dict
        del seqs
        del p_
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
    ########### end: de Bruijn sum graph construction ###########
    

    ########### start: node equilibrium distribution probability computation  ###########
    graph = nx.DiGraph()
    for nodes, weight in all_edge_weight_dict.items():
        graph.add_edge(
            nodes[0], nodes[1],
            weight=weight
        )
    # print(graph.edges.data("weight"))
    pgr = nx.pagerank(graph, alpha=args.pagerank_damping_factor)
    ########### end: node equilibrium distribution probability computation  ###########


    ########### start: training k-mer embedding vectors ###########
    clf = SeqGraph2Vec(
        p=args.p, 
        q=args.q, 
        workers=args.skip_gram_workers, 
        num_walks=args.num_walks, 
        walks_length=args.walks_length, 
        alpha=args.damping_factor_for_teleportation, 
        dimensions=args.kmer_vec_dimension, 
        pgr=pgr
        )
    clf.fit(path_to_edg_list_file=edge_list_path, path_to_embedding_file=args.path_to_kmer_embedding_file)
     ########### end: training k-mer embedding vectors ###########

