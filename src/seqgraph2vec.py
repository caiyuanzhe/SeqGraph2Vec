# -*- coding: utf-8 -*-
import sys
sys.path.extend(['.', '..'])
import numpy as np
import networkx as nx
from gensim.models import Word2Vec
from src import walker


class SeqGraph2Vec:
    """ Save a txt file recording all k-mers' related vectors. """

    def __init__(
        self,
        p: float = 1.0,
        q: float = 0.001,
        dimensions: int = 128,
        num_walks: int = 40,
        walks_length: int = 150,
        window: int = 10,
        min_count: int = 1,
        epochs: int = 1,
        workers: int = 4,
        alpha: float = 0.99,
        verbose: bool = True,
        pgr: dict = {},
    ):
        self.p = p
        self.q = q
        self.dimensions = dimensions
        self.num_walks = num_walks
        self.walks_length = walks_length
        self.window = window
        self.min_count = min_count
        self.epochs = epochs
        self.verbose = verbose
        self.workers = workers
        self.alpha = alpha
        self.pgr = pgr
        self.node_idmap = None
    

    def _read_graph(self, edge_list_path: str = 'networkfile.edg', extend: bool = False):
        walker_mode = getattr(walker, 'SparseOTF', None)
        graph = walker_mode(
            p=self.p,
            q=self.q,
            workers=self.workers,
            verbose=self.verbose,
            extend=extend,
            random_state=None
        )
        self.node_idmap = graph.read_edg(edge_list_path, weighted=True, directed=True)
        return graph

    def _simulate_walks(self, graph: nx.classes.graph.Graph) -> walker.SparseOTF:
        print('!!!Start to simulate random walks...')
        return graph.simulate_walks(self.num_walks, self.walks_length,  alpha=self.alpha, 
                                    pgr=self.pgr, node_idmap=self.node_idmap)

    def _learn_embeddings(
        self,
        walks: walker.SparseOTF,
        path_to_embeddings_file: str,
    ):
        print('!!!Start to learn kmer embeddings...')
        model = Word2Vec(
            walks,
            vector_size=self.dimensions,
            window=self.window,
            min_count=self.min_count,
            sg=1,
            workers=self.workers,
            epochs=self.epochs,
        )
        
        output_fp = path_to_embeddings_file
        if output_fp.endswith(".npz"):
            np.savez(output_fp, IDs=model.wv.index_to_key, data=model.wv.vectors)
        else:
            model.wv.save_word2vec_format(output_fp)

    def fit(
        self,
        path_to_edg_list_file: str,
        path_to_embedding_file: str,
    ):
        """ Get embeddings of k-mers fragmented from input sequences.

        Args:
            path_to_edg_list_file (str) : path to k-mers' edges list file.
            path_to_embeddings_file (str) : path to k-mers' embeddings file.
        """
        graph = self._read_graph(edge_list_path=path_to_edg_list_file)
        walks = self._simulate_walks(graph)
        self._learn_embeddings(walks, path_to_embedding_file)


