# -*- coding: utf-8 -*-
import sys
sys.path.extend(['.', '..'])
import numpy as np
import networkx as nx
from gensim.models import Word2Vec
from src import walker


class SeqGraph2Vec:
    """ Save a txt file recording all k-mers' related vectors.

    Args:
        p (float) : return parameter, optional (default = 1)
            The value less than 1 encourages returning back to
            previous vertex, and discourage for value grater than 1.
        q (float) : in-out parameter, optional (default = 0.001)
            The value less than 1 encourages walks to
            go "outward", and value greater than 1
            encourage walking within a localized neighborhood.
        dimensions (int) : dimensionality of the word vectors
            (default = 128).
        num_walks (int): number of walks starting from each node
            (default = 10).
        walks_length (int): length of walk
            (default = 80).
        window (int) : Maximum distance between the current and
            predicted k-mer within a sequence (default = 10).
        min_count (int) : Ignores all k-mers with total frequency
            lower than this (default = 1)
        epochs : Number of iterations (epochs) over the corpus
            (default = 1)
        workers (int) :  number of threads to be spawned for
            runing node2vec including walk generation and
            word2vec embedding, optional (default = 4)
        verbose (bool) : Whether or not to display walk generation
            progress (default = True).

    """

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

    def _simulate_walks(self, graph: nx.classes.graph.Graph, alpha: float = 0.99) -> walker.SparseOTF:
        print('!!!Start to simulate random walks...')
        return graph.simulate_walks(self.num_walks, self.walks_length,  alpha=alpha, 
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
        path_to_embeddings_file: str,
    ):
        """ Get embeddings of k-mers fragmented from input sequences.

        Args:
            graph (nx.classes.graph.Graph) : nx.DiGraph() object.
            mer (int) : sliding window length to fragment k-mers.
                slide only a single nucleotide.
            path_to_edg_list_file (str) : path to k-mers' edges list file.
            path_to_embeddings_file (str) : path to k-mers' embeddings file.
        """
        graph = self._read_graph(edge_list_path=path_to_edg_list_file)
        walks = self._simulate_walks(graph)
        self._learn_embeddings(walks, path_to_embeddings_file)


