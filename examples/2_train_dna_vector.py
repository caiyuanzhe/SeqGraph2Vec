# -*- coding: utf-8 -*-
import sys
sys.path.extend(['.', '..'])
import os
from gensim.models import KeyedVectors
from Bio import SeqIO
import re

from src.generators import seq2segs
from src.generators import seg2sentence
from src.generators import generate_query_segments
from src.generators import parse_seq

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


class SequenceEmbeddings:

    def __init__(
        self,
        mer: int,
        kmer2vec_file: str,
        seq_dir: str,
        ref_segment_length: int,
        query_segment_number: int,
        ref_segment_file: str,
        query_ref_segment_file: str,
        query_segment_file: str,
        sequence_vec_output_dir: str,

    ):
        """ segment embeddings """
        self.mer = mer  
        self.kmer2vec_file = kmer2vec_file
        self.seqs = parse_seq([seq_dir])
        self.seg_vec_output_dir = sequence_vec_output_dir

        self.ref_segment_length = ref_segment_length
        self.query_segment_number = query_segment_number

        self.ref_segment_file = ref_segment_file  
        self.query_ref_segment_file = query_ref_segment_file
        self.query_segment_file = query_segment_file

    def ref_segment_embeddings(self):

        # note that if an old one exists, we would not overwrite it.
        if not os.path.exists(self.ref_segment_file):
            seq2segs(
                self.seqs,
                self.ref_segment_length,
                self.ref_segment_file,
            )

        with open(self.ref_segment_file, 'r', encoding='utf-8') as fp:
            segs = [line.split('\n')[0] for line in fp.readlines()]
        sentences = seg2sentence(segs, self.mer)  # Tokenize

        from util.vectorizer import SeqVectorizer
        vecs = KeyedVectors.load_word2vec_format(self.kmer2vec_file)  # k-mer vectors

        clf = SeqVectorizer(vecs)
        clf.train(sentences)
        clf.save_embs_format(
            self.seg_vec_output_dir,
            f"{'ref_segment_vecs'}"
        )

    def query_segment_embeddings(self):

        if not os.path.exists(self.query_ref_segment_file) and not os.path.exists(self.query_segment_file):  # prevent overwriting
            generate_query_segments(
                self.ref_segment_file,
                self.ref_segment_length,
                self.query_segment_number,
                self.query_segment_file,
                self.query_ref_segment_file,
            )

        # load subsegments
        with open(self.query_segment_file, 'r', encoding='utf-8') as fp:
            query_segments = [line.split('\n')[0] for line in fp.readlines()]
        sentences = seg2sentence(query_segments, self.mer)  # Tokenize
        
        from util.vectorizer import SeqVectorizer
        vecs = KeyedVectors.load_word2vec_format(self.kmer2vec_file)  # k-mer2vec file
        
        clf = SeqVectorizer(vecs)
        clf.train(sentences, vector_size=128)
        clf.save_embs_format(
            self.seg_vec_output_dir,
            f"{'query_segment_vecs'}"
        )

    def train(self):
        """ Generate four types of files:
            1、One file of segments.
            2、One file of randomly extracted subsegments.
            3、One file of randomly extracted segments from which subsegments originally come.
            4、two files of vectors corresponding to segments/subsegments. 
        """
        self.ref_segment_embeddings()
        self.query_segment_embeddings()


if __name__ == '__main__':
    work_dir = '../data_dir/input/2_data_for_seq_search/'
    kmer_size = 8
    kmer2vec_file = work_dir+f"8mer-kmer-node2vec-embedding.txt"
    ref_segment_length = 150  # split original DNA sequences into numerous reference segments of length 150 bp 
    query_segment_number = 2  # extract [query_segment_number] query segments from reference segments. 

    # Note: if there already exists segment files, we won't overwrite it
    clf = SequenceEmbeddings(
        seq_dir=work_dir, 
        mer=kmer_size,  # set this, consistent with the length k of pre-trained k-mers
        kmer2vec_file=kmer2vec_file,  # set this
        ref_segment_length=ref_segment_length,  # set this
        query_segment_number=query_segment_number,  # set this
        ref_segment_file=work_dir+'ref_segment.txt',  # default name
        query_ref_segment_file=work_dir+'query_ref_segment.txt',  # default name
        query_segment_file=work_dir+'query_segment.txt',  # default name
        sequence_vec_output_dir=work_dir,
    )
    clf.train()

