# Fast similarity search for DNA sequences through de Bruijn sum graph embedding

Similarity search of DNA sequences is widely used in many genomic analyses. Today, next-generation sequencing (NGS) technologies are generating more and more DNA sequences. This requires more efficient sequence search methods that scale well to large sequence databases. Here, we present a new efficient DNA search algorithm, SeqGraph2Vec (+HNSW), that scales well to large data as shown in experiments, and provide Python implementation in this repository. 

## Installation

Install latest version with:
```bash
$ git clone https://github.com/caiyuanzhe/SeqGraph2Vec.git
$ cd SeqGraph2Vec
$ pip3 install -r requirements.txt
```

## Usage
- Guide to reproduce the experiments in our paper.

### Training k-mer embedding vectors

To train k-mer embedding vectors using *SeqGraph2Vec*, execute the following command from the *examples/* directory. 
Before running the command line, all input DNA sequence files should be placed in *input_data_dir* in advance. An example of input data file is "data_dir/input/small_data/toy_seqs.fna".
After running the command line, a k-mer embedding vector file is produced --- *path_to_kmer_embedding_file*.

```bash
cd examples
python3 1_train_kmer_vector.py --input_data_dir='../data_dir/input/1_data_for_kmer_vector_training/small_data/' --path_to_kmer_embedding_file='../data_dir/input/1_data_for_kmer_vector_training/small_data/kmer-embedding.txt' --kmer_size=8 --dataprocess_workers=8 --seq_file_num_to_load=8 --pagerank_damping_factor=0.85 --p=1.0 --q=0.001 --damping_factor_for_teleportation=0.99 --num_walks=40 --walks_length=150 --kmer_vec_dimension=128 --skip_gram_workers=8
```

### Computing DNA sequence embedding vectors

Once obtaining the k-mer embedding vectors, for a given input DNA sequence, the embedding vectors of its constituent k-mers are averaged to produce the sequence's embedding vector.
Before running the command line, several files including *the pretrained k-mer embedding vector file* (e.g. *path_to_kmer_embedding_file* in the above command line), and *the DNA input sequence files* that need to be vectorized should be placed in *work_dir* in advance. 
After running the command line, reference segments and query segments will be extracted from *the DNA input sequence files* as well as their corresponding sequence embedding vectors.

```bash
python3 2_train_dna_vector.py --work_dir='../data_dir/input/2_data_for_seq_search/' --path_to_kmer_embedding_file='../data_dir/input/2_data_for_seq_search/kmer-embedding.txt' --kmer_size=8 --ref_segment_length=150 --query_segment_number=2
```

### DNA similarity search as DNA vector search

DNA sequence search is transformed into a vector data search, which can be implemented very efficiently in vector search technology such as Faiss. 
Before running the command line, make sure that all the produced files in the previous step (i.e., Computing DNA sequence embedding vectors) is placed in *work_dir* in advance. 

```bash
python3 3_dna_vector_search.py --work_dir='../data_dir/input/2_data_for_seq_search/' --vertex_connection=100 --ef_search=2000 --ef_construction=128
```

### Options

Check out the full list of options available using:
```bash
python3 1_train_kmer_vector.py --help
python3 2_train_dna_vector.py --help
python3 3_dna_vector_search.py --help
```

Note: due to the Github storage limitation, currently only toy data set is provided in the repository (full data sets used in our paper can be downloaded from https://www.ncbi.nlm.nih.gov/).