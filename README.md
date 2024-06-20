# Official implementation of SeqGraph2Vec (+HNSW)

**This repository is an official Python implementation of**  ["Fast similarity search for DNA sequences through
de Bruijn sum graph embedding"]


## Usage
- To reproduce the experiments in our paper, run the following script:
```sh 
pip3 install -r requirements.txt  # Python Version: 3.8
cd examples 
```
```sh
python3 1_train_kmer_vector.py --input_data_dir='../data_dir/input/1_data_for_kmer_vector_training/small_data/' --path_to_kmer_embedding_file='../data_dir/input/1_data_for_kmer_vector_training/small_data/kmer-embedding.txt' --kmer_size=8 --dataprocess_workers=8 --seq_file_num_to_load=8 --pagerank_damping_factor=0.85 --p=1.0 --q=0.001 --damping_factor_for_teleportation=0.99 --num_walks=40 --walks_length=150 --kmer_vec_dimension=128 --skip_gram_workers=8
```
```sh
python3 2_train_dna_vector.py --work_dir='../data_dir/input/2_data_for_seq_search/' --kmer2vec_file='../data_dir/input/2_data_for_seq_search/kmer-embedding.txt' --kmer_size=8 --ref_segment_length=150 --query_segment_number=2
```
```sh
python3 3_dna_vector_search.py --work_dir='../data_dir/input/2_data_for_seq_search/' --vertex_connection=100 --ef_search=2000 --ef_construction=128
```

Note that: 
(i) due to the Github storage limitation, currently only toy data set is provided in the repository (full data sets used in our paper can be downloaded from https://www.ncbi.nlm.nih.gov/), (ii) please fine-tune the hyperparameters in "1_train_kmer_vector.py", "2_train_dna_vector.py" and "3_dna_vector_search.py" under "examples/ directory", according to own your data sets, and (iii) the sequence search result (i.e., search accuracy and speed) can be found in "data_dir/input/2_data_for_seq_search/faiss-result.log".
