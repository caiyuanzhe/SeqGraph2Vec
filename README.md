# Official implementation of SeqGraph2Vec (+HNSW)

**This repository is an official Python implementation of**  ["Fast similarity search for DNA sequences through
de Bruijn sum graph embedding"]

Similarity search of DNA sequences is widely used in many genomic analyses, such as pathogen detection, gene function annotation, and evolutionary relationship discovery. Today, next-generation sequencing (NGS) technologies are generating more and more DNA sequences. This requires more efficient sequence search methods that scale well to large sequence databases. Here, we present a new efficient DNA search algorithm that scales well to large data as shown in experiments. This algorithm involves three innovative techniques:
(i) de Bruijn sum graph, which is a natural representation of multiple DNA sequences, (ii) sampling from equilibrium distribution instead of traditional uniform distribution, and (iii) a technique to solve the sink difficulty in random walk sampling on the directed graph. 

A sequence corresponds to a path on de Bruijn sum graph, which generates a vector pooled from the path node embedding vectors. A query sequence similarly corresponds to a vector embedding. Thus, the similarity search becomes a vector data search, which can be implemented very efficiently in the vector database.

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

## Support
For support, please consider opening a GitHub issue and we will do our best to reply in a timely manner.
Alternatively, if you would like to keep the conversation private, feel free to contact Yuanzhe Cai at caiyuanzhe@sztu.edu.cn

## License


## Authors
Zhaochong Yu<sup>2</sup>, Zihang Yang a, Chris Ding b, Feijuan Huang c,* and Yuanzhe Cai a,*

a Shenzhen Technology University, Shenzhen 518118, China, b The Chinese University of Hong Kong, Shenzhen, Shenzhen 518172, China, c The Second People’s Hospital of Shenzhen, Shenzhen 518037, China.
## Funding
This research work was financially supported by Shenzhen Second People's Hospital COVID-19 Emergency Clinical Research Project (2023xgyj3357009), Shenzhen Science and Technology Program (20231127194506001), the Shenzhen Science and Technology Innovation Committee Funds (JSGG20220919091404008), the National Natural Science Foundation of China-81804154 and 2023 Special Fund for Science and Technology Innovation Strategy of Guangdong Province (Science and Technology Innovation Cultivation of College Students) (Pdjh2023b0468).
