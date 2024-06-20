# Official implementation of SeqGraph2Vec (+HNSW)

**This repository is an official Python implementation of**  ["Fast similarity search for DNA sequences through
de Bruijn sum graph embedding"]

## Requirements
- Python Version: 3.8
- arrow==1.2.2
- Bio==1.5.3
- faiss_cpu==1.7.3
- gensim==4.1.0
- icecream==2.1.1
- Logbook==1.5.3
- networkx==2.6.3
- numba==0.54.1
- numba_progress==0.0.2
- numpy==1.20.3
- nptyping==2.4.1
- pandas==1.5.2
- pecanpy==2.0.6
- prettytable==3.2.0
- psutil==5.9.1
- scikit_learn==1.2.1
- scipy==1.7.3
- setuptools==58.0.4
- tqdm==4.62.3

## Usage

- To reproduce the experiments in our paper, run the following script:
```sh
cd src
python3 setup.py build_ext -- inplace
```
```sh
cd examples 
sh examples/run.sh
```
Note that: 
(i) due to the Github storage limitation, currently only toy data set is provided in the repository (full data sets used in our paper can be downloaded from https://www.ncbi.nlm.nih.gov/), (ii) please fine-tune the hyperparameters in "1_train_kmer_vector.py", "2_train_dna_vector.py" and "3_dna_vector_search.py" under "examples/ directory", according to own your data sets, and (iii) the sequence search result (i.e., search accuracy and speed) can be found in "data_dir/input/2_data_for_seq_search/faiss-result.log".
