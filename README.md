# Official implementation of SeqGraph2Vec (+HNSW)

**This repository is an official Python implementation of**  ["Fast similarity search for DNA sequences through
de Bruijn sum graph embedding"]


## Usage

- To reproduce the experiments in our paper, run the following script:
```sh 
pip3 install -r requirements.txt  # Python Version: 3.8
```
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
