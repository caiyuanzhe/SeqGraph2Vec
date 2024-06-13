SeqGraph2Vec
====================================
SeqGraph2Vec is an open-source library to train fast and high-quality k-mer embedding from the k-mer graph. 
Within the k-mer embedding, we can easily compute the DNA sequence embedding. In our paper, the DNA sequence embedding is an average of the k-mer embedding.


------------------------------------

### Requirements
The codebase is implemented in Python 3.8 package versions. 

pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple


### Datasets
<p align="justify">
The code takes FASTA format files with file extension of **.fna**. Note that all training FASTA format files should be under the same input directory. A sample FASTA format file is included in the  `data_dir/input/` directory. </p>
<p align="justify">
Training the model is handled by the `src/cli.py` script which provides the following command line arguments.</p>
