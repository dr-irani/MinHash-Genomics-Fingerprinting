# CS447 Final Project:

## Repository Organization
    cs447-final-project/
      |- README.md
      |- minhash_combined.py --> source code for all MinHash implementations/modifications
      |- edit_distance_DP.py --> Needleman-Wunsch implementation (source: [Dr. Benjamin Langmead](https://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/tree/master/notebooks/))
      |- generate_data.py -->
      |- make_synth_data.py -->
      |- process_fastq.py -->
      |- profiler.py --> code for running time benchmarks for each MinHash implementation
      |- run_DP.py -->
      |- run_experiment.py --> 
      |- tune_parameters.py -->
      |- data/
         |- README.md
         |- ...
      |- notebooks/
         |- ...
      |- output/
         |- ...
      |- profile_out/
         |- ...
      |- synth_data/
         |- ...
      |- synth_out/
         |- ...
      |- tuning/
         |- ecoli/
            |- ...
         |- synth/
            |- ...
         |- ...


## MinHash Genomic Fingerprinting to Estimate Edit Distance

## Data Generation
Example commands of generating synthetic reads can be found in `/data/READMEmd`. 

## MinHash
All of our implementations for MinHash are contained in `minhash_combined.py` The following parameters are used as arguments when calling the program. 
``` 
-f1: "Path of sequence1 .txt file"
-f2: "Path of sequence2 .txt file"
-n: "Length of fingerprint"
-k: "Length of kmer"
-s: "Length of stride"
-m: "Toggle multi hash layers" (optional)
```

Here is an example command to compare files `file1` and `file2` using 128 hash functions, k-mer size 16, and stride length of 1.

```python minhash_combined.py -f1=file1 -f2=file2 -method=khash -n=128 -k=16 -s=1```

## Needleman-Wunsch
'run_DP.py' is used to compute the true edit distance. The script calls 'edit_distance_DP.py'

## Tuning Parameters
'tune_parameters.py' is used to find the optimal value of k and L. The script iterates through multiple combinations of k-mer lengths and size of fingerprint, and computes the Jaccard similarity for each of the MinHash implementation. It outputs Edit Distance vs. True Jaccard similiarity and the average error for Jaccard similarity for each of the different MinHash implementations. 

## Running Experiments
The 'run_experiment.py' is used to run the different experiments described in our paper.
