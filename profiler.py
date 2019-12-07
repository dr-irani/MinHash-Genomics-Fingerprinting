import cProfile
import os
import sys
import pstats
from collections import namedtuple
from minhash_combined import create_k_mer_set, containment_min_hash, calculate_jaccard, mash_distance, min_hash

SeqSet = namedtuple('SeqSet', 'set len')
num_hash = 128
kmer_len = 16
stride_len = 1

method = 'containment'
f1 = 'data/synth_1.txt'
f2 = 'data/synth_1_edited.txt'

with open(f1, 'r') as f:
    seq = f.read()
    set1 = SeqSet(create_k_mer_set(seq, kmer_len, stride_len), len(seq))

with open(f2, 'r') as f:
    seq = f.read()
    set2 = SeqSet(create_k_mer_set(seq, kmer_len, stride_len), len(seq))

def test():
    sets = sorted([set1, set2], key=lambda x: len(x.set), reverse=True)
    a = sets[0].set
    b = sets[1].set
    if method == 'containment':
        bloom_filter = containment_min_hash(a)
        jaccard = containment_min_hash(b, bloom_filter=bloom_filter)
    else:
        fp1, hash_fxns = min_hash(a, num_hash, method=method)
        fp2, hash_fxns = min_hash(b, num_hash, method=method, hash_fxns=hash_fxns)
        jaccard = calculate_jaccard(num_hash, fp1, fp2)
    mash_dist = mash_distance(jaccard, kmer_len)


cProfile.run('test()', 'OUTFILE')
p = pstats.Stats('OUTFILE')
p.sort_stats('tottime')
p.print_stats()