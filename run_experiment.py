import cProfile
import os
import sys
import pstats
from minhash_combined import create_k_mer_set, min_hash, calculate_jaccard, estimate_edit_distance
from edit_distance_DP import edDistDp

reads = []
filename = 'data/random_10000bp_50seq.txt'
with open(filename) as f:
	seq = f.readline()
	reads.append(seq)

kmer_list = [10, 20, 30, 40, 50]
sketch_list = [10, 50, 100, 500, 1000]

for kmer_len in kmer_list:
	for num_hash in sketch_list:

		for i in range(len(reads)):
			for j in range(i+1, len(reads)):
				test(reads[i], reads[j], kmer_len, num_hash)
		

def test(seq1, seq2, kmer_len, num_hash):
    """
    TODO: Call methods you want to profile here
    """

    true_edit_dist = edDistDp(seq1, seq2)

    stride_len = 1
    set1 = create_k_mer_set(seq1, kmer_len, stride_len)
    set2 = create_k_mer_set(seq1, kmer_len, stride_len)

    fp1, hash_fxns = min_hash(set1, num_hash, method='khash')
	fp2, hash_fxns = min_hash(set2, num_hash, method='khash', hash_fxns=hash_fxns)
	jaccard = calculate_jaccard(num_hash, fp1, fp2)
	est_edit_dist = estimate_edit_distance(jaccard, len_x, len_y)

	fp1, hash_fxns = min_hash(set1, num_hash, method='bottomk')
	fp2, hash_fxns = min_hash(set2, num_hash, method='bottomk', hash_fxns=hash_fxns)
	jaccard = calculate_jaccard(num_hash, fp1, fp2)
	est_edit_dist = estimate_edit_distance(jaccard, len_x, len_y)

	fp1, hash_fxns = min_hash(set1, num_hash, method='kpartition')
	fp2, hash_fxns = min_hash(set2, num_hash, method='kpartition', hash_fxns=hash_fxns)
	jaccard = calculate_jaccard(num_hash, fp1, fp2)
	est_edit_dist = estimate_edit_distance(jaccard, len_x, len_y)

cProfile.run('test()', 'OUTFILE')
p = pstats.Stats('OUTFILE')
p.sort_stats('tottime')
p.print_stats()