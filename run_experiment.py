#import cProfile
import os
import sys
import pstats
from minhash_combined import create_k_mer_set, min_hash, calculate_jaccard, estimate_edit_distance
from edit_distance_DP import edDistDp


def test(i, j, f, seq1, seq2, kmer_len, num_hash):
	"""
	TODO: Call methods you want to profile here
	"""

	true_edit_dist = edDistDp(seq1, seq2)

	stride_len = 1
	len_x = len(seq1)
	len_y = len(seq2)
	set1 = create_k_mer_set(seq1, kmer_len, stride_len)
	set2 = create_k_mer_set(seq1, kmer_len, stride_len)

	fp1, hash_fxns = min_hash(set1, num_hash, method='khash')
	fp2, hash_fxns = min_hash(set2, num_hash, method='khash', hash_fxns=hash_fxns)
	jaccard1 = calculate_jaccard(num_hash, fp1, fp2)
	est_edit_dist1 = estimate_edit_distance(jaccard1, len_x, len_y)

	fp1, hash_fxns = min_hash(set1, num_hash, method='bottomk')
	fp2, hash_fxns = min_hash(set2, num_hash, method='bottomk', hash_fxns=hash_fxns)
	jaccard2 = calculate_jaccard(num_hash, fp1, fp2)
	est_edit_dist2 = estimate_edit_distance(jaccard2, len_x, len_y)

	fp1, hash_fxns = min_hash(set1, num_hash, method='kpartition')
	fp2, hash_fxns = min_hash(set2, num_hash, method='kpartition', hash_fxns=hash_fxns)
	jaccard3 = calculate_jaccard(num_hash, fp1, fp2)
	est_edit_dist3 = estimate_edit_distance(jaccard3, len_x, len_y)

	f.write("(%d, %d)\t%d\t%f\t%d\t%f\t%f\t%d\t%f\t%d" % (i, j, true_edit_dist[0], jaccard1, est_edit_dist1[0], jaccard2, est_edit_dist2[0], jaccard3, est_edit_dist3[0]))

def main():

	reads = []
	filename = 'data/random_10000bp_50seq.txt'
	with open(filename) as f:
		seq = f.readline()
		while seq:
			reads.append(seq)
			seq = f.readline()

	kmer_list = [10, 20, 30, 40, 50]
	sketch_list = [10, 50, 100, 500, 1000]

	f = open("output/random_10000bp_50seq_output.txt", 'w')
	header = ["pairID", "ED", "J_khash", "ED_khash", "J_bottomk", "ED_bottomk", "J_kpart", "ED_kpart"]
	f.write('\t'.join([str(x) for x in header]))
	#import pdb; breakpoint()
	for kmer_len in kmer_list:
		for num_hash in sketch_list:
			for i in range(len(reads)):
				for j in range(i+1, len(reads)):
					print("read")
					test(i, j, f, reads[i], reads[j], kmer_len, num_hash)

if __name__ == '__main__':
	main()
	

#cProfile.run('test()', 'OUTFILE')
#p = pstats.Stats('OUTFILE')
#p.sort_stats('tottime')
#p.print_stats()