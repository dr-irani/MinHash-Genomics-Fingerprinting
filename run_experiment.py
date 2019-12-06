#import cProfile
import os
import sys
import pstats
from minhash_combined import create_k_mer_set, min_hash, calculate_jaccard, calculate_true_jaccard, containment_min_hash
from edit_distance_DP import edDistDp
from collections import namedtuple


def test(seq1, seq2, kmer_len, num_hash):
	"""
	TODO: Call methods you want to profile here
	"""

	#true_edit_dist = edDistDp(seq1, seq2)

	stride_len = 1
	SeqSet = namedtuple('SeqSet', 'set len')
	len_x = len(seq1)
	len_y = len(seq2)

	set1 = SeqSet(create_k_mer_set(seq1, kmer_len, stride_len), len(seq1))
	set2 = SeqSet(create_k_mer_set(seq2, kmer_len, stride_len), len(seq2))

	truejaccard = calculate_true_jaccard(set1.set, set2.set)

	fp1, hash_fxns = min_hash(set1.set, num_hash, method='khash')
	fp2, hash_fxns = min_hash(set2.set, num_hash, method='khash', hash_fxns=hash_fxns)
	jaccard1 = calculate_jaccard(num_hash, fp1, fp2)
	#est_edit_dist1 = estimate_edit_distance(jaccard1, len_x, len_y)

	fp1, hash_fxns = min_hash(set1.set, num_hash, method='bottomk')
	fp2, hash_fxns = min_hash(set2.set, num_hash, method='bottomk', hash_fxns=hash_fxns)
	jaccard2 = calculate_jaccard(num_hash, fp1, fp2)
	#est_edit_dist2 = estimate_edit_distance(jaccard2, len_x, len_y)

	fp1, hash_fxns = min_hash(set1.set, num_hash, method='kpartition')
	fp2, hash_fxns = min_hash(set2.set, num_hash, method='kpartition', hash_fxns=hash_fxns)
	jaccard3 = calculate_jaccard(num_hash, fp1, fp2)

	sets = sorted([set1, set2], key=lambda x: len(x.set), reverse=True)
	bloom_filter = containment_min_hash(sets[0].set)
	jaccard4 = containment_min_hash(sets[1].set, bloom_filter=bloom_filter)

	output = '{} \t {} \t {} \t {} \t {} \t {} \t {} \t'.format(len_x, len_y, truejaccard, jaccard1, jaccard2, jaccard3, jaccard4)
	return output
	
	#est_edit_dist3 = estimate_edit_distance(jaccard3, len_x, len_y)

	#f.write("(%d, %d)\t%d\t%f\t%d\t%f\t%f\t%d\t%f\t%d" % (i, j, true_edit_dist[0], jaccard1, est_edit_dist1[0], jaccard2, est_edit_dist2[0], jaccard3, est_edit_dist3[0]))

def synth_main():

	filenames1 = ['data/synth_{}.txt'.format(i) for i in range(1,11)]
	filenames2 = ['data/synth_{}_edited.txt'.format(i) for i in range(1, 11)]

	kmer_len = 16
	num_hash = 128

	for i in range(1, 11):
		fn1 = filenames1[i-1]
		fn2 = filenames2[i-1]
		with open(fn1) as f:
			seq1 = f.readline()
		with open(fn2) as f:
			seq2 = f.readline()

		fo_name = 'output/synth_{}_hash.txt'.format(i)
		fo = open(fo_name, 'w')
		output = test(seq1, seq2, kmer_len, num_hash)
		fo.write(output)
		print("%d completed" % i)
		fo.close

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
	synth_main()
	

#cProfile.run('test()', 'OUTFILE')
#p = pstats.Stats('OUTFILE')
#p.sort_stats('tottime')
#p.print_stats()