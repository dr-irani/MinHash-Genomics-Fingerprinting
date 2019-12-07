#import cProfile
import os
import sys
import pstats
from minhash_combined import create_k_mer_set, min_hash, calculate_jaccard, calculate_true_jaccard, containment_min_hash, estimate_edit_distance
from edit_distance_DP import edDistDp
from collections import namedtuple
import glob

def test(seq1, seq2, kmer_len, num_hash):

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

def ED_main():
	'''
	Combine all outputs (jaccard similarities and true edit distance) from our synthetic data.
	Calculate the estimated edit distance and output result in single file
	'''
	kmer_len = 16
	num_hash = 128

	hash_files = ['output/synth_{}_hash.txt'.format(i) for i in range(2, 11)]
	trueED_files = ['output/synth_{}_editdistance.txt'.format(i) for i in range(2, 11)]
	for i in range(1, 11):
		for j in range(5):
			hash_files.append('synth_output/synth_{}_{}_hash.txt'.format(i, j))
			trueED_files.append('synth_output/synth_{}_{}_editdistance.txt'.format(i, j))

	output_file = 'output/synthdata_estimated_editdistance.txt'
	fo = open(output_file, 'w')
	output_header = ['len_x', 'len_y', 'TJ', 'J1', 'J2', 'J3', 'J4', 'trueED', 
	'TED[0]', 'TED[1]', 'TED[2]', 'ED1[0]', 'ED1[1]', 'ED1[2]', 'ED2[0]', 'ED2[1]', 'ED2[2]', 
	'ED3[0]', 'ED3[1]', 'ED3[2]', 'ED4[0]', 'ED4[1]', 'ED4[2]']
	output = '\t'.join(output_header)
	fo.write(output + '\n')

	for i in range(len(hash_files)):
		hf_name = hash_files[i]
		ed_name = trueED_files[i]

		with open(hf_name) as f: 
			temp = f.readline()
			templist = temp.split()
			[len_x, len_y, TJ, J1, J2, J3, J4] = [float(i) for i in templist]
			len_x = int(len_x)
			len_y = int(len_y)

		with open(ed_name) as f: 
			temp = f.readline()
			templist = temp.split()
			[len_x1, len_y1, trueED] = [float(i) for i in templist]
			len_x1 = int(len_x1)
			len_y1 = int(len_y1)
			assert len_x == len_x1
			assert len_y == len_y1

		TED = estimate_edit_distance(TJ, len_x, len_y, kmer_len)
		ED1 = estimate_edit_distance(J1, len_x, len_y, kmer_len)
		ED2 = estimate_edit_distance(J2, len_x, len_y, kmer_len)
		ED3 = estimate_edit_distance(J3, len_x, len_y, kmer_len)
		ED4 = estimate_edit_distance(J4, len_x, len_y, kmer_len)

		output_list = [len_x, len_y, TJ, J1, J2, J3, J4, trueED, TED[0], TED[1], TED[2], 
		ED1[0], ED1[1], ED1[2], ED2[0], ED2[1], ED2[2], ED3[0], ED3[1], ED3[2], ED4[0], ED4[1], ED4[2]]
		output = '\t'.join(str(e) for e in output_list)
		fo.write(output + '\n')

	fo.close()

def ecoli_main():
	'''
	Get Jaccard similarity for the ecoli read pairs. 
	Calculate the estimated edit distance and output combined result in single file
	'''
	
	# Input and output files
	seq_file = 'data/ecoli_realdata.txt'
	ED_file = 'output/ecoli_true_editdistance.txt'
	hash_output_file = 'output/ecoli_hash.txt'
	ed_output_file = 'output/ecoli_estimated_editdistance.txt'

	# Parameters
	kmer_len = 16
	num_hash = 128

	# Calculate Jaccard similarity for all pairs, output to hash_output_file
	# f = open(seq_file, 'r')
	# fo = open(hash_output_file, 'w')

	# for i in range(0, 600, 2):
	# 	seq1 = f.readline()
	# 	seq2 = f.readline()
	# 	output = test(seq1, seq2, kmer_len, num_hash)
	# 	fo.write(output + '\n')
	# fo.close()
	# f.close()

	# Calculate estimated edit distance for all pairs, output to ed_output_file
	f_ed = open(ED_file, 'r')
	f_h = open(hash_output_file, 'r')

	fo = open(ed_output_file, 'w')
	output_header = ['len_x', 'len_y', 'TJ', 'J1', 'J2', 'J3', 'J4', 'trueED', 
	'TED[0]', 'TED[1]', 'TED[2]', 'ED1[0]', 'ED1[1]', 'ED1[2]', 'ED2[0]', 'ED2[1]', 'ED2[2]', 
	'ED3[0]', 'ED3[1]', 'ED3[2]', 'ED4[0]', 'ED4[1]', 'ED4[2]']
	output = '\t'.join(output_header)
	fo.write(output + '\n')

	f_ed.readline()  # get rid of header
	for i in range(0, 600, 2):
		temp = f_h.readline()
		templist = temp.split()
		[len_x, len_y, TJ, J1, J2, J3, J4] = [float(i) for i in templist]
		len_x = int(len_x)
		len_y = int(len_y)

		temp = f_ed.readline()
		templist = temp.split('\t')
		[len_x1, len_y1, trueED] = [float(i) for i in templist[1:]]
		len_x1 = int(len_x1)
		len_y1 = int(len_y1)
		assert len_x == len_x1
		assert len_y == len_y1

		TED = estimate_edit_distance(TJ, len_x, len_y, kmer_len)
		ED1 = estimate_edit_distance(J1, len_x, len_y, kmer_len)
		ED2 = estimate_edit_distance(J2, len_x, len_y, kmer_len)
		ED3 = estimate_edit_distance(J3, len_x, len_y, kmer_len)
		ED4 = estimate_edit_distance(J4, len_x, len_y, kmer_len)

		output_list = [len_x, len_y, TJ, J1, J2, J3, J4, trueED, TED[0], TED[1], TED[2], 
		ED1[0], ED1[1], ED1[2], ED2[0], ED2[1], ED2[2], ED3[0], ED3[1], ED3[2], ED4[0], ED4[1], ED4[2]]
		output = '\t'.join(str(e) for e in output_list)
		fo.write(output + '\n')

	fo.close()
	f_ed.close()
	f_h.close()

def synth_main():

	# For running the firt 10 files
	# filenames1 = ['data/synth_{}.txt'.format(i) for i in range(1,11)]
	# filenames2 = ['data/synth_{}_edited.txt'.format(i) for i in range(1, 11)]

	# kmer_len = 16
	# num_hash = 128

	# for i in range(1, 11):
	# 	fn1 = filenames1[i-1]
	# 	fn2 = filenames2[i-1]
	# 	with open(fn1) as f:
	# 		seq1 = f.readline()
	# 	with open(fn2) as f:
	# 		seq2 = f.readline()

	# 	fo_name = 'output/synth_{}_hash.txt'.format(i)
	# 	fo = open(fo_name, 'w')
	# 	output = test(seq1, seq2, kmer_len, num_hash)
	# 	fo.write(output)
	# 	print("%d completed" % i)
	# 	fo.close

	filenames1 = ['synth_data/synth_{}.txt'.format(i) for i in range(1,11)]

	kmer_len = 16
	num_hash = 128

	for i in range(1, 11):
		fn1 = filenames1[i-1]
		filenames2 = glob.glob('synth_data/synth_1_edited_*.txt')

		with open(fn1) as f1:
			seq1 = f1.readline()
			for j in range(5):
				fn2 = filenames2[j]
				with open(fn2) as f2:
					seq2 = f2.readline()
				fo_name = 'synth_output/synth_{}_{}_hash.txt'.format(i, j)
				fo = open(fo_name, 'w')
				output = test(seq1, seq2, kmer_len, num_hash)
				fo.write(output)
				print("%d, %d completed" % (i, j))
				fo.close()

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
	#synth_main()
	#ecoli_main()
	ED_main()
#cProfile.run('test()', 'OUTFILE')
#p = pstats.Stats('OUTFILE')
#p.sort_stats('tottime')
#p.print_stats()