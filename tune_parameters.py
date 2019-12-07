import os
import sys
from minhash_combined import create_k_mer_set, min_hash, calculate_jaccard, calculate_true_jaccard, containment_min_hash, estimate_edit_distance
from collections import namedtuple
import glob
import numpy as np

'''
Tune kmer_len and num_hash
'''

def get_jaccard(seq1, seq2, kmer_len, num_hash):


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

	fp1, hash_fxns = min_hash(set1.set, num_hash, method='bottomk')
	fp2, hash_fxns = min_hash(set2.set, num_hash, method='bottomk', hash_fxns=hash_fxns)
	jaccard2 = calculate_jaccard(num_hash, fp1, fp2)

	fp1, hash_fxns = min_hash(set1.set, num_hash, method='kpartition')
	fp2, hash_fxns = min_hash(set2.set, num_hash, method='kpartition', hash_fxns=hash_fxns)
	jaccard3 = calculate_jaccard(num_hash, fp1, fp2)

	sets = sorted([set1, set2], key=lambda x: len(x.set), reverse=True)
	bloom_filter = containment_min_hash(sets[0].set)
	jaccard4 = containment_min_hash(sets[1].set, bloom_filter=bloom_filter)

	output = '{} \t {} \t {} \t {} \t {} \t {} \t {} \t'.format(len_x, len_y, truejaccard, jaccard1, jaccard2, jaccard3, jaccard4)
	return output


def tune(reads1, reads2, kmer_len, num_hash, output_file):

	f = open(output_file, 'w')
	for i in range(len(reads1)):
		seq1 = reads1[i]
		seq2 = reads2[i]
		output = get_jaccard(seq1, seq2, kmer_len, num_hash)		
		f.write(output + '\n')

	f.close()

def create_readlists_ecoli():
	'''
	Create reads1, reads2 that contains the sequence pairs for Ecoli data
	'''
	filenames = 'data/ecoli_realdata.txt'
	f = open(filenames, 'r')

	reads1 = []
	reads2 = []

	for i in range(0, 600, 2):
		seq1 = f.readline()
		seq2 = f.readline()
		reads1.append(seq1)
		reads2.append(seq2)

	f.close()
		
	return reads1, reads2

def create_readlists():
	'''
	Create reads1, reads2 that contains the sequence pairs for simulated data
	'''
	filenames1 = ['data/synth_{}.txt'.format(i) for i in range(1,11)]
	filenames2 = ['data/synth_{}_edited.txt'.format(i) for i in range(1, 11)]

	reads1 = []
	reads2 = []

	for i in range(1, 11):
		fn1 = filenames1[i-1]
		fn2 = filenames2[i-1]

		with open(fn1) as f:
			seq1 = f.readline()
			reads1.append(seq1)
		with open(fn2) as f:
			seq2 = f.readline()
			reads2.append(seq2)

	filenames1 = ['synth_data/synth_{}.txt'.format(i) for i in range(1,11)]		
	for i in range(1, 11):
		
		temp = 'synth_data/synth_{}_edited_*.txt'.format(i)
		filenames2 = glob.glob(temp)

		fn1 = filenames1[i-1]
		with open(fn1) as f:
			seq1 = f.readline()

		for fn2 in filenames2:
			with open(fn2) as f:
				seq2 = f.readline()
			reads1.append(seq1)
			reads2.append(seq2)
		
	return reads1, reads2

def summary(kmer_len, num_hash):
	'''
	Summarize data and output for visualization.
	'''
	hash_file = 'tuning/synth/{}_mer_{}_numhash.txt'.format(kmer_len, num_hash)
	#hash_file = 'tuning/ecoli/{}_mer_{}_numhash.txt'.format(kmer_len, num_hash)

	trueJ = []
	i = 0
	err1 = 0
	err2 = 0
	err3 = 0
	err4 = 0
	with open(hash_file) as f:
		read = f.readline().split()
		#len_x = int(read[0])
		#len_y = int(read[1])

		TJ = float(read[2])
		err1 += np.abs(TJ - float(read[3])) # Cumulative error for L-Hash
		err2 += np.abs(TJ - float(read[4])) # Cumulative error for Bottom L
		err3 += np.abs(TJ - float(read[5])) # Cumulative error for L-partition
		err4 += np.abs(TJ - float(read[6])) # Cumulative error for Containment Hash
		i += 1
		trueJ.append(TJ)

	err1 = err1/i
	err2 = err2/i
	err3 = err3/i
	err4 = err4/i

	return trueJ, err1, err2, err3, err4

def main_summary():
	''' 
	Run summary() for each kmer, num_hash pair
	'''

	trueED_file = 'output/synthdata_true_editdistance_single.txt'
	#trueED_file = 'output/ecoli_true_editdistance_single.txt'

	with open(trueED_file) as f:
		trueED = [int(line) for line in iter(f.readline, r'')]

	kmer_list = [4, 8, 12, 16, 20]
	numhash_list = [32, 64, 128, 256, 512]

	trueJ_comb = np.zeros((len(trueED), 26))
	print(len(trueED))
	trueJ_comb[:, 0] = trueED
	err1 = [0] * 25
	err2 = [0] * 25
	err3 = [0] * 25
	err4 = [0] * 25
	i = 0
	for kmer_len in kmer_list:
		for num_hash in numhash_list:
			trueJ, err1[i], err2[i], err3[i], err4[i] = summary(kmer_len, num_hash)
			trueJ_comb[:, i+1] = np.asarray(trueJ).reshape((-1, 1))
			i += 1

	err1 = np.asarray(err1).reshape((-1, 1))
	err2 = np.asarray(err2).reshape((-1, 1))
	err3 = np.asarray(err3).reshape((-1, 1))
	err4 = np.asarray(err4).reshape((-1, 1))
	err = np.concatenate((err1, err2, err3, err4), axis = 1)

	np.savetxt("synth_trueJ.csv", trueJ_comb, delimiter = ',')
	np.savetxt("synth_hasherror.csv", err, delimiter = ',')

def main():

	#reads1, reads2 = create_readlists()
	reads1, reads2 = create_readlists_ecoli()

	kmer_list = [4, 8, 12, 16, 20]
	numhash_list = [32, 64, 128, 256, 512]

	for kmer_len in kmer_list:
		for num_hash in numhash_list:
			output_file = 'tuning/{}_mer_{}_numhash.txt'.format(kmer_len, num_hash)
			tune(reads1, reads2, kmer_len, num_hash, output_file)
			print("Finished %d mer %d numhash" % (kmer_len, num_hash))

if __name__ == '__main__':
	main_summary()