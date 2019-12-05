import random
import xxhash
import numpy as np
'''
Possible areas for improvement:
	- Loop ordering: Get k hash values at every kmer? (in other words, one iteration through read)
Stretch goals:
- If we aren't pushing time, maybe try multithreaded minhashing? :D
	I.e., if we have m threads and k hash functions, 
	each thread gets minhash for k/m functions
'''


def apply_hash(h, key):
	''' 
	Return h(key). 
	This function is a wrapper for xxhash functions with initialized seeds.
	Currently assume h is a xxhash.x32 object with initialized seed
	If we change choice of hash function later, it will be easier to change
	how we apply the hash (either through a function or an object) in this method
	'''
	h.update(key)
	val = h.intdigest() # TODO: What representation to return? (hex in str format?)
	h.reset()
	return val


def gen_k_hash_functions(k):
	'''
	Return k random 32-bit seeds for xxh32 (maybe experiment with 64 bit)
	Currently returns xxhash.x32 objects with initialized seeds.
	'''
	random_seeds = random.sample(range(0xfffffff),k) 
	return [xxhash.xxh32(seed=seed) for seed in random_seeds]


def min_hash(seq, hash_fxns, kmer_len, stride_len):
	'''
	Return minhash fingerprint of seq using k hash_fxns.
	Each hash key is a kmer of length kmer_len. 
	Stride length between kmers given by stride_len.
	'''
	fingerprint = [0]*len(hash_fxns)
	n = len(seq)
	for hash_index, h in enumerate(hash_fxns): 
		min_hval = float('inf')
		i = 0
		while i + kmer_len < n: 
			kmer = seq[i:i+kmer_len]
			curr_hval = apply_hash(h, kmer)
			min_hval = min(min_hval, curr_hval)
			i+=stride_len
		fingerprint[hash_index] = min_hval
	return fingerprint


def calculate_jaccard(n, m, k, f1, f2):
	if n == m:
		return np.sum(f1==f2) / k

if __name__ == '__main__':
	# TODO: Handle more elegant arg passing for hyperparameters
	# Place hyperparameters here
	num_hash_fxns = 13
	kmer_len = 20 
	stride_len = 10

	hash_fxns = gen_k_hash_functions(13)

	f = open('data/test.txt', 'r')
	seq = f.read()
	f.close()
	orig_fingerprint = min_hash(seq, hash_fxns, kmer_len, stride_len)
	# print(fingerprint)

	f = open('data/test.txt_edited_0i_0d_100s.txt', 'r')
	edit_seq = f.read()
	f.close()
	edit_fingerprint = min_hash(seq, hash_fxns, kmer_len, stride_len)


	print(orig_fingerprint)
	print(edit_fingerprint)
	print(np.count_nonzero(orig_fingerprint == edit_fingerprint))
	n,m = len(seq), len(edit_seq) 
	jaccard = calculate_jaccard(n,m,num_hash_fxns,orig_fingerprint,edit_fingerprint)
	est_edit_dist = int(jaccard / n)
	print(est_edit_dist)
    