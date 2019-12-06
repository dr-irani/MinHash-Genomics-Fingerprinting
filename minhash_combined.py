import random
import xxhash
import numpy as np
import argparse
'''
Possible areas for improvement:
	- Loop ordering: Get k hash values at every kmer? (in other words, one iteration through read)
Stretch goals:
- If we aren't pushing time, maybe try multithreaded minhashing? :D
	I.e., if we have m threads and k hash functions, 
	each thread gets minhash for k/m functions
'''

def create_k_mer_set(seq, kmer_len, stride_len):
	'''
	Return set of unique k-mers inside sequence.
	Stride length between kmers given by stride_len.
	'''
	kmer_set = set()
	i = 0
	while i + kmer_len < n: 
		kmer = seq[i:i+kmer_len]
		kmer_set.add(kmer)
		i += stride_len

	return kmer_set

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


def min_hash(set1, set2, num_hash, method):
	'''
	Return MinHash fingerprints of length num_hash of the two input sequences.
	3 methods: khash, bottomk, kpartition
	Each hash key is a kmer of length kmer_len.
	'''
	fingerprint1 = [0]*num_hash
	fingerprint2 = [0]*num_hash

	if method == "khash":
		# MinHash with k hash functions
		
		hash_fxns = gen_k_hash_functions(num_hash)
		for hash_index, h in enumerate(hash_fxns): 
			min_hval = float('inf')
			for kmer in set1:
				curr_hval = apply_hash(h, kmer)
				min_hval = min(min_hval, curr_hval)
			fingerprint1[hash_index] = min_hval

			min_hval = float('inf')
			for kmer in set2:
				curr_hval = apply_hash(h, kmer)
				min_hval = min(min_hval, curr_hval)
			fingerprint2[hash_index] = min_hval

	if method == "bottomk":
		# MinHash with bottom k hashes

		hash_fxns = gen_k_hash_functions(1)
		h = hash_fxns[0]
		
		hashes = []
		for kmer in set1:
			hashes.append(apply_hash(h, kmer))
		hashes.sort()
		fingerprint1 = hashes[:num_hash]

		hashes = []
		for kmer in set2:
			hashes.append(apply_hash(h, kmer))
		hashes.sort()
		fingerprint2 = hashes[:num_hash]

	if method == "kpartition":
		# MinHash with k partitions

		hash_fxns = gen_k_hash_functions(1)
		h = hash_fxns[0]

		# Literature suggested using the first few bits...but I'm going to try mod
		for kmer in set1:
			hval = apply_hash(h, kmer)
			i = hval % k
			fingerprint1[i] = min(fingerprint1[i], hval)

		for kmer in set2:
			hval = apply_hash(h, kmer)
			i = hval % k
			fingerprint1[i] = min(fingerprint2[i], hval)

	return fingerprint1, fingerprint2


def calculate_jaccard(n, m, k, f1, f2):
	if n == m:
		return np.sum(f1==f2) / k

if __name__ == '__main__':
	# TODO: Handle more elegant arg passing for hyperparameters
	parser = argparse.ArgumentParser()
	parser.add_argument("-f1", help="Path of sequence1 .txt file", required=True)
	parser.add_argument("-f2", help="Path of sequence2 .txt file", required=True)
	parser.add_argument("-n", help="Length of fingerprint", required=True, type=int)
	parser.add_argument("-k", help="Length of kmer", required=True, type=int)
	parser.add_argument("-s", help="Length of stride", required=True, type=int)
	args = parser.parse_args()

	# Place hyperparameters here
	num_hash = args.n
	kmer_len = args.k 
	stride_len = args.s

	fh1 = args.f1
	fh2 = args.f2

	f = open(fh1, 'r')
	seq = f.read()
	n = len(seq)
	f.close
	set1 = create_k_mer_set(seq, kmer_len, stride_len)

	f = open(fh1, 'r')
	seq = f.read()
	m = len(seq)
	f.close
	set2 = create_k_mer_set(seq, kmer_len, stride_len)

	fp1, fp2 = min_hash(set1, set2, num_hash, method='khash')
	jaccard = calculate_jaccard(n, m, num_hash, fp1, fp2)
	est_edit_dist = int(jaccard / n)

	print(fp1)
	print(fp2)
	print(np.count_nonzero(fp1 == fp2))
	print(est_edit_dist)
    