import random
import xxhash
import numpy as np
import argparse
import time


def get_args():
	parser = argparse.ArgumentParser()

	# Change according to how you are passing in your input files
	parser.add_argument("-f1", help="Path of sequence1 .txt file", required=True)
	parser.add_argument("-f2", help="Path of sequence2 .txt file", required=True)
	
	parser.add_argument("-n", help="Length of fingerprint", required=True, type=int)
	parser.add_argument("-k", help="Length of kmer", required=True, type=int)
	parser.add_argument("-s", help="Length of stride", required=True, type=int)
	args = parser.parse_args()
	return args


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
	while i + kmer_len < len(seq): 
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


def min_hash(seqset, num_hash, method, hash_fxns = None):
	'''
	Return MinHash fingerprint the input sequence and the hash function(s) used.
	3 methods: khash, bottomk, kpartition
	Each hash key is a kmer of length kmer_len.
	'''
	fingerprint = [0]*num_hash
	if hash_fxns == None:
		hash_fxns = gen_k_hash_functions(1)
	h = hash_fxns[0]


	if method == "khash":
		# MinHash with k hash functions
		
		if len(hash_fxns) != num_hash:	
			hash_fxns = gen_k_hash_functions(num_hash)
		
		for hash_index, h in enumerate(hash_fxns): 
			min_hval = float('inf')
			for kmer in seqset:
				curr_hval = apply_hash(h, kmer)
				min_hval = min(min_hval, curr_hval)
			fingerprint[hash_index] = min_hval

	if method == "bottomk":
		# MinHash with bottom k hashes		
		hashes = []
		for kmer in seqset:
			hashes.append(apply_hash(h, kmer))
		hashes.sort()
		fingerprint = hashes[:num_hash]

	if method == "kpartition":
		# MinHash with k partitions
		# Literature suggested using the first few bits...but I'm going to try mod
		for kmer in seqset:
			hval = apply_hash(h, kmer)
			i = hval % num_hash
			fingerprint[i] = min(fingerprint[i], hval)

	if method == "minimizerSeq":
		# MinHash that stores minimizer sequences
		pass

	return fingerprint, hash_fxns


def calculate_jaccard(k, f1, f2):
	''' Calculate Jaccard similarity '''
	s1 = set(f1)
	s2 = set(f2)
	union = list(s1.union(s2))
	union = set(union[:k])
	inter = s1.intersection(s2, union)

	return len(inter)/k

def estimate_edit_distance(jaccard, len_x, len_y):
	if len_x == len_y:
		return [int(jaccard * len_x)]
	else:
		alpha = min(len_x, len_y) / max(len_x, len_y)
		return [1 - alpha, (1+alpha) * (jaccard/(2-jaccard))]

def mash_distance(jaccard, kmer_len):
	# this seems to be pretty different from the actual edit distance though..
	# correlation yes, but very different scale
	
	return (-1/kmer_len) * np.log(2 * jaccard / (1 + jaccard))


def main():
	args = get_args()
	num_hash = args.n
	kmer_len = args.k 
	stride_len = args.s

	with open(args.f1, 'r') as f:
		seq = f.read()
		n = len(seq)
		set1 = create_k_mer_set(seq, kmer_len, stride_len)

	with open(args.f2, 'r') as f:
		seq = f.read()
		m = len(seq)
		set2 = create_k_mer_set(seq, kmer_len, stride_len)


	#start_time = time.time()
	fp1, hash_fxns = min_hash(set1, num_hash, method='khash')
	fp2, hash_fxns = min_hash(set2, num_hash, method='khash', hash_fxns=hash_fxns)
	jaccard = calculate_jaccard(num_hash, fp1, fp2)
	est_edit_dist = estimate_edit_distance(jaccard, len_x, len_y)
	mash_dist = mash_distance(jaccard, kmer_len)
	#elapsed_time = time.time() - start_time

	print(fp1)
	print(fp2)
	print(np.count_nonzero(fp1 == fp2))
	print(est_edit_dist)
	print(mash_dist)
	print(elapsed_time)

if __name__ == '__main__':
	main()
	
