import random
import xxhash
import numpy as np
import argparse
import time


class BloomFilter(object):
	def __init__(self, kmers, fp_prob=0.01):
		self.fp_prob = fp_prob
		self.num_kmers = kmers
		self.size = self.calculate_size()
		self.num_hash = self.calculate_num_hashes()
		self.filter = np.zeros(self.size)

	def calculate_size(self):
		return int(-(self.num_kmers * np.log(self.fp_prob)) / np.log(2) ** 2)

	def calculate_num_hashes(self):
		return int(self.size / self.num_kmers * np.log(2))

	def hash(self, kmer, i):
		return xxhash.xxh32(kmer, seed=i).intdigest() % self.size

	def add(self, kmer):
		for i in range(self.num_hash):
			idx = self.hash(kmer, i)
			self.filter[idx] == 1

	def find(self, kmer):
		for i in range(self.num_hash):
			idx = self.hash(kmer, i)
			if self.filter[idx] == 0: return False
		return True


def get_args():
	parser = argparse.ArgumentParser()
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


def min_hash(seqset, num_hash, method, hash_fxns=None, bloom_filter=None):
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
			i = hval % k
			fingerprint[i] = min(fingerprint[i], hval)
	return fingerprint, hash_fxns


def containment_min_hash(seqset, num_hash, bloom_filter=None):
	if bloom_filter is None:
		bloom_filter = BloomFilter(len(seqset))
		map(bloom_filter.add, seqset)
		return bloom_filter
		# for kmer in seqset:
		# 	bloom_filter.add(kmer)
	containment = 0
	for kmer in seqset:
		containment += 1 if bloom_filter.find(kmer) else 0
	containment -= np.round(bloom_filter.fp_prob * num_hash)
	containment /= num_hash
	return len(seqset) * containment / (len(seqset) + bloom_filter.num_kmers - len(seqset) * containment)



def calculate_containment_coeff(c, fp, num_hash, seqset):
	c -= np.round(fp_prob * num_hash)
	c = c / num_hash
	return len(seqset) * c / ()


def calculate_jaccard(n, m, k, f1, f2):
	if n == m:
		return np.sum(f1==f2) / k

def main():
	args = get_args()
	# Place hyperparameters here
	num_hash = args.n
	kmer_len = args.k 
	stride_len = args.s
	method = 'containment'

	with open(args.f1, 'r') as f:
		seq = f.read()
		n = len(seq)
		set1 = create_k_mer_set(seq, kmer_len, stride_len)

	with open(args.f2, 'r') as f:
		seq = f.read()
		m = len(seq)
		set2 = create_k_mer_set(seq, kmer_len, stride_len)


	start_time = time.time()
	if method == 'khash':
		sets = [set1, set2]
		sets.sort(key=len, reverse=True)
		bloom_filter = containment_min_hash(sets[0], num_hash)
		jaccard = containment_min_hash(sets[1], num_hash, bloom_filter=bloom_filter)
	else:
		fp1, hash_fxns = min_hash(set1, num_hash, method=method)
		fp2, hash_fxns = min_hash(set2, num_hash, method=method, hash_fxns=hash_fxns)
		jaccard = calculate_jaccard(n, m, num_hash, fp1, fp2)
		print(fp1)
		print(fp2)
		print(np.count_nonzero(fp1 == fp2))

	est_edit_dist = int(jaccard * n)
	elapsed_time = time.time() - start_time

	
	print(est_edit_dist)
	print(elapsed_time)


if __name__ == '__main__':
	main()
