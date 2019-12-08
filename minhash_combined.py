import random
import xxhash
import numpy as np
import argparse
import time
import math
from collections import namedtuple
from memory_profiler import profile

class BloomFilter(object):
	def __init__(self, kmers, fp_prob=0.01):
		self.fp_prob = fp_prob
		self.num_kmers = kmers
		self.size = self.calculate_size(kmers)
		self.num_hash = self.calculate_num_hashes(kmers)
		self.filter = np.zeros(self.size)

	def calculate_size(self, kmers):
		return int(-(kmers * np.log(self.fp_prob)) / np.log(2) ** 2)

	def calculate_num_hashes(self, kmers):
		return int(np.log(2) * self.size / kmers)

	def hash(self, kmer, i):
		return xxhash.xxh32(kmer, seed=i).intdigest() % self.size

	def add(self, kmer):
		for i in range(self.num_hash):
			idx = self.hash(kmer, i)
			self.filter[idx] = 1

	def find(self, kmer):
		for i in range(self.num_hash):
			idx = self.hash(kmer, i)
			if self.filter[idx] == 0: return False
		return True

	def union(self, kmers):
		union_filter = np.zeros(self.size)
		for kmer in kmers:
			for i in range(self.num_hash):
				idx = self.hash(kmer, i)
				union_filter[idx] = 1
		return np.bitwise_or(self.filter.astype(np.int32), union_filter.astype(np.int32))


def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f1", help="Path of sequence1 .txt file", required=True)
	parser.add_argument("-f2", help="Path of sequence2 .txt file", required=True)
	parser.add_argument("-n", help="Length of fingerprint", required=True, type=int)
	parser.add_argument("-k", help="Length of kmer", required=True, type=int)
	parser.add_argument("-s", help="Length of stride", required=True, type=int)
	parser.add_argument("-m", help="Use multi hash layers", required=False, default=False, type=bool)
	parser.add_argument("-method", help="Specify which MinHash implementation to use. Options include \"khash\", \"bottomk\", \"kpartition\", \"containment\"")

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


@profile
def min_hash(seqset, num_hash, method, hash_fxns=None):
	'''
	Return MinHash fingerprint the input sequence and the hash function(s) used.
	3 methods: khash, bottomk, kpartition
	Each hash key is a kmer of length kmer_len.
	'''
	fingerprint = [0]*num_hash
	hash_fxns = gen_k_hash_functions(1) if hash_fxns == None else hash_fxns

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
		h = hash_fxns[0]		
		hashes = []
		for kmer in seqset:
			hashes.append(apply_hash(h, kmer))
		hashes.sort()
		fingerprint = hashes[:num_hash]

	if method == "kpartition":
		# MinHash with k partitions
		# Literature suggested using the first few bits
		num_bit = int(math.log2(num_hash))
		h = hash_fxns[0]
		for kmer in seqset:
			hval = apply_hash(h, kmer)
			binary = bin(hval)
			i = int(binary[2:num_bit+2], 2)
			fingerprint[i] = min(fingerprint[i], hval)
	return fingerprint, hash_fxns

@profile
def containment_min_hash(seqset, bloom_filter=None):
	if bloom_filter is None:
		bloom_filter = BloomFilter(len(seqset))
		for kmer in seqset:
			bloom_filter.add(kmer)
		return bloom_filter

	union = -bloom_filter.size / bloom_filter.num_hash * np.log(1 - np.count_nonzero(bloom_filter.union(seqset)) / bloom_filter.size)
	intersection = len(seqset) + bloom_filter.num_kmers - union
	return intersection / union


def hash_kmers(seqset, seed=0):
	return set([xxhash.xxh32(kmer, seed=seed).hexdigest() for kmer in seqset])

def calculate_true_jaccard(s1, s2):
	union = s1.union(s2)
	inter = s1.intersection(s2)

	return len(inter)/len(union)

def calculate_jaccard(k, f1, f2):
	''' Calculate Jaccard similarity '''
	s1 = set(f1)
	s2 = set(f2)
	union = list(s1.union(s2))
	union = set(union[:k])
	inter = s1.intersection(s2, union)

	return len(inter)/k

def estimate_edit_distance(jaccard, len_x, len_y, kmer_len):
	'''
	Return 3 different estimates for edit distance.
	2 values from the bounds of the IEEE paper.
	1 value from our equation.
	'''
	edit_distance = [0] * 3
	max_len = max(len_x, len_y)
	min_len = min(len_x, len_y)
	
	alpha = min_len/max_len
	edit_distance[0] = (1 - alpha) * max_len
	edit_distance[1] = (1 + alpha) * (jaccard/(2-jaccard)) * max_len
	edit_distance[2] = (max_len - jaccard*min_len) / (kmer_len * (jaccard + 1))

	return edit_distance

def mash_distance(jaccard, kmer_len):
	# this seems to be pretty different from the actual edit distance
	# correlation yes, but very different scale
	return (-1/kmer_len) * np.log(2 * jaccard / (1 + jaccard))

@profile
def main():
	SeqSet = namedtuple('SeqSet', 'set len')
	args = get_args()
	num_hash = args.n
	kmer_len = args.k 
	stride_len = args.s
	method = args.methods

	with open(args.f1, 'r') as f:
		seq = f.read()
		set1 = SeqSet(create_k_mer_set(seq, kmer_len, stride_len), len(seq))

	with open(args.f2, 'r') as f:
		seq = f.read()
		set2 = SeqSet(create_k_mer_set(seq, kmer_len, stride_len), len(seq))

	start_time = time.time()
	sets = sorted([set1, set2], key=lambda x: len(x.set), reverse=True)
	set1 = hash_kmers(sets[0].set) if args.m else sets[0].set
	set2 = hash_kmers(sets[1].set) if args.m else sets[1].set

	if method == 'containment':
		bloom_filter = containment_min_hash(set1)
		jaccard = containment_min_hash(set2, bloom_filter=bloom_filter)
	else:
		fp1, hash_fxns = min_hash(set1, num_hash, method=method)
		fp2, hash_fxns = min_hash(set2, num_hash, method=method, hash_fxns=hash_fxns)
		jaccard = calculate_jaccard(num_hash, fp1, fp2)
	mash_dist = mash_distance(jaccard, kmer_len)


if __name__ == '__main__':
	main()
	
