import math
import numpy as np

def k_partition(seq, k, kmer_len, stride_len):

	'''
	Return k-partition minhash sketch fingerprint of seq.
	'''
	n = len(seq)
	p_size = math.floor(n/k)
	hashes = []

	for i in range(k):
		j = i * p_size
		min_value = float('inf')
		while j < (i+1) * p_size and j < n:
			if hash(seq[j:j+kmer_len]) < min_value:
				min_value = hash(seq[j:j+kmer_len])
				j += stride_len
		hashes.append(min_value)

	return hashes

def bottom_k(seq, k, kmer_len, stride_len):
	'''
	Return bottom-k minhash sketch fingerprint of seq.
	k-mers of kmer_len extracted from seq by taking stride
	lengths of stride_len.
	'''
	i = 0
	hashes = []
	n = len(seq)

	while i + kmer_len < n:
		hashes.append(hash(seq[i:i+kmer_len]))
		i += stride_len

	hashes.sort()

	return hashes[:k]

def estimate_edit_distance(jaccard, len_x, len_y):
	if len_x == len_y:
		return [jaccard * len_x]
	else:
		alpha = min(len_x, len_y) / max(len_x, len_y)
		return [1 - alpha, (1+alpha) * (jaccard/(2-jaccard))]

def mash_distance(jaccard, kmer_len):
	return (-1/kmer_len) * np.log(2 * jaccard / (1 + jaccard))
