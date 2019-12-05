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


def weighted_edit_distance(x, y, m, g):

	'''
	Calculate weighted edit distance between sequences x and y
	using matrix dynamic programming. Return weighted edit
	distance. For edit distance: m = 1, g = 1.
	'''
	D = np.zeros((len(x)+1, len(y)+1), dtype=int)
	D[0, 1:] = range(1, len(y)+1)
	D[1:, 0] = range(1, len(x)+1)

	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			delt = m if x[i-1] != y[j-1] else 0
			D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+g, D[i, j-1]+g)

	return D[len(x), len(y)]
