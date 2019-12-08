import math
import numpy as np

def k_partition(seq, k, kmer_len, stride_len):
	"""
    Return k-partition minhash sketch fingerprint of sequence.

    Parameters
    ----------
    seq : string sequence to perform bottom-k sketch on
    k : number of hash value representatives in k-partition sketch
	kmer_len : int the length of the k-mers of the sequence
	stride_len : int the stride length in extracting k-mers from the sequence

    Returns
    -------
    fingerprint : list of ints
        The k-partition fingerprint sketch of the sequence.
    """
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

def bottom_k(seq, k, kmer_len, stride_len=1):
	"""
    Return bottom-k minhash sketch fingerprint of sequence.

    Parameters
    ----------
    seq : string sequence to perform bottom-k sketch on
    k : number of hash value representatives in bottom-k sketch
	kmer_len : int the length of the k-mers of the sequence
	stride_len : int the stride length in extracting k-mers from the sequence

    Returns
    -------
    fingerprint : list of ints
        The bottom-k fingerprint sketch of the sequence.
    """
	i = 0
	hashes = []
	n = len(seq)

	while i + kmer_len < n:
		hashes.append(hash(seq[i:i+kmer_len]))
		i += stride_len

	hashes.sort()

	return hashes[:k]

def estimate_edit_distance(jaccard, len_x, len_y):
	"""
    Return the estimated edit distance between two seequences.

    Parameters
    ----------
    jaccard : double/float Jaccard similarity metric between the two sequences
    len_x : int length of sequence x
    len_y : int length of sequence y

    Returns
    -------
    edit_distance : list of int/float/double
        The estimated edit distance estimate or lower bound and upper bound on
        edit distance.
    """
	if len_x == len_y:
		return [jaccard * len_x]
	else:
		alpha = min(len_x, len_y) / max(len_x, len_y)
		return [1 - alpha, (1+alpha) * (jaccard/(2-jaccard))]

def mash_distance(jaccard, kmer_len):
	"""
    Calculates the mash distance between two sets of k-mers.

    Parameters
    ----------
    jaccard : double/float Jaccard similarity metric between the two sequences
	kmer_len : int the length of the k-mers of the sequence

    Returns
    -------
    mash distance : double/float
        The mash distance between two sets of k-mers.
    """
	return (-1/kmer_len) * np.log(2 * jaccard / (1 + jaccard))
